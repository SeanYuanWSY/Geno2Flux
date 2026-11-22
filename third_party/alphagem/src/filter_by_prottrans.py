import gc
import hashlib
import importlib.metadata
import json
import os
import tempfile
import tqdm
import cobra
from transformers import T5Tokenizer, T5EncoderModel
import torch
import pickle
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from src import findrefnames
from src.pipeline_contracts import protein_key, select_by_pca_distance


def _is_oom_error(error):
    message = str(error).lower()
    cuda_oom_type = getattr(torch.cuda, "OutOfMemoryError", None)
    if cuda_oom_type is not None and isinstance(error, cuda_oom_type):
        return True
    return (
        isinstance(error, RuntimeError)
        and "out of memory" in message
        and any(
            marker in message
            for marker in (
                "cuda",
                "cublas",
                "cudnn",
                "cufft",
                "cusolver",
                "cusparse",
            )
        )
    )


def _file_record(path):
    path = os.fspath(path)
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": os.path.getsize(path), "sha256": digest.hexdigest()}


def _write_json_atomic(path, payload):
    path = os.path.abspath(os.fspath(path))
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=os.path.dirname(path),
            prefix=".{}.".format(os.path.basename(path)),
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = handle.name
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


def load_model(
    model_path=None,
    local_files_only=False,
    *,
    force_cpu=False,
    return_audit=False,
):
    use_cuda = torch.cuda.is_available() and not force_cpu
    device = torch.device('cuda:0' if use_cuda else 'cpu')
    print(f"Loading model on {device}...")
    model_source = model_path or "Rostlab/prot_t5_xl_uniref50"
    audit = {
        "fallback_reason": None,
        "final_device": device.type,
        "final_dtype": "float16" if use_cuda else "float32",
        "model_source": str(model_source),
        "requested_device": "cpu" if force_cpu else (
            "cuda" if torch.cuda.is_available() else "cpu"
        ),
    }

    # Load the tokenizer
    tokenizer = T5Tokenizer.from_pretrained(
        model_source,
        do_lower_case=False,
        legacy=False,
        local_files_only=local_files_only,
    )

    # Load the model
    model_options = {"local_files_only": local_files_only}
    if use_cuda:
        # Materialize the 3B-parameter checkpoint directly in fp16. Loading
        # fp32 first and only then calling ``half()`` requires more memory than
        # an 8 GiB GPU before inference even begins.
        model_options["torch_dtype"] = torch.float16
    try:
        model = T5EncoderModel.from_pretrained(
            model_source,
            **model_options,
        )
    except (MemoryError, RuntimeError) as error:
        if not _is_oom_error(error):
            raise
        raise RuntimeError(
            "ProtTrans checkpoint materialization exhausted host memory"
        ) from error
    try:
        model = model.to(device)
    except (MemoryError, RuntimeError) as error:
        if not _is_oom_error(error) or not use_cuda:
            raise
        fallback_to_cpu = True
    else:
        fallback_to_cpu = False
    if fallback_to_cpu:
        del model
        gc.collect()
        torch.cuda.empty_cache()
        device = torch.device("cpu")
        model = T5EncoderModel.from_pretrained(
            model_source,
            local_files_only=local_files_only,
        )
        model = model.to(device)
        audit.update(
            {
                "fallback_reason": "cuda_model_transfer_oom",
                "final_device": "cpu",
                "final_dtype": "float32",
            }
        )

    result = (model, tokenizer, device)
    if return_audit:
        return (*result, audit)
    return result

def embedding(sequences, model, tokenizer, device):
    # prepare your protein sequences as a list
    sequence_examples = sequences

    # replace all rare/ambiguous amino acids by X and introduce white-space between all amino acids
    sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]

    # tokenize sequences and pad up to the longest sequence in the batch
    ids = tokenizer(sequence_examples, add_special_tokens=True, padding="longest")

    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)

    # generate embeddings
    with torch.no_grad():
        embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)
    
    embs=[]
    # Move to CPU once to avoid multiple device transfers if not needed immediately on GPU for other things
    # .float() might be needed if half precision was used
    last_hidden_states = embedding_repr.last_hidden_state.detach().cpu()
    
    for number, sequence in enumerate(sequences):
        # average over the sequence length (excluding padding if any, though attention mask handles computation, we need to correct length here)
        # Note: The original code used :len(sequence) which is correct for unpadded original length.
        # But wait, sequence_examples has spaces added. The tokenizer handles that.
        # The original code logic: embedding_repr.last_hidden_state[number,:len(sequence)]
        # However, len(sequence) is the length of the original amino acid string (e.g. "MVL...").
        # The tokenized input has spaces, so it's approx double length? No, tokenizer splits by space.
        # Let's double check T5Tokenizer. It treats "M V L" as 3 tokens?
        # ProtT5 tokenizer is character based mostly if space separated?
        # Actually better to stick to original logic: :len(sequence)
        # But wait, original code: 
        # sequence_examples = [" ".join(list...
        # ids = tokenizer(sequence_examples...
        # If I have "M V" (len 3 chars), tokenizer might produce tokens corresponding to M and V.
        # If strict preservation of original logic is desired:
        embs.append(last_hidden_states[number, :len(sequence)].mean(dim=0).numpy())
            
    return embs


def generate_uniform_embeddings(
    sequences,
    *,
    model_path=None,
    local_files_only=False,
    batch_size=8,
):
    model, tokenizer, device, audit = load_model(
        model_path=model_path,
        local_files_only=local_files_only,
        return_audit=True,
    )
    if device.type != "cuda":
        embeddings = [
            embedding([sequence], model, tokenizer, device)[0]
            for sequence in tqdm.tqdm(
                sequences,
                desc="CPU Processing",
            )
        ]
        return embeddings, model, tokenizer, device, audit

    embeddings = [None] * len(sequences)
    fallback_required = False
    for batch_start in tqdm.tqdm(
        range(0, len(sequences), batch_size),
        desc="GPU Processing",
    ):
        batch_end = min(batch_start + batch_size, len(sequences))
        batch_sequences = sequences[batch_start:batch_end]
        batch_indices = list(range(batch_start, batch_end))
        try:
            values = embedding(
                batch_sequences,
                model,
                tokenizer,
                device,
            )
            for index, value in zip(batch_indices, values):
                embeddings[index] = value
        except (MemoryError, RuntimeError) as error:
            if not _is_oom_error(error):
                raise
            batch_oom = True
        else:
            batch_oom = False
        if batch_oom:
            torch.cuda.empty_cache()
            for index in batch_indices:
                try:
                    embeddings[index] = embedding(
                        [sequences[index]],
                        model,
                        tokenizer,
                        device,
                    )[0]
                except (MemoryError, RuntimeError) as single_error:
                    if not _is_oom_error(single_error):
                        raise
                    single_oom = True
                else:
                    single_oom = False
                if single_oom:
                    fallback_required = True
                    break
        if fallback_required:
            break

    if not fallback_required:
        if any(value is None for value in embeddings):
            raise RuntimeError(
                "ProtTrans GPU embedding pass left incomplete outputs"
            )
        return embeddings, model, tokenizer, device, audit

    del model
    gc.collect()
    torch.cuda.empty_cache()
    model, tokenizer, device, cpu_audit = load_model(
        model_path=model_path,
        local_files_only=local_files_only,
        force_cpu=True,
        return_audit=True,
    )
    embeddings = [
        embedding([sequence], model, tokenizer, device)[0]
        for sequence in tqdm.tqdm(
            sequences,
            desc="CPU full recompute",
        )
    ]
    cpu_audit["fallback_reason"] = (
        "cuda_inference_oom_full_cpu_recompute"
    )
    cpu_audit["requested_device"] = "cuda"
    return embeddings, model, tokenizer, device, cpu_audit

def findtargethomos(data,dataindex,distance=1.2):
    return select_by_pca_distance(data, dataindex, distance)

def findtargethomos_by_level(data, dataindex, level_map, low_distance=1.2, high_distance=None):
    del level_map, high_distance
    return select_by_pca_distance(data, dataindex, low_distance)

def filter_by_prottrans(
    name,
    refname,
    distance,
    high_distance=None,
    *,
    model_path=None,
    local_files_only=False,
    write_audit=False,
):
    del high_distance
    batch_size=8
    prottrans_model_path = model_path
    
    homologs = pd.read_excel(f'working/{name}/matrix_homolog{name}_preforprottrans.xlsx')
    if 'confidence_level' not in homologs.columns:
        homologs['confidence_level'] = 'high'
    groupedhomologs = {}
    for refgene, group in homologs.groupby('refmodelgene'):
        groupedhomologs[refgene] = list(zip(group['tarmodelgene'], group['confidence_level']))
    model_paths = {
        "yeast": ("sbml", "models/yeast-GEM.xml"),
        "ecoli": ("json", "models/iML1515.json"),
        "strco": ("sbml", "models/Sco-GEM.xml"),
        "human": ("sbml", "models/Human-GEM.xml"),
        "synechocystis": ("sbml", "models/iSynCJ816.xml"),
    }
    try:
        model_format, reference_model_path = model_paths[refname]
    except KeyError as error:
        raise ValueError(f"Unsupported ProtTrans reference: {refname}") from error
    refmodel = (
        cobra.io.load_json_model(reference_model_path)
        if model_format == "json"
        else cobra.io.read_sbml_model(reference_model_path)
    )
    genes=[gene.id for gene in refmodel.genes]
    final_homologs=[]
    sortedgenes={}
    for i in groupedhomologs.keys():
        if i in genes:
            sortedgenes[i]=groupedhomologs[i]
        else:
            for gene, _level in groupedhomologs[i]:
                final_homologs.append([i,gene])
    indexes=[]
    genes2embedding=[]
    for i in sortedgenes.keys():
        level_map = {}
        unique_genes = []
        for gene, level in sortedgenes[i]:
            if gene not in level_map:
                unique_genes.append(gene)
                level_map[gene] = level
            elif level == 'high':
                level_map[gene] = 'high'
        if len(unique_genes)>1:
            genes2embedding.append((i, unique_genes, level_map))
        else:
            final_homologs.append([i, unique_genes[0]])
    indexseq=[]
    records1={
        protein_key(record.id): str(record.seq)
        for record in SeqIO.parse(
            f'./data_available/{refname}.fasta', 'fasta'
        )
    }
    records2={
        protein_key(record.id): str(record.seq)
        for record in SeqIO.parse(
            f'./working/{name}/{name}.fasta', 'fasta'
        )
    }
    findrefnames.predata(refname)
    for genes in genes2embedding:
        indexes.append(findrefnames.findmodelname(genes[0]))
        indexes=indexes+genes[1]
    indexes=list(dict.fromkeys(indexes))
    missing_sequences = []
    for i in indexes:
        if i in records1:
            indexseq.append(records1[i])
        elif i in records2:
            indexseq.append(records2[i])
        else:
            missing_sequences.append(i)
    if missing_sequences:
        raise KeyError(
            "ProtTrans sequences missing for IDs: "
            + ", ".join(missing_sequences[:20])
        )
    audit_path = f"working/{name}/prottrans_execution.json"
    try:
        (
            indexemb,
            model,
            tokenizer,
            device,
            execution_audit,
        ) = generate_uniform_embeddings(
            indexseq,
            model_path=prottrans_model_path,
            local_files_only=local_files_only,
            batch_size=batch_size,
        )
    except Exception as error:
        if write_audit:
            _write_json_atomic(
                audit_path,
                {
                    "error_message": str(error)[:1000],
                    "error_type": type(error).__name__,
                    "model_source": str(
                        prottrans_model_path
                        or "Rostlab/prot_t5_xl_uniref50"
                    ),
                    "status": "failed",
                    "uniform_embedding_precision": True,
                },
            )
        raise
    
    embeddings={}
    for i in range(len(indexes)):
        embeddings[indexes[i]]=indexemb[i]
    del indexemb,indexseq,indexes
    embedding_path = 'working/{name}/{name}.emb'.format(name=name)
    with open(embedding_path, 'wb') as handle:
        # noinspection PyTypeChecker
        pickle.dump(embeddings, handle)
    # with open('working/{name}/{name}.emb'.format(name=name), 'rb') as handle:
    #     embeddings=pickle.load(handle)
    for genes in genes2embedding:
        reference_id = findrefnames.findmodelname(genes[0])
        if reference_id not in embeddings:
            raise KeyError(f"ProtTrans reference embedding missing: {reference_id}")
        data=[embeddings[reference_id]]
        dataindex=[genes[0]]
        for gene in genes[1]:
            if gene not in embeddings:
                raise KeyError(f"ProtTrans target embedding missing: {gene}")
            data.append(embeddings[gene])
            dataindex.append(gene)
        dataindex=np.array(dataindex)
        homoss=findtargethomos(
            data,
            dataindex,
            distance=distance,
        )[1:]
        for homomomo in homoss:
            final_homologs.append([genes[0],homomomo])
    del records1, records2
    matrix_path = f'working/{name}/matrix_homolog{name}.xlsx'
    pd.DataFrame(final_homologs).to_excel(matrix_path)
    if write_audit:
        execution_audit.update(
            {
                "embedding_count": len(embeddings),
                "embedding_output": _file_record(embedding_path),
                "group_count": len(genes2embedding),
                "matrix_output": _file_record(matrix_path),
                "status": "completed",
                "torch_version": importlib.metadata.version("torch"),
                "transformers_version": importlib.metadata.version(
                    "transformers"
                ),
                "uniform_embedding_precision": True,
            }
        )
        _write_json_atomic(audit_path, execution_audit)
    
    # Clean up memory
    if 'model' in locals():
        del model
    if 'tokenizer' in locals():
        del tokenizer
    gc.collect()
    torch.cuda.empty_cache()
