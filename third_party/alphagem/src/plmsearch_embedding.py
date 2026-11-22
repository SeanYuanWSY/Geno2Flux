"""Checkpointed ESM-1b embedding generation with provenance sidecars."""

import argparse
import json
import pickle
from pathlib import Path

import torch
import torch.serialization as serialization
from esm import pretrained
from esm.data import FastaBatchedDataset

from src.pipeline_contracts import (
    assert_partial_embedding_keys,
    file_sha256,
    write_embedding_provenance,
)


serialization.add_safe_globals([argparse.Namespace])


def _atomic_pickle(payload, destination):
    destination = Path(destination)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(destination)


def _atomic_json(payload, destination):
    destination = Path(destination)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(destination)


def _resolve_model_file(model_identifier):
    supplied = Path(model_identifier)
    candidates = [supplied, Path(f"{model_identifier}.pt")]
    if not supplied.is_absolute():
        candidates.append(
            Path(torch.hub.get_dir())
            / "checkpoints"
            / f"{supplied.stem}.pt"
        )
    return next((path.resolve() for path in candidates if path.is_file()), None)


def _generate(
    fasta,
    embedding_result,
    *,
    esm_model_path,
    nogpu,
    token_budget,
    checkpoint_every,
    embedding_semantics,
):
    fasta = Path(fasta)
    embedding_result = Path(embedding_result)
    partial = embedding_result.with_suffix(".partial.pkl")
    partial_metadata = partial.with_suffix(".meta.json")
    metadata = embedding_result.with_name(
        embedding_result.stem + ".meta.json"
    )

    pooling_names = {
        "corrected": "mean_residue_excluding_bos_eos",
        "legacy_repo": (
            "legacy_truncated_slice_includes_eos_for_long_sequences"
        ),
    }
    try:
        pooling = pooling_names[embedding_semantics]
    except KeyError as error:
        raise ValueError(
            f"Unknown embedding semantics: {embedding_semantics}"
        ) from error

    model, alphabet = pretrained.load_model_and_alphabet(esm_model_path)
    model.eval()
    model_file = _resolve_model_file(esm_model_path)
    cache_context = {
        "schema_version": 1,
        "fasta_sha256": file_sha256(fasta),
        "model_identifier": Path(esm_model_path).stem,
        "model_sha256": (
            file_sha256(model_file) if model_file is not None else None
        ),
        "representation_layer": 33,
        "pooling": pooling,
        "max_residues": 1022,
    }
    use_gpu = torch.cuda.is_available() and not nogpu
    if use_gpu:
        model = model.cuda()
        print("Transferred ESM-1b to GPU")

    dataset = FastaBatchedDataset.from_file(str(fasta))
    batches = dataset.get_batch_indices(
        token_budget,
        extra_toks_per_seq=1,
    )
    loader = torch.utils.data.DataLoader(
        dataset,
        collate_fn=alphabet.get_batch_converter(1022),
        batch_sampler=batches,
    )
    embeddings = {}
    if partial.is_file():
        if not partial_metadata.is_file():
            raise ValueError(
                f"Partial embedding provenance is missing: {partial_metadata}"
            )
        observed_context = json.loads(
            partial_metadata.read_text(encoding="utf-8")
        )
        expected_context = {
            **cache_context,
            "embedding_sha256": file_sha256(partial),
        }
        if observed_context != expected_context:
            raise ValueError(
                "Partial embedding cache does not match the current "
                "FASTA/model/pooling contract"
            )
        with partial.open("rb") as handle:
            embeddings = pickle.load(handle)
        assert_partial_embedding_keys(embeddings, fasta)
        print(f"Resuming {len(embeddings)} cached proteins from {partial}")

    print(f"Read {fasta} with {len(dataset)} sequences")
    with torch.no_grad():
        for batch_index, (labels, sequences, tokens) in enumerate(loader):
            if all(label in embeddings for label in labels):
                continue
            print(
                f"Processing {batch_index + 1}/{len(batches)} "
                f"({tokens.size(0)} sequences)"
            )
            if use_gpu:
                tokens = tokens.to(device="cuda", non_blocking=True)
            representations = model(
                tokens,
                repr_layers=[33],
                return_contacts=False,
            )["representations"][33]
            for index, label in enumerate(labels):
                if embedding_semantics == "corrected":
                    residue_stop = min(len(sequences[index]), 1022) + 1
                else:
                    residue_stop = len(sequences[index]) + 1
                embeddings[label] = (
                    representations[index, 1:residue_stop]
                    .mean(0)
                    .clone()
                    .cpu()
                )
            if (batch_index + 1) % checkpoint_every == 0:
                _atomic_pickle(embeddings, partial)
                _atomic_json(
                    {
                        **cache_context,
                        "embedding_sha256": file_sha256(partial),
                    },
                    partial_metadata,
                )
                if use_gpu:
                    torch.cuda.empty_cache()

    if len(embeddings) != len(dataset):
        raise RuntimeError(
            f"Embedding generation incomplete: {len(embeddings)}/{len(dataset)}"
        )
    _atomic_pickle(embeddings, embedding_result)
    partial.unlink(missing_ok=True)
    partial_metadata.unlink(missing_ok=True)
    write_embedding_provenance(
        embedding_result,
        fasta,
        metadata,
        model_identifier=esm_model_path,
        model_path=model_file,
        token_budget=token_budget,
        generator_path=Path(__file__),
        pooling=pooling,
    )
    return embedding_result


def embedding_generate(
    name,
    esm_model_path="esm1b_t33_650M_UR50S",
    nogpu=False,
    token_budget=16384,
    checkpoint_every=25,
    embedding_semantics="corrected",
):
    return _generate(
        f"./working/{name}/{name}.fasta",
        f"./working/{name}/{name}_embedding.pkl",
        esm_model_path=esm_model_path,
        nogpu=nogpu,
        token_budget=token_budget,
        checkpoint_every=checkpoint_every,
        embedding_semantics=embedding_semantics,
    )


def embedding_generate_base(
    name,
    esm_model_path="esm1b_t33_650M_UR50S",
    nogpu=False,
    token_budget=16384,
    checkpoint_every=25,
    embedding_semantics="corrected",
):
    return _generate(
        f"./data_available/{name}.fasta",
        f"./working/{name}_embedding.pkl",
        esm_model_path=esm_model_path,
        nogpu=nogpu,
        token_budget=token_budget,
        checkpoint_every=checkpoint_every,
        embedding_semantics=embedding_semantics,
    )
