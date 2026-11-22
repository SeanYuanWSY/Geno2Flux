from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from src.pipeline_contracts import normalize_alphagem_target_id


def _normalized_gene_id(record_id, *, strict=False):
    if strict:
        return normalize_alphagem_target_id(record_id)

    value = str(record_id)
    try:
        value = value.split("SerialNumber")[1]
    except IndexError:
        pass
    try:
        gene_id = value.split("|")[1]
    except IndexError:
        gene_id = value
    return f"sp|{gene_id}|{gene_id}"


def rename_genes_in_fasta(input_fasta, output_fasta, *, strict=False):
    records = []
    normalized_to_source = {}
    with Path(input_fasta).open(encoding="utf-8") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            source_id = record.id
            normalized_id = _normalized_gene_id(source_id, strict=strict)
            if strict and normalized_id in normalized_to_source:
                raise ValueError(
                    "Formal FASTA normalization would collapse distinct records "
                    f"{normalized_to_source[normalized_id]!r} and {source_id!r} "
                    f"to {normalized_id!r}"
                )
            normalized_to_source[normalized_id] = source_id

            sequence = str(record.seq)
            if strict and "*" in sequence[:-1]:
                raise ValueError(
                    "Formal FASTA normalization found an internal stop codon in "
                    f"{source_id!r}"
                )
            if strict:
                sequence = sequence.removesuffix("*")
                if not sequence:
                    raise ValueError(
                        "Formal FASTA normalization produced an empty sequence "
                        f"for {source_id!r}"
                    )
            else:
                sequence = sequence.replace("*", "")

            record.id = normalized_id
            record.name = ""
            record.description = ""
            record.seq = Seq(sequence)
            records.append(record)

    if strict and not records:
        raise ValueError("Formal FASTA normalization found no input records")

    output = Path(output_fasta)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as output_handle:
        written = SeqIO.write(records, output_handle, "fasta")
    if strict and written != len(records):
        raise RuntimeError(
            "Formal FASTA normalization did not write every validated record: "
            f"expected={len(records)}, written={written}"
        )


def handle(fasta, name, *, strict=False):
    rename_genes_in_fasta(
        fasta,
        f"./working/{name}/{name}.fasta",
        strict=strict,
    )
