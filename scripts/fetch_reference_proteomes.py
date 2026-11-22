#!/usr/bin/env python3
"""Download NCBI RefSeq or UniProt reference proteomes for homology-guided gene prediction.

Creates ``combined.faa``, a reversible ``protein_id_map.tsv``, and ``manifest.json``
with SHA-256 checksums so evidence-guided predictors (such as funannotate) can use
auditable protein homology inputs.

Examples:
    # NCBI RefSeq assembly accession:
    uv run --frozen python scripts/fetch_reference_proteomes.py \
        --accession SACCE=GCF_000146045.2 \
        --output /work/ref_proteomes

    # Custom UniProt proteome IDs (e.g. UP000002311):
    uv run --frozen python scripts/fetch_reference_proteomes.py \
        --uniprot-proteome YEAST=UP000002311 \
        --output /work/ref_proteomes
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen
import zipfile

from Bio import SeqIO
from Bio.Seq import Seq


ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def _parse_kv(items: list[str] | None) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in items or []:
        if "=" in item:
            token, acc = item.split("=", 1)
        else:
            token = acc = item
        token = token.strip()
        acc = acc.strip()
        if not token or not acc:
            raise ValueError(f"Invalid source specification: {item!r}")
        if token in out:
            raise ValueError(f"Duplicate source token: {token!r}")
        out[token] = acc
    return out


def _ingest_fasta_bytes(
    *,
    fasta_bytes: bytes,
    token: str,
    source_id: str,
    source_file_name: str,
    records: list,
    mapping: list[dict[str, str]],
) -> tuple[int, str]:
    source_sha = hashlib.sha256(fasta_bytes).hexdigest()
    seen: set[str] = set()
    count = 0
    for r in SeqIO.parse(io.StringIO(fasta_bytes.decode("utf-8")), "fasta"):
        sequence = str(r.seq).upper()
        if sequence.endswith("*"):
            sequence = sequence[:-1]
        if not sequence or (set(sequence) - ALLOWED_AA) or r.id in seen:
            raise ValueError(f"Invalid reference protein: {source_id}/{r.id}")
        seen.add(r.id)
        original = r.id
        prefixed = f"{token}|{original}"
        mapping.append(
            {
                "prefixed_id": prefixed,
                "source_token": token,
                "source_accession": source_id,
                "original_id": original,
                "original_description": r.description,
                "source_file": source_file_name,
                "source_sha256": source_sha,
            }
        )
        r.id = prefixed
        r.description = ""
        r.seq = Seq(sequence)
        records.append(r)
        count += 1
    if count == 0:
        raise ValueError(f"Empty reference proteome for {source_id}")
    return count, source_sha


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", "--output-dir", dest="output", type=Path, required=True)
    p.add_argument(
        "--accession",
        action="append",
        help="NCBI RefSeq/GenBank assembly accession (TOKEN=GCF_xxx or GCF_xxx); may be repeated",
    )
    p.add_argument(
        "--uniprot-proteome",
        action="append",
        help="UniProt proteome ID (TOKEN=UP00000xxxx or UP00000xxxx); may be repeated",
    )
    a = p.parse_args(argv)
    refseq_sources = _parse_kv(a.accession)
    uniprot_sources = _parse_kv(a.uniprot_proteome)
    if not refseq_sources and not uniprot_sources:
        p.error("specify at least one --accession or --uniprot-proteome")

    a.output.mkdir(parents=True, exist_ok=False)
    records: list = []
    mapping: list[dict[str, str]] = []
    manifest: list[dict[str, object]] = []

    for token, accession in refseq_sources.items():
        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            + accession
            + "/download?include_annotation_type=PROT_FASTA"
        )
        with urlopen(Request(url, headers={"User-Agent": "Geno2Flux/0.1"}), timeout=180) as response:
            body = response.read()
        archive = a.output / f"{accession}.zip"
        archive.write_bytes(body)
        with zipfile.ZipFile(io.BytesIO(body)) as z:
            names = [n for n in z.namelist() if n.endswith("/protein.faa") and accession in n]
            if len(names) != 1:
                raise ValueError(f"NCBI archive lacks an unambiguous protein FASTA: {accession}")
            fasta = z.read(names[0])
        source = a.output / f"{accession}.faa"
        source.write_bytes(fasta)
        count, source_sha = _ingest_fasta_bytes(
            fasta_bytes=fasta,
            token=token,
            source_id=accession,
            source_file_name=source.name,
            records=records,
            mapping=mapping,
        )
        manifest.append(
            {
                "provider": "ncbi_datasets",
                "accession": accession,
                "url": url,
                "records": count,
                "sha256": source_sha,
                "archive_sha256": hashlib.sha256(body).hexdigest(),
            }
        )

    for token, up_id in uniprot_sources.items():
        url = (
            "https://rest.uniprot.org/uniprotkb/stream"
            f"?compressed=true&format=fasta&query=(proteome:{up_id})"
        )
        with urlopen(Request(url, headers={"User-Agent": "Geno2Flux/0.1"}), timeout=180) as response:
            gz_body = response.read()
        archive = a.output / f"{up_id}.fasta.gz"
        archive.write_bytes(gz_body)
        fasta = gzip.decompress(gz_body)
        source = a.output / f"{up_id}.faa"
        source.write_bytes(fasta)
        count, source_sha = _ingest_fasta_bytes(
            fasta_bytes=fasta,
            token=token,
            source_id=up_id,
            source_file_name=source.name,
            records=records,
            mapping=mapping,
        )
        manifest.append(
            {
                "provider": "uniprot",
                "proteome_id": up_id,
                "url": url,
                "records": count,
                "sha256": source_sha,
                "archive_sha256": hashlib.sha256(gz_body).hexdigest(),
            }
        )

    if len({r.id for r in records}) != len(records):
        raise ValueError("Duplicate merged protein identifiers")
    SeqIO.write(records, a.output / "combined.faa", "fasta")
    with (a.output / "protein_id_map.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, list(mapping[0]), delimiter="\t")
        w.writeheader()
        w.writerows(mapping)
    combined_hash = hashlib.sha256((a.output / "combined.faa").read_bytes()).hexdigest()
    (a.output / "manifest.json").write_text(
        json.dumps({"sources": manifest, "combined_sha256": combined_hash}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(combined_hash)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
