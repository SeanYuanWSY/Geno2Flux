#!/usr/bin/env python3
"""Create a reversible, sequence-preserving FASTA ID normalization copy."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path


SAFE_ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


@dataclass(frozen=True)
class Record:
    original_id: str
    original_header: str
    sequence: str


def read_fasta(path: Path) -> list[Record]:
    records: list[Record] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    records.append(make_record(header, chunks, path))
                header = line[1:]
                if not header or not header.split()[0]:
                    raise ValueError(f"empty FASTA header at {path}:{line_number}")
                chunks = []
            elif line.strip():
                if header is None:
                    raise ValueError(f"sequence before first header at {path}:{line_number}")
                chunks.append(line.strip())
    if header is not None:
        records.append(make_record(header, chunks, path))
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    identifiers = [record.original_id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"duplicate FASTA IDs in {path}")
    return records


def make_record(header: str, chunks: list[str], path: Path) -> Record:
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"empty FASTA sequence for {header!r} in {path}")
    return Record(header.split()[0], header, sequence)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--map", dest="mapping", type=Path, required=True)
    parser.add_argument("--prefix", default="PROT")
    args = parser.parse_args()

    resolved_paths = [args.input.resolve(strict=False), args.output.resolve(strict=False), args.mapping.resolve(strict=False)]
    if len(set(resolved_paths)) != len(resolved_paths):
        raise ValueError("input, output, and mapping paths must be pairwise distinct")
    if not SAFE_ID_RE.fullmatch(args.prefix):
        raise ValueError("prefix must use only safe GFF/FASTA identifier characters")
    if args.output.exists() or args.mapping.exists():
        raise FileExistsError("output and mapping paths must not already exist")

    records = read_fasta(args.input)
    width = max(6, len(str(len(records))))
    safe_ids = [f"{args.prefix}{index:0{width}d}" for index in range(1, len(records) + 1)]
    if len(safe_ids) != len(set(safe_ids)) or not all(SAFE_ID_RE.fullmatch(value) for value in safe_ids):
        raise ValueError("generated safe IDs are invalid or non-unique")
    original_ids = {record.original_id for record in records}
    if original_ids & set(safe_ids):
        raise ValueError("safe and original FASTA ID namespaces overlap")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.mapping.parent.mkdir(parents=True, exist_ok=True)
    output_temp: Path | None = None
    mapping_temp: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="ascii", newline="\n", dir=args.output.parent,
            prefix=f".{args.output.name}.", delete=False,
        ) as output_handle, tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=args.mapping.parent,
            prefix=f".{args.mapping.name}.", delete=False,
        ) as map_handle:
            output_temp = Path(output_handle.name)
            mapping_temp = Path(map_handle.name)
            fields = ["order", "original_id", "safe_id", "length", "sequence_sha256", "original_header"]
            writer = csv.DictWriter(map_handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for order, (record, safe_id) in enumerate(zip(records, safe_ids), 1):
                output_handle.write(f">{safe_id}\n")
                for start in range(0, len(record.sequence), 60):
                    output_handle.write(record.sequence[start : start + 60] + "\n")
                writer.writerow({
                    "order": order,
                    "original_id": record.original_id,
                    "safe_id": safe_id,
                    "length": len(record.sequence),
                    "sequence_sha256": sha256_text(record.sequence),
                    "original_header": record.original_header,
                })

        normalized = read_fasta(output_temp)
        if len(normalized) != len(records) or len(safe_ids) != len(records):
            raise AssertionError("record count changed during normalization")
        for original, safe, expected_safe_id in zip(records, normalized, safe_ids):
            if safe.original_id != expected_safe_id or safe.sequence != original.sequence:
                raise AssertionError(f"sequence/order mismatch after normalization for {original.original_id}")
        os.replace(mapping_temp, args.mapping)
        mapping_temp = None
        os.replace(output_temp, args.output)
        output_temp = None
    finally:
        for path in (output_temp, mapping_temp):
            if path is not None and path.exists():
                path.unlink()


if __name__ == "__main__":
    main()
