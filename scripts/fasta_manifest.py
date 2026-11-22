#!/usr/bin/env python3
"""Create auditable FASTA QC manifests without third-party dependencies."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


DNA_ALPHABET = frozenset("ACGTURYSWKMBDHVN")
PROTEIN_ALPHABET = frozenset("ABCDEFGHIKLMNPQRSTVWXYZJUO*")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_fasta(path: Path, sequence_type: str) -> dict[str, Any]:
    allowed = DNA_ALPHABET if sequence_type == "dna" else PROTEIN_ALPHABET
    ids: list[str] = []
    lengths: list[int] = []
    invalid_records: list[dict[str, Any]] = []
    internal_stop_ids: list[str] = []
    current_id: str | None = None
    current_sequence: list[str] = []

    def finish_record() -> None:
        if current_id is None:
            return
        sequence = "".join(current_sequence).upper()
        ids.append(current_id)
        lengths.append(len(sequence))
        invalid = sorted(set(sequence) - allowed)
        if invalid:
            invalid_records.append({"id": current_id, "characters": "".join(invalid)})
        if sequence_type == "protein" and "*" in sequence[:-1]:
            internal_stop_ids.append(current_id)

    with path.open("rt", encoding="utf-8", errors="strict") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                current_id = header.split()[0]
                current_sequence = []
            elif current_id is None:
                if line.strip():
                    raise ValueError(f"Sequence before first header at line {line_number}")
            else:
                current_sequence.append(line)
    finish_record()

    duplicate_ids = sorted(identifier for identifier, n in Counter(ids).items() if n > 1)
    empty_ids = [identifier for identifier, length in zip(ids, lengths) if length == 0]
    empty_file_issue = int(not ids or sum(lengths) == 0)
    blocking_issue_count = (
        empty_file_issue
        + len(duplicate_ids)
        + len(empty_ids)
        + len(invalid_records)
        + len(internal_stop_ids)
    )
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "sequence_type": sequence_type,
        "records": len(ids),
        "unique_ids": len(set(ids)),
        "duplicate_ids": duplicate_ids,
        "empty_sequence_ids": empty_ids,
        "total_sequence_length": sum(lengths),
        "minimum_length": min(lengths) if lengths else 0,
        "maximum_length": max(lengths) if lengths else 0,
        "invalid_records": invalid_records,
        "internal_stop_ids": internal_stop_ids,
        "blocking_issue_count": blocking_issue_count,
        "empty_fasta": bool(empty_file_issue),
        "valid": blocking_issue_count == 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--type", choices=("dna", "protein"), required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--no-fail",
        action="store_true",
        help="Write the manifest but return zero even when blocking issues exist.",
    )
    args = parser.parse_args()

    result = inspect_fasta(args.fasta, args.type)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0 if result["valid"] or args.no_fail else 2


if __name__ == "__main__":
    raise SystemExit(main())
