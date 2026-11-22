#!/usr/bin/env python3
"""Verify that a soft-masked FASTA differs from its source only by case."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


class SoftmaskValidationError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(chunks)))
                fields = line[1:].split()
                identifier = fields[0] if fields else ""
                if not identifier or any(identifier == prior for prior, _ in records):
                    raise SoftmaskValidationError(
                        f"invalid or duplicate FASTA ID at line {line_number}: {path}"
                    )
                chunks = []
            elif line:
                if identifier is None or any(character.isspace() for character in line):
                    raise SoftmaskValidationError(f"invalid FASTA structure at line {line_number}")
                chunks.append(line)
    if identifier is not None:
        records.append((identifier, "".join(chunks)))
    if not records or any(not sequence for _, sequence in records):
        raise SoftmaskValidationError(f"empty FASTA record set: {path}")
    return records


def validate(
    original_path: Path,
    masked_path: Path,
    per_scaffold_path: Path,
    intervals_path: Path,
    summary_path: Path,
) -> dict[str, object]:
    original = _read_fasta(original_path)
    masked = _read_fasta(masked_path)
    if len(original) != len(masked):
        raise SoftmaskValidationError("record count differs between original and masked FASTA")
    original_ids = [identifier for identifier, _ in original]
    masked_ids = [identifier for identifier, _ in masked]
    if original_ids != masked_ids:
        raise SoftmaskValidationError(
            "record IDs or order differ between original and masked FASTA"
        )

    rows: list[tuple[str, int, int, float]] = []
    intervals: list[tuple[str, int, int]] = []
    total_bases = 0
    total_masked = 0
    for index in range(len(original)):
        identifier, source = original[index]
        _, result = masked[index]
        if len(source) != len(result):
            raise SoftmaskValidationError(f"sequence length differs for {identifier}")
        if source.upper() != result.upper():
            raise SoftmaskValidationError(
                f"sequence content differs beyond letter case for {identifier}"
            )
        masked_count = sum(character.islower() for character in result)
        rows.append((identifier, len(result), masked_count, masked_count / len(result)))
        total_bases += len(result)
        total_masked += masked_count
        start: int | None = None
        for index, character in enumerate(result):
            if character.islower() and start is None:
                start = index
            elif not character.islower() and start is not None:
                intervals.append((identifier, start, index))
                start = None
        if start is not None:
            intervals.append((identifier, start, len(result)))

    per_scaffold_path.parent.mkdir(parents=True, exist_ok=True)
    with per_scaffold_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("scaffold\tlength\tmasked_bases\tmasked_fraction\n")
        for identifier, length, masked_count, fraction in rows:
            handle.write(f"{identifier}\t{length}\t{masked_count}\t{fraction:.12f}\n")
    with intervals_path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, start, end in intervals:
            handle.write(f"{identifier}\t{start}\t{end}\n")

    summary: dict[str, object] = {
        "status": "PASS_CASE_ONLY_DIFFERENCE",
        "original_path": str(original_path),
        "original_sha256": _sha256(original_path),
        "masked_path": str(masked_path),
        "masked_sha256": _sha256(masked_path),
        "records": len(original),
        "total_bases": total_bases,
        "masked_bases": total_masked,
        "masked_fraction": total_masked / total_bases,
        "interval_count": len(intervals),
        "ids_equal": True,
        "order_equal": True,
        "lengths_equal": True,
        "uppercase_sequences_equal": True,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", required=True, type=Path)
    parser.add_argument("--masked", required=True, type=Path)
    parser.add_argument("--per-scaffold", required=True, type=Path)
    parser.add_argument("--intervals", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    result = validate(
        args.original,
        args.masked,
        args.per_scaffold,
        args.intervals,
        args.summary,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
