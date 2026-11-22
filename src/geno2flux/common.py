"""Shared validation and local provenance."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def save_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False, default=str)
        + "\n"
    )

def table(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def finite(value: object) -> float | None:
    try:
        n = float(value)
    except (ValueError, TypeError):
        return None
    return n if math.isfinite(n) else None

def resource_check(root: Path, manifest: Path) -> dict:
    records = json.loads(manifest.read_text())
    if not isinstance(records, dict) or not records:
        raise ValueError("Resource lock must be a nonempty path-to-SHA256 object")
    resolved = root.resolve()
    for name, expected in records.items():
        p = (root / name).resolve()
        if not p.is_relative_to(resolved):
            raise ValueError("Resource lock path escapes the resource root")
        if not p.is_file() or sha256(p) != expected:
            raise ValueError("Resource fingerprint mismatch: " + name)
    return records
