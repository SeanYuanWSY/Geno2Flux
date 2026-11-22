import filecmp
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO


STRUCTURE_COLUMNS = ("Structure", "Structure Data")


def _file_record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def _write_json_atomic(path, payload):
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _alphagem_record_fields(record):
    fields = record.id.split("|")
    if len(fields) < 3:
        raise ValueError(
            "AlphaGEM FASTA identifiers must use the form "
            "'prefix|Entry|Entry Name'; got {!r}".format(record.id)
        )
    return fields[1], fields[2]


def _structure_column(table):
    matches = [column for column in STRUCTURE_COLUMNS if column in table.columns]
    if not matches:
        raise ValueError(
            "Structure mapping must contain 'Structure' (README format) or "
            "'Structure Data' (legacy format)"
        )
    return matches[0]


def _safe_entry(entry):
    entry = str(entry).strip()
    if not entry or entry.lower() == "nan":
        raise ValueError("Structure mapping contains an empty Entry value")
    if Path(entry).name != entry or entry in {".", ".."}:
        raise ValueError("Unsafe Entry value in structure mapping: {!r}".format(entry))
    return entry


def _find_structure_file(structure_dir, structure_value):
    value = str(structure_value).strip()
    if not value or value.lower() == "nan":
        raise ValueError("Structure mapping contains an empty structure value")

    basename = Path(value).name
    candidates = [basename]
    if not basename.lower().endswith(".pdb"):
        candidates.extend(
            [
                "{}.pdb".format(basename),
                "AF-{}-F1-model_v4.pdb".format(basename),
            ]
        )

    for candidate_name in candidates:
        candidate = structure_dir / candidate_name
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Mapped structure {!r} was not found in {}".format(
            value,
            structure_dir,
        )
    )


def _validate_staged_copy(source, destination):
    if not destination.is_file():
        raise FileNotFoundError("Staged structure was not created: {}".format(destination))
    if not filecmp.cmp(str(source), str(destination), shallow=False):
        raise IOError(
            "Staged structure differs from its input: {} -> {}".format(
                source,
                destination,
            )
        )


def generare_list_without_structure(name):
    rows = []
    genes = "working/{}/{}.fasta".format(name, name)
    records = SeqIO.parse(genes, 'fasta')
    for record in records:
        entry, entry_name = _alphagem_record_fields(record)
        rows.append(
            {
                "Entry": entry,
                "Entry Name": entry_name,
                "Structure Data": entry,
            }
        )
    table = pd.DataFrame(rows)
    table.to_excel("working/{}/{}.xlsx".format(name, name), index=False)


def generare_list_with_structure(
    name,
    listfile='',
    structurefile='',
    *,
    allow_missing_structures=False,
):
    """Prepare a validated, immutable-input structure staging directory.

    The returned absolute path is the only target-structure path that structure
    alignment and pLDDT filtering should consume.
    """
    if not listfile:
        rows = []
        genes = "working/{}/{}.fasta".format(name, name)
        records = SeqIO.parse(genes, 'fasta')
        for record in records:
            entry, entry_name = _alphagem_record_fields(record)
            rows.append(
                {
                    "Entry": entry,
                    "Entry Name": entry_name,
                    "Structure Data": entry,
                }
            )
        table = pd.DataFrame(rows)
    else:
        list_path = Path(listfile).expanduser()
        if not list_path.is_file():
            raise FileNotFoundError(
                "Structure mapping file not found: {}".format(list_path)
            )
        table = pd.read_excel(str(list_path))

    missing_columns = [
        column for column in ("Entry", "Entry Name") if column not in table.columns
    ]
    if missing_columns:
        raise ValueError(
            "Structure mapping is missing required columns: {}".format(
                ", ".join(missing_columns)
            )
        )
    structure_column = _structure_column(table)

    structure_dir = Path(structurefile).expanduser().resolve()
    if not structure_dir.is_dir():
        raise NotADirectoryError(
            "Structure input directory not found: {}".format(structure_dir)
        )

    stage_dir = (
        Path.cwd() / "working" / name / "structures"
    ).resolve()
    stage_dir.mkdir(parents=True, exist_ok=True)

    expected_names = set()
    seen_entries = set()
    staged_pairs = []
    staged_entries = []
    missing_structures = []
    mapping_row_count = 0
    for row_index, row in table.iterrows():
        mapping_row_count += 1
        entry = _safe_entry(row["Entry"])
        destination_name = "AF-{}-F1-model_v4.pdb".format(entry)
        if entry in seen_entries:
            raise ValueError(
                "Duplicate Entry {!r} in structure mapping (row {})".format(
                    entry,
                    row_index,
                )
            )
        seen_entries.add(entry)
        structure_value = str(row[structure_column]).strip()
        try:
            source = _find_structure_file(structure_dir, structure_value)
        except FileNotFoundError:
            if not allow_missing_structures:
                raise
            missing_structures.append(
                {
                    "entry": entry,
                    "row_index": int(row_index),
                    "structure_value": structure_value,
                }
            )
            continue
        destination = stage_dir / destination_name
        expected_names.add(destination_name)
        staged_pairs.append((source, destination))
        staged_entries.append((entry, source, destination))

    unexpected = {
        path.name for path in stage_dir.glob("*.pdb")
    } - expected_names
    if unexpected:
        raise ValueError(
            "Structure staging directory contains files not present in the "
            "current mapping: {}".format(", ".join(sorted(unexpected)))
        )

    for source, destination in staged_pairs:
        if source.resolve() != destination.resolve():
            shutil.copy2(str(source), str(destination))
        _validate_staged_copy(source, destination)

    output = Path.cwd() / "working" / name / "{}.xlsx".format(name)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_excel(str(output), index=False)

    staged_names = {path.name for path in stage_dir.glob("*.pdb")}
    if staged_names != expected_names:
        raise RuntimeError(
            "Structure staging validation failed: expected {}, found {}".format(
                len(expected_names),
                len(staged_names),
            )
        )
    if allow_missing_structures:
        audit_path = (
            Path.cwd()
            / "working"
            / name
            / "structure_staging_audit.json"
        )
        _write_json_atomic(
            audit_path,
            {
                "allow_missing_structures": True,
                "mapping_row_count": mapping_row_count,
                "missing_structure_count": len(missing_structures),
                "missing_structures": missing_structures,
                "staged_structure_count": len(staged_entries),
                "staged_structures": [
                    {
                        "entry": entry,
                        "file": destination.name,
                        **_file_record(destination),
                    }
                    for entry, _source, destination in staged_entries
                ],
            },
        )
    return str(stage_dir)

