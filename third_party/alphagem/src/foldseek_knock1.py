"""Run Foldseek with explicit score semantics and fresh formal outputs."""

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def _file_record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _copy_atomic(source, destination):
    source = Path(source)
    destination = Path(destination)
    temporary = None
    try:
        with source.open("rb") as input_handle, tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as output_handle:
            temporary = Path(output_handle.name)
            shutil.copyfileobj(input_handle, output_handle)
            output_handle.flush()
            os.fsync(output_handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _write_excel_atomic(table, destination):
    destination = Path(destination)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".xlsx",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
        table.to_excel(temporary, index=False)
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _write_json_atomic(payload, destination):
    destination = Path(destination)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _reject_unsafe_existing_output(path, label):
    path = Path(path)
    if path.is_symlink() or (path.exists() and not path.is_file()):
        raise RuntimeError(
            f"Formal Foldseek {label} must be a regular file: {path}"
        )


def _assert_formal_output_path(path, root, label):
    root = Path(root).resolve(strict=True)
    path = Path(path)
    try:
        relative = path.relative_to(root)
    except ValueError as error:
        raise RuntimeError(
            f"Formal Foldseek {label} escapes the repository: {path}"
        ) from error
    current = root
    for component in relative.parts:
        current = current / component
        if current.is_symlink():
            raise RuntimeError(
                f"Formal Foldseek {label} path does not permit symlinks: "
                f"{current}"
            )
    resolved = path.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise RuntimeError(
            f"Formal Foldseek {label} resolves outside the repository: "
            f"{path} -> {resolved}"
        ) from error


def foldseekfind(
    path_taryeast_structure="",
    name="",
    refname="",
    *,
    allow_score_fallback=True,
    require_fresh_output=False,
):
    root = Path.cwd().resolve(strict=True)
    target_structures = (
        Path(path_taryeast_structure).resolve()
        if path_taryeast_structure
        else (root / "struct_data" / "taryeast" / name).resolve()
    )
    reference_structures = root / "struct_data" / refname
    if not target_structures.is_dir():
        raise NotADirectoryError(
            f"Target structure directory not found: {target_structures}"
        )
    if not reference_structures.is_dir():
        raise NotADirectoryError(
            f"Reference structure directory not found: {reference_structures}"
        )

    working = root / "working" / name
    if require_fresh_output:
        _assert_formal_output_path(
            working,
            root,
            "working directory",
        )
    working.mkdir(parents=True, exist_ok=True)
    if require_fresh_output:
        _assert_formal_output_path(
            working,
            root,
            "working directory",
        )
        if not working.is_dir():
            raise RuntimeError(
                f"Formal Foldseek working path is not a directory: {working}"
            )
    raw_output = working / f"aln{name}.csv"
    matrix_output = working / f"matrix_foldseek_{name}.xlsx"
    audit_output = working / "foldseek_execution.json"
    if require_fresh_output:
        for path, label in (
            (raw_output, "raw output"),
            (matrix_output, "matrix output"),
            (audit_output, "execution audit"),
        ):
            _assert_formal_output_path(path, root, label)
            _reject_unsafe_existing_output(path, label)
            if path.exists():
                path.unlink()

    started = time.time()
    started_at = datetime.now(timezone.utc).isoformat()
    score_fields = (
        ("alntmscore", "qtmscore")
        if allow_score_fallback
        else ("alntmscore",)
    )
    successful_field = None
    successful_command = None
    last_error = None

    with tempfile.TemporaryDirectory(
        dir=working,
        prefix=".foldseek-run-",
    ) as run_directory:
        run_root = Path(run_directory)
        database_prefix = run_root / "foldseek_db" / "target"
        database_prefix.parent.mkdir(parents=True)
        createdb_command = [
            "foldseek",
            "createdb",
            str(target_structures),
            str(database_prefix),
        ]
        subprocess.run(createdb_command, check=True)

        for score_field in score_fields:
            with tempfile.TemporaryDirectory(
                dir=run_root,
                prefix=f"{score_field}-",
            ) as attempt_directory:
                attempt_root = Path(attempt_directory)
                attempt_output = attempt_root / "alignments.tsv"
                command = [
                    "foldseek",
                    "easy-search",
                    str(reference_structures),
                    str(database_prefix),
                    str(attempt_output),
                    str(attempt_root / "tmp"),
                    "--format-output",
                    f"query,target,{score_field},prob,qcov,tcov",
                    "--alignment-type",
                    "2",
                    "-e",
                    "0.00001",
                ]
                try:
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as error:
                    if not allow_score_fallback:
                        raise
                    last_error = error
                    continue
                if (
                    not attempt_output.is_file()
                    or attempt_output.is_symlink()
                    or attempt_output.stat().st_size == 0
                ):
                    raise RuntimeError(
                        "Foldseek completed without a fresh non-empty regular "
                        f"alignment output: {attempt_output}"
                    )
                alignments = pd.read_csv(
                    attempt_output,
                    sep="\t",
                    names=["ref", "tar", "tms", "pro", "rcov", "tcov"],
                )
                _copy_atomic(attempt_output, raw_output)
                successful_field = score_field
                successful_command = command
                break

        if successful_field is None:
            raise RuntimeError(
                "foldseek easy-search failed for both supported TM-score fields"
            ) from last_error

    alignments = alignments[alignments["pro"] >= 1].copy()
    alignments["ref"] = alignments["ref"].map(
        lambda value: value[3 : value.rfind("-F1-")]
    )
    alignments["tar"] = alignments["tar"].map(
        lambda value: value[3 : value.rfind("-F1-")]
    )
    alignments["pro"] = alignments["tms"].astype(float)
    alignments["tms"] = 1.0
    if require_fresh_output:
        _write_excel_atomic(alignments, matrix_output)
        audit = {
            "schema_version": 1,
            "started_at_utc": started_at,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "allow_score_fallback": allow_score_fallback,
            "score_field": successful_field,
            "createdb_command": createdb_command,
            "easy_search_command": successful_command,
            "raw_output": _file_record(raw_output),
            "matrix_output": _file_record(matrix_output),
        }
        _write_json_atomic(audit, audit_output)
    else:
        alignments.to_excel(matrix_output, index=False)
    print(
        "foldseek easy-search succeeded with score field: "
        f"{successful_field}"
    )
    print(f"Foldseek cost {time.time() - started:.2f} seconds")
