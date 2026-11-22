import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from src.pipeline_contracts import protein_key
def _clean_path():
    return Path.cwd() / "tools" / "CLEAN" / "app"


def _assert_local_path(path, root, *, label, require_file=False):
    root = Path(root).resolve(strict=True)
    path = Path(path)
    if not path.is_absolute():
        path = root / path
    try:
        relative = path.relative_to(root)
    except ValueError as error:
        raise RuntimeError(f"{label} escapes the repository: {path}") from error
    current = root
    for component in relative.parts:
        current = current / component
        if current.is_symlink():
            raise RuntimeError(
                f"{label} path does not permit symlinks: {current}"
            )
    resolved = path.resolve(strict=require_file)
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise RuntimeError(
            f"{label} resolves outside the repository: {path} -> {resolved}"
        ) from error
    if require_file and (
        not path.is_file()
        or path.is_symlink()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"{label} is missing, empty, or not a regular file: {path}"
        )
    return path


def _file_record(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "bytes": Path(path).stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _runtime_file_record(path):
    path = Path(path)
    metadata = os.stat(path, follow_symlinks=False)
    return {
        **_file_record(path),
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "mtime_ns": metadata.st_mtime_ns,
        "ctime_ns": metadata.st_ctime_ns,
    }


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


def _copy_file_atomic(source, destination):
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


def _write_fasta_atomic(records, destination):
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
            written = SeqIO.write(records, handle, "fasta")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        return written
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def get_fasta(name, *, require_safe_paths=False):
    repo_root = Path.cwd().resolve(strict=True)
    clean_path = _clean_path()
    with open(f"working/{name}/{name}.fasta", encoding="utf-8") as fasta_file:
        fasta = list(SeqIO.parse(fasta_file, "fasta"))
    tary = pd.read_excel(f"working/{name}/{name}.xlsx")
    homo = pd.read_excel(f"working/{name}/matrix_homolog{name}.xlsx")
    target_ids = [protein_key(record.id) for record in fasta]
    table_ids = [str(tary.iat[index, 0]) for index in range(len(tary.index))]
    homolog_ids = {
        str(homo.iat[index, 2])
        for index in range(len(homo.index))
        if not pd.isna(homo.iat[index, 2])
    }
    if require_safe_paths:
        if len(target_ids) != len(set(target_ids)):
            raise RuntimeError(
                "Formal CLEAN membership input has duplicate target FASTA IDs"
            )
        if len(table_ids) != len(set(table_ids)) or set(table_ids) != set(
            target_ids
        ):
            raise RuntimeError(
                "Formal CLEAN membership table does not exactly cover the "
                "normalized target FASTA"
            )
        unknown_homologs = sorted(homolog_ids - set(target_ids))
        if unknown_homologs:
            raise RuntimeError(
                "Formal CLEAN homolog matrix references unknown target IDs: "
                f"{unknown_homologs[:20]}"
            )
    left = set(table_ids) - homolog_ids
    records_left = [
        record for record in fasta if protein_key(record.id) in left
    ]
    if require_safe_paths and len(records_left) != len(left):
        raise RuntimeError(
            "Formal CLEAN non-homolog selection is not one-to-one with the "
            "validated target artifacts"
        )
    working_input = Path("working") / name / f"{name}_homoleft.fasta"
    if require_safe_paths:
        _assert_local_path(
            working_input.parent,
            repo_root,
            label="Formal CLEAN working directory",
        )
        working_input.parent.mkdir(parents=True, exist_ok=True)
        _assert_local_path(
            working_input.parent,
            repo_root,
            label="Formal CLEAN working directory",
        )
        if not working_input.parent.is_dir():
            raise RuntimeError(
                "Formal CLEAN working directory is not a directory: "
                f"{working_input.parent}"
            )
        _assert_local_path(
            working_input,
            repo_root,
            label="Formal CLEAN working input FASTA",
        )
        if working_input.is_symlink():
            raise RuntimeError(
                "Formal CLEAN working input FASTA does not permit symlinks: "
                f"{working_input}"
            )
        if working_input.exists() and not working_input.is_file():
            raise RuntimeError(
                "Formal CLEAN working input FASTA must be a regular file "
                f"when it exists: {working_input}"
            )
        working_count = _write_fasta_atomic(records_left, working_input)
        if working_count != len(records_left):
            raise RuntimeError(
                "Formal CLEAN input staging did not preserve every selected "
                "non-homolog record"
            )
    else:
        clean_input = (
            clean_path / "data" / "inputs" / f"{name}_homoleft.fasta"
        )
        clean_input.parent.mkdir(parents=True, exist_ok=True)
        with clean_input.open("w", encoding="utf-8") as fasta_out:
            SeqIO.write(records_left, fasta_out, "fasta")
        with working_input.open("w", encoding="utf-8") as fasta_out:
            SeqIO.write(records_left, fasta_out, "fasta")


def _prepare_isolated_clean_runtime_impl(
    repo_root,
    clean_path,
    working,
    input_fasta,
    name,
    esm_model_path,
    esm_token_budget,
    failed_runtimes,
):
    working = Path(working).resolve(strict=True)
    runtime = Path(
        tempfile.mkdtemp(
            prefix="clean-runtime-",
            dir=working,
        )
    ).resolve(strict=True)
    failed_runtimes.append(runtime)
    _assert_local_path(
        runtime,
        repo_root,
        label="Formal CLEAN isolated runtime",
    )
    if runtime.is_symlink() or not runtime.is_dir():
        raise RuntimeError(
            f"Formal CLEAN isolated runtime is invalid: {runtime}"
        )
    if (
        isinstance(esm_token_budget, bool)
        or not isinstance(esm_token_budget, int)
        or esm_token_budget <= 0
    ):
        raise ValueError(
            "Formal CLEAN ESM token budget must be a positive integer"
        )
    runtime_input = (
        runtime / "data" / "inputs" / f"{name}_homoleft.fasta"
    )
    runtime_input.parent.mkdir(parents=True)
    (runtime / "data" / "esm_data").mkdir()
    (runtime / "results" / "inputs").mkdir(parents=True)
    _copy_file_atomic(input_fasta, runtime_input)
    esm_model_path = Path(esm_model_path)
    if esm_model_path.name != "esm1b_t33_650M_UR50S.pt":
        raise RuntimeError(
            "Formal CLEAN requires the canonical ESM-1b checkpoint filename: "
            "esm1b_t33_650M_UR50S.pt"
        )
    if esm_model_path.is_symlink():
        raise RuntimeError(
            "Formal CLEAN ESM-1b model weights do not permit symlinks: "
            f"{esm_model_path}"
        )
    esm_model_path = esm_model_path.resolve(strict=True)
    esm_regression_path = esm_model_path.with_name(
        esm_model_path.stem + "-contact-regression.pt"
    )
    if esm_regression_path.is_symlink():
        raise RuntimeError(
            "Formal CLEAN ESM-1b contact-regression weights do not permit "
            f"symlinks: {esm_regression_path}"
        )
    esm_regression_path = esm_regression_path.resolve(strict=True)
    for path, label in (
        (esm_model_path, "Formal CLEAN ESM-1b model weights"),
        (
            esm_regression_path,
            "Formal CLEAN ESM-1b contact-regression weights",
        ),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"{label} is missing or empty: {path}")
    checkpoint_dir = esm_model_path.parent
    if (
        checkpoint_dir.name != "checkpoints"
        or checkpoint_dir.parent.name != "hub"
    ):
        raise RuntimeError(
            "Formal CLEAN ESM weights must be inside a pinned "
            "TORCH_HOME/hub/checkpoints directory"
        )
    torch_home = checkpoint_dir.parent.parent
    weight_records = {
        "source_model": {
            "observed_path": str(esm_model_path),
            **_runtime_file_record(esm_model_path),
        },
        "source_regression": {
            "observed_path": str(esm_regression_path),
            **_runtime_file_record(esm_regression_path),
        },
    }

    static_links = []

    def link(relative, target):
        target = Path(target)
        _assert_local_path(
            target,
            repo_root,
            label=f"Formal CLEAN static dependency {relative}",
            require_file=target.is_file(),
        )
        if not target.exists() or target.is_symlink():
            raise RuntimeError(
                "Formal CLEAN static dependency is missing or a symlink: "
                f"{target}"
            )
        destination = runtime / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.symlink_to(
            target,
            target_is_directory=target.is_dir(),
        )
        if not destination.is_symlink() or destination.resolve() != (
            target.resolve()
        ):
            raise RuntimeError(
                f"Formal CLEAN static link is invalid: {destination}"
            )
        static_links.append(
            {
                "path": str(destination),
                "target": str(target.resolve()),
            }
        )

    link(Path("src"), clean_path / "src")
    canonical_extract = clean_path / "esm" / "scripts" / "extract.py"
    _assert_local_path(
        canonical_extract,
        repo_root,
        label="Formal CLEAN canonical ESM extraction script",
        require_file=True,
    )
    source_wrapper = repo_root / "src" / "clean_esm_extract_wrapper.py"
    _assert_local_path(
        source_wrapper,
        repo_root,
        label="Formal CLEAN ESM extraction wrapper source",
        require_file=True,
    )
    runtime_extract = runtime / "esm" / "scripts" / "extract.py"
    runtime_extract.parent.mkdir(parents=True)
    _copy_file_atomic(source_wrapper, runtime_extract)
    esm_extraction = {
        "token_budget": esm_token_budget,
        "canonical_extract": {
            "observed_path": str(canonical_extract),
            **_runtime_file_record(canonical_extract),
        },
        "source_wrapper": {
            "observed_path": str(source_wrapper),
            **_runtime_file_record(source_wrapper),
        },
        "runtime_wrapper": {
            "observed_path": str(runtime_extract),
            **_runtime_file_record(runtime_extract),
        },
    }
    canonical_data = clean_path / "data"
    for child in sorted(canonical_data.iterdir()):
        if child.name in {"inputs", "esm_data"}:
            continue
        link(Path("data") / child.name, child)
    return (
        runtime,
        runtime_input,
        static_links,
        weight_records,
        esm_extraction,
        torch_home,
    )


def _remove_failed_clean_runtimes(runtimes):
    failures = []
    for runtime in reversed(runtimes):
        try:
            shutil.rmtree(runtime)
        except OSError as error:
            failures.append(f"{runtime}: {error}")
    if failures:
        raise RuntimeError(
            "Formal CLEAN failed and its isolated runtime could not be "
            "removed: " + "; ".join(failures)
        )


def _prepare_isolated_clean_runtime(
    repo_root,
    clean_path,
    working,
    input_fasta,
    name,
    esm_model_path,
    esm_token_budget,
):
    failed_runtimes = []
    try:
        prepared = _prepare_isolated_clean_runtime_impl(
            repo_root,
            clean_path,
            working,
            input_fasta,
            name,
            esm_model_path,
            esm_token_budget,
            failed_runtimes,
        )
    except BaseException as original_error:
        try:
            _remove_failed_clean_runtimes(failed_runtimes)
        except BaseException as cleanup_error:
            raise original_error from cleanup_error
        raise
    failed_runtimes.clear()
    return prepared


def _verify_static_links(static_links, clean_path):
    clean_path = Path(clean_path).resolve(strict=True)
    for record in static_links:
        path = Path(record["path"])
        target = Path(record["target"])
        if (
            not path.is_symlink()
            or path.resolve(strict=True) != target
            or target != target.resolve(strict=True)
        ):
            raise RuntimeError(
                f"Formal CLEAN static link changed during execution: {path}"
            )
        try:
            target.relative_to(clean_path)
        except ValueError as error:
            raise RuntimeError(
                "Formal CLEAN static link escapes the canonical application: "
                f"{path} -> {target}"
            ) from error


def _clean_impl(
    name,
    *,
    require_fresh_output=False,
    esm_model_path=None,
    esm_token_budget=4096,
    failed_runtimes,
):
    repo_root = Path.cwd().resolve(strict=True)
    clean_path = _clean_path()
    script = clean_path / "CLEAN_infer_fasta.py"
    if require_fresh_output:
        _assert_local_path(
            clean_path,
            repo_root,
            label="Formal CLEAN application root",
        )
        _assert_local_path(
            script,
            repo_root,
            label="Formal CLEAN entry point",
            require_file=True,
        )
    elif not script.is_file():
        raise FileNotFoundError(f"CLEAN entry point not found: {script}")
    working = Path("working") / name
    working_input = working / f"{name}_homoleft.fasta"
    canonical_input = (
        clean_path / "data" / "inputs" / f"{name}_homoleft.fasta"
    )
    input_fasta = working_input if require_fresh_output else canonical_input
    if require_fresh_output:
        _assert_local_path(
            input_fasta,
            repo_root,
            label="Formal CLEAN input FASTA",
            require_file=True,
        )
    elif not input_fasta.is_file() or input_fasta.stat().st_size == 0:
        raise FileNotFoundError(f"CLEAN input FASTA is missing or empty: {input_fasta}")
    canonical_result = (
        clean_path
        / "results"
        / "inputs"
        / f"{name}_homoleft_maxsep.csv"
    )
    destination = working / f"{name}_homoleft_maxsep.csv"
    audit_path = working / "clean_execution.json"
    runtime = None
    static_links = []
    weight_records = None
    esm_extraction = None
    if require_fresh_output:
        if esm_model_path is None:
            raise ValueError(
                "Formal CLEAN execution requires a local ESM-1b weight file"
            )
        (
            runtime,
            runtime_input,
            static_links,
            weight_records,
            esm_extraction,
            torch_home,
        ) = (
            _prepare_isolated_clean_runtime(
                repo_root,
                clean_path,
                working,
                input_fasta,
                name,
                esm_model_path,
                esm_token_budget,
            )
        )
        failed_runtimes.append(runtime)
        input_fasta = runtime_input
        result = (
            runtime
            / "results"
            / "inputs"
            / f"{name}_homoleft_maxsep.csv"
        )
    else:
        result = canonical_result
    if require_fresh_output:
        for path, label in (
            (result, "Formal CLEAN result"),
            (destination, "Formal CLEAN working result"),
            (audit_path, "Formal CLEAN execution audit"),
        ):
            _assert_local_path(path, repo_root, label=label)
            if path.is_symlink():
                raise RuntimeError(f"{label} does not permit symlinks: {path}")
            if path.exists() and not path.is_file():
                raise RuntimeError(
                    f"{label} must be a regular file when it exists: {path}"
                )
    removed_preexisting_result = False
    command = [
        sys.executable,
        str(script),
        "--fasta_data",
        f"{name}_homoleft",
    ]
    started_at = datetime.now(timezone.utc).isoformat()
    entrypoint_record = _runtime_file_record(script)
    input_record = _runtime_file_record(input_fasta)
    working_input_record = (
        _file_record(working_input) if require_fresh_output else None
    )
    environment = None
    nested_python_record = None
    command_cwd = clean_path
    if require_fresh_output:
        environment = os.environ.copy()
        environment["PYTHONPATH"] = os.pathsep.join(
            (
                str(clean_path / "src"),
                str(clean_path / "esm"),
            )
        )
        environment["PYTHONNOUSERSITE"] = "1"
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        environment["TORCH_HOME"] = str(torch_home)
        environment["ALPHAGEM_CLEAN_ESM_TOKEN_BUDGET"] = str(
            esm_token_budget
        )
        environment["ALPHAGEM_CLEAN_ESM_EXTRACT"] = (
            esm_extraction["canonical_extract"]["observed_path"]
        )
        for variable in (
            "HTTP_PROXY",
            "HTTPS_PROXY",
            "http_proxy",
            "https_proxy",
        ):
            environment[variable] = "http://127.0.0.1:1"
        environment["NO_PROXY"] = ""
        environment["no_proxy"] = ""
        environment["PATH"] = os.pathsep.join(
            [
                str(Path(sys.executable).parent),
                environment.get("PATH", ""),
            ]
        )
        nested_python = shutil.which("python", path=environment["PATH"])
        nested_python_path = (
            Path(nested_python).absolute()
            if nested_python is not None
            else None
        )
        if (
            nested_python_path is None
            or nested_python_path.parent
            != Path(sys.executable).absolute().parent
            or nested_python_path.resolve()
            != Path(sys.executable).resolve()
        ):
            raise RuntimeError(
                "Formal CLEAN nested ESM extraction cannot bind `python` "
                f"to the active interpreter: {nested_python!r}"
            )
        nested_python_record = {
            "invocation_path": str(nested_python_path),
            "resolved_path": str(nested_python_path.resolve()),
            "sys_prefix": sys.prefix,
            **_runtime_file_record(nested_python_path.resolve()),
        }
        command_cwd = runtime
    subprocess.run(
        command,
        cwd=command_cwd,
        check=True,
        **({"env": environment} if environment is not None else {}),
    )
    if require_fresh_output:
        _assert_local_path(
            script,
            repo_root,
            label="Formal CLEAN entry point",
            require_file=True,
        )
        _assert_local_path(
            input_fasta,
            repo_root,
            label="Formal CLEAN input FASTA",
            require_file=True,
        )
        if (
            _runtime_file_record(script) != entrypoint_record
            or _runtime_file_record(input_fasta) != input_record
        ):
            raise RuntimeError(
                "Formal CLEAN entry point or input changed while CLEAN was "
                "executing"
            )
        _verify_static_links(static_links, clean_path)
        for record in (
            esm_extraction["canonical_extract"],
            esm_extraction["source_wrapper"],
            esm_extraction["runtime_wrapper"],
        ):
            path = Path(record["observed_path"])
            if _runtime_file_record(path) != {
                key: value
                for key, value in record.items()
                if key != "observed_path"
            }:
                raise RuntimeError(
                    "Formal CLEAN ESM extraction controller changed during "
                    f"execution: {path}"
                )
        for record in weight_records.values():
            path = Path(record["observed_path"])
            if _runtime_file_record(path) != {
                key: value
                for key, value in record.items()
                if key != "observed_path"
            }:
                raise RuntimeError(
                    "Formal CLEAN ESM weight changed during execution: "
                    f"{path}"
                )
    if not result.is_file() or result.is_symlink() or result.stat().st_size == 0:
        raise RuntimeError(
            f"CLEAN completed without a fresh non-empty regular result: {result}"
        )
    if require_fresh_output:
        _assert_local_path(
            result,
            repo_root,
            label="Formal CLEAN result",
            require_file=True,
        )
        result_record = _runtime_file_record(result)
        _copy_file_atomic(result, destination)
        if _runtime_file_record(result) != result_record:
            raise RuntimeError(
                "Formal CLEAN result changed while it was copied to the "
                "working directory"
            )
    else:
        shutil.copy2(result, destination)
    if require_fresh_output:
        audit = {
            "schema_version": 3,
            "command": command,
            "command_cwd": str(command_cwd),
            "started_at_utc": started_at,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "removed_preexisting_result": removed_preexisting_result,
            "isolated_runtime": True,
            "runtime_directory": str(runtime),
            "static_links": static_links,
            "nested_python": nested_python_record,
            "esm_weights": weight_records,
            "esm_weight_storage": "direct_pinned_source",
            "esm_extraction": esm_extraction,
            "torch_home": str(torch_home),
            "network_policy": "loopback_refuse_proxy",
            "pythonpath": environment["PYTHONPATH"],
            "python_no_user_site": environment["PYTHONNOUSERSITE"],
            "clean_entrypoint": entrypoint_record,
            "input_fasta": input_record,
            "working_input": working_input_record,
            "clean_result": result_record,
            "working_result": _file_record(destination),
        }
        if any(
            audit["clean_result"][field] != audit["working_result"][field]
            for field in ("bytes", "sha256")
        ):
            raise RuntimeError(
                "Formal CLEAN working copy differs from the newly generated result"
            )
        _write_json_atomic(audit_path, audit)
        failed_runtimes.remove(runtime)


def clean(
    name,
    *,
    require_fresh_output=False,
    esm_model_path=None,
    esm_token_budget=4096,
):
    failed_runtimes = []
    try:
        return _clean_impl(
            name,
            require_fresh_output=require_fresh_output,
            esm_model_path=esm_model_path,
            esm_token_budget=esm_token_budget,
            failed_runtimes=failed_runtimes,
        )
    except BaseException as original_error:
        try:
            _remove_failed_clean_runtimes(failed_runtimes)
        except BaseException as cleanup_error:
            raise original_error from cleanup_error
        raise


def clean_result(name):
    get_fasta(name)
    clean(name)
