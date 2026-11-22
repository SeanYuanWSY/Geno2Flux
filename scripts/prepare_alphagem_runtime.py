#!/usr/bin/env python3
"""Install the repository-tracked AlphaGEM source into an existing runtime root.

This copies ordinary files and commits a local source-only Git snapshot. It
never downloads resources, creates symlinks, overwrites an existing file, or
starts an AlphaGEM experiment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BUNDLE_ROOT = REPOSITORY_ROOT / "third_party" / "alphagem"
MANIFEST_PATH = BUNDLE_ROOT / "SOURCE_MANIFEST.json"
RUNTIME_GITIGNORE = """# Generated runtime: only the frozen source snapshot is committed.
/data_available/
/rhea/
/tools/
/plmsearchtools/
/working/
/temp/
/struct_data/
/var/
/models/*
__pycache__/
*.pyc
.DS_Store
"""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_real_runtime_path(runtime_root: Path) -> None:
    if any(path.is_symlink() for path in (runtime_root, *runtime_root.parents)):
        raise ValueError("Runtime path contains a symlink")
    if runtime_root.resolve(strict=False).is_relative_to(REPOSITORY_ROOT.resolve()):
        raise ValueError("Runtime root must be outside the source repository")


def source_records() -> tuple[dict, list[dict]]:
    if not BUNDLE_ROOT.is_dir() or BUNDLE_ROOT.is_symlink():
        raise ValueError("Bundled AlphaGEM source is missing or is a symlink")
    if any(path.is_symlink() for path in BUNDLE_ROOT.rglob("*")):
        raise ValueError("Bundled AlphaGEM source contains a symlink")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or not isinstance(
        manifest.get("files"), list
    ):
        raise ValueError("Invalid AlphaGEM source manifest")
    records = manifest["files"]
    if not records:
        raise ValueError("Empty AlphaGEM source manifest")
    seen: set[str] = set()
    for item in records:
        relative = Path(item["path"])
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or not relative.parts
            or relative.as_posix() in seen
        ):
            raise ValueError(f"Unsafe or duplicate source path: {relative}")
        seen.add(relative.as_posix())
        source = BUNDLE_ROOT / relative
        if (
            not source.is_file()
            or source.stat().st_size != item["bytes"]
            or sha256(source) != item["sha256"]
        ):
            raise ValueError(f"Bundled AlphaGEM source differs: {relative}")
    return manifest, records


def inspect_runtime(runtime_root: Path, records: list[dict]) -> dict:
    require_real_runtime_path(runtime_root)
    missing = []
    verified = []
    for item in records:
        relative = Path(item["path"])
        destination = runtime_root / relative
        if any(
            (runtime_root / Path(*relative.parts[:index])).is_symlink()
            for index in range(1, len(relative.parts))
        ):
            raise ValueError(f"Runtime path contains a symlink: {relative}")
        if destination.is_symlink():
            raise ValueError(f"Runtime source file is a symlink: {relative}")
        if not destination.exists():
            missing.append(relative.as_posix())
        elif not destination.is_file() or sha256(destination) != item["sha256"]:
            raise ValueError(f"Runtime source file conflicts: {relative}")
        else:
            verified.append(relative.as_posix())
    return {"verified": verified, "missing": missing}


def _git(runtime_root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", *arguments],
        cwd=runtime_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _write_control_file(path: Path, content: str) -> None:
    if path.is_symlink():
        raise ValueError(f"Runtime control file is a symlink: {path.name}")
    if path.exists():
        if not path.is_file() or path.read_text(encoding="utf-8") != content:
            raise ValueError(f"Runtime control file conflicts: {path.name}")
        return
    try:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(content)
    except FileExistsError as error:
        raise ValueError(f"Runtime control file appeared: {path.name}") from error


def runtime_source_commit(runtime_root: Path, records: list[dict]) -> str:
    """Require a clean, committed local Git identity for the formal profile."""
    require_real_runtime_path(runtime_root)
    if (runtime_root / ".git").is_symlink():
        raise ValueError("Runtime Git metadata is a symlink")
    if not (runtime_root / ".git").is_dir():
        raise ValueError("Runtime source has no Git commit; use --prepare first")
    if (runtime_root / ".gitignore").is_symlink() or (
        runtime_root / "SOURCE_IDENTITY.json"
    ).is_symlink():
        raise ValueError("Runtime source identity metadata contains a symlink")
    if (runtime_root / ".gitignore").read_text(encoding="utf-8") != RUNTIME_GITIGNORE:
        raise ValueError("Runtime source ignore rules differ")
    identity = json.loads((runtime_root / "SOURCE_IDENTITY.json").read_text(encoding="utf-8"))
    expected = {
        "schema_version": 1,
        "upstream_commit": json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["source_commit"],
        "source_manifest_sha256": sha256(MANIFEST_PATH),
    }
    if identity != expected:
        raise ValueError("Runtime source identity differs from bundled snapshot")
    if inspect_runtime(runtime_root, records)["missing"]:
        raise ValueError("Runtime source is incomplete")
    if Path(_git(runtime_root, "rev-parse", "--show-toplevel")).resolve() != runtime_root.resolve():
        raise ValueError("Runtime Git repository is not rooted at runtime root")
    commit = _git(runtime_root, "rev-parse", "HEAD")
    if _git(runtime_root, "status", "--porcelain=v1", "--untracked-files=all"):
        raise ValueError("Runtime source Git repository is dirty")
    tracked = set(_git(runtime_root, "ls-files").splitlines())
    required = {item["path"] for item in records} | {".gitignore", "SOURCE_IDENTITY.json"}
    if not required.issubset(tracked):
        raise ValueError("Runtime Git commit does not track the complete source snapshot")
    return commit


def ensure_runtime_source_commit(runtime_root: Path, records: list[dict]) -> str:
    """Create a local source-only commit, never a remote or a symlink."""
    if (runtime_root / ".git").is_symlink():
        raise ValueError("Runtime Git metadata is a symlink")
    git_directory = runtime_root / ".git"
    if git_directory.exists() and not (runtime_root / "SOURCE_IDENTITY.json").is_file():
        raise ValueError("Existing runtime Git repository lacks source identity")
    identity = {
        "schema_version": 1,
        "upstream_commit": json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["source_commit"],
        "source_manifest_sha256": sha256(MANIFEST_PATH),
    }
    _write_control_file(runtime_root / ".gitignore", RUNTIME_GITIGNORE)
    _write_control_file(
        runtime_root / "SOURCE_IDENTITY.json",
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
    )
    if not git_directory.exists():
        _git(runtime_root, "init", "-q")
    elif not git_directory.is_dir():
        raise ValueError("Runtime Git metadata is not a directory")
    try:
        _git(runtime_root, "rev-parse", "--verify", "HEAD")
    except subprocess.CalledProcessError:
        _git(
            runtime_root,
            "add", "-f", "--", ".gitignore", "SOURCE_IDENTITY.json",
            *(item["path"] for item in records),
        )
        _git(
            runtime_root,
            "-c", "user.name=Geno2Flux Runtime",
            "-c", "user.email=runtime@localhost.invalid",
            "commit", "-q", "-m", "Freeze AlphaGEM source snapshot for formal runs",
        )
    return runtime_source_commit(runtime_root, records)


def install_runtime(runtime_root: Path, records: list[dict]) -> dict:
    require_real_runtime_path(runtime_root)
    runtime_root.mkdir(parents=True, exist_ok=True)
    before = inspect_runtime(runtime_root, records)
    installed = []
    for item in records:
        relative = Path(item["path"])
        if relative.as_posix() not in before["missing"]:
            continue
        destination = runtime_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.parent.is_symlink():
            raise ValueError(f"Runtime parent is a symlink: {relative}")
        with tempfile.NamedTemporaryFile(
            prefix=f".{destination.name}.incoming-",
            dir=destination.parent,
            delete=False,
        ) as temporary:
            temp_path = Path(temporary.name)
        try:
            shutil.copyfile(BUNDLE_ROOT / relative, temp_path)
            shutil.copymode(BUNDLE_ROOT / relative, temp_path)
            if sha256(temp_path) != item["sha256"]:
                raise ValueError(f"Copied AlphaGEM source differs: {relative}")
            os.link(temp_path, destination)  # atomic create; never overwrite
            installed.append(relative.as_posix())
        finally:
            temp_path.unlink(missing_ok=True)
    after = inspect_runtime(runtime_root, records)
    if after["missing"]:
        raise RuntimeError("Runtime source installation is incomplete")
    commit = ensure_runtime_source_commit(runtime_root, records)
    return {"installed": installed, "verified": after["verified"], "runtime_source_commit": commit}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify or copy bundled AlphaGEM source without symlinks."
    )
    parser.add_argument("--runtime-root", required=True, type=Path)
    parser.add_argument(
        "--prepare",
        action="store_true",
        help="Copy verified source files; without this flag, inspect only.",
    )
    args = parser.parse_args(argv)
    manifest, records = source_records()
    # Canonicalize command-line aliases (macOS maps /tmp to /private/tmp)
    # before applying the destination no-symlink safety check.
    runtime_root = args.runtime_root.expanduser().resolve(strict=False)
    require_real_runtime_path(runtime_root)
    result = (
        install_runtime(runtime_root, records)
        if args.prepare
        else inspect_runtime(runtime_root, records)
    )
    if not args.prepare and not result["missing"]:
        result["runtime_source_commit"] = runtime_source_commit(runtime_root, records)
    print(
        json.dumps(
            {
                "source_commit": manifest["source_commit"],
                "source_files": len(records),
                "runtime_root": str(runtime_root),
                "prepared": args.prepare,
                **result,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if not result.get("missing") else 2


if __name__ == "__main__":
    raise SystemExit(main())
