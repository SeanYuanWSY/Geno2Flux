#!/usr/bin/env python3
"""One-click fetch + preflight for Geno2Flux external tools, data, and envs.

Profiles compose: ``gapfill`` implies ``screening`` and ``alphagem`` implies
``gapfill``; ``all`` expands to every profile. Downloads are resumable and
every artifact with a published checksum is verified before it is moved into
place. ``--check`` performs a real readiness audit (imports, interpreters,
binaries, files) instead of fetching.

Examples:
    # Preview the full plan without touching disk or network:
    bash run.sh fetch-resources --profile all --dry-run

    # DIAMOND into tools/bin/ (Linux x86_64 tarball; conda/brew hints on macOS):
    bash run.sh fetch-resources --profile screening

    # AlphaGEM data + ESM-1b weights into an external runtime root:
    bash run.sh fetch-resources --profile alphagem \
        --runtime-root /work/alphagem-runtime --with-eggnog

    # ANNEVO source + fungal weights:
    bash run.sh fetch-resources --profile annotation --tool-root /work/annotation-tools

    # Real readiness audit of every stage:
    bash run.sh check-env --runtime-root /work/alphagem-runtime \
        --alphagem-python /path/to/AlphaGEM-env/bin/python
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import posixpath
import shlex
import shutil
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BIN_DIR = REPO_ROOT / "tools" / "bin"
DEEPLOC_VENV = REPO_ROOT / "tools" / "deeploc"
SUPPORTED_PLATFORMS = (
    "linux-x86_64",
    "linux-aarch64",
    "darwin-aarch64",
    "darwin-x86_64",
)

# ---------------------------------------------------------------------------
# Official, live-verified download specifications. Do not substitute mirrors.
# ---------------------------------------------------------------------------

DIAMOND_VERSION = "2.1.10"
DIAMOND_URLS = {
    # No official macOS/aarch64 GitHub release asset exists.
    "linux-x86_64": (
        "https://github.com/bbuchfink/diamond/releases/download/"
        "v2.1.10/diamond-linux64.tar.gz"
    ),
}
DIAMOND_MANUAL_HINT = (
    "No official DIAMOND binary release exists for this platform. Install "
    "DIAMOND 2.1.10 with either:\n"
    "    conda install -c conda-forge -c bioconda diamond=2.1.10\n"
    "    brew install diamond\n"
    "then re-run `bash run.sh check-env --profile screening`."
)

BLAST_PLUS = {
    "version": "2.16.0",
    "base": "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0",
    "archives": {
        "linux-x86_64": {
            "file": "ncbi-blast-2.16.0+-x64-linux.tar.gz",
            "md5": "48f66c9e01ea5136e381b2bf6fc62036",
        },
        "linux-aarch64": {
            "file": "ncbi-blast-2.16.0+-aarch64-linux.tar.gz",
            "md5": "d5ae1faf81bb85ad872d6efa4a0079e9",
        },
        "darwin-aarch64": {
            "file": "ncbi-blast-2.16.0+-aarch64-macosx.tar.gz",
            "md5": "13ce616ef80afd7e3e4c5f44ec57a521",
        },
        "darwin-x86_64": {
            "file": "ncbi-blast-2.16.0+-x64-macosx.tar.gz",
            "md5": "7888e1ccf54b5add921867b7ec0cfd69",
        },
    },
    "binaries": ("tblastn", "makeblastdb", "blastp"),
}

ALPHAGEM_DATA_SOURCE = {
    "url": "https://ndownloader.figshare.com/files/59612702",
    "doi": "https://doi.org/10.6084/m9.figshare.30630524",
    "bytes": 13073797016,
    "md5": "73aaaec107c5faec9f4c642706de5b59",
    "sha256": "1ef308718a647ce339fe1423aac573489567f3c0d3e2f155e42fa1a954bd5875",
}

ESM_WEIGHTS = {
    "esm1b_t33_650M_UR50S.pt": {
        "url": "https://dl.fbaipublicfiles.com/fair-esm/models/"
        "esm1b_t33_650M_UR50S.pt",
        "bytes": 7828576466,
    },
    "esm1b_t33_650M_UR50S-contact-regression.pt": {
        "url": "https://dl.fbaipublicfiles.com/fair-esm/regression/"
        "esm1b_t33_650M_UR50S-contact-regression.pt",
        "bytes": 3174,
        "md5": "e7fe626dfd516fb6824bd1d30192bdb1",
    },
}

# eggNOG 5.0.2 contracts mirrored from
# third_party/alphagem/src/paper_tool_contracts.py
EGGNOG = {
    "db_gz_url": "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz",
    "db_sha256": "4d6857c129423889edcd607995b28d97944787a588224f31dd7f745fe9ef1e4f",
    "db_bytes": 41370988544,
    "proteomes_url": "http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa",
    "proteomes_sha256": "efc3685ebcbfb1d73359231dbb823cd4cf8bfc80c1ca1b98c7117d1c4674b8f2",
    "taxid_url": "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv",
    "taxid_sha256": "c3839fada0046a1c83d2eda6625a2c806931e3e982b84c43b7147dbe8474a125",
    "selected_fasta_sha256": (
        "0f64c4114fbcd4bbd78bb7493854226c0baeb5f5dfc9fe6768da8e41a0d7da21"
    ),
    "fungi_dmnd_sha256": (
        "b9c737213a3713af370dbc03e9d798ea8ef9a576754d5707d7a7c75ec9546c16"
    ),
    "fungi_dmnd_bytes": 879076733,
    "selection_contract": (
        "eggnog-mapper 2.1.12 create_dbs.py: retain taxids whose Named "
        "Lineage (4th column, comma-separated) contains the exact token "
        "'Fungi', then retain e5.proteomes.faa records whose header taxid "
        "precedes the first dot"
    ),
}

# 9-file eggNOG static data bundle mirrored from
# third_party/alphagem/src/paper_tool_contracts.py
# (PAPER_TOOL_BUNDLE_CONTRACTS["eggnog_static_data"]).
EGGNOG_STATIC_DATA = {
    "path": "tools/eggnog_mapper/data",
    "excluded_relative_roots": ("eggnog.db", "eggnog.db.gz", "fungi.dmnd"),
    "file_count": 9,
    "bytes": 16224797454,
    "sha256": (
        "2128f9ae1e5cab076409822c04fc5b83f9c75d977513c87f92a098b6a1edf5f7"
    ),
}


PROTT5_REPO = "Rostlab/prot_t5_xl_uniref50"
PROTT5_URL = f"https://huggingface.co/{PROTT5_REPO}"

ANNEVO = {
    "archive_url": "https://github.com/xjtu-omics/ANNEVO/archive/refs/heads/main.tar.gz",
    "weights_url": "https://raw.githubusercontent.com/xjtu-omics/ANNEVO/main/saved_model/ANNEVO_Fungi.pt",
    "weights_bytes": 58562884,
    "repo_url": "https://github.com/xjtu-omics/ANNEVO",
}

HELIXER_FUNGI_URL = (
    "https://zenodo.org/records/10836346/files/fungi_v0.3_a_0100.h5?download=1"
)

DEEPLOC = {
    "2.1": {
        "page": "https://services.healthtech.dtu.dk/services/DeepLoc-2.1/",
        "form": (
            "https://services.healthtech.dtu.dk/cgi-bin/sw_request"
            "?software=deeploc&version=2.1&packageversion=2.1&platform=All"
        ),
    },
    "2.0": {
        "page": "https://services.healthtech.dtu.dk/services/DeepLoc-2.0/",
        "form": (
            "https://services.healthtech.dtu.dk/cgi-bin/sw_request"
            "?software=deeploc&version=2.0&packageversion=2.0&platform=All"
        ),
    },
}

ALL_PROFILES = ("screening", "gapfill", "alphagem", "annotation", "deeploc")
# Ordered for the --check summary table.
CHECK_ORDER = ("core", "screening", "gapfill", "deeploc", "alphagem", "annotation")

_CORE_IMPORTS = ("cobra", "highspy", "scipy", "Bio", "pandas", "numpy", "openpyxl")
_ALPHAGEM_IMPORTS = ("torch", "esm", "sklearn", "cobra", "Bio", "pandas", "CLEAN|clean")
_DEEPEC_IMPORTS = ("torch", "transformers")
_ANNEVO_IMPORTS = ("torch", "h5py", "Bio")

# Runtime files required for a READY alphagem check: (relpath, origin hint).
_ALPHAGEM_REQUIRED_FILES: tuple[tuple[str, str], ...] = (
    ("AlphaGEM.py", "bundled source (fetch-resources --profile alphagem)"),
    ("models/yeast-GEM.xml", "bundled source"),
    ("models/reaction_pool_final.json", "bundled source"),
    ("data_available/yeast.fasta", "Figshare data_source.tar.xz"),
    ("data_available/yeast.xlsx", "Figshare data_source.tar.xz"),
    ("rhea/rhea2uniprot_sprot.tsv", "Figshare data_source.tar.xz"),
    ("rhea/rhea_protembedding1.pkl", "Figshare data_source.tar.xz"),
    ("rhea/rhea_protembedding2.pkl", "Figshare data_source.tar.xz"),
    ("rhea/rhea_protembedding3.pkl", "Figshare data_source.tar.xz"),
    ("rhea/rhea_protembedding4.pkl", "Figshare data_source.tar.xz"),
    ("rhea/rhea_protembedding5.pkl", "Figshare data_source.tar.xz"),
    ("rhea/rhea_protembedding6.pkl", "Figshare data_source.tar.xz"),
    ("rhea/rhea.rdf", "Figshare data_source.tar.xz"),
    ("plmsearchtools/model/plmsearch.sav", "Figshare data_source.tar.xz"),
    ("tools/CLEAN/app/CLEAN_infer_fasta.py", "Figshare data_source.tar.xz"),
    ("tools/DeepProZyme", "Figshare data_source.tar.xz"),
    ("tools/eggnog_mapper/emapper.py", "Figshare data_source.tar.xz"),
    ("tools/eggnog_mapper/data/eggnog.db", "--with-eggnog download"),
    ("tools/eggnog_mapper/data/fungi.dmnd", "--with-eggnog diamond makedb"),
)

_PROBE_IMPORTS_SNIPPET = """\
import importlib.util, sys
missing = []
for spec in sys.argv[1:]:
    if all(importlib.util.find_spec(name) is None for name in spec.split('|')):
        missing.append(spec)
print('MISSING:' + ','.join(missing) if missing else 'OK')
sys.exit(1 if missing else 0)
"""

_PROBE_HIGHS_MILP_SNIPPET = """\
import sys
try:
    import highspy
except ImportError:
    print('MISSING:highspy'); sys.exit(1)
h = highspy.Highs()
h.setOptionValue('output_flag', False)
lp = highspy.HighsLp()
lp.num_col_ = 2; lp.num_row_ = 1
lp.col_cost_ = [-1.0, 0.0]
lp.col_lower_ = [0.0, 0.0]; lp.col_upper_ = [1.0, 1.0]
lp.row_lower_ = [1.0]; lp.row_upper_ = [highspy.kHighsInf]
lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
lp.a_matrix_.start_ = [0, 1, 2]; lp.a_matrix_.index_ = [0, 0]
lp.a_matrix_.value_ = [1.0, 1.0]
lp.integrality_ = [highspy.HighsVarType.kInteger, highspy.HighsVarType.kContinuous]
assert h.passModel(lp) == highspy.HighsStatus.kOk
assert h.run() == highspy.HighsStatus.kOk
assert h.getInfo().objective_function_value == -1.0
print('MILP_OK')
"""


class FetchError(RuntimeError):
    """Download, checksum, extraction, or install failure."""


# ---------------------------------------------------------------------------
# Platform + hashing helpers
# ---------------------------------------------------------------------------


def detect_platform(override: str | None = None) -> str:
    """Return ``<os>-<arch>`` or raise ValueError for unsupported combos."""
    arch_aliases = {"arm64": "aarch64", "amd64": "x86_64"}
    if override:
        raw = override.strip().lower()
        if "-" in raw:
            os_part, arch_part = raw.split("-", 1)
            key = f"{os_part}-{arch_aliases.get(arch_part, arch_part)}"
        else:
            key = raw
    else:
        sysname = sys.platform
        if sysname.startswith("linux"):
            os_part = "linux"
        elif sysname == "darwin":
            os_part = "darwin"
        else:
            raise ValueError(f"Unsupported operating system: {sysname}")
        machine = os.uname().machine.lower()
        if machine in ("x86_64", "amd64"):
            arch = "x86_64"
        elif machine in ("aarch64", "arm64"):
            arch = "aarch64"
        else:
            raise ValueError(f"Unsupported CPU architecture: {machine}")
        key = f"{os_part}-{arch}"
    if key not in SUPPORTED_PLATFORMS:
        raise ValueError(
            f"Unsupported platform {key!r}; expected one of "
            + ", ".join(SUPPORTED_PLATFORMS)
        )
    return key


def hash_file(path: Path, *algorithms: str) -> dict[str, str]:
    """Compute all requested digests of *path* in a single streaming pass."""
    hashers = {algo: hashlib.new(algo) for algo in algorithms}
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            for hasher in hashers.values():
                hasher.update(block)
    return {algo: h.hexdigest() for algo, h in hashers.items()}


def verify_file(
    path: Path,
    *,
    size: int | None = None,
    sha256: str | None = None,
    md5: str | None = None,
) -> tuple[bool, str]:
    """Check size and/or checksums. ``(False, reason)`` on any mismatch."""
    if not path.is_file():
        return False, "not a regular file"
    if size is not None and path.stat().st_size != size:
        return False, (
            f"size mismatch (got {path.stat().st_size}, expected {size})"
        )
    wanted = [a for a, v in (("sha256", sha256), ("md5", md5)) if v]
    if wanted:
        digests = hash_file(path, *wanted)
        if sha256 is not None and digests["sha256"] != sha256:
            return False, f"sha256 mismatch ({digests['sha256']})"
        if md5 is not None and digests["md5"] != md5:
            return False, f"md5 mismatch ({digests['md5']})"
    return True, "ok"


def _verify_spec(path: Path, spec: dict) -> tuple[bool, str]:
    return verify_file(
        path,
        size=spec.get("bytes"),
        sha256=spec.get("sha256"),
        md5=spec.get("md5"),
    )


# ---------------------------------------------------------------------------
# Resumable downloader (curl preferred, urllib fallback)
# ---------------------------------------------------------------------------


def _curl_download(url: str, tmp: Path) -> None:
    cmd = [
        "curl", "-fsSL",
        "--retry", "3",
        "--retry-delay", "2",
        "--continue-at", "-",
        "-o", str(tmp),
        url,
    ]
    result = subprocess.run(cmd)
    if result.returncode == 0:
        return
    if tmp.exists() and tmp.stat().st_size > 0:
        # Server may have rejected the range request: restart cleanly once.
        tmp.unlink(missing_ok=True)
        result = subprocess.run(cmd)
    if result.returncode != 0:
        raise FetchError(f"curl exited {result.returncode} for {url}")


def _urllib_download(url: str, tmp: Path) -> None:
    have = tmp.stat().st_size if tmp.exists() else 0
    request = urllib.request.Request(url)
    if have:
        request.add_header("Range", f"bytes={have}-")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            if have and response.status != 206:
                # Server ignored the Range header: start over.
                have = 0
            mode = "ab" if have else "wb"
            with tmp.open(mode) as out:
                shutil.copyfileobj(response, out, 1024 * 1024)
    except urllib.error.URLError as exc:
        raise FetchError(f"download failed for {url}: {exc}") from exc


def download_file(
    url: str,
    dest: Path,
    *,
    expected_bytes: int | None = None,
    sha256: str | None = None,
    md5: str | None = None,
    force: bool = False,
) -> str:
    """Fetch *url* into *dest* with resume support and verify-before-move.

    A ``<dest>.part`` partial survives transport failures so a later run
    resumes where it stopped (``curl --continue-at -`` or a Range request).
    Checksums are verified before the file is moved into place; a verified
    partial failing checksum is removed rather than silently overwriting.

    Returns ``"already_present"``, ``"downloaded"``, or ``"overwritten"``.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        ok, why = verify_file(dest, size=expected_bytes, sha256=sha256, md5=md5)
        if ok and not force:
            return "already_present"
        if not ok and not force:
            raise FetchError(
                f"{dest} exists but fails verification ({why}); "
                "delete it or re-run with --force"
            )

    tmp = dest.with_name(dest.name + ".part")
    try:
        if shutil.which("curl"):
            _curl_download(url, tmp)
        else:
            _urllib_download(url, tmp)
    except FetchError:
        raise
    except Exception as exc:  # transport failure: keep the partial for resume
        raise FetchError(f"download failed for {url}: {exc}") from exc

    ok, why = verify_file(tmp, size=expected_bytes, sha256=sha256, md5=md5)
    if not ok:
        tmp.unlink(missing_ok=True)
        raise FetchError(f"verification failed for {url}: {why}")
    existed = dest.exists()
    os.replace(tmp, dest)
    return "overwritten" if existed else "downloaded"


# ---------------------------------------------------------------------------
# Safe tar extraction (rejects traversal, absolute paths, escaping links)
# ---------------------------------------------------------------------------


def _member_relpath(name: str, strip_components: int) -> str | None:
    norm = posixpath.normpath(name)
    if norm.startswith("/"):
        raise FetchError(f"absolute path in archive: {name!r}")
    parts = [p for p in norm.split("/") if p not in ("", ".")]
    if any(p == ".." for p in parts):
        raise FetchError(f"path traversal in archive: {name!r}")
    parts = parts[strip_components:]
    if not parts:
        return None
    return "/".join(parts)

def _check_within(dest_dir: Path, target: Path) -> None:
    resolved = target.resolve(strict=False)
    if not resolved.is_relative_to(dest_dir):
        raise FetchError(f"archive member escapes destination: {target}")


def safe_extract_tar(
    archive: Path,
    dest_dir: Path,
    *,
    strip_components: int = 0,
    force: bool = False,
) -> dict:
    """Extract a .tar / .tar.gz / .tar.xz defensively into *dest_dir*.

    Rejects ``..`` traversal, absolute member paths, device/FIFO members,
    and any symlink/hardlink whose target escapes the archive root. Existing
    files are never overwritten unless ``force`` is set.
    """
    archive = Path(archive)
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_resolved = dest_dir.resolve()
    with tarfile.open(archive, "r:*") as tar:
        planned: list[tuple[tarfile.TarInfo, str]] = []
        for member in tar.getmembers():
            rel = _member_relpath(member.name, strip_components)
            if rel is None:
                continue
            if member.isdev() or member.isfifo():
                raise FetchError(
                    f"unsupported member type {member.type!r}: {member.name!r}"
                )
            if member.issym() or member.islnk():
                link_base = posixpath.dirname(member.name)
                target_name = posixpath.normpath(
                    posixpath.join(link_base, member.linkname)
                )
                link_rel = _member_relpath(target_name, strip_components)
                if link_rel is None:
                    raise FetchError(
                        f"link target escapes archive root: {member.name!r} "
                        f"-> {member.linkname!r}"
                    )
            planned.append((member, rel))
        conflicts = [rel for _, rel in planned if (dest_dir / rel).exists()]
        if conflicts and not force:
            raise FetchError(
                f"refusing to overwrite existing paths in {dest_dir}: "
                f"{conflicts[:5]}{' ...' if len(conflicts) > 5 else ''} "
                "(re-run with --force to replace)"
            )
        written: list[str] = []
        for member, rel in planned:
            target = dest_dir / rel
            _check_within(dest_resolved, target)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            if member.issym() or member.islnk():
                if target.exists() or target.is_symlink():
                    target.unlink()
                if member.issym():
                    target.symlink_to(member.linkname)
                else:
                    link_rel = _member_relpath(
                        posixpath.normpath(
                            posixpath.join(
                                posixpath.dirname(member.name), member.linkname
                            )
                        ),
                        strip_components,
                    )
                    os.link(dest_dir / link_rel, target)
                written.append(rel)
                continue
            source = tar.extractfile(member)
            if source is None:
                raise FetchError(f"cannot read archive member {member.name!r}")
            with source, target.open("wb") as out:
                shutil.copyfileobj(source, out, 1024 * 1024)
            target.chmod(member.mode & 0o777 or 0o644)
            written.append(rel)
    return {"extracted": written, "files": len(written)}


# ---------------------------------------------------------------------------
# Binary resolution across drop-in dirs, env override, and PATH
# ---------------------------------------------------------------------------


def binary_search_dirs(bin_dir: Path | None = None) -> list[Path]:
    dirs: list[Path] = []
    if bin_dir is not None:
        dirs.append(Path(bin_dir))
    env_dir = os.environ.get("G2F_TOOLS_BIN")
    if env_dir:
        dirs.append(Path(env_dir).expanduser())
    if DEFAULT_BIN_DIR.is_dir():
        dirs.append(DEFAULT_BIN_DIR)
    return dirs


def resolve_binary(name: str, bin_dir: Path | None = None) -> Path | None:
    """Find an executable *name* in drop-in dirs, G2F_TOOLS_BIN, then PATH."""
    for directory in binary_search_dirs(bin_dir):
        candidate = directory / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    found = shutil.which(name)
    return Path(found) if found else None


def run_command(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    """subprocess.run wrapper kept small so tests can mock it."""
    return subprocess.run(cmd, **kwargs)


# ---------------------------------------------------------------------------
# eggNOG fungi.dmnd selection contract (eggnog-mapper 2.1.12 create_dbs.py)
# ---------------------------------------------------------------------------


def select_fungi_taxids(taxid_info_path: Path) -> set[str]:
    """Taxids whose 'Named Lineage' column contains the exact token Fungi."""
    taxids: set[str] = set()
    with open(taxid_info_path, "r", encoding="utf-8") as table:
        for index, line in enumerate(table):
            if index == 0 or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            named_lineage = {token.strip() for token in fields[3].split(",")}
            if "Fungi" in named_lineage:
                taxids.add(fields[0].strip())
    if not taxids:
        raise FetchError(f"no Fungi taxids parsed from {taxid_info_path}")
    return taxids


def filter_proteome_fasta(
    proteomes_path: Path, taxids: set[str], out_path: Path
) -> dict:
    """Keep records whose header token before the first '.' is in *taxids*."""
    hasher = hashlib.sha256()
    selected_records = 0
    selected_lines = discarded_lines = 0
    accept = False
    with open(proteomes_path, "r", encoding="utf-8") as src, open(
        out_path, "w", encoding="utf-8", newline=""
    ) as out:
        for line in src:
            stripped = line.strip()
            if line.startswith(">"):
                taxid = line[1:].split(".")[0]
                accept = taxid in taxids
                if accept:
                    selected_records += 1
            if accept:
                out.write(stripped + "\n")
                hasher.update((stripped + "\n").encode())
                selected_lines += 1
            elif not line.startswith(">"):
                discarded_lines += 1
    return {
        "selected_records": selected_records,
        "selected_lines": selected_lines,
        "discarded_lines": discarded_lines,
        "sha256": hasher.hexdigest(),
    }


# ---------------------------------------------------------------------------
# Python interpreter probes (subprocess-isolated so heavy deps never leak in)
# ---------------------------------------------------------------------------


def probe_python_packages(
    python_bin: str, modules: tuple[str, ...], timeout: int = 120
) -> dict:
    """Run a ``python -c`` import probe; report missing modules and status.

    Module names may contain ``|`` alternatives (``CLEAN|clean``) where any
    one satisfying the spec counts as present.
    """
    expanded = Path(python_bin).expanduser()
    resolved = shutil.which(python_bin) or (
        str(expanded) if expanded.is_file() else None
    )
    result = {
        "interpreter": resolved or python_bin,
        "status": "unknown",
        "missing": [],
        "detail": "",
    }
    if resolved is None:
        result.update(status="error", detail="interpreter not found")
        return result
    try:
        proc = run_command(
            [resolved, "-c", _PROBE_IMPORTS_SNIPPET, *modules],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        result.update(status="error", detail=str(exc))
        return result
    output = ((proc.stdout or "") + " " + (proc.stderr or "")).strip()
    if proc.returncode == 0 and "OK" in (proc.stdout or ""):
        result["status"] = "ok"
    elif "MISSING:" in (proc.stdout or ""):
        missing_line = proc.stdout.split("MISSING:", 1)[1].splitlines()[0]
        result.update(
            status="missing",
            missing=[m for m in missing_line.split(",") if m],
            detail=output[:400],
        )
    else:
        result.update(status="error", detail=output[:400])
    return result


def probe_highs_milp(python_bin: str | None = None, timeout: int = 120) -> dict:
    """Verify the HiGHS LP/MILP backend end-to-end in a fresh interpreter."""
    binary = python_bin or sys.executable
    try:
        proc = run_command(
            [binary, "-c", _PROBE_HIGHS_MILP_SNIPPET],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"status": "error", "detail": str(exc), "interpreter": binary}
    if proc.returncode == 0 and "MILP_OK" in (proc.stdout or ""):
        return {"status": "ok", "interpreter": binary}
    detail = ((proc.stdout or "") + " " + (proc.stderr or "")).strip()
    return {"status": "error", "detail": detail[:400], "interpreter": binary}


# ---------------------------------------------------------------------------
# Profile expansion
# ---------------------------------------------------------------------------


def expand_profiles(profiles: list[str] | tuple[str, ...]) -> list[str]:
    """Expand ``all`` and implications into an ordered unique list.

    ``alphagem`` implies ``gapfill``, which in turn implies ``screening``.
    """
    requested = {p.strip().lower() for p in profiles if p and p.strip()}
    if not requested:
        requested = {"screening"}
    unknown = requested - set(ALL_PROFILES) - {"all"}
    if unknown:
        raise ValueError(
            f"unknown profile(s) {sorted(unknown)}; choose from "
            + ", ".join((*ALL_PROFILES, "all"))
        )
    selected: set[str] = set(ALL_PROFILES) if "all" in requested else set()
    for name in requested & set(ALL_PROFILES):
        selected.add(name)
        if name == "alphagem":
            selected.add("gapfill")
        if name in ("alphagem", "gapfill"):
            selected.add("screening")
    return [name for name in ALL_PROFILES if name in selected]


# ---------------------------------------------------------------------------
# Execution context + step bookkeeping
# ---------------------------------------------------------------------------


@dataclass
class Context:
    args: argparse.Namespace
    platform: str
    bin_dir: Path
    runtime_root: Path | None
    tool_root: Path | None
    weights_dir: Path | None
    dry_run: bool = False
    force: bool = False
    steps: list[dict] = field(default_factory=list)

    def run_step(self, profile: str, action: str, fn, status: str = "planned", **fields) -> dict:
        """Record a step; in dry-run mark it *status*, else execute *fn*."""
        step = {"profile": profile, "action": action, **fields}
        self.steps.append(step)
        if self.dry_run:
            step["status"] = status
        else:
            fn(step)
        return step


def _spec_fields(spec: dict) -> dict:
    out = {"url": spec.get("url")}
    for key in ("bytes", "md5", "sha256"):
        if spec.get(key):
            out[key] = spec[key]
    return out


def _fetch_step(ctx: Context, profile: str, url: str, dest: Path, spec: dict, label: str):
    def do(step: dict) -> None:
        step["status"] = download_file(
            url,
            dest,
            expected_bytes=spec.get("bytes"),
            sha256=spec.get("sha256"),
            md5=spec.get("md5"),
            force=ctx.force,
        )

    return ctx.run_step(
        profile, "download", do, label=label, dest=str(dest), **_spec_fields(spec)
    )


# ---------------------------------------------------------------------------
# Profile: screening (DIAMOND)
# ---------------------------------------------------------------------------


def run_screening(ctx: Context) -> None:
    dest = ctx.bin_dir / "diamond"
    if ctx.platform in DIAMOND_URLS:
        url = DIAMOND_URLS[ctx.platform]
        tarball = ctx.bin_dir / "_downloads" / "diamond-linux64.tar.gz"
        existing = resolve_binary("diamond", ctx.bin_dir)
        if existing is not None and not ctx.force:
            ctx.run_step(
                "screening",
                "install_diamond",
                lambda s: s.update(status="already_installed", binary=str(existing)),
                status="already_installed",
                url=url,
                dest=str(dest),
                note="diamond already on PATH or in a drop-in directory",
            )
            return

        def install(step: dict) -> None:
            ctx.bin_dir.mkdir(parents=True, exist_ok=True)
            download_file(url, tarball, force=ctx.force)
            extracted = False
            with tarfile.open(tarball, "r:*") as tar:
                for member in tar.getmembers():
                    rel = _member_relpath(member.name, 0)
                    if rel is None or not member.isfile():
                        continue
                    if posixpath.basename(rel) != "diamond":
                        continue
                    target = ctx.bin_dir / "diamond"
                    _check_within(ctx.bin_dir.resolve(), target)
                    src = tar.extractfile(member)
                    assert src is not None
                    with src, target.open("wb") as out:
                        shutil.copyfileobj(src, out)
                    target.chmod(0o755)
                    extracted = True
            if not extracted:
                raise FetchError("diamond executable not found inside tarball")
            step["status"] = "installed"
            step["binary"] = str(dest)

        ctx.run_step(
            "screening",
            "install_diamond",
            install,
            url=url,
            tarball=str(tarball),
            dest=str(dest),
            note="upstream publishes no checksum; binary verified by extraction",
        )
        return

    # No official binary release for this platform.
    existing = resolve_binary("diamond", ctx.bin_dir)
    if existing is not None:
        ctx.run_step(
            "screening",
            "install_diamond",
            lambda s: s.update(status="already_installed", binary=str(existing)),
            status="already_installed",
            note="diamond already discoverable",
        )
        return
    ctx.run_step(
        "screening",
        "install_diamond",
        lambda s: s.update(status="manual_action_required"),
        status="manual_action_required",
        note=DIAMOND_MANUAL_HINT,
        platform=ctx.platform,
    )


# ---------------------------------------------------------------------------
# Profile: gapfill (NCBI BLAST+)
# ---------------------------------------------------------------------------


def run_gapfill(ctx: Context) -> None:
    spec = BLAST_PLUS["archives"][ctx.platform]
    url = f"{BLAST_PLUS['base']}/{spec['file']}"
    tarball = ctx.bin_dir / "_downloads" / spec["file"]

    present = {b: resolve_binary(b, ctx.bin_dir) for b in BLAST_PLUS["binaries"]}
    if all(present.values()) and not ctx.force:
        ctx.run_step(
            "gapfill",
            "install_blast_plus",
            lambda s: s.update(
                status="already_installed",
                binaries={k: str(v) for k, v in present.items()},
            ),
            status="already_installed",
            note="tblastn/makeblastdb/blastp already discoverable",
        )
        return

    def install(step: dict) -> None:
        ctx.bin_dir.mkdir(parents=True, exist_ok=True)
        download_file(url, tarball, md5=spec["md5"], force=ctx.force)
        installed: list[str] = []
        with tarfile.open(tarball, "r:*") as tar:
            for member in tar.getmembers():
                rel = _member_relpath(member.name, 0)
                if rel is None or not member.isfile():
                    continue
                parts = rel.split("/")
                if len(parts) == 3 and parts[1] == "bin":
                    target = ctx.bin_dir / parts[2]
                    _check_within(ctx.bin_dir.resolve(), target)
                    src = tar.extractfile(member)
                    assert src is not None
                    with src, target.open("wb") as out:
                        shutil.copyfileobj(src, out)
                    target.chmod(0o755)
                    installed.append(parts[2])
        missing = [b for b in BLAST_PLUS["binaries"] if b not in installed]
        if missing:
            raise FetchError(f"BLAST+ tarball lacked expected binaries: {missing}")
        step["status"] = "installed"
        step["binaries"] = sorted(installed)

    ctx.run_step(
        "gapfill",
        "install_blast_plus",
        install,
        url=url,
        tarball=str(tarball),
        dest=str(ctx.bin_dir),
        md5=spec["md5"],
        note="extracts */bin/* executables into --bin-dir",
    )


# ---------------------------------------------------------------------------
# Profile: alphagem (Figshare data source + ESM weights + optional extras)
# ---------------------------------------------------------------------------


def _prepare_runtime_module():
    try:
        from scripts import prepare_alphagem_runtime as mod  # type: ignore
    except ImportError:
        import prepare_alphagem_runtime as mod  # type: ignore
    return mod


def _alphagem_config_snippet(ctx: Context) -> dict:
    return {
        "alphagem_runtime_root": (
            str(ctx.runtime_root)
            if ctx.runtime_root
            else "/path/to/alphagem-runtime"
        ),
        "reference_fasta": (
            str(ctx.runtime_root / "data_available" / "yeast.fasta")
            if ctx.runtime_root
            else "/path/to/alphagem-runtime/data_available/yeast.fasta"
        ),
        "reference_xlsx": (
            str(ctx.runtime_root / "data_available" / "yeast.xlsx")
            if ctx.runtime_root
            else "/path/to/alphagem-runtime/data_available/yeast.xlsx"
        ),
        "alphagem_python": "<separate AlphaGEM conda env python (AlphaGEM.yml)>",
        "alphagem_mode": "plmsearch",
        "deepec_python": "<deepecenvironment.yml env python>",
        "esm_model_path": (
            str(ctx.weights_dir / "esm" / "esm1b_t33_650M_UR50S.pt")
            if ctx.weights_dir
            else "<weights-dir>/esm/esm1b_t33_650M_UR50S.pt"
        ),
        "prottrans_model_path": (
            (
                str(ctx.weights_dir / "prot_t5_xl_uniref50")
                if ctx.weights_dir
                else "<weights-dir>/prot_t5_xl_uniref50"
            )
            if bool(getattr(ctx.args, "with_prott5", False))
            else None
        ),
        "structure_maplist": None,
        "structure_dir": None,
    }


def _alphagem_env_notice() -> str:
    return (
        "AlphaGEM steps 1-3 need a SEPARATE Python environment: create it "
        "from third_party/alphagem/AlphaGEM.yml (conda env create -f "
        "AlphaGEM.yml) plus third_party/alphagem/deepecenvironment.yml for "
        'DeepEC, then set "alphagem_python"/G2F_ALPHAGEM_PYTHON and '
        '"deepec_python" to those interpreters. The uv-managed Geno2Flux '
        "env intentionally lacks torch/esm. Verify with `bash run.sh "
        "check-env --runtime-root ... --alphagem-python ...`."
    )


def run_alphagem(ctx: Context) -> None:
    if ctx.runtime_root is None:
        if ctx.dry_run:
            ctx.run_step(
                "alphagem",
                "plan_runtime",
                lambda s: s.update(status="needs_runtime_root"),
                status="needs_runtime_root",
                note="pass --runtime-root (outside the repo) or --dest to see the full plan",
            )
            return
        raise FetchError(
            "--profile alphagem requires --runtime-root (or --dest) outside the repository"
        )
    runtime = ctx.runtime_root
    weights = ctx.weights_dir or (runtime / "var" / "weights")

    def stage_source(step: dict) -> None:
        mod = _prepare_runtime_module()
        mod.require_real_runtime_path(runtime)
        _manifest, records = mod.source_records()
        result = mod.install_runtime(runtime, records)
        step["status"] = "installed" if result["installed"] else "already_present"
        step["installed_files"] = len(result["installed"])
        step["runtime_source_commit"] = result.get("runtime_source_commit")

    ctx.run_step(
        "alphagem",
        "stage_alphagem_source",
        stage_source,
        dest=str(runtime),
        note="copies bundled third_party/alphagem snapshot via prepare_alphagem_runtime",
    )

    archive = runtime / "var" / "downloads" / "data_source.tar.xz"
    extract_marker = runtime / "var" / "downloads" / "data_source.extracted.json"
    ds = ALPHAGEM_DATA_SOURCE
    _fetch_step(ctx, "alphagem", ds["url"], archive, ds, "data_source.tar.xz")

    def extract_data_source(step: dict) -> None:
        ok, why = _verify_spec(archive, ds)
        if not ok:
            raise FetchError(f"{archive} fails verification ({why})")
        expected_dirs = ("data_available", "rhea", "plmsearchtools", "tools")
        if (
            not ctx.force
            and extract_marker.is_file()
            and all((runtime / d).is_dir() for d in expected_dirs)
        ):
            try:
                marker_data = json.loads(extract_marker.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                marker_data = {}
            if marker_data.get("sha256") == ds.get("sha256"):
                step["status"] = "already_extracted"
                step["files"] = marker_data.get("files")
                return
        info = safe_extract_tar(archive, runtime, force=ctx.force)
        extract_marker.parent.mkdir(parents=True, exist_ok=True)
        extract_marker.write_text(
            json.dumps({"sha256": ds.get("sha256"), "files": info["files"]}, indent=2) + "\n",
            encoding="utf-8",
        )
        step["status"] = "extracted"
        step["files"] = info["files"]
    ctx.run_step(
        "alphagem",
        "extract_data_source",
        extract_data_source,
        dest=str(runtime),
        note="places data_available/, rhea/, plmsearchtools/, tools/ under runtime root",
    )

    esm_dir = weights / "esm"
    for name, spec in ESM_WEIGHTS.items():
        _fetch_step(ctx, "alphagem", spec["url"], esm_dir / name, spec, name)

    if ctx.args.with_eggnog:
        _run_eggnog(ctx, runtime)
    else:
        ctx.run_step(
            "alphagem",
            "eggnog",
            lambda s: s.update(status="skipped"),
            status="skipped",
            note="required for fresh AlphaGEM runs (~41 GiB): re-run with --with-eggnog unless eggnog.db + fungi.dmnd are already staged",
        )

    if ctx.args.with_prott5:
        prott5_dir = weights / "prot_t5_xl_uniref50"

        def prott5(step: dict) -> None:
            hf_cli = shutil.which("huggingface-cli") or shutil.which("hf")
            cmd = [
                hf_cli or "huggingface-cli",
                "download",
                PROTT5_REPO,
                "--local-dir",
                str(prott5_dir),
            ]
            step["command"] = " ".join(cmd)
            if hf_cli is None:
                step["status"] = "manual_action_required"
                return
            prott5_dir.mkdir(parents=True, exist_ok=True)
            proc = run_command(cmd)
            if proc.returncode != 0:
                raise FetchError(f"huggingface-cli download exited {proc.returncode}")
            step["status"] = "downloaded"
            step["dest"] = str(prott5_dir)

        ctx.run_step(
            "alphagem",
            "download_prott5",
            prott5,
            url=PROTT5_URL,
            dest=str(prott5_dir),
            note="only needed for alphagem_mode=structure_alignment",
        )

    ctx.run_step(
        "alphagem",
        "config_snippet",
        lambda s: s.update(status="info"),
        status="info",
        config=_alphagem_config_snippet(ctx),
        notice=_alphagem_env_notice(),
    )


def _run_eggnog(ctx: Context, runtime: Path) -> None:
    data_dir = runtime / "tools" / "eggnog_mapper" / "data"
    cache_dir = runtime / "var" / "downloads" / "eggnog"
    gz_path = cache_dir / "eggnog.db.gz"
    db_path = data_dir / "eggnog.db"

    def fetch_db(step: dict) -> None:
        step["status"] = download_file(EGGNOG["db_gz_url"], gz_path, force=ctx.force)

    ctx.run_step(
        "alphagem",
        "download_eggnog_db_gz",
        fetch_db,
        url=EGGNOG["db_gz_url"],
        dest=str(gz_path),
        note="no published hash for the .gz; the decompressed db sha256 is verified next",
    )

    def gunzip_db(step: dict) -> None:
        ok, _why = verify_file(
            db_path, size=EGGNOG["db_bytes"], sha256=EGGNOG["db_sha256"]
        )
        if ok and not ctx.force:
            step["status"] = "already_present"
            return
        data_dir.mkdir(parents=True, exist_ok=True)
        tmp = db_path.with_name(".eggnog.db.partial")
        try:
            with gzip.open(gz_path, "rb") as src, tmp.open("wb") as out:
                shutil.copyfileobj(src, out, 4 * 1024 * 1024)
            ok, why = verify_file(
                tmp, size=EGGNOG["db_bytes"], sha256=EGGNOG["db_sha256"]
            )
            if not ok:
                tmp.unlink(missing_ok=True)
                raise FetchError(f"decompressed eggnog.db failed verification: {why}")
            os.replace(tmp, db_path)
        finally:
            tmp.unlink(missing_ok=True)
        step["status"] = "verified"
        step["sha256"] = EGGNOG["db_sha256"]

    ctx.run_step(
        "alphagem",
        "decompress_eggnog_db",
        gunzip_db,
        dest=str(db_path),
        sha256=EGGNOG["db_sha256"],
        bytes=EGGNOG["db_bytes"],
    )

    proteomes = cache_dir / "e5.proteomes.faa"
    taxid_info = cache_dir / "e5.taxid_info.tsv"
    _fetch_step(
        ctx,
        "alphagem",
        EGGNOG["proteomes_url"],
        proteomes,
        {"url": EGGNOG["proteomes_url"], "sha256": EGGNOG["proteomes_sha256"]},
        "e5.proteomes.faa",
    )
    _fetch_step(
        ctx,
        "alphagem",
        EGGNOG["taxid_url"],
        taxid_info,
        {"url": EGGNOG["taxid_url"], "sha256": EGGNOG["taxid_sha256"]},
        "e5.taxid_info.tsv",
    )

    fungi_fasta = cache_dir / "e5.proteomes.fungi.faa"

    def build_dmnd(step: dict) -> None:
        dmnd = data_dir / "fungi.dmnd"
        if dmnd.is_file() and not ctx.force:
            ok, _why = verify_file(
                dmnd,
                size=EGGNOG["fungi_dmnd_bytes"],
                sha256=EGGNOG["fungi_dmnd_sha256"],
            )
            if ok:
                step["status"] = "already_present"
                step["dmnd_bytes"] = EGGNOG["fungi_dmnd_bytes"]
                step["dmnd_sha256"] = EGGNOG["fungi_dmnd_sha256"]
                return
        cache_dir.mkdir(parents=True, exist_ok=True)
        data_dir.mkdir(parents=True, exist_ok=True)
        taxids = select_fungi_taxids(taxid_info)
        step["fungi_taxids"] = len(taxids)
        stats = filter_proteome_fasta(proteomes, taxids, fungi_fasta)
        step.update(
            selected_records=stats["selected_records"],
            selected_fasta_sha256=stats["sha256"],
        )
        if stats["sha256"] != EGGNOG["selected_fasta_sha256"]:
            raise FetchError(
                "fungi FASTA selection sha256 "
                f"{stats['sha256']} != contract "
                f"{EGGNOG['selected_fasta_sha256']} — refusing to build "
                "fungi.dmnd from a selection that deviates from the verified "
                "contract"
            )
        diamond = resolve_binary("diamond", ctx.bin_dir)
        if diamond is None:
            raise FetchError(
                "diamond not found for `diamond makedb`; run --profile screening first"
            )
        proc = run_command(
            [
                str(diamond), "makedb",
                "--in", str(fungi_fasta),
                "-d", str(dmnd.with_suffix("")),
            ]
        )
        if proc.returncode != 0 or not dmnd.is_file():
            raise FetchError(f"diamond makedb exited {proc.returncode}")
        actual_bytes = dmnd.stat().st_size
        digest = hash_file(dmnd, "sha256")["sha256"]
        step["dmnd_bytes"] = actual_bytes
        step["dmnd_sha256"] = digest
        if (
            actual_bytes != EGGNOG["fungi_dmnd_bytes"]
            or digest != EGGNOG["fungi_dmnd_sha256"]
        ):
            step["status"] = "contract_mismatch"
            step["warning"] = (
                "fungi.dmnd does not match the verified paper contract: got "
                f"{actual_bytes} bytes / sha256 {digest}, expected "
                f"{EGGNOG['fungi_dmnd_bytes']} bytes / sha256 "
                f"{EGGNOG['fungi_dmnd_sha256']}. The unified pipeline invokes "
                "AlphaGEM.py --profile paper, which hard-stops on this "
                "mismatch; rebuilding needs Linux x86_64 DIAMOND v2.1.10 or "
                "a prebuilt contract-matching fungi.dmnd."
            )
        else:
            step["status"] = "built"

    ctx.run_step(
        "alphagem",
        "build_fungi_dmnd",
        build_dmnd,
        dest=str(data_dir / "fungi.dmnd"),
        contract_bytes=EGGNOG["fungi_dmnd_bytes"],
        contract_sha256=EGGNOG["fungi_dmnd_sha256"],
        selection_contract=EGGNOG["selection_contract"],
    )


# ---------------------------------------------------------------------------
# Profile: annotation (ANNEVO + optional Helixer)
# ---------------------------------------------------------------------------


def run_annotation(ctx: Context) -> None:
    if ctx.tool_root is None:
        if ctx.dry_run:
            ctx.run_step(
                "annotation",
                "plan_tool_root",
                lambda s: s.update(status="needs_tool_root"),
                status="needs_tool_root",
                note="pass --tool-root or --dest to see the full plan",
            )
            return
        raise FetchError("--profile annotation requires --tool-root (or --dest)")
    tool_root = ctx.tool_root
    annevo_dir = tool_root / "tools" / "annevo" / "annevo-current"
    archive = tool_root / "_downloads" / "annevo-main.tar.gz"

    def fetch_archive(step: dict) -> None:
        step["status"] = download_file(
            ANNEVO["archive_url"], archive, force=ctx.force
        )

    ctx.run_step(
        "annotation",
        "download_annevo_source",
        fetch_archive,
        url=ANNEVO["archive_url"],
        dest=str(archive),
        note="GitHub main-branch archive; no published checksum",
    )

    def extract_archive(step: dict) -> None:
        if (
            not ctx.force
            and (annevo_dir / "prediction.py").is_file()
            and (annevo_dir / "decoding.py").is_file()
        ):
            step["status"] = "already_extracted"
            return
        info = safe_extract_tar(
            archive, annevo_dir, strip_components=1, force=ctx.force
        )
        step["status"] = "extracted"
        step["files"] = info["files"]
    ctx.run_step(
        "annotation",
        "extract_annevo_source",
        extract_archive,
        dest=str(annevo_dir),
        note="strips top-level ANNEVO-main/ so prediction.py sits at the root",
    )

    weights_dest = tool_root / "weights" / "annevo" / "ANNEVO_Fungi.pt"
    _fetch_step(
        ctx,
        "annotation",
        ANNEVO["weights_url"],
        weights_dest,
        {"url": ANNEVO["weights_url"], "bytes": ANNEVO["weights_bytes"]},
        "ANNEVO_Fungi.pt",
    )

    if ctx.args.with_helixer:
        _fetch_step(
            ctx,
            "annotation",
            HELIXER_FUNGI_URL,
            tool_root / "weights" / "helixer" / "fungi_v0.3_a_0100.h5",
            {"url": HELIXER_FUNGI_URL},
            "fungi_v0.3_a_0100.h5",
        )

    ctx.run_step(
        "annotation",
        "env_notice",
        lambda s: s.update(status="info"),
        status="info",
        note=(
            "ANNEVO runs in its own PyTorch environment. Create one at "
            f"{tool_root}/envs/annevo (python -m venv + pip install torch "
            "h5py biopython + ANNEVO requirements), or pass --tool-python / "
            "--annevo-python. `predict_proteome.py --tool gff` needs no "
            "external predictor."
        ),
    )


# ---------------------------------------------------------------------------
# Profile: deeploc (DTU license-gated; local tarball install only)
# ---------------------------------------------------------------------------


def _resolve_deeploc() -> tuple[list[str] | None, str]:
    """Delegate to geno2flux.localization; inline fallback if deps missing."""
    try:
        from geno2flux.localization import resolve_deeploc_command

        return resolve_deeploc_command()
    except FileNotFoundError as exc:
        return None, f"misconfigured: {exc}"
    except ImportError:
        pass
    # Fallback mirroring resolve_deeploc_command's documented search order.
    env_val = os.environ.get("G2F_DEEPLOC_BIN", "").strip()
    if env_val:
        tokens = shlex.split(env_val)
        if tokens:
            head = Path(tokens[0]).expanduser()
            if head.is_file() and os.access(head, os.X_OK):
                return [str(head.resolve()), *tokens[1:]], "env_G2F_DEEPLOC_BIN"
            found = shutil.which(tokens[0])
            if found:
                return [found, *tokens[1:]], "env_G2F_DEEPLOC_BIN"
            return None, "G2F_DEEPLOC_BIN points to a missing binary"
    for slot in (DEEPLOC_VENV / "bin" / "deeploc2", DEEPLOC_VENV / "deeploc2"):
        if slot.is_file() and os.access(slot, os.X_OK):
            return [str(slot)], "repo_tools_deeploc_slot"
    found = shutil.which("deeploc2")
    if found:
        return [found], "system_PATH"
    return None, "not_found"


def run_deeploc(ctx: Context) -> None:
    tarball = getattr(ctx.args, "deeploc_tarball", None)
    if not tarball:
        ctx.run_step(
            "deeploc",
            "license_request",
            lambda s: s.update(status="manual_action_required"),
            status="manual_action_required",
            note=(
                "DeepLoc is license-gated by DTU Health Tech — it is NEVER "
                "downloaded automatically. Request the academic package, then "
                "re-run with --deeploc-tarball /path/to/DeepLoc-2.x.tar.gz"
            ),
            links={
                "deeploc_2_1_page": DEEPLOC["2.1"]["page"],
                "deeploc_2_1_request": DEEPLOC["2.1"]["form"],
                "deeploc_2_0_page": DEEPLOC["2.0"]["page"],
                "deeploc_2_0_request": DEEPLOC["2.0"]["form"],
            },
        )
        return
    tarball = Path(tarball).expanduser()
    venv = DEEPLOC_VENV
    deeploc_bin = venv / "bin" / "deeploc2"

    def install(step: dict) -> None:
        if not tarball.is_file():
            raise FetchError(f"--deeploc-tarball not found: {tarball}")
        pip = venv / "bin" / "pip"
        if not pip.exists():
            proc = run_command(
                [sys.executable, "-m", "venv", str(venv)],
                capture_output=True,
                text=True,
            )
            if proc.returncode != 0:
                raise FetchError(
                    f"venv creation failed: {proc.stderr.strip()[:300]}"
                )
        proc = run_command(
            [str(pip), "install", str(tarball)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise FetchError(f"pip install failed: {proc.stderr.strip()[:400]}")
        if not (deeploc_bin.is_file() and os.access(deeploc_bin, os.X_OK)):
            raise FetchError(f"{deeploc_bin} missing/non-executable after install")
        cmd, source = _resolve_deeploc()
        if cmd is None:
            raise FetchError(
                "deeploc2 installed but resolve_deeploc_command() cannot detect it"
            )
        step["status"] = "installed"
        step["command"] = cmd
        step["detected_via"] = source

    ctx.run_step(
        "deeploc",
        "install_deeploc_tarball",
        install,
        tarball=str(tarball),
        venv=str(venv),
        note="python -m venv tools/deeploc && tools/deeploc/bin/pip install <tarball>",
    )


# ---------------------------------------------------------------------------
# --check: real readiness audit
# ---------------------------------------------------------------------------


def _check_core(report: dict) -> None:
    probe = probe_python_packages(sys.executable, _CORE_IMPORTS)
    milp = probe_highs_milp(sys.executable)
    ready = probe["status"] == "ok" and milp["status"] == "ok"
    report.update(
        interpreter=probe["interpreter"],
        missing_python_packages=list(probe["missing"]),
        highs_milp=milp,
        status="READY" if ready else "BLOCKED",
    )
    if not ready:
        report["remediation"] = (
            "run `bash run.sh setup` (uv sync --frozen --extra test) so the "
            "core env gains cobra/highspy/scipy/biopython/pandas/openpyxl"
        )


def _check_screening(report: dict, ctx: Context) -> None:
    binary = resolve_binary("diamond", ctx.bin_dir)
    report["searched"] = [str(d) for d in binary_search_dirs(ctx.bin_dir)] + ["PATH"]
    if binary is None:
        if ctx.platform in DIAMOND_URLS:
            remediation = "run `bash run.sh fetch-resources --profile screening`"
        else:
            remediation = DIAMOND_MANUAL_HINT
        report.update(
            status="BLOCKED",
            missing_binaries=["diamond"],
            remediation=remediation,
        )
        return
    proc = run_command([str(binary), "version"], capture_output=True, text=True)
    report.update(
        binary=str(binary),
        version_output=(proc.stdout or "").strip().splitlines()[:1],
    )
    if proc.returncode != 0:
        report.update(
            status="BLOCKED",
            detail="`diamond version` failed",
            remediation="re-install via fetch-resources --profile screening --force",
        )
    else:
        report["status"] = "READY"


def _check_gapfill(report: dict, ctx: Context) -> None:
    binaries: dict[str, str | None] = {}
    missing: list[str] = []
    for name in ("tblastn", "makeblastdb"):
        found = resolve_binary(name, ctx.bin_dir)
        binaries[name] = str(found) if found else None
        if found is None:
            missing.append(name)
    report.update(
        binaries=binaries,
        missing_binaries=missing,
        searched=[str(d) for d in binary_search_dirs(ctx.bin_dir)] + ["PATH"],
    )
    if missing:
        report.update(
            status="BLOCKED",
            remediation="run `bash run.sh fetch-resources --profile gapfill`",
        )
        return
    report["version_output"] = {}
    for name, path in binaries.items():
        proc = run_command([path, "-version"], capture_output=True, text=True)
        if proc.returncode != 0:
            report.update(
                status="BLOCKED",
                detail=f"`{name} -version` failed",
                remediation="re-install via fetch-resources --profile gapfill --force",
            )
            return
        report["version_output"][name] = (proc.stdout or "").strip().splitlines()[:1]
    milp = probe_highs_milp(sys.executable)
    report["highs_milp"] = milp
    if milp["status"] != "ok":
        report.update(
            status="BLOCKED",
            detail="HiGHS MILP backend failed",
            remediation="run `bash run.sh setup` to install highspy",
        )
        return
    report["status"] = "READY"


def _check_deeploc(report: dict) -> None:
    cmd, source = _resolve_deeploc()
    report.update(
        detection=source,
        links={
            "deeploc_2_1_request": DEEPLOC["2.1"]["form"],
            "deeploc_2_0_request": DEEPLOC["2.0"]["form"],
        },
    )
    if cmd is not None:
        report.update(status="READY", command=cmd)
        return
    report.update(
        status="FALLBACK_AVAILABLE",
        command=None,
        note=(
            "deeploc2 not detected (checked config, G2F_DEEPLOC_BIN, "
            "tools/deeploc/, PATH). The pipeline still runs with default "
            "compartments unless require_deeploc=true. To enable, request "
            "the DTU package and run `bash run.sh fetch-resources --profile "
            "deeploc --deeploc-tarball /path/to/DeepLoc-2.x.tar.gz`."
        ),
    )


def _verify_eggnog_static_data(runtime: Path) -> dict | None:
    """Verify the 9-file eggNOG static data bundle against EGGNOG_STATIC_DATA.

    Mirrors the aggregate-identity algorithm in
    third_party/alphagem/src/pipeline_runner.py (``_paper_bundle_record``):
    sorted ``<file_sha256>  <repo-relative-posix-path>`` lines hashed with
    SHA-256, plus file_count and total bytes. Returns a contract-mismatch
    dict on any deviation, ``None`` when the bundle verifies.
    """
    spec = EGGNOG_STATIC_DATA
    bundle_rel = spec["path"]
    root = runtime / bundle_rel

    def mismatch(reason: str, observed: dict | None = None) -> dict:
        return {
            "path": bundle_rel,
            "reason": reason,
            "expected_file_count": spec["file_count"],
            "expected_bytes": spec["bytes"],
            "expected_sha256": spec["sha256"],
            "observed_file_count": (observed or {}).get("file_count"),
            "observed_bytes": (observed or {}).get("bytes"),
            "observed_sha256": (observed or {}).get("sha256"),
        }

    if not root.is_dir():
        return mismatch(f"eggnog static data directory is missing: {root}")

    excluded = tuple(
        Path(relative).parts for relative in spec["excluded_relative_roots"]
    )

    def is_excluded(relative: Path) -> bool:
        return any(
            relative.parts[: len(prefix)] == prefix for prefix in excluded
        )

    entries: list[tuple[str, str, int]] = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if "__pycache__" in relative.parts or relative.suffix == ".pyc":
            return mismatch(
                f"eggnog static data contains forbidden Python bytecode: {path}"
            )
        if is_excluded(relative):
            continue
        if path.is_symlink():
            return mismatch(
                f"eggnog static data contains a symlink: {path}"
            )
        if path.is_dir():
            continue
        if not path.is_file():
            return mismatch(
                f"eggnog static data contains a non-regular entry: {path}"
            )
        entries.append(
            (
                path.relative_to(runtime).as_posix(),
                hash_file(path, "sha256")["sha256"],
                path.stat().st_size,
            )
        )
    if not entries:
        return mismatch("eggnog static data contains no regular files")
    entries.sort()
    aggregate = hashlib.sha256()
    for relative, file_sha256, _size in entries:
        aggregate.update(f"{file_sha256}  {relative}\n".encode("utf-8"))
    observed = {
        "file_count": len(entries),
        "bytes": sum(size for _rel, _sha, size in entries),
        "sha256": aggregate.hexdigest(),
    }
    expected = {
        key: spec[key] for key in ("file_count", "bytes", "sha256")
    }
    if observed != expected:
        return mismatch(
            "eggnog_static_data bundle deviates from the verified paper "
            f"contract: {json.dumps({'observed': observed})}",
            observed,
        )
    return None


def _check_alphagem(report: dict, ctx: Context) -> None:
    runtime = ctx.runtime_root
    missing_files: list[dict] = []
    if runtime is None:
        missing_files = [
            {"path": rel, "hint": hint} for rel, hint in _ALPHAGEM_REQUIRED_FILES
        ]
        report["runtime_files_note"] = (
            "pass --runtime-root to audit runtime + ESM weight files"
        )
    else:
        for rel, hint in _ALPHAGEM_REQUIRED_FILES:
            if not (runtime / rel).exists():
                missing_files.append({"path": rel, "hint": hint})
        esm_arg = getattr(ctx.args, "esm_model_path", None)
        if esm_arg:
            model = Path(esm_arg).expanduser()
            regression = model.with_name(
                model.name.replace(".pt", "-contact-regression.pt")
            )
            candidates = [model, regression]
        else:
            base = ctx.weights_dir or (runtime / "var" / "weights")
            candidates = [base / "esm" / name for name in ESM_WEIGHTS]
        for path in candidates:
            if not path.is_file():
                missing_files.append(
                    {"path": str(path), "hint": "ESM-1b weights download"}
                )

    contract_mismatches: list[dict] = []
    if runtime is not None:
        dmnd_path = runtime / "tools" / "eggnog_mapper" / "data" / "fungi.dmnd"
        if dmnd_path.is_file():
            ok, why = verify_file(
                dmnd_path,
                size=EGGNOG["fungi_dmnd_bytes"],
                sha256=EGGNOG["fungi_dmnd_sha256"],
            )
            if not ok:
                contract_mismatches.append(
                    {
                        "path": "tools/eggnog_mapper/data/fungi.dmnd",
                        "reason": why,
                        "expected_bytes": EGGNOG["fungi_dmnd_bytes"],
                        "expected_sha256": EGGNOG["fungi_dmnd_sha256"],
                    }
                )
        db_path = runtime / "tools" / "eggnog_mapper" / "data" / "eggnog.db"
        if db_path.is_file():
            ok, why = verify_file(
                db_path,
                size=EGGNOG["db_bytes"],
                sha256=EGGNOG["db_sha256"],
            )
            if not ok:
                contract_mismatches.append(
                    {
                        "path": "tools/eggnog_mapper/data/eggnog.db",
                        "reason": why,
                        "expected_bytes": EGGNOG["db_bytes"],
                        "expected_sha256": EGGNOG["db_sha256"],
                    }
                )
        for stray_name in (
            "e5.proteomes.faa",
            "e5.taxid_info.tsv",
            "e5.proteomes.fungi.faa",
        ):
            stray_path = runtime / "tools" / "eggnog_mapper" / "data" / stray_name
            if stray_path.exists():
                contract_mismatches.append(
                    {
                        "path": f"tools/eggnog_mapper/data/{stray_name}",
                        "reason": (
                            "unexpected intermediate file inside eggnog_static_data "
                            "(AlphaGEM --profile paper rejects extra files in "
                            "tools/eggnog_mapper/data; move or delete it)"
                        ),
                    }
                )
        static_data_mismatch = _verify_eggnog_static_data(runtime)
        if static_data_mismatch is not None:
            contract_mismatches.append(static_data_mismatch)


    missing_binaries = [
        b for b in ("blastp", "diamond") if resolve_binary(b, ctx.bin_dir) is None
    ]

    interpreter = (
        getattr(ctx.args, "alphagem_python", None)
        or os.environ.get("G2F_ALPHAGEM_PYTHON")
        or sys.executable
    )
    alpha_probe = probe_python_packages(str(interpreter), _ALPHAGEM_IMPORTS)

    deepec = getattr(ctx.args, "deepec_python", None)
    deepec_probe = (
        probe_python_packages(str(deepec), _DEEPEC_IMPORTS) if deepec else None
    )
    deepec_warning = None if deepec else (
        "deepec_python not provided; DeepEC-based EC annotation needs it. "
        "Set --deepec-python or the deepec_python config key to audit it."
    )

    dirty_git_entries: list[str] = []
    if runtime is not None and (runtime / ".git").exists():
        git_proc = run_command(
            ["git", "-C", str(runtime), "status", "--porcelain"],
            capture_output=True,
            text=True,
        )
        if git_proc.returncode == 0:
            dirty_git_entries = [
                line for line in (git_proc.stdout or "").splitlines() if line.strip()
            ]

    ready = (
        not missing_files
        and not missing_binaries
        and not dirty_git_entries
        and not contract_mismatches
        and alpha_probe["status"] == "ok"
        and (deepec_probe is None or deepec_probe["status"] == "ok")
    )
    report.update(
        runtime_root=str(runtime) if runtime else None,
        weights_dir=str(ctx.weights_dir) if ctx.weights_dir else None,
        missing_files=missing_files,
        missing_binaries=missing_binaries,
        contract_mismatches=contract_mismatches,
        dirty_git_entries=dirty_git_entries,
        missing_python_packages=list(alpha_probe["missing"]),
        interpreter_used=alpha_probe["interpreter"],
        alphagem_probe=alpha_probe,
        deepec_probe=deepec_probe,
    )
    if deepec_warning:
        report["deepec_warning"] = deepec_warning
    if ready:
        report["status"] = "READY"
    else:
        report["status"] = "BLOCKED"
        report["remediation"] = [
            "data: bash run.sh fetch-resources --profile alphagem "
            "--runtime-root <dir> --with-eggnog",
            "interpreter: conda env create -f third_party/alphagem/AlphaGEM.yml "
            'then set "alphagem_python" or G2F_ALPHAGEM_PYTHON',
            "binaries: bash run.sh fetch-resources --profile gapfill",
        ]
        if contract_mismatches:
            report["remediation"].append(
                "contract: fungi.dmnd/eggnog.db/eggnog_static_data deviate "
                "from the verified paper contract; rebuild with Linux x86_64 "
                "DIAMOND v2.1.10 via "
                "fetch-resources --profile alphagem --runtime-root <dir> "
                "--with-eggnog --force"
            )


def _check_annotation(report: dict, ctx: Context) -> None:
    tool_root = ctx.tool_root
    required = (
        ("tools/annevo/annevo-current/prediction.py", "ANNEVO source archive"),
        ("weights/annevo/ANNEVO_Fungi.pt", "ANNEVO fungal weights"),
    )
    missing_files: list[dict] = []
    annevo_py = getattr(ctx.args, "annevo_python", None)
    if tool_root is None:
        missing_files = [{"path": rel, "hint": hint} for rel, hint in required]
        report["tool_root_note"] = "pass --tool-root to audit annotation files"
    else:
        for rel, hint in required:
            if not (tool_root / rel).exists():
                missing_files.append({"path": rel, "hint": hint})
        default_env = tool_root / "envs" / "annevo" / "bin" / "python"
        if annevo_py is None and default_env.is_file():
            annevo_py = str(default_env)
    probe = (
        probe_python_packages(str(annevo_py), _ANNEVO_IMPORTS)
        if annevo_py
        else {"status": "not_configured", "interpreter": None, "missing": []}
    )
    ready = not missing_files and probe["status"] == "ok"
    report.update(
        tool_root=str(tool_root) if tool_root else None,
        missing_files=missing_files,
        annevo_probe=probe,
        gff_passthrough_ready=True,
        gff_note="`predict_proteome.py --tool gff` needs no external predictors",
        status="READY" if ready else "BLOCKED",
    )
    if not ready:
        report["remediation"] = (
            "run `bash run.sh fetch-resources --profile annotation --tool-root "
            "<dir>` and create a PyTorch env at <tool-root>/envs/annevo (or "
            "pass --annevo-python); --tool gff remains available meanwhile"
        )


def run_preflight_checks(ctx: Context, profiles: list[str]) -> dict:
    """Audit core + requested profiles; every check probes real artifacts."""
    checks: dict[str, dict] = {name: {"status": "SKIPPED"} for name in CHECK_ORDER}
    _check_core(checks["core"])
    runners = {
        "screening": _check_screening,
        "gapfill": _check_gapfill,
        "deeploc": lambda r, _c: _check_deeploc(r),
        "alphagem": _check_alphagem,
        "annotation": _check_annotation,
    }
    for name in profiles:
        runners[name](checks[name], ctx)
    ok_states = {"READY", "FALLBACK_AVAILABLE"}
    overall_ok = checks["core"]["status"] in ok_states and all(
        checks[name]["status"] in ok_states for name in profiles
    )
    return {
        "overall": "OK" if overall_ok else "INCOMPLETE",
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Output rendering
# ---------------------------------------------------------------------------


def _human_bytes(value: int | None) -> str:
    if value is None:
        return "?"
    size = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{int(size)} B" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TiB"


def _print_step_extras(step: dict, indent: str = "      ") -> None:
    for key in ("note", "warning", "notice"):
        if step.get(key):
            for line in str(step[key]).splitlines():
                print(f"{indent}{key}: {line}")
    if step.get("command"):
        print(f"{indent}command: {step['command']}")
    if step.get("config"):
        print(f"{indent}config snippet:")
        print(json.dumps(step["config"], indent=8))
    if step.get("links"):
        for key, url in step["links"].items():
            print(f"{indent}{key}: {url}")


def print_plan(ctx: Context) -> None:
    print(f"Plan for platform {ctx.platform} (dry-run, nothing touched):\n")
    for index, step in enumerate(ctx.steps, 1):
        bits = [f"[{index:>2}] ({step['profile']}) {step['action']}"]
        if step.get("status"):
            bits.append(f"[{step['status']}]")
        if step.get("url"):
            bits.append(f"url={step['url']}")
        if step.get("dest"):
            bits.append(f"-> {step['dest']}")
        if step.get("bytes"):
            bits.append(_human_bytes(step["bytes"]))
        bits.extend(
            f"{k}={step[k][:12]}…" for k in ("sha256", "md5") if step.get(k)
        )
        print("  " + " | ".join(bits))
        _print_step_extras(step)
    print()


def print_steps_summary(ctx: Context) -> None:
    print(f"\nResults (platform {ctx.platform}):\n")
    for step in ctx.steps:
        status = step.get("status", "?")
        extra = []
        for key in ("binary", "dest", "dmnd_sha256", "detail"):
            if step.get(key):
                extra.append(f"{key}={step[key]}")
        print(
            f"  ({step['profile']}) {step['action']}: {status}"
            + ("  " + " ".join(extra) if extra else "")
        )
        _print_step_extras(step)
    manual = [s for s in ctx.steps if s.get("status") == "manual_action_required"]
    if manual:
        print(
            f"\n{len(manual)} step(s) need manual action (see notes above); "
            "all attempted downloads and checksums succeeded."
        )
    mismatched = [
        s for s in ctx.steps if s.get("status") == "contract_mismatch"
    ]
    if mismatched:
        actions = ", ".join(s["action"] for s in mismatched)
        print(
            f"ERROR: contract verification failed for {actions}; fresh "
            "reconstruction is not ready (AlphaGEM.py --profile paper "
            "hard-stops on this mismatch).",
            file=sys.stderr,
        )


def print_check_summary(result: dict) -> None:
    print(f"\nPreflight: {result['overall']}\n")
    name_w = max(len(n) for n in CHECK_ORDER)
    for name in CHECK_ORDER:
        check = result["checks"][name]
        print(f"  {name:<{name_w}}  {check['status']}")
        for key in (
            "missing_binaries",
            "missing_python_packages",
            "dirty_git_entries",
            "interpreter_used",
            "runtime_files_note",
            "tool_root_note",
            "note",
            "deepec_warning",
            "gff_note",
            "detail",
        ):
            value = check.get(key)
            if not value:
                continue
            text = value if isinstance(value, str) else json.dumps(value)
            for line in str(text).splitlines():
                print(f"    {'':<{name_w}}  {key}: {line}")
        links = check.get("links")
        if isinstance(links, dict) and check.get("status") != "READY":
            for lkey, lurl in links.items():
                print(f"    {'':<{name_w}}  {lkey}: {lurl}")
        for item in check.get("missing_files", [])[:10]:
            print(
                f"    {'':<{name_w}}  missing file: {item['path']} ({item['hint']})"
            )
        for item in check.get("contract_mismatches", []):
            expected = ""
            if item.get("expected_file_count") is not None:
                expected += f"; expected {item['expected_file_count']} files"
            if item.get("expected_bytes") is not None:
                expected += f"; expected {item['expected_bytes']} bytes"
            if item.get("expected_sha256"):
                expected += f" / sha256 {item['expected_sha256']}"
            print(
                f"    {'':<{name_w}}  contract mismatch: {item['path']} "
                f"({item['reason']}{expected})"
            )
        remediation = check.get("remediation")
        if remediation:
            lines = remediation if isinstance(remediation, list) else [remediation]
            for line in lines:
                for sub in str(line).splitlines():
                    print(f"    {'':<{name_w}}  fix: {sub}")
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--profile",
        nargs="+",
        default=None,
        metavar="PROFILE",
        help=(
            "one or more of screening gapfill alphagem annotation deeploc all "
            "(default: screening when fetching, all for --check)"
        ),
    )
    parser.add_argument(
        "--dest",
        type=Path,
        help="umbrella dir: derives --runtime-root=<dest>/alphagem-runtime and "
        "--tool-root=<dest>/annotation-tools when those flags are omitted",
    )
    parser.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR)
    parser.add_argument(
        "--runtime-root", type=Path, help="AlphaGEM runtime dir (outside repo)"
    )
    parser.add_argument("--tool-root", type=Path, help="annotation tools root")
    parser.add_argument(
        "--weights-dir",
        type=Path,
        help="ESM/ProtT5 weights dir (default <runtime-root>/var/weights)",
    )
    parser.add_argument(
        "--with-eggnog",
        action="store_true",
        help="fetch eggnog.db + build fungi.dmnd (~41 GiB decompressed)",
    )
    parser.add_argument(
        "--with-prott5",
        action="store_true",
        help="fetch ProtT5 XL weights (structure_alignment mode only)",
    )
    parser.add_argument(
        "--with-helixer", action="store_true", help="fetch Helixer fungal weights"
    )
    parser.add_argument(
        "--deeploc-tarball",
        type=Path,
        help="licensed DTU DeepLoc-2.x.tar.gz to install into tools/deeploc",
    )
    parser.add_argument("--alphagem-python", help="AlphaGEM env interpreter to probe")
    parser.add_argument("--deepec-python", help="DeepEC env interpreter to probe")
    parser.add_argument("--annevo-python", help="ANNEVO env interpreter to probe")
    parser.add_argument(
        "--esm-model-path", type=Path, help="explicit ESM-1b .pt for --check"
    )
    parser.add_argument("--dry-run", action="store_true", help="print plan only")
    parser.add_argument(
        "--check", action="store_true", help="readiness audit instead of fetching"
    )
    parser.add_argument("--json", action="store_true", help="machine-readable output")
    parser.add_argument(
        "--force", action="store_true", help="re-download/overwrite verified files"
    )
    parser.add_argument("--platform", help="override detected platform (for --dry-run)")
    return parser


def _resolve_roots(
    args: argparse.Namespace,
) -> tuple[Path | None, Path | None, Path | None]:
    runtime = args.runtime_root
    tool_root = args.tool_root
    if args.dest:
        dest = args.dest.expanduser()
        if runtime is None:
            runtime = dest / "alphagem-runtime"
        if tool_root is None:
            tool_root = dest / "annotation-tools"
    runtime = runtime.expanduser().resolve(strict=False) if runtime else None
    tool_root = tool_root.expanduser().resolve(strict=False) if tool_root else None
    weights = args.weights_dir or (runtime / "var" / "weights" if runtime else None)
    weights = weights.expanduser().resolve(strict=False) if weights else None
    return runtime, tool_root, weights


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        platform = detect_platform(args.platform)
    except ValueError as exc:
        parser.error(str(exc))

    profiles_arg = args.profile or (["all"] if args.check else ["screening"])
    try:
        profiles = expand_profiles(profiles_arg)
    except ValueError as exc:
        parser.error(str(exc))

    runtime_root, tool_root, weights_dir = _resolve_roots(args)
    ctx = Context(
        args=args,
        platform=platform,
        bin_dir=args.bin_dir.expanduser().resolve(strict=False),
        runtime_root=runtime_root,
        tool_root=tool_root,
        weights_dir=weights_dir,
        dry_run=args.dry_run,
        force=args.force,
    )

    if args.check:
        result = run_preflight_checks(ctx, profiles)
        if args.json:
            print(json.dumps(result, indent=2, default=str))
        else:
            print_check_summary(result)
        return 0 if result["overall"] == "OK" else 2

    if not args.dry_run:
        try:
            if "alphagem" in profiles:
                if ctx.runtime_root is None:
                    raise FetchError(
                        "--profile alphagem requires --runtime-root (or --dest) outside the repository"
                    )
                _prepare_runtime_module().require_real_runtime_path(ctx.runtime_root)
            if "annotation" in profiles and ctx.tool_root is None:
                raise FetchError(
                    "--profile annotation requires --tool-root (or --dest)"
                )
            tarball = getattr(ctx.args, "deeploc_tarball", None)
            if "deeploc" in profiles and tarball is not None:
                if not Path(tarball).expanduser().is_file():
                    raise FetchError(f"DeepLoc tarball not found: {tarball}")
        except (FetchError, ValueError, RuntimeError) as exc:
            if args.json:
                print(json.dumps({"error": str(exc), "steps": []}, indent=2))
            else:
                print(f"ERROR: {exc}", file=sys.stderr)
            return 2

    for name in profiles:
        try:
            {
                "screening": run_screening,
                "gapfill": run_gapfill,
                "alphagem": run_alphagem,
                "annotation": run_annotation,
                "deeploc": run_deeploc,
            }[name](ctx)
        except (FetchError, ValueError, RuntimeError) as exc:
            if args.json:
                print(
                    json.dumps(
                        {"error": str(exc), "steps": ctx.steps},
                        indent=2,
                        default=str,
                    )
                )
            else:
                print_steps_summary(ctx)
                print(f"ERROR: {exc}", file=sys.stderr)
            return 2

    payload = {
        "platform": platform,
        "profiles": profiles,
        "dry_run": args.dry_run,
        "steps": ctx.steps,
    }
    has_mismatch = any(
        s.get("status") == "contract_mismatch" for s in ctx.steps
    )
    if has_mismatch:
        payload["contract_mismatch"] = True
    if args.json:
        print(json.dumps(payload, indent=2, default=str))
    elif args.dry_run:
        print_plan(ctx)
    else:
        print_steps_summary(ctx)
    return 2 if has_mismatch else 0


if __name__ == "__main__":
    raise SystemExit(main())
