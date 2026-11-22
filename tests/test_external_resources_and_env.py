"""Deterministic tests for fetch_external_resources.py and pipeline env/interpreter wiring."""
from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import os
import stat
import subprocess
import sys
import tarfile
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import geno2flux.pipeline as pipeline
from geno2flux.pipeline import (
    PipelineError,
    resolve_alphagem_python,
    resolve_tool_executable,
)
import scripts.fetch_external_resources as fetcher


def _make_executable(path: Path, body: str = "#!/bin/sh\necho 'mock 2.16.0'\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _make_tar_gz(dest: Path, entries: dict[str, bytes]) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(dest, "w:gz") as tar:
        for name, data in entries.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            info.mode = 0o755
            tar.addfile(info, io.BytesIO(data))
    return hashlib.md5(dest.read_bytes()).hexdigest()


_EGGNOG_STATIC_NAMES = (
    "eggnog_mapper_default.tax.db",
    "eggnog_mapper_default.tax.db.md5",
    "eggnog_mapper_default.tax.db.pkl",
    "eggnog_mapper_default.tax.db.pkl.md5",
    "eggnog_mapper_default.tax.db.tar.gz",
    "eggnog_mapper_default.tax.db.tar.gz.md5",
    "genepred_linkage.sqlite",
    "og2level.tsv",
    "eggnog4.functional_categories.txt",
)


def _write_eggnog_static_data(runtime: Path) -> dict:
    """Populate the 9-file eggnog_static_data bundle and return a patched
    EGGNOG_STATIC_DATA contract matching what was written."""
    data_dir = runtime / fetcher.EGGNOG_STATIC_DATA["path"]
    data_dir.mkdir(parents=True, exist_ok=True)
    for index, name in enumerate(_EGGNOG_STATIC_NAMES):
        (data_dir / name).write_bytes(
            f"static-eggnog-data-{index}\n".encode("utf-8")
        )
    entries = []
    excluded = tuple(
        Path(rel).parts
        for rel in fetcher.EGGNOG_STATIC_DATA["excluded_relative_roots"]
    )
    for path in sorted(data_dir.rglob("*")):
        relative = path.relative_to(data_dir)
        if any(
            relative.parts[: len(prefix)] == prefix for prefix in excluded
        ):
            continue
        if path.is_file():
            entries.append(
                (
                    path.relative_to(runtime).as_posix(),
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    path.stat().st_size,
                )
            )
    entries.sort()
    aggregate = hashlib.sha256()
    for relative, file_sha256, _size in entries:
        aggregate.update(f"{file_sha256}  {relative}\n".encode("utf-8"))
    return {
        **fetcher.EGGNOG_STATIC_DATA,
        "file_count": len(entries),
        "bytes": sum(size for _rel, _sha, size in entries),
        "sha256": aggregate.hexdigest(),
    }


def test_expand_profiles_and_detect_platform():
    assert fetcher.expand_profiles(["screening"]) == ["screening"]
    assert fetcher.expand_profiles(["gapfill"]) == ["screening", "gapfill"]
    assert fetcher.expand_profiles(["alphagem"]) == ["screening", "gapfill", "alphagem"]
    assert fetcher.expand_profiles(["all"]) == [
        "screening",
        "gapfill",
        "alphagem",
        "annotation",
        "deeploc",
    ]
    with pytest.raises(ValueError):
        fetcher.expand_profiles(["unknown-profile"])

    assert fetcher.detect_platform("linux-x86_64") == "linux-x86_64"
    assert fetcher.detect_platform("linux-amd64") == "linux-x86_64"
    assert fetcher.detect_platform("linux-arm64") == "linux-aarch64"
    assert fetcher.detect_platform("darwin-aarch64") == "darwin-aarch64"
    assert fetcher.detect_platform("darwin-" + "arm64") == "darwin-aarch64"
    with pytest.raises(ValueError):
        fetcher.detect_platform("windows-x86_64")


def test_upfront_validation_blocks_downloads_on_invalid_runtime_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Passing an inside-repo runtime root must fail BEFORE downloading BLAST+/DIAMOND."""
    download_called = False

    def fail_if_called(*args, **kwargs):
        nonlocal download_called
        download_called = True
        raise AssertionError("download_file must not run before argument validation")

    monkeypatch.setattr(fetcher, "download_file", fail_if_called)
    inside_repo = fetcher.REPO_ROOT / "_invalid_inside_runtime"
    rc = fetcher.main(
        [
            "--profile",
            "alphagem",
            "--bin-dir",
            str(tmp_path / "bin"),
            "--runtime-root",
            str(inside_repo),
            "--platform",
            "linux-x86_64",
        ]
    )
    assert rc == 2
    assert not download_called

    rc_ann = fetcher.main(
        [
            "--profile",
            "annotation",
            "--bin-dir",
            str(tmp_path / "bin"),
            "--platform",
            "linux-x86_64",
        ]
    )
    assert rc_ann == 2
    assert not download_called


def test_safe_extract_tar_blocks_traversal_and_escaping_symlinks(tmp_path: Path):
    # 1. Path traversal (..)
    bad_rel = tmp_path / "bad_rel.tar.gz"
    with tarfile.open(bad_rel, "w:gz") as tar:
        data = b"evil"
        info = tarfile.TarInfo(name="../escape.txt")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with pytest.raises(fetcher.FetchError):
        fetcher.safe_extract_tar(bad_rel, tmp_path / "out_rel")

    # 2. Absolute path
    bad_abs = tmp_path / "bad_abs.tar.gz"
    with tarfile.open(bad_abs, "w:gz") as tar:
        data = b"evil"
        info = tarfile.TarInfo(name="/tmp/evil_abs.txt")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with pytest.raises(fetcher.FetchError):
        fetcher.safe_extract_tar(bad_abs, tmp_path / "out_abs")

    # 3. Escaping symlink
    bad_link = tmp_path / "bad_link.tar.gz"
    with tarfile.open(bad_link, "w:gz") as tar:
        info = tarfile.TarInfo(name="link_out")
        info.type = tarfile.SYMTYPE
        info.linkname = "../../etc/passwd"
        tar.addfile(info)
    with pytest.raises(fetcher.FetchError):
        fetcher.safe_extract_tar(bad_link, tmp_path / "out_link")

    # 4. Valid archive with strip_components=1
    good = tmp_path / "good.tar.gz"
    _make_tar_gz(good, {"ANNEVO-main/prediction.py": b"print('ok')\n"})
    dest = tmp_path / "out_good"
    res = fetcher.safe_extract_tar(good, dest, strip_components=1)
    assert res["files"] == 1
    assert (dest / "prediction.py").read_text(encoding="utf-8") == "print('ok')\n"


def test_download_file_verifies_checksum_before_promoting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    payload = b"verified-payload-bytes"
    expected_sha = hashlib.sha256(payload).hexdigest()
    expected_md5 = hashlib.md5(payload).hexdigest()

    def fake_curl(url: str, tmp: Path) -> None:
        tmp.write_bytes(payload)

    monkeypatch.setattr(fetcher.shutil, "which", lambda cmd: "/usr/bin/curl" if cmd == "curl" else None)
    monkeypatch.setattr(fetcher, "_curl_download", fake_curl)

    dest = tmp_path / "artifact.bin"
    # Mismatched sha256 must fail and leave dest absent
    with pytest.raises(fetcher.FetchError, match="verification failed"):
        fetcher.download_file(
            "https://example.invalid/artifact.bin",
            dest,
            sha256="0" * 64,
        )
    assert not dest.exists()

    # Matching checksum succeeds; second call skips download
    status1 = fetcher.download_file(
        "https://example.invalid/artifact.bin",
        dest,
        expected_bytes=len(payload),
        sha256=expected_sha,
        md5=expected_md5,
    )
    assert status1 == "downloaded"
    assert dest.read_bytes() == payload

    status2 = fetcher.download_file(
        "https://example.invalid/artifact.bin",
        dest,
        expected_bytes=len(payload),
        sha256=expected_sha,
        md5=expected_md5,
    )
    assert status2 == "already_present"


def test_eggnog_fungi_selection_contract(tmp_path: Path):
    taxid_tsv = tmp_path / "e5.taxid_info.tsv"
    taxid_tsv.write_text(
        "# Taxid\tSci.Name\tRank\tNamed Lineage\n"
        "4932\tSaccharomyces cerevisiae\tspecies\tEukaryota,Opisthokonta,Fungi,Ascomycota\n"
        "9606\tHomo sapiens\tspecies\tEukaryota,Opisthokonta,Metazoa\n"
        "9999\tPseudo organism\tspecies\tEukaryota,NotFungi,Other\n",
        encoding="utf-8",
    )
    taxids = fetcher.select_fungi_taxids(taxid_tsv)
    assert taxids == {"4932"}

    proteomes = tmp_path / "e5.proteomes.faa"
    proteomes.write_text(
        ">4932.YAL001C\nMSEQONE\n>9606.ENSP0001\nMSEQTWO\n>9999.FAKE\nMSEQTHREE\n",
        encoding="utf-8",
    )
    out_faa = tmp_path / "fungi.faa"
    stats = fetcher.filter_proteome_fasta(proteomes, taxids, out_faa)
    assert stats["selected_records"] == 1
    assert out_faa.read_text(encoding="utf-8") == ">4932.YAL001C\nMSEQONE\n"


def test_fetch_gapfill_and_annotation_profiles_install_expected_layout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    bin_dir = tmp_path / "bin"
    tool_root = tmp_path / "annotation-tools"

    diamond_tar = tmp_path / "src_diamond.tar.gz"
    _make_tar_gz(diamond_tar, {"diamond": b"#!/bin/sh\necho 'diamond version 2.1.10'\n"})

    blast_tar = tmp_path / "src_blast.tar.gz"
    blast_md5 = _make_tar_gz(
        blast_tar,
        {
            "ncbi-blast-2.16.0+/bin/tblastn": b"#!/bin/sh\necho 'tblastn: 2.16.0+'\n",
            "ncbi-blast-2.16.0+/bin/makeblastdb": b"#!/bin/sh\necho 'makeblastdb: 2.16.0+'\n",
            "ncbi-blast-2.16.0+/bin/blastp": b"#!/bin/sh\necho 'blastp: 2.16.0+'\n",
        },
    )

    annevo_tar = tmp_path / "src_annevo.tar.gz"
    _make_tar_gz(
        annevo_tar,
        {
            "ANNEVO-main/prediction.py": b"# annevo prediction\n",
            "ANNEVO-main/decoding.py": b"# annevo decoding\n",
        },
    )
    annevo_pt_bytes = b"MOCK_ANNEVO_WEIGHTS"

    patched_blast = {
        **fetcher.BLAST_PLUS,
        "archives": {
            **fetcher.BLAST_PLUS["archives"],
            "linux-x86_64": {
                "file": "ncbi-blast-2.16.0+-x64-linux.tar.gz",
                "md5": blast_md5,
            },
        },
    }
    monkeypatch.setattr(fetcher, "BLAST_PLUS", patched_blast)
    monkeypatch.setattr(
        fetcher,
        "ANNEVO",
        {**fetcher.ANNEVO, "weights_bytes": len(annevo_pt_bytes)},
    )
    # Isolate binary resolution from host PATH
    monkeypatch.setattr(fetcher, "resolve_binary", lambda name, bdir=None: None)

    def fake_download(url: str, dest: Path, **kwargs) -> str:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if "diamond" in url:
            dest.write_bytes(diamond_tar.read_bytes())
        elif "ncbi-blast" in url:
            dest.write_bytes(blast_tar.read_bytes())
        elif "ANNEVO/archive" in url:
            dest.write_bytes(annevo_tar.read_bytes())
        elif "ANNEVO_Fungi.pt" in url:
            dest.write_bytes(annevo_pt_bytes)
        else:
            raise AssertionError(f"Unexpected URL in offline test: {url}")
        return "downloaded"

    monkeypatch.setattr(fetcher, "download_file", fake_download)

    rc = fetcher.main(
        [
            "--profile",
            "gapfill",
            "annotation",
            "deeploc",
            "--bin-dir",
            str(bin_dir),
            "--tool-root",
            str(tool_root),
            "--platform",
            "linux-x86_64",
            "--force",
        ]
    )
    assert rc == 0
    for binary in ("diamond", "tblastn", "makeblastdb", "blastp"):
        p = bin_dir / binary
        assert p.is_file() and os.access(p, os.X_OK), f"missing executable {binary}"

    # Verify paths match scripts/predict_proteome.py --tool-root contract verbatim
    assert (tool_root / "tools" / "annevo" / "annevo-current" / "prediction.py").is_file()
    assert (tool_root / "tools" / "annevo" / "annevo-current" / "decoding.py").is_file()
    assert (tool_root / "weights" / "annevo" / "ANNEVO_Fungi.pt").read_bytes() == annevo_pt_bytes


def test_preflight_check_blocks_alphagem_on_uv_python_and_passes_with_alphagem_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
):
    """Regression test: even when ALL AlphaGEM data + binaries exist,
    --check MUST mark alphagem BLOCKED if the interpreter lacks torch/esm/CLEAN,
    and READY once a valid alphagem_python interpreter is supplied."""
    bin_dir = tmp_path / "bin"
    for name in ("diamond", "blastp", "tblastn", "makeblastdb"):
        _make_executable(bin_dir / name, f"#!/bin/sh\necho '{name} 2.16.0'\n")

    runtime = tmp_path / "alphagem-runtime"
    for rel, _hint in fetcher._ALPHAGEM_REQUIRED_FILES:
        target = runtime / rel
        if rel == "tools/DeepProZyme":
            target.mkdir(parents=True, exist_ok=True)
            (target / "run_deepec.py").write_text("# mock\n", encoding="utf-8")
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("ok\n", encoding="utf-8")
    esm_dir = runtime / "var" / "weights" / "esm"
    esm_dir.mkdir(parents=True, exist_ok=True)
    for fname in fetcher.ESM_WEIGHTS:
        (esm_dir / fname).write_bytes(b"weights")
    mock_ok_bytes = b"ok\n"
    monkeypatch.setattr(
        fetcher,
        "EGGNOG",
        {
            **fetcher.EGGNOG,
            "db_bytes": len(mock_ok_bytes),
            "db_sha256": hashlib.sha256(mock_ok_bytes).hexdigest(),
            "fungi_dmnd_bytes": len(mock_ok_bytes),
            "fungi_dmnd_sha256": hashlib.sha256(mock_ok_bytes).hexdigest(),
        },
    )
    static_contract = _write_eggnog_static_data(runtime)
    monkeypatch.setattr(fetcher, "EGGNOG_STATIC_DATA", static_contract)

    # 1. Under default uv Python (sys.executable), AlphaGEM MUST be BLOCKED for missing torch/esm
    monkeypatch.delenv("G2F_ALPHAGEM_PYTHON", raising=False)
    rc_blocked = fetcher.main(
        [
            "--check",
            "--profile",
            "alphagem",
            "--bin-dir",
            str(bin_dir),
            "--runtime-root",
            str(runtime),
            "--json",
        ]
    )
    out_blocked = json.loads(capsys.readouterr().out)
    assert rc_blocked == 2
    assert out_blocked["overall"] == "INCOMPLETE"
    assert out_blocked["checks"]["core"]["status"] == "READY"
    assert out_blocked["checks"]["screening"]["status"] == "READY"
    assert out_blocked["checks"]["gapfill"]["status"] == "READY"
    ag = out_blocked["checks"]["alphagem"]
    assert ag["status"] == "BLOCKED"
    assert ag["missing_files"] == []
    assert ag["missing_binaries"] == []
    assert ag["contract_mismatches"] == []
    assert "torch" in ag["missing_python_packages"]
    assert "esm" in ag["missing_python_packages"]

    # 2. With a mock AlphaGEM Python interpreter that satisfies all imports -> READY
    fake_ag_py = _make_executable(
        tmp_path / "envs" / "AlphaGEM" / "bin" / "python",
        "#!/bin/sh\necho 'OK'\nexit 0\n",
    )
    fake_deepec_py = _make_executable(
        tmp_path / "envs" / "deepec" / "bin" / "python",
        "#!/bin/sh\necho 'OK'\nexit 0\n",
    )
    rc_ready = fetcher.main(
        [
            "--check",
            "--profile",
            "alphagem",
            "--bin-dir",
            str(bin_dir),
            "--runtime-root",
            str(runtime),
            "--alphagem-python",
            str(fake_ag_py),
            "--deepec-python",
            str(fake_deepec_py),
            "--json",
        ]
    )
    out_ready = json.loads(capsys.readouterr().out)
    assert rc_ready == 0
    assert out_ready["overall"] == "OK"
    assert out_ready["checks"]["alphagem"]["status"] == "READY"
    assert out_ready["checks"]["alphagem"]["interpreter_used"] == str(fake_ag_py)
    # 3. eggnog_static_data contract violations MUST block check-env even when
    # interpreter + binaries + eggnog.db + fungi.dmnd are otherwise fine.
    check_argv = [
        "--check",
        "--profile",
        "alphagem",
        "--bin-dir",
        str(bin_dir),
        "--alphagem-python",
        str(fake_ag_py),
        "--deepec-python",
        str(fake_deepec_py),
        "--json",
    ]

    def build_contract_runtime(case: str, *, write_static: bool) -> Path:
        case_rt = tmp_path / case
        data_dir = case_rt / "tools" / "eggnog_mapper" / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        (data_dir / "eggnog.db").write_bytes(mock_ok_bytes)
        (data_dir / "fungi.dmnd").write_bytes(mock_ok_bytes)
        if write_static:
            _write_eggnog_static_data(case_rt)
        return case_rt

    def run_check(case_rt: Path) -> tuple[int, dict]:
        rc = fetcher.main(check_argv + ["--runtime-root", str(case_rt)])
        return rc, json.loads(capsys.readouterr().out)

    # (i) only eggnog.db + fungi.dmnd, no static data files -> BLOCKED
    rt_missing = build_contract_runtime("rt_missing_static", write_static=False)
    rc, payload = run_check(rt_missing)
    assert rc == 2
    ag_missing = payload["checks"]["alphagem"]
    assert ag_missing["status"] == "BLOCKED"
    assert any(
        m["path"] == "tools/eggnog_mapper/data"
        and m.get("expected_file_count") == 9
        for m in ag_missing["contract_mismatches"]
    )

    # (ii) unexpected extra files inside the bundle -> BLOCKED
    rt_extra = build_contract_runtime("rt_extra_files", write_static=True)
    extra_dir = rt_extra / "tools" / "eggnog_mapper" / "data"
    (extra_dir / "e5.proteomes.faa").write_bytes(b"stray")
    (extra_dir / "unexpected.dat").write_bytes(b"unexpected")
    rc, payload = run_check(rt_extra)
    assert rc == 2
    ag_extra = payload["checks"]["alphagem"]
    assert ag_extra["status"] == "BLOCKED"
    assert any(
        m["path"] == "tools/eggnog_mapper/data"
        for m in ag_extra["contract_mismatches"]
    )
    fetcher.print_check_summary(payload)
    assert "eggnog_mapper" in capsys.readouterr().out

    # (iii) same-length byte mutation of eggnog.db -> BLOCKED via db sha256
    rt_mutated = build_contract_runtime("rt_mutated_db", write_static=True)
    mutated = bytes(b ^ 0xFF for b in mock_ok_bytes)
    (rt_mutated / "tools" / "eggnog_mapper" / "data" / "eggnog.db").write_bytes(
        mutated
    )
    rc, payload = run_check(rt_mutated)
    assert rc == 2
    ag_mutated = payload["checks"]["alphagem"]
    assert ag_mutated["status"] == "BLOCKED"
    assert any(
        m["path"] == "tools/eggnog_mapper/data/eggnog.db"
        for m in ag_mutated["contract_mismatches"]
    )



def test_deeploc_check_returns_fallback_available_with_dtu_links(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
):
    bin_dir = tmp_path / "bin"
    _make_executable(bin_dir / "diamond", "#!/bin/sh\necho 'diamond version 2.1.10'\n")
    monkeypatch.setattr(fetcher, "_resolve_deeploc", lambda: (None, "not_found"))

    rc = fetcher.main(
        [
            "--check",
            "--profile",
            "screening",
            "deeploc",
            "--bin-dir",
            str(bin_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert rc == 0
    assert payload["overall"] == "OK"
    dl = payload["checks"]["deeploc"]
    assert dl["status"] == "FALLBACK_AVAILABLE"
    assert "services.healthtech.dtu.dk" in dl["links"]["deeploc_2_1_request"]
    assert "services.healthtech.dtu.dk" in dl["links"]["deeploc_2_0_request"]


def test_resolve_alphagem_python_and_tool_executable(tmp_path: Path):
    cfg_py = _make_executable(tmp_path / "cfg_env" / "bin" / "python", "#!/bin/sh\nexit 0\n")
    env_py = _make_executable(tmp_path / "env_env" / "bin" / "python", "#!/bin/sh\nexit 0\n")

    # 1. Explicit config beats G2F_ALPHAGEM_PYTHON
    py, src = resolve_alphagem_python(
        str(cfg_py),
        base_dir=tmp_path,
        env={"G2F_ALPHAGEM_PYTHON": str(env_py)},
    )
    assert py == str(cfg_py.resolve())
    assert src == "configured"

    # 2. G2F_ALPHAGEM_PYTHON beats sys.executable
    py, src = resolve_alphagem_python(
        None,
        base_dir=tmp_path,
        env={"G2F_ALPHAGEM_PYTHON": str(env_py)},
    )
    assert py == str(env_py.resolve())
    assert src == "env_G2F_ALPHAGEM_PYTHON"

    # 3. Fallback to sys.executable
    py, src = resolve_alphagem_python(None, base_dir=tmp_path, env={})
    assert py == sys.executable
    assert src == "current_python"

    # 4. Explicit missing path raises PipelineError
    with pytest.raises(PipelineError):
        resolve_alphagem_python(str(tmp_path / "no_such_python"), base_dir=tmp_path, env={})

    # 5. Tool binary auto-discovery: extra_bin_dirs -> G2F_TOOLS_BIN -> repo/tools/bin -> unchanged
    fake_repo = tmp_path / "repo"
    repo_diamond = _make_executable(fake_repo / "tools" / "bin" / "diamond")
    extra_diamond = _make_executable(tmp_path / "extra_bin" / "diamond")
    env_diamond = _make_executable(tmp_path / "env_bin" / "diamond")

    assert (
        resolve_tool_executable(
            "diamond",
            extra_bin_dirs=[tmp_path / "extra_bin"],
            env={"G2F_TOOLS_BIN": str(tmp_path / "env_bin")},
            repo_root=fake_repo,
        )
        == str(extra_diamond.resolve())
    )
    assert (
        resolve_tool_executable(
            "diamond",
            env={"G2F_TOOLS_BIN": str(tmp_path / "env_bin")},
            repo_root=fake_repo,
        )
        == str(env_diamond.resolve())
    )
    assert (
        resolve_tool_executable("diamond", env={}, repo_root=fake_repo)
        == str(repo_diamond.resolve())
    )
    # Missing drop-in keeps bare name unchanged for PATH lookup / test mocks
    assert resolve_tool_executable("diamond", env={}, repo_root=tmp_path / "empty") == "diamond"


def _write_maplist_xlsx(path: Path) -> Path:
    """Write a real AlphaGEM structure maplist workbook (pd.read_excel format)."""
    import pandas as pd

    pd.DataFrame(
        [
            {
                "Entry": "P001",
                "Entry Name": "P001_YEAST",
                "Structure Data": "P001",
            }
        ]
    ).to_excel(str(path), index=False)
    return path


def test_pipeline_fresh_steps_use_alphagem_python_and_hint_on_default_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from tests.test_pipeline_contracts import (
        ExternalBoundaryStop,
        _fake_runtime,
        _light_stubs,
        _tiny_model_pair,
        _write,
    )

    rt = _fake_runtime(tmp_path)
    custom_py = _make_executable(tmp_path / "alphagem_env" / "bin" / "python", "#!/bin/sh\nexit 0\n")
    prott5 = tmp_path / "prot_t5_snapshot"
    prott5.mkdir()
    maplist = _write_maplist_xlsx(tmp_path / "maplist.xlsx")
    struct_dir = tmp_path / "structures"
    struct_dir.mkdir()
    (struct_dir / "p1.pdb").write_text("END\n", encoding="utf-8")

    cfg = _write(
        tmp_path,
        proteome_sha=True,
        runtime_root=str(rt),
        extra={
            "alphagem_python": str(custom_py),
            "alphagem_mode": "structure_alignment",
            "prottrans_model_path": str(prott5),
            "structure_maplist": str(maplist),
            "structure_dir": str(struct_dir),
        },
    )
    cfg_data = json.loads(cfg.read_text(encoding="utf-8"))
    cfg_data.pop("withgaps_model", None)
    cfg.write_text(json.dumps(cfg_data), encoding="utf-8")
    out = tmp_path / "out_custom_py"
    calls: list[list[str]] = []

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        name = cmd[cmd.index("--name") + 1]
        yml = Path(kw["cwd"]) / "working" / name / f"{name}-GEM_withgaps.yml"
        yml.parent.mkdir(parents=True, exist_ok=True)
        model, _ = _tiny_model_pair()
        import cobra

        cobra.io.save_yaml_model(model, str(yml))
        return MagicMock(returncode=0)

    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        stack.enter_context(patch.object(pipeline.subprocess, "run", side_effect=fake_run))
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(cfg, out, localization_order_override="pre-gapfill")

    assert calls and calls[0][0] == str(custom_py.resolve())
    assert calls[0][calls[0].index("--mode") + 1] == "structure_alignment"
    assert calls[0][calls[0].index("--prottrans-model-path") + 1] == str(prott5)
    assert calls[0][calls[0].index("--maplist") + 1] == str(maplist)
    assert calls[0][calls[0].index("--structure") + 1] == str(struct_dir)
    manifest = json.loads((out / "pipeline_manifest.json").read_text(encoding="utf-8"))
    step13 = next(s for s in manifest["steps"] if s["step"] == "step1_3_withgaps")
    assert step13["detail"]["alphagem_python"] == str(custom_py.resolve())
    assert step13["detail"]["alphagem_python_source"] == "configured"

    # When alphagem_python is omitted and AlphaGEM.py fails under default uv Python,
    # pipeline raises PipelineError with actionable alphagem_python / check-env hint.
    monkeypatch.delenv("G2F_ALPHAGEM_PYTHON", raising=False)
    case_fail = tmp_path / "case_fail"
    case_fail.mkdir()
    cfg_default = _write(
        case_fail,
        proteome_sha=True,
        runtime_root=str(rt),
    )
    cfg_default_data = json.loads(cfg_default.read_text(encoding="utf-8"))
    cfg_default_data.pop("withgaps_model", None)
    cfg_default.write_text(json.dumps(cfg_default_data), encoding="utf-8")

    def failing_run(cmd, **kw):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd)

    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        stack.enter_context(patch.object(pipeline.subprocess, "run", side_effect=failing_run))
        with pytest.raises(PipelineError, match="alphagem_python"):
            pipeline.run_pipeline(
                cfg_default,
                tmp_path / "out_fail_hint",
                localization_order_override="pre-gapfill",
            )


def test_alphagem_fetch_keeps_runtime_git_clean_and_rerunnable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Fetching --profile alphagem must place download caches & default weights under
    the gitignored var/ dir so formal AlphaGEM git cleanliness checks pass, and a
    second run without --force must succeed idempotently."""
    import scripts.prepare_alphagem_runtime as prep

    bin_dir = tmp_path / "bin"
    for name in ("diamond", "tblastn", "makeblastdb", "blastp"):
        _make_executable(bin_dir / name, f"#!/bin/sh\necho '{name} 2.16.0'\n")

    runtime = tmp_path / "alphagem-runtime"
    ds_tar = tmp_path / "src_data_source.tar.xz"
    with tarfile.open(ds_tar, "w:xz") as tar:
        entries = {
            "data_available/yeast.fasta": b">Y1\nMAAA\n",
            "data_available/yeast.xlsx": b"xlsx",
            "rhea/rhea2uniprot_sprot.tsv": b"ID\tRHEA\n",
            "rhea/rhea.rdf": b"<rdf/>\n",
            "plmsearchtools/model/plmsearch.sav": b"sav",
            "tools/CLEAN/app/CLEAN_infer_fasta.py": b"# clean\n",
            "tools/DeepProZyme/run.py": b"# deepec\n",
            "tools/eggnog_mapper/emapper.py": b"# emapper\n",
        }
        for i in range(1, 7):
            entries[f"rhea/rhea_protembedding{i}.pkl"] = b"pkl"
        for rel, data in entries.items():
            info = tarfile.TarInfo(name=rel)
            info.size = len(data)
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(data))
    ds_bytes = ds_tar.read_bytes()
    ds_sha = hashlib.sha256(ds_bytes).hexdigest()
    ds_md5 = hashlib.md5(ds_bytes).hexdigest()

    esm_pt = b"MOCK_ESM1B_PT"
    esm_reg = b"MOCK_ESM1B_REG"
    esm_reg_md5 = hashlib.md5(esm_reg).hexdigest()

    monkeypatch.setattr(
        fetcher,
        "ALPHAGEM_DATA_SOURCE",
        {
            **fetcher.ALPHAGEM_DATA_SOURCE,
            "bytes": len(ds_bytes),
            "sha256": ds_sha,
            "md5": ds_md5,
        },
    )
    monkeypatch.setattr(
        fetcher,
        "ESM_WEIGHTS",
        {
            "esm1b_t33_650M_UR50S.pt": {
                "url": fetcher.ESM_WEIGHTS["esm1b_t33_650M_UR50S.pt"]["url"],
                "bytes": len(esm_pt),
            },
            "esm1b_t33_650M_UR50S-contact-regression.pt": {
                "url": fetcher.ESM_WEIGHTS["esm1b_t33_650M_UR50S-contact-regression.pt"]["url"],
                "bytes": len(esm_reg),
                "md5": esm_reg_md5,
            },
        },
    )

    def fake_curl(url: str, tmp: Path) -> None:
        if "figshare" in url:
            tmp.write_bytes(ds_bytes)
        elif url.endswith("esm1b_t33_650M_UR50S.pt"):
            tmp.write_bytes(esm_pt)
        elif url.endswith("esm1b_t33_650M_UR50S-contact-regression.pt"):
            tmp.write_bytes(esm_reg)
        else:
            raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(fetcher.shutil, "which", lambda cmd: "/usr/bin/curl" if cmd == "curl" else None)
    monkeypatch.setattr(fetcher, "_curl_download", fake_curl)

    argv = [
        "--profile",
        "alphagem",
        "--bin-dir",
        str(bin_dir),
        "--runtime-root",
        str(runtime),
        "--platform",
        "linux-x86_64",
    ]
    # First run: stages source, downloads into var/, extracts cleanly
    assert fetcher.main(argv) == 0
    assert (runtime / "var" / "downloads" / "data_source.tar.xz").is_file()
    assert (runtime / "var" / "weights" / "esm" / "esm1b_t33_650M_UR50S.pt").is_file()
    assert not (runtime / "_downloads").exists()
    assert not (runtime / "weights").exists()

    # Runtime Git tree MUST be clean for AlphaGEM --profile paper
    git_status = subprocess.run(
        ["git", "-C", str(runtime), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    assert git_status == "", f"runtime git tree is dirty: {git_status!r}"
    _manifest, records = prep.source_records()
    commit = prep.runtime_source_commit(runtime, records)
    assert len(commit) == 40

    # Second run (without --force) must be cleanly idempotent
    assert fetcher.main(argv) == 0
    git_status_after = subprocess.run(
        ["git", "-C", str(runtime), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    assert git_status_after == ""


def test_structure_alignment_fresh_steps_fail_fast_before_runtime_copy(tmp_path: Path):
    """alphagem_mode=structure_alignment with fresh steps 1-3 must fail fast
    before copying the AlphaGEM runtime when required inputs are missing/empty,
    while reusing withgaps_model must not require them."""
    from tests.test_pipeline_contracts import (
        ExternalBoundaryStop,
        _fake_runtime,
        _light_stubs,
        _write,
    )

    rt = _fake_runtime(tmp_path)
    # 1. Missing all three required keys -> PipelineError before runtime view copy
    c1_dir = tmp_path / "case_missing"
    c1_dir.mkdir()
    cfg_missing = _write(
        c1_dir,
        proteome_sha=True,
        runtime_root=str(rt),
        extra={"alphagem_mode": "structure_alignment"},
    )
    data = json.loads(cfg_missing.read_text(encoding="utf-8"))
    data.pop("withgaps_model", None)
    cfg_missing.write_text(json.dumps(data), encoding="utf-8")

    out1 = tmp_path / "out_missing"
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        mock_sub = stack.enter_context(patch.object(pipeline.subprocess, "run"))
        with pytest.raises(PipelineError) as exc_info:
            pipeline.run_pipeline(cfg_missing, out1, localization_order_override="pre-gapfill")
    msg = str(exc_info.value)
    assert "structure_maplist" in msg
    assert "structure_dir" in msg
    assert "prottrans_model_path" in msg
    assert not (out1 / "alphagem_runtime_view").exists()
    mock_sub.assert_not_called()

    # 2. structure_dir present but contains no *.pdb -> PipelineError before copy
    c2_dir = tmp_path / "case_empty_pdb"
    c2_dir.mkdir()
    maplist = _write_maplist_xlsx(c2_dir / "maplist.xlsx")
    empty_struct = c2_dir / "empty_pdbs"
    empty_struct.mkdir()
    prott5 = c2_dir / "prott5"
    prott5.mkdir()
    cfg_empty = _write(
        c2_dir,
        proteome_sha=True,
        runtime_root=str(rt),
        extra={
            "alphagem_mode": "structure_alignment",
            "structure_maplist": str(maplist),
            "structure_dir": str(empty_struct),
            "prottrans_model_path": str(prott5),
        },
    )
    data2 = json.loads(cfg_empty.read_text(encoding="utf-8"))
    data2.pop("withgaps_model", None)
    cfg_empty.write_text(json.dumps(data2), encoding="utf-8")

    out2 = tmp_path / "out_empty_pdb"
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        mock_sub2 = stack.enter_context(patch.object(pipeline.subprocess, "run"))
        with pytest.raises(PipelineError, match=r"\*\.pdb"):
            pipeline.run_pipeline(cfg_empty, out2, localization_order_override="pre-gapfill")
    assert not (out2 / "alphagem_runtime_view").exists()
    mock_sub2.assert_not_called()

    # 3. Reusing withgaps_model bypasses fresh steps 1-3 -> does NOT require structure keys
    c3_dir = tmp_path / "case_reuse"
    c3_dir.mkdir()
    cfg_reuse = _write(
        c3_dir,
        proteome_sha=True,
        withgaps_sha=True,
        runtime_root=str(rt),
        extra={"alphagem_mode": "structure_alignment"},
    )
    out3 = tmp_path / "out_reuse"
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(cfg_reuse, out3, localization_order_override="pre-gapfill")


def test_fungi_dmnd_contract_mismatch_blocks_fetch_and_check_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
):
    """fetch-resources --with-eggnog and check-env must both flag fungi.dmnd
    contract mismatches (exit 2 / BLOCKED) and pass once contract bytes+sha256 match."""
    import gzip

    bin_dir = tmp_path / "bin"
    for name in ("tblastn", "makeblastdb", "blastp"):
        _make_executable(bin_dir / name, f"#!/bin/sh\necho '{name} 2.16.0'\n")

    expected_dmnd = b"VERIFIED_LINUX_DIAMOND_2_1_10_FUNGI_DB"
    wrong_dmnd = b"WRONG_PLATFORM_OR_VERSION_DMND"
    dmnd_mode = {"payload": wrong_dmnd}

    # Mock diamond binary: writes dmnd_mode['payload'] to `< -d arg >.dmnd`
    diamond_bin = _make_executable(
        bin_dir / "diamond",
        "#!/bin/sh\necho 'diamond version 2.1.10'\n",
    )

    real_run_command = fetcher.run_command

    def fake_run_command(cmd, **kwargs):
        if len(cmd) >= 2 and Path(cmd[0]).name == "diamond" and cmd[1] == "makedb":
            out_prefix = Path(cmd[cmd.index("-d") + 1])
            dmnd_file = out_prefix.with_suffix(".dmnd")
            dmnd_file.parent.mkdir(parents=True, exist_ok=True)
            dmnd_file.write_bytes(dmnd_mode["payload"])
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
        return real_run_command(cmd, **kwargs)

    monkeypatch.setattr(fetcher, "run_command", fake_run_command)

    runtime = tmp_path / "alphagem-runtime"
    ds_tar = tmp_path / "src_data_source.tar.xz"
    with tarfile.open(ds_tar, "w:xz") as tar:
        entries = {
            "data_available/yeast.fasta": b">Y1\nMAAA\n",
            "data_available/yeast.xlsx": b"xlsx",
            "rhea/rhea2uniprot_sprot.tsv": b"ID\tRHEA\n",
            "rhea/rhea.rdf": b"<rdf/>\n",
            "plmsearchtools/model/plmsearch.sav": b"sav",
            "tools/CLEAN/app/CLEAN_infer_fasta.py": b"# clean\n",
            "tools/DeepProZyme/run.py": b"# deepec\n",
            "tools/eggnog_mapper/emapper.py": b"# emapper\n",
        }
        for i in range(1, 7):
            entries[f"rhea/rhea_protembedding{i}.pkl"] = b"pkl"
        for rel, data in entries.items():
            info = tarfile.TarInfo(name=rel)
            info.size = len(data)
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(data))
    ds_bytes = ds_tar.read_bytes()
    esm_pt = b"MOCK_ESM1B_PT"
    esm_reg = b"MOCK_ESM1B_REG"

    eggnog_db_raw = b"MOCK_EGGNOG_SQLITE_DB"
    eggnog_db_gz = gzip.compress(eggnog_db_raw)
    taxid_bytes = (
        b"# Taxid\tSci.Name\tRank\tNamed Lineage\n"
        b"4932\tSaccharomyces cerevisiae\tspecies\tEukaryota,Opisthokonta,Fungi,Ascomycota\n"
    )
    proteomes_bytes = b">4932.YAL001C\nMSEQONE\n"
    selected_bytes = b">4932.YAL001C\nMSEQONE\n"

    monkeypatch.setattr(
        fetcher,
        "ALPHAGEM_DATA_SOURCE",
        {
            **fetcher.ALPHAGEM_DATA_SOURCE,
            "bytes": len(ds_bytes),
            "sha256": hashlib.sha256(ds_bytes).hexdigest(),
            "md5": hashlib.md5(ds_bytes).hexdigest(),
        },
    )
    monkeypatch.setattr(
        fetcher,
        "ESM_WEIGHTS",
        {
            "esm1b_t33_650M_UR50S.pt": {
                "url": fetcher.ESM_WEIGHTS["esm1b_t33_650M_UR50S.pt"]["url"],
                "bytes": len(esm_pt),
            },
            "esm1b_t33_650M_UR50S-contact-regression.pt": {
                "url": fetcher.ESM_WEIGHTS["esm1b_t33_650M_UR50S-contact-regression.pt"]["url"],
                "bytes": len(esm_reg),
                "md5": hashlib.md5(esm_reg).hexdigest(),
            },
        },
    )
    monkeypatch.setattr(
        fetcher,
        "EGGNOG",
        {
            **fetcher.EGGNOG,
            "db_bytes": len(eggnog_db_raw),
            "db_sha256": hashlib.sha256(eggnog_db_raw).hexdigest(),
            "proteomes_sha256": hashlib.sha256(proteomes_bytes).hexdigest(),
            "taxid_sha256": hashlib.sha256(taxid_bytes).hexdigest(),
            "selected_fasta_sha256": hashlib.sha256(selected_bytes).hexdigest(),
            "fungi_dmnd_bytes": len(expected_dmnd),
            "fungi_dmnd_sha256": hashlib.sha256(expected_dmnd).hexdigest(),
        },
    )

    def fake_curl(url: str, tmp: Path) -> None:
        if "figshare" in url:
            tmp.write_bytes(ds_bytes)
        elif url.endswith("esm1b_t33_650M_UR50S.pt"):
            tmp.write_bytes(esm_pt)
        elif url.endswith("esm1b_t33_650M_UR50S-contact-regression.pt"):
            tmp.write_bytes(esm_reg)
        elif url.endswith("eggnog.db.gz"):
            tmp.write_bytes(eggnog_db_gz)
        elif url.endswith("e5.proteomes.faa"):
            tmp.write_bytes(proteomes_bytes)
        elif url.endswith("e5.taxid_info.tsv"):
            tmp.write_bytes(taxid_bytes)
        else:
            raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(fetcher.shutil, "which", lambda cmd: "/usr/bin/curl" if cmd == "curl" else None)
    monkeypatch.setattr(fetcher, "_curl_download", fake_curl)

    fetch_argv = [
        "--profile",
        "alphagem",
        "--with-eggnog",
        "--bin-dir",
        str(bin_dir),
        "--runtime-root",
        str(runtime),
        "--platform",
        "linux-x86_64",
        "--json",
    ]
    # 1. Mismatched fungi.dmnd -> fetch exits 2 and reports contract_mismatch
    rc_bad = fetcher.main(fetch_argv)
    payload_bad = json.loads(capsys.readouterr().out)
    assert rc_bad == 2
    assert payload_bad["contract_mismatch"] is True
    dmnd_step = next(s for s in payload_bad["steps"] if s["action"] == "build_fungi_dmnd")
    assert dmnd_step["status"] == "contract_mismatch"
    assert "AlphaGEM.py --profile paper" in dmnd_step["warning"]

    # 2. check-env against the mismatched fungi.dmnd -> BLOCKED with contract_mismatches
    fake_ag_py = _make_executable(
        tmp_path / "envs" / "AlphaGEM" / "bin" / "python",
        "#!/bin/sh\necho 'OK'\nexit 0\n",
    )
    fake_deepec_py = _make_executable(
        tmp_path / "envs" / "deepec" / "bin" / "python",
        "#!/bin/sh\necho 'OK'\nexit 0\n",
    )
    check_argv = [
        "--check",
        "--profile",
        "alphagem",
        "--bin-dir",
        str(bin_dir),
        "--runtime-root",
        str(runtime),
        "--alphagem-python",
        str(fake_ag_py),
        "--deepec-python",
        str(fake_deepec_py),
        "--json",
    ]
    rc_check_bad = fetcher.main(check_argv)
    check_bad = json.loads(capsys.readouterr().out)
    assert rc_check_bad == 2
    ag_bad = check_bad["checks"]["alphagem"]
    assert ag_bad["status"] == "BLOCKED"
    assert any(m["path"] == "tools/eggnog_mapper/data/fungi.dmnd" for m in ag_bad["contract_mismatches"])

    # 3. Rebuild with contract-matching diamond output -> exits 0, rerun skips, check-env READY
    dmnd_mode["payload"] = expected_dmnd
    rc_good = fetcher.main(fetch_argv)
    payload_good = json.loads(capsys.readouterr().out)
    assert rc_good == 0
    assert "contract_mismatch" not in payload_good
    dmnd_good = next(s for s in payload_good["steps"] if s["action"] == "build_fungi_dmnd")
    assert dmnd_good["status"] == "built"
    assert dmnd_good["dmnd_bytes"] == len(expected_dmnd)

    rc_rerun = fetcher.main(fetch_argv)
    payload_rerun = json.loads(capsys.readouterr().out)
    assert rc_rerun == 0
    dmnd_rerun = next(s for s in payload_rerun["steps"] if s["action"] == "build_fungi_dmnd")
    assert dmnd_rerun["status"] == "already_present"

    # The 9-file eggnog_static_data bundle is staged alongside the fetched db
    # and the rebuilt fungi.dmnd (as a real prepared runtime would have it).
    monkeypatch.setattr(
        fetcher, "EGGNOG_STATIC_DATA", _write_eggnog_static_data(runtime)
    )

    rc_check_ok = fetcher.main(check_argv)
    check_ok = json.loads(capsys.readouterr().out)
    assert rc_check_ok == 0
    assert check_ok["checks"]["alphagem"]["status"] == "READY"
    assert check_ok["checks"]["alphagem"]["contract_mismatches"] == []
    # Verify eggNOG intermediate files stayed in var/downloads/eggnog/, NOT tools/eggnog_mapper/data/
    eggnog_data = runtime / "tools" / "eggnog_mapper" / "data"
    assert sorted(p.name for p in eggnog_data.iterdir()) == sorted(
        ["eggnog.db", "fungi.dmnd", *list(_EGGNOG_STATIC_NAMES)]
    )
    eggnog_cache = runtime / "var" / "downloads" / "eggnog"
    for cached in ("eggnog.db.gz", "e5.proteomes.faa", "e5.taxid_info.tsv", "e5.proteomes.fungi.faa"):
        assert (eggnog_cache / cached).is_file()

    # A stray intermediate file inside tools/eggnog_mapper/data/ must block check-env
    # (AlphaGEM --profile paper validates eggnog_static_data file_count/bytes/sha256)
    stray = eggnog_data / "e5.proteomes.faa"
    stray.write_bytes(b"stray")
    (eggnog_data / "unexpected.dat").write_bytes(b"unexpected")
    rc_stray = fetcher.main(check_argv)
    check_stray = json.loads(capsys.readouterr().out)
    assert rc_stray == 2
    ag_stray = check_stray["checks"]["alphagem"]
    assert ag_stray["status"] == "BLOCKED"
    assert any(
        m["path"] == "tools/eggnog_mapper/data/e5.proteomes.faa"
        for m in ag_stray["contract_mismatches"]
    )
    # The aggregate eggnog_static_data bundle contract must also flag extras
    assert any(
        m["path"] == "tools/eggnog_mapper/data"
        and m.get("expected_file_count") == 9
        for m in ag_stray["contract_mismatches"]
    )
    # Human-readable summary rendering must handle stray + bundle entries
    fetcher.print_check_summary(check_stray)
    rendered = capsys.readouterr().out
    assert "tools/eggnog_mapper/data/e5.proteomes.faa" in rendered
    stray.unlink()
    (eggnog_data / "unexpected.dat").unlink()

    # Missing static data files (only eggnog.db + fungi.dmnd) must block check-env
    for static_name in _EGGNOG_STATIC_NAMES:
        (eggnog_data / static_name).unlink()
    rc_missing = fetcher.main(check_argv)
    check_missing = json.loads(capsys.readouterr().out)
    assert rc_missing == 2
    ag_missing = check_missing["checks"]["alphagem"]
    assert ag_missing["status"] == "BLOCKED"
    assert any(
        m["path"] == "tools/eggnog_mapper/data"
        and m.get("expected_file_count") == 9
        for m in ag_missing["contract_mismatches"]
    )
    _write_eggnog_static_data(runtime)

    # A same-length byte mutation of eggnog.db must be caught by db sha256
    db_path = eggnog_data / "eggnog.db"
    original_db = db_path.read_bytes()
    db_path.write_bytes(b"X" * len(original_db))
    rc_mutated = fetcher.main(check_argv)
    check_mutated = json.loads(capsys.readouterr().out)
    assert rc_mutated == 2
    ag_mutated = check_mutated["checks"]["alphagem"]
    assert ag_mutated["status"] == "BLOCKED"
    assert any(
        m["path"] == "tools/eggnog_mapper/data/eggnog.db"
        for m in ag_mutated["contract_mismatches"]
    )
    db_path.write_bytes(original_db)
    assert diamond_bin.is_file()


def test_alphagem_config_snippet_and_plmsearch_flag_gating(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    """Default --profile alphagem config snippet must emit
    prottrans_model_path=null, the emitted config must pass fresh plmsearch
    validation, and AlphaGEM.py must never see structure-alignment flags in
    plmsearch mode."""
    from tests.test_pipeline_contracts import (
        ExternalBoundaryStop,
        _fake_runtime,
        _light_stubs,
        _tiny_model_pair,
        _write,
    )

    bin_dir = tmp_path / "bin"
    _make_executable(bin_dir / "diamond", "#!/bin/sh\necho 'diamond 2.1.10'\n")

    # Default dry-run snippet: no --with-prott5 -> prottrans_model_path is null
    rt = _fake_runtime(tmp_path)
    rc = fetcher.main(
        [
            "--profile", "alphagem", "--dry-run", "--json",
            "--bin-dir", str(bin_dir),
            "--runtime-root", str(rt),
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert rc == 0
    snippet_step = next(s for s in payload["steps"] if s["action"] == "config_snippet")
    snippet = snippet_step["config"]
    assert snippet["prottrans_model_path"] is None

    # With --with-prott5 the snippet must point at prot_t5_xl_uniref50
    rc = fetcher.main(
        [
            "--profile", "alphagem", "--dry-run", "--json",
            "--bin-dir", str(bin_dir),
            "--runtime-root", str(rt),
            "--with-prott5",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert rc == 0
    snippet_p5 = next(s for s in payload["steps"] if s["action"] == "config_snippet")[
        "config"
    ]
    assert snippet_p5["prottrans_model_path"].endswith("prot_t5_xl_uniref50")

    # The default snippet merged into a pipeline config must pass fresh
    # plmsearch validation and never forward structure-alignment flags.
    custom_py = _make_executable(
        tmp_path / "alphagem_env" / "bin" / "python", "#!/bin/sh\nexit 0\n"
    )
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        runtime_root=str(rt),
        extra={**snippet, "alphagem_python": str(custom_py)},
    )
    cfg_data = json.loads(cfg.read_text(encoding="utf-8"))
    cfg_data.pop("withgaps_model", None)
    # The snippet's deepec_python is a documentation placeholder, not a path
    cfg_data.pop("deepec_python", None)
    cfg.write_text(json.dumps(cfg_data), encoding="utf-8")

    calls: list[list[str]] = []

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        name = cmd[cmd.index("--name") + 1]
        yml = Path(kw["cwd"]) / "working" / name / f"{name}-GEM_withgaps.yml"
        yml.parent.mkdir(parents=True, exist_ok=True)
        model, _ = _tiny_model_pair()
        import cobra

        cobra.io.save_yaml_model(model, str(yml))
        return MagicMock(returncode=0)

    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        stack.enter_context(
            patch.object(pipeline.subprocess, "run", side_effect=fake_run)
        )
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, tmp_path / "out_plmsearch",
                localization_order_override="pre-gapfill",
            )

    assert calls
    assert calls[0][calls[0].index("--mode") + 1] == "plmsearch"
    for forbidden in ("--prottrans-model-path", "--maplist", "--structure"):
        assert forbidden not in calls[0]


def test_stale_structure_paths_do_not_block_plmsearch_or_skip(
    tmp_path: Path, synthetic_screening_case: Path
):
    """Stale/nonexistent prottrans_model_path, structure_dir, and
    structure_maplist must not fail plmsearch or localization_order='skip'
    runs; plmsearch must never forward those flags."""
    from tests.test_pipeline_contracts import (
        ExternalBoundaryStop,
        _fake_runtime,
        _light_stubs,
        _tiny_model_pair,
        _write,
    )

    stale = {
        "prottrans_model_path": str(tmp_path / "missing-prott5"),
        "structure_dir": str(tmp_path / "missing-structures"),
        "structure_maplist": str(tmp_path / "missing-maplist.xlsx"),
    }

    # 1. Fresh plmsearch run with stale structure keys must not fail
    rt = _fake_runtime(tmp_path)
    custom_py = _make_executable(
        tmp_path / "alphagem_env" / "bin" / "python", "#!/bin/sh\nexit 0\n"
    )
    c1_dir = tmp_path / "case_plmsearch_stale"
    c1_dir.mkdir()
    cfg = _write(
        c1_dir,
        proteome_sha=True,
        runtime_root=str(rt),
        extra={**stale, "alphagem_python": str(custom_py)},
    )
    cfg_data = json.loads(cfg.read_text(encoding="utf-8"))
    cfg_data.pop("withgaps_model", None)
    cfg.write_text(json.dumps(cfg_data), encoding="utf-8")

    calls: list[list[str]] = []

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        name = cmd[cmd.index("--name") + 1]
        yml = Path(kw["cwd"]) / "working" / name / f"{name}-GEM_withgaps.yml"
        yml.parent.mkdir(parents=True, exist_ok=True)
        model, _ = _tiny_model_pair()
        import cobra

        cobra.io.save_yaml_model(model, str(yml))
        return MagicMock(returncode=0)

    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        stack.enter_context(
            patch.object(pipeline.subprocess, "run", side_effect=fake_run)
        )
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, tmp_path / "out_plmsearch_stale",
                localization_order_override="pre-gapfill",
            )
    assert calls
    assert calls[0][calls[0].index("--mode") + 1] == "plmsearch"
    for forbidden in ("--prottrans-model-path", "--maplist", "--structure"):
        assert forbidden not in calls[0]

    # 2. Existing-GEM localization_order='skip' run must not validate the
    # stale structure paths either.
    cfg = json.loads((synthetic_screening_case / "config.json").read_text())
    cfg["localization_order"] = "skip"
    cfg.pop("target_module", None)
    cfg["target_reaction"] = "EX_lac_e"
    cfg.update(stale)
    cfg_path = synthetic_screening_case / "stale_structure_skip.json"
    cfg_path.write_text(json.dumps(cfg, indent=2))
    out_dir = tmp_path / "stale_skip_out"
    result = pipeline.run_pipeline(cfg_path, out_dir)
    assert result["paths"]["screening"]["target_reaction"] == "EX_lac_e"


def test_eggnog_constants_mirror_paper_contracts():
    """fetcher.EGGNOG and fetcher.EGGNOG_STATIC_DATA must match
    third_party/alphagem/src/paper_tool_contracts.py exactly."""
    spec = importlib.util.spec_from_file_location(
        "paper_tool_contracts",
        Path(fetcher.__file__).resolve().parents[1]
        / "third_party"
        / "alphagem"
        / "src"
        / "paper_tool_contracts.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    db = module.PAPER_EGGNOG_DATABASE_CONTRACT
    assert fetcher.EGGNOG["db_bytes"] == db["bytes"]
    assert fetcher.EGGNOG["db_sha256"] == db["sha256"]
    assert fetcher.EGGNOG["db_gz_url"] == db["official_download_url"]

    fungi = module.PAPER_EGGNOG_TAXON_DATABASE_CONTRACTS["yeast"]
    assert fetcher.EGGNOG["fungi_dmnd_bytes"] == fungi["bytes"]
    assert fetcher.EGGNOG["fungi_dmnd_sha256"] == fungi["sha256"]
    assert fetcher.EGGNOG["proteomes_url"] == fungi["source_urls"]["proteomes"]
    assert fetcher.EGGNOG["taxid_url"] == fungi["source_urls"]["taxid_info"]
    assert (
        fetcher.EGGNOG["proteomes_sha256"]
        == fungi["source_sha256"]["proteomes"]
    )
    assert (
        fetcher.EGGNOG["taxid_sha256"]
        == fungi["source_sha256"]["taxid_info"]
    )
    assert (
        fetcher.EGGNOG["selected_fasta_sha256"]
        == fungi["source_sha256"]["selected_uncompressed_fasta"]
    )

    bundle = module.PAPER_TOOL_BUNDLE_CONTRACTS["eggnog_static_data"]
    assert (
        tuple(fetcher.EGGNOG_STATIC_DATA["excluded_relative_roots"])
        == tuple(bundle["excluded_relative_roots"])
    )
    for key in ("path", "file_count", "bytes", "sha256"):
        assert fetcher.EGGNOG_STATIC_DATA[key] == bundle[key], (
            f"EGGNOG_STATIC_DATA[{key!r}] deviates from paper contract"
        )
