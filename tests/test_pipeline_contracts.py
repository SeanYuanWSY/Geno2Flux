"""General contract tests for reuse identity and runtime roots.

Real boundaries only — wrong-file reuse, per-file directory verification,
read-only external runtimes, and compartment semantics on synthetic
models. No source-string or mock-echo assertions.
"""

from __future__ import annotations

import hashlib
import json
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from geno2flux import pipeline
from geno2flux.knockout import is_pathway_reaction


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _inputs(tmp: Path) -> Path:
    d = tmp / "inputs"
    d.mkdir(exist_ok=True)
    return d


def _write(
    tmp: Path,
    *,
    extra: dict | None = None,
    proteome_sha: str | None | bool = None,
    withgaps_sha: str | None | bool = None,
    runtime_root: str | None = None,
    loc_sha: str | None | bool = True,
) -> Path:
    d = _inputs(tmp)
    proteome = d / "proteome.fasta"
    proteome.write_text(">p1\nMKTVL")
    withgaps = d / "withgaps.yml"
    import cobra
    wg_model, _ = _tiny_model_pair()
    cobra.io.save_yaml_model(wg_model, str(withgaps))
    loc = d / "loc.tsv"
    loc.write_text("protein_id\tlocalizations\np1\tNucleus\n")
    for n in ("ref.fasta", "ref.xlsx"):
        (d / n).write_text("x")
    (d / "self_hits.tsv").write_text(
        "p1\tp1\t100.0\t5\t1e-10\t100.0\t5\t5\n"
    )
    if proteome_sha is True:
        proteome_sha = _sha(proteome)
    if withgaps_sha is True:
        withgaps_sha = _sha(withgaps)
    if loc_sha is True:
        loc_sha = _sha(loc)
    cfg = {
        "name": "fungal_annevo_v1",
        "refname": "yeast",
        "solver": "highs",
        "proteome_fasta": str(proteome),
        "withgaps_model": str(withgaps),
        "localization_tsv": str(loc),
        "localization_sha256": loc_sha if isinstance(loc_sha, str) else None,
        "reference_fasta": str(d / "ref.fasta"),
        "reference_xlsx": str(d / "ref.xlsx"),
        "biomass": "R_target",
        "diamond_self_hits": str(d / "self_hits.tsv"),
        "diamond_self_hits_sha256": _sha(d / "self_hits.tsv"),
    }
    if proteome_sha is not None:
        cfg["proteome_sha256"] = proteome_sha
    if withgaps_sha is not None:
        cfg["withgaps_sha256"] = withgaps_sha
    if runtime_root is not None:
        cfg["alphagem_runtime_root"] = runtime_root
    cfg.update(extra or {})
    cp = tmp / "config.json"
    cp.write_text(json.dumps(cfg))
    return cp


class ExternalBoundaryStop(RuntimeError):
    """Named stop at a heavy external-call boundary this env can't honour."""


def _tiny_model_pair():
    import cobra

    def build(tag):
        m = cobra.Model(f"tiny_{tag}")
        met = cobra.Metabolite(f"a_{tag}_c", compartment="c")
        rx = cobra.Reaction(f"R_{tag}", lower_bound=0, upper_bound=1000)
        rx.add_metabolites({met: 1})
        m.add_reactions([rx])
        return m

    return build("target"), build("template")


def _light_stubs():
    """Honest boundary set.

    Everything before the boundary runs REAL: identity gates, runtime
    view construction, real file I/O, localization table load, step-4
    compartment refinement on real tiny models. The heavy external
    consumers this environment cannot honour raise ExternalBoundaryStop
    at their named boundary; tests assert exactly that sentinel — no
    no-op writers, no swallowed exceptions, no mock-echo success.
    """
    return [
        patch.object(pipeline, "_use_default_solver", MagicMock()),
        patch.object(
            pipeline, "_prepare_alphagem_runtime",
            side_effect=_prepare_synthetic_runtime,
        ),
        patch.object(pipeline, "_ensure_data_available", MagicMock()),
        patch.object(
            pipeline, "_run_gapfill",
            MagicMock(side_effect=ExternalBoundaryStop(
                "step5 external AlphaGEM/HiGHS gap-fill"
            )),
        ),
    ]


def _prepare_synthetic_runtime(root, manifest):
    """Vendored-path runtime stub: REAL tiny template where step 4 reads."""
    import cobra
    root.mkdir(parents=True, exist_ok=True)
    (root / "models").mkdir(exist_ok=True)
    tpl, _ = _tiny_model_pair()
    cobra.io.write_sbml_model(tpl, str(root / "models" / "yeast-GEM.xml"))
    (root / "models" / "reaction_pool_final.json").write_text("{}")


# ---------------------------------------------------------------------------
# Reuse identity gate
# ---------------------------------------------------------------------------


def test_wrong_file_reuse_fails_before_heavy_work(tmp_path: Path):
    """Expected hash of DIFFERENT content -> SHA mismatch, nothing runs."""
    other = tmp_path / "other.bin"
    other.write_bytes(b"not the proteome")
    cfg = _write(tmp_path, proteome_sha=_sha(other))
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        stack.enter_context(
            patch.object(
                pipeline, "_record_file",
                MagicMock(side_effect=lambda p: _record_file_stub()),
            )
        )
        with pytest.raises(pipeline.PipelineError, match="SHA256 mismatch"):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )
    # failed BEFORE any heavy work: no runtime copy was even prepared
    assert not (tmp_path / "out" / "alphagem_runtime").exists()


def test_missing_expected_hash_fails_by_default(tmp_path: Path):
    cfg = _write(tmp_path)  # no proteome_sha256
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(
            pipeline.PipelineError, match="proteome_fasta.*proteome_sha256"
        ):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )
    assert not (tmp_path / "out" / "alphagem_runtime").exists()


def test_opted_in_reuse_records_unverified_in_manifest(tmp_path: Path):
    cfg = _write(tmp_path, extra={"allow_unverified_reuse": True})
    out = tmp_path / "out"
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, out, localization_order_override="pre-gapfill"
            )
    m = json.loads((out / "pipeline_manifest.json").read_text())
    step0 = next(s for s in m["steps"] if s["step"] == "step0_proteome")
    assert step0["inputs"]["identity"]["identity"] == "unverified_opted_in"
    assert step0["inputs"]["identity"]["expected_sha256"] is None


def test_verified_reuse_records_verified(tmp_path: Path):
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        withgaps_sha=True,
    )
    out = tmp_path / "out"
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, out, localization_order_override="pre-gapfill"
            )
    m = json.loads((out / "pipeline_manifest.json").read_text())
    step0 = next(s for s in m["steps"] if s["step"] == "step0_proteome")
    assert step0["inputs"]["identity"]["identity"] == "verified"


def test_pre_gapfill_dir_requires_per_file_map(tmp_path: Path):
    d = _inputs(tmp_path)
    pdir = tmp_path / "pre"
    pdir.mkdir()
    (pdir / "fungal_annevo_v1-GEM_withgaps.yml").write_text("id: x\n")
    (pdir / "extra.json").write_text("{}\n")
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        extra={
            "withgaps_model": None,
            "pre_gapfill_dir": str(pdir),
            "pre_gapfill_dir_sha256_map": {
                "fungal_annevo_v1-GEM_withgaps.yml": _sha(
                    pdir / "fungal_annevo_v1-GEM_withgaps.yml"
                ),
                # extra.json deliberately NOT declared
            },
        },
    )
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(
            pipeline.PipelineError, match="extra.json.*sha256_map"
        ):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )


def test_pre_gapfill_dir_declared_but_absent_fails(tmp_path: Path):
    pdir = tmp_path / "pre"
    pdir.mkdir()
    (pdir / "fungal_annevo_v1-GEM_withgaps.yml").write_text("id: x\n")
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        extra={
            "withgaps_model": None,
            "pre_gapfill_dir": str(pdir),
            "pre_gapfill_dir_sha256_map": {
                "fungal_annevo_v1-GEM_withgaps.yml": _sha(
                    pdir / "fungal_annevo_v1-GEM_withgaps.yml"
                ),
                "ghost.txt": "0" * 64,
            },
        },
    )
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(pipeline.PipelineError, match="ghost.txt"):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )


def test_pre_gapfill_dir_wrong_file_hash_fails(tmp_path: Path):
    pdir = tmp_path / "pre"
    pdir.mkdir()
    (pdir / "fungal_annevo_v1-GEM_withgaps.yml").write_text("id: x\n")
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        extra={
            "withgaps_model": None,
            "pre_gapfill_dir": str(pdir),
            "pre_gapfill_dir_sha256_map": {
                "fungal_annevo_v1-GEM_withgaps.yml": "f" * 64,  # wrong
            },
        },
    )
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(pipeline.PipelineError, match="SHA256 mismatch"):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )


# ---------------------------------------------------------------------------
# External runtime root (read-only contract)
# ---------------------------------------------------------------------------


def _fake_runtime(tmp: Path) -> Path:
    rt = tmp / "prepared_runtime"
    (rt / "src").mkdir(parents=True)
    (rt / "models").mkdir()
    (rt / "data_available").mkdir()
    (rt / "AlphaGEM.py").write_text("# prepared\n")
    (rt / "tools").mkdir()
    (rt / "tools" / "rhea_marker.txt").write_text("rhea\n")
    (rt / "SOURCE_IDENTITY.json").write_text("{}\n")
    (rt / "src" / "probe_mod.py").write_text("VALUE = 1\n")
    import cobra
    tpl, _ = _tiny_model_pair()
    cobra.io.write_sbml_model(tpl, str(rt / "models" / "yeast-GEM.xml"))
    (rt / "models" / "reaction_pool_final.json").write_text("{}")
    (rt / "data_available" / "yeast.fasta").write_text(">y\nMK")
    (rt / "data_available" / "yeast.xlsx").write_text("xlsx")
    return rt


def _tree_state(root: Path) -> dict[str, tuple[float, int, str]]:
    out: dict[str, tuple[float, int, str]] = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            st = p.stat()
            out[str(p.relative_to(root))] = (
                st.st_mtime_ns, st.st_size, _sha(p)
            )
    return out


def test_external_runtime_missing_entries_fails_with_list(tmp_path: Path):
    rt = _fake_runtime(tmp_path)
    (rt / "models" / "reaction_pool_final.json").unlink()
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        withgaps_sha=True,
        runtime_root=str(rt),
    )
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(
            pipeline.PipelineError, match="reaction_pool_final"
        ):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )


def test_external_runtime_missing_references_fail_read_only(tmp_path: Path):
    rt = _fake_runtime(tmp_path)
    (rt / "data_available" / "yeast.xlsx").unlink()
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        withgaps_sha=True,
        runtime_root=str(rt),
    )
    before = _tree_state(rt)
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(
            pipeline.PipelineError, match="yeast.xlsx"
        ):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )
    assert _tree_state(rt) == before  # nothing written, nothing healed
    assert not (tmp_path / "out" / "alphagem_runtime_view").exists()


def test_external_runtime_fresh_view_fully_writable(tmp_path):
    """Fresh steps 1-3 run against the VIEW, never the read-only tree.

    Real boundaries: the subprocess consumer's cwd is the view; the view
    is a COMPLETE per-run copy (entry file included, so the entry's own
    __file__ resolves REPO_ROOT inside the view); the original tree is
    byte-identical afterwards. Upstream tool EXECUTION is stubbed at the
    named boundary (recorded honestly in the manifest as not validated).
    """
    rt = _fake_runtime(tmp_path)
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        runtime_root=str(rt),
    )
    del_cfg = json.loads(cfg.read_text())
    del_cfg.pop("withgaps_model", None)  # force fresh steps 1-3
    cfg.write_text(json.dumps(del_cfg))
    out = tmp_path / "out"
    calls = []

    def fake_run(cmd, **kw):
        calls.append((list(cmd), kw.get("cwd")))
        name = cmd[cmd.index("--name") + 1]
        yml = Path(kw["cwd"]) / "working" / name / f"{name}-GEM_withgaps.yml"
        yml.parent.mkdir(parents=True, exist_ok=True)
        model, _ = _tiny_model_pair()
        import cobra
        cobra.io.save_yaml_model(model, str(yml))
        return MagicMock(returncode=0)

    before = _tree_state(rt)
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        stack.enter_context(
            patch.object(pipeline.subprocess, "run", side_effect=fake_run)
        )
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, out, localization_order_override="pre-gapfill"
            )
    view = out / "alphagem_runtime_view"
    # the upstream consumer ran INSIDE the view
    assert calls, "steps 1-3 subprocess must run"
    cmd0, cwd0 = calls[0]
    assert Path(cwd0) == view
    assert str(view / "AlphaGEM.py") in cmd0
    # fresh mode: COMPLETE per-run copy — the real entry resolves its
    # REPO_ROOT from __file__, so AlphaGEM.py must be a real copy whose
    # parent IS the view (a symlink would chdir/write back into the
    # original tree); tool siblings and source identity ride along
    assert (view / "models").is_dir() and not (view / "models").is_symlink()
    assert (view / "src").is_dir() and not (view / "src").is_symlink()
    assert (view / "AlphaGEM.py").is_file() and not (
        view / "AlphaGEM.py"
    ).is_symlink()
    assert (view / "AlphaGEM.py").resolve().parent == view.resolve()
    assert (view / "tools" / "rhea_marker.txt").is_file()
    assert (view / "SOURCE_IDENTITY.json").is_file()
    assert _tree_state(rt) == before
    m = json.loads((out / "pipeline_manifest.json").read_text())
    rp = next(s for s in m["steps"] if s["step"] == "runtime_prepare")
    det = rp["detail"] if "detail" in rp else rp["outputs"]
    assert det["fresh_steps_1_3"].startswith("view is a COMPLETE")
    assert "not validated" in det["fresh_steps_1_3"]


def test_external_runtime_used_read_only_with_isolated_working(tmp_path):
    rt = _fake_runtime(tmp_path)
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        withgaps_sha=True,
        runtime_root=str(rt),
    )
    out = tmp_path / "out"
    before = _tree_state(rt)
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, out, localization_order_override="pre-gapfill"
            )
    # runtime tree untouched
    assert _tree_state(rt) == before
    # per-run VIEW under the output carries all relative writes
    view = out / "alphagem_runtime_view"
    assert (view / "working" / "fungal_annevo_v1").is_dir()
    assert (view / "data_available" / "yeast.fasta").is_file()
    assert not (rt / "working").exists()
    # models/AlphaGEM.py stay symlinked; src is a per-run COPY because
    # in-process imports would write __pycache__ through a symlink
    assert (view / "models").is_symlink()
    assert not (view / "src").is_symlink()
    # REAL import through the view: bytecode lands in the view copy,
    # never in the original tree
    import sys as _sys
    import importlib
    _sys.path.insert(0, str(view / "src"))
    try:
        mod = importlib.import_module("probe_mod")  # REAL import
        assert mod.VALUE == 1
    finally:
        _sys.path.pop(0)
        _sys.modules.pop("probe_mod", None)
        _sys.modules.pop("src", None)  # namespace pkg would shadow others
    # bytecode cache landed in the per-run COPY, never in the original
    assert any(
        p.suffix == ".pyc" for p in (view / "src" / "__pycache__").glob("probe_mod.*")
    )
    assert not (rt / "src" / "__pycache__").exists()

    # Same-name package precedence: poison the cache with the VENDORED
    # bundle's src, then prefer_runtime must re-resolve into the view —
    # the real gap-fill import path cannot silently bypass the view.
    import sys as _s2
    import importlib as _il
    _path_snapshot = list(_s2.path)
    _src_mods_snapshot = {
        k: v for k, v in _s2.modules.items()
        if k == "src" or k.startswith("src.")
    }
    bundle_src = Path(pipeline.__file__).resolve().parents[2] / "third_party" / "alphagem"
    _s2.path.insert(0, str(bundle_src))
    try:
        _il.import_module("src")  # resolves to the vendored bundle
        vendored_pkg = _s2.modules["src"]
        assert not list(vendored_pkg.__path__)[0].startswith(str(view))
    finally:
        _s2.path.remove(str(bundle_src))
    pipeline._repo_imports(view, prefer_runtime=True)
    view_pkg = _il.import_module("src")  # re-import resolves fresh
    assert list(view_pkg.__path__)[0].startswith(str(view))
    # and the vendored bundle is no longer on sys.path ahead of the view
    assert str(bundle_src) not in _s2.path
    # restore EXACT sys.path/src-module state (prefer_runtime removed the
    # vendored bundle path; other suite tests rely on it)
    _s2.path[:] = _path_snapshot
    for k in [k for k in list(_s2.modules) if k == "src" or k.startswith("src.")]:
        if k not in _src_mods_snapshot:
            del _s2.modules[k]
    # relative-write containment: writing through the view never
    # touches the original tree
    (view / "data_available" / "runtime_write_probe.txt").write_text("x")
    (view / "working" / "fungal_annevo_v1" / "probe.txt").write_text("x")
    assert not (rt / "data_available" / "runtime_write_probe.txt").exists()
    assert _tree_state(rt) == before
    m = json.loads((out / "pipeline_manifest.json").read_text())
    prep = next(
        s for s in m["steps"] if s["step"] == "runtime_prepare"
    )
    assert (
        prep["detail"]["runtime_source"]
        == "external_configured_runtime_read_only"
    )


# ---------------------------------------------------------------------------
# Malformed expected hashes (falsy/malformed can never count as verified)
# ---------------------------------------------------------------------------


def test_falsy_or_malformed_expected_hash_rejected(tmp_path):
    d = _inputs(tmp_path)
    for bad in ("", "0" * 63, "zz" * 32, 0, False):
        cfg = _write(tmp_path)
        c = json.loads(cfg.read_text())
        c["proteome_sha256"] = bad
        cfg.write_text(json.dumps(c))
        with pytest.raises(pipeline.PipelineError, match="malformed expected"):
            pipeline.run_pipeline(
                cfg, tmp_path / f"out_{abs(hash(str(bad)))}",
                localization_order_override="pre-gapfill",
            )


def test_pre_gapfill_dir_falsy_map_value_rejected(tmp_path):
    pdir = tmp_path / "pre"
    pdir.mkdir()
    (pdir / "fungal_annevo_v1-GEM_withgaps.yml").write_text("id: x\n")
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        extra={
            "withgaps_model": None,
            "pre_gapfill_dir": str(pdir),
            "pre_gapfill_dir_sha256_map": {
                "fungal_annevo_v1-GEM_withgaps.yml": "",
            },
        },
    )
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(
            pipeline.PipelineError, match="malformed expected hash"
        ):
            pipeline.run_pipeline(
                cfg, tmp_path / "out",
                localization_order_override="pre-gapfill",
            )


# ---------------------------------------------------------------------------
# Output placement guards
# ---------------------------------------------------------------------------


def test_output_inside_external_runtime_rejected_zero_writes(tmp_path):
    rt = _fake_runtime(tmp_path)
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path, proteome_sha=True, withgaps_sha=True, runtime_root=str(rt)
    )
    before = _tree_state(rt)
    for out_path in (rt / "out", rt):  # under it, and equal to it
        with pytest.raises(
            pipeline.PipelineError, match="inside the read-only external"
        ):
            pipeline.run_pipeline(
                cfg, out_path, localization_order_override="pre-gapfill"
            )
    assert _tree_state(rt) == before


def test_output_inside_runtime_via_symlink_rejected(tmp_path):
    rt = _fake_runtime(tmp_path)
    link = tmp_path / "link_to_rt"
    link.symlink_to(rt, target_is_directory=True)
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path, proteome_sha=True, withgaps_sha=True, runtime_root=str(rt)
    )
    before = _tree_state(rt)
    with pytest.raises(
        pipeline.PipelineError, match="inside the read-only external"
    ):
        pipeline.run_pipeline(
            cfg, link / "out", localization_order_override="pre-gapfill"
        )
    assert _tree_state(rt) == before


def test_configured_but_missing_withgaps_refuses_not_fresh(tmp_path):
    """External runtime + withgaps_model pointing at a MISSING file must
    fail before the view is built — never silently re-dispatch as fresh
    steps 1-3 (which would run the entry against symlinked resources)."""
    rt = _fake_runtime(tmp_path)
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path, proteome_sha=True, withgaps_sha=True, runtime_root=str(rt)
    )
    c = json.loads(cfg.read_text())
    c["withgaps_model"] = str(d / "does_not_exist.yml")
    cfg.write_text(json.dumps(c))
    out = tmp_path / "out"
    before = _tree_state(rt)
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with patch.object(
            pipeline.subprocess, "run",
            side_effect=AssertionError("no subprocess may run"),
        ):
            with pytest.raises(
                pipeline.PipelineError, match="refusing to guess"
            ):
                pipeline.run_pipeline(
                    cfg, out, localization_order_override="pre-gapfill"
                )
    assert _tree_state(rt) == before
    assert not (out / "alphagem_runtime_view").exists()


def test_configured_pre_gapfill_wrong_type_refuses(tmp_path):
    rt = _fake_runtime(tmp_path)
    d = _inputs(tmp_path)
    not_a_dir = d / "actually_a_file.tsv"
    not_a_dir.write_text("x\n")
    cfg = _write(
        tmp_path,
        proteome_sha=True,
        runtime_root=str(rt),
        extra={
            "withgaps_model": None,
            "pre_gapfill_dir": str(not_a_dir),
        },
    )
    out = tmp_path / "out"
    before = _tree_state(rt)
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with patch.object(
            pipeline.subprocess, "run",
            side_effect=AssertionError("no subprocess may run"),
        ):
            with pytest.raises(
                pipeline.PipelineError, match="not a directory.*refusing"
            ):
                pipeline.run_pipeline(
                    cfg, out, localization_order_override="pre-gapfill"
                )
    assert _tree_state(rt) == before
    assert not (out / "alphagem_runtime_view").exists()


def test_second_run_same_output_refused_before_heavy_work(tmp_path):
    rt = _fake_runtime(tmp_path)
    d = _inputs(tmp_path)
    cfg = _write(
        tmp_path, proteome_sha=True, withgaps_sha=True, runtime_root=str(rt)
    )
    out = tmp_path / "out"
    with ExitStack() as stack:
        for p in _light_stubs():
            stack.enter_context(p)
        with pytest.raises(ExternalBoundaryStop, match="step5"):
            pipeline.run_pipeline(
                cfg, out, localization_order_override="pre-gapfill"
            )
    with pytest.raises(pipeline.PipelineError, match="previous run"):
        pipeline.run_pipeline(
            cfg, out, localization_order_override="pre-gapfill"
        )


# ---------------------------------------------------------------------------
# Import-state isolation between runtime roots (external -> vendored)
# ---------------------------------------------------------------------------


def test_vendored_run_after_external_purges_view_modules(tmp_path):
    """Real module import: after an external-view run, the vendored mode
    must re-resolve src.* to the vendored bundle, not keep executing the
    stale view copy."""
    import sys
    import importlib

    rt = _fake_runtime(tmp_path)
    (rt / "src" / "__init__.py").write_text("")
    (rt / "src" / "gapfilling.py").write_text("RUNTIME_MARKER = 'view'\n")
    view = tmp_path / "out" / "alphagem_runtime_view"
    view.mkdir(parents=True)
    import shutil as _sh
    _sh.copytree(rt / "src", view / "src")

    snapshot_paths = list(sys.path)
    snapshot_src = {
        k: v for k, v in sys.modules.items()
        if k == "src" or k.startswith("src.")
    }
    try:
        # external run state: view's src imported and cached
        pipeline._repo_imports(view, prefer_runtime=True)
        mod_view = importlib.import_module("src.gapfilling")
        assert mod_view.RUNTIME_MARKER == "view"
        # vendored run: must purge the view residue and re-resolve into
        # the vendored bundle (REAL import, not a sys.path assertion)
        vendored = pipeline._bundle_root()
        pipeline._repo_imports(vendored, prefer_runtime=False)
        mod_vendored = importlib.import_module("src.gapfilling")
        origin = Path(mod_vendored.__file__).resolve()
        assert origin.is_relative_to(Path(vendored).resolve())
        assert not origin.is_relative_to(view.resolve())
    finally:
        sys.path[:] = snapshot_paths
        for k in [k for k in list(sys.modules)
                  if k == "src" or k.startswith("src.")]:
            if k not in snapshot_src:
                del sys.modules[k]


# ---------------------------------------------------------------------------
# Knockout compartment semantics (synthetic model)
# ---------------------------------------------------------------------------


def _model_with(rxn):
    import cobra

    m = cobra.Model("synthetic")
    m.add_reactions([rxn])
    return m


def test_compartment_attribute_not_id_suffix_drives_filter():
    """IDs may lie about compartments; Metabolite.compartment must win."""
    import cobra

    # Same "_c" id suffix, DIFFERENT explicit compartments -> multi-compartment
    a = cobra.Metabolite("same_c", compartment="c")
    b = cobra.Metabolite("also_c", compartment="m")
    rxn = cobra.Reaction("R1", lower_bound=-1000, upper_bound=1000)
    rxn.add_metabolites({a: -1, b: 1})
    rxn.subsystem = "other"
    _model_with(rxn)
    assert rxn not in rxn.model.boundary
    assert is_pathway_reaction(rxn) is True  # suffix heuristic missed this

    # Different suffixes, SAME compartment -> NOT multi-compartment
    c = cobra.Metabolite("x_c", compartment="c")
    d = cobra.Metabolite("y_e", compartment="c")
    rxn2 = cobra.Reaction("R2", lower_bound=-1000, upper_bound=1000)
    rxn2.add_metabolites({c: -1, d: 1})
    rxn2.subsystem = "other"
    _model_with(rxn2)
    assert is_pathway_reaction(rxn2) is False  # suffix heuristic lied here


@pytest.mark.parametrize("key", [
    "proteome_fasta", "localization_tsv", "gapfilled_model",
    "diamond_self_hits", "id_map", "essential_genes_tsv",
])
def test_configured_missing_file_fails_before_runtime(tmp_path, key):
    cfg = _write(tmp_path, proteome_sha=True, withgaps_sha=True,
                 extra={key: str(tmp_path / "missing")})
    with patch.object(pipeline, "_prepare_alphagem_runtime") as prepare:
        with pytest.raises(pipeline.PipelineError, match=f"{key}.*readable file"):
            pipeline.run_pipeline(cfg, tmp_path / "out")
    prepare.assert_not_called()


@pytest.mark.parametrize("value", ["false", "true", 1, None])
def test_unverified_reuse_requires_json_boolean(tmp_path, value):
    cfg = _write(tmp_path, extra={"allow_unverified_reuse": value})
    with patch.object(pipeline, "_prepare_alphagem_runtime") as prepare:
        with pytest.raises(pipeline.PipelineError, match="JSON boolean"):
            pipeline.run_pipeline(cfg, tmp_path / "out")
    prepare.assert_not_called()


@pytest.mark.parametrize("contents, threshold", [
    ("gene_id\tloss\nA\t0.9\n", 0.5),
    ("gene_id\tloss\tstatus\nA\t0.9\toptimal\n", 2),
])
def test_invalid_essential_config_fails_before_runtime(tmp_path, contents, threshold):
    essential = tmp_path / "growth.tsv"
    essential.write_text(contents)
    cfg = _write(tmp_path, proteome_sha=True, withgaps_sha=True,
                 extra={"essential_genes_tsv": str(essential),
                        "essential_ko_loss_threshold": threshold})
    with patch.object(pipeline, "_prepare_alphagem_runtime") as prepare:
        with pytest.raises(pipeline.PipelineError):
            pipeline.run_pipeline(cfg, tmp_path / "out")
    prepare.assert_not_called()
