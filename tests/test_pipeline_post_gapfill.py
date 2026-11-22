"""Behavior tests for localization_order=post-gapfill wiring.

These tests exercise ``run_pipeline`` end-to-end with all expensive externals
(AlphaGEM runtime, cobra I/O, gapfill MILP, refinement, FSEOF, paralog/KO)
mocked out. They assert:

1. ``post-gapfill`` *without* ``gapfilled_model``: ``_run_gapfill`` is invoked
   EXACTLY ONCE to build the localization base, with ``out_dir`` ending in
   ``_post_gapfill_base``, and BEFORE ``refine_compartments`` runs.
2. ``post-gapfill`` *with* ``gapfilled_model``: ``_run_gapfill`` is NEVER
   invoked (the archived gapfilled model is the localization base directly).
3. ``pre-gapfill``: ``refine_compartments`` runs BEFORE ``_run_gapfill``
   (order regression guard).

The call-order assertions use a shared recorder so we don't rely on code
inspection or keyword matching.
"""

from __future__ import annotations

import json
from contextlib import ExitStack
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from geno2flux import pipeline
from geno2flux.pipeline import run_pipeline


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fake_model() -> MagicMock:
    """Minimal cobra.Model stand-in: reactions + genes iterables."""
    m = MagicMock()
    m.reactions = []
    m.genes = []
    m.medium = {}
    return m


def _write_minimal_config(tmp_path: Path, **overrides: Any) -> Path:
    """Write a pipeline config pointing at stub inputs in tmp_path.

    Required keys the pipeline actually reads:
      proteome_fasta, withgaps_model, localization_tsv,
      reference_fasta, reference_xlsx (runtime bundling)
    """
    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    for name in (
        "proteome.fasta",
        "withgaps.yml",
        "loc.tsv",
        "ref.fasta",
        "ref.xlsx",
    ):
        (inputs / name).write_text("x")

    cfg: dict[str, Any] = {
        "name": "test",
        "refname": "yeast",
        "solver": "highs",
        "target_reaction": "TARGET",
        "allow_unverified_reuse": True,
        "proteome_fasta": str(inputs / "proteome.fasta"),
        "withgaps_model": str(inputs / "withgaps.yml"),
        "localization_tsv": str(inputs / "loc.tsv"),
        "reference_fasta": str(inputs / "ref.fasta"),
        "reference_xlsx": str(inputs / "ref.xlsx"),
    }
    cfg.update(overrides)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(cfg))
    return config_path


@pytest.fixture
def stub_pipeline(tmp_path: Path):
    """Patch every heavy external in ``pipeline``; yield the call recorder.

    Yields ``(calls, mock_run_gapfill, mock_refine)``. ``calls`` is a list of
    step-name strings appended in the order the mocked steps run, so tests can
    assert relative ordering without caring about absolute timing.
    """
    calls: list[str] = []

    fake_run_gapfill = MagicMock(
        side_effect=lambda **kw: (
            calls.append("gapfill"),
            Path(kw["out_dir"]) / "final.xml",
        )[1]
    )
    fake_refine = MagicMock(
        side_effect=lambda *a, **kw: (
            calls.append("refine"),
            {"pulled_reaction_count": 0, "relocated": []},
        )[1]
    )
    fake_load_models = MagicMock(return_value=(_fake_model(), _fake_model()))

    fake_fseof_result = {
        "reaction_rows": [],
        "maximum_product": 0.0,
        "raw_candidate_count": 0,
        "robust_candidate_count": 0,
        "scans": [],
    }

    patches = [
        patch.object(pipeline, "_prepare_alphagem_runtime", return_value=None),
        patch.object(pipeline, "_ensure_data_available", return_value=None),
        patch.object(pipeline, "_check_sha", return_value=None),
        patch.object(pipeline, "_use_default_solver", return_value=None),
        patch.object(pipeline, "load_localization_table", return_value={}),
        patch.object(pipeline, "load_models_for_refinement", fake_load_models),
        patch.object(pipeline, "refine_compartments", fake_refine),
        patch.object(pipeline, "_run_gapfill", fake_run_gapfill),
        patch.object(pipeline, "mount_target_module", return_value={"product_reaction": "TARGET"}),
        patch.object(pipeline, "_attach_solver", return_value=None),
        patch.object(
            pipeline, "run_published_fseof", return_value=fake_fseof_result
        ),
        patch.object(pipeline, "_load_protein_map", return_value=[]),
        patch.object(pipeline, "_gene_report", return_value=[]),
        patch.object(pipeline, "_write_tsv", return_value=None),
        patch.object(pipeline, "run_diamond_self_alignment", return_value=None),
        patch.object(pipeline, "load_self_hits", return_value=[]),
        patch.object(pipeline, "build_families", return_value=[]),
        patch.object(pipeline, "rank_paralog_families", return_value=[]),
        patch.object(pipeline, "write_family_tables", return_value=None),
        patch.object(
            pipeline,
            "_record_file",
            side_effect=lambda p: {
                "path": str(p),
                "bytes": 0,
                "sha256": "0" * 64,
            },
        ),
        patch.object(
            pipeline, "run_ko_suite", return_value={"top_reaction_loss": None}
        ),
        patch.object(
            pipeline.cobra.io,
            "read_sbml_model",
            side_effect=lambda p: _fake_model(),
        ),
        patch.object(pipeline.cobra.io, "write_sbml_model", return_value=None),
        patch.object(pipeline.cobra.io, "save_yaml_model", return_value=None),
        patch.object(pipeline.cobra.io, "save_json_model", return_value=None),
        patch.object(pipeline.shutil, "copy2", return_value=None),
        patch.object(pipeline, "save_json", return_value=None),
    ]

    with ExitStack() as stack:
        for p in patches:
            stack.enter_context(p)
        yield calls, fake_run_gapfill, fake_refine


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_post_gapfill_calls_gapfill_once_before_refine(
    stub_pipeline, tmp_path: Path
):
    """post-gapfill, no archived model: gapfill the raw withgaps first."""
    calls, fake_run_gapfill, fake_refine = stub_pipeline
    config_path = _write_minimal_config(tmp_path)
    out = tmp_path / "out"

    run_pipeline(
        config_path, out, localization_order_override="post-gapfill"
    )

    # Exactly one gapfill call.
    assert fake_run_gapfill.call_count == 1, (
        f"expected 1 _run_gapfill call, got {fake_run_gapfill.call_count}"
    )
    # It targeted the dedicated _post_gapfill_base directory.
    kwargs = fake_run_gapfill.call_args.kwargs
    assert str(kwargs["out_dir"]).endswith("_post_gapfill_base"), (
        f"_run_gapfill out_dir should end in _post_gapfill_base, "
        f"got {kwargs['out_dir']}"
    )
    # It ran BEFORE compartment refinement.
    assert calls == ["gapfill", "refine"], (
        f"expected gapfill->refine order, got {calls}"
    )


def test_post_gapfill_with_archived_model_skips_gapfill(
    stub_pipeline, tmp_path: Path
):
    """post-gapfill + gapfilled_model: never call _run_gapfill."""
    calls, fake_run_gapfill, fake_refine = stub_pipeline
    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    archived = inputs / "archived_gapfilled.xml"
    archived.write_text("x")
    config_path = _write_minimal_config(
        tmp_path, gapfilled_model=str(archived)
    )
    out = tmp_path / "out"

    run_pipeline(
        config_path, out, localization_order_override="post-gapfill"
    )

    assert fake_run_gapfill.call_count == 0, (
        f"expected 0 _run_gapfill calls, got {fake_run_gapfill.call_count}"
    )
    assert calls == ["refine"], (
        f"expected only refine to run, got {calls}"
    )


def test_pre_gapfill_refines_before_gapfill(
    stub_pipeline, tmp_path: Path
):
    """pre-gapfill order: refine runs first, then gapfill (regression)."""
    calls, fake_run_gapfill, fake_refine = stub_pipeline
    config_path = _write_minimal_config(tmp_path)
    out = tmp_path / "out"

    run_pipeline(
        config_path, out, localization_order_override="pre-gapfill"
    )

    assert fake_run_gapfill.call_count == 1, (
        f"expected 1 _run_gapfill call, got {fake_run_gapfill.call_count}"
    )
    kwargs = fake_run_gapfill.call_args.kwargs
    assert str(kwargs["out_dir"]).endswith("pre-gapfill"), (
        f"pre-gapfill _run_gapfill should write to path_dir, "
        f"got {kwargs['out_dir']}"
    )
    assert calls == ["refine", "gapfill"], (
        f"expected refine->gapfill order, got {calls}"
    )


def test_screening_keeps_gapfilled_model_medium(stub_pipeline, tmp_path):
    """Reusing a model must not replace its medium with the template's."""
    import cobra

    def nutrient_model(bound):
        model = cobra.Model("medium_contract")
        nutrient = cobra.Metabolite("nutrient_e", compartment="e")
        model.add_metabolites([nutrient])
        model.add_boundary(nutrient, type="exchange", reaction_id="EX_nutrient")
        model.medium = {"EX_nutrient": bound}
        return model

    archived = tmp_path / "archived.xml"
    archived.write_text("archived input boundary")
    config_path = _write_minimal_config(tmp_path, gapfilled_model=str(archived))
    candidate = nutrient_model(17)
    reference = nutrient_model(1)

    class ScreeningBoundary(RuntimeError):
        pass

    def observe_medium(model, *args, **kwargs):
        assert model.medium == {"EX_nutrient": 17}
        raise ScreeningBoundary("FSEOF external boundary")

    with patch.object(
        pipeline.cobra.io, "read_sbml_model",
        side_effect=lambda path: reference if "/models/" in path else candidate,
    ), patch.object(pipeline, "run_published_fseof", side_effect=observe_medium):
        with pytest.raises(ScreeningBoundary, match="FSEOF"):
            run_pipeline(config_path, tmp_path / "out",
                         localization_order_override="post-gapfill")


def test_checkpoint_hook_restored_after_gapfill_failure(tmp_path):
    import sys
    from types import ModuleType
    import cobra

    runtime = tmp_path / "runtime"
    runtime.mkdir()
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    model = cobra.Model("checkpoint")
    original = MagicMock(return_value={})
    module = ModuleType("src.gapfilling")
    module._prune_unsupported_reactions = original

    def fail_gapfill(*args, **kwargs):
        module._prune_unsupported_reactions(model, [])
        raise RuntimeError("gapfill failed")

    module.gapfill = fail_gapfill
    package = ModuleType("src")
    package.gapfilling = module
    with patch.dict(sys.modules, {"src": package, "src.gapfilling": module}), \
         patch.object(pipeline, "_repo_imports"), \
         patch.object(pipeline, "_use_default_solver"):
        with pytest.raises(RuntimeError, match="gapfill failed"):
            pipeline._run_gapfill(
                name="toy", refname="yeast", runtime_root=runtime,
                working=runtime / "working" / "toy", solver="highs",
                cfg={"audit_checkpoint_dir": "checkpoints"}, cdir=config_dir,
                out_dir=tmp_path / "out", tag="pre",
            )
    assert module._prune_unsupported_reactions is original
    assert (config_dir / "checkpoints" / "toy_model2_pre_audit.xml").is_file()
    assert not (runtime / "checkpoints").exists()
