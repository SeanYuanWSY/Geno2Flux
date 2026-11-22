"""Tests for knockout prediction module."""

import csv
from types import SimpleNamespace
from unittest.mock import MagicMock
import cobra
import pytest
from pathlib import Path

from geno2flux.knockout import (
    run_single_reaction_ko,
    run_single_gene_ko,
    run_family_ko,
    run_ko_suite,
    is_pathway_reaction,
    _loss,
)


@pytest.fixture
def toy_model():
    """Minimal toy model: A -> B -> C -> product."""
    m = cobra.Model("toy")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    c = cobra.Metabolite("c_c", compartment="c")
    p = cobra.Metabolite("p_c", compartment="c")
    m.add_metabolites([a, b, c, p])
    r1 = cobra.Reaction("R1")
    r1.add_metabolites({a: -1, b: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("R2")
    r2.add_metabolites({b: -1, c: 1})
    r2.gene_reaction_rule = "g2"
    r3 = cobra.Reaction("R3")
    r3.add_metabolites({c: -1, p: 1})
    r3.gene_reaction_rule = "g3"
    dm = cobra.Reaction("DM_p")
    dm.add_metabolites({p: -1})
    ex = cobra.Reaction("EX_a")
    ex.add_metabolites({a: -1})
    ex.bounds = (-10, 10)
    m.add_reactions([r1, r2, r3, dm, ex])
    m.objective = "DM_p"
    return m


class TestLoss:
    def test_optimal(self):
        assert _loss(0.5, 1.0, "optimal") == pytest.approx(0.5)

    def test_infeasible(self):
        assert _loss(None, 1.0, "infeasible") == 1.0

    def test_zero_obj(self):
        assert _loss(0.0, 1.0, "optimal") == 1.0

    def test_tiny_obj(self):
        assert _loss(1e-13, 1.0, "optimal") == 1.0


class TestSingleReactionKO:
    def test_essential_reaction(self, toy_model):
        opt, nonopt = run_single_reaction_ko(toy_model, "DM_p", Path("/tmp"))
        r1 = [r for r in opt if r["reaction_id"] == "R1"][0]
        assert r1["loss"] == 1.0
        assert r1["status"] == "optimal"
        assert all(r["status"] == "optimal" for r in opt)
        assert all(r["status"] != "optimal" for r in nonopt)

    def test_ranking_order(self, toy_model):
        opt, _ = run_single_reaction_ko(toy_model, "DM_p", Path("/tmp"))
        losses = [r["loss"] for r in opt]
        assert losses == sorted(losses, reverse=True)

    def test_nonoptimal_excluded_from_rank(self, toy_model):
        """Non-optimal rows must not appear in the ranked table."""
        # Gene KO: g3 is essential for R3→DM_p
        opt, nonopt = run_single_gene_ko(toy_model, "DM_p", Path("/tmp"))
        g3_opt = [r for r in opt if r["gene_id"] == "g3"]
        g3_nonopt = [r for r in nonopt if r["gene_id"] == "g3"]
        # g3 KO should give loss=1.0 (either optimal-with-0 or nonoptimal)
        if g3_opt:
            assert g3_opt[0]["loss"] == 1.0
        if g3_nonopt:
            assert g3_nonopt[0]["status"] != "optimal"
        # g3 must be in exactly one list
        assert len(g3_opt) + len(g3_nonopt) == 1


class TestSingleGeneKO:
    def test_essential_gene(self, toy_model):
        opt, _ = run_single_gene_ko(toy_model, "DM_p", Path("/tmp"))
        g1 = [r for r in opt if r["gene_id"] == "g1"][0]
        assert g1["loss"] == 1.0
        assert all(r["status"] == "optimal" for r in opt)


class TestFamilyKO:
    def test_family_ko(self, toy_model):
        families = {"FAM1": ["g1", "g2"], "FAM2": ["g3"]}
        opt, _ = run_family_ko(toy_model, "DM_p", families, Path("/tmp"))
        assert len(opt) == 2
        assert all(r["loss"] == 1.0 for r in opt)
        assert all(r["status"] == "optimal" for r in opt)


class TestPathwayFilter:
    def test_boundary(self, toy_model):
        ex = toy_model.reactions.get_by_id("EX_a")
        assert is_pathway_reaction(ex)  # boundary

    def test_internal_not_filtered(self, toy_model):
        r1 = toy_model.reactions.get_by_id("R1")
        assert not is_pathway_reaction(r1)


class TestFamilyKOConsistency:
    """Family KO input must match default tier exactly (no cross-tier mixing)."""

    def test_family_ids_match_default_tier(self):
        default_fams = [
            {"family_id": "FAM00001", "members": ["g1", "g2"]},
            {"family_id": "FAM00002", "members": ["g3"]},
        ]
        strict_fams = [{"family_id": "FAM00001", "members": ["g1"]}]
        loose_fams = [{"family_id": "FAM00001", "members": ["g1", "g2", "g4"]}]
        all_families = {"default": default_fams, "strict": strict_fams, "loose": loose_fams}

        ko_families = {f["family_id"]: f["members"] for f in all_families["default"]}
        for fam in all_families["default"]:
            fid = fam["family_id"]
            assert fid in ko_families
            assert sorted(ko_families[fid]) == sorted(fam["members"])
        assert len(ko_families) == len(default_fams)
        assert sorted(ko_families["FAM00001"]) == ["g1", "g2"]


def _run_ko(kind, model, out_dir):
    if kind == "reaction":
        return run_single_reaction_ko(model, "DM_p", out_dir)
    if kind == "gene":
        return run_single_gene_ko(model, "DM_p", out_dir)
    if kind == "family":
        return run_family_ko(model, "DM_p", {"FAM1": ["g1", "g2"]}, out_dir)
    return run_ko_suite(model, "DM_p", out_dir=out_dir)


def _model_state(model):
    # Compare the linear objective and its constant, ignoring zero terms
    # introduced by the solver when a model context restores the objective.
    objective_coefficients = {}
    for term, coefficient in model.objective.expression.as_coefficients_dict().items():
        value = float(coefficient)
        if value != 0.0:
            key = ("variable", term.name) if hasattr(term, "name") else ("constant",)
            objective_coefficients[key] = value
    return (
        objective_coefficients,
        model.objective_direction,
        [(r.id, r.bounds) for r in model.reactions],
        [(g.id, g.functional) for g in model.genes],
    )


@pytest.mark.parametrize("kind", ["reaction", "gene", "family", "suite"])
def test_knockout_restores_caller_objective_bounds_and_genes(kind, toy_model, tmp_path):
    toy_model.objective = "R2"
    toy_model.objective_direction = "min"
    before = _model_state(toy_model)
    result = _run_ko(kind, toy_model, tmp_path)
    assert _model_state(toy_model) == before
    if kind == "suite":
        assert result["base_objective"] == pytest.approx(10.0)


@pytest.mark.parametrize("kind", ["reaction", "gene", "family", "suite"])
@pytest.mark.parametrize(
    "status,objective",
    [
        ("optimal", 0.0),
        ("optimal", -1.0),
        ("optimal", None),
        ("optimal", float("nan")),
        ("optimal", float("inf")),
        ("infeasible", None),
        ("time_limit", 10.0),
    ],
)
def test_undefined_baseline_cannot_produce_knockout_rankings(
    kind, status, objective, toy_model, tmp_path, monkeypatch
):
    toy_model.objective = "R2"
    toy_model.objective_direction = "min"
    before = _model_state(toy_model)
    optimize = MagicMock(return_value=SimpleNamespace(status=status, objective_value=objective))
    monkeypatch.setattr(toy_model, "optimize", optimize)
    with pytest.raises(ValueError, match="baseline"):
        _run_ko(kind, toy_model, tmp_path / "ko")
    assert optimize.call_count == 1
    assert _model_state(toy_model) == before
    assert not (tmp_path / "ko").exists()


@pytest.mark.parametrize("kind", ["reaction", "gene", "family", "suite"])
def test_knockout_solve_exception_restores_caller_state(
    kind, toy_model, tmp_path, monkeypatch
):
    toy_model.objective = "R2"
    toy_model.objective_direction = "min"
    before = _model_state(toy_model)
    baseline = SimpleNamespace(status="optimal", objective_value=10.0)
    baseline_calls = 2 if kind == "suite" else 1
    optimize = MagicMock(side_effect=[baseline] * baseline_calls + [RuntimeError("solve failed")])
    monkeypatch.setattr(toy_model, "optimize", optimize)
    with pytest.raises(RuntimeError, match="solve failed"):
        _run_ko(kind, toy_model, tmp_path)
    assert _model_state(toy_model) == before


def test_or_rule_requires_simultaneous_family_knockout(toy_model, tmp_path):
    toy_model.reactions.get_by_id("R1").gene_reaction_rule = "g1 or g1_copy"
    opt, nonopt = run_single_gene_ko(toy_model, "DM_p", tmp_path)
    rows = {r["gene_id"]: r for r in opt}
    assert not nonopt
    assert rows["g1"]["loss"] == pytest.approx(0.0)
    assert rows["g1_copy"]["loss"] == pytest.approx(0.0)
    family_opt, family_nonopt = run_family_ko(
        toy_model, "DM_p", {"FAM_OR": ["g1", "g1_copy"]}, tmp_path
    )
    assert not family_nonopt
    assert family_opt[0]["loss"] == pytest.approx(1.0)
    assert all(g.functional for g in toy_model.genes)
    assert toy_model.reactions.get_by_id("R1").bounds == (0.0, 1000.0)


def test_nonoptimal_solve_is_supplement_only(toy_model, tmp_path, monkeypatch):
    baseline = SimpleNamespace(status="optimal", objective_value=10.0)
    optimize = MagicMock(side_effect=[
        baseline,
        SimpleNamespace(status="time_limit", objective_value=9.0),
        SimpleNamespace(status="optimal", objective_value=0.0),
        SimpleNamespace(status="optimal", objective_value=5.0),
    ])
    monkeypatch.setattr(toy_model, "optimize", optimize)
    opt, nonopt = run_single_gene_ko(toy_model, "DM_p", tmp_path)
    assert [r["ko_rank"] for r in opt] == [1, 2]
    assert [r["loss"] for r in opt] == [1.0, 0.5]
    assert len(nonopt) == 1
    assert nonopt[0]["status"] == "time_limit"
    assert nonopt[0]["obj"] is None
    assert "ko_rank" not in nonopt[0]


@pytest.mark.parametrize("families", [None, {}])
def test_suite_writes_empty_family_tables(toy_model, tmp_path, families):
    result = run_ko_suite(toy_model, "DM_p", families=families, out_dir=tmp_path)
    assert result["family_ko"] == 0
    assert result["family_ko_nonoptimal"] == 0
    for name in ("family_ko.tsv", "family_ko_nonoptimal.tsv"):
        with (tmp_path / name).open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            assert "family_id" in reader.fieldnames
            assert list(reader) == []


@pytest.mark.parametrize("base", [0.0, -1.0, float("nan"), float("inf")])
def test_loss_rejects_undefined_denominator(base):
    with pytest.raises(ValueError, match="baseline"):
        _loss(1.0, base, "optimal")


@pytest.mark.parametrize("objective", [float("nan"), float("inf")])
def test_loss_rejects_nonfinite_optimal_objective(objective):
    with pytest.raises(ValueError, match="knockout objective"):
        _loss(objective, 1.0, "optimal")


def test_loss_keeps_unclamped_improvement():
    assert _loss(1.5, 1.0, "optimal") == pytest.approx(-0.5)
