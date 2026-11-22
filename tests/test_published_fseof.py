"""Behavioral checks of the source-defined ordinary-GEM FSEOF port."""

import json
from pathlib import Path
import sys

import cobra
import pytest

from geno2flux.published_fseof import (
    FSEOFConfig,
    PublishedFSEOFError,
    classify_fseof_trajectory,
    run_published_fseof,
    validate_flux_witness,
)


def metabolic_toy(reverse=False):
    model = cobra.Model("fseof_toy")
    model.solver = "glpk"
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    uptake = cobra.Reaction("UP")
    uptake.add_metabolites({a: 1})
    uptake.bounds = (0, 10)
    biomass = cobra.Reaction("BIO")
    biomass.add_metabolites({a: -1})
    biomass.bounds = (0, 100)
    conversion = cobra.Reaction("REV" if reverse else "CONVERT")
    conversion.add_metabolites({a: 1, b: -1} if reverse else {a: -1, b: 1})
    conversion.bounds = (-100, 100) if reverse else (0, 100)
    conversion.gene_reaction_rule = "g1 or g2"
    product = cobra.Reaction("PRODUCT")
    product.add_metabolites({b: -1})
    product.bounds = (0, 100)
    model.add_reactions([uptake, biomass, conversion, product])
    model.objective = "BIO"
    return model


@pytest.mark.parametrize(
    "values,raw,robust,direction",
    [
        ([1, 2, 3], True, True, 1),
        ([-1, -2, -3], True, True, -1),
        ([0, 1, 2], False, False, 0),
        ([1, 1, 2], False, False, 1),
        ([1, 3, 2], False, False, 1),
        ([-1, -2, 1], False, False, -1),
        ([1, 1 + 1e-8, 2], True, False, 1),
    ],
)
def test_source_candidate_rule(values, raw, robust, direction):
    result = classify_fseof_trajectory(values, [1, 2, 4], 1e-6)
    assert result["raw_candidate"] is raw
    assert result["robust_candidate"] is robust
    assert result["direction"] == direction
    assert result["slope"] == pytest.approx(abs(values[-1] - values[0]) / 3)


def test_reverse_flux_candidate_and_point_optimal_growth():
    result = run_published_fseof(
        metabolic_toy(reverse=True), "BIO", "PRODUCT", source_identity={"commit": "test-fixture"}
    )
    assert result["maximum_product"] == pytest.approx(10)
    assert result["enforced_product_lower_bounds"] == pytest.approx([0.9 * i for i in range(1, 11)])
    rows = {row["reaction_id"]: row for row in result["reaction_rows"]}
    assert rows["REV"]["raw_candidate"]
    assert rows["REV"]["robust_candidate"]
    assert rows["REV"]["direction"] == -1
    assert rows["REV"]["slope"] == pytest.approx(1)
    assert rows["REV"]["gene_ids"] == ["g1", "g2"]
    assert not rows["BIO"]["candidate"]
    for scan in result["scans"]:
        assert scan["point_optimal_growth"] == pytest.approx(10 - scan["product_lower_bound"])
        assert scan["growth_lower_bound"] == pytest.approx(scan["point_optimal_growth"] * 0.999999)
        assert scan["representative_growth"] >= scan["growth_lower_bound"] - 1e-6
        assert scan["secondary_abs_flux_objective"] > result["maximum_product"]
    assert result["metadata"]["secondary_growth_fraction"] == 0.999999
    assert result["metadata"]["source_identity"]["commit"] == "test-fixture"
    json.dumps(result, allow_nan=False)


def test_native_highspy_fseof_witnesses_are_strict_json():
    pytest.importorskip("highspy")
    from geno2flux import highspy_backend

    interface = highspy_backend._load_interface_module()
    interface.register_highs_interface()
    model = metabolic_toy()
    model.solver = "highs"
    model.solver.configuration.tolerances.feasibility = 1e-9
    model.solver.configuration.continuous_lp_backend = "linprog"
    with highspy_backend.scoped_highspy_backend():
        highspy_backend.install_highspy_backend()
        result = run_published_fseof(
            model, "BIO", "PRODUCT", config=FSEOFConfig(iterations=2),
            source_identity={"fixture": "native-highspy-json"},
        )
    decoded = json.loads(json.dumps(result, allow_nan=False))
    assert decoded["maximum_product"] == pytest.approx(10)
    witnesses = [decoded["maximum_product_witness"]]
    for scan in decoded["scans"]:
        witnesses.extend([scan["growth_optimum_witness"], scan["representative_witness"]])
    for witness in witnesses:
        receipt = witness["audit"]["backend_receipt"]
        assert receipt["backend"] == "highspy.Highs native"
        assert receipt["status"] == "optimal"
        assert receipt["primals_from_backend"] is True


def test_product_is_lower_bound_not_equality():
    model = cobra.Model("growth_coproduces_product")
    model.solver = "glpk"
    a, b, c = [cobra.Metabolite(name, compartment="c") for name in ("a", "b", "c")]
    uptake = cobra.Reaction("UP")
    uptake.add_metabolites({a: 1})
    uptake.bounds = (0, 10)
    coproduce = cobra.Reaction("COPRODUCE")
    coproduce.add_metabolites({a: -1, b: 1, c: 1})
    product = cobra.Reaction("PRODUCT")
    product.add_metabolites({b: -1})
    biomass = cobra.Reaction("BIO")
    biomass.add_metabolites({c: -1})
    model.add_reactions([uptake, coproduce, product, biomass])
    model.objective = "BIO"
    result = run_published_fseof(model, "BIO", "PRODUCT")
    # Equality would restrict growth to the scan's enforced level.  The
    # published lower bound retains growth near 10 and excess product.
    assert result["scans"][0]["product_lower_bound"] == pytest.approx(0.9)
    assert all(scan["point_optimal_growth"] == pytest.approx(10) for scan in result["scans"])
    assert all(scan["product_flux"] > scan["product_lower_bound"] for scan in result["scans"])


def test_real_extra_constraints_and_auxiliary_variables_are_retained():
    model = metabolic_toy()
    auxiliary = model.problem.Variable("capacity_proxy", lb=0, ub=6)
    model.add_cons_vars([auxiliary])
    model.add_cons_vars([model.problem.Constraint(
        model.reactions.UP.flux_expression - auxiliary, ub=0, name="extra_capacity_link"
    )])
    result = run_published_fseof(model, "BIO", "PRODUCT")
    assert result["maximum_product"] == pytest.approx(6)
    for scan in result["scans"]:
        for label in ("growth_optimum_witness", "representative_witness"):
            witness = scan[label]
            assert witness["constraint_values"]["extra_capacity_link"] <= 1e-6
            assert witness["primal_values"]["capacity_proxy"] <= 6 + 1e-6
            assert witness["fluxes"]["UP"] <= 6 + 1e-6
            assert witness["audit"]["solver_constraint_count"] >= len(model.constraints)


def test_forged_auxiliary_assignment_is_rejected_by_whole_lp_audit():
    model = metabolic_toy()
    auxiliary = model.problem.Variable("capacity_proxy", lb=0, ub=10)
    model.add_cons_vars([auxiliary])
    model.add_cons_vars([model.problem.Constraint(
        model.reactions.UP.flux_expression - auxiliary, ub=0, name="extra_capacity_link"
    )])
    solution = model.optimize()
    primals = dict(model.solver.primal_values)
    primals["capacity_proxy"] = 0
    with pytest.raises(PublishedFSEOFError, match="Invalid flux witness"):
        validate_flux_witness(model, solution, 1e-6, primal_values=primals)


def test_nonzero_metabolite_row_bounds_are_not_replaced_with_zero():
    model = metabolic_toy()
    model.constraints["a_c"].ub = 2
    model.constraints["a_c"].lb = 2
    result = run_published_fseof(model, "BIO", "PRODUCT")
    assert result["maximum_product"] == pytest.approx(8)
    assert result["maximum_product_witness"]["constraint_values"]["a_c"] == pytest.approx(2)


def test_model_is_restored_after_success():
    model = metabolic_toy()
    model.reactions.PRODUCT.lower_bound = 2
    constraints = [constraint.name for constraint in model.constraints]
    objective = str(model.objective.expression)
    result = run_published_fseof(model, "BIO", "PRODUCT")
    assert result["product_original_bounds"] == [2, 100]
    # RAVEN assigns the bound, rather than max(old_bound, scan_bound).
    assert result["scans"][0]["product_flux"] == pytest.approx(0.9)
    assert model.reactions.PRODUCT.bounds == (2, 100)
    assert [constraint.name for constraint in model.constraints] == constraints
    assert str(model.objective.expression) == objective
    assert model.objective_direction == "max"


def test_infeasible_model_is_not_a_zero_candidate_result():
    model = metabolic_toy()
    model.reactions.UP.upper_bound = 0
    model.reactions.PRODUCT.lower_bound = 1
    objective = str(model.objective.expression)
    with pytest.raises(PublishedFSEOFError, match="Non-optimal|failed"):
        run_published_fseof(model, "BIO", "PRODUCT")
    assert str(model.objective.expression) == objective
    assert model.reactions.PRODUCT.lower_bound == 1


def test_unbounded_product_is_rejected():
    model = metabolic_toy()
    for reaction in model.reactions:
        reaction.upper_bound = float("inf")
    with pytest.raises(PublishedFSEOFError, match="failed|Non-optimal"):
        run_published_fseof(model, "BIO", "PRODUCT")


def test_missing_or_nonfinite_primals_are_rejected():
    model = metabolic_toy()
    solution = model.optimize()
    primals = dict(model.solver.primal_values)
    first = next(iter(primals))
    del primals[first]
    with pytest.raises(PublishedFSEOFError, match="Incomplete"):
        validate_flux_witness(model, solution, 1e-6, primal_values=primals)
    primals = dict(model.solver.primal_values)
    primals[first] = float("nan")
    with pytest.raises(PublishedFSEOFError, match="Nonfinite"):
        validate_flux_witness(model, solution, 1e-6, primal_values=primals)


def test_forged_reaction_flux_assignment_is_rejected():
    model = metabolic_toy()
    solution = model.optimize()
    solution.fluxes["BIO"] = 0
    with pytest.raises(PublishedFSEOFError, match="Invalid flux witness"):
        validate_flux_witness(model, solution, 1e-6)


def test_integer_model_is_not_silently_called_ordinary_fseof():
    model = metabolic_toy()
    switch = model.problem.Variable("switch", type="binary")
    model.add_cons_vars([switch])
    with pytest.raises(PublishedFSEOFError, match="continuous LP"):
        run_published_fseof(model, "BIO", "PRODUCT")


def test_hiGHS_glpk_parity_on_unique_path_toy():
    bundle_src = Path(__file__).resolve().parents[1] / "third_party" / "alphagem" / "src"
    if str(bundle_src) not in sys.path:
        sys.path.insert(0, str(bundle_src))
    import optlang_highs_interface as highs

    previous = cobra.util.solver.solvers.get("highs")
    highs.register_highs_interface()
    try:
        glpk_model = metabolic_toy(reverse=True)
        highs_model = glpk_model.copy()
        highs_model.solver = "highs"
        glpk_result = run_published_fseof(glpk_model, "BIO", "PRODUCT")
        highs_result = run_published_fseof(highs_model, "BIO", "PRODUCT")
        assert highs_result["maximum_product"] == pytest.approx(glpk_result["maximum_product"])
        for highs_row, glpk_row in zip(highs_result["flux_matrix"], glpk_result["flux_matrix"]):
            assert highs_row == pytest.approx(glpk_row, abs=1e-6)
        for highs_row, glpk_row in zip(highs_result["reaction_rows"], glpk_result["reaction_rows"]):
            assert highs_row["robust_candidate"] == glpk_row["robust_candidate"]
    finally:
        if previous is None:
            cobra.util.solver.solvers.pop("highs", None)
        else:
            cobra.util.solver.solvers["highs"] = previous


@pytest.mark.parametrize("kwargs", [
    {"iterations": 1}, {"iterations": True}, {"maximum_fraction": 1},
    {"maximum_fraction": float("nan")}, {"feasibility_tolerance": 0},
])
def test_invalid_scan_config_is_rejected(kwargs):
    with pytest.raises(ValueError):
        FSEOFConfig(**kwargs)


def test_nonfinite_source_identity_is_rejected():
    with pytest.raises(ValueError, match="source_identity"):
        run_published_fseof(metabolic_toy(), "BIO", "PRODUCT", source_identity={"value": float("nan")})
