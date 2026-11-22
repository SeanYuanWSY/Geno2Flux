"""Tests for the bundled optlang-compatible HiGHS interface (scipy milp)."""

import sys
from pathlib import Path

import cobra
import pytest


BUNDLE_ROOT = Path(__file__).resolve().parents[1] / "third_party" / "alphagem"
for entry in (str(BUNDLE_ROOT), str(BUNDLE_ROOT / "src")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

import optlang_highs_interface as highs  # noqa: E402


@pytest.fixture(autouse=True)
def restore_solver():
    original = cobra.Configuration().solver
    yield
    cobra.Configuration().solver = original
    cobra.util.solver.solvers.pop("highs", None)


@pytest.fixture
def registered():
    assert highs.register_highs_interface() == "highs"
    return highs


def _toy_metabolic_model():
    model = cobra.Model("toy")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    uptake = cobra.Reaction("EX_a")
    uptake.add_metabolites({a: -1})
    uptake.bounds = (-1000, 1000)
    convert = cobra.Reaction("R1")
    convert.add_metabolites({a: -1, b: 1})
    convert.bounds = (0, 1000)
    biomass = cobra.Reaction("BIO")
    biomass.add_metabolites({b: -1})
    biomass.bounds = (0, 1000)
    model.add_reactions([uptake, convert, biomass])
    model.objective = "BIO"
    return model


def test_registration_makes_highs_a_cobra_solver(registered):
    assert cobra.util.solver.solvers["highs"] is highs
    cobra.Configuration().solver = "highs"
    model = _toy_metabolic_model()
    assert model.solver.interface is highs


def test_lp_parity_with_glpk_on_metabolic_model(registered):
    objectives = {}
    for solver in ("glpk", "highs"):
        cobra.Configuration().solver = solver
        model = _toy_metabolic_model()
        solution = model.optimize()
        assert solution.status == "optimal"
        objectives[solver] = solution.objective_value
    assert abs(objectives["glpk"] - objectives["highs"]) < 1e-7


def test_binary_milp_forcing_and_primal(registered):
    cobra.Configuration().solver = "highs"
    model = cobra.Model("milp")
    interface = model.problem
    keep = interface.Variable("keep", type="binary")
    flux = interface.Variable("flux", lb=0.0, ub=10.0)
    model.add_cons_vars([keep, flux])
    model.add_cons_vars([
        interface.Constraint(flux - 10.0 * keep, ub=0, name="link"),
        interface.Constraint(flux, lb=5.0, name="floor"),
    ])
    model.objective = interface.Objective(keep, direction="min")
    solution = model.optimize()
    assert solution.status == "optimal"
    assert abs(solution.objective_value - 1.0) < 1e-7
    assert keep.primal > 0.5
    assert abs(flux.primal - 5.0) < 1e-6


def test_infeasible_and_unbounded_status_mapping(registered):
    cobra.Configuration().solver = "highs"
    infeasible = cobra.Model("infeasible")
    interface = infeasible.problem
    x = interface.Variable("x", lb=0.0, ub=1.0)
    infeasible.add_cons_vars([
        x,
        interface.Constraint(x, lb=5.0, name="impossible"),
    ])
    infeasible.objective = interface.Objective(x, direction="max")
    assert infeasible.optimize().status == "infeasible"

    unbounded = cobra.Model("unbounded")
    y = interface.Variable("y", lb=0.0)
    unbounded.add_cons_vars([y])
    unbounded.objective = interface.Objective(y, direction="max")
    # optlang-level mapping is honest; cobra's optimize() raises on unbounded
    # regardless of backend (shared cobra behavior).
    assert unbounded.solver.optimize() == "unbounded"
    with pytest.raises(cobra.exceptions.OptimizationError, match="unbounded"):
        unbounded.optimize()


def test_limit_statuses_map_honestly(registered, monkeypatch):
    cobra.Configuration().solver = "highs"
    model = cobra.Model("limit")
    interface = model.problem
    x = interface.Variable("x", lb=0.0, ub=1.0)
    model.add_cons_vars([x])
    model.objective = interface.Objective(x, direction="max")

    class FakeResult:
        status = 1
        success = False
        message = "Iteration or time limit reached"
        x = None
        fun = None
        mip_node_count = None
        mip_dual_bound = None
        mip_gap = None

    monkeypatch.setattr(
        highs.scipy.optimize, "milp", lambda *a, **k: FakeResult()
    )
    # optlang-level status mapping; cobra's full-Solution path needs primals,
    # which a limit result does not provide (same for GLPK).
    assert model.solver.optimize() == "iteration_limit"
    model.solver.configuration.timeout = 5.0
    assert model.solver.optimize() == "time_limit"
    receipt = model.solver.problem.last_result_receipt
    assert receipt["raw_status"] == 1
    assert receipt["success"] is False
    assert receipt["forwarded_options"]["time_limit"] == 5.0


def test_context_manager_rolls_back_added_constraints(registered):
    cobra.Configuration().solver = "highs"
    model = _toy_metabolic_model()
    before_variables = [v.name for v in model.solver.variables]
    before_constraints = [c.name for c in model.solver.constraints]
    before_objective = model.optimize().objective_value
    with model:
        interface = model.problem
        keep = interface.Variable("temp_keep", type="binary")
        model.add_cons_vars([
            keep,
            interface.Constraint(keep, lb=1, name="temp_force"),
        ])
        assert "temp_force" in [c.name for c in model.solver.constraints]
        assert model.optimize().status == "optimal"
    assert [v.name for v in model.solver.variables] == before_variables
    assert [c.name for c in model.solver.constraints] == before_constraints
    assert abs(model.optimize().objective_value - before_objective) < 1e-7


def test_solver_switch_clones_glpk_model_into_highs(registered):
    cobra.Configuration().solver = "glpk"
    model = _toy_metabolic_model()
    glpk_objective = model.optimize().objective_value
    model.solver = "highs"
    assert model.solver.interface is highs
    solution = model.optimize()
    assert solution.status == "optimal"
    assert abs(solution.objective_value - glpk_objective) < 1e-7


def test_model_copy_roundtrips_through_highs(registered):
    cobra.Configuration().solver = "highs"
    model = _toy_metabolic_model()
    copy = model.copy()
    assert copy.solver.interface is highs
    solution = copy.optimize()
    assert solution.status == "optimal"
    assert abs(solution.objective_value - model.optimize().objective_value) < 1e-7


def test_continuous_linprog_options_copy_and_rollback(registered, monkeypatch):
    model = _toy_metabolic_model(); model.solver = "highs"
    model.solver.configuration.continuous_lp_backend = "linprog"
    model.solver.configuration.tolerances.feasibility = 1e-9
    calls = []
    original = highs.scipy.optimize.linprog
    def record(*args, **kwargs):
        calls.append(kwargs.copy()); return original(*args, **kwargs)
    monkeypatch.setattr(highs.scipy.optimize, "linprog", record)
    assert model.optimize().objective_value == pytest.approx(1000)
    assert calls[-1]["options"]["primal_feasibility_tolerance"] == 1e-9
    assert calls[-1]["options"]["dual_feasibility_tolerance"] == 1e-9
    assert "mip_rel_gap" not in calls[-1]["options"]
    clone = model.copy()
    assert clone.solver.configuration.continuous_lp_backend == "linprog"
    assert clone.solver.configuration.tolerances.feasibility == 1e-9
    with clone:
        clone.reactions.R1.upper_bound = 5
        assert clone.optimize().objective_value == pytest.approx(5)
    assert clone.optimize().objective_value == pytest.approx(1000)
    receipt = clone.solver.problem.last_result_receipt
    assert receipt["backend"] == "scipy.optimize.linprog"


def test_continuous_linprog_translates_ranges_equalities_free_and_negative_variables(registered):
    model = cobra.Model("range"); model.solver = "highs"
    model.solver.configuration.continuous_lp_backend = "linprog"
    p = model.problem
    x, y, z = p.Variable("x", lb=None, ub=None), p.Variable("y", lb=-10, ub=-1), p.Variable("z")
    model.add_cons_vars([x, y, z])
    model.add_cons_vars([p.Constraint(x, lb=-3, ub=-2, name="range"),
                        p.Constraint(x+y, lb=-5, ub=-5, name="equal"),
                        p.Constraint(z, ub=0, name="upper"),
                        p.Constraint(z, lb=0, name="lower"),
                        p.Constraint(x+y+z, name="free")])
    model.objective = p.Objective(x, direction="max")
    assert model.optimize().objective_value == pytest.approx(-2)
    model.objective.direction = "min"
    assert model.optimize().objective_value == pytest.approx(-3)
    assert y.primal == pytest.approx(-2)


def test_linprog_selection_preserves_binary_MILP_route(registered, monkeypatch):
    model = cobra.Model("integer"); model.solver = "highs"
    model.solver.configuration.continuous_lp_backend = "linprog"
    p = model.problem; x = p.Variable("x", type="binary")
    model.add_cons_vars([x]); model.objective = p.Objective(x, direction="max")
    monkeypatch.setattr(highs.scipy.optimize, "linprog", lambda *a, **k: pytest.fail("integer sent to LP"))
    assert model.optimize().objective_value == pytest.approx(1)
    assert model.solver.problem.last_result_receipt["backend"] == "scipy.optimize.milp"


def test_continuous_linprog_honest_nonoptimal_status(registered):
    model = cobra.Model("infeasibleLP"); model.solver = "highs"
    model.solver.configuration.continuous_lp_backend = "linprog"
    p = model.problem; x = p.Variable("x", lb=0, ub=1)
    model.add_cons_vars([x, p.Constraint(x, lb=2)]); model.objective = p.Objective(x)
    assert model.solver.optimize() == "infeasible"
    assert not model.solver.problem.last_result_receipt["success"]


def test_published_fseof_with_real_continuous_LP_precision(registered):
    from geno2flux.pipeline import _attach_solver
    from geno2flux.published_fseof import run_published_fseof

    m = _toy_metabolic_model()
    product = cobra.Reaction("PRODUCT")
    product.add_metabolites({m.metabolites.b_c: -1})
    product.bounds = (0, 1000)
    m.add_reactions([product])
    _attach_solver(m, "highs")
    result = run_published_fseof(m, "BIO", "PRODUCT")
    assert len(result["scans"]) == 10
    assert all(
        s["representative_witness"]["audit"]["backend_receipt"]["backend"] == "scipy.optimize.linprog"
        for s in result["scans"]
    )
    m.solver = "glpk"
    m.solver.configuration.tolerances.feasibility = 1e-9
    if hasattr(m.solver.configuration.tolerances, "optimality"):
        m.solver.configuration.tolerances.optimality = 1e-9
    elif getattr(m.solver.configuration, "_smcp", None) is not None:
        m.solver.configuration._smcp.tol_dj = 1e-9
    other = run_published_fseof(m, "BIO", "PRODUCT")
    assert all(
        s["representative_witness"]["audit"]["actual_glpk_simplex_tolerances"]["tol_dj"] == 1e-9
        for s in other["scans"]
    )


def test_gapfiller_selects_same_reaction_on_both_backends(registered):
    from cobra.flux_analysis.gapfilling import GapFiller

    def build_model_missing_link():
        model = cobra.Model("toy_gap")
        a = cobra.Metabolite("a_c", compartment="c")
        b = cobra.Metabolite("b_c", compartment="c")
        c = cobra.Metabolite("c_c", compartment="c")
        uptake = cobra.Reaction("EX_a")
        uptake.add_metabolites({a: -1})
        uptake.bounds = (-1000, 1000)
        convert = cobra.Reaction("R1")
        convert.add_metabolites({a: -1, b: 1})
        convert.bounds = (0, 1000)
        biomass = cobra.Reaction("BIO")
        biomass.add_metabolites({c: -1})
        biomass.bounds = (0, 1000)
        model.add_reactions([uptake, convert, biomass])
        model.objective = "BIO"
        return model

    def build_universal():
        universal = cobra.Model("universal")
        b2 = cobra.Metabolite("b_c", compartment="c")
        c2 = cobra.Metabolite("c_c", compartment="c")
        r2 = cobra.Reaction("R2")
        r2.add_metabolites({b2: -1, c2: 1})
        r2.bounds = (0, 1000)
        universal.add_reactions([r2])
        return universal

    results = {}
    for solver in ("glpk", "highs"):
        cobra.Configuration().solver = solver
        model = build_model_missing_link()
        gapfiller = GapFiller(
            model, build_universal(), lower_bound=0.05,
            demand_reactions=False, integer_threshold=1e-9,
        )
        filled = gapfiller.fill()
        results[solver] = {
            "selected": sorted(reaction.id for reaction in filled[0]),
            "objective": float(gapfiller.model.solver.objective.value),
            "status": model.solver.status,
        }
        if solver == "highs":
            receipt = gapfiller.model.solver.problem.last_result_receipt
            results[solver]["receipt"] = receipt
    assert results["glpk"]["selected"] == results["highs"]["selected"] == ["R2"]
    assert results["glpk"]["status"] == results["highs"]["status"] == "optimal"
    assert abs(results["glpk"]["objective"] - results["highs"]["objective"]) < 1e-7
    assert results["highs"]["receipt"]["raw_status"] == 0





def test_prune_unsupported_reactions_runs_identically_on_highs(registered):
    gapfilling = pytest.importorskip("src.gapfilling")

    def build():
        model = cobra.Model("prune")
        a = cobra.Metabolite("a_c", compartment="c")
        b = cobra.Metabolite("b_c", compartment="c")
        c = cobra.Metabolite("c_c", compartment="c")
        uptake = cobra.Reaction("EX_a")
        uptake.add_metabolites({a: -1})
        uptake.bounds = (-1000, 1000)
        convert = cobra.Reaction("R1")
        convert.add_metabolites({a: -1, b: 1})
        convert.bounds = (0, 1000)
        gap_one = cobra.Reaction("GAP1")
        gap_one.add_metabolites({b: -1, c: 1})
        gap_one.bounds = (0, 1000)
        gap_two = cobra.Reaction("GAP2")
        gap_two.add_metabolites({b: -1, c: 1})
        gap_two.bounds = (0, 1000)
        biomass = cobra.Reaction("BIO")
        biomass.add_metabolites({c: -1})
        biomass.bounds = (0, 1000)
        model.add_reactions([uptake, convert, gap_one, gap_two, biomass])
        model.objective = "BIO"
        return model

    outcomes = {}
    for solver in ("glpk", "highs"):
        cobra.Configuration().solver = solver
        model = build()
        outcomes[solver] = gapfilling._prune_unsupported_reactions(
            model, ["GAP1", "GAP2"]
        )
    for solver in ("glpk", "highs"):
        outcome = outcomes[solver]
        # Exactly one of the two interchangeable gap reactions is retained.
        assert outcome["minimum_retained_reaction_count"] == 1
        assert len(outcome["retained_algorithmic_reactions"]) == 1
        assert outcome["minimum_set_identity_is_unique"] is False
        assert outcome["final_growth_status"] == "optimal"
    assert (
        outcomes["glpk"]["retained_algorithmic_reactions"]
        == outcomes["highs"]["retained_algorithmic_reactions"]
    )
    assert outcomes["highs"]["growth_retention_rule"]["solver"]["interface"] == (
        "optlang_highs_interface"
    )


def test_rename_middle_variable_keeps_column_mapping(registered):
    """Renaming a mounted variable must not scramble problem columns."""
    import optlang.glpk_interface

    outcomes = {}
    for solver, interface in (
        ("glpk", optlang.glpk_interface),
        ("highs", highs),
    ):
        model = interface.Model(name="rename")
        x = interface.Variable("x", lb=0.0, ub=10.0)
        y = interface.Variable("y", lb=0.0, ub=10.0)
        z = interface.Variable("z", lb=0.0, ub=10.0)
        model.add([x, y, z])
        model.add([
            interface.Constraint(2.0 * x + 3.0 * y, ub=12.0, name="c1"),
            interface.Constraint(z, lb=4.0, name="c2"),
        ])
        # Non-uniform objective: a scrambled column mapping changes the value.
        model.objective = interface.Objective(x + 2.0 * y + z, direction="max")
        # Rename the middle variable after the problem is fully assembled.
        y.name = "w"
        status = model.optimize()
        outcomes[solver] = (
            status,
            model.objective.value,
            {variable.name: variable.primal for variable in model.variables},
        )
    assert outcomes["glpk"][0] == outcomes["highs"][0] == "optimal"
    assert abs(outcomes["glpk"][1] - outcomes["highs"][1]) < 1e-7
    # Both backends renamed y -> w; the highs column mapping must not
    # scramble: every renamed/untouched variable matches glpk's value.
    for key in ("x", "w", "z"):
        assert abs(outcomes["highs"][2][key] - outcomes["glpk"][2][key]) < 1e-7


def test_solver_switch_then_infeasible_does_not_crash(registered):
    """glpk->highs clone injects presolve='auto'; a non-optimal solve must not
    trip the validating presolve setter via optlang's auto-retry."""
    cobra.Configuration().solver = "glpk"
    model = _toy_metabolic_model()
    model.solver = "highs"
    assert model.solver.interface is highs
    interface = model.problem
    blocker = interface.Constraint(
        model.reactions.BIO.forward_variable, lb=2000.0, name="blocker"
    )
    model.add_cons_vars([blocker])
    solution = model.optimize()
    assert solution.status == "infeasible"
    assert model.solver.configuration.presolve is True
    assert model.solver.configuration.presolve_auto_mapped is True


def test_pipeline_runner_configures_highs_solver(tmp_path):
    from src.pipeline_runner import PipelineRunner

    runner = PipelineRunner(tmp_path)
    runner._configure_solver("highs")
    assert cobra.Configuration().solver is highs
    assert cobra.util.solver.solvers["highs"] is highs
    runner._configure_solver("glpk")
    assert cobra.Configuration().solver != "highs"
