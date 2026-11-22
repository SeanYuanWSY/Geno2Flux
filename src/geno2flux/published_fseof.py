# SPDX-License-Identifier: GPL-3.0-only
"""Ordinary-GEM port of the published RAVEN FSEOF implementation.

FSEOF prioritizes flux-increase candidates; it does not simulate expression
changes or estimate a gene's intervention benefit.  The caller supplies a
prepared model, medium, product reaction, and the frozen source identities.
No medium, product chemistry, enzyme constraint, or growth floor is invented
here.  In particular, the 0.999999 fraction below is the secondary-solve
numerical relaxation in RAVEN ``solveLP.m``, not a biological growth fraction.

Primary implementation sources:
https://github.com/SysBioChalmers/RAVEN/blob/main/core/FSEOF.m
https://github.com/SysBioChalmers/RAVEN/blob/main/solver/solveLP.m

The exact strict comparison is retained as ``raw_candidate``.  A separate,
explicit engineering-tolerance flag is reported as ``robust_candidate``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from typing import Any, Mapping, Sequence

import cobra


RAVEN_GROWTH_OPTIMALITY_FRACTION = 0.999999
FSEOF_SOURCE_URL = (
    "https://raw.githubusercontent.com/SysBioChalmers/RAVEN/main/core/FSEOF.m"
)
SOLVELP_SOURCE_URL = (
    "https://raw.githubusercontent.com/SysBioChalmers/RAVEN/main/solver/solveLP.m"
)


class PublishedFSEOFError(RuntimeError):
    """A failed solve or invalid witness cannot be converted to a zero score."""


@dataclass(frozen=True)
class FSEOFConfig:
    """Published scan settings and separately named numerical tolerances."""

    iterations: int = 10
    maximum_fraction: float = 0.9
    feasibility_tolerance: float = 1e-6
    monotonicity_tolerance: float = 1e-6

    def __post_init__(self) -> None:
        if isinstance(self.iterations, bool) or not isinstance(self.iterations, int):
            raise ValueError("iterations must be an integer")
        if self.iterations < 2:
            raise ValueError("At least two iterations are needed for the published slope")
        if not math.isfinite(self.maximum_fraction) or not 0 < self.maximum_fraction < 1:
            raise ValueError("maximum_fraction must be strictly between zero and one")
        for name in ("feasibility_tolerance", "monotonicity_tolerance"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")


def classify_fseof_trajectory(
    fluxes: Sequence[float], enforced_levels: Sequence[float], tolerance: float
) -> dict[str, Any]:
    """Apply FSEOF.m's first-round direction and strict successive increases.

    ``slope`` is RAVEN's endpoint slope for *all* reactions, including those
    that fail its candidate rule.  It is meaningful for prioritization only
    when accompanied by the candidate flag.
    """
    values = [float(value) for value in fluxes]
    levels = [float(value) for value in enforced_levels]
    if len(values) < 2 or len(values) != len(levels):
        raise ValueError("A trajectory needs at least two matching scan levels")
    if not all(math.isfinite(value) for value in values + levels):
        raise ValueError("Nonfinite trajectory")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if not all(right > left for left, right in zip(levels, levels[1:])):
        raise ValueError("Enforced levels must strictly increase")
    direction = 1 if values[0] > 0 else -1 if values[0] < 0 else 0
    increments = [direction * (right - left) for left, right in zip(values, values[1:])]
    raw = bool(direction and all(change > 0 for change in increments))
    robust = bool(
        raw and abs(values[0]) > tolerance
        and all(change > tolerance for change in increments)
    )
    return {
        "direction": direction,
        "raw_candidate": raw,
        "robust_candidate": robust,
        "numerically_fragile_candidate": bool(raw and not robust),
        "slope": abs(values[-1] - values[0]) / (levels[-1] - levels[0]),
        "initial_flux": values[0],
        "final_flux": values[-1],
        "minimum_directional_increment": min(increments) if direction else None,
    }


def _finite_float(value: Any, description: str) -> float:
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise PublishedFSEOFError(f"Missing numeric {description}") from exc
    if not math.isfinite(number):
        raise PublishedFSEOFError(f"Nonfinite {description}")
    return number


def _linear_value(expression: Any, primals: Mapping[str, float]) -> float:
    """Recompute a linear expression, rather than trusting solver row primals."""
    value = 0.0
    for term, coefficient in expression.as_coefficients_dict().items():
        coefficient = _finite_float(coefficient, "linear coefficient")
        if term == 1:
            value += coefficient
        elif term.is_Symbol:
            if term.name not in primals:
                raise PublishedFSEOFError(f"Missing variable {term.name}")
            value += coefficient * primals[term.name]
        else:
            raise PublishedFSEOFError("FSEOF requires a linear continuous model")
    return _finite_float(value, "recomputed expression")


def _bound_violation(value: float, lower: Any, upper: Any) -> float:
    result = 0.0
    if lower is not None:
        lower = float(lower)
        if math.isnan(lower):
            raise PublishedFSEOFError("NaN lower bound")
        result = max(result, lower - value)
    if upper is not None:
        upper = float(upper)
        if math.isnan(upper):
            raise PublishedFSEOFError("NaN upper bound")
        result = max(result, value - upper)
    return result


def validate_flux_witness(
    model: cobra.Model,
    solution: Any,
    tolerance: float,
    *,
    primal_values: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Validate every active solver row and variable while that LP is intact.

    Net reaction fluxes alone are insufficient: extra optlang constraints
    can contain auxiliary variables or individually split fluxes.  The full
    primal assignment is therefore recorded and checked, along with actual
    metabolite-row bounds, reaction bounds and the current objective.
    ``status == optimal`` certifies solver-reported optimality; no independent
    dual certificate is claimed by this function.
    """
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if solution.status != "optimal" or model.solver.status != "optimal":
        raise PublishedFSEOFError(f"Non-optimal solve: {solution.status}")
    reaction_ids = [reaction.id for reaction in model.reactions]
    if list(solution.fluxes.index) != reaction_ids:
        raise PublishedFSEOFError("Solution reaction identities/order do not match the model")
    fluxes = {
        rid: _finite_float(solution.fluxes[rid], f"flux {rid}") for rid in reaction_ids
    }
    raw_primals = model.solver.primal_values if primal_values is None else primal_values
    expected_variables = {variable.name for variable in model.solver.variables}
    if set(raw_primals) != expected_variables:
        raise PublishedFSEOFError("Incomplete or unexpected solver primal variables")
    primals = {
        name: _finite_float(value, f"primal {name}") for name, value in raw_primals.items()
    }
    variable_error = max(
        (_bound_violation(primals[v.name], v.lb, v.ub) for v in model.solver.variables),
        default=0.0,
    )
    constraint_values: dict[str, float] = {}
    constraint_error = 0.0
    worst_constraint = None
    for constraint in model.solver.constraints:
        value = _linear_value(constraint.expression, primals)
        constraint_values[constraint.name] = value
        error = _bound_violation(value, constraint.lb, constraint.ub)
        if error > constraint_error:
            constraint_error, worst_constraint = error, constraint.name
    reaction_error = 0.0
    primal_flux_error = 0.0
    balances = {metabolite.id: 0.0 for metabolite in model.metabolites}
    for reaction in model.reactions:
        flux = fluxes[reaction.id]
        reaction_error = max(
            reaction_error, _bound_violation(flux, reaction.lower_bound, reaction.upper_bound)
        )
        net_primal = primals[reaction.forward_variable.name] - primals[reaction.reverse_variable.name]
        primal_flux_error = max(primal_flux_error, abs(flux - net_primal))
        for metabolite, coefficient in reaction.metabolites.items():
            balances[metabolite.id] += float(coefficient) * flux
    mass_error = 0.0
    for metabolite in model.metabolites:
        balance = _finite_float(balances[metabolite.id], f"balance {metabolite.id}")
        constraint = model.constraints[metabolite.id]
        mass_error = max(mass_error, _bound_violation(balance, constraint.lb, constraint.ub))
    objective = _finite_float(solution.objective_value, "solution objective")
    computed_objective = _linear_value(model.objective.expression, primals)
    objective_error = abs(objective - computed_objective)
    audit = {
        "maximum_variable_bound_violation": variable_error,
        "maximum_solver_constraint_violation": constraint_error,
        "worst_solver_constraint": worst_constraint,
        "maximum_reaction_bound_violation": reaction_error,
        "maximum_metabolite_bound_violation": mass_error,
        "maximum_flux_primal_disagreement": primal_flux_error,
        "objective_recalculation_error": objective_error,
        "absolute_feasibility_tolerance": float(tolerance),
        "solver_constraint_count": len(model.solver.constraints),
        "solver_variable_count": len(model.solver.variables),
        "solver_reported_status": "optimal",
        "optimality_certificate": "solver status; independent dual certificate not computed",
        "solver_interface": model.solver.interface.__name__,
        "configured_feasibility_tolerance": model.solver.configuration.tolerances.feasibility,
        "backend_receipt": getattr(model.solver.problem, "last_result_receipt", None),
    }
    simplex = getattr(model.solver.configuration, "_smcp", None)
    if simplex is not None:
        audit["actual_glpk_simplex_tolerances"] = {
            "tol_bnd": float(simplex.tol_bnd), "tol_dj": float(simplex.tol_dj),
            "tol_piv": float(simplex.tol_piv),
        }
    if max(variable_error, constraint_error, reaction_error, mass_error,
           primal_flux_error, objective_error) > tolerance:
        raise PublishedFSEOFError("Invalid flux witness: " + repr(audit))
    return {
        "fluxes": fluxes,
        "primal_values": primals,
        "constraint_values": constraint_values,
        "objective_value": objective,
        "audit": audit,
    }


def _optimal_witness(model: cobra.Model, tolerance: float, stage: str) -> dict[str, Any]:
    try:
        solution = model.optimize()
        return validate_flux_witness(model, solution, tolerance)
    except (cobra.exceptions.OptimizationError, ValueError) as exc:
        raise PublishedFSEOFError(f"{stage} failed: {exc}") from exc


def run_published_fseof(
    model: cobra.Model,
    biomass_reaction: str,
    product_reaction: str,
    *,
    config: FSEOFConfig | Mapping[str, Any] | None = None,
    source_identity: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return JSON-serializable scan witnesses and reaction-level candidates.

    The supplied model is restored on return or failure.  Existing extra
    constraints, variable bounds, and nonzero metabolite-row bounds are
    retained.  The published assignment of the target reaction's *lower*
    bound is exact: the caller's upper bound is retained, but any preexisting
    lower bound is replaced during the scan and restored afterwards.

    Source identities are provenance supplied by the caller, not verified
    or fabricated here.  ``candidate`` aliases RAVEN's strict raw rule, and
    ``robust_candidate`` is a separately reported numerical sensitivity flag.
    """
    if config is None:
        config = FSEOFConfig()
    elif isinstance(config, Mapping):
        config = FSEOFConfig(**config)
    elif not isinstance(config, FSEOFConfig):
        raise TypeError("config must be FSEOFConfig or a matching mapping")
    if biomass_reaction == product_reaction:
        raise ValueError("Biomass and product reactions must differ")
    model.reactions.get_by_id(biomass_reaction)
    model.reactions.get_by_id(product_reaction)
    if any(variable.type != "continuous" for variable in model.solver.variables):
        raise PublishedFSEOFError("Ordinary RAVEN FSEOF requires continuous LP variables")
    source_identity = dict(source_identity or {})
    try:
        json.dumps(source_identity, allow_nan=False)
    except (ValueError, TypeError) as exc:
        raise ValueError("source_identity must contain finite JSON-serializable values") from exc
    reaction_ids = [reaction.id for reaction in model.reactions]
    tolerance = config.feasibility_tolerance
    original_product_bounds = list(model.reactions.get_by_id(product_reaction).bounds)
    original_constraint_names = [constraint.name for constraint in model.solver.constraints]
    scans: list[dict[str, Any]] = []
    with model:
        product = model.reactions.get_by_id(product_reaction)
        biomass = model.reactions.get_by_id(biomass_reaction)
        model.objective = product
        model.objective_direction = "max"
        maximum_witness = _optimal_witness(model, tolerance, "Maximum product LP")
        maximum_product = maximum_witness["fluxes"][product_reaction]
        if maximum_product <= tolerance:
            raise PublishedFSEOFError("No positive product maximum above feasibility tolerance")
        scaled_maximum = config.maximum_fraction * maximum_product
        levels = [step * scaled_maximum / config.iterations
                  for step in range(1, config.iterations + 1)]
        model.objective = biomass
        model.objective_direction = "max"
        for index, lower_bound in enumerate(levels, 1):
            product.lower_bound = lower_bound
            growth_witness = _optimal_witness(model, tolerance, f"Scan {index} growth LP")
            maximum_growth = growth_witness["fluxes"][biomass_reaction]
            if maximum_growth < -tolerance:
                raise PublishedFSEOFError(f"Scan {index} has negative optimal biomass flux")
            # The frozen solveLP.m lines 176-179 allow 1e-6 relative
            # objective loss (line numbers refer to the caller's snapshot).
            # The negative branch is retained for numerical negatives near 0.
            growth_floor = maximum_growth * (
                RAVEN_GROWTH_OPTIMALITY_FRACTION if maximum_growth > 0 else 1.000001
            )
            with model:
                floor = model.problem.Constraint(
                    biomass.flux_expression, lb=growth_floor,
                    name=f"published_fseof_growth_floor_{index}",
                )
                model.add_cons_vars([floor])
                model.objective = model.problem.Objective(0, direction="min")
                model.objective.set_linear_coefficients({
                    variable: 1.0
                    for reaction in model.reactions
                    for variable in (reaction.forward_variable, reaction.reverse_variable)
                })
                representative = _optimal_witness(model, tolerance, f"Scan {index} abs-flux LP")
                growth = representative["fluxes"][biomass_reaction]
                if growth < growth_floor - tolerance or growth > maximum_growth + tolerance:
                    raise PublishedFSEOFError(f"Scan {index} violates point-optimal growth range")
                representative["audit"]["point_optimal_growth"] = maximum_growth
                representative["audit"]["required_growth_lower_bound"] = growth_floor
                representative["audit"]["representative_growth"] = growth
                scans.append({
                    "step": index,
                    "product_lower_bound": lower_bound,
                    "product_flux": representative["fluxes"][product_reaction],
                    "point_optimal_growth": maximum_growth,
                    "representative_growth": growth,
                    "growth_lower_bound": growth_floor,
                    "secondary_abs_flux_objective": representative["objective_value"],
                    "growth_optimum_witness": growth_witness,
                    "representative_witness": representative,
                })
        flux_matrix = [
            [scan["representative_witness"]["fluxes"][rid] for scan in scans]
            for rid in reaction_ids
        ]
        reaction_rows: list[dict[str, Any]] = []
        for reaction, fluxes in zip(model.reactions, flux_matrix):
            classification = classify_fseof_trajectory(
                fluxes, levels, config.monotonicity_tolerance
            )
            reaction_rows.append({
                "reaction_id": reaction.id,
                "reaction_name": reaction.name,
                "gpr": reaction.gene_reaction_rule,
                "gene_ids": sorted(gene.id for gene in reaction.genes),
                "boundary": bool(reaction.boundary),
                "reversible": bool(reaction.reversibility),
                "raven_reported_direction": (
                    classification["direction"] if reaction.reversibility else 0
                ),
                "candidate": classification["raw_candidate"],
                **classification,
            })
        result = {
            "schema_version": 1,
            "method": "RAVEN ordinary-GEM FSEOF Python port",
            "claim": "flux-increase candidate prioritization; not causal gene overexpression benefit",
            "config": asdict(config),
            "biomass_reaction": biomass_reaction,
            "product_reaction": product_reaction,
            "product_original_bounds": original_product_bounds,
            "maximum_product": maximum_product,
            "maximum_product_witness": maximum_witness,
            "maximum_product_witness_kind": "primary-product LP only, not secondary minimum-flux representative",
            "scaled_maximum_product": scaled_maximum,
            "enforced_product_lower_bounds": levels,
            "reaction_ids": reaction_ids,
            "flux_matrix": flux_matrix,
            "scans": scans,
            "reaction_rows": reaction_rows,
            "raw_candidate_count": sum(row["raw_candidate"] for row in reaction_rows),
            "robust_candidate_count": sum(row["robust_candidate"] for row in reaction_rows),
            "metadata": {
                "source_identity": source_identity,
                "source_urls": {"FSEOF.m": FSEOF_SOURCE_URL, "solveLP.m": SOLVELP_SOURCE_URL},
                "source_identity_verified_here": False,
                "scan_rule": "i * maximum_fraction * maximum_product / iterations, i=1..iterations",
                "product_constraint": "lower bound, not equality",
                "growth_rule": "maximize biomass at each scan; no global WT growth floor",
                "secondary_objective": "minimize sum of forward and reverse reaction flux variables",
                "secondary_growth_fraction": RAVEN_GROWTH_OPTIMALITY_FRACTION,
                "secondary_growth_fraction_origin": "frozen RAVEN solveLP.m lines 176-179",
                "maximum_product_port_adaptation": "primary product optimum supplies the same target-maximum value; its witness is not the source's secondary minimum-flux representative",
                "candidate_rule": "first flux nonzero; every following directional increment strictly positive",
                "slope_rule": "absolute endpoint flux difference / enforced-product lower-bound span",
                "numerical_adaptation": "robust_candidate additionally requires initial magnitude and every increment > monotonicity_tolerance; raw rule retained",
                "original_constraint_names": original_constraint_names,
                "solver_interface": model.solver.interface.__name__,
                "cobra_version": cobra.__version__,
            },
        }
    # Prevent JSON NaN/Infinity from turning an unsuccessful stage into an
    # apparently completed, portable scientific result.
    try:
        json.dumps(result, allow_nan=False)
    except (ValueError, TypeError) as exc:
        raise PublishedFSEOFError("Output is not finite JSON data") from exc
    return result
