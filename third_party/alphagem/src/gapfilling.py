from cobra.flux_analysis import gapfilling
import cobra
import csv
from Bio import SeqIO
from src import findrefnames
import pickle
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from sympy import Symbol, false, true
from sympy.logic import simplify_logic
from cobra.core.gene import GPR
from cobra.util.solver import linear_reaction_coefficients

from src.pipeline_contracts import (
    best_gapfill_hits,
    file_sha256,
    gapfill_media_plan,
    map_tblastn_hits_to_target_proteins,
    load_coordinate_translation_exceptions,
    protein_key,
    validate_gapfill_cds,
    validate_gapfill_genome,
)


def _finite_objective_or_none(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        objective = float(value)
    except (TypeError, ValueError):
        return None
    return objective if math.isfinite(objective) else None


def _require_optimal_finite_objective(
    solution,
    context,
    *,
    minimum_exclusive=None,
):
    status = getattr(solution, "status", None)
    objective = _finite_objective_or_none(
        getattr(solution, "objective_value", None)
    )
    if status != "optimal" or objective is None:
        raise RuntimeError(
            f"{context} requires an optimal solution with a finite objective: "
            f"status={status}, objective="
            f"{getattr(solution, 'objective_value', None)}"
        )
    if (
        minimum_exclusive is not None
        and objective <= float(minimum_exclusive)
    ):
        raise RuntimeError(
            f"{context} did not meet the required growth threshold: "
            f"status={status}, objective={objective}, "
            f"required>{float(minimum_exclusive)}"
        )
    return objective


def _solution_record(solution):
    return {
        "status": getattr(solution, "status", None),
        "objective_value": _finite_objective_or_none(
            getattr(solution, "objective_value", None)
        ),
    }


def _enforce_essential_reactions(
    model,
    reference_model,
    reaction_ids,
    *,
    require_all=False,
):
    records = []
    for reaction_id in sorted(set(reaction_ids)):
        try:
            reference_reaction = reference_model.reactions.get_by_id(
                reaction_id
            )
        except KeyError as error:
            if require_all:
                raise KeyError(
                    "Essential reaction is absent from reference model: "
                    f"{reaction_id}"
                ) from error
            records.append(
                {
                    "reaction_id": reaction_id,
                    "status": "missing_from_reference",
                    "final_present": reaction_id in model.reactions,
                }
            )
            continue
        if reaction_id in model.reactions:
            records.append(
                {
                    "reaction_id": reaction_id,
                    "status": "already_present",
                    "final_present": True,
                }
            )
            continue
        try:
            model.add_reactions([reference_reaction.copy()])
        except Exception as error:
            if require_all:
                raise RuntimeError(
                    "Failed to restore required essential reaction: "
                    f"{reaction_id}"
                ) from error
            records.append(
                {
                    "reaction_id": reaction_id,
                    "status": "add_failed",
                    "final_present": False,
                }
            )
            continue
        final_present = reaction_id in model.reactions
        if require_all and not final_present:
            raise RuntimeError(
                "Essential reaction addition returned without restoring "
                f"{reaction_id}"
            )
        records.append(
            {
                "reaction_id": reaction_id,
                "status": "added" if final_present else "add_failed",
                "final_present": final_present,
            }
        )
    return records


def _model_equivalence_issues(left, right):
    issues = []
    collections = (
        ("reaction", left.reactions, right.reactions),
        ("metabolite", left.metabolites, right.metabolites),
        ("gene", left.genes, right.genes),
    )
    for label, left_items, right_items in collections:
        left_ids = {item.id for item in left_items}
        right_ids = {item.id for item in right_items}
        if left_ids != right_ids:
            issues.append(
                f"{label} IDs differ: left_only="
                f"{sorted(left_ids - right_ids)[:20]}, right_only="
                f"{sorted(right_ids - left_ids)[:20]}"
            )
    shared_reactions = {
        reaction.id for reaction in left.reactions
    } & {reaction.id for reaction in right.reactions}
    for reaction_id in sorted(shared_reactions):
        left_reaction = left.reactions.get_by_id(reaction_id)
        right_reaction = right.reactions.get_by_id(reaction_id)
        left_stoichiometry = {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in left_reaction.metabolites.items()
        }
        right_stoichiometry = {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in right_reaction.metabolites.items()
        }
        if left_stoichiometry != right_stoichiometry:
            issues.append(f"{reaction_id} stoichiometry differs")
        if (
            not math.isclose(
                float(left_reaction.lower_bound),
                float(right_reaction.lower_bound),
                rel_tol=0.0,
                abs_tol=1e-9,
            )
            or not math.isclose(
                float(left_reaction.upper_bound),
                float(right_reaction.upper_bound),
                rel_tol=0.0,
                abs_tol=1e-9,
            )
        ):
            issues.append(f"{reaction_id} bounds differ")
        if (
            left_reaction.gene_reaction_rule
            != right_reaction.gene_reaction_rule
        ):
            issues.append(f"{reaction_id} GPR differs")
    left_objective = {
        reaction.id: float(coefficient)
        for reaction, coefficient in linear_reaction_coefficients(
            left
        ).items()
    }
    right_objective = {
        reaction.id: float(coefficient)
        for reaction, coefficient in linear_reaction_coefficients(
            right
        ).items()
    }
    if left_objective != right_objective:
        issues.append("objective coefficients differ")
    if left.objective.direction != right.objective.direction:
        issues.append("objective directions differ")
    return issues


def _read_blast_rows(path):
    rows = []
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return rows
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not row:
                continue
            if len(row) != 12:
                raise ValueError(
                    f"Malformed BLAST row at {path}:{line_number}; "
                    f"expected 12 fields, observed {len(row)}"
                )
            rows.append(
                {
                    "query": row[0],
                    "subject": row[1],
                    "evalue": float(row[10]),
                    "bitscore": float(row[11]),
                    "sstart": int(row[8]),
                    "send": int(row[9]),
                }
            )
    return rows


def _replace_gpr_gene(rule, old_gene, new_gene):
    pattern = re.compile(
        rf"(?<![A-Za-z0-9_.:-]){re.escape(old_gene)}"
        rf"(?![A-Za-z0-9_.:-])"
    )
    return pattern.sub(new_gene, rule)


def _rewrite_gpr_for_supported_hits(reaction, gene_hits):
    """Prune unsupported Boolean branches and replace supported genes."""

    if not reaction.gene_reaction_rule or not reaction.genes:
        return None
    expression = reaction.gpr.as_symbolic()
    substitutions = {}
    for symbol in expression.free_symbols:
        hit = gene_hits.get(str(symbol))
        substitutions[symbol] = (
            Symbol(protein_key(hit["subject"])) if hit is not None else false
        )
    reduced = simplify_logic(expression.subs(substitutions), force=True)
    if reduced in {false, true}:
        return None
    return GPR.from_symbolic(reduced).to_string()


def _apply_paper_reaction_support(
    model,
    reaction_ids,
    query_to_model_gene,
    hits,
):
    gene_hits = {
        query_to_model_gene[query]: hit
        for query, hit in hits.items()
        if query in query_to_model_gene
    }
    supported = []
    unsupported = []
    reaction_audit = {}
    for reaction_id in sorted(set(reaction_ids)):
        reaction = model.reactions.get_by_id(reaction_id)
        original_rule = reaction.gene_reaction_rule
        reference_genes = sorted(gene.id for gene in reaction.genes)
        rewritten_rule = _rewrite_gpr_for_supported_hits(reaction, gene_hits)
        if rewritten_rule is None:
            unsupported.append(reaction_id)
            status = "unsupported"
        else:
            reaction.gene_reaction_rule = rewritten_rule
            supported.append(reaction_id)
            status = "sequence_supported"
        reaction_audit[reaction_id] = {
            "status": status,
            "original_gpr": original_rule,
            "final_gpr": (
                reaction.gene_reaction_rule
                if status == "sequence_supported"
                else None
            ),
            "reference_genes": reference_genes,
            "supported_reference_genes": sorted(
                gene_id for gene_id in reference_genes if gene_id in gene_hits
            ),
        }
    return supported, unsupported, reaction_audit


def _prune_unsupported_reactions(
    model,
    reaction_ids,
    *,
    objective_fraction=0.5,
    objective_cap=0.1,
):
    """Retain the smallest unsupported-reaction set needed to preserve growth."""

    baseline = model.optimize()
    baseline_objective = _require_optimal_finite_objective(
        baseline,
        "unsupported-reaction pruning baseline",
    )
    threshold = objective_fraction * min(
        baseline_objective,
        objective_cap,
    )
    try:
        feasibility_tolerance = float(
            model.solver.configuration.tolerances.feasibility
        )
    except (AttributeError, TypeError, ValueError) as error:
        raise RuntimeError(
            "The selected solver does not expose a numeric feasibility "
            "tolerance required for auditable strict growth retention"
        ) from error
    pruning_solver = _solver_tolerance_snapshot(model)
    growth_margin = max(
        100.0 * feasibility_tolerance,
        1e-6,
        abs(threshold) * 1e-6,
    )
    growth_floor = threshold + growth_margin
    candidate_ids = [
        reaction_id
        for reaction_id in sorted(set(reaction_ids))
        if reaction_id in model.reactions
    ]
    forced_flux_candidates = [
        reaction_id
        for reaction_id in candidate_ids
        if not (
            model.reactions.get_by_id(reaction_id).lower_bound
            <= 0
            <= model.reactions.get_by_id(reaction_id).upper_bound
        )
    ]
    if forced_flux_candidates:
        raise RuntimeError(
            "Unsupported-reaction minimization requires every candidate to "
            "permit zero flux so reaction addition is monotonic; forced-flux "
            f"candidates={forced_flux_candidates[:20]}"
        )
    individual_trials = {}
    for reaction_id in candidate_ids:
        with model:
            model.remove_reactions([reaction_id], remove_orphans=True)
            solution = model.optimize()
            objective = (
                float(solution.objective_value)
                if solution.objective_value is not None
                else None
            )
        individual_trials[reaction_id] = {
            "status_after_individual_removal": solution.status,
            "objective_after_individual_removal": objective,
            "individually_growth_required": not (
                solution.status == "optimal"
                and objective is not None
                and objective + feasibility_tolerance >= growth_floor
            ),
        }

    if not candidate_ids:
        return {
            "growth_retention_rule": {
                "objective_fraction": objective_fraction,
                "objective_cap": objective_cap,
                "absolute_threshold": threshold,
                "solver_feasibility_tolerance": feasibility_tolerance,
                "solver": pruning_solver,
                "strict_growth_margin": growth_margin,
                "required_growth_floor": growth_floor,
                "selection": (
                    "minimum_cardinality_unsupported_reaction_set_milp"
                ),
                "source": (
                    "operationalization_of_paper_seriously_weaken_using_"
                    "released_code_thresholds"
                ),
            },
            "removed_unsupported_reactions": [],
            "retained_algorithmic_reactions": [],
            "minimum_retained_reaction_count": 0,
            "minimum_set_identity_is_unique": True,
            "final_growth_status": baseline.status,
            "final_growth_objective": baseline_objective,
            "removal_trials": {},
        }

    optimization_model = model.copy()
    interface = optimization_model.problem
    binary_variables = {}
    linking_constraints = []
    for index, reaction_id in enumerate(candidate_ids):
        reaction = optimization_model.reactions.get_by_id(reaction_id)
        keep = interface.Variable(
            f"alphagem_keep_gap_{index}",
            type="binary",
        )
        binary_variables[reaction_id] = keep
        for direction, flux_variable in (
            ("forward", reaction.forward_variable),
            ("reverse", reaction.reverse_variable),
        ):
            capacity = flux_variable.ub
            if capacity is None or not math.isfinite(float(capacity)):
                raise RuntimeError(
                    "Paper gap-reaction minimization requires finite reaction "
                    f"bounds: {reaction_id} {direction} upper bound={capacity}"
                )
            linking_constraints.append(
                interface.Constraint(
                    flux_variable - float(capacity) * keep,
                    ub=0,
                    name=f"alphagem_link_gap_{index}_{direction}",
                )
            )
    optimization_model.add_cons_vars(
        list(binary_variables.values()) + linking_constraints
    )
    growth_constraint = interface.Constraint(
        optimization_model.objective.expression,
        lb=growth_floor,
        name="alphagem_gapfill_growth_floor",
    )
    optimization_model.add_cons_vars([growth_constraint])
    optimization_model.objective = interface.Objective(
        sum(binary_variables.values(), 0),
        direction="min",
    )
    validation_cut_counter = [0]

    def optimize_with_direct_growth_validation():
        seen_selected_sets = set()
        while True:
            solution = optimization_model.optimize()
            if solution.status == "infeasible":
                return solution, None
            if solution.status != "optimal":
                raise RuntimeError(
                    "Unsupported-reaction MILP returned a non-decisive solver "
                    f"status: {solution.status}"
                )
            selected_ids = {
                reaction_id
                for reaction_id in candidate_ids
                if float(binary_variables[reaction_id].primal or 0.0) > 0.5
            }
            selected_key = tuple(sorted(selected_ids))
            if selected_key in seen_selected_sets:
                raise RuntimeError(
                    "Unsupported-reaction MILP repeated a rejected binary "
                    f"assignment despite a no-good cut: {selected_key}"
                )
            seen_selected_sets.add(selected_key)
            with model:
                excluded_ids = [
                    reaction_id
                    for reaction_id in candidate_ids
                    if reaction_id not in selected_ids
                ]
                if excluded_ids:
                    model.remove_reactions(
                        excluded_ids,
                        remove_orphans=True,
                    )
                direct_solution = model.optimize()
                direct_objective = _finite_objective_or_none(
                    direct_solution.objective_value
                )
            if direct_solution.status not in {"optimal", "infeasible"}:
                raise RuntimeError(
                    "Direct FBA validation returned a non-decisive solver "
                    f"status: {direct_solution.status}"
                )
            if (
                direct_solution.status == "optimal"
                and direct_objective is None
            ):
                raise RuntimeError(
                    "Direct FBA validation returned an optimal status without "
                    "an objective value"
                )
            directly_valid = (
                direct_solution.status == "optimal"
                and direct_objective is not None
                and direct_objective + feasibility_tolerance >= growth_floor
                and direct_objective > threshold
            )
            if directly_valid:
                return solution, {
                    "retained_ids": selected_ids,
                    "status": direct_solution.status,
                    "objective": direct_objective,
                }
            excluded_variables = [
                binary_variables[reaction_id]
                for reaction_id in candidate_ids
                if reaction_id not in selected_ids
            ]
            if not excluded_variables:
                raise RuntimeError(
                    "Even retaining every unsupported reaction failed the "
                    "direct growth-floor validation"
                )
            cut_index = validation_cut_counter[0]
            validation_cut_counter[0] += 1
            optimization_model.add_cons_vars(
                [
                    interface.Constraint(
                        sum(excluded_variables, 0),
                        lb=1,
                        name=f"alphagem_direct_growth_cut_{cut_index}",
                    )
                ]
            )

    minimum_solution, minimum_direct = (
        optimize_with_direct_growth_validation()
    )
    if minimum_solution.status != "optimal":
        raise RuntimeError(
            "Could not solve the minimum unsupported-reaction retention "
            f"problem: {minimum_solution.status}"
        )
    minimum_count = int(
        round(
            _require_optimal_finite_objective(
                minimum_solution,
                "minimum unsupported-reaction retention MILP",
            )
        )
    )
    cardinality_expression = sum(binary_variables.values(), 0)
    minimum_cardinality = interface.Constraint(
        cardinality_expression,
        lb=minimum_count,
        ub=minimum_count,
        name="alphagem_fix_minimum_gap_cardinality",
    )
    optimization_model.add_cons_vars([minimum_cardinality])
    identity_audit = {}
    for index, reaction_id in enumerate(candidate_ids):
        with optimization_model:
            optimization_model.add_cons_vars(
                [
                    interface.Constraint(
                        binary_variables[reaction_id],
                        ub=0,
                        name=f"alphagem_test_gap_absent_{index}",
                    ),
                ]
            )
            absent_solution, absent_direct = (
                optimize_with_direct_growth_validation()
            )
        with optimization_model:
            optimization_model.add_cons_vars(
                [
                    interface.Constraint(
                        binary_variables[reaction_id],
                        lb=1,
                        name=f"alphagem_test_gap_present_{index}",
                    ),
                ]
            )
            present_solution, present_direct = (
                optimize_with_direct_growth_validation()
            )
        can_be_absent = (
            absent_solution.status == "optimal"
            and absent_direct is not None
        )
        can_be_present = (
            present_solution.status == "optimal"
            and present_direct is not None
        )
        if not can_be_absent and not can_be_present:
            raise RuntimeError(
                "Minimum-set identity audit found no feasible state for "
                f"reaction {reaction_id}"
            )
        identity_audit[reaction_id] = {
            "present_in_some_minimum_set": can_be_present,
            "present_in_all_minimum_sets": not can_be_absent,
            "absent_in_some_minimum_set": can_be_absent,
            "absent_in_all_minimum_sets": not can_be_present,
            "minimum_set_membership_interchangeable": (
                can_be_absent and can_be_present
            ),
        }

    canonical_constraints = []
    for index, reaction_id in enumerate(candidate_ids):
        with optimization_model:
            optimization_model.add_cons_vars(
                [
                    interface.Constraint(
                        binary_variables[reaction_id],
                        ub=0,
                        name=f"alphagem_try_canonical_absent_{index}",
                    )
                ]
            )
            canonical_trial, canonical_direct = (
                optimize_with_direct_growth_validation()
            )
        canonical_value = (
            0
            if canonical_trial.status == "optimal"
            and canonical_direct is not None
            else 1
        )
        fixed = interface.Constraint(
            binary_variables[reaction_id],
            lb=canonical_value,
            ub=canonical_value,
            name=f"alphagem_canonical_gap_{index}",
        )
        optimization_model.add_cons_vars([fixed])
        canonical_constraints.append(fixed)
    canonical_solution, canonical_direct = (
        optimize_with_direct_growth_validation()
    )
    if canonical_solution.status != "optimal":
        raise RuntimeError(
            "Could not construct the deterministic canonical minimum "
            f"unsupported-reaction set: {canonical_solution.status}"
        )
    if canonical_direct is None:
        raise RuntimeError(
            "Canonical minimum unsupported-reaction set was not directly "
            "validated against the original model"
        )
    retained = [
        reaction_id
        for reaction_id in candidate_ids
        if float(binary_variables[reaction_id].primal or 0.0) > 0.5
    ]
    if len(retained) != minimum_count:
        raise RuntimeError(
            "MILP returned an inconsistent unsupported-reaction count: "
            f"objective={minimum_count}, selected={len(retained)}"
        )
    removed = [
        reaction_id
        for reaction_id in candidate_ids
        if reaction_id not in retained
    ]
    if removed:
        model.remove_reactions(removed, remove_orphans=True)
    for reaction_id in retained:
        model.reactions.get_by_id(reaction_id).gene_reaction_rule = ""
    final_solution = model.optimize()
    final_objective = _finite_objective_or_none(
        final_solution.objective_value
    )
    if (
        final_solution.status != "optimal"
        or final_objective is None
        or final_objective + feasibility_tolerance < growth_floor
        or final_objective <= threshold
    ):
        raise RuntimeError(
            "Minimum unsupported-reaction set failed the post-removal growth "
            f"check: status={final_solution.status}, objective={final_objective}"
        )

    decisions = {}
    for reaction_id in candidate_ids:
        membership = identity_audit[reaction_id]
        if reaction_id not in retained:
            decision = (
                "removed_canonical_interchangeable_member"
                if membership["minimum_set_membership_interchangeable"]
                else "removed_from_all_minimum_sets"
            )
        elif membership["present_in_all_minimum_sets"]:
            decision = "retained_in_all_minimum_sets"
        else:
            decision = "retained_canonical_interchangeable_member"
        decisions[reaction_id] = {
            **individual_trials[reaction_id],
            **membership,
            "retention_threshold": threshold,
            "decision": decision,
        }
    return {
        "growth_retention_rule": {
            "objective_fraction": objective_fraction,
            "objective_cap": objective_cap,
            "absolute_threshold": threshold,
            "solver_feasibility_tolerance": feasibility_tolerance,
            "solver": pruning_solver,
            "strict_growth_margin": growth_margin,
            "required_growth_floor": growth_floor,
            "selection": "minimum_cardinality_unsupported_reaction_set_milp",
            "direct_fba_validation_cuts": validation_cut_counter[0],
            "canonical_tie_break": (
                "sorted_reaction_ids_prefer_absent_when_feasible"
            ),
            "source": (
                "operationalization_of_paper_seriously_weaken_using_"
                "released_code_thresholds"
            ),
        },
        "removed_unsupported_reactions": removed,
        "retained_algorithmic_reactions": retained,
        "minimum_retained_reaction_count": minimum_count,
        "minimum_set_identity_is_unique": not any(
            value["minimum_set_membership_interchangeable"]
            for value in identity_audit.values()
        ),
        "final_growth_status": final_solution.status,
        "final_growth_objective": final_objective,
        "removal_trials": decisions,
    }


def read_gapgenes(
    model,
    genes,
    solution,
    refname,
    name,
    *,
    bitscore_threshold,
    bitscore_operator,
    search_backend="legacy_auto",
    coordinate_rows=None,
    tie_break="strict",
    unmapped_hit_policy="error",
):
    findrefnames.predata(refname)
    rows = _read_blast_rows(f"working/{name}/outfile{name}.csv")
    tie_audit = []
    unmapped_audit = []
    if search_backend == "tblastn_genome":
        hits = map_tblastn_hits_to_target_proteins(
            rows,
            coordinate_rows or [],
            threshold=bitscore_threshold,
            operator=bitscore_operator,
            tie_break=tie_break,
            tie_audit=tie_audit,
            unmapped_hit_policy=unmapped_hit_policy,
            unmapped_audit=unmapped_audit,
        )
    else:
        hits = best_gapfill_hits(
            rows,
            threshold=bitscore_threshold,
            operator=bitscore_operator,
        )
    nonfindgenes = []
    mapped = 0
    query_to_model_gene = {}
    for gene in genes:
        source_gene = findrefnames.findentryname(protein_key(gene))
        if search_backend == "tblastn_genome" and not source_gene:
            raise KeyError(
                f"Reference protein is absent from the model-name mapping: {gene}"
            )
        query_to_model_gene[gene] = source_gene
        hit = hits.get(gene)
        if hit is None:
            nonfindgenes.append(source_gene)
            continue
        if search_backend == "tblastn_genome":
            mapped += 1
            continue
        new_gene = protein_key(hit["subject"])
        for reaction_id in solution:
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.gene_reaction_rule = _replace_gpr_gene(
                reaction.gene_reaction_rule,
                source_gene,
                new_gene,
            )
        mapped += 1
    support_audit = {}
    unsupported_reactions = []
    if search_backend == "tblastn_genome":
        (
            supported_reactions,
            unsupported_reactions,
            reaction_audit,
        ) = _apply_paper_reaction_support(
            model,
            solution,
            query_to_model_gene,
            hits,
        )
        support_audit = {
            "sequence_supported_reactions": supported_reactions,
            "sequence_unsupported_reactions": unsupported_reactions,
            "reaction_sequence_support": reaction_audit,
            "selected_sequence_hits": {
                query: {
                    key: value
                    for key, value in hit.items()
                    if key
                    in {
                        "subject",
                        "genome_subject",
                        "sstart",
                        "send",
                        "evalue",
                        "bitscore",
                        "coordinate_evidence",
                        "coordinate_overlap_bp",
                    }
                }
                for query, hit in sorted(hits.items())
            },
        }
    return model, nonfindgenes, unsupported_reactions, {
        "query_count": len(genes),
        "raw_hit_count": len(rows),
        "retained_query_count": len(hits),
        "mapped_query_count": mapped,
        "unmapped_query_count": len(nonfindgenes),
        "hit_tie_break": tie_break,
        "hit_tie_breaks": tie_audit,
        "hit_unmapped_policy": unmapped_hit_policy,
        "hit_unmapped": unmapped_audit,
        **support_audit,
    }


def _run_gapfill_sequence_search(
    *,
    name,
    refname,
    search_backend,
    gapfill_cds,
    gapfill_genome,
    gapfill_coordinate_map,
    gapfill_coordinate_segments,
    gapfill_coordinate_contract,
):
    working = Path("working") / name
    query = working / f"{refname}_gap.fasta"
    output = working / f"outfile{name}.csv"
    database = working / f"{name}_gapfill_db"

    if search_backend == "tblastn_genome":
        if (
            gapfill_genome is None
            or gapfill_coordinate_map is None
            or gapfill_coordinate_segments is None
        ):
            raise ValueError(
                "tblastn_genome requires a raw genome, coordinate envelope "
                "map, and CDS-segment map"
            )
        genome_path = Path(gapfill_genome)
        contract = validate_gapfill_genome(
            genome_path,
            working / f"{name}.fasta",
            Path(gapfill_coordinate_map),
            Path(gapfill_coordinate_segments),
            translation_exceptions=(
                load_coordinate_translation_exceptions(
                    Path(gapfill_coordinate_contract)
                )
                if gapfill_coordinate_contract is not None
                else None
            ),
        )
        makeblastdb_command = [
            "makeblastdb",
            "-in",
            str(genome_path),
            "-dbtype",
            "nucl",
            "-input_type",
            "fasta",
            "-out",
            str(database),
        ]
        tblastn_command = [
            "tblastn",
            "-query",
            str(query),
            "-db",
            str(database),
            "-out",
            str(output),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend "
            "sstart send evalue bitscore",
            "-evalue",
            "1e-1",
            "-max_target_seqs",
            "1000000",
        ]
        subprocess.run(makeblastdb_command, check=True)
        subprocess.run(tblastn_command, check=True)
        execution = {
            "backend": "tblastn_genome",
            "commands": [makeblastdb_command, tblastn_command],
            "inputs": {
                "query": {
                    "path": str(query.resolve()),
                    "sha256": file_sha256(query),
                },
                "raw_genome": {
                    "path": str(genome_path.resolve()),
                    "sha256": file_sha256(genome_path),
                },
                "coordinate_map": {
                    "path": str(Path(gapfill_coordinate_map).resolve()),
                    "sha256": file_sha256(Path(gapfill_coordinate_map)),
                },
                "coordinate_segments": {
                    "path": str(
                        Path(gapfill_coordinate_segments).resolve()
                    ),
                    "sha256": file_sha256(
                        Path(gapfill_coordinate_segments)
                    ),
                },
            },
            "raw_output": {
                "path": str(output.resolve()),
                "bytes": output.stat().st_size,
                "sha256": file_sha256(output),
            },
        }
        execution_path = working / "gapfill_sequence_search.json"
        execution_path.write_text(
            json.dumps(execution, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return "tblastn", contract["segment_rows"], execution

    if search_backend == "tblastn_cds":
        if gapfill_cds is None:
            raise ValueError("tblastn_cds requires an explicit nucleotide CDS FASTA")
        cds_path = Path(gapfill_cds)
        validate_gapfill_cds(cds_path, working / f"{name}.fasta")
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(cds_path),
                "-dbtype",
                "nucl",
                "-input_type",
                "fasta",
                "-out",
                str(database),
            ],
            check=True,
        )
        subprocess.run(
            [
                "tblastn",
                "-query",
                str(query),
                "-db",
                str(database),
                "-out",
                str(output),
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore",
                "-evalue",
                "1e-1",
                "-max_target_seqs",
                "50",
            ],
            check=True,
        )
        return "tblastn", None, None

    if search_backend != "legacy_auto":
        raise ValueError(f"Unsupported gap-fill search backend: {search_backend}")
    legacy_nucleotide = working / f"{name}_non_anno.fasta"
    if legacy_nucleotide.is_file():
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(legacy_nucleotide),
                "-dbtype",
                "nucl",
                "-input_type",
                "fasta",
                "-out",
                str(database),
            ],
            check=True,
        )
        subprocess.run(
            [
                "tblastn",
                "-query",
                str(query),
                "-db",
                str(database),
                "-out",
                str(output),
                "-outfmt",
                "6",
                "-evalue",
                "1e-1",
                "-max_target_seqs",
                "1",
            ],
            check=True,
        )
        return "tblastn", None, None

    subprocess.run(
        [
            "diamond",
            "makedb",
            "--in",
            str(working / f"{name}.fasta"),
            "--db",
            str(database),
        ],
        check=True,
    )
    subprocess.run(
        [
            "diamond",
            "blastp",
            "--query",
            str(query),
            "--db",
            str(database),
            "--out",
            str(output),
            "--outfmt",
            "6",
        ],
        check=True,
    )
    return "diamond_blastp", None, None

def _export_final_models(
    name,
    model,
    *,
    require_optimal_solution=False,
    minimum_objective=1e-5,
    audit_factory=None,
    normalize_model=True,
):
    working = Path("working") / name
    canonical_xml = working / f"{name}-GEM.xml"
    canonical_json = working / f"{name}-GEM.json"
    canonical_audit = working / "gapfill_audit.json"
    if require_optimal_solution:
        for canonical in (
            canonical_xml,
            canonical_json,
            canonical_audit,
        ):
            if canonical.is_symlink() or (
                canonical.exists() and not canonical.is_file()
            ):
                raise RuntimeError(
                    "Formal gap-fill model output must be a regular file: "
                    f"{canonical}"
                )
    with tempfile.TemporaryDirectory(
        dir=working,
        prefix=".gapfill-export-",
    ) as temporary_directory:
        temporary_root = Path(temporary_directory)
        temporary_xml = temporary_root / canonical_xml.name
        temporary_json = temporary_root / canonical_json.name
        cobra.io.write_sbml_model(model, temporary_xml)
        cobra.io.save_json_model(model, temporary_json)
        if normalize_model:
            exported_model = cobra.io.read_sbml_model(temporary_xml)
            exported_model.id = f"{name}"
            exported_model.name = f"{name}-GEM"
            for gene in exported_model.genes:
                gene.annotation["uniprot"] = gene.id
            cobra.io.write_sbml_model(exported_model, temporary_xml)
            cobra.io.save_json_model(exported_model, temporary_json)
        final_xml_model = cobra.io.read_sbml_model(temporary_xml)
        final_json_model = cobra.io.load_json_model(temporary_json)
        format_issues = _model_equivalence_issues(
            final_xml_model,
            final_json_model,
        )
        if format_issues:
            raise RuntimeError(
                "Final XML and JSON models are not equivalent: "
                + "; ".join(format_issues[:20])
            )
        final_xml_solution = final_xml_model.optimize()
        final_json_solution = final_json_model.optimize()
        if require_optimal_solution:
            for output_format, solution in (
                ("XML", final_xml_solution),
                ("JSON", final_json_solution),
            ):
                _require_optimal_finite_objective(
                    solution,
                    f"round-tripped final gap-fill {output_format} model",
                    minimum_exclusive=minimum_objective,
                )
        artifacts = [
            (temporary_xml, canonical_xml),
            (temporary_json, canonical_json),
        ]
        if audit_factory is not None:
            audit_payload = audit_factory(
                final_xml_model,
                final_xml_solution,
            )
            if not isinstance(audit_payload, dict):
                raise RuntimeError(
                    "Gap-fill audit factory must return a dictionary"
                )
            audit_payload = dict(audit_payload)
            audit_payload["format_validation"] = {
                "formats_equivalent": True,
                "xml": _solution_record(final_xml_solution),
                "json": _solution_record(final_json_solution),
            }
            audit_payload["model_outputs"] = {
                "xml": {
                    "bytes": temporary_xml.stat().st_size,
                    "sha256": file_sha256(temporary_xml),
                },
                "json": {
                    "bytes": temporary_json.stat().st_size,
                    "sha256": file_sha256(temporary_json),
                },
            }
            temporary_audit = temporary_root / canonical_audit.name
            with temporary_audit.open("w", encoding="utf-8") as handle:
                json.dump(audit_payload, handle, indent=2, sort_keys=True)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            artifacts.append((temporary_audit, canonical_audit))
        backups = {}
        for _, canonical in artifacts:
            if canonical.exists():
                backup = temporary_root / f"prior-{canonical.name}"
                shutil.copy2(canonical, backup)
                backups[canonical] = backup
        published = []
        try:
            for temporary, canonical in artifacts:
                os.replace(temporary, canonical)
                published.append(canonical)
        except Exception:
            for canonical in reversed(published):
                backup = backups.get(canonical)
                if backup is not None:
                    os.replace(backup, canonical)
                elif canonical.exists():
                    canonical.unlink()
            raise
    return final_xml_model, final_xml_solution


def _solver_tolerance_snapshot(model):
    observed = {}
    for name in ("feasibility", "integrality", "optimality"):
        try:
            observed[name] = float(
                getattr(model.solver.configuration.tolerances, name)
            )
        except (AttributeError, TypeError, ValueError):
            observed[name] = None
    return {
        "interface": model.solver.interface.__name__,
        "tolerances": observed,
    }


def _set_solver_tolerances(model, requested, *, require_all=False):
    unsupported = []
    for name, value in requested.items():
        try:
            setattr(model.solver.configuration.tolerances, name, value)
        except (AttributeError, TypeError, ValueError):
            unsupported.append(name)
    record = {
        **_solver_tolerance_snapshot(model),
        "requested_tolerances": dict(requested),
        "unsupported_tolerances": unsupported,
    }
    mismatched = {
        name: {
            "requested": float(value),
            "observed": record["tolerances"].get(name),
        }
        for name, value in requested.items()
        if (
            name not in unsupported
            and (
                record["tolerances"].get(name) is None
                or not math.isclose(
                    record["tolerances"][name],
                    float(value),
                    rel_tol=1e-12,
                    abs_tol=0.0,
                )
            )
        )
    }
    record["mismatched_tolerances"] = mismatched
    if require_all and (unsupported or mismatched):
        raise RuntimeError(
            "Selected solver cannot enforce all formal gap-fill tolerances: "
            f"interface={record['interface']}, unsupported={unsupported}, "
            f"mismatched={mismatched}"
        )
    return record


def _gapfill_tolerance_request(model):
    """Return the enforceable formal tolerance contract for this backend.

    The paper pins COBRApy but does not declare a solver or a common set of
    backend tolerance names. GLPK/optlang exposes feasibility and integrality
    but no optimality setting, so requesting the latter would make the formal
    run fail for a parameter that the selected backend cannot represent.
    """
    interface = str(model.solver.interface.__name__).lower()
    requested = {
        "feasibility": 1e-9,
        "integrality": 1e-9,
    }
    if "glpk" not in interface:
        requested["optimality"] = 1e-9
    return requested


def gapfill(
    name,
    refname,
    refmodelname,
    grothmedium='min',
    esstinal=None,
    strategy="selected_medium",
    search_backend="legacy_auto",
    gapfill_cds=None,
    gapfill_genome=None,
    gapfill_coordinate_map=None,
    gapfill_coordinate_segments=None,
    gapfill_coordinate_contract=None,
    bitscore_threshold=200.0,
    bitscore_operator="ge",
    tie_break="strict",
    unmapped_hit_policy="error",
    require_solver_tolerances=False,
    require_optimal_solutions=False,
    require_essential_reactions=False,
):
    if esstinal is None:
        esstinal=[]
    if require_optimal_solutions:
        prior_audit = Path("working") / name / "gapfill_audit.json"
        if prior_audit.is_symlink() or (
            prior_audit.exists() and not prior_audit.is_file()
        ):
            raise RuntimeError(
                "Formal gap-fill audit path must be a regular file: "
                f"{prior_audit}"
            )
    media_plan = gapfill_media_plan(
        strategy,
        grothmedium,
        refname not in {"human", "strco", "synechocystis"},
    )
    sequence_audit = {
        "search_backend": search_backend,
        "search_tool": None,
        "search_execution": None,
        "bitscore_threshold": float(bitscore_threshold),
        "bitscore_operator": bitscore_operator,
        "sequence_search_performed": False,
        "sequence_input": {
            key: value
            for key, value in {
                "raw_genome": (
                    {
                        "path": str(Path(gapfill_genome).resolve()),
                        "input_kind": "raw_genome",
                        "sha256": file_sha256(Path(gapfill_genome)),
                    }
                    if gapfill_genome is not None
                    else None
                ),
                "coordinate_map": (
                    {
                        "path": str(Path(gapfill_coordinate_map).resolve()),
                        "input_kind": "protein_to_genome_coordinate_map",
                        "sha256": file_sha256(Path(gapfill_coordinate_map)),
                    }
                    if gapfill_coordinate_map is not None
                    else None
                ),
                "coordinate_segments": (
                    {
                        "path": str(
                            Path(gapfill_coordinate_segments).resolve()
                        ),
                        "input_kind": "protein_to_genome_cds_segments",
                        "sha256": file_sha256(
                            Path(gapfill_coordinate_segments)
                        ),
                    }
                    if gapfill_coordinate_segments is not None
                    else None
                ),
                "coordinate_contract": (
                    {
                        "path": str(
                            Path(gapfill_coordinate_contract).resolve()
                        ),
                        "input_kind": "coordinate_recovery_contract",
                        "sha256": file_sha256(
                            Path(gapfill_coordinate_contract)
                        ),
                    }
                    if gapfill_coordinate_contract is not None
                    else None
                ),
                "nucleotide_cds_proxy": (
                    {
                        "path": str(Path(gapfill_cds).resolve()),
                        "input_kind": "nucleotide_cds",
                        "sha256": file_sha256(Path(gapfill_cds)),
                    }
                    if gapfill_cds is not None
                    else None
                ),
            }.items()
            if value is not None
        } or None,
        "query_count": 0,
        "raw_hit_count": 0,
        "retained_query_count": 0,
        "mapped_query_count": 0,
        "unmapped_query_count": 0,
        "original_model_solver": None,
        "gapfiller_solvers": [],
        "solver_check_policy": (
            "optimal_finite_objective_required"
            if require_optimal_solutions
            else "released_code_compatibility"
        ),
        "solver_checks": {},
    }
    refmodel=cobra.io.read_sbml_model(f'models/{refmodelname}')
    # reactionpool=cobra.io.load_json_model('./models/reaction_pool_final.json')
    findrefnames.predata(refname)
    model = cobra.io.load_yaml_model(f'working/{name}/{name}-GEM_withgaps.yml')
    sequence_audit["original_model_solver"] = _solver_tolerance_snapshot(
        model
    )
    if refname == 'human':
        model.objective=model.reactions.MAR00021
        model.reactions.MAR00021.bounds=[0,1000]
    minimal_medium = refmodel.medium.copy()
    if strategy == "selected_medium" and grothmedium == "full":
        full_medium_path = Path(f"data_available/{refname}_full_medium.pkl")
        if not full_medium_path.is_file():
            raise FileNotFoundError(
                f"Full medium requested but not found: {full_medium_path}"
            )
        with full_medium_path.open("rb") as file:
            selected_medium = pickle.load(file)
        refmodel.medium = selected_medium
        model.medium = selected_medium
    else:
        selected_medium = minimal_medium
        model.medium = selected_medium
    initial_solution = model.optimize()
    sequence_audit["solver_checks"]["initial_model"] = _solution_record(
        initial_solution
    )
    initial_objective = _finite_objective_or_none(
        initial_solution.objective_value
    )
    if require_optimal_solutions:
        initial_objective = _require_optimal_finite_objective(
            initial_solution,
            "initial target model",
        )
    if (
        initial_solution.status == "optimal"
        and initial_objective is not None
        and initial_objective > 0.001
        and (not require_essential_reactions or not esstinal)
    ):
        def early_audit_factory(_exported_model, exported_solution):
            exported_objective = _finite_objective_or_none(
                exported_solution.objective_value
            )
            sequence_audit["solver_checks"]["final_roundtrip"] = (
                _solution_record(exported_solution)
            )
            return {
                "gapfill_exit_reason": "initial_model_already_grows",
                "strategy": strategy,
                "growth_medium": grothmedium,
                "media_plan": list(media_plan),
                "solver_status": exported_solution.status,
                "objective_value": exported_objective,
                "candidate_reactions": [],
                "essential_reactions": list(esstinal),
                "net_added_reactions": [],
                **sequence_audit,
            }

        _export_final_models(
            name,
            model,
            require_optimal_solution=require_optimal_solutions,
            minimum_objective=0.001,
            audit_factory=early_audit_factory,
            normalize_model=require_optimal_solutions,
        )
        return 0
    # gap = gapfilling.GapFiller(model, universal=refmodel, integer_threshold=1e-12, demand_reactions=False,lower_bound=0.05)
    # for j in gap.indicators:
    #     i=gap.model.variables.get(key=j.rxn_id)
    #     if i._get_primal() >0:
    #         solution.append(i.name)
    # refmodel_reactions=[rrr.id for rrr in refmodel.reactions]
    # for reac in model.reactions:
    #     if reac.id not in refmodel_reactions:
    #       refmodel.add_reactions([reac])
    # fba=cobra.flux_analysis.pfba(refmodel)
    # model2 = cobra.io.load_yaml_model(f'working/{name}/{name}-GEM_withgaps.yml')
    model2=model.copy()
    pre_reaction_ids = {reaction.id for reaction in model2.reactions}
    targenes=[]
    allgenes=[x for x in SeqIO.parse(open(f'data_available/{refname}.fasta'),'fasta')]#打开fasta文件
    allgenes_by_id = {protein_key(record.id): record for record in allgenes}
    reacgaps=[]
    # for solution1 in solution:
    #     if fba.fluxes.get(solution1)>=1e-10 or fba.fluxes.get(solution1)<=-1e-10:
    #         reacgaps.append(gap.model.reactions.get_by_id(solution1))
    reference_solution = refmodel.optimize()
    sequence_audit["solver_checks"]["reference_model"] = _solution_record(
        reference_solution
    )
    objectbound = reference_solution.objective_value
    if require_optimal_solutions:
        objectbound = _require_optimal_finite_objective(
            reference_solution,
            "reference model gap-fill bound",
            minimum_exclusive=0.0,
        )
    gapfiller = gapfilling.GapFiller(model2, refmodel, integer_threshold=1e-09, demand_reactions=False,lower_bound=min(objectbound*0.8,0.09))
    sequence_audit["gapfiller_solvers"].append(
        {
            "medium": grothmedium,
            **_set_solver_tolerances(
                gapfiller.model,
                _gapfill_tolerance_request(gapfiller.model),
                require_all=require_solver_tolerances,
            ),
        }
    )
    try:
        reacgaps= gapfiller.fill()[0]
    except Exception as error:
        if strategy == "selected_medium":
            raise RuntimeError(
                f"Gap-filling failed on {grothmedium} medium"
            ) from error
        reacgaps=[]
        fba=gapfiller.model.optimize()
        for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            if fba.fluxes[i.id] != 0:
                try:
                    reacgaps.append(refmodel.reactions.get_by_id(i.id))
                except:
                    continue
    for r in esstinal:
        try:
            reacgaps.append(refmodel.reactions.get_by_id(r))
        except KeyError as error:
            if strategy == "selected_medium":
                raise KeyError(
                    f"Essential reaction is absent from reference model: {r}"
                ) from error
            print('Can not find essential reaction {r}',r)
            continue
    for reac in reacgaps:
        for genes in reac.genes:
            reference_entry = findrefnames.findmodelname(genes.id)
            if not reference_entry or reference_entry not in allgenes_by_id:
                if search_backend == "tblastn_genome":
                    raise KeyError(
                        "Gap reaction gene cannot be mapped to a reference "
                        f"protein sequence: reaction={reac.id}, gene={genes.id}"
                    )
                continue
            targenes.append(allgenes_by_id[reference_entry])
    medium = minimal_medium
    if (
        media_plan == ("min", "full")
    ):
        with open(f'data_available/{refname}_full_medium.pkl', 'rb') as file:
            fullmedium = pickle.load(file)
        medium=refmodel.medium.copy()
        refmodel.medium = fullmedium
        objectbound = refmodel.optimize().objective_value
        fba2=cobra.flux_analysis.pfba(refmodel)
        reacgaps2 = []
        model2.medium=fullmedium
        # for solution1 in solution:
        #     if fba2.fluxes.get(solution1) >= 1e-10 or fba2.fluxes.get(solution1) <= -1e-10:
        #         reacgaps2.append(gap.model.reactions.get_by_id(solution1))
        gapfiller = gapfilling.GapFiller(model2, refmodel, integer_threshold=1e-09, demand_reactions=False,lower_bound=min(objectbound*0.5,0.05))
        sequence_audit["gapfiller_solvers"].append(
            {
                "medium": "full",
                **_set_solver_tolerances(
                    gapfiller.model,
                    _gapfill_tolerance_request(gapfiller.model),
                    require_all=require_solver_tolerances,
                ),
            }
        )
        try:
            reacgaps2 = gapfiller.fill()[0]
        except:
          reacgaps2 = []
          fba = gapfiller.model.optimize()
          for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            if fba.fluxes[i.id] != 0:
                try:
                    reacgaps2.append(refmodel.reactions.get_by_id(i.id))
                except:
                    continue
        if refname=='yeast':
           reacgaps2.append(refmodel.reactions.get_by_id('rhea_21624_c'))
        for reac in reacgaps2:
            for genes in reac.genes:
                reference_entry = findrefnames.findmodelname(genes.id)
                if not reference_entry or reference_entry not in allgenes_by_id:
                    continue
                targenes.append(allgenes_by_id[reference_entry])
        #pd.DataFrame(reacgaps).to_excel(f'working/{name}/reacgapsMM.xlsx')
        #pd.DataFrame(reacgaps2).to_excel(f'working/{name}/reacgaps2FULL.xlsx')
        reacgaps=reacgaps+reacgaps2
    reacgaps=[
        refmodel.reactions.get_by_id(reaction_id)
        for reaction_id in sorted({reaction.id for reaction in reacgaps})
    ]
    model2.add_reactions(reacgaps)
    solution_after_add = model2.optimize()
    sequence_audit["solver_checks"]["after_candidate_addition"] = (
        _solution_record(solution_after_add)
    )
    print(solution_after_add)
    solution_after_add_objective = _finite_objective_or_none(
        solution_after_add.objective_value
    )
    if require_optimal_solutions:
        solution_after_add_objective = _require_optimal_finite_objective(
            solution_after_add,
            "target model after candidate gap-reaction addition",
            minimum_exclusive=1e-5,
        )
    if (
        solution_after_add.status != "optimal"
        or solution_after_add_objective is None
        or solution_after_add_objective <= 1e-5
    ):
        if strategy == "selected_medium":
            raise RuntimeError(
                "Selected-medium gap-fill solution did not restore growth: "
                f"status={solution_after_add.status}, "
                f"objective={solution_after_add.objective_value}"
            )
        gapnow=[]
        for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            with model2:
                model2.add_reactions([i])
                if model2.optimize().objective_value>1e-4:
                    gapnow.append(i)
        # print(gapnow)
        for reac in gapnow:
            for genes in reac.genes:
                reference_entry = findrefnames.findmodelname(genes.id)
                if not reference_entry or reference_entry not in allgenes_by_id:
                    continue
                targenes.append(allgenes_by_id[reference_entry])
        model2.add_reactions(gapnow)
        print(f'Adding {len(gapnow)} genes in gapfilling')
    nonfindgenes=[]
    targenes = list({record.id: record for record in targenes}.values())
    if len(targenes)!=0:
       with open(f'./working/{name}/{refname}_gap.fasta','w') as fastafile:
           SeqIO.write(targenes,fastafile,'fasta')
       (
           sequence_audit["search_tool"],
           coordinate_rows,
           sequence_audit["search_execution"],
       ) = _run_gapfill_sequence_search(
           name=name,
           refname=refname,
           search_backend=search_backend,
           gapfill_cds=gapfill_cds,
           gapfill_genome=gapfill_genome,
           gapfill_coordinate_map=gapfill_coordinate_map,
           gapfill_coordinate_segments=gapfill_coordinate_segments,
           gapfill_coordinate_contract=gapfill_coordinate_contract,
       )
       sequence_audit["sequence_search_performed"] = True
       model2, nonfindgenes, unsupported_reactions, hit_audit = read_gapgenes(
           model=model2,
           genes=[x.id for x in targenes],
           solution=[reac.id for reac in reacgaps],
           refname=refname,
           name=name,
           bitscore_threshold=bitscore_threshold,
           bitscore_operator=bitscore_operator,
           search_backend=search_backend,
           coordinate_rows=coordinate_rows,
           tie_break=tie_break,
           unmapped_hit_policy=unmapped_hit_policy,
       )
       sequence_audit.update(hit_audit)
    else:
       unsupported_reactions = (
           [reaction.id for reaction in reacgaps]
           if search_backend == "tblastn_genome"
           else []
       )
    # print([reac.id for reac in reacgaps])
       # for met in model.metabolites:
       #     if met.compartment == '':
       #         model.metabolites.get_by_id(met.id).compartment = 'c'
       # for re in model.reactions:
       #     for cp in re.compartments:
       #         if cp == '':
       #             model.reactions.get_by_id(re.id).compartments.remove('')
       #             model.reactions.get_by_id(re.id).compartments.add('c')
    if strategy == "legacy_dual_medium" and grothmedium=='min':
        model2.medium=medium
    #changed
    # with model2:
    #     nonfindreactions=[r.id for r in cobra.manipulation.delete.knock_out_model_genes(model2, nonfindgenes)]
    #     pd.DataFrame(nonfindreactions).to_excel(f'working/{name}/reacgaps_non_gene.xlsx')
    paper_reaction_pruning = None
    if search_backend == "tblastn_genome":
        paper_reaction_pruning = _prune_unsupported_reactions(
            model2,
            unsupported_reactions,
        )
        sequence_audit.update(paper_reaction_pruning)
    genes_to_filter = (
        nonfindgenes if search_backend != "tblastn_genome" else []
    )
    filter1 = None
    if genes_to_filter:
        filter_solution = model2.optimize()
        sequence_audit["solver_checks"]["gene_filter_baseline"] = (
            _solution_record(filter_solution)
        )
        filter_objective = _finite_objective_or_none(
            filter_solution.objective_value
        )
        if require_optimal_solutions:
            filter_objective = _require_optimal_finite_objective(
                filter_solution,
                "gap-fill gene-filter baseline",
                minimum_exclusive=0.0,
            )
        if filter_objective is None:
            raise RuntimeError(
                "Gap-fill gene-filter baseline has no finite objective: "
                f"status={filter_solution.status}, "
                f"objective={filter_solution.objective_value}"
            )
        filter1 = min(filter_objective, 0.1)
    num_nonreactions=[]
    for gene in genes_to_filter:
        try:
            geneid=model2.genes.get_by_id(gene)
        except:
            continue
        with model2:
            cobra.manipulation.delete.knock_out_model_genes(model2,[geneid])
            try:
              if model2.optimize().objective_value<=0.5*filter1:#changed
                num_nonreactions=num_nonreactions+[r.id for r in geneid.reactions]
                continue
            except:
                num_nonreactions=num_nonreactions+[r.id for r in geneid.reactions]
                continue
        # print('delete gene ',gene)
        cobra.manipulation.delete.remove_genes(model2,[geneid])
    model2.compartments=refmodel.compartments
    essential_reaction_enforcement = _enforce_essential_reactions(
        model2,
        refmodel,
        esstinal,
        require_all=require_essential_reactions,
    )
    sequence_audit["essential_reaction_enforcement"] = (
        essential_reaction_enforcement
    )
    # if refname=='ecoli':
    #     cobra.io.save_json_model(model2,f'./working/{name}/{name}-GEM.json')
    #     model2 = cobra.io.load_json_model(f'./working/{name}/{name}-GEM.json')
    #     model2.id = f'{name}'
    #     model2.name = f'{name}-GEM'
    #     for i in model2.reactions:
    #         keys_to_remove = [k for k, v in i.annotation.items() if v in [[], '']]
    #         # 删除这些键
    #         for k in keys_to_remove:
    #             del i.annotation[k]
    #     cobra.io.save_json_model(model2,f'./working/{name}/{name}-GEM.json')
    #     cobra.io.write_sbml_model(model2, f'./working/{name}/{name}-GEM.xml')
    # else:
    for i in model2.reactions:
        keys_to_remove = [k for k, v in i.annotation.items() if v in [[], '']]
        # 删除这些键
        for k in keys_to_remove:
                del i.annotation[k]
    pre_export_solution = model2.optimize()
    sequence_audit["solver_checks"]["pre_export_model"] = _solution_record(
        pre_export_solution
    )
    if require_optimal_solutions:
        _require_optimal_finite_objective(
            pre_export_solution,
            "final gap-fill model before export",
            minimum_exclusive=1e-5,
        )
    print(pre_export_solution.objective_value)
    def completed_audit_factory(exported_model, exported_solution):
        sequence_audit["solver_checks"]["final_roundtrip"] = _solution_record(
            exported_solution
        )
        final_objective = _finite_objective_or_none(
            exported_solution.objective_value
        )
        if require_optimal_solutions:
            final_objective = _require_optimal_finite_objective(
                exported_solution,
                "round-tripped final gap-fill model",
                minimum_exclusive=1e-5,
            )
        if require_essential_reactions:
            missing_essential = sorted(
                reaction_id
                for reaction_id in set(esstinal)
                if reaction_id not in exported_model.reactions
            )
            if missing_essential:
                raise RuntimeError(
                    "Round-tripped final model is missing required essential "
                    f"reactions: {missing_essential}"
                )
        net_added_reactions = sorted(
            {reaction.id for reaction in exported_model.reactions}
            - pre_reaction_ids
        )
        return {
            "gapfill_exit_reason": "gapfill_completed",
            "strategy": strategy,
            "growth_medium": grothmedium,
            "media_plan": list(media_plan),
            "solver": str(exported_model.solver.interface),
            "solver_status": exported_solution.status,
            "objective_value": final_objective,
            "candidate_reactions": sorted(
                {reaction.id for reaction in reacgaps}
            ),
            "essential_reactions": sorted(set(esstinal)),
            "net_added_reactions": net_added_reactions,
            **sequence_audit,
        }

    model2, final_solution = _export_final_models(
        name,
        model2,
        require_optimal_solution=require_optimal_solutions,
        minimum_objective=1e-5,
        audit_factory=completed_audit_factory,
    )
    print(f'Number of all gapfilling reactions: {len(reacgaps)}')
    print(f'Number of non-bio genes: {len(set([r.id for r in reacgaps])&set(num_nonreactions))}')
