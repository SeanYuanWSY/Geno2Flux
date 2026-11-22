import ast
import cobra
import json
import math
import pandas as pd
import pickle
from pathlib import Path
from tqdm import tqdm

from cobra.core import GPR

from src.pipeline_contracts import reaction_is_mass_and_charge_balanced


def check_reactions(reac,model):
    del model
    try:
        imbalance = reac.check_mass_balance()
    except (TypeError, ValueError) as error:
        print(f"Rejecting {reac.id}: balance could not be checked: {error}")
        return 0
    if not reaction_is_mass_and_charge_balanced(imbalance):
        print(f"Rejecting {reac.id}: mass/charge imbalance {imbalance}")
        return 0
    return 1


def postmerge_metabolite_identity_conflicts(reaction, model):
    """Find same-ID metabolites that COBRA would replace on model insertion."""
    conflicts = []
    for candidate_metabolite in reaction.metabolites:
        try:
            model_metabolite = model.metabolites.get_by_id(candidate_metabolite.id)
        except KeyError:
            continue
        candidate_identity = {
            "formula": candidate_metabolite.formula,
            "charge": candidate_metabolite.charge,
            "compartment": candidate_metabolite.compartment,
        }
        model_identity = {
            "formula": model_metabolite.formula,
            "charge": model_metabolite.charge,
            "compartment": model_metabolite.compartment,
        }
        incomplete = any(
            value is None or value == ""
            for identity in (candidate_identity, model_identity)
            for value in identity.values()
        )
        if incomplete or candidate_identity != model_identity:
            conflicts.append({
                "metabolite_id": candidate_metabolite.id,
                "conflict_kind": (
                    "unverifiable_identity" if incomplete else "identity_mismatch"
                ),
                "candidate_identity": candidate_identity,
                "model_identity": model_identity,
            })
    return conflicts


ALREADY_PRESENT_POLICIES = (
    "author_skip",
    "merge_verified_identical_or",
    "merge_verified_identical_or_multisource",
    "merge_threeway_or_multisource",
)
CONFLICT_POLICIES = (
    "reject",
    "suffix_new_metabolite",
    "remap_by_model_identity",
)
# A candidate gene joins a merged GPR only when at least this many of the
# annotation columns anno_A/anno_B/anno_C/anno_D are non-empty for it.
MULTISOURCE_MIN_SOURCES = 2
MULTISOURCE_ANNOTATION_COLUMNS = ("anno_A", "anno_B", "anno_C", "anno_D")
# Coverage ratios above this share of input candidates are listed as
# escalation signals. They are signals only and never raise.
ESCALATION_RATIO_THRESHOLD = 0.10


def _author_name_key(name):
    """Author metdict equivalence: lowercased name with parentheses stripped."""
    return (name or "").lower().replace("(", "").replace(")", "")


def _source_name(source_label, metabolite_id, source_metabolites,
                 pool_raw_names):
    """Source-side metabolite name for the informational name comparison.

    Pool-sourced candidates use the raw frozen-pool JSON name (the in-memory
    pool reaction object may already carry the model registry name after a
    same-stage add unified the metabolite). Reference-sourced candidates use
    the reference reaction's metabolite object names.
    """
    if source_label == "frozen_reaction_pool" and pool_raw_names is not None:
        return pool_raw_names.get(metabolite_id)
    metabolite = source_metabolites.get(metabolite_id)
    return metabolite.name if metabolite is not None else None


def _normalized_metabolite_names(name):
    """Author metdict keys: lowercased name and its parenthesis-free form."""
    lowered = (name or "").lower()
    return {lowered, lowered.replace("(", "").replace(")", "")}


def _metabolite_identity(metabolite):
    return {
        "formula": metabolite.formula,
        "charge": metabolite.charge,
        "compartment": metabolite.compartment,
    }


def _identity_incomplete(identity):
    return any(value is None or value == "" for value in identity.values())


def _identity_matches_in_model(
    model, candidate_metabolite, exclude_id, allowed_ids=None
):
    """Model metabolites name- and identity-equal to the candidate side.

    allowed_ids pins the search to metabolites present at stage start so the
    result cannot depend on how many candidates ran before this one.
    """
    candidate_identity = _metabolite_identity(candidate_metabolite)
    if _identity_incomplete(candidate_identity):
        return []
    candidate_names = _normalized_metabolite_names(candidate_metabolite.name)
    if candidate_names == {""}:
        return []
    matches = []
    for metabolite in model.metabolites:
        if metabolite.id == exclude_id:
            continue
        if allowed_ids is not None and metabolite.id not in allowed_ids:
            continue
        if not _normalized_metabolite_names(metabolite.name) & candidate_names:
            continue
        if _metabolite_identity(metabolite) == candidate_identity:
            matches.append(metabolite)
    return matches


def _plan_conflict_resolution(
    reaction, model, conflicts, conflict_policy, stage_start_metabolite_ids=None
):
    """Plan all-or-nothing conflict repairs; None means fall back to reject."""
    candidate_metabolites = {m.id: m for m in reaction.metabolites}
    operations = []
    planned_ids = set()
    for conflict in conflicts:
        metabolite_id = conflict["metabolite_id"]
        candidate_metabolite = candidate_metabolites[metabolite_id]
        if conflict_policy == "suffix_new_metabolite":
            candidate_identity = conflict["candidate_identity"]
            if _identity_incomplete(candidate_identity):
                return None, {
                    "suffix_fallback": "incomplete_candidate_identity",
                    "fallback_metabolite_id": metabolite_id,
                }
            new_id = f"{metabolite_id}__pool"
            if new_id in planned_ids or new_id in candidate_metabolites:
                raise RuntimeError(f"Suffix metabolite ID collides: {new_id}")
            existing_suffix = None
            try:
                existing_suffix = model.metabolites.get_by_id(new_id)
            except KeyError:
                pass
            else:
                # A suffix metabolite from an earlier candidate of the same
                # pool is chemically identical and is reused; the same ID
                # with a divergent identity is a hard collision.
                if _metabolite_identity(existing_suffix) != candidate_identity:
                    raise RuntimeError(
                        f"Suffix metabolite ID collides: {new_id}"
                    )
            planned_ids.add(new_id)
            operations.append({
                "kind": "suffix",
                "metabolite_id": metabolite_id,
                "new_metabolite_id": new_id,
                "candidate_metabolite": candidate_metabolite,
                "existing_metabolite": existing_suffix,
            })
        elif conflict_policy == "remap_by_model_identity":
            matches = _identity_matches_in_model(
                model,
                candidate_metabolite,
                exclude_id=metabolite_id,
                allowed_ids=stage_start_metabolite_ids,
            )
            if len(matches) != 1:
                return None, {
                    "remap_fallback": "no_unique_match",
                    "fallback_metabolite_id": metabolite_id,
                    "match_count": len(matches),
                }
            operations.append({
                "kind": "remap",
                "metabolite_id": metabolite_id,
                "candidate_metabolite": candidate_metabolite,
                "target_metabolite": matches[0],
            })
        else:
            raise ValueError(f"Unknown conflict policy: {conflict_policy}")
    return operations, None


def _apply_conflict_resolution(reaction, operations, record):
    """Rewrite the candidate onto pool-suffixed or identity-matched metabolites."""
    suffix_records = []
    remap_records = []
    for operation in operations:
        candidate_metabolite = operation["candidate_metabolite"]
        stoich = reaction.metabolites[candidate_metabolite]
        if operation["kind"] == "suffix":
            suffix_metabolite = operation["existing_metabolite"]
            if suffix_metabolite is None:
                suffix_metabolite = cobra.Metabolite(
                    operation["new_metabolite_id"],
                    formula=candidate_metabolite.formula,
                    charge=candidate_metabolite.charge,
                    compartment=candidate_metabolite.compartment,
                    name=candidate_metabolite.name,
                )
            reaction.add_metabolites({suffix_metabolite: stoich})
            reaction.add_metabolites({candidate_metabolite: -1 * stoich})
            suffix_records.append({
                "metabolite_id": operation["metabolite_id"],
                "new_metabolite_id": operation["new_metabolite_id"],
                "name": candidate_metabolite.name,
                "reused_existing_pool_metabolite": (
                    operation["existing_metabolite"] is not None
                ),
                **_metabolite_identity(candidate_metabolite),
            })
        else:
            target = operation["target_metabolite"]
            reaction.add_metabolites({target: stoich})
            reaction.add_metabolites({candidate_metabolite: -1 * stoich})
            remap_records.append({
                "metabolite_id": operation["metabolite_id"],
                "remapped_to": target.id,
                "match_basis": "normalized_name_and_identity",
                **_metabolite_identity(candidate_metabolite),
            })
    if suffix_records:
        record["addition_kind"] = "suffix_metabolite"
        record["suffix_metabolites"] = suffix_records
    if remap_records:
        record["addition_kind"] = "remapped_metabolite"
        record["remapped_metabolites"] = remap_records


def admit_candidate_reaction(
    reaction, model, record, conflict_policy="reject",
    stage_start_metabolite_ids=None,
):
    """Conservatively reject a candidate before COBRA merges its metabolites."""
    if conflict_policy not in CONFLICT_POLICIES:
        raise ValueError(f"Unknown conflict policy: {conflict_policy}")
    # Preserve the author's first gate and its status for candidates that are
    # already imbalanced in their original source chemistry.
    if check_reactions(reaction, model) == 0:
        record["status"] = "rejected_mass_or_charge_balance"
        return
    conflicts = postmerge_metabolite_identity_conflicts(reaction, model)
    if conflicts and conflict_policy != "reject":
        operations, fallback = _plan_conflict_resolution(
            reaction, model, conflicts, conflict_policy,
            stage_start_metabolite_ids=stage_start_metabolite_ids,
        )
        if operations is None:
            record["status"] = "rejected_postmerge_metabolite_identity_conflict"
            record["rejection_reason"] = "postmerge_metabolite_identity_conflict"
            record["identity_conflicts"] = conflicts
            record.update(fallback)
            return
        _apply_conflict_resolution(reaction, operations, record)
        conflicts = postmerge_metabolite_identity_conflicts(reaction, model)
        if conflicts:
            raise RuntimeError(
                f"Conflict resolution left conflicts on {reaction.id}"
            )
    if conflicts:
        record["status"] = "rejected_postmerge_metabolite_identity_conflict"
        record["rejection_reason"] = "postmerge_metabolite_identity_conflict"
        record["identity_conflicts"] = conflicts
        return
    model.add_reactions([reaction])
    try:
        final_imbalance = model.reactions.get_by_id(reaction.id).check_mass_balance()
    except (TypeError, ValueError) as error:
        raise RuntimeError(
            f"Added reaction {reaction.id} has unverifiable final balance"
        ) from error
    if not reaction_is_mass_and_charge_balanced(final_imbalance):
        raise RuntimeError(
            f"Added reaction {reaction.id} became imbalanced in the model: "
            f"{final_imbalance}"
        )
    record["status"] = "added"


def _gpr_genes_and_complexity(rule):
    """Gene set of a GPR rule; flags complex (and-containing) rules.

    cobra normalizes a malformed rule to an empty GPR (body None) with a
    warning instead of raising, so the observable outcome is an empty gene
    set; the caller's try/except only guards versions that do raise.
    """
    if not isinstance(rule, str) or not rule.strip():
        return set(), False
    gpr = GPR.from_string(rule)
    if gpr.body is None:
        return set(), False
    has_complex = any(
        isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And)
        for node in ast.walk(gpr.body)
    )
    return set(gpr.genes), has_complex


def _verify_existing_matches_source(existing, source):
    """Machine-check that a model reaction and its candidate source agree."""
    existing_metabolites = {m.id: m for m in existing.metabolites}
    source_metabolites = {m.id: m for m in source.metabolites}
    if set(existing_metabolites) != set(source_metabolites):
        return "stoichiometry_or_identity_differ", {
            "model_only_metabolites": sorted(
                set(existing_metabolites) - set(source_metabolites)
            ),
            "source_only_metabolites": sorted(
                set(source_metabolites) - set(existing_metabolites)
            ),
        }
    existing_coefficients = {
        metabolite_id: float(existing.metabolites[metabolite])
        for metabolite_id, metabolite in existing_metabolites.items()
    }
    source_coefficients = {
        metabolite_id: float(source.metabolites[metabolite])
        for metabolite_id, metabolite in source_metabolites.items()
    }
    coefficient_differences = sorted(
        metabolite_id
        for metabolite_id, coefficient in existing_coefficients.items()
        if not math.isclose(
            coefficient,
            source_coefficients[metabolite_id],
            rel_tol=0.0,
            abs_tol=1e-9,
        )
    )
    if coefficient_differences:
        return "stoichiometry_or_identity_differ", {
            "coefficient_differences": coefficient_differences,
        }
    incomplete = []
    identity_differences = []
    for metabolite_id in sorted(existing_metabolites):
        existing_identity = _metabolite_identity(existing_metabolites[metabolite_id])
        source_identity = _metabolite_identity(source_metabolites[metabolite_id])
        if _identity_incomplete(existing_identity) or _identity_incomplete(
            source_identity
        ):
            incomplete.append(metabolite_id)
        elif existing_identity != source_identity:
            identity_differences.append({
                "metabolite_id": metabolite_id,
                "model_identity": existing_identity,
                "source_identity": source_identity,
            })
    if incomplete:
        return "unverifiable_identity", {"incomplete_metabolites": incomplete}
    if identity_differences:
        return "stoichiometry_or_identity_differ", {
            "identity_differences": identity_differences,
        }
    if not (
        math.isclose(
            existing.lower_bound, source.lower_bound, rel_tol=0.0, abs_tol=1e-9
        )
        and math.isclose(
            existing.upper_bound, source.upper_bound, rel_tol=0.0, abs_tol=1e-9
        )
    ):
        return "direction_or_bounds_differ", {
            "model_bounds": [existing.lower_bound, existing.upper_bound],
            "source_bounds": [source.lower_bound, source.upper_bound],
        }
    return None, None


def _nadp_family(name):
    """Heuristic cofactor family from a metabolite name ('nadp' vs 'nad')."""
    lowered = (name or "").lower()
    if "nadp" in lowered:
        return "nadp"
    if "nad" in lowered:
        return "nad"
    return None


def _has_cofactor_name_pair(names_left, names_right):
    """Heuristic NAD(H)<->NADP(H) swap detection by metabolite names.

    This is a name-based heuristic only: a differing metabolite pair counts
    when one side's name contains 'nadp' and the other side's contains 'nad'
    without 'p'. It does not prove cofactor chemistry.
    """
    families_left = {_nadp_family(name) for name in names_left}
    families_right = {_nadp_family(name) for name in names_right}
    families_left.discard(None)
    families_right.discard(None)
    return bool(families_left) and bool(families_right) and (
        families_left != families_right
    )


def _classify_reaction_divergence(left, right):
    """Extended divergence classification between two reaction versions.

    Returns None when the two versions agree under the same criteria as
    _verify_existing_matches_source, else {"classification", "detail"}.
    The classification extends the two-way conflict vocabulary with
    reversed_stoichiometry (identical metabolite sets, all coefficients
    sign-flipped) and the heuristic cofactor_pair label.
    """
    left_metabolites = {m.id: m for m in left.metabolites}
    right_metabolites = {m.id: m for m in right.metabolites}
    if set(left_metabolites) != set(right_metabolites):
        only_left = sorted(set(left_metabolites) - set(right_metabolites))
        only_right = sorted(set(right_metabolites) - set(left_metabolites))
        detail = {
            "left_only_metabolites": only_left,
            "right_only_metabolites": only_right,
        }
        if _has_cofactor_name_pair(
            [left_metabolites[m].name for m in only_left],
            [right_metabolites[m].name for m in only_right],
        ):
            return {"classification": "cofactor_pair", "detail": detail}
        return {"classification": "stoichiometry_or_identity_differ", "detail": detail}
    coefficient_differences = sorted(
        metabolite_id
        for metabolite_id, metabolite in left_metabolites.items()
        if not math.isclose(
            float(left.metabolites[metabolite]),
            float(right.metabolites[right_metabolites[metabolite_id]]),
            rel_tol=0.0,
            abs_tol=1e-9,
        )
    )
    if coefficient_differences:
        reversed_all = all(
            math.isclose(
                float(left.metabolites[left_metabolites[metabolite_id]]),
                -float(right.metabolites[right_metabolites[metabolite_id]]),
                rel_tol=0.0,
                abs_tol=1e-9,
            )
            for metabolite_id in left_metabolites
        ) and any(
            not math.isclose(
                float(left.metabolites[left_metabolites[metabolite_id]]),
                0.0,
                rel_tol=0.0,
                abs_tol=1e-9,
            )
            for metabolite_id in left_metabolites
        )
        detail = {"coefficient_differences": coefficient_differences}
        if reversed_all:
            return {"classification": "reversed_stoichiometry", "detail": detail}
        return {"classification": "stoichiometry_or_identity_differ", "detail": detail}
    incomplete = []
    identity_differences = []
    for metabolite_id in sorted(left_metabolites):
        left_identity = _metabolite_identity(left_metabolites[metabolite_id])
        right_identity = _metabolite_identity(right_metabolites[metabolite_id])
        if _identity_incomplete(left_identity) or _identity_incomplete(
            right_identity
        ):
            incomplete.append(metabolite_id)
        elif left_identity != right_identity:
            identity_differences.append({
                "metabolite_id": metabolite_id,
                "left_identity": left_identity,
                "right_identity": right_identity,
            })
    if incomplete:
        return {
            "classification": "unverifiable_identity",
            "detail": {"incomplete_metabolites": incomplete},
        }
    if identity_differences:
        if _has_cofactor_name_pair(
            [left_metabolites[row["metabolite_id"]].name for row in identity_differences],
            [right_metabolites[row["metabolite_id"]].name for row in identity_differences],
        ):
            return {
                "classification": "cofactor_pair",
                "detail": {"identity_differences": identity_differences},
            }
        return {
            "classification": "stoichiometry_or_identity_differ",
            "detail": {"identity_differences": identity_differences},
        }
    if not (
        math.isclose(
            left.lower_bound, right.lower_bound, rel_tol=0.0, abs_tol=1e-9
        )
        and math.isclose(
            left.upper_bound, right.upper_bound, rel_tol=0.0, abs_tol=1e-9
        )
    ):
        return {
            "classification": "direction_or_bounds_differ",
            "detail": {
                "left_bounds": [left.lower_bound, left.upper_bound],
                "right_bounds": [right.lower_bound, right.upper_bound],
            },
        }
    return None


def _threeway_gate(existing, template_reaction, pool_reaction):
    """Three-way identity gate; returns a divergence detail dict or None.

    Compared pairwise in the fixed order template_pool, model_template,
    model_pool; the first divergent pair sets the headline classification.
    """
    pairs = (
        ("template_pool", template_reaction, pool_reaction),
        ("model_template", existing, template_reaction),
        ("model_pool", existing, pool_reaction),
    )
    comparisons = {}
    divergent_pairs = []
    for label, left, right in pairs:
        divergence = _classify_reaction_divergence(left, right)
        comparisons[label] = divergence
        if divergence is not None:
            divergent_pairs.append(label)
    if not divergent_pairs:
        return None
    return {
        "classification": comparisons[divergent_pairs[0]]["classification"],
        "divergent_pairs": divergent_pairs,
        "pairs": comparisons,
    }


def _rhea_key(value):
    """Normalize a spreadsheet rhea cell to the candidate ID's number string."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)) or int(value) != value:
            return None
        return str(int(value))
    text = str(value).strip()
    return text or None


def _multisource_gene_votes(name):
    """Per-(gene, rhea) counts of distinct non-empty annotation sources.

    The multiple_score table is keyed per (gene, rhea) pair, so only evidence
    rows pointing at the candidate's own rhea number count. With the current
    profile (ensemble_score 0.4 while the largest single-tool weight in
    calculate_score is anno_B's 0.3), every candidate that reached the
    gpr_score table already has >= 2 non-empty sources at the pair level, so
    this filter is a redundant safeguard there, not a stricter screen; a
    missing evidence file still blocks every merge.
    Returns None when the evidence file is absent.
    """
    path = Path("working") / name / f"{name}multiple_score.xlsx"
    if not path.is_file():
        return None
    table = pd.read_excel(path)
    if not {"reaction", "rhea"} <= set(table.columns):
        return {}
    present_columns = [
        column
        for column in MULTISOURCE_ANNOTATION_COLUMNS
        if column in table.columns
    ]
    sources = {}
    for _, row in table.iterrows():
        gene = row["reaction"]
        rhea = _rhea_key(row["rhea"])
        if not isinstance(gene, str) or not gene.strip() or rhea is None:
            continue
        pair_sources = sources.setdefault((gene, rhea), set())
        for column in present_columns:
            value = row[column]
            if pd.notna(value) and str(value).strip():
                pair_sources.add(column)
    return {pair: len(columns) for pair, columns in sources.items()}


def _merge_already_present(
    existing, row, record, refmodel, refids, reactionpools,
    already_present_policy, name, rhea, model, pool_raw_names,
):
    """Verify an already-present candidate and merge its GPR when identical."""
    reaction_id = record["reaction_id"]
    in_refmodel = reaction_id in refids
    try:
        pool_version = reactionpools.reactions.get_by_id(reaction_id)
    except KeyError:
        pool_version = None
    if not in_refmodel and pool_version is None:
        record["status"] = "already_present_conflict"
        record["conflict_kind"] = "unverifiable_identity"
        record["conflict_detail"] = {"reason": "candidate_source_absent"}
        return
    if in_refmodel:
        source_reaction = refmodel.reactions.get_by_id(reaction_id)
        source_label = "yeast9_reference"
    else:
        source_reaction = pool_version
        source_label = "frozen_reaction_pool"
    record["candidate_source"] = source_label
    threeway_check = None
    if already_present_policy == "merge_threeway_or_multisource":
        if in_refmodel and pool_version is not None:
            # Mechanistically conservative gate: the existing model reaction,
            # the Yeast9 reference version, and the frozen-pool (Rhea-derived)
            # version must all be identical before any GPR is merged, so a
            # direction or cofactor conflict between the two templates can
            # never be silently merged away.
            gate = _threeway_gate(
                existing,
                refmodel.reactions.get_by_id(reaction_id),
                pool_version,
            )
            if gate is not None:
                record["status"] = "already_present_conflict"
                record["conflict_kind"] = "template_pool_divergence"
                record["conflict_detail"] = gate
                return
            threeway_check = "model_template_pool_identical"
        elif in_refmodel:
            # No third version to check against: fall back to the two-way
            # verification against the reference model.
            record["pool_coverage"] = "not_in_pool"
    conflict_kind, detail = _verify_existing_matches_source(existing, source_reaction)
    candidate_genes = set()
    existing_genes = set()
    if conflict_kind is None:
        candidate_value = row.get("reaction")
        candidate_rule = candidate_value if isinstance(candidate_value, str) else ""
        try:
            existing_genes, existing_complex = _gpr_genes_and_complexity(
                existing.gene_reaction_rule
            )
            candidate_genes, candidate_complex = _gpr_genes_and_complexity(
                candidate_rule
            )
        except (SyntaxError, TypeError, ValueError):
            conflict_kind = "unverifiable_identity"
            detail = {"reason": "gpr_parse_failure"}
        else:
            if existing_complex or candidate_complex:
                conflict_kind = "complex_gpr_not_merged"
                detail = {
                    "existing_rule": existing.gene_reaction_rule,
                    "candidate_rule": candidate_rule,
                }
    if conflict_kind is not None:
        record["status"] = "already_present_conflict"
        record["conflict_kind"] = conflict_kind
        record["conflict_detail"] = detail
        return
    # Metabolite names are informational only: formula/charge/compartment
    # already pin the passing identity, and real-world name spellings differ.
    # A name mismatch is recorded but never blocks a merge; formula-equal or
    # name-equal isomers remain a residual risk either way. Comparison uses
    # the author's metdict equivalence key (lowercase, parentheses stripped),
    # the same relation the add path uses when it swaps metabolites by name;
    # raw spelling differences like 'H(+)' vs 'H+' are not mismatches.
    # The model side reads the REGISTRY name (what the validator recomputes
    # from the saved model). The pool side reads the RAW frozen JSON name:
    # cobra add_reactions re-points a same-stage added reaction's metabolite
    # references to the registry objects, so the in-memory pool reaction then
    # carries registry names and would hide real differences (server case:
    # registry 'ammonium' vs frozen pool 'NH4(+)' on rhea_45048_c).
    source_metabolites = {m.id: m for m in source_reaction.metabolites}
    name_mismatch_metabolites = sorted(
        metabolite.id
        for metabolite in existing.metabolites
        if model.metabolites.get_by_id(metabolite.id).name
        and _source_name(source_label, metabolite.id, source_metabolites,
                         pool_raw_names)
        and _author_name_key(model.metabolites.get_by_id(metabolite.id).name)
        != _author_name_key(
            _source_name(source_label, metabolite.id, source_metabolites,
                         pool_raw_names)
        )
    )
    verification = {
        "candidate_source": source_label,
        "stoichiometry": "metabolite_ids_and_signed_coefficients_equal",
        "metabolite_identity": "formula_charge_compartment_complete_and_equal",
        "flux_bounds": "lower_and_upper_bounds_equal",
        "metabolite_names": "informational_only_not_a_pass_condition",
        "name_mismatch": bool(name_mismatch_metabolites),
        "name_mismatch_metabolites": name_mismatch_metabolites,
        "shared_metabolite_count": len(existing.metabolites),
        "existing_gene_count": len(existing_genes),
        "candidate_gene_count": len(candidate_genes),
        "residual_risks": (
            "same-formula isomers and same-name isomers cannot be excluded "
            "by this check"
        ),
    }
    if threeway_check is not None:
        verification["three_way_check"] = threeway_check
    if already_present_policy in {
        "merge_verified_identical_or_multisource",
        "merge_threeway_or_multisource",
    }:
        votes = _multisource_gene_votes(name)
        eligible = sorted(
            gene
            for gene in candidate_genes
            if votes is not None
            and votes.get((gene, rhea), 0) >= MULTISOURCE_MIN_SOURCES
        )
        record["skipped_insufficient_evidence"] = sorted(
            set(candidate_genes) - set(eligible)
        )
        verification["min_annotation_sources"] = MULTISOURCE_MIN_SOURCES
        verification["evidence_granularity"] = "gene_rhea_pair"
        if not eligible:
            record["status"] = "already_present_no_merge"
            record["no_merge_reason"] = "all_candidate_genes_below_min_sources"
            record["verification"] = verification
            return
        candidate_genes = set(eligible)
    merged_genes = sorted(candidate_genes - existing_genes)
    existing.gene_reaction_rule = " or ".join(sorted(existing_genes | candidate_genes))
    record["status"] = "already_present_gpr_merged"
    record["merge_kind"] = (
        "into_empty_gpr" if not existing_genes else "into_existing_gpr"
    )
    record["merged_genes"] = merged_genes
    record["verification"] = verification


def _mark_blocked_by_construction(candidate_audit, saved_model_path):
    """Flag resolved additions whose new metabolites appear in one reaction.

    Reloads the saved stage model so the flag reflects the persisted state.
    A metabolite is "new" for a candidate when this stage created it (a
    non-reused __pool suffix) or when it entered the model with this very
    reaction (recorded as unmatched). A resolved addition with no new
    metabolites is never marked.
    """
    resolved = [
        record
        for record in candidate_audit
        if record.get("addition_kind")
        in {"suffix_metabolite", "remapped_metabolite"}
    ]
    if not resolved:
        return
    saved_model = cobra.io.load_yaml_model(saved_model_path)
    use_count = {}
    for reaction in saved_model.reactions:
        for metabolite in reaction.metabolites:
            use_count[metabolite.id] = use_count.get(metabolite.id, 0) + 1
    for record in resolved:
        new_ids = [
            entry["new_metabolite_id"]
            for entry in record.get("suffix_metabolites", [])
            if not entry.get("reused_existing_pool_metabolite")
        ]
        new_ids += record.get("unmatched_metabolites", [])
        record["blocked_by_construction"] = bool(new_ids) and all(
            use_count.get(metabolite_id, 0) == 1 for metabolite_id in new_ids
        )


def _build_reaction_pool_coverage(
    candidate_audit, reactionpool, already_present_policy, conflict_policy
):
    """Defect-type coverage and escalation signals for one stage run."""
    input_count = len(candidate_audit)
    status_counts = {}
    for record in candidate_audit:
        status = record["status"]
        status_counts[status] = status_counts.get(status, 0) + 1

    def ratio(count):
        return (count / input_count) if input_count else 0.0

    suffix_added = sum(
        1
        for record in candidate_audit
        if record.get("addition_kind") == "suffix_metabolite"
    )
    remapped_added = sum(
        1
        for record in candidate_audit
        if record.get("addition_kind") == "remapped_metabolite"
    )
    added_blocked = sum(
        1 for record in candidate_audit if record.get("blocked_by_construction")
    )
    identity_rejected = status_counts.get(
        "rejected_postmerge_metabolite_identity_conflict", 0
    )
    already_present_breakdown = {
        "author_skip": status_counts.get("already_present", 0),
        "gpr_merged": status_counts.get("already_present_gpr_merged", 0),
        "no_merge": status_counts.get("already_present_no_merge", 0),
        "conflict": status_counts.get("already_present_conflict", 0),
    }
    conflict_kinds = {}
    for record in candidate_audit:
        if record["status"] == "already_present_conflict":
            kind = record.get("conflict_kind")
            conflict_kinds[kind] = conflict_kinds.get(kind, 0) + 1
    defects = {
        "missing_from_reaction_pool": {
            "count": status_counts.get("missing_from_reaction_pool", 0),
        },
        "rejected_mass_or_charge_balance": {
            "count": status_counts.get("rejected_mass_or_charge_balance", 0),
        },
        "identity_conflict": {
            "count": identity_rejected + suffix_added + remapped_added,
            "breakdown": {
                "rejected": identity_rejected,
                "suffix_metabolite_added": suffix_added,
                "remapped_metabolite_added": remapped_added,
                "added_blocked_by_construction": added_blocked,
            },
        },
        "already_present": {
            "count": sum(already_present_breakdown.values()),
            "breakdown": already_present_breakdown,
        },
    }
    for defect in defects.values():
        defect["ratio"] = ratio(defect["count"])
    escalation_signals = [
        {
            "type": defect_type,
            "count": defect["count"],
            "denominator": input_count,
            "ratio": defect["ratio"],
        }
        for defect_type, defect in defects.items()
        if defect["ratio"] > ESCALATION_RATIO_THRESHOLD
    ]
    return {
        "input_candidate_count": input_count,
        "policies": {
            "already_present_policy": already_present_policy,
            "conflict_policy": conflict_policy,
        },
        "status_counts": status_counts,
        "status_ratios": {
            status: ratio(count) for status, count in sorted(status_counts.items())
        },
        "defects": defects,
        "already_present_conflict_kinds": conflict_kinds,
        "escalation_signals": escalation_signals,
        "missing_from_reaction_pool_ids": [
            record["reaction_id"]
            for record in candidate_audit
            if record["status"] == "missing_from_reaction_pool"
        ],
        # Plain pool size statistic. The pool's vintage must never be inferred
        # from its maximum reaction ID; IDs are not dated version markers.
        "reaction_pool": {"reaction_count": len(reactionpool.reactions)},
    }


def model_reaction(name,refmodel,strict_identity=False,already_present_policy="author_skip",conflict_policy="reject"):
    if already_present_policy not in ALREADY_PRESENT_POLICIES:
        raise ValueError(
            f"Unknown already-present policy: {already_present_policy}"
        )
    if conflict_policy not in CONFLICT_POLICIES:
        raise ValueError(f"Unknown conflict policy: {conflict_policy}")
    if not strict_identity and (
        already_present_policy != "author_skip" or conflict_policy != "reject"
    ):
        raise ValueError(
            "Non-default reaction-pool policies require strict_identity=True"
        )
    refmodel=cobra.io.read_sbml_model(f'models/{refmodel}')
    refids=[r.id for r in refmodel.reactions]
    model = cobra.io.load_yaml_model(f'./working/{name}/model2{name}.yml')
    fba_begin=model.optimize().objective_value
    gpr_e = pd.read_excel(f'./working/{name}/{name}gpr_score.xlsx')
    pklopen=open('data_available/general2specific.pkl', 'rb')
    general2spe=pickle.load(pklopen)
    reactionpool=cobra.io.load_json_model('models/reaction_pool_final.json')
    with open('models/reaction_pool_final.json', encoding='utf-8') as pool_stream:
        pool_raw = json.load(pool_stream)
    pool_raw_names = {
        met["id"]: met.get("name")
        for met in pool_raw.get("metabolites", [])
        if isinstance(met, dict) and isinstance(met.get("id"), str)
    }
    pklopen.close()
    Reactions=[]
    reactionpools=reactionpool.copy()
    metdict={metss.name.lower():metss for metss in model.metabolites if metss.compartment=='c'}
    for metss in model.metabolites:
        if metss.compartment=='c':
            metdict[metss.name.lower().replace('(','').replace(')','')]=metss
    # Stage-start metabolite snapshot: remap matching must not depend on how
    # many candidates already ran (keeps draft-model revalidation exact).
    stage_start_metabolite_ids = {metss.id for metss in model.metabolites}
    candidate_audit = []
    for index, row in tqdm(gpr_e.iterrows()):
        reactions=[]
        reactions.append(str(row['rhea']))
        for rea in reactions:
            reaction_id = 'rhea_'+str(rea)+'_c'
            record = {
                "input_row": int(index) if isinstance(index, int) else str(index),
                "reaction_id": reaction_id,
                "status": None,
                "unmatched_metabolites": [],
            }
            try:
                existing_reaction = model.reactions.get_by_id(reaction_id)#已经添加过了
            except KeyError:
                existing_reaction = None
            if existing_reaction is not None:
                if already_present_policy == "author_skip":
                    record["status"] = "already_present"
                    candidate_audit.append(record)
                    continue
                _merge_already_present(
                    existing_reaction, row, record, refmodel, refids,
                    reactionpools, already_present_policy, name,
                    _rhea_key(row["rhea"]), model, pool_raw_names,
                )
                candidate_audit.append(record)
                continue
            already =0
            try:
               reaction_add=reactionpools.reactions.get_by_id(reaction_id)
            except KeyError:
                print(reaction_id+' does not exist in reaction pool')
                record["status"] = "missing_from_reaction_pool"
                candidate_audit.append(record)
                continue
            if (reaction_add.reactants==[] or reaction_add.products==[]) and reaction_add.id not in refids:
                print('reaction error')
                record["status"] = "incomplete_stoichiometry"
                candidate_audit.append(record)
                continue
            if reaction_add.id in refids:
                reaction_add=refmodel.reactions.get_by_id(reaction_add.id)
                already=1
            if strict_identity:
                record["candidate_source"] = (
                    "yeast9_reference" if already else "frozen_reaction_pool"
                )
            reaction_add.gene_reaction_rule=row['reaction']
            reaction_add.annotation['from']='NonHomo'
            #TODO: DEBUG REACTION ERROR OCCUR
            a=list(reaction_add.metabolites.keys())
            for met in a:
                 if already==1:
                     continue
                 try:
                     model.metabolites.get_by_id(met)
                     continue
                 except KeyError:
                     pass
                 metsearch=met.name.lower() if met.name.lower() in list(metdict.keys()) else met.name.lower().replace('(','').replace(')','')
                 try:
                     if metdict[metsearch].elements==met.elements and metdict[
                         metsearch].charge==met.charge and metdict[metsearch].compartment==met.compartment:
                         stoich=reaction_add.metabolites[met]
                         reaction_add.add_metabolites({metdict[metsearch]:stoich})
                         reaction_add.add_metabolites({met:-1*stoich})
                         continue
                     # try:
                     #   if (metdict[metsearch].elements['H'] - met.elements['H'] == metdict[
                     #     metsearch].charge - met.charge and metdict[metsearch].elements.pop(
                     #     'H') == met.elements.pop('H')) and metdict[metsearch].compartment==met.compartment:
                     #     stoich = reaction_add.metabolites[met]
                     #     reaction_add.add_metabolites({metdict[metsearch]: reaction_add.metabolites[stoich],
                     #                                   met: -1 * reaction_add.metabolites[stoich]})
                     # except:
                     #   pass
                 except KeyError:
                     record["unmatched_metabolites"].append(met.id)
            if strict_identity:
                admit_candidate_reaction(
                    reaction_add, model, record, conflict_policy=conflict_policy,
                    stage_start_metabolite_ids=stage_start_metabolite_ids,
                )
            elif check_reactions(reaction_add, model) != 0:
                model.add_reactions([reaction_add])
                record["status"] = "added"
            else:
                record["status"] = "rejected_mass_or_charge_balance"
            candidate_audit.append(record)
    model.id=f'{name}-GEM'
    for i in model.reactions:
        keys_to_remove = [k for k, v in i.annotation.items() if v in [[], '']]
        # 删除这些键
        for k in keys_to_remove:
            del i.annotation[k]
    cobra.io.save_yaml_model(model, f'./working/{name}/{name}-GEM_withgaps.yml')
    _mark_blocked_by_construction(
        candidate_audit, f'./working/{name}/{name}-GEM_withgaps.yml'
    )
    audit_path = Path("working") / name / "reaction_pool_audit.json"
    status_counts = {}
    for record in candidate_audit:
        status = record["status"]
        status_counts[status] = status_counts.get(status, 0) + 1
    audit_path.write_text(
        json.dumps(
            {
                "input_candidate_count": len(candidate_audit),
                "status_counts": status_counts,
                "candidates": candidate_audit,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    coverage_path = Path("working") / name / "reaction_pool_coverage.json"
    coverage_path.write_text(
        json.dumps(
            _build_reaction_pool_coverage(
                candidate_audit,
                reactionpool,
                already_present_policy,
                conflict_policy,
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
