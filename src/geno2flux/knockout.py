"""Knockout prediction: single-reaction, single-gene, and paralog-family KO."""

from __future__ import annotations

import csv
import re
from math import isfinite
from pathlib import Path
from typing import Any

import cobra


def _loss(obj: float | None, base: float, status: str) -> float:
    """Compute loss = 1 - v_KO/v_base (pre-registered formula, no clamping).

    Non-optimal status, missing or ~zero objective → loss=1.0 (target dead).
    Negative values are kept verbatim: a KO that improves production is real.
    """
    if not isfinite(base) or base <= 1e-12:
        raise ValueError("KO loss requires a finite positive baseline objective")
    if status != "optimal":
        return 1.0
    if obj is not None and not isfinite(obj):
        raise ValueError("KO loss requires a finite knockout objective")
    if obj is None or obj <= 1e-12:
        return 1.0
    return 1.0 - obj / base


def _target_baseline(model: cobra.Model, target_rxn_id: str) -> float:
    """Set the target objective and reject an undefined loss denominator.

    Call inside a model context so both success and failure restore the
    caller's objective and direction.
    """
    model.objective = target_rxn_id
    model.objective_direction = "max"
    sol = model.optimize()
    base = sol.objective_value
    if sol.status != "optimal" or base is None or not isfinite(base) or base <= 1e-12:
        raise ValueError(
            f"KO target {target_rxn_id!r} requires an optimal finite positive "
            f"baseline (status={sol.status!r}, objective={base!r})"
        )
    return float(base)


def run_single_reaction_ko(
    model: cobra.Model,
    target_rxn_id: str,
    out_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Rank reaction knockouts against a positive target baseline; restore model."""
    rows = []
    with model:
        W = _target_baseline(model, target_rxn_id)
        for rxn in model.reactions:
            rid = rxn.id
            if rid == target_rxn_id:
                continue
            with model:
                rxn.knock_out()
                sol = model.optimize()
                obj = sol.objective_value if sol.status == "optimal" else None
                loss = _loss(obj, W, sol.status)
                rows.append({
                    "reaction_id": rid,
                    "reaction_name": rxn.name,
                    "subsystem": rxn.subsystem or "",
                    "obj": obj,
                    "loss": loss,
                    "status": sol.status,
                })
    # Only rank optimal solutions; non-optimal go to separate list.
    optimal_rows = [r for r in rows if r["status"] == "optimal"]
    nonoptimal_rows = [r for r in rows if r["status"] != "optimal"]
    optimal_rows.sort(key=lambda r: (-r["loss"], r["reaction_id"]))
    for i, r in enumerate(optimal_rows):
        r["ko_rank"] = i + 1
    nonoptimal_rows.sort(key=lambda r: r["reaction_id"])
    return optimal_rows, nonoptimal_rows


def run_single_gene_ko(
    model: cobra.Model,
    target_rxn_id: str,
    out_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Rank gene knockouts against a positive target baseline; restore model."""
    rows = []
    with model:
        W = _target_baseline(model, target_rxn_id)
        for g in model.genes:
            gid = g.id
            with model:
                g.knock_out()
                sol = model.optimize()
                obj = sol.objective_value if sol.status == "optimal" else None
                loss = _loss(obj, W, sol.status)
                rows.append({
                    "gene_id": gid,
                    "gene_name": g.name,
                    "obj": obj,
                    "loss": loss,
                    "status": sol.status,
                })
    optimal_rows = [r for r in rows if r["status"] == "optimal"]
    nonoptimal_rows = [r for r in rows if r["status"] != "optimal"]
    optimal_rows.sort(key=lambda r: (-r["loss"], r["gene_id"]))
    for i, r in enumerate(optimal_rows):
        r["ko_rank"] = i + 1
    nonoptimal_rows.sort(key=lambda r: r["gene_id"])
    return optimal_rows, nonoptimal_rows


def run_family_ko(
    model: cobra.Model,
    target_rxn_id: str,
    families: dict[str, list[str]],
    out_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Knock out each family's model genes together and restore the model."""
    rows = []
    with model:
        W = _target_baseline(model, target_rxn_id)
        model_gene_ids = {g.id for g in model.genes}
        for fam_id, members in sorted(families.items()):
            with model:
                for gid in members:
                    if gid in model_gene_ids:
                        model.genes.get_by_id(gid).knock_out()
                sol = model.optimize()
                obj = sol.objective_value if sol.status == "optimal" else None
                loss = _loss(obj, W, sol.status)
                rows.append({
                    "family_id": fam_id,
                    "members": ";".join(sorted(members)),
                    "size": len(members),
                    "obj": obj,
                    "loss": loss,
                    "status": sol.status,
                })
    optimal_rows = [r for r in rows if r["status"] == "optimal"]
    nonoptimal_rows = [r for r in rows if r["status"] != "optimal"]
    optimal_rows.sort(key=lambda r: (-r["loss"], r["family_id"]))
    for i, r in enumerate(optimal_rows):
        r["ko_rank"] = i + 1
    nonoptimal_rows.sort(key=lambda r: r["family_id"])
    return optimal_rows, nonoptimal_rows


def is_pathway_reaction(rxn: cobra.Reaction) -> bool:
    """Filter: boundary OR metabolites span ≥2 compartments OR respiratory subsystem."""
    if rxn in rxn.model.boundary:
        return True
    comps = {m.compartment or "" for m in rxn.metabolites}
    if len(comps) >= 2:
        return True
    sub = (rxn.subsystem or "").lower()
    if re.search(
        r"oxidative phosphorylation|respiratory chain|electron transport|atp synthase",
        sub,
    ):
        return True
    return False


def write_ko_tsv(path: Path, rows: list[dict[str, Any]], id_col: str,
                 nonoptimal: bool = False, filtered: bool = False) -> None:
    """Write KO results to TSV.
    nonoptimal=True writes the nonoptimal supplement (no rank column).
    filtered=True uses filtered_rank as primary rank for pathway-filtered tables.
    """
    if filtered:
        fieldnames = {
            "reaction_id": ["filtered_rank", "ko_rank", "reaction_id", "reaction_name", "subsystem", "obj", "loss", "status"],
            "gene_id": ["filtered_rank", "ko_rank", "gene_id", "gene_name", "obj", "loss", "status"],
            "family_id": ["filtered_rank", "ko_rank", "family_id", "members", "size", "obj", "loss", "status"],
        }[id_col]
    else:
        fieldnames = {
            "reaction_id": ["ko_rank", "reaction_id", "reaction_name", "subsystem", "obj", "loss", "status"],
            "gene_id": ["ko_rank", "gene_id", "gene_name", "obj", "loss", "status"],
            "family_id": ["ko_rank", "family_id", "members", "size", "obj", "loss", "status"],
        }[id_col]
    if nonoptimal:
        fieldnames = [f for f in fieldnames if f not in ("ko_rank", "filtered_rank")]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run_ko_suite(
    model: cobra.Model,
    target_rxn_id: str,
    families: dict[str, list[str]] | None = None,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    """Run all KO types, preserving the caller's model and writing all tables."""
    with model:
        base_objective = _target_baseline(model, target_rxn_id)
    out_dir = out_dir or Path("ko")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Single reaction KO: returns (optimal, nonoptimal)
    rxn_opt, rxn_nonopt = run_single_reaction_ko(model, target_rxn_id, out_dir)
    write_ko_tsv(out_dir / "single_reaction_ko.tsv", rxn_opt, "reaction_id")
    write_ko_tsv(out_dir / "single_reaction_ko_nonoptimal.tsv", rxn_nonopt, "reaction_id", nonoptimal=True)

    # Single gene KO
    gene_opt, gene_nonopt = run_single_gene_ko(model, target_rxn_id, out_dir)
    write_ko_tsv(out_dir / "single_gene_ko.tsv", gene_opt, "gene_id")
    write_ko_tsv(out_dir / "single_gene_ko_nonoptimal.tsv", gene_nonopt, "gene_id", nonoptimal=True)

    # Filtered reaction KO (post-hoc analysis) — only optimal rows
    filtered = [r for r in rxn_opt if not is_pathway_reaction(model.reactions.get_by_id(r["reaction_id"]))]
    filtered.sort(key=lambda r: (-r["loss"], r["reaction_id"]))
    for i, r in enumerate(filtered):
        r["filtered_rank"] = i + 1
    write_ko_tsv(out_dir / "single_reaction_ko_filtered.tsv", filtered, "reaction_id", filtered=True)

    # Family KO
    fam_opt, fam_nonopt = [], []
    if families:
        fam_opt, fam_nonopt = run_family_ko(model, target_rxn_id, families, out_dir)
    write_ko_tsv(out_dir / "family_ko.tsv", fam_opt, "family_id")
    write_ko_tsv(out_dir / "family_ko_nonoptimal.tsv", fam_nonopt, "family_id", nonoptimal=True)

    return {
        "target": target_rxn_id,
        "base_objective": base_objective,
        "single_reaction_ko": len(rxn_opt),
        "single_reaction_ko_nonoptimal": len(rxn_nonopt),
        "single_gene_ko": len(gene_opt),
        "single_gene_ko_nonoptimal": len(gene_nonopt),
        "family_ko": len(fam_opt),
        "family_ko_nonoptimal": len(fam_nonopt),
        "top_reaction_loss": rxn_opt[0]["loss"] if rxn_opt else None,
        "top_gene_loss": gene_opt[0]["loss"] if gene_opt else None,
        "top_family_loss": fam_opt[0]["loss"] if fam_opt else None,
    }
