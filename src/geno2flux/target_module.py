"""Configurable target reaction & product-module mounting for metabolic screening.

The caller must choose one of two target modes; no biological product
or strain is assumed by the toolkit:

1. ``existing_reaction``: use a reaction already present in the model (for
   example an exchange, sink, or biosynthetic reaction) as the screening target.
2. ``custom`` (mapping or JSON spec file): declare expected host metabolites,
   new bookkeeping metabolites, internal reactions (enforcing mass and charge
   balance by default), and a boundary/demand reaction to optimize.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Mapping

import cobra





def load_target_spec(
    spec_or_path: Mapping[str, Any] | Path | str | None,
    *,
    target_reaction: str | None = None,
    base_dir: Path | None = None,
) -> dict[str, Any]:
    """Normalize a target specification from a mapping or JSON file path."""
    if isinstance(spec_or_path, (str, Path)):
        candidate_str = str(spec_or_path).strip()
        path = Path(candidate_str)
        if not path.is_absolute() and base_dir is not None:
            path = (base_dir / path).resolve()
        if not path.is_file():
            raise ValueError(f"Target module specification file not found: {path}")
        loaded = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(loaded, dict):
            raise ValueError(f"Target module JSON must be an object: {path}")
        return load_target_spec(loaded, target_reaction=target_reaction, base_dir=path.parent)

    if isinstance(spec_or_path, Mapping):
        spec = dict(spec_or_path)
        mode = spec.get("mode")
        if target_reaction and "target_reaction" not in spec:
            spec["target_reaction"] = target_reaction
        return spec

    if target_reaction:
        return {"mode": "existing_reaction", "target_reaction": target_reaction}

    raise ValueError(
        "A target is required: set 'target_reaction' for an existing model reaction "
        "or provide a custom 'target_module' specification"
    )


def mount_target_module(
    model: cobra.Model,
    spec_or_path: Mapping[str, Any] | Path | str | None = None,
    *,
    target_reaction: str | None = None,
    base_dir: Path | None = None,
) -> dict[str, Any]:
    """Validate or mount the configured target module in-place and return a receipt."""
    spec = load_target_spec(
        spec_or_path, target_reaction=target_reaction, base_dir=base_dir
    )
    mode = spec.get("mode")
    if mode is None:
        mode = "custom" if ("reactions" in spec or "new_metabolites" in spec) else "existing_reaction"

    if mode == "existing_reaction":
        rxn_id = spec.get("target_reaction")
        if not rxn_id or not isinstance(rxn_id, str):
            raise ValueError("existing_reaction mode requires a non-empty 'target_reaction'")
        if rxn_id not in model.reactions:
            raise ValueError(
                f"Configured target_reaction {rxn_id!r} is absent from the model"
            )
        return {
            "mode": "existing_reaction",
            "product_reaction": rxn_id,
            "added_metabolites": [],
            "added_reactions": [],
        }

    if mode != "custom":
        raise ValueError(
            f"Unsupported target module mode {mode!r}; expected 'existing_reaction' or 'custom'"
        )

    rxn_id = spec.get("target_reaction")
    if not rxn_id or not isinstance(rxn_id, str):
        raise ValueError("Custom target module requires a non-empty 'target_reaction'")

    for prefix in spec.get("forbid_reaction_prefixes", []):
        if any(r.id.startswith(prefix) for r in model.reactions):
            raise ValueError(
                f"Model already contains a reaction starting with forbidden prefix {prefix!r}"
            )
    if rxn_id in model.reactions:
        raise ValueError(f"Target reaction {rxn_id!r} already exists in the model")

    for mid, expected in (spec.get("expected_metabolites") or {}).items():
        if mid not in model.metabolites:
            raise ValueError(f"Missing expected host metabolite {mid}")
        met = model.metabolites.get_by_id(mid)
        exp_formula = expected.get("formula")
        exp_charge = expected.get("charge")
        exp_comp = expected.get("compartment")
        if exp_formula is not None and met.formula != exp_formula:
            raise ValueError(
                f"Product chemistry contract mismatch on {mid}: formula {met.formula!r} != {exp_formula!r}"
            )
        if exp_charge is not None and met.charge != exp_charge:
            raise ValueError(
                f"Product chemistry contract mismatch on {mid}: charge {met.charge!r} != {exp_charge!r}"
            )
        if exp_comp is not None and met.compartment != exp_comp:
            raise ValueError(
                f"Product chemistry contract mismatch on {mid}: compartment {met.compartment!r} != {exp_comp!r}"
            )

    added_metabolites: list[str] = []
    new_met_objs: list[cobra.Metabolite] = []
    for met_spec in spec.get("new_metabolites", []):
        mid = met_spec["id"]
        if mid in model.metabolites:
            raise ValueError(f"Product bookkeeping metabolite already exists: {mid}")
        met = cobra.Metabolite(
            mid,
            formula=met_spec.get("formula"),
            charge=met_spec.get("charge"),
            compartment=met_spec.get("compartment", "c"),
            name=met_spec.get("name", mid),
        )
        new_met_objs.append(met)
        added_metabolites.append(mid)
    if new_met_objs:
        model.add_metabolites(new_met_objs)

    added_reactions: list[str] = []
    new_rxn_objs: list[cobra.Reaction] = []
    balance_checks: dict[str, dict[str, float]] = {}
    for rxn_spec in spec.get("reactions", []):
        rid = rxn_spec["id"]
        if rid in model.reactions:
            raise ValueError(f"Product module reaction already exists: {rid}")
        rxn = cobra.Reaction(
            rid,
            name=rxn_spec.get("name", rid),
            subsystem=rxn_spec.get("subsystem", ""),
            lower_bound=float(rxn_spec.get("lower_bound", 0.0)),
            upper_bound=float(rxn_spec.get("upper_bound", 1000.0)),
        )
        stoich = rxn_spec.get("stoichiometry") or {}
        if not stoich:
            raise ValueError(f"Reaction {rid!r} has empty stoichiometry")
        met_map: dict[cobra.Metabolite, float] = {}
        for mid, coeff in stoich.items():
            if mid not in model.metabolites:
                raise ValueError(f"Reaction {rid!r} references missing metabolite {mid!r}")
            met_map[model.metabolites.get_by_id(mid)] = float(coeff)
        rxn.add_metabolites(met_map)
        if rxn_spec.get("gene_reaction_rule"):
            rxn.gene_reaction_rule = str(rxn_spec["gene_reaction_rule"])
        for k, v in (rxn_spec.get("annotation") or {}).items():
            rxn.annotation[str(k)] = str(v)

        # Internal (>=2 metabolites) reactions check mass/charge balance by default;
        # single-metabolite demand/sink reactions default to False unless requested.
        default_check = len(met_map) >= 2
        if bool(rxn_spec.get("check_mass_balance", default_check)):
            imbalance = rxn.check_mass_balance()
            if imbalance:
                raise ValueError(f"Unbalanced target module reaction {rid!r}: {imbalance}")
            balance_checks[rid] = imbalance

        new_rxn_objs.append(rxn)
        added_reactions.append(rid)

    if new_rxn_objs:
        model.add_reactions(new_rxn_objs)

    if rxn_id not in model.reactions:
        raise ValueError(
            f"Target reaction {rxn_id!r} was not added by custom target module and is absent from model"
        )

    receipt: dict[str, Any] = {
        "mode": "custom",
        "product_reaction": rxn_id,
        "added_metabolites": added_metabolites,
        "added_reactions": added_reactions,
        "mass_charge_balance": balance_checks,
    }
    if "metadata" in spec and isinstance(spec["metadata"], Mapping):
        receipt.update(spec["metadata"])
    return receipt




def _known_carbon(metabolite: cobra.Metabolite) -> float | None:
    formula = metabolite.formula
    if not formula or not re.fullmatch(r"(?:[A-Z][a-z]?(?:[0-9]+(?:\.[0-9]+)?)?)+", formula):
        return None
    if any(symbol in {"R", "X"} for symbol in re.findall(r"[A-Z][a-z]?", formula)):
        return None
    elements = metabolite.elements
    if elements is None:
        return None
    return float(elements.get("C", 0))


def carbon_witness(
    model: cobra.Model,
    fluxes: Mapping[str, float],
    *,
    product_reaction: str,
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Count boundary carbon and flag active internal known-carbon defects."""
    carbon_in = carbon_out = 0.0
    unknown_boundary = []
    boundary_rows = []
    active_known_carbon_defects = []
    active_unknown = []
    for reaction in model.reactions:
        value = float(fluxes[reaction.id])
        if not math.isfinite(value):
            raise ValueError("Nonfinite carbon witness")
        if abs(value) <= tolerance:
            continue
        missing = [m.id for m in reaction.metabolites if _known_carbon(m) is None]
        if reaction.boundary:
            for metabolite, coefficient in reaction.metabolites.items():
                amount = float(coefficient) * value
                carbons = _known_carbon(metabolite)
                if carbons is None:
                    unknown_boundary.append(
                        {"reaction": reaction.id, "metabolite": metabolite.id, "amount": amount}
                    )
                    continue
                carbon = carbons * amount
                carbon_in += max(0.0, carbon)
                carbon_out += max(0.0, -carbon)
                boundary_rows.append(
                    {
                        "reaction": reaction.id,
                        "metabolite": metabolite.id,
                        "flux": value,
                        "carbon_into_model": carbon,
                    }
                )
        elif missing:
            active_unknown.append(
                {
                    "reaction": reaction.id,
                    "missing_formula_metabolites": missing,
                    "flux": value,
                    "name": reaction.name,
                }
            )
        else:
            residual = sum(
                float(coefficient) * (_known_carbon(met) or 0.0)
                for met, coefficient in reaction.metabolites.items()
            )
            if abs(residual) > tolerance:
                active_known_carbon_defects.append(
                    {
                        "reaction": reaction.id,
                        "carbon_residual": residual,
                        "flux": value,
                        "carbon_creation_rate": residual * value,
                    }
                )
    product_flux = float(fluxes.get(product_reaction, 0.0))
    product_rxn = model.reactions.get_by_id(product_reaction)
    product_carbons = sum(
        abs(float(coeff)) * (_known_carbon(met) or 0.0)
        for met, coeff in product_rxn.metabolites.items()
        if float(coeff) < 0
    )
    product_c_out = product_carbons * product_flux
    return {
        "known_boundary_carbon_input": carbon_in,
        "known_boundary_carbon_output": carbon_out,
        "product_carbon_output": product_c_out,
        "product_within_known_input": product_c_out <= carbon_in + tolerance,
        "known_output_within_known_input": carbon_out <= carbon_in + tolerance,
        "boundary_rows": boundary_rows,
        "active_unknown_formula_boundaries": unknown_boundary,
        "active_internal_known_carbon_defects": active_known_carbon_defects,
        "active_internal_unknown_formula_reactions": active_unknown,
        "global_elemental_conservation_proven": False,
        "note": "Mixed-unit pseudocomponents or missing formulas are reported separately and never treated as zero carbon",
    }
