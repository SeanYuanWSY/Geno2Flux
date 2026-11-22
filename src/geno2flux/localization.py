"""Automatic subcellular-localization compartment refinement for a GEM.

DeepLoc predictions are compared against the reaction-compartment families
present in the yeast template (``yeast-GEM.xml``).  When a single-compartment
reaction family exists in the template in >= 2 compartments, and the target
model already carries that family in some compartments ``a_comps`` but lacks a
template compartment ``mc``, the corresponding template reaction is pulled in
automatically -- provided a complete enzyme alternative from the target
reaction family's GPR is predicted by DeepLoc to localize to ``mc``.
AND subunits remain required in both pulled and pruned reaction rules.

No reaction ID is ever hard-coded; the rule is purely signature- and
localization-driven.
"""
from __future__ import annotations

import ast
import csv
import os
import re
import shlex
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import cobra

from .paralog_grouping import protein_key

# DeepLoc compartment label -> yeast-GEM compartment code(s).
DL_TO_COMP: dict[str, set[str]] = {
    "Cytoplasm": {"c"},
    "Mitochondrion": {"m", "mm"},
    "Peroxisome": {"p"},
    "Nucleus": {"n"},
    "Endoplasmic reticulum": {"er", "erm"},
    "Golgi apparatus": {"g", "gm"},
    "Lysosome/Vacuole": {"v", "vm"},
    "Cell membrane": {"ce"},
    "Extracellular": {"e", "ce"},
    "Plastid": set(),
}

_REPO_ROOT = Path(__file__).resolve().parents[2]


def default_deeploc_slot_paths(repo_root: Path | None = None) -> tuple[Path, Path]:
    """Return the two standard plug-and-play DeepLoc binary slots in the repository."""
    root = (repo_root or _REPO_ROOT).resolve()
    return (
        root / "tools" / "deeploc" / "bin" / "deeploc2",
        root / "tools" / "deeploc" / "deeploc2",
    )


def _is_executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def resolve_deeploc_command(
    configured_cmd: Sequence[str] | str | None = None,
    *,
    env: Mapping[str, str] | None = None,
    repo_root: Path | None = None,
) -> tuple[list[str] | None, str]:
    """Locate an executable DeepLoc CLI command or return ``(None, 'not_found')``.

    Search order:
    1. Explicit ``configured_cmd`` (e.g. ``deeploc_cmd`` in JSON config or
       ``--deeploc-bin`` when non-default). If explicitly configured by the
       caller, its leading executable must exist or ``FileNotFoundError`` is raised.
    2. ``G2F_DEEPLOC_BIN`` environment variable (command string or binary path).
       If set to a missing binary, ``FileNotFoundError`` is raised.
    3. Reserved drop-in slot ``<repo>/tools/deeploc/bin/deeploc2`` or
       ``<repo>/tools/deeploc/deeploc2``.
    4. ``deeploc2`` on ``PATH``.
    """
    environ = os.environ if env is None else env

    if configured_cmd is not None:
        if isinstance(configured_cmd, str):
            tokens = shlex.split(configured_cmd)
        else:
            tokens = [str(t) for t in configured_cmd]
        if tokens:
            head = tokens[0]
            head_path = Path(head).expanduser()
            if _is_executable(head_path):
                return [str(head_path.resolve()), *tokens[1:]], "configured_cmd"
            found = shutil.which(head, path=environ.get("PATH"))
            if found:
                return [found, *tokens[1:]], "configured_cmd"
            raise FileNotFoundError(
                f"Configured DeepLoc executable was not found or is not executable: {head!r}"
            )

    env_val = (environ.get("G2F_DEEPLOC_BIN") or "").strip()
    if env_val:
        tokens = shlex.split(env_val)
        if tokens:
            head = tokens[0]
            head_path = Path(head).expanduser()
            if _is_executable(head_path):
                return [str(head_path.resolve()), *tokens[1:]], "env_G2F_DEEPLOC_BIN"
            found = shutil.which(head, path=environ.get("PATH"))
            if found:
                return [found, *tokens[1:]], "env_G2F_DEEPLOC_BIN"
            raise FileNotFoundError(
                f"G2F_DEEPLOC_BIN points to a missing or non-executable binary: {env_val!r}"
            )

    for slot in default_deeploc_slot_paths(repo_root):
        if _is_executable(slot):
            return [str(slot)], "repo_tools_deeploc_slot"

    found_path = shutil.which("deeploc2", path=environ.get("PATH"))
    if found_path:
        return [found_path], "system_PATH"

    return None, "not_found"

# Recognised compartment suffixes on metabolite IDs.
_COMP_SUFFIX = re.compile(r"_(c|m|p|e|n|er|g|v|ce|erm|lp|vm|mm|gm)$")


def strip_comp(metabolite_id: str) -> str:
    """Drop the trailing compartment suffix from a metabolite ID."""
    return _COMP_SUFFIX.sub("", metabolite_id)


def rxn_chem_sig(rxn: cobra.Reaction) -> tuple:
    """Compartment-free stoichiometric signature of a reaction.

    Two reactions that catalyse the same chemistry in different compartments
    share the same signature.
    """
    return tuple(
        sorted(
            (strip_comp(met.id), round(coeff, 6))
            for met, coeff in rxn.metabolites.items()
        )
    )


def _rxn_compartments(rxn: cobra.Reaction) -> set[str]:
    comps: set[str] = set()
    for met in rxn.metabolites:
        match = _COMP_SUFFIX.search(met.id)
        if match:
            comps.add(match.group(1))
        else:
            # Fall back to the declared compartment attribute when the
            # metabolite ID carries no recognisable suffix.
            comps.add(met.compartment or "")
    return comps


def load_localization_table(path: Path) -> dict[str, set[str]]:
    """Map protein_id -> set of yeast-GEM compartment codes predicted by DeepLoc.

    Accepts the merged ``localization_tsv`` (``deeploc_localizations`` column,
    ``|``-separated DeepLoc labels) or a raw ``deeploc_csv`` with a
    ``Localizations``/``localization`` column. FASTA identifiers use the same
    accession key as AlphaGEM and the FSEOF proteome report.
    """
    mapping: dict[str, set[str]] = {}
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Localization table has no header: {path}")
        fields = {f.lower(): f for f in reader.fieldnames}
        loc_col = None
        for candidate in (
            "deeploc_localizations",
            "localizations",
            "localization",
            "deeploc_localization",
        ):
            if candidate in fields:
                loc_col = fields[candidate]
                break
        id_col = None
        for candidate in ("protein_id", "protein", "id", "gene_id", "seq_id"):
            if candidate in fields:
                id_col = fields[candidate]
                break
        if loc_col is None or id_col is None:
            raise ValueError(
                f"Localization table lacks protein_id/localizations columns: {path}"
            )
        for row in reader:
            pid = (row.get(id_col) or "").strip()
            if not pid:
                continue
            pid = protein_key(pid)
            comps: set[str] = set()
            for label in re.split(r"[|;]", row.get(loc_col) or ""):
                label = label.strip()
                if label:
                    comps |= DL_TO_COMP.get(label, set())
            mapping[pid] = comps
    if not mapping:
        raise ValueError(f"Localization table yielded no protein mappings: {path}")
    return mapping


def _single_compartment_sig_families(
    model: cobra.Model,
) -> dict[tuple, dict[str, list[cobra.Reaction]]]:
    """Group single-compartment GPR-bearing reactions by chemical signature.

    Returns ``{sig: {compartment: [reaction, ...]}}`` keeping only families that
    span >= 2 compartments.
    """
    families: dict[tuple, dict[str, list[cobra.Reaction]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for rxn in model.reactions:
        comps = _rxn_compartments(rxn)
        if len(comps) != 1:
            continue  # only single-compartment reactions are family anchors
        if not rxn.gene_reaction_rule.strip():
            continue  # template anchor must carry a GPR
        families[rxn_chem_sig(rxn)][next(iter(comps))].append(rxn)
    return {sig: comps for sig, comps in families.items() if len(comps) >= 2}


def _model_sig_index(
    model: cobra.Model,
) -> dict[tuple, dict[str, set[str]]]:
    """Index target single-compartment reactions: ``sig -> comp -> {gene ids}``."""
    index: dict[tuple, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for rxn in model.reactions:
        comps = _rxn_compartments(rxn)
        if len(comps) != 1:
            continue
        comp = next(iter(comps))
        index[rxn_chem_sig(rxn)][comp] |= {gene.id for gene in rxn.genes}
    return index


def _restricted_gpr(node: ast.AST | None, allowed: set[str]) -> str | None:
    """Keep complete GPR alternatives using only the allowed genes.

    A missing AND subunit removes the whole complex alternative; OR branches
    can survive independently. None means no functional alternative remains.
    """
    if node is None:
        return None
    if isinstance(node, ast.Name):
        return node.id if node.id in allowed else None
    if isinstance(node, ast.BoolOp):
        children = [_restricted_gpr(child, allowed) for child in node.values]
        if isinstance(node.op, ast.And):
            if any(child is None for child in children):
                return None
            operator = " and "
        elif isinstance(node.op, ast.Or):
            operator = " or "
        else:
            raise ValueError(f"Unsupported GPR operator: {type(node.op).__name__}")
        kept = sorted({child for child in children if child is not None})
        if not kept:
            return None
        if len(kept) == 1:
            return kept[0]
        return operator.join(f"({child})" for child in kept)
    raise ValueError(f"Unsupported GPR node: {type(node).__name__}")


def refine_compartments(
    model: cobra.Model,
    template: cobra.Model,
    localization: Mapping[str, set[str]],
    *,
    prune_relocated_genes: bool = False,
) -> dict[str, Any]:
    """Pull DeepLoc-supported missing-compartment reactions into ``model``.

    Returns an audit ledger (also JSON-serializable) describing every pulled
    reaction.  ``model`` is modified in place.
    """
    template_families = _single_compartment_sig_families(template)
    model_index = _model_sig_index(model)
    # Snapshot rules before optional pruning so a gene supported in multiple
    # missing compartments can seed each one with the original complex logic.
    source_rules: dict[tuple, list[cobra.core.GPR]] = defaultdict(list)
    for rxn in model.reactions:
        if len(_rxn_compartments(rxn)) == 1 and rxn.gene_reaction_rule.strip():
            source_rules[rxn_chem_sig(rxn)].append(rxn.gpr.copy())
    existing_ids = {rxn.id for rxn in model.reactions}
    orig_metabolites = {met.id for met in model.metabolites}
    pulled: list[dict[str, Any]] = []
    for sig, template_comps in template_families.items():
        if sig not in model_index:
            continue
        present_comps = model_index[sig]
        a_comps = set(present_comps)
        family_genes: set[str] = set().union(*present_comps.values())
        if not family_genes:
            continue
        for mc in sorted(template_comps):
            if mc in a_comps:
                continue  # compartment already present in target
            dl_genes = sorted(
                gene for gene in family_genes if mc in localization.get(gene, set())
            )
            if not dl_genes:
                continue  # no DeepLoc evidence for the missing compartment
            allowed = set(dl_genes)
            rules = sorted({
                rule for gpr in source_rules[sig]
                if (rule := _restricted_gpr(gpr.body, allowed)) is not None
            })
            if not rules:
                continue  # localization of one complex subunit is insufficient
            supported_rule = " or ".join(f"({rule})" for rule in rules)
            for template_rxn in sorted(
                template_comps[mc], key=lambda r: r.id
            ):
                if template_rxn.id in existing_ids:
                    continue
                new_rxn = template_rxn.copy()
                new_rxn.gene_reaction_rule = supported_rule
                trigger_genes = sorted(gene.id for gene in new_rxn.genes)
                # Computed BEFORE add_reactions: adding the reaction would
                # otherwise register its metabolites and always yield True.
                substrates_ok = all(
                    met.id in orig_metabolites for met in new_rxn.metabolites
                )
                model.add_reactions([new_rxn])
                existing_ids.add(new_rxn.id)
                pulled.append(
                    {
                        "reaction_id": new_rxn.id,
                        "reaction_name": new_rxn.name,
                        "source_compartments": sorted(a_comps),
                        "target_compartment": mc,
                        "trigger_genes": trigger_genes,
                        "target_substrates_present": substrates_ok,
                        "template_source_reaction": template_rxn.id,
                    }
                )
                if prune_relocated_genes:
                    _prune_genes(model_index, sig, mc, trigger_genes, model, localization)

    return {
        "pulled_reaction_count": len(pulled),
        "pulled_reactions": pulled,
        "prune_relocated_genes": prune_relocated_genes,
        "final_reaction_count": len(model.reactions),
    }


def _prune_genes(
    model_index: dict[tuple, dict[str, set[str]]],
    sig: tuple,
    moved_to: str,
    moved_genes: list[str],
    model: cobra.Model,
    localization: Mapping[str, set[str]],
) -> None:
    """Prune relocated alternatives only when a complete resident GPR survives.

    Genes predicted in both source and destination remain in the source rule.
    """
    resident = model_index[sig]
    for comp in resident:
        if comp == moved_to:
            continue
        moved = {gene for gene in moved_genes if comp not in localization.get(gene, set())}
        if not moved:
            continue
        for rxn in model.reactions:
            if rxn_chem_sig(rxn) != sig or _rxn_compartments(rxn) != {comp}:
                continue
            if not rxn.gene_reaction_rule.strip():
                continue
            allowed = {gene.id for gene in rxn.genes} - moved
            rule = _restricted_gpr(rxn.gpr.body, allowed)
            if rule is not None and allowed != {gene.id for gene in rxn.genes}:
                rxn.gene_reaction_rule = rule


def load_models_for_refinement(
    target_path: Path, template_path: Path
) -> tuple[cobra.Model, cobra.Model]:
    """Load target (YAML or SBML) and template (SBML) models."""
    target = _load_model(target_path)
    template = cobra.io.read_sbml_model(str(template_path))
    return target, template


def _load_model(path: Path) -> cobra.Model:
    suffix = path.suffix.lower()
    if suffix in {".yml", ".yaml"}:
        return cobra.io.load_yaml_model(str(path))
    if suffix == ".json":
        return cobra.io.load_json_model(str(path))
    return cobra.io.read_sbml_model(str(path))
