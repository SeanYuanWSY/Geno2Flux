"""Synthetic-fixture tests for automatic compartment refinement."""

from __future__ import annotations

import csv
from pathlib import Path

import cobra
import pytest

from geno2flux.localization import (
    DL_TO_COMP,
    load_localization_table,
    refine_compartments,
    rxn_chem_sig,
    strip_comp,
)


def _met(mid, compartment):
    return cobra.Metabolite(mid, compartment=compartment)


def _rxn(rid, left, right, gpr="", lower=-1000.0, upper=1000.0):
    r = cobra.Reaction(rid, lower_bound=lower, upper_bound=upper)
    r.add_metabolites({m: c for m, c in list(left) + list(right)})
    r.gene_reaction_rule = gpr
    return r


def _mdh_template():
    """Yeast template carrying MDH (rhea_21432 family) in c/m/p."""
    ref = cobra.Model("yeast")
    rxns = []
    for comp, gene in (("c", "G_c"), ("m", "G_m"), ("p", "G_p")):
        rxns.append(
            _rxn(
                f"rhea_21432_{comp}",
                [
                    (_met(f"chebi15589_{comp}", comp), -1),
                    (_met(f"chebi57540_{comp}", comp), -1),
                ],
                [
                    (_met(f"chebi16452_{comp}", comp), 1),
                    (_met(f"chebi24636_{comp}", comp), 1),
                    (_met(f"chebi57945_{comp}", comp), 1),
                ],
                gpr=gene,
            )
        )
    ref.add_reactions(rxns)
    return ref


def _target_with_mito_mdh():
    """Target model that only has the mitochondrial MDH member."""
    model = cobra.Model("target")
    rxn = _rxn(
        "rhea_21432_m",
        [(_met("chebi15589_m", "m"), -1), (_met("chebi57540_m", "m"), -1)],
        [
            (_met("chebi16452_m", "m"), 1),
            (_met("chebi24636_m", "m"), 1),
            (_met("chebi57945_m", "m"), 1),
        ],
        gpr="ANV010759 or ANV010130 or ANV008085 or ANV005681",
    )
    model.add_reactions([rxn])
    return model


def test_strip_comp_suffix():
    assert strip_comp("chebi15361_c") == "chebi15361"
    assert strip_comp("chebi16004_erm") == "chebi16004"
    assert strip_comp("plainmet") == "plainmet"


def test_rxn_chem_sig_compartment_free():
    a = _rxn("x_c", [(_met("m1_c", "c"), -1)], [(_met("m2_c", "c"), 1)])
    b = _rxn("x_m", [(_met("m1_m", "m"), -1)], [(_met("m2_m", "m"), 1)])
    assert rxn_chem_sig(a) == rxn_chem_sig(b)


def test_pulls_cytosolic_mdh_with_deeploc_evidence():
    ref = _mdh_template()
    target = _target_with_mito_mdh()
    # DeepLoc predicts ANV010759 + ANV010130 as Cytoplasm -> comp c.
    localization = {
        "ANV010759": {"c"},
        "ANV010130": {"c"},
        "ANV008085": {"m"},
        "ANV005681": {"m"},
    }
    audit = refine_compartments(target, ref, localization)
    pulled_ids = {p["reaction_id"] for p in audit["pulled_reactions"]}
    # c pulled (cytosolic genes present); p pulled only if a peroxisomal gene
    # existed in the family -- it does not, so exactly one pull.
    assert "rhea_21432_c" in pulled_ids
    assert "rhea_21432_p" not in pulled_ids
    rxn = target.reactions.get_by_id("rhea_21432_c")
    assert set(g.id for g in rxn.genes) == {"ANV010759", "ANV010130"}
    assert audit["pulled_reaction_count"] == len(pulled_ids)


def test_no_pull_without_deeploc_evidence():
    ref = _mdh_template()
    target = _target_with_mito_mdh()
    localization = {g: {"m"} for g in ("ANV010759", "ANV010130", "ANV008085", "ANV005681")}
    audit = refine_compartments(target, ref, localization)
    assert audit["pulled_reaction_count"] == 0


def test_no_pull_when_compartment_already_present():
    ref = _mdh_template()
    target = _target_with_mito_mdh()
    # add the cytosolic member already so c is not missing
    target.add_reactions(
        [
            _rxn(
                "rhea_21432_c",
                [(_met("chebi15589_c", "c"), -1), (_met("chebi57540_c", "c"), -1)],
                [
                    (_met("chebi16452_c", "c"), 1),
                    (_met("chebi24636_c", "c"), 1),
                    (_met("chebi57945_c", "c"), 1),
                ],
                gpr="ANV010759",
            )
        ]
    )
    localization = {"ANV010759": {"c"}, "ANV010130": {"c"}}
    audit = refine_compartments(target, ref, localization)
    assert all(p["reaction_id"] != "rhea_21432_c" for p in audit["pulled_reactions"])


def test_load_localization_table(tmp_path):
    tsv = tmp_path / "loc.tsv"
    with tsv.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["protein_id", "deeploc_localizations"])
        w.writerow(["ANV1", "Cytoplasm|Nucleus"])
        w.writerow(["ANV2", "Mitochondrion"])
        w.writerow(["ANV3", "Extracellular"])
    table = load_localization_table(tsv)
    assert table["ANV1"] == {"c", "n"}
    assert table["ANV2"] == {"m", "mm"}
    assert table["ANV3"] == {"e", "ce"}


def test_prune_relocated_genes_removes_moved_gene():
    ref = _mdh_template()
    target = _target_with_mito_mdh()
    localization = {"ANV010759": {"c"}, "ANV010130": {"c"}, "ANV008085": {"m"}, "ANV005681": {"m"}}
    refine_compartments(target, ref, localization, prune_relocated_genes=True)
    mito = target.reactions.get_by_id("rhea_21432_m")
    assert "ANV010759" not in mito.gene_reaction_rule
    assert "ANV010130" not in mito.gene_reaction_rule
    assert "ANV008085" in mito.gene_reaction_rule


def test_partial_complex_localization_does_not_pull_a_reaction():
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "A and B"
    audit = refine_compartments(target, _mdh_template(), {"A": {"c"}, "B": {"m"}})
    assert audit["pulled_reaction_count"] == 0
    assert "rhea_21432_c" not in target.reactions


def test_pulled_complex_preserves_required_subunits_and_alternatives():
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "(A and B) or C"
    audit = refine_compartments(
        target, _mdh_template(), {"A": {"c"}, "B": {"c"}, "C": {"c"}},
    )
    pulled = target.reactions.get_by_id("rhea_21432_c")
    assert pulled.gpr.eval({"C"})
    assert not pulled.gpr.eval({"A", "C"})
    assert not pulled.gpr.eval({"B", "C"})
    assert pulled.gpr.eval({"A", "B"})  # the C isoenzyme remains sufficient
    assert audit["pulled_reactions"][0]["trigger_genes"] == ["A", "B", "C"]


def test_incomplete_complex_does_not_leak_subunits_into_pulled_gpr():
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "(A and B) or C"
    audit = refine_compartments(
        target, _mdh_template(), {"A": {"c"}, "B": {"m"}, "C": {"c"}},
        prune_relocated_genes=True,
    )
    pulled = target.reactions.get_by_id("rhea_21432_c")
    resident = target.reactions.get_by_id("rhea_21432_m")
    assert {gene.id for gene in pulled.genes} == {"C"}
    assert not pulled.gpr.eval({"C"})
    assert audit["pulled_reactions"][0]["trigger_genes"] == ["C"]
    assert {gene.id for gene in resident.genes} == {"A", "B"}
    assert not resident.gpr.eval({"A"})
    assert not resident.gpr.eval({"B"})


def test_pruning_keeps_resident_complex_boolean_logic():
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "(A and B) or (C and D)"
    refine_compartments(
        target, _mdh_template(),
        {"A": {"c"}, "B": {"c"}, "C": {"m"}, "D": {"m"}},
        prune_relocated_genes=True,
    )
    resident = target.reactions.get_by_id("rhea_21432_m")
    assert {gene.id for gene in resident.genes} == {"C", "D"}
    assert resident.gpr.eval(set())
    assert not resident.gpr.eval({"C"})
    assert not resident.gpr.eval({"D"})
    pulled = target.reactions.get_by_id("rhea_21432_c")
    assert {gene.id for gene in pulled.genes} == {"A", "B"}
    assert not pulled.gpr.eval({"A"})
    assert not pulled.gpr.eval({"B"})


def test_pruning_retains_genes_predicted_in_source_and_destination():
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "A or B"
    refine_compartments(
        target, _mdh_template(), {"A": {"c", "m"}, "B": {"m"}},
        prune_relocated_genes=True,
    )
    resident = target.reactions.get_by_id("rhea_21432_m")
    assert {gene.id for gene in resident.genes} == {"A", "B"}
    assert resident.gpr.eval({"B"})


def test_pruning_does_not_erase_evidence_for_second_missing_compartment():
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "(A and B) or C"
    audit = refine_compartments(
        target, _mdh_template(), {"A": {"c", "p"}, "B": {"c", "p"}, "C": {"m"}},
        prune_relocated_genes=True,
    )
    assert audit["pulled_reaction_count"] == 2
    for compartment in ("c", "p"):
        pulled = target.reactions.get_by_id(f"rhea_21432_{compartment}")
        assert {gene.id for gene in pulled.genes} == {"A", "B"}
        assert not pulled.gpr.eval({"A"})
        assert not pulled.gpr.eval({"B"})
    assert target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule == "C"


def test_raw_deeploc_fasta_identifier_matches_target_gene(tmp_path):
    path = tmp_path / "deeploc.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Protein_ID", "Localizations"])
        writer.writerow(["sp|P1|ENTRY_NAME", "Cytoplasm"])
    target = _target_with_mito_mdh()
    target.reactions.get_by_id("rhea_21432_m").gene_reaction_rule = "P1"
    localization = load_localization_table(path)
    assert localization == {"P1": {"c"}}
    audit = refine_compartments(target, _mdh_template(), localization)
    assert audit["pulled_reaction_count"] == 1
    assert audit["pulled_reactions"][0]["trigger_genes"] == ["P1"]
