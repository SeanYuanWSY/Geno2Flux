"""Tests for caller-selected target reactions and custom product modules."""
from __future__ import annotations

import json
from pathlib import Path

import cobra
import pytest

from geno2flux.target_module import (
    carbon_witness,
    load_target_spec,
    mount_target_module,
)


CUSTOM_SPEC = {
    "mode": "custom",
    "target_reaction": "DM_target_unit_c",
    "expected_metabolites": {
        "source_c": {"formula": "C2H4O2", "charge": 0, "compartment": "c"},
    },
    "new_metabolites": [
        {
            "id": "target_unit_c",
            "formula": "C2H4O2",
            "charge": 0,
            "compartment": "c",
        }
    ],
    "reactions": [
        {
            "id": "TARGET_UNIT_SYN",
            "stoichiometry": {"source_c": -1.0, "target_unit_c": 1.0},
        },
        {
            "id": "DM_target_unit_c",
            "stoichiometry": {"target_unit_c": -1.0},
        },
    ],
}


def _host_model() -> cobra.Model:
    m = cobra.Model("host")
    m.add_metabolites([
        cobra.Metabolite("source_c", formula="C2H4O2", charge=0, compartment="c"),
        cobra.Metabolite("product_e", formula="C2H4O2", charge=0, compartment="e"),
    ])
    ex_product = cobra.Reaction("EX_product_e", lower_bound=0.0, upper_bound=1000.0)
    ex_product.add_metabolites({m.metabolites.get_by_id("product_e"): -1.0})
    m.add_reactions([ex_product])
    return m


def test_existing_reaction_mode_requires_no_new_reactions():
    m = _host_model()
    before_rxns = len(m.reactions)
    receipt = mount_target_module(m, target_reaction="EX_product_e")
    assert receipt["mode"] == "existing_reaction"
    assert receipt["product_reaction"] == "EX_product_e"
    assert len(m.reactions) == before_rxns


def test_existing_reaction_mode_rejects_missing_reaction():
    m = _host_model()
    with pytest.raises(ValueError, match="absent from the model"):
        mount_target_module(m, target_reaction="DOES_NOT_EXIST")


def test_custom_module_accepts_inline_mapping_and_json_file(tmp_path: Path):
    m1 = _host_model()
    r1 = mount_target_module(m1, CUSTOM_SPEC)
    assert r1["mode"] == "custom"
    assert r1["product_reaction"] == "DM_target_unit_c"
    assert not m1.reactions.get_by_id("TARGET_UNIT_SYN").check_mass_balance()

    spec_file = tmp_path / "custom.json"
    spec_file.write_text(json.dumps(CUSTOM_SPEC), encoding="utf-8")
    m2 = _host_model()
    r2 = mount_target_module(m2, spec_file)
    assert r2["product_reaction"] == "DM_target_unit_c"
    assert load_target_spec(spec_file)["reactions"] == CUSTOM_SPEC["reactions"]


def test_missing_target_does_not_assume_any_product():
    with pytest.raises(ValueError, match="A target is required"):
        load_target_spec(None)
    with pytest.raises(ValueError, match="A target is required"):
        mount_target_module(_host_model())


def test_custom_module_enforces_mass_and_charge_balance(tmp_path: Path):
    spec_path = tmp_path / "custom.json"
    spec_path.write_text(json.dumps(CUSTOM_SPEC), encoding="utf-8")
    m = _host_model()
    receipt = mount_target_module(m, spec_path)
    assert receipt["mode"] == "custom"
    assert receipt["product_reaction"] == "DM_target_unit_c"
    assert "TARGET_UNIT_SYN" in m.reactions

    unbalanced = dict(CUSTOM_SPEC)
    unbalanced["new_metabolites"] = [
        {"id": "target_unit_c", "formula": "C3H4O2", "charge": 0, "compartment": "c"}
    ]
    m_bad = _host_model()
    with pytest.raises(ValueError, match="Unbalanced target module reaction"):
        mount_target_module(m_bad, unbalanced)


def test_carbon_witness_works_for_arbitrary_product_reaction():
    m = _host_model()
    fluxes = {r.id: 0.0 for r in m.reactions}
    fluxes["EX_product_e"] = 2.5
    witness = carbon_witness(m, fluxes, product_reaction="EX_product_e")
    assert witness["product_carbon_output"] == pytest.approx(5.0)
