"""Tests for carbon-accounting helpers on caller-provided targets."""
import cobra
import pytest

from geno2flux.target_module import carbon_witness


def model():
    m = cobra.Model("chemistry")
    product = cobra.Metabolite("product_c", formula="C2H4O2", charge=0, compartment="c")
    sink = cobra.Reaction("sink_product")
    sink.add_metabolites({product: -1})
    m.add_reactions([sink])
    return m


def test_carbon_audit_detects_known_carbon_creation_and_unknown_formula():
    m = model()
    a = cobra.Metabolite("a_c", formula="C", compartment="c")
    b = cobra.Metabolite("b_c", formula="C2", compartment="c")
    broken = cobra.Reaction("broken")
    broken.add_metabolites({a: -1, b: 1})
    u = cobra.Metabolite("unknown_c", compartment="c")
    unknown = cobra.Reaction("unknown")
    unknown.add_metabolites({u: -1, a: 1})
    m.add_reactions([broken, unknown])
    fluxes = {r.id: 0 for r in m.reactions}
    fluxes.update(broken=2, unknown=3)
    q = carbon_witness(m, fluxes, product_reaction="sink_product")
    assert q["active_internal_known_carbon_defects"][0]["carbon_creation_rate"] == 2
    assert q["active_internal_unknown_formula_reactions"][0]["reaction"] == "unknown"
    assert not q["global_elemental_conservation_proven"]


@pytest.mark.parametrize("formula", ["R", "C11H16O13PR2", "C*", "X", "carbon"])
def test_symbolic_and_invalid_formulas_are_unknown_not_zero_carbon(formula):
    m = model()
    a = cobra.Metabolite("u_c", formula=formula, compartment="c")
    boundary = cobra.Reaction("sink_u")
    boundary.add_metabolites({a: -1})
    m.add_reactions([boundary])
    flux = {r.id: 0 for r in m.reactions}
    flux["sink_u"] = -1
    q = carbon_witness(m, flux, product_reaction="sink_product")
    assert q["active_unknown_formula_boundaries"]
    assert q["known_boundary_carbon_input"] == 0
