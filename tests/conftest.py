"""Shared pytest fixtures.

Sets the cobra default solver to a working backend so that ``cobra.Model()`` /
``model.copy()`` inside localization code never instantiate the machine's
possibly-broken default (e.g. expired gurobi) solver. HiGHS is preferred when
registered (as under the vendored AlphaGEM ``optlang_highs_interface``);
otherwise falls back to a listed free solver.

Also provides ``synthetic_screening_case(tmp_path)`` to materialize an
in-memory synthetic metabolic model + coding sequences inside pytest's
temporary directory without keeping any toy model directory in the repository.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from Bio.Data import CodonTable
import cobra
import pytest
from optlang import available_solvers

# Make the top-level ``scripts`` package importable for legacy test modules.
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))


@pytest.fixture(autouse=True)
def _default_solver() -> None:
    """Set the cobra default solver to a working backend for the test."""
    for name in ("highs", "glpk", "scipy"):
        if available_solvers.get(name.upper(), False):
            try:
                cobra.Configuration().solver = name
            except cobra.exceptions.SolverNotFound:
                continue
            break


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _formula(counts: list[tuple[str, int]]) -> str:
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in counts)


def _synth_aa(seed: int, *, mut: int = 0) -> str:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    chars = ["M"]
    for i in range(47):
        idx = (seed * 7 + i * 3 + (1 if (mut and i % 9 == 0) else 0)) % len(alphabet)
        chars.append(alphabet[idx])
    return "".join(chars)


def _build_synthetic_models() -> tuple[cobra.Model, cobra.Model]:
    """Construct balanced base & template models programmatically in memory."""
    f_atp = _formula([("C", 10), ("H", 12), ("N", 5), ("O", 13), ("P", 3)])
    f_adp = _formula([("C", 10), ("H", 12), ("N", 5), ("O", 10), ("P", 2)])
    f_amp = _formula([("C", 10), ("H", 12), ("N", 5), ("O", 7), ("P", 1)])
    f_nadh = _formula([("C", 21), ("H", 27), ("N", 7), ("O", 14), ("P", 2)])
    f_nad = _formula([("C", 21), ("H", 26), ("N", 7), ("O", 14), ("P", 2)])

    specs = [
        ("glc_e", "D-glucose [e]", "e", "C6H12O6", 0),
        ("glc_c", "D-glucose [c]", "c", "C6H12O6", 0),
        ("pyr_c", "Pyruvate [c]", "c", "C3H3O3", -1),
        ("co2_e", "CO2 [e]", "e", "CO2", 0),
        ("co2_c", "CO2 [c]", "c", "CO2", 0),
        ("hco3_c", "Bicarbonate [c]", "c", "CHO3", -1),
        ("h2o_e", "H2O [e]", "e", "H2O", 0),
        ("h2o_c", "H2O [c]", "c", "H2O", 0),
        ("h_e", "H+ [e]", "e", "H", 1),
        ("h_c", "H+ [c]", "c", "H", 1),
        ("o2_e", "O2 [e]", "e", "O2", 0),
        ("o2_c", "O2 [c]", "c", "O2", 0),
        ("oaa_c", "Oxaloacetate [c]", "c", "C4H2O5", -2),
        ("oaa_m", "Oxaloacetate [m]", "m", "C4H2O5", -2),
        ("nadh_c", "NADH [c]", "c", f_nadh, -2),
        ("nad_c", "NAD+ [c]", "c", f_nad, -1),
        ("nadh_m", "NADH [m]", "m", f_nadh, -2),
        ("nad_m", "NAD+ [m]", "m", f_nad, -1),
        ("h_m", "H+ [m]", "m", "H", 1),
        ("fum_c", "Fumarate [c]", "c", "C4H2O4", -2),
        ("chebi15589_c", "(S)-malate [c]", "c", "C4H4O5", -2),
        ("chebi15589_m", "(S)-malate [m]", "m", "C4H4O5", -2),
        ("chebi30616_c", "ATP [c]", "c", f_atp, -4),
        ("adp_c", "ADP [c]", "c", f_adp, -3),
        ("pi_c", "Phosphate [c]", "c", "HO4P", -2),
        ("chebi456215_c", "AMP [c]", "c", f_amp, -2),
        ("chebi33019_c", "Diphosphate [c]", "c", "HO7P2", -3),
        ("lac_c", "Lactate [c]", "c", "C3H5O3", -1),
        ("lac_e", "Lactate [e]", "e", "C3H5O3", -1),
        ("biomass_c", "Biomass precursor [c]", "c", "C3H3O3", -1),
    ]
    mets = {
        mid: cobra.Metabolite(mid, name=name, compartment=comp, formula=formula, charge=charge)
        for mid, name, comp, formula, charge in specs
    }

    def rxn(rid: str, name: str, stoich: dict[str, float], lb: float, ub: float, gpr: str = "") -> cobra.Reaction:
        r = cobra.Reaction(rid, name=name, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({mets[k]: v for k, v in stoich.items()})
        if gpr:
            r.gene_reaction_rule = gpr
        return r

    base = cobra.Model("synthetic_draft")
    base.add_reactions(
        [
            rxn("EX_glc_e", "Glucose exchange", {"glc_e": -1}, -10.0, 1000.0),
            rxn("EX_o2_e", "O2 exchange", {"o2_e": -1}, -1000.0, 1000.0),
            rxn("EX_co2_e", "CO2 exchange", {"co2_e": -1}, -1000.0, 1000.0),
            rxn("EX_h2o_e", "H2O exchange", {"h2o_e": -1}, -1000.0, 1000.0),
            rxn("EX_h_e", "H+ exchange", {"h_e": -1}, -1000.0, 1000.0),
            rxn("EX_lac_e", "Lactate exchange", {"lac_e": -1}, 0.0, 1000.0),
            rxn("GLCt", "Glucose transport", {"glc_e": -1, "glc_c": 1}, 0.0, 1000.0, "G_GLC"),
            rxn("O2t", "O2 transport", {"o2_e": -1, "o2_c": 1}, -1000.0, 1000.0),
            rxn("CO2t", "CO2 transport", {"co2_e": -1, "co2_c": 1}, -1000.0, 1000.0),
            rxn("H2Ot", "H2O transport", {"h2o_e": -1, "h2o_c": 1}, -1000.0, 1000.0),
            rxn("Ht", "H+ transport", {"h_e": -1, "h_c": 1}, -1000.0, 1000.0),
            rxn(
                "GLYC",
                "Glycolysis lumped",
                {
                    "glc_c": -1,
                    "nad_c": -2,
                    "adp_c": -2,
                    "pi_c": -2,
                    "pyr_c": 2,
                    "nadh_c": 2,
                    "chebi30616_c": 2,
                    "h2o_c": 2,
                    "h_c": 2,
                },
                0.0,
                1000.0,
                "G_GLC",
            ),
            rxn(
                "RESP",
                "Respiratory NADH oxidation",
                {
                    "nadh_c": -2,
                    "o2_c": -1,
                    "adp_c": -2,
                    "pi_c": -2,
                    "h_c": -4,
                    "nad_c": 2,
                    "chebi30616_c": 2,
                    "h2o_c": 4,
                },
                0.0,
                1000.0,
            ),
            rxn(
                "ATPM",
                "ATP turnover",
                {"chebi30616_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1, "h_c": 1},
                0.0,
                1000.0,
            ),
            rxn(
                "HCO3E",
                "Carbonic anhydrase",
                {"co2_c": -1, "h2o_c": -1, "hco3_c": 1, "h_c": 1},
                -1000.0,
                1000.0,
            ),
            rxn(
                "PYC_c",
                "Pyruvate carboxylase",
                {
                    "pyr_c": -1,
                    "hco3_c": -1,
                    "chebi30616_c": -1,
                    "oaa_c": 1,
                    "adp_c": 1,
                    "pi_c": 1,
                    "h_c": 1,
                },
                0.0,
                1000.0,
                "G_PYC1 or G_PYC2",
            ),
            rxn(
                "FUM_c",
                "Fumarase (cytosolic, reversible)",
                {"fum_c": -1, "h2o_c": -1, "chebi15589_c": 1},
                -1000.0,
                1000.0,
            ),
            rxn(
                "MDH_m",
                "Malate dehydrogenase (mitochondrial)",
                {"oaa_m": -1, "nadh_m": -1, "h_m": -1, "chebi15589_m": 1, "nad_m": 1},
                -1000.0,
                1000.0,
                "G_MDH1 or G_MDH2",
            ),
            rxn(
                "ADK1",
                "Adenylate kinase",
                {"chebi456215_c": -1, "chebi30616_c": -1, "adp_c": 2},
                -1000.0,
                1000.0,
                "G_ATP",
            ),
            rxn(
                "PPA",
                "Inorganic diphosphatase",
                {"chebi33019_c": -1, "h2o_c": -1, "pi_c": 2, "h_c": 1},
                0.0,
                1000.0,
            ),
            rxn(
                "LDH_c",
                "Lactate dehydrogenase",
                {"pyr_c": -1, "nadh_c": -1, "h_c": -1, "lac_c": 1, "nad_c": 1},
                0.0,
                1000.0,
                "G_FERM",
            ),
            rxn(
                "LACt",
                "Lactate export",
                {"lac_c": -1, "lac_e": 1},
                0.0,
                1000.0,
                "G_FERM",
            ),
            rxn(
                "BIOMASS_SYN",
                "Biomass synthesis",
                {
                    "pyr_c": -1,
                    "chebi30616_c": -1,
                    "h2o_c": -1,
                    "biomass_c": 1,
                    "adp_c": 1,
                    "pi_c": 1,
                    "h_c": 1,
                },
                0.0,
                1000.0,
                "G_BIO",
            ),
            rxn(
                "r_4041",
                "Biomass growth reaction",
                {"biomass_c": -1},
                0.0,
                1000.0,
                "G_BIO",
            ),
        ]
    )
    base.objective = "r_4041"

    template = base.copy()
    mdh_c = cobra.Reaction("MDH_c", name="Malate dehydrogenase (cytosolic)", lower_bound=-1000.0, upper_bound=1000.0)
    mdh_c.add_metabolites(
        {
            template.metabolites.get_by_id("oaa_c"): -1,
            template.metabolites.get_by_id("nadh_c"): -1,
            template.metabolites.get_by_id("h_c"): -1,
            template.metabolites.get_by_id("chebi15589_c"): 1,
            template.metabolites.get_by_id("nad_c"): 1,
        }
    )
    mdh_c.gene_reaction_rule = "TPL_MDH_C"
    template.add_reactions([mdh_c])
    return base, template


def build_synthetic_screening_dir(case_dir: Path) -> Path:
    """Write a complete synthetic fixture directory under ``case_dir`` for tests."""
    case_dir.mkdir(parents=True, exist_ok=True)
    proteins = {
        "G_ATP": _synth_aa(1),
        "G_BIO": _synth_aa(2),
        "G_FERM": _synth_aa(3),
        "G_GLC": _synth_aa(4),
        "G_MDH1": _synth_aa(5, mut=0),
        "G_MDH2": _synth_aa(5, mut=1),
        "G_PYC1": _synth_aa(6, mut=0),
        "G_PYC2": _synth_aa(6, mut=1),
    }
    back_table = CodonTable.unambiguous_dna_by_id[1].back_table
    flank = "A" * 40 + "C" * 20 + "G" * 20
    genome_chunks = [flank]
    gff_lines = ["##gff-version 3"]
    cursor = len(flank)
    for pid, seq in proteins.items():
        cds = "".join(back_table[aa] for aa in seq) + "TAA"
        start = cursor + 1
        end = cursor + len(cds)
        gff_lines.extend(
            [
                f"chr1\tTest\tgene\t{start}\t{end}\t.\t+\t.\tID=gene_{pid}",
                f"chr1\tTest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={pid};Parent=gene_{pid}",
                f"chr1\tTest\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds_{pid};Parent={pid}",
            ]
        )
        genome_chunks.append(cds)
        genome_chunks.append(flank)
        cursor = end + len(flank)

    genome_seq = "".join(genome_chunks)
    gff_lines.insert(1, f"##sequence-region chr1 1 {len(genome_seq)}")

    (case_dir / "genome.fasta").write_text(f">chr1\n{genome_seq}\n", encoding="utf-8")
    (case_dir / "annotation.gff3").write_text("\n".join(gff_lines) + "\n", encoding="utf-8")
    (case_dir / "proteome.fasta").write_text(
        "".join(f">{pid}\n{seq}\n" for pid, seq in proteins.items()),
        encoding="utf-8",
    )
    (case_dir / "localization.tsv").write_text(
        "protein_id\tdeeploc_localizations\n"
        "G_ATP\tCytoplasm\n"
        "G_BIO\tCytoplasm\n"
        "G_FERM\tCytoplasm\n"
        "G_GLC\tCytoplasm\n"
        "G_MDH1\tCytoplasm|Mitochondrion\n"
        "G_MDH2\tMitochondrion\n"
        "G_PYC1\tCytoplasm\n"
        "G_PYC2\tCytoplasm\n",
        encoding="utf-8",
    )
    L = len(next(iter(proteins.values())))
    (case_dir / "self_hits.tsv").write_text(
        f"G_ATP\tG_ATP\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_BIO\tG_BIO\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_FERM\tG_FERM\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_GLC\tG_GLC\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_MDH1\tG_MDH1\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_MDH1\tG_MDH2\t85.0\t{L}\t1e-50\t180\t{L}\t{L}\n"
        f"G_MDH2\tG_MDH2\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_MDH2\tG_MDH1\t85.0\t{L}\t1e-50\t180\t{L}\t{L}\n"
        f"G_PYC1\tG_PYC1\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_PYC1\tG_PYC2\t88.0\t{L}\t1e-55\t190\t{L}\t{L}\n"
        f"G_PYC2\tG_PYC2\t100\t{L}\t1e-60\t200\t{L}\t{L}\n"
        f"G_PYC2\tG_PYC1\t88.0\t{L}\t1e-55\t190\t{L}\t{L}\n",
        encoding="utf-8",
    )
    (case_dir / "target_module.json").write_text(
        json.dumps(
            {
                "mode": "custom",
                "target_reaction": "DM_target_unit_c",
                "expected_metabolites": {
                    "chebi15589_c": {"formula": "C4H4O5", "charge": -2, "compartment": "c"},
                    "chebi30616_c": {"formula": _formula([("C", 10), ("H", 12), ("N", 5), ("O", 13), ("P", 3)]), "charge": -4, "compartment": "c"},
                    "chebi456215_c": {"formula": _formula([("C", 10), ("H", 12), ("N", 5), ("O", 7), ("P", 1)]), "charge": -2, "compartment": "c"},
                    "chebi33019_c": {"formula": "HO7P2", "charge": -3, "compartment": "c"},
                },
                "new_metabolites": [
                    {"id": "target_unit_c", "formula": "C4H3O4", "charge": -1, "compartment": "c"},
                ],
                "reactions": [
                    {
                        "id": "TARGET_UNIT_SYN",
                        "stoichiometry": {
                            "chebi15589_c": -1.0,
                            "chebi30616_c": -1.0,
                            "chebi456215_c": 1.0,
                            "chebi33019_c": 1.0,
                            "target_unit_c": 1.0,
                        },
                    },
                    {
                        "id": "DM_target_unit_c",
                        "check_mass_balance": False,
                        "stoichiometry": {"target_unit_c": -1.0},
                    },
                ],
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    base, template = _build_synthetic_models()
    cobra.io.write_sbml_model(base, str(case_dir / "base_model.xml"))
    cobra.io.write_sbml_model(template, str(case_dir / "template_model.xml"))

    cfg = {
        "schema_version": 1,
        "name": "synthetic_case",
        "solver": "highs",
        "localization_order": "post-gapfill",
        "biomass": "r_4041",
        "target_module": "target_module.json",
        "compute_growth_essentiality": True,
        "essential_ko_loss_threshold": 0.5,
        "proteome_fasta": "proteome.fasta",
        "proteome_sha256": _sha256(case_dir / "proteome.fasta"),
        "localization_tsv": "localization.tsv",
        "localization_sha256": _sha256(case_dir / "localization.tsv"),
        "gapfilled_model": "base_model.xml",
        "gapfilled_model_sha256": _sha256(case_dir / "base_model.xml"),
        "template_model": "template_model.xml",
        "template_model_sha256": _sha256(case_dir / "template_model.xml"),
        "diamond_self_hits": "self_hits.tsv",
        "diamond_self_hits_sha256": _sha256(case_dir / "self_hits.tsv"),
        "scan": {
            "iterations": 10,
            "maximum_fraction": 0.9,
            "feasibility_tolerance": 1e-06,
            "monotonicity_tolerance": 1e-06,
        },
    }
    (case_dir / "config.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    return case_dir


@pytest.fixture
def synthetic_screening_case(tmp_path: Path) -> Path:
    return build_synthetic_screening_dir(tmp_path / "_synthetic_case")
