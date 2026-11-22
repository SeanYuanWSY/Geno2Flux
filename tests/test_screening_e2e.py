"""End-to-end integration tests using programmatic temporary fixtures.

Exercises both the unified pipeline (`geno2flux pipeline`) and the individual
stage CLI (`scripts/run_stage.py`) on real metabolic models, real HiGHS solves,
real compartment refinement, real FSEOF, real paralog clustering, real KO
suites, and real growth-essentiality filtering—without shipping any toy example
directory in the repository.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from geno2flux.pipeline import run_pipeline
import scripts.run_stage as run_stage


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_unified_pipeline_end_to_end(tmp_path: Path, synthetic_screening_case: Path):
    out_dir = tmp_path / "e2e_run"
    result = run_pipeline(synthetic_screening_case / "config.json", out_dir)
    assert "post-gapfill" in result["paths"]
    summary = result["paths"]["post-gapfill"]

    # 1. Localization pulled cytosolic MDH_c via G_MDH1 (Cytoplasm|Mitochondrion)
    assert summary["pulled"] == 1
    loc_audit = json.loads(
        (out_dir / "post-gapfill" / "localization_refinement_audit.json").read_text()
    )
    assert loc_audit["pulled_reactions"][0]["reaction_id"] == "MDH_c"
    assert loc_audit["pulled_reactions"][0]["trigger_genes"] == ["G_MDH1"]

    # 2. Essentiality filter dropped growth-essential G_GLC & G_BIO (loss == 1.0 >= 0.5)
    ef = summary["essential_filter"]
    assert ef["essentiality_source"] == "growth_ko_table_threshold"
    assert ef["essential_genes"] == 2

    # 3. Paralog compensation demonstrated:
    # Single-gene KO of G_PYC1 or G_PYC2 has loss == 0.0, while joint family KO has loss == 1.0!
    sg_ko = {
        r["gene_id"]: float(r["loss"])
        for r in _read_tsv(out_dir / "post-gapfill" / "knockout" / "single_gene_ko_essential_filtered.tsv")
    }
    assert "G_GLC" not in sg_ko
    assert "G_BIO" not in sg_ko
    assert sg_ko["G_PYC1"] == 0.0
    assert sg_ko["G_PYC2"] == 0.0
    assert sg_ko["G_MDH1"] == 1.0

    fam_ko = _read_tsv(out_dir / "post-gapfill" / "knockout" / "family_ko_essential_filtered.tsv")
    pyc_fam = next(r for r in fam_ko if r["members"] == "G_PYC1;G_PYC2")
    assert float(pyc_fam["loss"]) == 1.0

    # 4. FSEOF identifies PYC, MDH1, and ADK1 as top non-essential overexpression candidates
    fseof_filtered = _read_tsv(
        out_dir / "post-gapfill" / "fseof" / "all_genes_fseof_essential_filtered.tsv"
    )
    top_candidates = {r["gene_id"] for r in fseof_filtered if r["filtered_rank"] != ""}
    assert top_candidates == {"G_ATP", "G_MDH1", "G_PYC1", "G_PYC2"}


def test_screening_only_skip_localization_with_existing_reaction_target(
    tmp_path: Path, synthetic_screening_case: Path
):
    """Users with an existing SBML model and an endogenous target reaction can run
    FSEOF + paralogs + KO + essentiality filtering directly via localization_order='skip'."""
    cfg = json.loads((synthetic_screening_case / "config.json").read_text())
    cfg["localization_order"] = "skip"
    cfg.pop("target_module", None)
    cfg["target_reaction"] = "EX_lac_e"
    cfg_path = synthetic_screening_case / "lactate_screen.json"
    cfg_path.write_text(json.dumps(cfg, indent=2))

    out_dir = tmp_path / "lactate_out"
    result = run_pipeline(cfg_path, out_dir)
    summary = result["paths"]["screening"]
    assert summary["target_reaction"] == "EX_lac_e"
    fseof_rows = _read_tsv(out_dir / "screening" / "fseof" / "all_genes_fseof_essential_filtered.tsv")
    g_ferm = next(r for r in fseof_rows if r["gene_id"] == "G_FERM")
    assert g_ferm["filtered_rank"] == "1"


def test_run_stage_cli_end_to_end(tmp_path: Path, synthetic_screening_case: Path):
    """Exercise localize -> target -> fseof -> paralogs -> knockout -> filter via scripts/run_stage.py."""
    stage_dir = tmp_path / "stages"
    inputs_dir = synthetic_screening_case

    run_stage.main(
        [
            "localize",
            "--workdir", str(stage_dir),
            "--inputs", str(inputs_dir),
            "--draft-model", str(inputs_dir / "base_model.xml"),
            "--template-model", str(inputs_dir / "template_model.xml"),
            "--localization-table", str(inputs_dir / "localization.tsv"),
        ]
    )
    assert (stage_dir / "localized.yml").is_file()

    run_stage.main(
        [
            "target",
            "--workdir", str(stage_dir),
            "--inputs", str(inputs_dir),
            "--source-model", str(stage_dir / "localized.yml"),
            "--target-spec", str(inputs_dir / "target_module.json"),
        ]
    )
    assert (stage_dir / "product.xml").is_file()

    run_stage.main(
        [
            "fseof",
            "--workdir", str(stage_dir),
            "--inputs", str(inputs_dir),
            "--growth-reaction", "r_4041",
        ]
    )
    assert (stage_dir / "all_genes_fseof.tsv").is_file()

    run_stage.main(
        [
            "paralogs",
            "--workdir", str(stage_dir),
            "--inputs", str(inputs_dir),
            "--self-hits", str(inputs_dir / "self_hits.tsv"),
        ]
    )
    assert (stage_dir / "families.tsv").is_file()

    run_stage.main(
        [
            "knockout",
            "--workdir", str(stage_dir),
            "--inputs", str(inputs_dir),
        ]
    )
    assert (stage_dir / "knockout" / "family_ko.tsv").is_file()

    run_stage.main(
        [
            "filter",
            "--workdir", str(stage_dir),
            "--inputs", str(inputs_dir),
            "--growth-reaction", "r_4041",
        ]
    )
    assert (stage_dir / "all_genes_fseof_essential_filtered.tsv").is_file()
    assert (stage_dir / "families_essential_filtered.tsv").is_file()
    assert (stage_dir / "knockout" / "family_ko_essential_filtered.tsv").is_file()
