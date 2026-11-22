"""End-to-end tests for DeepLoc slot/auto-fallback and genome-to-proteome entry."""
from __future__ import annotations

import json
import os
import stat
from pathlib import Path
from unittest.mock import patch

import pytest

from geno2flux.localization import resolve_deeploc_command
from geno2flux.pipeline import PipelineError, run_pipeline
import scripts.fetch_reference_proteomes as fetch_ref
import scripts.predict_proteome as predict_proteome
import scripts.run_stage as run_stage


def _write_mock_deeploc(bin_path: Path) -> Path:
    """Create an executable mock deeploc2 script that emits DeepLoc CSV output."""
    bin_path.parent.mkdir(parents=True, exist_ok=True)
    bin_path.write_text(
        "#!/usr/bin/env python3\n"
        "import argparse, csv\n"
        "from pathlib import Path\n"
        "p = argparse.ArgumentParser()\n"
        "p.add_argument('-f', required=True)\n"
        "p.add_argument('-o', required=True)\n"
        "a = p.parse_args()\n"
        "out = Path(a.o)\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "with (out / 'results_mock.csv').open('w', newline='', encoding='utf-8') as fh:\n"
        "    w = csv.writer(fh)\n"
        "    w = csv.writer(fh)\n"
        "    w.writerow(['Protein_ID', 'Localizations'])\n"
        "    w.writerow(['G_MDH1', 'Cytoplasm|Mitochondrion'])\n"
        "    w.writerow(['G_MDH2', 'Mitochondrion'])\n",
        encoding="utf-8",
    )
    bin_path.chmod(bin_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return bin_path


def test_resolve_deeploc_command_priority(tmp_path: Path):
    fake_repo = tmp_path / "repo"
    slot_bin = _write_mock_deeploc(fake_repo / "tools" / "deeploc" / "bin" / "deeploc2")
    env_bin = _write_mock_deeploc(tmp_path / "env_bin" / "deeploc2")
    cfg_bin = _write_mock_deeploc(tmp_path / "cfg_bin" / "deeploc2")

    # 1. Explicit configured_cmd wins over everything
    cmd, src = resolve_deeploc_command(
        [str(cfg_bin)],
        env={"G2F_DEEPLOC_BIN": str(env_bin), "PATH": ""},
        repo_root=fake_repo,
    )
    assert cmd == [str(cfg_bin.resolve())]
    assert src == "configured_cmd"

    # 2. G2F_DEEPLOC_BIN wins over repo slot
    cmd, src = resolve_deeploc_command(
        None,
        env={"G2F_DEEPLOC_BIN": str(env_bin), "PATH": ""},
        repo_root=fake_repo,
    )
    assert cmd == [str(env_bin.resolve())]
    assert src == "env_G2F_DEEPLOC_BIN"

    # 3. Repo slot tools/deeploc/bin/deeploc2 auto-detected when env/PATH empty
    cmd, src = resolve_deeploc_command(
        None,
        env={"PATH": ""},
        repo_root=fake_repo,
    )
    assert cmd == [str(slot_bin.resolve())]
    assert src == "repo_tools_deeploc_slot"

    # 4. Returns (None, 'not_found') cleanly when nothing is installed
    empty_repo = tmp_path / "empty_repo"
    empty_repo.mkdir()
    cmd, src = resolve_deeploc_command(
        None,
        env={"PATH": ""},
        repo_root=empty_repo,
    )
    assert cmd is None
    assert src == "not_found"


def test_pipeline_auto_fallback_when_deeploc_absent(
    tmp_path: Path, synthetic_screening_case: Path, monkeypatch: pytest.MonkeyPatch
):
    """When neither localization_tsv nor DeepLoc is available, pipeline falls back
    gracefully to default compartments and completes FSEOF + paralogs + KO."""
    monkeypatch.delenv("G2F_DEEPLOC_BIN", raising=False)
    cfg = json.loads((synthetic_screening_case / "config.json").read_text(encoding="utf-8"))
    cfg.pop("localization_tsv", None)
    cfg.pop("localization_sha256", None)
    cfg.pop("target_module", None)
    cfg["target_reaction"] = "EX_lac_e"
    cfg_path = synthetic_screening_case / "no_deeploc_config.json"
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")

    out_dir = tmp_path / "out_fallback"
    with patch("geno2flux.pipeline.resolve_deeploc_command", return_value=(None, "not_found")):
        result = run_pipeline(cfg_path, out_dir)

    summary = result["paths"]["post-gapfill"]
    assert summary["pulled"] == 0
    audit = json.loads(
        (out_dir / "post-gapfill" / "localization_refinement_audit.json").read_text(encoding="utf-8")
    )
    assert audit["deeploc_status"] == "fallback_default_compartments"
    assert audit["fallback"] == "deeploc_unavailable_default_compartments"

    manifest = json.loads((out_dir / "pipeline_manifest.json").read_text(encoding="utf-8"))
    step4 = next(s for s in manifest["steps"] if s["step"] == "step4_localization[post-gapfill]")
    assert step4["detail"]["deeploc_status"] == "fallback_default_compartments"
    assert step4["inputs"]["localization_tsv"] is None

    # Enforce strict failure when require_deeploc=true
    cfg["require_deeploc"] = True
    strict_cfg = synthetic_screening_case / "strict_deeploc_config.json"
    strict_cfg.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    with patch("geno2flux.pipeline.resolve_deeploc_command", return_value=(None, "not_found")):
        with pytest.raises(PipelineError, match="require_deeploc=true"):
            run_pipeline(strict_cfg, tmp_path / "out_strict_fail")


def test_pipeline_auto_detects_deeploc_via_env_or_slot(
    tmp_path: Path, synthetic_screening_case: Path, monkeypatch: pytest.MonkeyPatch
):
    """Dropping a DeepLoc binary into G2F_DEEPLOC_BIN (or tools/deeploc/) runs
    compartment refinement automatically without localization_tsv in config."""
    mock_dl = _write_mock_deeploc(tmp_path / "tools" / "deeploc" / "bin" / "deeploc2")
    monkeypatch.setenv("G2F_DEEPLOC_BIN", str(mock_dl))

    cfg = json.loads((synthetic_screening_case / "config.json").read_text(encoding="utf-8"))
    cfg.pop("localization_tsv", None)
    cfg.pop("localization_sha256", None)
    cfg_path = synthetic_screening_case / "auto_deeploc_config.json"
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")

    out_dir = tmp_path / "out_auto_dl"
    result = run_pipeline(cfg_path, out_dir)
    summary = result["paths"]["post-gapfill"]
    assert summary["pulled"] == 1
    audit = json.loads(
        (out_dir / "post-gapfill" / "localization_refinement_audit.json").read_text(encoding="utf-8")
    )
    assert audit["deeploc_status"] == "predicted_with_deeploc"
    assert audit["deeploc_source"] == "env_G2F_DEEPLOC_BIN"
    assert audit["pulled_reactions"][0]["reaction_id"] == "MDH_c"


def test_run_stage_localize_fallback_when_deeploc_absent(
    tmp_path: Path, synthetic_screening_case: Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.delenv("G2F_DEEPLOC_BIN", raising=False)
    stage_dir = tmp_path / "stage_loc_fallback"
    with patch("geno2flux.localization.resolve_deeploc_command", return_value=(None, "not_found")):
        run_stage.main(
            [
                "localize",
                "--workdir", str(stage_dir),
                "--inputs", str(synthetic_screening_case),
                "--draft-model", str(synthetic_screening_case / "base_model.xml"),
                "--template-model", str(synthetic_screening_case / "template_model.xml"),
            ]
        )
    assert (stage_dir / "localized.yml").is_file()
    audit = json.loads((stage_dir / "localization_audit.json").read_text(encoding="utf-8"))
    assert audit["deeploc_status"] == "fallback_default_compartments"
    assert audit["pulled_reaction_count"] == 0


def test_predict_proteome_builtin_splicer_without_gffread(
    tmp_path: Path, synthetic_screening_case: Path
):
    """predict_proteome.py standardizes genome + GFF3 -> proteome + AlphaGEM contracts
    using the built-in Python splicer even when gffread is absent."""
    out_dir = tmp_path / "pred_gff"
    rc = predict_proteome.main(
        [
            "--tool", "gff",
            "--genome", str(synthetic_screening_case / "genome.fasta"),
            "--gff", str(synthetic_screening_case / "annotation.gff3"),
            "--case-id", "strain_from_genome",
            "--keep-original-ids",
            "--gffread", "builtin",
            "--output", str(out_dir),
        ]
    )
    assert rc == 0
    receipt = json.loads((out_dir / "predict_proteome_receipt.json").read_text(encoding="utf-8"))
    assert receipt["status"] == "completed"
    assert receipt["target_proteins"] == 8
    assert receipt["standardization"]["gffread_path"] == "builtin_python"
    assert (out_dir / "inputs" / "proteome.fasta").read_bytes() == (
        synthetic_screening_case / "proteome.fasta"
    ).read_bytes()
    assert (out_dir / "inputs" / "target-input-contract.json").is_file()
    assert (out_dir / "inputs" / "gapfill-coordinate-contract.json").is_file()
    assert (out_dir / "03_annotation_qc" / "annotation_qc.json").is_file()


def test_predict_proteome_custom_predictor_command(
    tmp_path: Path, synthetic_screening_case: Path
):
    """predict_proteome.py --tool custom invokes an external gene caller and builds contracts."""
    mock_caller = tmp_path / "mock_caller.py"
    mock_caller.write_text(
        "#!/usr/bin/env python3\n"
        "import shutil, sys\n"
        f"shutil.copy2({str(synthetic_screening_case / 'annotation.gff3')!r}, sys.argv[2])\n",
        encoding="utf-8",
    )
    mock_caller.chmod(0o755)

    out_dir = tmp_path / "pred_custom"
    rc = predict_proteome.main(
        [
            "--tool", "custom",
            "--predictor-cmd", f"{os.fspath(mock_caller)} {{genome}} {{gff}}",
            "--genome", str(synthetic_screening_case / "genome.fasta"),
            "--case-id", "custom_pred",
            "--prefix", "GENE",
            "--gffread", "builtin",
            "--output", str(out_dir),
        ]
    )
    assert rc == 0
    receipt = json.loads((out_dir / "predict_proteome_receipt.json").read_text(encoding="utf-8"))
    assert receipt["target_proteins"] == 8
    assert (out_dir / "inputs" / "proteome.fasta").read_text(encoding="utf-8").count(">GENE") == 8


def test_fetch_reference_proteomes_requires_explicit_sources(tmp_path):
    output = tmp_path / "reference_proteomes"
    with pytest.raises(SystemExit) as exc:
        fetch_ref.main(["--output", str(output)])
    assert exc.value.code == 2
    assert not output.exists()


def test_fetch_reference_proteomes_ingest_helper():
    records: list = []
    mapping: list[dict[str, str]] = []
    fasta_bytes = b">sp|P001|PYC1_YEAST Pyruvate carboxylase\nMKTVL*\n"
    count, sha = fetch_ref._ingest_fasta_bytes(
        fasta_bytes=fasta_bytes,
        token="YEAST",
        source_id="UP000002311",
        source_file_name="UP000002311.faa",
        records=records,
        mapping=mapping,
    )
    assert count == 1
    assert len(sha) == 64
    assert records[0].id == "YEAST|sp|P001|PYC1_YEAST"
    assert str(records[0].seq) == "MKTVL"
    assert mapping[0]["source_accession"] == "UP000002311"
