"""Tests for the generalized AlphaGEM target & coordinate contract builder."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import scripts.build_coordinate_contract as contract_builder
from scripts.build_coordinate_contract import build


def _fixture(tmp_path: Path, *, mismatch: bool = False):
    genome = tmp_path / "genome.fna"
    genome.write_text(
        ">scaffold,plus\nATGAAATAA\n>scaffold,minus\nTTATTTCAT\n",
        encoding="utf-8",
    )
    gff = tmp_path / "prediction.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "scaffold,plus\tANNEVO\tmRNA\t1\t9\t.\t+\t.\tID=tx,plus;Parent=gene,plus\n"
        "scaffold,plus\tANNEVO\tCDS\t1\t9\t.\t+\t0\tID=cds,plus;Parent=tx,plus\n"
        "scaffold,minus\tANNEVO\tmRNA\t1\t9\t.\t-\t.\tID=tx,minus;Parent=gene,minus\n"
        "scaffold,minus\tANNEVO\tCDS\t1\t9\t.\t-\t0\tID=cds,minus;Parent=tx,minus\n",
        encoding="utf-8",
    )
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(
        ">tx,plus\n" + ("MA" if mismatch else "MK") + "\n>tx,minus\nMK\n",
        encoding="utf-8",
    )
    return genome, gff, proteins


def test_build_coordinate_contract_normalizes_comma_ids_and_validates(tmp_path: Path):
    genome, gff, proteins = _fixture(tmp_path)
    out = tmp_path / "contracts"
    receipt = build(
        genome=genome,
        gff=gff,
        proteins=proteins,
        case_id="fungal_sample_01",
        prefix="FUNGI",
        source_label="annevo",
        output=out,
    )
    assert receipt["status"] == "strict_contract_validated"
    assert receipt["target_proteins"] == 2
    assert receipt["coordinate_segment_records"] == 2
    assert (out / "proteome.fasta").read_text().count(">FUNGI") == 2
    assert (out / "genome.fasta").is_file()
    assert json.loads((out / "target-input-contract.json").read_text())["source_fasta"]["record_count"] == 2
    with (out / "coordinate-map.tsv").open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["refseq_protein_id"] for r in rows} == {"tx,plus", "tx,minus"}
    assert {r["strand"] for r in rows} == {"+", "-"}


def test_build_coordinate_contract_allows_representative_isoform_subset(tmp_path: Path):
    """Longest-isoform protein FASTAs legitimately omit alternative transcripts in the GFF3."""
    genome, gff, _ = _fixture(tmp_path)
    subset_faa = tmp_path / "longest_isoform.faa"
    subset_faa.write_text(">tx,plus\nMK\n", encoding="utf-8")
    out = tmp_path / "subset_contracts"
    receipt = build(
        genome=genome,
        gff=gff,
        proteins=subset_faa,
        case_id="isoform_subset",
        output=out,
    )
    assert receipt["target_proteins"] == 1
    assert receipt["non_selected_gff_transcripts"] == 1


def test_build_coordinate_contract_rejects_translation_mismatch_cleanly(tmp_path: Path):
    genome, gff, proteins = _fixture(tmp_path, mismatch=True)
    out = tmp_path / "failed_contracts"
    with pytest.raises(ValueError, match="translation differs"):
        build(genome=genome, gff=gff, proteins=proteins, case_id="bad_case", output=out)
    assert not out.exists()


def test_build_coordinate_contract_cli_keep_original_ids(
    tmp_path: Path, synthetic_screening_case: Path
):
    out = tmp_path / "cli_contracts"
    assert contract_builder.main(
        [
            "--genome", str(synthetic_screening_case / "genome.fasta"),
            "--gff", str(synthetic_screening_case / "annotation.gff3"),
            "--proteins", str(synthetic_screening_case / "proteome.fasta"),
            "--case-id", "cli_demo",
            "--keep-original-ids",
            "--output", str(out),
        ]
    ) == 0
    assert (out / "target-input-contract.json").is_file()
    assert (out / "gapfill-coordinate-contract.json").is_file()
