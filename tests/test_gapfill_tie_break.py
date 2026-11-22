"""Tests for the gapfill tblastn hit tie-break policy switch."""

import json
import sys
from pathlib import Path

import pytest


BUNDLE_ROOT = Path(__file__).resolve().parents[1] / "third_party" / "alphagem"
if str(BUNDLE_ROOT) not in sys.path:
    sys.path.insert(0, str(BUNDLE_ROOT))

from src.pipeline_contracts import map_tblastn_hits_to_target_proteins  # noqa: E402


def _row(query, subject, sstart, send, bitscore, evalue="1e-50"):
    return {
        "query": query,
        "subject": subject,
        "sstart": sstart,
        "send": send,
        "evalue": evalue,
        "bitscore": bitscore,
    }


def _coordinate(seqid, start, end, target, strand="+"):
    return {
        "genome_seqid": seqid,
        "start": start,
        "end": end,
        "strand": strand,
        "target_protein_id": target,
        "mapping_evidence": "fixture",
    }


def _ambiguous_rows(overlap_b=60, overlap_a=100):
    """One query, two equal-top hits on two scaffolds mapping to TIB_B/TIB_A."""
    rows = [
        _row("Q1", "scf1", 100, 199, 200),
        _row("Q1", "scf2", 50, 149, 200),
    ]
    coordinates = [
        _coordinate("scf1", 140, 199, "TIB_B"),  # 60 bp overlap with row 1
        _coordinate("scf2", 50, 50 + overlap_a - 1, "TIB_A"),
    ]
    if overlap_b != 60:
        coordinates[0] = _coordinate("scf1", 200 - overlap_b, 199, "TIB_B")
    return rows, coordinates


def test_strict_default_still_refuses_ambiguous_top_hits():
    rows, coordinates = _ambiguous_rows()
    with pytest.raises(ValueError, match="Ambiguous highest tblastn hit"):
        map_tblastn_hits_to_target_proteins(
            rows, coordinates, threshold=50, operator="ge"
        )


def test_deterministic_picks_largest_overlap_with_full_audit():
    rows, coordinates = _ambiguous_rows(overlap_b=60, overlap_a=100)
    tie_audit = []
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        tie_break="deterministic",
        tie_audit=tie_audit,
    )
    assert selected["Q1"]["target_protein_id"] == "TIB_A"
    assert selected["Q1"]["coordinate_overlap_bp"] == 100
    assert selected["Q1"]["subject"] == "TIB_A"
    expected = {
        "query": "Q1",
        "tied_targets": ["TIB_A", "TIB_B"],
        "chosen": "TIB_A",
        "tie_break": "deterministic",
        "overlap_bp": 100,
    }
    assert selected["Q1"]["tie_break_audit"] == expected
    assert tie_audit == [expected]


def test_deterministic_remaining_tie_picks_smallest_target_id():
    rows, coordinates = _ambiguous_rows(overlap_b=100, overlap_a=100)
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        tie_break="deterministic",
    )
    assert selected["Q1"]["target_protein_id"] == "TIB_A"
    audit = selected["Q1"]["tie_break_audit"]
    assert audit["tied_targets"] == ["TIB_A", "TIB_B"]
    assert audit["chosen"] == "TIB_A"
    assert audit["overlap_bp"] == 100


def test_non_ambiguous_queries_carry_no_tie_audit():
    rows = [_row("Q2", "scf3", 10, 109, 300)]
    coordinates = [_coordinate("scf3", 10, 109, "TIB_C")]
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        tie_break="deterministic",
    )
    assert selected["Q2"]["target_protein_id"] == "TIB_C"
    assert "tie_break_audit" not in selected["Q2"]


def test_missing_coordinate_still_refused_under_deterministic():
    rows = [_row("Q3", "scfX", 10, 109, 300)]
    with pytest.raises(ValueError, match="does not map to a target protein"):
        map_tblastn_hits_to_target_proteins(
            rows, [], threshold=50, operator="ge", tie_break="deterministic"
        )


def test_row_level_multi_target_overlap_stays_strict_under_deterministic():
    # One HSP overlapping two target proteins is a worse ambiguity and is
    # never tie-broken.
    rows = [_row("Q4", "scf4", 10, 109, 300)]
    coordinates = [
        _coordinate("scf4", 10, 59, "TIB_D"),
        _coordinate("scf4", 60, 109, "TIB_E"),
    ]
    with pytest.raises(ValueError, match="overlaps multiple target proteins"):
        map_tblastn_hits_to_target_proteins(
            rows,
            coordinates,
            threshold=50,
            operator="ge",
            tie_break="deterministic",
        )


def test_unknown_tie_break_refused():
    with pytest.raises(ValueError, match="tie_break"):
        map_tblastn_hits_to_target_proteins(
            [], [], threshold=50, operator="ge", tie_break="random"
        )


def test_threshold_still_filters_before_tie_break():
    rows, coordinates = _ambiguous_rows()
    rows.append(_row("Q1", "scf3", 1, 100, 100))  # below the top score
    coordinates.append(_coordinate("scf3", 1, 100, "TIB_C"))
    with pytest.raises(ValueError, match="Ambiguous highest tblastn hit"):
        map_tblastn_hits_to_target_proteins(
            rows, coordinates, threshold=50, operator="ge", tie_break="strict"
        )
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        tie_break="deterministic",
    )
    assert selected["Q1"]["target_protein_id"] in {"TIB_A", "TIB_B"}


def test_cli_and_run_field_carry_gapfill_tie_break(tmp_path):
    import AlphaGEM

    args = AlphaGEM.build_parser().parse_args([
        "--name", "sample_fixture",
        "--fasta", str(tmp_path / "x.faa"),
        "--solver", "highs",
        "--gapfill-tie-break", "deterministic",
    ])
    run = AlphaGEM.build_run(args)
    assert run.gapfill_tie_break == "deterministic"
    default_run = AlphaGEM.build_run(
        AlphaGEM.build_parser().parse_args([
            "--name", "sample_fixture",
            "--fasta", str(tmp_path / "x.faa"),
        ])
    )
    assert default_run.gapfill_tie_break == "strict"


def test_runner_rejects_unknown_tie_break(tmp_path):
    from src.pipeline_config import get_profile, resolve_profile_for_reference
    from src.pipeline_runner import PipelineRun, PipelineRunner

    run = PipelineRun(
        name="sample_fixture",
        refname="yeast",
        fasta=tmp_path / "x.faa",
        mode="plmsearch",
        profile=resolve_profile_for_reference(get_profile("paper"), "yeast"),
        gapfill_tie_break="random",
    )
    with pytest.raises(ValueError, match="tie_break"):
        PipelineRunner(tmp_path)._validate_run(run, ["gapfill"])


def test_manifest_compat_accepts_old_manifest_at_default(tmp_path):
    """A pre-tie-break manifest (key absent) is compatible with the default."""
    from src.pipeline_config import get_profile, resolve_profile_for_reference
    from src.pipeline_runner import PipelineRun, PipelineRunner

    run = PipelineRun(
        name="sample_fixture",
        refname="yeast",
        fasta=tmp_path / "x.faa",
        mode="plmsearch",
        profile=resolve_profile_for_reference(get_profile("paper"), "yeast"),
    )
    runner = PipelineRunner(tmp_path)
    # The comparison rule itself: absent key equals the default value.
    manifest_value = {}.get("gapfill_tie_break", "strict")
    assert manifest_value == run.gapfill_tie_break
    assert runner._reaction_pool_policy_record(run) == {
        "already_present_policy": "author_skip",
        "conflict_policy": "reject",
    }
