"""Tests for the gapfill unmapped-hit policy switch."""

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


def test_error_policy_default_keeps_unmapped_hit_fatal():
    rows = [_row("Q1", "scf88", 11785, 12561, 231)]
    coordinates = [_coordinate("scf88", 1, 100, "TIB_A")]
    with pytest.raises(ValueError, match="does not map to a target protein"):
        map_tblastn_hits_to_target_proteins(
            rows, coordinates, threshold=50, operator="ge"
        )


def test_skip_policy_all_unmapped_leaves_query_absent_with_audit():
    rows = [
        _row("Q1", "scf88", 11785, 12561, 231),
        _row("Q1", "scf31", 300, 450, 231),
    ]
    coordinates = [
        _coordinate("scf88", 1, 100, "TIB_A"),
        _coordinate("scf31", 900, 999, "TIB_B"),
    ]
    audit = []
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        unmapped_hit_policy="skip_with_audit",
        unmapped_audit=audit,
    )
    assert "Q1" not in selected
    # Audit rows follow the function's deterministic hit ranking order
    # (bitscore, evalue, subject, coordinates), not the raw file order.
    assert audit[:2] == [
        {
            "query": "Q1",
            "subject": "scf31",
            "sstart": 300,
            "send": 450,
            "strand": "+",
            "bitscore": 231.0,
            "evalue": 1e-50,
            "reason": "no_overlapping_target_coordinate",
        },
        {
            "query": "Q1",
            "subject": "scf88",
            "sstart": 11785,
            "send": 12561,
            "strand": "+",
            "bitscore": 231.0,
            "evalue": 1e-50,
            "reason": "no_overlapping_target_coordinate",
        },
    ]
    # The query-level record makes the fully unmapped outcome explicit.
    assert audit[2]["query"] == "Q1"
    assert audit[2]["reason"] == "all_equal_top_hits_unmapped"


def test_skip_policy_mapped_row_wins_and_skip_is_recorded():
    rows = [
        _row("Q1", "scf88", 11785, 12561, 231),  # unmapped
        _row("Q1", "scf1", 50, 149, 231),        # maps to TIB_A
    ]
    coordinates = [_coordinate("scf1", 50, 149, "TIB_A")]
    audit = []
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        unmapped_hit_policy="skip_with_audit",
        unmapped_audit=audit,
    )
    assert selected["Q1"]["target_protein_id"] == "TIB_A"
    assert selected["Q1"]["subject"] == "TIB_A"
    assert [record["reason"] for record in audit] == [
        "no_overlapping_target_coordinate"
    ]
    assert audit[0]["subject"] == "scf88"
    # Mixed outcome: no query-level all-unmapped record.
    assert not any(
        record["reason"] == "all_equal_top_hits_unmapped" for record in audit
    )


def test_skip_policy_strand_mismatch_is_still_unmapped():
    rows = [_row("Q1", "scf1", 50, 149, 231)]
    # Coordinates exist only on the opposite strand.
    coordinates = [_coordinate("scf1", 50, 149, "TIB_A", strand="-")]
    audit = []
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        unmapped_hit_policy="skip_with_audit",
        unmapped_audit=audit,
    )
    assert "Q1" not in selected
    assert audit[-1]["reason"] == "all_equal_top_hits_unmapped"


def test_single_hsp_multi_target_ambiguity_stays_fatal_under_skip():
    rows = [_row("Q1", "scf1", 10, 109, 300)]
    coordinates = [
        _coordinate("scf1", 10, 59, "TIB_A"),
        _coordinate("scf1", 60, 109, "TIB_B"),
    ]
    with pytest.raises(ValueError, match="overlaps multiple target proteins"):
        map_tblastn_hits_to_target_proteins(
            rows,
            coordinates,
            threshold=50,
            operator="ge",
            unmapped_hit_policy="skip_with_audit",
        )


def test_empty_coordinate_table_skip_vs_error():
    rows = [_row("Q1", "scf1", 10, 109, 300)]
    with pytest.raises(ValueError, match="does not map to a target protein"):
        map_tblastn_hits_to_target_proteins(
            rows, [], threshold=50, operator="ge"
        )
    audit = []
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        [],
        threshold=50,
        operator="ge",
        unmapped_hit_policy="skip_with_audit",
        unmapped_audit=audit,
    )
    assert "Q1" not in selected
    assert audit[-1]["reason"] == "all_equal_top_hits_unmapped"


def test_unknown_unmapped_hit_policy_refused():
    with pytest.raises(ValueError, match="unmapped_hit_policy"):
        map_tblastn_hits_to_target_proteins(
            [], [], threshold=50, operator="ge", unmapped_hit_policy="guess"
        )


def test_skip_composes_with_deterministic_tie_break():
    # One unmapped row plus two equal-top rows mapping to two targets: the
    # skip happens first, then the tie is broken deterministically.
    rows = [
        _row("Q1", "scf88", 11785, 12561, 200),  # unmapped
        _row("Q1", "scf1", 100, 199, 200),
        _row("Q1", "scf2", 50, 149, 200),
    ]
    coordinates = [
        _coordinate("scf1", 140, 199, "TIB_B"),
        _coordinate("scf2", 50, 149, "TIB_A"),
    ]
    unmapped_audit = []
    selected = map_tblastn_hits_to_target_proteins(
        rows,
        coordinates,
        threshold=50,
        operator="ge",
        tie_break="deterministic",
        unmapped_hit_policy="skip_with_audit",
        unmapped_audit=unmapped_audit,
    )
    assert selected["Q1"]["target_protein_id"] == "TIB_A"
    assert selected["Q1"]["tie_break_audit"]["chosen"] == "TIB_A"
    assert [r["reason"] for r in unmapped_audit] == [
        "no_overlapping_target_coordinate"
    ]


def test_cli_and_run_field_carry_unmapped_hit_policy(tmp_path):
    import AlphaGEM

    args = AlphaGEM.build_parser().parse_args([
        "--name", "sample_fixture",
        "--fasta", str(tmp_path / "x.faa"),
        "--solver", "highs",
        "--gapfill-unmapped-hit-policy", "skip_with_audit",
    ])
    run = AlphaGEM.build_run(args)
    assert run.gapfill_unmapped_hit_policy == "skip_with_audit"
    default_run = AlphaGEM.build_run(
        AlphaGEM.build_parser().parse_args([
            "--name", "sample_fixture",
            "--fasta", str(tmp_path / "x.faa"),
        ])
    )
    assert default_run.gapfill_unmapped_hit_policy == "error"


def test_runner_rejects_unknown_unmapped_hit_policy(tmp_path):
    from src.pipeline_config import get_profile, resolve_profile_for_reference
    from src.pipeline_runner import PipelineRun, PipelineRunner

    run = PipelineRun(
        name="sample_fixture",
        refname="yeast",
        fasta=tmp_path / "x.faa",
        mode="plmsearch",
        profile=resolve_profile_for_reference(get_profile("paper"), "yeast"),
        gapfill_unmapped_hit_policy="drop",
    )
    with pytest.raises(ValueError, match="unmapped_hit_policy"):
        PipelineRunner(tmp_path)._validate_run(run, ["gapfill"])
