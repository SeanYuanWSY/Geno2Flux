"""Tests for paralog grouping and the unified pipeline glue (manifest,
config resolution, CLI wiring) using synthetic fixtures."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from geno2flux.cli import main as cli_main
from geno2flux.paralog_grouping import (
    TIERS,
    build_families,
    load_self_hits,
    protein_key,
    rank_paralog_families,
)


def _hit(q, s, pident, cov=1.0):
    return {
        "qseqid": q,
        "sseqid": s,
        "pident": pident,
        "length": 100,
        "evalue": 1e-50,
        "bitscore": 500.0,
        "qlen": 100,
        "slen": 100,
        "cov": cov,
    }


def test_protein_key_strips_sp_prefix():
    assert protein_key("sp|ANV000001|ANV000001") == "ANV000001"
    assert protein_key("ANV000001") == "ANV000001"


def test_build_families_directional_edge():
    hits = [
        _hit("A", "B", 90),
        _hit("B", "A", 90),
        _hit("A", "C", 60, 0.5),  # below default
        _hit("D", "D", 100),      # self hit never reaches here (filtered upstream)
    ]
    fams = build_families(hits, TIERS["default"])
    assert sorted(map(sorted, fams)) == [["A", "B"]]


def test_build_families_connected_component_transitive():
    hits = [
        _hit("A", "B", 90), _hit("B", "A", 90),
        _hit("B", "C", 90), _hit("C", "B", 90),
    ]
    fams = build_families(hits, TIERS["default"])
    assert sorted(map(sorted, fams)) == [["A", "B", "C"]]


def test_strict_tier_splits_family():
    hits = [
        _hit("A", "B", 75, 0.85),
        _hit("B", "A", 75, 0.85),
    ]
    assert build_families(hits, TIERS["default"]) == [["A", "B"]]
    assert build_families(hits, TIERS["strict"]) == []


def test_rank_paralog_families_includes_singletons():
    scores = {"A": 2.0, "B": 1.0, "C": None, "D": 2.0, "SOLO": 5.0}
    fams = [["A", "B"], ["C"], ["D"]]
    rows = rank_paralog_families(
        scores, fams, all_protein_ids=["A", "B", "C", "D", "SOLO"]
    )
    fam_solo = next(r for r in rows if r["members"] == ["SOLO"])
    assert fam_solo["size"] == 1
    assert fam_solo["max_score"] == 5.0
    assert fam_solo["max_score_rank"] == 1
    fam_ab = next(r for r in rows if set(r["members"]) == {"A", "B"})
    fam_d = next(r for r in rows if r["members"] == ["D"])
    assert fam_ab["max_score"] == 2.0
    assert fam_ab["mean_score"] == 1.5
    assert fam_ab["max_score_rank"] == 2  # SOLO (5.0) outranks
    assert fam_d["max_score_rank"] == 2


def test_load_self_hits(tmp_path):
    tsv = tmp_path / "hits.tsv"
    tsv.write_text(
        "sp|ANV1|ANV1\tsp|ANV2|ANV2\t90.0\t90\t1e-50\t500\t100\t100\n"
        "sp|ANV1|ANV1\tsp|ANV1|ANV1\t100\t100\t0\t1000\t100\t100\n"  # self, dropped
        "sp|A|A\tsp|B|B\t90\t50\t1e-1\t300\t100\t100\n"               # evalue>1e-10
    )
    hits = load_self_hits(tsv)
    assert len(hits) == 1
    assert hits[0]["qseqid"] == "ANV1"
    assert hits[0]["cov"] == 0.9  # length 90 / max(100,100)


def test_cli_has_pipeline_subcommand(capsys):
    with pytest.raises(SystemExit) as exc:
        cli_main()
    assert exc.value.code != 0
    err = capsys.readouterr().err
    assert "pipeline" in err

def test_pipeline_result_json_serializes_paths(tmp_path):
    """Regression: run_pipeline result dict must survive save_json even when a
    Path object leaks into it."""
    import json as _json
    from geno2flux.common import save_json

    results = {"paths": {"pre-gapfill": {"final_model": Path("/x/m.xml")}},
               "manifest": Path("/x/pipeline_manifest.json")}
    out = tmp_path / "pipeline_result.json"
    save_json(out, results)
    parsed = json.loads(out.read_text())
    assert parsed["manifest"] == "/x/pipeline_manifest.json"
    assert parsed["paths"]["pre-gapfill"]["final_model"] == "/x/m.xml"

def test_family_ids_stable_under_input_shuffle():
    """Shuffling hit order must not change FAM ids or membership."""
    hits = [
        _hit("P1", "P2", 90), _hit("P3", "P4", 95), _hit("P2", "P5", 90),
        _hit("Z1", "Z2", 99), _hit("A1", "A2", 85),
    ]
    fams_a = build_families(hits, TIERS["default"])
    rows_a = rank_paralog_families({}, fams_a, all_protein_ids=["S1", "S2"])
    shuffled = list(reversed(hits))
    fams_b = build_families(shuffled, TIERS["default"])
    rows_b = rank_paralog_families({}, fams_b, all_protein_ids=["S2", "S1"])
    assert fams_a == fams_b
    assert [r["family_id"] for r in rows_a] == [r["family_id"] for r in rows_b]
    assert [tuple(r["members"]) for r in rows_a] == [
        tuple(r["members"]) for r in rows_b
    ]
