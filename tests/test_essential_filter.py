"""Essentiality-filtered ranking tests (synthetic tables, real logic).

Covers the growth-vs-product rule, the inclusive 0.5 boundary, optimal-
row restriction, member->family propagation, original-rank preservation
and contiguous filtered_rank. No mock echoes: every assertion checks
computed outputs of the real functions on synthetic input tables.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from geno2flux.essential_filter import (
    apply_essential_filter,
    essential_family_ids,
    load_essential_gene_ids,
    validate_threshold,
    write_filtered_tsv,
)


def _write_tsv(path: Path, header: list[str], rows: list[list]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _growth_table(tmp: Path) -> Path:
    # boundary probe: EXACTLY 0.5 is essential; 0.4999 is not;
    # non-optimal rows never are; unparseable loss never is.
    return _write_tsv(
        tmp / "single_gene_ko_growth.tsv",
        ["ko_rank", "gene_id", "gene_name", "obj", "loss", "status"],
        [
            [1, "G_BOUNDARY", "b", 0.5, "0.5", "optimal"],        # in
            [2, "G_JUST_BELOW", "j", 0.5001, "0.4999", "optimal"],  # out
            [3, "G_HIGH", "h", 0.1, "0.9", "optimal"],            # in
            [4, "G_NONOPT", "n", "", "1.0", "infeasible"],        # out
            [5, "G_BADLOSS", "x", 0.7, "not-a-number", "optimal"],  # out
        ],
    )


def test_growth_table_threshold_boundary_inclusive(tmp_path: Path):
    record = load_essential_gene_ids(_growth_table(tmp_path), 0.5)
    assert record["mode"] == "growth_ko_table_threshold"
    assert record["ids"] == {"G_BOUNDARY", "G_HIGH"}



def test_half_ko_table_missing_status_is_rejected(tmp_path: Path):
    """gene_id+loss WITHOUT status is a config error, never an id list
    (reading it as a list would mark EVERY row essential)."""
    bad = _write_tsv(
        tmp_path / "half_ko.tsv",
        ["ko_rank", "gene_id", "loss"],
        [[1, "G1", "0.9"], [2, "G2", "0.1"]],
    )
    with pytest.raises(ValueError, match="refusing to guess"):
        load_essential_gene_ids(bad, 0.5)

def test_explicit_id_list_bypasses_threshold(tmp_path: Path):
    table = _write_tsv(
        tmp_path / "essential_ids.tsv", ["gene_id"], [["G_A"], ["G_B"]]
    )
    record = load_essential_gene_ids(table, 0.5)
    assert record["mode"] == "explicit_id_list"
    assert record["ids"] == {"G_A", "G_B"}
    assert record["threshold"] is None


def test_threshold_validation_rejects_out_of_range():
    for bad in (-0.01, 1.01, "zero", None):
        with pytest.raises(ValueError):
            validate_threshold(bad)
    assert validate_threshold(0.5) == 0.5
    assert validate_threshold("0.25") == 0.25


def test_family_essential_any_member(tmp_path: Path):
    essential = load_essential_gene_ids(_growth_table(tmp_path), 0.5)["ids"]
    families = [
        {"family_id": "FAM00001", "members": ["G_OK1", "G_OK2"]},
        {"family_id": "FAM00002", "members": ["G_OK1", "G_HIGH"]},
        {"family_id": "FAM00003", "members": ["G_BOUNDARY"]},
        # read-back shape: members as ';'-joined string
        {"family_id": "FAM00004", "members": "G_OK2;G_OK3"},
    ]
    assert essential_family_ids(families, essential) == {
        "FAM00002",
        "FAM00003",
    }


def test_filter_preserves_ranks_and_contiguous_filtered_rank():
    essential = {"G_DROP"}  # e.g. raw_rank 2
    rows = [
        {"gene_id": "A", "raw_rank": 1, "robust_rank": 2},
        {"gene_id": "G_DROP", "raw_rank": 2, "robust_rank": 1},
        {"gene_id": "C", "raw_rank": 3, "robust_rank": 3},
    ]
    result = apply_essential_filter(rows, essential, "gene_id", "raw_rank")
    kept = result["kept"]
    assert [r["gene_id"] for r in kept] == ["A", "C"]
    assert [r["filtered_rank"] for r in kept] == [1, 2]
    # original ranks verbatim, untouched by the filter
    assert kept[0]["raw_rank"] == 1 and kept[1]["raw_rank"] == 3
    assert [r["gene_id"] for r in result["dropped"]] == ["G_DROP"]
    assert "filtered_rank" not in result["dropped"][0]


def test_string_ranks_sort_numerically():
    # read-back TSV rows carry string ranks: 2 < 10 must hold numerically
    rows = [
        {"gene_id": "B", "ko_rank": "10"},
        {"gene_id": "A", "ko_rank": "2"},
        {"gene_id": "C", "ko_rank": "1"},
    ]
    kept = apply_essential_filter(rows, set(), "gene_id", "ko_rank")["kept"]
    assert [r["gene_id"] for r in kept] == ["C", "A", "B"]
    assert [r["filtered_rank"] for r in kept] == [1, 2, 3]


def test_unranked_rows_sort_last_and_get_blank_rank():
    """Empty source rank => EMPTY filtered_rank (never a fake bottom
    rank); filtered_row_index still orients the row; order unchanged."""
    rows = [
        {"gene_id": "Z", "raw_rank": None},
        {"gene_id": "A", "raw_rank": None},
        {"gene_id": "M", "raw_rank": 1},
    ]
    kept = apply_essential_filter(rows, set(), "gene_id", "raw_rank")["kept"]
    assert [r["gene_id"] for r in kept] == ["M", "A", "Z"]
    assert [r["filtered_rank"] for r in kept] == [1, "", ""]
    assert [r["filtered_row_index"] for r in kept] == [1, 2, 3]
    # no numeric leakage into unranked rows
    assert all(r["filtered_rank"] == "" for r in kept if r["raw_rank"] is None)


def test_gene_v2_order_robust_primary_raw_tail():
    """v2 gene order: robust_rank ascending, empty-robust tail by
    raw_rank. Empty robust_rank => EMPTY filtered_rank (no fake bottom
    rank); row order and filtered_row_index still cover every row."""
    rows = [
        {"gene_id": "R2", "raw_rank": 2, "robust_rank": 2},
        {"gene_id": "R1", "raw_rank": 1, "robust_rank": 3},
        {"gene_id": "T1", "raw_rank": 5, "robust_rank": None},
        {"gene_id": "T2", "raw_rank": 4, "robust_rank": None},
    ]
    kept = apply_essential_filter(
        rows, set(), "gene_id", "raw_rank", order="gene_v2"
    )["kept"]
    assert [r["gene_id"] for r in kept] == ["R2", "R1", "T2", "T1"]
    assert [r["filtered_rank"] for r in kept] == [1, 2, "", ""]
    assert [r["filtered_row_index"] for r in kept] == [1, 2, 3, 4]


def test_family_positional_columns_match_pk_ground_truth():
    """filtered_rank_max/mean are DENSE 1..n positions (pK formal
    comparison ground truth): re-sort by the respective source rank
    ascending (id tie-break), rows lacking that rank in the tail."""
    from geno2flux.essential_filter import positional_filtered_rank

    families = [
        {"family_id": "F1", "max_score": 3.0, "mean_score": 2.0,
         "max_score_rank": 1, "mean_score_rank": 2},
        {"family_id": "F2", "max_score": 3.0, "mean_score": 1.0,
         "max_score_rank": 1, "mean_score_rank": 3},   # ties F1 on max
        {"family_id": "F3", "max_score": 1.0, "mean_score": 3.0,
         "max_score_rank": 3, "mean_score_rank": 1},   # essential, dropped
        {"family_id": "F4", "max_score": None, "mean_score": None,
         "max_score_rank": None, "mean_score_rank": None},
        {"family_id": "F5", "max_score": 2.0, "mean_score": 2.0,
         "max_score_rank": 2, "mean_score_rank": 2},
    ]
    result = apply_essential_filter(
        families, {"F3"}, "family_id", "max_score_rank"
    )
    kept = result["kept"]
    positional_filtered_rank(kept, "max_score_rank", "family_id",
                             "filtered_rank_max")
    positional_filtered_rank(kept, "mean_score_rank", "family_id",
                             "filtered_rank_mean")
    by_id = {r["family_id"]: r for r in kept}
    # max order: F1,F2 (rank 1, id tie-break) -> 1,2; F5 (2) -> 3;
    # F4 unscored -> BLANK (never a fake bottom rank)
    assert by_id["F1"]["filtered_rank_max"] == 1
    assert by_id["F2"]["filtered_rank_max"] == 2
    assert by_id["F5"]["filtered_rank_max"] == 3
    assert by_id["F4"]["filtered_rank_max"] == ""
    # mean judged independently against mean_score_rank: F4 blank there
    # too; scored families (PYC-style) keep dense numeric ranks
    assert by_id["F1"]["filtered_rank_mean"] == 1
    assert by_id["F5"]["filtered_rank_mean"] == 2
    assert by_id["F2"]["filtered_rank_mean"] == 3
    assert by_id["F4"]["filtered_rank_mean"] == ""
    # contiguous positional filtered_rank (by max) dense over ranked only
    assert [r["filtered_rank"] for r in kept] == [1, 2, 3, ""]


def test_filtered_tsv_appends_rank_column_last(tmp_path: Path):
    rows = [{"gene_id": "A", "raw_rank": 1, "filtered_rank": 1,
             "filtered_row_index": 1}]
    out = tmp_path / "f.tsv"
    write_filtered_tsv(out, rows, ["gene_id", "raw_rank"])
    text = out.read_text()
    assert (
        text.splitlines()[0]
        == "gene_id\traw_rank\tfiltered_rank\tfiltered_row_index"
    )
    assert text.splitlines()[1] == "A\t1\t1\t1"


def test_step9b_wiring_real_dirs_real_manifest(tmp_path: Path):
    """The pipeline step-9b helper end-to-end on real synthetic TSVs.

    Real files in, real Manifest, four filtered TSVs + detail JSON +
    manifest entry out; also proves the not-configured no-op summary.
    """
    from geno2flux import pipeline as P

    out = tmp_path / "out"
    fseof_dir = out / "pre-gapfill" / "fseof"
    para_dir = out / "pre-gapfill" / "paralog"
    ko_dir = out / "pre-gapfill" / "knockout"

    gene_rows = [
        {"gene_id": "G1", "in_model": True, "raw_score": 3.0,
         "robust_score": 3.0, "raw_candidate_reactions": "R1",
         "robust_candidate_reactions": "R1",
         "all_associated_reactions": "R1", "raw_rank": 1,
         "robust_rank": 1},
        {"gene_id": "G2", "in_model": True, "raw_score": 2.0,
         "robust_score": 2.0, "raw_candidate_reactions": "R2",
         "robust_candidate_reactions": "R2",
         "all_associated_reactions": "R2", "raw_rank": 2,
         "robust_rank": 2},
    ]
    write_filtered_tsv(fseof_dir / "all_genes_fseof.tsv",
                       [dict(r) for r in gene_rows],
                       list(gene_rows[0].keys()))
    family_rows = [
        {"family_id": "FAM00001", "size": 1, "n_scored": 1,
         "max_score": 3.0, "mean_score": 3.0, "max_score_rank": 1,
         "mean_score_rank": 1, "members": ["G1"]},
        {"family_id": "FAM00002", "size": 1, "n_scored": 1,
         "max_score": 2.0, "mean_score": 2.0, "max_score_rank": 2,
         "mean_score_rank": 2, "members": ["G2"]},
    ]
    write_filtered_tsv(
        para_dir / "paralog_family_rankings.tsv",
        [dict(r, members=";".join(r["members"])) for r in family_rows],
        ["family_id", "size", "n_scored", "max_score", "mean_score",
         "max_score_rank", "mean_score_rank", "members"],
    )
    _write_tsv(
        ko_dir / "single_gene_ko.tsv",
        ["ko_rank", "gene_id", "gene_name", "obj", "loss", "status"],
        [[1, "G1", "a", 0.2, 0.8, "optimal"],
         [2, "G2", "b", 0.9, 0.1, "optimal"]],
    )
    _write_tsv(
        ko_dir / "family_ko.tsv",
        ["ko_rank", "family_id", "members", "size", "obj", "loss",
         "status"],
        [[1, "FAM00001", "G1", 1, 0.2, 0.8, "optimal"],
         [2, "FAM00002", "G2", 1, 0.9, 0.1, "optimal"]],
    )
    growth = _write_tsv(
        tmp_path / "single_gene_ko_growth.tsv",
        ["ko_rank", "gene_id", "gene_name", "obj", "loss", "status"],
        [[1, "G1", "a", 0.4, 0.6, "optimal"],   # essential (>= 0.5)
        [2, "G2", "b", 0.9, 0.1, "optimal"]],
    )

    manifest = P.Manifest(out)
    summary = P._write_essential_filtered_outputs(
        fseof_dir=fseof_dir, para_dir=para_dir, ko_dir=ko_dir,
        gene_rows=gene_rows, family_rows=family_rows,
        essential_tsv=growth, cfg={"essential_ko_loss_threshold": 0.5},
        manifest=manifest, od="pre-gapfill",
    )
    assert summary["essential_genes"] == 1
    assert summary["essential_families"] == 1  # FAM00001 via member G1
    assert summary["fseof_genes_kept"] == 1
    assert summary["gene_ko_kept"] == 1

    gene_f = (fseof_dir / "all_genes_fseof_essential_filtered.tsv").read_text()
    lines = gene_f.splitlines()
    assert lines[0].endswith("filtered_row_index")
    assert lines[1].startswith("G2") and lines[1].endswith("1\t1")
    assert "G1" not in gene_f
    fam_ko = (ko_dir / "family_ko_essential_filtered.tsv").read_text()
    assert "FAM00001" not in fam_ko and "FAM00002" in fam_ko
    # family ranking filtered: three rank columns, tie-aware pair present
    fam_f = (para_dir
             / "paralog_family_rankings_essential_filtered.tsv").read_text()
    fam_header = fam_f.splitlines()[0].split("\t")
    assert fam_header[-4:] == [
        "filtered_rank_max", "filtered_rank_mean",
        "filtered_rank", "filtered_row_index",
    ]
    # members serialized like the source table (';'-joined), never list text
    assert "['" not in fam_f
    kept_line = next(ln for ln in fam_f.splitlines() if "FAM00002" in ln)
    assert kept_line.split("\t")[fam_header.index("members")] == "G2"
    # header has no duplicate columns
    assert len(fam_header) == len(set(fam_header))
    assert "FAM00001" not in fam_f
    # growth table filtered: essential row dropped, survivor keeps ko_rank
    growth_f = (ko_dir
                / "single_gene_ko_growth_essential_filtered.tsv").read_text()
    assert "G1\t" not in growth_f and "\tG2\t" in growth_f

    detail = json.loads((ko_dir / "essential_filter_detail.json").read_text())
    assert detail["essentiality_source"] == "growth_ko_table_threshold"
    assert detail["growth_ko_kept"] == 1
    m = json.loads((out / "pipeline_manifest.json").read_text())
    entry = next(s for s in m["steps"]
                 if s["step"] == "step9b_essential_filter[pre-gapfill]")
    assert set(entry["outputs"]) == {
        "fseof_gene_filtered", "family_ranking_filtered",
        "single_gene_ko_filtered", "family_ko_filtered",
        "growth_ko_filtered", "detail",
    }

    # not configured -> no-op summary, no outputs
    out2 = tmp_path / "out2"
    s2 = P._write_essential_filtered_outputs(
        fseof_dir=out2, para_dir=out2, ko_dir=out2, gene_rows=[],
        family_rows=[], essential_tsv=None, cfg={}, manifest=None,
        od="pre-gapfill",
    )
    assert s2 == {"essentiality_source": "not_configured"}


def test_product_ko_table_is_not_an_essentiality_source(tmp_path: Path):
    """The rule's sharpest edge: a gene with PRODUCT-loss 1.0 but growth
    loss below threshold stays in the filtered gene-KO table."""
    essential = load_essential_gene_ids(_growth_table(tmp_path), 0.5)["ids"]
    product_ko_rows = [
        # G_PRODUCTLETHAL is absent from the growth table entirely
        {"gene_id": "G_PRODUCTLETHAL", "ko_rank": 1, "loss": "1.0"},
        {"gene_id": "G_HIGH", "ko_rank": 2, "loss": "0.1"},
    ]
    kept = apply_essential_filter(
        product_ko_rows, essential, "gene_id", "ko_rank"
    )["kept"]
    assert [r["gene_id"] for r in kept] == ["G_PRODUCTLETHAL"]


def test_gene_v2_tail_keeps_raw_candidates_before_unscored_rows():
    rows = [
        {"gene_id": "UNSCORED", "raw_rank": None, "robust_rank": None},
        {"gene_id": "RAW_2", "raw_rank": "2", "robust_rank": ""},
        {"gene_id": "ROBUST", "raw_rank": "3", "robust_rank": "1"},
        {"gene_id": "RAW_1", "raw_rank": "1", "robust_rank": ""},
    ]
    kept = apply_essential_filter(
        rows, set(), "gene_id", "raw_rank", order="gene_v2"
    )["kept"]
    assert [row["gene_id"] for row in kept] == [
        "ROBUST", "RAW_1", "RAW_2", "UNSCORED"
    ]
    assert [row["filtered_rank"] for row in kept] == [1, "", "", ""]


@pytest.mark.parametrize("invalid_rank", ["nan", "inf", "-inf", float("nan")])
def test_nonfinite_source_ranks_never_become_candidates(invalid_rank):
    from geno2flux.essential_filter import positional_filtered_rank

    rows = [
        {"gene_id": "BAD", "ko_rank": invalid_rank},
        {"gene_id": "GOOD", "ko_rank": "1"},
    ]
    kept = apply_essential_filter(rows, set(), "gene_id", "ko_rank")["kept"]
    positional_filtered_rank(kept, "ko_rank", "gene_id", "alternate_rank")
    assert [row["gene_id"] for row in kept] == ["GOOD", "BAD"]
    assert [row["filtered_rank"] for row in kept] == [1, ""]
    assert [row["alternate_rank"] for row in kept] == [1, ""]


def test_growth_table_ignores_nonfinite_losses(tmp_path):
    path = _write_tsv(
        tmp_path / "growth.tsv", ["gene_id", "loss", "status"],
        [["NAN", "nan", "optimal"], ["INF", "inf", "optimal"],
         ["ESSENTIAL", "0.5", "optimal"]],
    )
    assert load_essential_gene_ids(path, 0.5)["ids"] == {"ESSENTIAL"}


def test_filtered_writer_preserves_additional_source_columns(tmp_path):
    from geno2flux.essential_filter import read_tsv_rows

    rows = [{"gene_id": "A", "ko_rank": 1, "source_receipt": "run-1",
             "filtered_rank": 1, "filtered_row_index": 1}]
    path = tmp_path / "filtered.tsv"
    write_filtered_tsv(path, rows, ["gene_id", "ko_rank"])
    fields, saved = read_tsv_rows(path)
    assert fields == ["gene_id", "ko_rank", "source_receipt",
                      "filtered_rank", "filtered_row_index"]
    assert saved[0]["source_receipt"] == "run-1"


def test_paralog_hit_keys_match_fseof_and_essential_gene_keys(tmp_path):
    from geno2flux.paralog_grouping import (
        TIERS, build_families, load_self_hits, protein_key, rank_paralog_families,
    )

    hits_path = tmp_path / "hits.tsv"
    hits_path.write_text(
        "sp|P1|ENTRY_A\tsp|P2|ENTRY_B\t90\t90\t1e-50\t500\t100\t100\n"
    )
    hits = load_self_hits(hits_path)
    scores = {protein_key("sp|P1|ENTRY_A"): 3.0,
              protein_key("sp|P2|ENTRY_B"): 1.0}
    families = rank_paralog_families(
        scores, build_families(hits, TIERS["default"]), all_protein_ids=scores,
    )
    assert len(families) == 1
    assert families[0]["members"] == ["P1", "P2"]
    assert families[0]["n_scored"] == 2
    assert families[0]["mean_score"] == 2.0
    assert essential_family_ids(families, {"P1"}) == {families[0]["family_id"]}


@pytest.fixture
def recompute_script():
    import runpy

    return runpy.run_path(str(
        Path(__file__).resolve().parents[1]
        / "scripts" / "recompute_essential_filtered.py"
    ))


def _comparison_tables(root, targets, *, empty=False):
    for relative, id_col, rank_cols in targets:
        fields = [id_col, "source_score", *rank_cols]
        rows = [] if empty else [["A", "3.0", *([1] * len(rank_cols))]]
        _write_tsv(root / relative, fields, rows)


def test_recompute_comparison_accepts_empty_survivor_tables(tmp_path, recompute_script):
    mine, theirs = tmp_path / "mine", tmp_path / "theirs"
    for root in (mine, theirs):
        _comparison_tables(root, recompute_script["_COMPARE_TARGETS"], empty=True)
    assert recompute_script["_compare"](mine, theirs) == 0


@pytest.mark.parametrize("side", ["mine", "theirs"])
def test_recompute_comparison_requires_rank_headers_on_both_sides(
    tmp_path, recompute_script, side,
):
    mine, theirs = tmp_path / "mine", tmp_path / "theirs"
    targets = recompute_script["_COMPARE_TARGETS"]
    for root in (mine, theirs):
        _comparison_tables(root, targets, empty=True)
    relative, id_col, _ = targets[0]
    _write_tsv(tmp_path / side / relative, [id_col, "source_score"], [])
    assert recompute_script["_compare"](mine, theirs) == 1


def test_recompute_comparison_detects_changed_source_values(tmp_path, recompute_script):
    mine, theirs = tmp_path / "mine", tmp_path / "theirs"
    targets = recompute_script["_COMPARE_TARGETS"]
    for root in (mine, theirs):
        _comparison_tables(root, targets)
    assert recompute_script["_compare"](mine, theirs) == 0
    relative, id_col, rank_cols = targets[0]
    _write_tsv(theirs / relative, [id_col, "source_score", *rank_cols],
               [["A", "4.0", 1]])
    assert recompute_script["_compare"](mine, theirs) == 1


def test_recompute_preserves_empty_source_table_headers(
    tmp_path, recompute_script, monkeypatch,
):
    from geno2flux.essential_filter import read_tsv_rows

    source, output = tmp_path / "source", tmp_path / "output"
    tables = [
        ("fseof/all_genes_fseof.tsv", ["gene_id", "raw_rank", "robust_rank"]),
        ("paralog/paralog_family_rankings.tsv",
         ["family_id", "members", "max_score_rank", "mean_score_rank"]),
        ("knockout/single_gene_ko.tsv", ["gene_id", "ko_rank"]),
        ("knockout/family_ko.tsv", ["family_id", "ko_rank"]),
    ]
    for relative, fields in tables:
        _write_tsv(source / relative, [*fields, "source_receipt"], [])
    growth = _write_tsv(
        tmp_path / "growth.tsv", ["gene_id", "loss", "status", "source_receipt"], [],
    )
    monkeypatch.setattr("sys.argv", [
        "recompute", "--input-dir", str(source), "--output-dir", str(output),
        "--growth-tsv", str(growth),
    ])
    assert recompute_script["main"]() == 0
    for relative, _, _ in recompute_script["_COMPARE_TARGETS"]:
        fields, rows = read_tsv_rows(output / relative)
        assert "source_receipt" in fields
        assert not rows
