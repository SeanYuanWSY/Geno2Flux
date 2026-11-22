#!/usr/bin/env python3
"""Recompute essential-filtered target tables from a completed screening output directory.

Usage (from the repository root):
    uv run --frozen --no-sync python scripts/recompute_essential_filtered.py \
        --input-dir <completed_order_dir> \
        --growth-tsv <single_gene_ko_growth.tsv> \
        --output-dir <recompute_dir> \
        [--compare-dir <reference_dir>] [--threshold 0.5]

Applies the exact Step 9b rules:
- essential genes: growth-KO optimal rows with loss >= threshold (inclusive)
- essential families: any member essential
- gene filtered_rank: robust_rank ascending, empty-robust tail by raw_rank
- family: positional filtered_rank + positional filtered_rank_max/mean
- growth table filtered of the essential genes it identified
--compare-dir checks all five filtered outputs by ID and shared source columns,
including the required filtered rank columns and empty survivor tables.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from geno2flux.essential_filter import (  # noqa: E402
    apply_essential_filter,
    essential_family_ids,
    load_essential_gene_ids,
    positional_filtered_rank,
    read_tsv_rows,
    validate_threshold,
    write_filter_manifest_json,
    write_filtered_tsv,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--growth-tsv", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--threshold", default="0.5")
    ap.add_argument("--compare-dir", type=Path, default=None)
    args = ap.parse_args()

    threshold = validate_threshold(args.threshold)
    src = args.input_dir
    record = load_essential_gene_ids(args.growth_tsv, threshold)
    ess_genes = record["ids"]
    if record["mode"] != "growth_ko_table_threshold":
        raise ValueError("--growth-tsv requires a growth-KO table with gene_id/loss/status")

    fam_fields, fam_rank_rows = read_tsv_rows(
        src / "paralog" / "paralog_family_rankings.tsv"
    )
    gene_fields, gene_rows = read_tsv_rows(src / "fseof" / "all_genes_fseof.tsv")

    ess_fams = essential_family_ids(fam_rank_rows, ess_genes)
    out = args.output_dir
    (out / "fseof").mkdir(parents=True, exist_ok=True)
    (out / "paralog").mkdir(parents=True, exist_ok=True)
    (out / "knockout").mkdir(parents=True, exist_ok=True)

    gres = apply_essential_filter(
        gene_rows, ess_genes, "gene_id", "raw_rank", order="gene_v2"
    )
    write_filtered_tsv(
        out / "fseof" / "all_genes_fseof_essential_filtered.tsv",
        gres["kept"], gene_fields,
    )
    fres = apply_essential_filter(
        fam_rank_rows, ess_fams, "family_id", "max_score_rank"
    )
    positional_filtered_rank(
        fres["kept"], "max_score_rank", "family_id", "filtered_rank_max"
    )
    positional_filtered_rank(
        fres["kept"], "mean_score_rank", "family_id", "filtered_rank_mean"
    )
    write_filtered_tsv(
        out / "paralog" / "paralog_family_rankings_essential_filtered.tsv",
        fres["kept"], fam_fields + ["filtered_rank_max", "filtered_rank_mean"],
    )
    growth_fields, growth_rows = read_tsv_rows(args.growth_tsv)
    gg = apply_essential_filter(
        [r for r in growth_rows if r.get("status") == "optimal"],
        ess_genes, "gene_id", "ko_rank",
    )
    write_filtered_tsv(
        out / "knockout" / "single_gene_ko_growth_essential_filtered.tsv",
        gg["kept"], growth_fields,
    )

    gene_ko_fields, gene_ko_rows = read_tsv_rows(
        src / "knockout" / "single_gene_ko.tsv"
    )
    kg = apply_essential_filter(gene_ko_rows, ess_genes, "gene_id", "ko_rank")
    write_filtered_tsv(
        out / "knockout" / "single_gene_ko_essential_filtered.tsv",
        kg["kept"], gene_ko_fields,
    )
    fam_ko_fields, fam_ko_rows = read_tsv_rows(src / "knockout" / "family_ko.tsv")
    kf = apply_essential_filter(fam_ko_rows, ess_fams, "family_id", "ko_rank")
    write_filtered_tsv(
        out / "knockout" / "family_ko_essential_filtered.tsv",
        kf["kept"], fam_ko_fields,
    )
    write_filter_manifest_json(
        out / "essential_filter_detail.json",
        {
            "essentiality_source": record["mode"],
            "essentiality_input": str(args.growth_tsv),
            "threshold": record["threshold"],
            "essential_genes": len(ess_genes),
            "essential_families": len(ess_fams),
            "fseof_genes_kept": len(gres["kept"]),
            "families_kept": len(fres["kept"]),
            "gene_ko_kept": len(kg["kept"]),
            "family_ko_kept": len(kf["kept"]),
        },
    )
    print(
        f"recomputed: essential_genes={len(ess_genes)} "
        f"essential_families={len(ess_fams)} "
        f"fseof_kept={len(gres['kept'])} families_kept={len(fres['kept'])} "
        f"growth_ko_kept={len(gg['kept'])}"
    )

    if args.compare_dir:
        return _compare(out, args.compare_dir)
    return 0


_COMPARE_TARGETS = [
    (
        "fseof/all_genes_fseof_essential_filtered.tsv",
        "gene_id",
        ["filtered_rank"],
    ),
    (
        "paralog/paralog_family_rankings_essential_filtered.tsv",
        "family_id",
        ["filtered_rank_max", "filtered_rank_mean"],
    ),
    (
        "knockout/single_gene_ko_growth_essential_filtered.tsv",
        "gene_id",
        ["filtered_rank"],
    ),
    (
        "knockout/single_gene_ko_essential_filtered.tsv",
        "gene_id",
        ["filtered_rank"],
    ),
    (
        "knockout/family_ko_essential_filtered.tsv",
        "family_id",
        ["filtered_rank"],
    ),
]


def _compare(mine: Path, theirs: Path) -> int:
    """Compare ID sets, ID uniqueness, and every shared column across both directories."""
    failures = 0
    for rel, id_col, rank_cols in _COMPARE_TARGETS:
        a_path, b_path = mine / rel, theirs / rel
        if not b_path.exists():
            print(f"MISSING in compare dir: {rel}")
            failures += 1
            continue
        a_fields, a_rows = read_tsv_rows(a_path)
        b_fields, b_rows = read_tsv_rows(b_path)
        missing_columns = False
        for label, fields in (("mine", a_fields), ("theirs", b_fields)):
            missing = sorted({id_col, *rank_cols} - set(fields))
            if missing:
                print(f"DIFF {rel} ({label}): missing required columns {missing}")
                failures += 1
                missing_columns = True
        if missing_columns:
            continue
        for label, rows in (("mine", a_rows), ("theirs", b_rows)):
            ids = [r[id_col] for r in rows]
            if len(ids) != len(set(ids)):
                dupes = sorted({i for i in ids if ids.count(i) > 1})
                print(f"DIFF {rel} ({label}): duplicate {id_col}: {dupes[:5]}")
                failures += 1
        a_map = {r[id_col]: r for r in a_rows}
        b_map = {r[id_col]: r for r in b_rows}
        only_a = sorted(set(a_map) - set(b_map))
        only_b = sorted(set(b_map) - set(a_map))
        if only_a or only_b or len(a_rows) != len(b_rows):
            print(
                f"DIFF {rel}: only_mine={only_a[:5]} "
                f"only_theirs={only_b[:5]} rows {len(a_rows)} vs {len(b_rows)}"
            )
            failures += 1
            continue
        shared_cols = sorted((set(a_fields) & set(b_fields)) - {id_col})
        row_diffs = 0
        for key in sorted(a_map):
            for col in shared_cols:
                va = (a_map[key].get(col) or "").strip()
                vb = (b_map[key].get(col) or "").strip()
                if va != vb and not (va in ("", "None") and vb in ("", "None")):
                    print(f"DIFF {rel} {key}: {col} {va!r} vs {vb!r}")
                    row_diffs += 1
        if row_diffs:
            failures += 1
        else:
            print(
                f"OK {rel}: {len(a_rows)} rows identical on "
                f"{len(shared_cols)} shared columns (incl. {rank_cols})"
            )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
