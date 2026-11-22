"""Growth-essentiality filtering for target prioritization rankings.

ESSENTIALITY SOURCE — growth, not product. The pipeline's own
``single_gene_ko`` tables optimize the user-selected PRODUCT objective;
using their loss as "essential" would delete the wrong genes. Essential
genes therefore come from a GROWTH-objective single-gene KO table
(biomass objective, same medium, same solver; e.g.
``single_gene_ko_growth.tsv``), supplied by config:

- ``essential_genes_tsv``: path to that table (columns ``gene_id``,
  ``loss``, ``status``; status ``optimal`` rows with
  ``loss >= essential_ko_loss_threshold`` are essential, boundary
  INCLUSIVE — 0.5 exactly counts). A bare id list (header without
  ``loss``/``status``) is accepted as an explicit essential set.
- ``essential_ko_loss_threshold``: default 0.5 (at least half of the
  base biomass flux lost). Config-validated to [0, 1].

Derived rules: essential FAMILY = any member essential. Filtered
outputs (FSEOF gene ranking, family ranking, single-gene KO, family KO)
drop essential entries, keep every original column and rank verbatim,
and append a contiguous ``filtered_rank`` over the survivors. Reaction
tables are untouched: their ``single_reaction_ko_filtered`` is a PATHWAY
filter (boundary / multi-compartment / respiratory), a different axis —
see README. No gene, reaction or family list is hardcoded anywhere.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable


def validate_threshold(value: Any) -> float:
    """Threshold must be a finite fraction in [0, 1]; never silently coerced."""
    try:
        threshold = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"essential_ko_loss_threshold must be a number, got {value!r}"
        ) from exc
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(
            f"essential_ko_loss_threshold must lie in [0, 1], got {threshold}"
        )
    return threshold


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, list(reader)


def load_essential_gene_ids(path: Path, threshold: float) -> dict[str, Any]:
    """Essential gene ids from a growth-KO table or a bare id list.

    Accepts exactly two shapes and refuses everything else loudly:
    - a growth-KO table with ALL of ``gene_id``/``loss``/``status``
      (threshold derivation, optimal rows, inclusive boundary);
    - a SINGLE-column table (an explicit essential id list; the column
      name is irrelevant, threshold bypassed).
    A multi-column table missing any KO column is a CONFIG ERROR — a
    half-KO table must never be silently read as an all-rows id list
    (that would drop every gene from the filtered rankings).
    """
    fields, rows = read_tsv_rows(path)
    has_ko_columns = {"gene_id", "loss", "status"} <= set(fields)
    if has_ko_columns:
        essential: set[str] = set()
        for row in rows:
            if row.get("status") != "optimal":
                continue
            try:
                loss = float(row.get("loss"))
            except (TypeError, ValueError):
                continue
            if math.isfinite(loss) and loss >= threshold:
                essential.add(row["gene_id"])
        return {
            "ids": essential,
            "mode": "growth_ko_table_threshold",
            "source": str(path),
            "threshold": threshold,
        }
    if len(fields) == 1:
        id_col = fields[0]
        return {
            "ids": {value for row in rows if (value := (row.get(id_col) or "").strip())},
            "mode": "explicit_id_list",
            "source": str(path),
            "threshold": None,
        }
    raise ValueError(
        f"essential_genes_tsv is neither a growth-KO table "
        f"(needs gene_id/loss/status columns; got {fields}) nor a "
        f"single-column id list - refusing to guess"
    )


def essential_family_ids(
    families: Iterable[dict[str, Any]], essential_genes: set[str]
) -> set[str]:
    """A family is essential when ANY member gene is essential."""
    essential: set[str] = set()
    for fam in families:
        members = fam.get("members", [])
        if isinstance(members, str):
            members = [m for m in members.split(";") if m]
        if any(m in essential_genes for m in members):
            essential.add(fam["family_id"])
    return essential


def _numeric_rank(value: Any) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _rank_key(row: dict[str, Any], primary_rank: str, id_key: str):
    numeric = _numeric_rank(row.get(primary_rank))
    return (
        numeric is None,  # ranked rows first, unranked last
        numeric if numeric is not None else 0.0,
        str(row.get(id_key, "")),
    )


def _gene_key(row: dict[str, Any], id_key: str):
    """Order robust candidates first, raw-only candidates next, unscored last."""
    robust = _numeric_rank(row.get("robust_rank"))
    if robust is not None:
        return (0, robust, str(row.get(id_key, "")))
    raw = _numeric_rank(row.get("raw_rank"))
    return (
        1 if raw is not None else 2,
        raw if raw is not None else 0.0,
        str(row.get(id_key, "")),
    )


def _has_rank(row: dict[str, Any], rank_key: str) -> bool:
    return _numeric_rank(row.get(rank_key)) is not None


def apply_essential_filter(
    rows: list[dict[str, Any]],
    essential: set[str],
    id_key: str,
    primary_rank: str,
    *,
    order: str = "primary",
) -> dict[str, Any]:
    """Split rows into kept/dropped; kept get a dense filtered_rank.

    Original rank columns preserved verbatim. ``order="primary"``:
    rows sort by the primary rank ascending (numeric-aware; read-back
    TSV rows carry string ranks) with the id as tie-break; rows without
    a primary rank sort last. ``order="gene_v2"``: pre-registered v2
    gene order (robust_rank primary, empty-robust tail by raw_rank).

    A row whose SOURCE rank is empty gets an EMPTY filtered_rank — a
    non-candidate must never look like a bottom-ranked candidate.
    Every kept row also gets ``filtered_row_index`` (its 1..n position
    in the output file): that is a ROW INDEX for orientation, NOT a
    rank, and must never be quoted as one. Dropped rows are returned
    unmodified.
    """
    kept = [dict(r) for r in rows if r.get(id_key) not in essential]
    dropped = [r for r in rows if r.get(id_key) in essential]
    if order == "gene_v2":
        kept.sort(key=lambda r: _gene_key(r, id_key))
    else:
        kept.sort(key=lambda r: _rank_key(r, primary_rank, id_key))
    next_rank = 1
    for index, row in enumerate(kept, start=1):
        row["filtered_row_index"] = index  # position, NOT a rank
        if _has_rank(row, "robust_rank" if order == "gene_v2"
                     else primary_rank):
            row["filtered_rank"] = next_rank
            next_rank += 1
        else:
            row["filtered_rank"] = ""
    return {"kept": kept, "dropped": dropped}


def positional_filtered_rank(
    kept: list[dict[str, Any]],
    primary_rank: str,
    id_key: str,
    rank_name: str,
) -> None:
    """Attach a POSITIONAL contiguous rank under ``rank_name``.

    Dense 1..m positions of a re-sort by ``primary_rank`` ascending
    (id tie-break); rows whose ``primary_rank`` is EMPTY get an EMPTY
    value — each column is judged independently against its own source
    rank (a family unscored on mean must not borrow a mean rank).
    """
    ordered = sorted(
        kept, key=lambda r: _rank_key(r, primary_rank, id_key)
    )
    next_rank = 1
    for row in ordered:
        if _has_rank(row, primary_rank):
            row[rank_name] = next_rank
            next_rank += 1
        else:
            row[rank_name] = ""

def write_filtered_tsv(
    path: Path, rows: list[dict[str, Any]], fieldnames: list[str]
) -> None:
    """Filtered rows with filtered_rank (rank; EMPTY for unranked) and
    filtered_row_index (row position, NOT a rank) as the LAST columns.

    List values (e.g. in-memory ``members``) are ';'-joined to match
    the source table format, and appended columns are skipped when the
    caller already included them (header stays duplicate-free).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    appended = ("filtered_rank", "filtered_row_index")
    fields = list(dict.fromkeys(col for col in fieldnames if col not in appended))
    # Retain additional source columns, including columns unknown to the
    # canonical schemas used by the pipeline's in-memory callers.
    for row in rows:
        for col in row:
            if col not in fields and col not in appended:
                fields.append(col)
    fields.extend(appended)
    normalized = [
        {
            k: (";".join(v) if isinstance(v, list) else v)
            for k, v in row.items()
        }
        for row in rows
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(normalized)


def write_filter_manifest_json(path: Path, detail: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(detail, indent=2, allow_nan=False, sort_keys=True) + "\n"
    )


FSEOF_GENE_FIELDS = [
    "gene_id",
    "in_model",
    "raw_score",
    "robust_score",
    "raw_candidate_reactions",
    "robust_candidate_reactions",
    "all_associated_reactions",
    "raw_rank",
    "robust_rank",
]

FAMILY_FIELDS = [
    "family_id",
    "size",
    "n_scored",
    "max_score",
    "mean_score",
    "max_score_rank",
    "mean_score_rank",
    "members",
]

KO_GENE_FIELDS = ["ko_rank", "gene_id", "gene_name", "obj", "loss", "status"]

KO_FAMILY_FIELDS = [
    "ko_rank",
    "family_id",
    "members",
    "size",
    "obj",
    "loss",
    "status",
]
