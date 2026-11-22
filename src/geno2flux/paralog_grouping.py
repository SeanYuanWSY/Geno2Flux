"""Paralog family grouping and competition ranking for FSEOF gene scores.

An all-vs-all DIAMOND self-alignment of the proteome yields pairwise hits.
Proteins are merged into families as connected components under three
pre-registered threshold tiers, matching the audited reference implementation
``scripts/a01_paralog_ranks.py``.  Each family (multi-member components plus
every singleton protein in the proteome) receives ``max`` and ``mean``
aggregate FSEOF scores plus competition ranks.
"""
from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping

# Pre-registered tiers: (min_pident_percent, min_coverage_fraction).
TIERS: dict[str, tuple[float, float]] = {
    "default": (70.0, 0.80),
    "strict": (80.0, 0.90),
    "loose": (50.0, 0.70),
}

E_CUT = 1e-10

# DIAMOND blastp outfmt 6 field order used here. ``length`` is the HSP length
# and coverage is recomputed as ``length / max(qlen, slen)`` so no dependence on
# the ``qcovhsp``/``scovhsp`` conventions (which are 0..100 percents).
BLAST_OUTFMT_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "evalue",
    "bitscore", "qlen", "slen",
)


def protein_key(protein_id: str) -> str:
    """Normalize a FASTA/BLAST identifier to its bare accession.

    ``sp|P12345|ENTRY_NAME`` -> ``P12345`` (the second ``|``-field),
    matching AlphaGEM's protein keys and the FSEOF proteome loader.
    Plain IDs pass through.
    """
    value = str(protein_id).split()[0]
    parts = value.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return value




def run_diamond_self_alignment(
    proteome_fasta: Path,
    out_tsv: Path,
    *,
    diamond_bin: str = "diamond",
    threads: int = 8,
    evalue: str = "1e-10",
    max_target_seqs: int = 500,
    workdir: Path | None = None,
) -> Path:
    """Run ``diamond makedb`` + ``diamond blastp`` self-alignment.

    Emits the 8-column layout ``qseqid sseqid pident length evalue bitscore
    qlen slen`` expected by :func:`load_self_hits`.
    """
    wd = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="paralog-"))
    wd.mkdir(parents=True, exist_ok=True)
    db = wd / "self.dmnd"
    subprocess.run(
        [diamond_bin, "makedb", "--in", str(proteome_fasta), "--db", str(db)],
        check=True,
    )
    subprocess.run(
        [
            diamond_bin, "blastp",
            "--query", str(proteome_fasta),
            "--db", str(db),
            "--out", str(out_tsv),
            "--outfmt", "6", *BLAST_OUTFMT_FIELDS,
            "--evalue", evalue,
            "--max-target-seqs", str(max_target_seqs),
            "--threads", str(threads),
        ],
        check=True,
    )
    return out_tsv


def load_self_hits(tsv_path: Path) -> list[dict[str, Any]]:
    """Parse a DIAMOND self-alignment TSV into hit dicts.

    The first eight columns must follow ``BLAST_OUTFMT_FIELDS``; trailing
    columns are ignored. Coverage is computed from ``length/qlen/slen``.
    Rows failing ``qseqid != sseqid`` or ``evalue <= E_CUT`` are dropped.
    """
    hits: list[dict[str, Any]] = []
    with tsv_path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 8:
                continue
            q = protein_key(row[0])
            s = protein_key(row[1])
            pid = float(row[2])
            length = int(float(row[3]))
            ev = float(row[4])
            qlen = int(float(row[6]))
            slen = int(float(row[7]))
            if q == s or ev > E_CUT:
                continue
            cov = length / max(qlen, slen)
            hits.append(
                {
                    "qseqid": q,
                    "sseqid": s,
                    "pident": pid,
                    "length": length,
                    "evalue": ev,
                    "bitscore": float(row[5]),
                    "qlen": qlen,
                    "slen": slen,
                    "cov": cov,
                }
            )
    return hits


def _components(nodes: set[str], edges: list[tuple[str, str]]) -> list[list[str]]:
    """Union-find connected components over undirected edges."""
    parent = {n: n for n in nodes}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    groups: dict[str, list[str]] = {}
    for n in nodes:
        groups.setdefault(find(n), []).append(n)
    return list(groups.values())


def build_families(
    hits: Iterable[Mapping[str, Any]], tier: tuple[float, float] | Mapping[str, float]
) -> list[list[str]]:
    """Multi-member connected components over hits passing the tier.

    Any single-direction hit row with ``pident >= min_pid`` and
    ``cov >= min_cov`` creates an undirected edge (union-find then merges
    transitively), exactly as the reference ``cluster(pairs, min_pid, min_cov)``.
    Returns only families with >= 2 members.
    """
    if isinstance(tier, Mapping):
        min_pid = tier["min_pident"]
        min_cov = tier["min_cov"]
    else:
        min_pid, min_cov = tier
    nodes: set[str] = set()
    edges: list[tuple[str, str]] = []
    for hit in hits:
        if hit["pident"] >= min_pid and hit["cov"] >= min_cov:
            q, s = hit["qseqid"], hit["sseqid"]
            nodes.add(q)
            nodes.add(s)
            edges.append((q, s))
    groups = _components(nodes, edges)
    # Deterministic order: sort each family, then sort families by their
    # member tuple so FAM numbering is stable across runs/input order.
    return sorted((sorted(g) for g in groups if len(g) >= 2))


def rank_paralog_families(
    gene_scores: Mapping[str, float | None],
    families: list[list[str]],
    *,
    all_protein_ids: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    """Aggregate per-family scores and emit competition-ranked rows.

    ``families`` is the multi-member component list from :func:`build_families`.
    Every protein in ``all_protein_ids`` not already inside a multi-member
    family is appended as a singleton ``size=1`` family so the full proteome is
    ranked.  ``gene_scores`` maps protein_id -> raw FSEOF score (or ``None``).
    """
    in_group: set[str] = set()
    for members in families:
        in_group.update(members)
    groups = sorted((sorted(g) for g in families))
    if all_protein_ids is not None:
        singletons = sorted(m for m in set(all_protein_ids) if m not in in_group)
        groups.extend([m] for m in singletons)

    rows: list[dict[str, Any]] = []
    for index, members in enumerate(groups):
        scored = [gene_scores.get(m) for m in members]
        scored = [s for s in scored if s is not None]
        rows.append(
            {
                "family_id": f"FAM{index:05d}",
                "members": members,
                "size": len(members),
                "n_scored": len(scored),
                "max_score": max(scored) if scored else None,
                "mean_score": (sum(scored) / len(scored)) if scored else None,
            }
        )
    for agg in ("max_score", "mean_score"):
        present = [r for r in rows if r[agg] is not None]
        for row in rows:
            if row[agg] is None:
                row[f"{agg}_rank"] = None
            else:
                row[f"{agg}_rank"] = 1 + sum(
                    1 for r in present if r[agg] > row[agg]
                )
    return rows


def write_family_tables(
    families: list[dict[str, Any]], json_path: Path, tsv_path: Path
) -> None:
    """Emit machine-readable JSON and flat TSV family rankings."""
    import json

    json_path.write_text(
        json.dumps({"families": families}, indent=2, allow_nan=False) + "\n"
    )
    with tsv_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "family_id",
                "size",
                "n_scored",
                "max_score",
                "mean_score",
                "max_score_rank",
                "mean_score_rank",
                "members",
            ]
        )
        for row in families:
            writer.writerow(
                [
                    row["family_id"],
                    row["size"],
                    row["n_scored"],
                    row["max_score"],
                    row["mean_score"],
                    row["max_score_rank"],
                    row["mean_score_rank"],
                    ";".join(row["members"]),
                ]
            )
