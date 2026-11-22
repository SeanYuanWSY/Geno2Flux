#!/usr/bin/env python3
"""Create deterministic annotation and proteome QC tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")


@dataclass(frozen=True)
class Feature:
    seqid: str
    kind: str
    start: int
    end: int
    strand: str
    attrs: dict[str, tuple[str, ...]]

    @property
    def identifier(self) -> str | None:
        values = self.attrs.get("ID", ())
        return values[0] if len(values) == 1 else None

    @property
    def parents(self) -> tuple[str, ...]:
        return self.attrs.get("Parent", ())


def parse_attrs(raw: str) -> dict[str, tuple[str, ...]]:
    if raw in {"", "."}:
        return {}
    parts = raw.split(";")
    if parts[-1] == "":
        parts = parts[:-1]
    attrs: dict[str, tuple[str, ...]] = {}
    for part in parts:
        if not part or "=" not in part:
            raise ValueError(f"malformed GFF3 attributes: {raw}")
        key, value = part.split("=", 1)
        key = unquote(key)
        if key in attrs:
            raise ValueError(f"duplicate GFF3 attribute: {key}")
        attrs[key] = tuple(unquote(item) for item in value.split(","))
    return attrs


def read_gff(path: Path) -> list[Feature]:
    features: list[Feature] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 columns")
            features.append(
                Feature(
                    seqid=fields[0],
                    kind=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    attrs=parse_attrs(fields[8]),
                )
            )
    return features


def read_fasta(path: Path) -> tuple[list[tuple[str, str]], Counter[str]]:
    records: list[tuple[str, str]] = []
    counts: Counter[str] = Counter()
    current_id: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks).upper()))
                current_id = line[1:].split()[0]
                counts[current_id] += 1
                chunks = []
            elif current_id is None:
                raise ValueError(f"sequence before FASTA header in {path}")
            else:
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, "".join(chunks).upper()))
    return records, counts


def read_mapping(
    path: Path,
) -> tuple[set[str], set[str], dict[str, str], dict[str, str]]:
    representative_proteins: set[str] = set()
    transcript_to_gene: dict[str, str] = {}
    transcript_to_protein: dict[str, str] = {}
    representative_counts: Counter[str] = Counter()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            transcript_id = row["transcript_id"]
            if transcript_id in transcript_to_gene:
                raise ValueError(f"duplicate transcript in mapping: {transcript_id}")
            protein_id = row["protein_id"]
            if protein_id in transcript_to_protein.values():
                raise ValueError(f"duplicate protein in mapping: {protein_id}")
            transcript_to_gene[transcript_id] = row["gene_id"]
            transcript_to_protein[transcript_id] = protein_id
            if row["representative"].lower() == "true":
                representative_proteins.add(protein_id)
                representative_counts[row["gene_id"]] += 1
    invalid_genes = sorted(
        gene_id for gene_id in set(transcript_to_gene.values())
        if representative_counts[gene_id] != 1
    )
    if invalid_genes:
        raise ValueError(
            "each gene must have exactly one representative: "
            f"{invalid_genes[:5]}"
        )
    representative_transcripts = {
        transcript_id
        for transcript_id, protein_id in transcript_to_protein.items()
        if protein_id in representative_proteins
    }
    return (
        representative_transcripts,
        representative_proteins,
        transcript_to_gene,
        transcript_to_protein,
    )


def read_genome_lengths(path: Path) -> dict[str, int]:
    records, counts = read_fasta(path)
    duplicate_ids = [identifier for identifier, count in counts.items() if count != 1]
    if duplicate_ids:
        raise ValueError(f"duplicate genome IDs: {duplicate_ids[:5]}")
    return {identifier: len(sequence) for identifier, sequence in records}


def read_bed(path: Path | None) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if path is None:
        return intervals
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected BED3")
            start_zero = int(fields[1])
            end_zero = int(fields[2])
            if start_zero < 0 or end_zero <= start_zero:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            intervals[fields[0]].append((start_zero + 1, end_zero))
    for seqid in intervals:
        intervals[seqid].sort()
        merged: list[tuple[int, int]] = []
        for left, right in intervals[seqid]:
            if merged and left <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], right))
            else:
                merged.append((left, right))
        intervals[seqid] = merged
    return intervals


def unique_features(features: list[Feature], kind: str) -> dict[str, Feature]:
    result: dict[str, Feature] = {}
    for feature in features:
        if feature.kind != kind:
            continue
        identifier = feature.identifier
        if identifier is None:
            raise ValueError(f"{kind} must have exactly one ID")
        if identifier in result:
            raise ValueError(f"duplicate {kind} ID: {identifier}")
        result[identifier] = feature
    return result


def excluded_counts(path: Path) -> tuple[int, int]:
    reason_rows = 0
    model_keys: set[tuple[str, str]] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2):
            reason_rows += 1
            key = next(
                (
                    (field, row.get(field, ""))
                    for field in ("transcript_id", "gene_id", "protein_id")
                    if row.get(field, "")
                ),
                ("line", str(line_number)),
            )
            model_keys.add(key)
    return len(model_keys), reason_rows


def protein_length_and_internal_stop(sequence: str) -> tuple[int, bool]:
    translated = sequence[:-1] if sequence.endswith("*") else sequence
    return len(translated), "*" in translated


def protein_rows(
    method: str, records: list[tuple[str, str]]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for identifier, sequence in records:
        protein_length, internal_stop = protein_length_and_internal_stop(sequence)
        rows.append(
            {
                "method": method,
                "protein_id": identifier,
                "length": protein_length,
                "internal_stop": str(internal_stop).lower(),
                "invalid_amino_acid": str(bool(set(sequence) - VALID_AA)).lower(),
            }
        )
    return rows


def overlap_bases(start: int, end: int, intervals: list[tuple[int, int]]) -> int:
    total = 0
    for left, right in intervals:
        if right < start:
            continue
        if left > end:
            break
        total += max(0, min(end, right) - max(start, left) + 1)
    return total


def percentile(values: list[int], fraction: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def n50(values: list[int]) -> int:
    threshold = sum(values) / 2
    cumulative = 0
    for value in sorted(values, reverse=True):
        cumulative += value
        if cumulative >= threshold:
            return value
    return 0


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_qc(
    *,
    method: str,
    genome_path: Path,
    gff_path: Path,
    all_proteins_path: Path,
    representative_proteins_path: Path,
    mapping_path: Path,
    excluded_path: Path,
    output_dir: Path,
    repeat_bed: Path | None = None,
    edge_bp: int = 500,
) -> dict[str, object]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    genome_lengths = read_genome_lengths(genome_path)
    features = read_gff(gff_path)
    (
        representatives,
        representative_protein_ids,
        transcript_to_gene,
        transcript_to_protein,
    ) = read_mapping(mapping_path)
    repeats = read_bed(repeat_bed)

    genes = unique_features(features, "gene")
    transcript_features = [
        feature for feature in features if feature.kind in {"mRNA", "transcript"}
    ]
    transcripts: dict[str, Feature] = {}
    for transcript in transcript_features:
        transcript_id = transcript.identifier
        if transcript_id is None:
            raise ValueError("transcript must have exactly one ID")
        if transcript_id in transcripts:
            raise ValueError(f"duplicate transcript ID: {transcript_id}")
        transcripts[transcript_id] = transcript
    if set(transcript_to_gene) != set(transcripts):
        raise ValueError("mapping transcript IDs do not equal GFF3 transcript IDs")
    if set(transcript_to_gene.values()) != set(genes):
        raise ValueError("mapping gene IDs do not equal GFF3 gene IDs")
    for transcript_id, gene_id in transcript_to_gene.items():
        if transcripts[transcript_id].parents != (gene_id,):
            raise ValueError(f"mapping/GFF3 parent mismatch: {transcript_id}")
    cds_by_transcript: dict[str, list[Feature]] = defaultdict(list)
    exon_by_transcript: dict[str, list[Feature]] = defaultdict(list)
    feature_counts = Counter(feature.kind for feature in features)
    for feature in features:
        if feature.kind == "CDS":
            for parent in feature.parents:
                cds_by_transcript[parent].append(feature)
        elif feature.kind == "exon":
            for parent in feature.parents:
                exon_by_transcript[parent].append(feature)

    structure_rows: list[dict[str, object]] = []
    for transcript_id in sorted(representatives):
        transcript = transcripts[transcript_id]
        gene_id = transcript_to_gene[transcript_id]
        gene = genes[gene_id]
        cds = sorted(cds_by_transcript[transcript_id], key=lambda item: item.start)
        exons = sorted(exon_by_transcript.get(transcript_id, cds), key=lambda item: item.start)
        introns = [max(0, right.start - left.end - 1) for left, right in zip(exons, exons[1:])]
        gene_span = gene.end - gene.start + 1
        repeat_bp = overlap_bases(gene.start, gene.end, repeats.get(gene.seqid, []))
        scaffold_length = genome_lengths[gene.seqid]
        structure_rows.append(
            {
                "method": method,
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "scaffold": gene.seqid,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "gene_span": gene_span,
                "transcript_span": transcript.end - transcript.start + 1,
                "cds_length": sum(item.end - item.start + 1 for item in cds),
                "cds_segment_count": len(cds),
                "exon_count": len(exons),
                "intron_count": len(introns),
                "intron_total_length": sum(introns),
                "intron_mean_length": statistics.mean(introns) if introns else 0,
                "single_exon": str(len(exons) == 1).lower(),
                "scaffold_edge": str(
                    gene.start <= edge_bp
                    or scaffold_length - gene.end < edge_bp
                ).lower(),
                "repeat_overlap_bp": repeat_bp,
                "repeat_overlap_fraction": repeat_bp / gene_span,
            }
        )

    structure_fields = list(structure_rows[0]) if structure_rows else []
    write_tsv(output_dir / "gene_model_structure.tsv", structure_rows, structure_fields)

    all_records, all_counts = read_fasta(all_proteins_path)
    representative_records, representative_counts = read_fasta(representative_proteins_path)
    all_protein_ids = set(all_counts)
    expected_all_proteins = set(transcript_to_protein.values())
    if all_protein_ids != expected_all_proteins:
        raise ValueError("all-isoform FASTA IDs do not equal mapping protein IDs")
    if set(representative_counts) != representative_protein_ids:
        raise ValueError("representative FASTA IDs do not equal mapping representatives")
    length_rows = protein_rows(method, representative_records)
    write_tsv(
        output_dir / "protein_lengths.tsv",
        length_rows,
        ["method", "protein_id", "length", "internal_stop", "invalid_amino_acid"],
    )
    lengths = [int(row["length"]) for row in length_rows]
    excluded_count, excluded_reason_rows = excluded_counts(excluded_path)
    duplicate_all = sum(count - 1 for count in all_counts.values() if count > 1)
    duplicate_representative = sum(
        count - 1 for count in representative_counts.values() if count > 1
    )
    summary = {
        "method": method,
        "genes": len(genes),
        "transcripts": len(transcripts),
        "cds_features": feature_counts["CDS"],
        "all_isoform_proteins": len(all_records),
        "representative_proteins": len(representative_records),
        "unique_all_isoform_protein_ids": len(all_counts),
        "duplicate_all_isoform_protein_records": duplicate_all,
        "duplicate_representative_protein_records": duplicate_representative,
        "empty_representative_proteins": sum(not sequence for _, sequence in representative_records),
        "invalid_amino_acid_representative_proteins": sum(
            row["invalid_amino_acid"] == "true" for row in length_rows
        ),
        "internal_stop_representative_proteins": sum(
            row["internal_stop"] == "true" for row in length_rows
        ),
        "excluded_models": excluded_count,
        "excluded_reason_rows": excluded_reason_rows,
        "total_amino_acids": sum(lengths),
        "minimum_protein_length": min(lengths, default=0),
        "maximum_protein_length": max(lengths, default=0),
        "mean_protein_length": statistics.mean(lengths) if lengths else math.nan,
        "median_protein_length": statistics.median(lengths) if lengths else math.nan,
        "q1_protein_length": percentile(lengths, 0.25),
        "q3_protein_length": percentile(lengths, 0.75),
        "protein_n50": n50(lengths),
        "proteins_lt_50": sum(value < 50 for value in lengths),
        "proteins_lt_100": sum(value < 100 for value in lengths),
        "proteins_gt_1000": sum(value > 1000 for value in lengths),
        "proteins_gt_2000": sum(value > 2000 for value in lengths),
        "single_exon_genes": sum(row["single_exon"] == "true" for row in structure_rows),
        "scaffold_edge_genes": sum(row["scaffold_edge"] == "true" for row in structure_rows),
        "repeat_overlapping_genes": sum(
            float(row["repeat_overlap_fraction"]) > 0 for row in structure_rows
        ),
        "plus_strand_genes": sum(row["strand"] == "+" for row in structure_rows),
        "minus_strand_genes": sum(row["strand"] == "-" for row in structure_rows),
    }
    summary_rows = [
        {"method": method, "metric": key, "value": value}
        for key, value in summary.items()
        if key != "method"
    ]
    write_tsv(output_dir / "basic_metrics.tsv", summary_rows, ["method", "metric", "value"])
    (output_dir / "annotation_qc.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--all-proteins", type=Path, required=True)
    parser.add_argument("--representative-proteins", type=Path, required=True)
    parser.add_argument("--map", dest="mapping", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--repeat-bed", type=Path)
    parser.add_argument("--edge-bp", type=int, default=500)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    build_qc(
        method=args.method,
        genome_path=args.genome,
        gff_path=args.gff3,
        all_proteins_path=args.all_proteins,
        representative_proteins_path=args.representative_proteins,
        mapping_path=args.mapping,
        excluded_path=args.excluded,
        output_dir=args.output_dir,
        repeat_bed=args.repeat_bed,
        edge_bp=args.edge_bp,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
