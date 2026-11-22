#!/usr/bin/env python3
"""Standardize a coding GFF3 with one strict, method-neutral extraction path.

The script uses gffread for sequence extraction, then independently translates
the extracted CDS with translation table 1. Invalid models are reported with a
closed reason-code vocabulary instead of being silently discarded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from urllib.parse import quote, unquote

REASON_CODES = {
    "MISSING_PARENT",
    "DUPLICATE_ID",
    "INVALID_COORDINATE",
    "INVALID_STRAND",
    "INVALID_PHASE",
    "EMPTY_CDS",
    "CDS_NOT_MULTIPLE_OF_3",
    "INTERNAL_STOP",
    "EMPTY_PROTEIN",
    "INVALID_AMINO_ACID",
    "TRANSLATION_MISMATCH",
    "UNRESOLVED_GFF_GRAPH",
}
TRANSCRIPT_TYPES = {"mRNA", "transcript"}
TRANSCRIPT_CHILD_TYPES = {"CDS", "exon", "five_prime_UTR", "three_prime_UTR", "UTR"}
MULTI_VALUE_ATTRS = {"Alias", "Dbxref", "Note", "Ontology_term"}
DNA = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
AA = set("ACDEFGHIKLMNPQRSTVWYX*")
CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class StandardizationError(RuntimeError):
    pass


@dataclass(frozen=True)
class Feature:
    line_number: int
    seqid: str
    source: str
    kind: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict[str, tuple[str, ...]]

    @property
    def identifier(self) -> str | None:
        values = self.attrs.get("ID", ())
        return values[0] if len(values) == 1 else None

    @property
    def parents(self) -> tuple[str, ...]:
        return self.attrs.get("Parent", ())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_attrs(raw: str) -> dict[str, tuple[str, ...]]:
    attrs: dict[str, tuple[str, ...]] = {}
    if raw == ".":
        return attrs
    if raw == "":
        raise StandardizationError("malformed GFF3 attributes: empty field")
    items = raw.split(";")
    if items[-1] == "":
        items = items[:-1]
    for item in items:
        if not item or "=" not in item:
            raise StandardizationError(f"malformed GFF3 attributes: {raw}")
        key_raw, value_raw = item.split("=", 1)
        key = unquote(key_raw)
        if key in attrs:
            raise StandardizationError(f"duplicate GFF3 attribute key: {key}")
        # Parent is resolved after the complete entity-ID set is known. This
        # preserves Helixer IDs that contain literal commas while still
        # supporting standards-compliant comma-separated Parent lists. Other
        # standard multi-value fields follow normal GFF3 comma semantics.
        if key == "Parent":
            values = (value_raw,)
        elif key in MULTI_VALUE_ATTRS:
            values = tuple(unquote(value) for value in value_raw.split(","))
        else:
            values = (unquote(value_raw),)
        if not key or not values or any(not value for value in values):
            raise StandardizationError(f"empty GFF3 attribute: {raw}")
        attrs[key] = values
    return attrs


def _format_attrs(attrs: dict[str, tuple[str, ...]]) -> str:
    if not attrs:
        return "."
    safe = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:^*$@!+_?-|"
    return ";".join(
        f"{quote(key, safe=safe)}={','.join(quote(value, safe=safe) for value in values)}"
        for key, values in attrs.items()
    )


def _read_fasta(path: Path, *, decode_ids: bool = False) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                fields = line[1:].split()
                current = fields[0] if fields else ""
                if decode_ids:
                    current = unquote(current)
                if not current or current in records:
                    raise StandardizationError(f"invalid or duplicate FASTA ID in {path}")
                chunks = []
            elif line:
                if current is None or any(character.isspace() for character in line):
                    raise StandardizationError(f"invalid FASTA structure in {path}")
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    if not records:
        raise StandardizationError(f"empty FASTA: {path}")
    return records


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")


def _parse_gff(path: Path) -> tuple[list[str], list[Feature]]:
    directives: list[str] = []
    features: list[Feature] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                directives.append(line)
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise StandardizationError(f"GFF3 line {line_number} does not have 9 columns")
            seqid, source, kind, start_raw, end_raw, score, strand, phase, attrs_raw = fields
            try:
                start, end = int(start_raw), int(end_raw)
            except ValueError as exc:
                raise StandardizationError(
                    f"GFF3 line {line_number} has invalid coordinates"
                ) from exc
            features.append(
                Feature(
                    line_number,
                    seqid,
                    source,
                    kind,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    _parse_attrs(attrs_raw),
                )
            )
    if not features:
        raise StandardizationError(f"GFF3 contains no features: {path}")
    entity_ids = {
        feature.identifier
        for feature in features
        if feature.kind in {"gene", *TRANSCRIPT_TYPES} and feature.identifier
    }
    resolved: list[Feature] = []
    for feature in features:
        parents = feature.attrs.get("Parent")
        if not parents:
            resolved.append(feature)
            continue
        raw_parent = parents[0]
        decoded_parent = unquote(raw_parent)
        split_parents = tuple(unquote(parent) for parent in raw_parent.split(","))
        if decoded_parent in entity_ids:
            new_parents = (decoded_parent,)
        elif len(split_parents) > 1:
            new_parents = split_parents
        else:
            new_parents = (decoded_parent,)
        attrs = dict(feature.attrs)
        attrs["Parent"] = new_parents
        resolved.append(replace(feature, attrs=attrs))
    return directives, resolved


def _reason(
    exclusions: dict[str, dict[str, str]], transcript_id: str, code: str, details: str
) -> None:
    if code not in REASON_CODES:
        raise AssertionError(code)
    exclusions[transcript_id].setdefault(code, details)


def _feature_line(feature: Feature, parents: tuple[str, ...] | None = None) -> str:
    attrs = dict(feature.attrs)
    if parents is not None:
        attrs["Parent"] = parents
    return "\t".join(
        (
            feature.seqid,
            feature.source,
            feature.kind,
            str(feature.start),
            str(feature.end),
            feature.score,
            feature.strand,
            feature.phase,
            _format_attrs(attrs),
        )
    )


def _write_subset_gff(
    path: Path,
    assembly: dict[str, str],
    features: list[Feature],
    transcript_to_gene: dict[str, str],
    kept_transcripts: set[str],
) -> None:
    kept_genes = {transcript_to_gene[transcript] for transcript in kept_transcripts}
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        for seqid, sequence in assembly.items():
            handle.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
        for feature in features:
            identifier = feature.identifier
            if feature.kind == "gene":
                if identifier in kept_genes:
                    handle.write(_feature_line(feature) + "\n")
                continue
            if feature.kind in TRANSCRIPT_TYPES:
                if identifier in kept_transcripts:
                    handle.write(_feature_line(feature) + "\n")
                continue
            parents = tuple(parent for parent in feature.parents if parent in kept_transcripts)
            if parents:
                handle.write(_feature_line(feature, parents) + "\n")


def _run_gffread(
    gffread: Path,
    genome: Path,
    input_gff: Path,
    output_gff: Path,
    transcripts: Path,
    cds: Path,
    proteins: Path,
    stderr_path: Path,
) -> None:
    command = [
        str(gffread),
        str(input_gff),
        "-g",
        str(genome),
        "-F",
        "--keep-exon-attrs",
        "--keep-genes",
        "-O",
        "-P",
        "-S",
        "-o",
        str(output_gff),
        "-w",
        str(transcripts),
        "-x",
        str(cds),
        "-y",
        str(proteins),
        "-E",
    ]
    with stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
        completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=stderr, check=False)
    if completed.returncode != 0:
        raise StandardizationError(
            f"gffread failed with exit code {completed.returncode}; see {stderr_path}"
        )


def _canonicalize_gffread_gff(path: Path) -> None:
    """Remove gffread's volatile command comment while preserving directives."""
    lines = path.read_text(encoding="utf-8").splitlines()
    stable = [line for line in lines if not line.startswith("# ")]
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(stable) + "\n")


def _canonicalize_log(path: Path, replacements: dict[str, str]) -> None:
    text = path.read_text(encoding="utf-8")
    for source, replacement in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(source, replacement)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _extract_spliced_python(
    assembly: dict[str, str],
    features: list[Feature],
    kept_transcripts: set[str],
) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Pure-Python fallback when external gffread is not installed."""
    transcripts_by_id = {
        f.identifier: f for f in features if f.kind in TRANSCRIPT_TYPES and f.identifier
    }
    cds_by_parent: dict[str, list[Feature]] = defaultdict(list)
    exon_by_parent: dict[str, list[Feature]] = defaultdict(list)
    for f in features:
        if f.kind == "CDS":
            for p in f.parents:
                cds_by_parent[p].append(f)
        elif f.kind == "exon":
            for p in f.parents:
                exon_by_parent[p].append(f)

    tx_seqs: dict[str, str] = {}
    cds_seqs: dict[str, str] = {}
    prot_seqs: dict[str, str] = {}
    for tid in sorted(kept_transcripts):
        tx = transcripts_by_id[tid]
        chrom = assembly[tx.seqid]
        exons = exon_by_parent.get(tid) or cds_by_parent.get(tid, [])
        ordered_exons = sorted(exons, key=lambda x: x.start, reverse=(tx.strand == "-"))
        tx_parts = [
            chrom[e.start - 1 : e.end] if tx.strand == "+" else _revcomp(chrom[e.start - 1 : e.end])
            for e in ordered_exons
        ]
        tx_seqs[tid] = "".join(tx_parts)

        cds_list = sorted(cds_by_parent.get(tid, []), key=lambda x: x.start, reverse=(tx.strand == "-"))
        cds_parts = [
            chrom[c.start - 1 : c.end] if tx.strand == "+" else _revcomp(chrom[c.start - 1 : c.end])
            for c in cds_list
        ]
        first_phase = int(cds_list[0].phase) if cds_list and cds_list[0].phase in {"0", "1", "2"} else 0
        spliced_cds = "".join(cds_parts)[first_phase:]
        cds_seqs[tid] = spliced_cds
        prot_seqs[tid] = _translate(spliced_cds) if len(spliced_cds) % 3 == 0 and spliced_cds else ""
    return tx_seqs, cds_seqs, prot_seqs

def _translate(cds: str) -> str:
    amino: list[str] = []
    for index in range(0, len(cds), 3):
        codon = cds[index : index + 3].upper()
        amino.append(CODONS[codon] if set(codon) <= set("ACGT") else "X")
    return "".join(amino)


def _is_partial(feature: Feature, translated: str) -> bool:
    markers = {"partial", "start_range", "end_range"}
    marked = any(key in feature.attrs for key in markers)
    return marked or not translated.startswith("M") or not translated.endswith("*")


def _validate_structure(
    assembly: dict[str, str], features: list[Feature]
) -> tuple[
    dict[str, Feature],
    dict[str, str],
    set[str],
    dict[str, dict[str, str]],
]:
    gene_occurrences: dict[str, list[Feature]] = defaultdict(list)
    transcript_occurrences: dict[str, list[Feature]] = defaultdict(list)
    cds_by_parent: dict[str, list[Feature]] = defaultdict(list)
    child_by_parent: dict[str, list[Feature]] = defaultdict(list)
    orphan_children: list[Feature] = []
    exclusions: dict[str, dict[str, str]] = defaultdict(dict)
    for feature in features:
        identifier = feature.identifier
        if feature.kind == "gene":
            if identifier:
                gene_occurrences[identifier].append(feature)
            else:
                _reason(
                    exclusions,
                    f"line_{feature.line_number}",
                    "UNRESOLVED_GFF_GRAPH",
                    "gene feature has no unique ID",
                )
        elif feature.kind in TRANSCRIPT_TYPES:
            if identifier:
                transcript_occurrences[identifier].append(feature)
            else:
                _reason(
                    exclusions,
                    f"line_{feature.line_number}",
                    "UNRESOLVED_GFF_GRAPH",
                    "transcript feature has no unique ID",
                )
        if feature.kind == "CDS":
            for parent in feature.parents:
                cds_by_parent[parent].append(feature)
        if feature.kind not in {"gene", *TRANSCRIPT_TYPES}:
            if not feature.parents and feature.kind in TRANSCRIPT_CHILD_TYPES:
                orphan_children.append(feature)
            for parent in feature.parents:
                child_by_parent[parent].append(feature)

    genes = {identifier: occurrences[0] for identifier, occurrences in gene_occurrences.items()}
    transcripts = {
        identifier: occurrences[0]
        for identifier, occurrences in transcript_occurrences.items()
    }
    duplicate_ids = {
        identifier
        for identifier in set(gene_occurrences) | set(transcript_occurrences)
        if len(gene_occurrences.get(identifier, ()))
        + len(transcript_occurrences.get(identifier, ()))
        > 1
    }
    transcript_to_gene: dict[str, str] = {}
    coding_transcripts = set(transcripts) | {
        parent
        for feature in features
        if feature.kind in TRANSCRIPT_CHILD_TYPES
        for parent in feature.parents
    }
    for orphan in orphan_children:
        _reason(
            exclusions,
            f"line_{orphan.line_number}",
            "MISSING_PARENT",
            f"{orphan.kind} feature has no Parent",
        )
    for gene_id in sorted(genes):
        children = [
            transcript_id
            for transcript_id, transcript in transcripts.items()
            if transcript.parents == (gene_id,)
        ]
        if gene_id in duplicate_ids:
            for transcript_id in children:
                _reason(
                    exclusions,
                    transcript_id,
                    "DUPLICATE_ID",
                    "parent gene ID occurs more than once",
                )
        if not children:
            marker = f"gene::{gene_id}"
            transcript_to_gene[marker] = gene_id
            _reason(
                exclusions,
                marker,
                "DUPLICATE_ID" if gene_id in duplicate_ids else "UNRESOLVED_GFF_GRAPH",
                "gene ID occurs more than once"
                if gene_id in duplicate_ids
                else "gene has no transcript child",
            )
    for transcript_id in sorted(coding_transcripts):
        transcript = transcripts.get(transcript_id)
        if transcript is None:
            _reason(exclusions, transcript_id, "MISSING_PARENT", "child Parent is not a transcript")
            continue
        if transcript_id in duplicate_ids:
            _reason(
                exclusions,
                transcript_id,
                "DUPLICATE_ID",
                "transcript ID occurs more than once",
            )
        if len(transcript.parents) != 1 or transcript.parents[0] not in genes:
            _reason(
                exclusions,
                transcript_id,
                "MISSING_PARENT",
                "transcript Parent is not one gene",
            )
            continue
        gene_id = transcript.parents[0]
        transcript_to_gene[transcript_id] = gene_id
        gene = genes[gene_id]
        if gene_id in duplicate_ids:
            _reason(
                exclusions,
                transcript_id,
                "DUPLICATE_ID",
                "parent gene ID occurs more than once",
            )
        if (
            transcript.strand not in {"+", "-"}
            or gene.strand not in {"+", "-"}
            or gene.strand != transcript.strand
        ):
            _reason(
                exclusions,
                transcript_id,
                "INVALID_STRAND",
                "gene/transcript strand is invalid",
            )
        if (
            transcript.seqid not in assembly
            or gene.seqid not in assembly
            or gene.seqid != transcript.seqid
            or gene.start < 1
            or gene.end < gene.start
            or gene.end > len(assembly.get(gene.seqid, ""))
            or transcript.start < 1
            or transcript.end < transcript.start
            or transcript.end > len(assembly.get(transcript.seqid, ""))
            or gene.start > transcript.start
            or gene.end < transcript.end
        ):
            _reason(
                exclusions,
                transcript_id,
                "INVALID_COORDINATE",
                "gene/transcript coordinates are invalid",
            )
        for child in child_by_parent.get(transcript_id, []):
            if (
                child.seqid != transcript.seqid
                or child.start < transcript.start
                or child.end > transcript.end
                or child.start < 1
                or child.end < child.start
                or child.end > len(assembly.get(child.seqid, ""))
            ):
                _reason(
                    exclusions,
                    transcript_id,
                    "INVALID_COORDINATE",
                    f"{child.kind} coordinates are invalid",
                )
            if child.strand not in {"+", "-"} or child.strand != transcript.strand:
                _reason(
                    exclusions,
                    transcript_id,
                    "INVALID_STRAND",
                    f"{child.kind} strand is invalid",
                )
        segments = cds_by_parent.get(transcript_id, [])
        if not segments:
            _reason(exclusions, transcript_id, "EMPTY_CDS", "transcript has no CDS feature")
            continue
        ordered = sorted(segments, key=lambda item: item.start, reverse=transcript.strand == "-")
        cumulative = 0
        first_phase = 0
        for index, segment in enumerate(ordered):
            if segment.phase not in {"0", "1", "2"}:
                _reason(exclusions, transcript_id, "INVALID_PHASE", "CDS phase is not 0, 1, or 2")
                continue
            phase = int(segment.phase)
            if index == 0:
                first_phase = phase
            expected = first_phase if index == 0 else (3 - cumulative % 3) % 3
            if phase != expected:
                _reason(exclusions, transcript_id, "INVALID_PHASE", "CDS phase continuity failed")
            cumulative += segment.end - segment.start + 1 - (first_phase if index == 0 else 0)
        genomic_order = sorted(segments, key=lambda item: item.start)
        if any(
            genomic_order[index].end >= genomic_order[index + 1].start
            for index in range(len(genomic_order) - 1)
        ):
            _reason(exclusions, transcript_id, "INVALID_COORDINATE", "CDS segments overlap")

    structurally_valid = {
        transcript_id
        for transcript_id in transcript_to_gene
        if transcript_id not in exclusions
    }
    return transcripts, transcript_to_gene, structurally_valid, exclusions


def _resolve_gffread(gffread: Path | str | None) -> Path | None:
    if gffread is None or str(gffread).strip() in {"", "builtin", "python"}:
        return None
    raw = str(gffread).strip()
    if raw == "auto":
        found = shutil.which("gffread")
        return Path(found) if found else None
    candidate = Path(raw)
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate
    found = shutil.which(raw)
    if found:
        return Path(found)
    raise StandardizationError(f"Configured gffread executable not found: {gffread}")


def standardize(
    genome_path: Path,
    gff_path: Path,
    method: str,
    output_dir: Path,
    gffread: Path | str | None = "auto",
    *,
    prefix: str | None = None,
) -> dict[str, object]:
    if not re.fullmatch(r"[a-z][a-z0-9_]*", method):
        raise StandardizationError("method must match [a-z][a-z0-9_]*")
    if prefix is None:
        prefix = f"standardized.{method}"
    elif not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", prefix):
        raise StandardizationError("prefix must match [A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise StandardizationError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_gffread = _resolve_gffread(gffread)
    source_hash_before = _sha256(genome_path)
    source_fai = Path(f"{genome_path}.fai")
    source_fai_existed = source_fai.exists()
    assembly = _read_fasta(genome_path)
    if any(not sequence or set(sequence) - DNA for sequence in assembly.values()):
        raise StandardizationError("assembly contains empty records or unsupported DNA symbols")
    _, features = _parse_gff(gff_path)
    transcripts, transcript_to_gene, structural, exclusions = _validate_structure(
        assembly, features
    )

    with tempfile.TemporaryDirectory(prefix="geno2flux-standardize-", dir=output_dir) as tmp_raw:
        tmp = Path(tmp_raw)
        genome_copy = tmp / "genome.fasta"
        shutil.copyfile(genome_path, genome_copy)
        if _sha256(genome_copy) != source_hash_before:
            raise StandardizationError("run-local assembly copy does not match the source")
        candidate_gff = tmp / "candidate.gff3"
        _write_subset_gff(candidate_gff, assembly, features, transcript_to_gene, structural)
        if resolved_gffread is not None:
            normalized_gff = tmp / "candidate.normalized.gff3"
            candidate_tx = tmp / "candidate.transcripts.fna"
            candidate_cds = tmp / "candidate.cds.fna"
            candidate_proteins = tmp / "candidate.proteins.faa"
            candidate_stderr = output_dir / "gffread.candidate.stderr.txt"
            _run_gffread(
                resolved_gffread,
                genome_copy,
                candidate_gff,
                normalized_gff,
                candidate_tx,
                candidate_cds,
                candidate_proteins,
                candidate_stderr,
            )
            _canonicalize_log(
                candidate_stderr,
                {str(tmp): "<TMP>", str(output_dir): "<OUTPUT_DIR>"},
            )
            _read_fasta(candidate_tx, decode_ids=True)
            extracted_cds = _read_fasta(candidate_cds, decode_ids=True)
            extracted_proteins = _read_fasta(candidate_proteins, decode_ids=True)
        else:
            _, extracted_cds, extracted_proteins = _extract_spliced_python(
                assembly, features, structural
            )
        for transcript_id in structural:
            if transcript_id not in extracted_cds or transcript_id not in extracted_proteins:
                _reason(
                    exclusions,
                    transcript_id,
                    "UNRESOLVED_GFF_GRAPH",
                    "gffread did not emit CDS and protein for structurally valid transcript",
                )
                continue
            cds = extracted_cds[transcript_id]
            protein = extracted_proteins[transcript_id]
            if not cds:
                _reason(exclusions, transcript_id, "EMPTY_CDS", "extracted CDS is empty")
                continue
            if len(cds) % 3:
                _reason(
                    exclusions,
                    transcript_id,
                    "CDS_NOT_MULTIPLE_OF_3",
                    f"extracted CDS length is {len(cds)}",
                )
                continue
            if not protein:
                _reason(exclusions, transcript_id, "EMPTY_PROTEIN", "extracted protein is empty")
                continue
            if set(protein) - AA:
                _reason(
                    exclusions,
                    transcript_id,
                    "INVALID_AMINO_ACID",
                    "extracted protein has unsupported symbols",
                )
                continue
            translated = _translate(cds)
            normalized_protein = protein[:-1] if protein.endswith("*") else protein
            normalized_translation = translated[:-1] if translated.endswith("*") else translated
            if "*" in normalized_protein or "*" in normalized_translation:
                _reason(
                    exclusions,
                    transcript_id,
                    "INTERNAL_STOP",
                    "translation has an internal stop",
                )
                continue
            if normalized_protein != normalized_translation:
                _reason(
                    exclusions,
                    transcript_id,
                    "TRANSLATION_MISMATCH",
                    "gffread protein differs from independent table-1 translation",
                )

        valid = {transcript_id for transcript_id in structural if transcript_id not in exclusions}
        if not valid:
            raise StandardizationError("no valid coding transcript remains after standardization")
        final_candidate = tmp / "final_candidate.gff3"
        _write_subset_gff(final_candidate, assembly, features, transcript_to_gene, valid)
        final_gff = output_dir / f"{prefix}.standardized.gff3"
        if resolved_gffread is not None:
            final_tx_raw = tmp / "final.transcripts.fna"
            final_cds_raw = tmp / "final.cds.fna"
            final_proteins_raw = tmp / "final.proteins.faa"
            final_stderr = output_dir / "gffread.final.stderr.txt"
            _run_gffread(
                resolved_gffread,
                genome_copy,
                final_candidate,
                final_gff,
                final_tx_raw,
                final_cds_raw,
                final_proteins_raw,
                final_stderr,
            )
            _canonicalize_gffread_gff(final_gff)
            _canonicalize_log(
                final_stderr,
                {str(tmp): "<TMP>", str(output_dir): "<OUTPUT_DIR>"},
            )
            final_tx = _read_fasta(final_tx_raw, decode_ids=True)
            final_cds = _read_fasta(final_cds_raw, decode_ids=True)
            final_proteins_with_stops = _read_fasta(final_proteins_raw, decode_ids=True)
        else:
            shutil.copy2(final_candidate, final_gff)
            final_tx, final_cds, final_proteins_with_stops = _extract_spliced_python(
                assembly, features, valid
            )

    if set(final_cds) != valid or set(final_proteins_with_stops) != valid:
        raise StandardizationError(
            "final gffread output does not exactly match the valid transcript set"
        )
    if set(final_tx) != valid:
        raise StandardizationError(
            "final transcript FASTA does not exactly match the valid transcript set"
        )
    if _sha256(genome_path) != source_hash_before:
        raise StandardizationError("source assembly changed during standardization")
    if not source_fai_existed and source_fai.exists():
        raise StandardizationError("gffread created an index beside the source assembly")

    proteins = {
        identifier: sequence[:-1] if sequence.endswith("*") else sequence
        for identifier, sequence in final_proteins_with_stops.items()
    }
    complete: dict[str, bool] = {}
    partial: dict[str, bool] = {}
    for transcript_id in valid:
        translated = _translate(final_cds[transcript_id])
        complete[transcript_id] = translated.startswith("M") and translated.endswith("*")
        partial[transcript_id] = _is_partial(transcripts[transcript_id], translated)
        if proteins[transcript_id] != (translated[:-1] if translated.endswith("*") else translated):
            raise StandardizationError(f"final translation mismatch for {transcript_id}")

    representatives: dict[str, str] = {}
    by_gene: dict[str, list[str]] = defaultdict(list)
    for transcript_id in valid:
        by_gene[transcript_to_gene[transcript_id]].append(transcript_id)
    for gene_id, choices in by_gene.items():
        representatives[gene_id] = sorted(
            choices,
            key=lambda transcript_id: (
                -int(complete[transcript_id]),
                int(partial[transcript_id]),
                -len(proteins[transcript_id]),
                transcript_id,
            ),
        )[0]

    ordered_valid = sorted(valid)
    _write_fasta(
        output_dir / f"{prefix}.cds.fna",
        [(identifier, final_cds[identifier]) for identifier in ordered_valid],
    )
    _write_fasta(
        output_dir / f"{prefix}.transcripts.fna",
        [(identifier, final_tx[identifier]) for identifier in ordered_valid],
    )
    _write_fasta(
        output_dir / f"{prefix}.proteins.all_isoforms.faa",
        [(identifier, proteins[identifier]) for identifier in ordered_valid],
    )
    representative_ids = sorted(representatives.values())
    _write_fasta(
        output_dir / f"{prefix}.proteins.longest_isoform.faa",
        [(identifier, proteins[identifier]) for identifier in representative_ids],
    )

    mapping_path = output_dir / f"{prefix}.gene_transcript_protein_map.tsv"
    with mapping_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "method\tgene_id\ttranscript_id\tprotein_id\trepresentative\t"
            "complete_translation\tpartial\tcds_length\tprotein_length\n"
        )
        for transcript_id in ordered_valid:
            gene_id = transcript_to_gene[transcript_id]
            handle.write(
                f"{method}\t{gene_id}\t{transcript_id}\t{transcript_id}\t"
                f"{str(representatives[gene_id] == transcript_id).lower()}\t"
                f"{str(complete[transcript_id]).lower()}\t{str(partial[transcript_id]).lower()}\t"
                f"{len(final_cds[transcript_id])}\t{len(proteins[transcript_id])}\n"
            )

    excluded_path = output_dir / f"{prefix}.excluded_models.tsv"
    with excluded_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("method\tgene_id\ttranscript_id\tprotein_id\treason_code\tdetails\n")
        for transcript_id in sorted(exclusions):
            gene_id = transcript_to_gene.get(transcript_id, "")
            for code, details in sorted(exclusions[transcript_id].items()):
                handle.write(
                    f"{method}\t{gene_id}\t{transcript_id}\t{transcript_id}\t{code}\t"
                    f"{details.replace(chr(9), ' ')}\n"
                )

    outputs = sorted(path for path in output_dir.iterdir() if path.is_file())
    summary: dict[str, object] = {
        "method": method,
        "assembly_sha256": source_hash_before,
        "input_gff3_sha256": _sha256(gff_path),
        "gffread_path": str(resolved_gffread) if resolved_gffread is not None else "builtin_python",
        "gffread_version": (
            subprocess.run(
                [str(resolved_gffread), "--version"], capture_output=True, text=True, check=True
            ).stdout.strip()
            if resolved_gffread is not None
            else "builtin_python_table1"
        ),
        "input_gene_count": sum(feature.kind == "gene" for feature in features),
        "input_transcript_count": sum(feature.kind in TRANSCRIPT_TYPES for feature in features),
        "input_cds_feature_count": sum(feature.kind == "CDS" for feature in features),
        "included_gene_count": len(representatives),
        "included_transcript_count": len(valid),
        "representative_protein_count": len(representatives),
        "excluded_transcript_count": len(exclusions),
        "translation_match_count": len(valid),
        "output_sha256": {path.name: _sha256(path) for path in outputs},
    }
    summary_path = output_dir / f"{prefix}.standardization_qc.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--gff3", required=True, type=Path)
    parser.add_argument("--method", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--gffread",
        default="auto",
        help="Path to gffread executable, 'auto' (use gffread if in PATH else built-in Python splicer), or 'builtin'",
    )
    parser.add_argument("--prefix", help="Output filename prefix (default: standardized.<method>)")
    args = parser.parse_args()
    summary = standardize(
        args.genome,
        args.gff3,
        args.method,
        args.output_dir,
        args.gffread,
        prefix=args.prefix,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
