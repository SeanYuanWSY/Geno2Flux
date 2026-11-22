#!/usr/bin/env python3
"""Build AlphaGEM target-input and gap-fill coordinate contracts from genome, GFF3, and proteome FASTA.

Validates every accepted protein in ``--proteins`` against the exact coding sequence
spliced from ``--genome`` according to ``--gff`` (translation table 1), emits safe
sequential protein IDs with a reversible ``id_map.tsv`` (or keeps original IDs when
``--keep-original-ids`` is set and all IDs are already safe), and validates the
resulting contracts with the bundled AlphaGEM contract validators.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
ALPHAGEM_CONTRACTS = ROOT / "third_party" / "alphagem" / "src" / "pipeline_contracts.py"
ALPHAGEM_SOURCE_MANIFEST = ROOT / "third_party" / "alphagem" / "SOURCE_MANIFEST.json"
NUCLEOTIDES = set("ACGTRYSWKMBDHVN")
SAFE_PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,31}$")
SAFE_PROTEIN_ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, int | str]:
    return {"bytes": path.stat().st_size, "sha256": sha256(path)}


def contracts_module():
    """Load the AlphaGEM coordinate validator bundled in this repository."""
    manifest = json.loads(ALPHAGEM_SOURCE_MANIFEST.read_text(encoding="utf-8"))
    expected = next(
        row for row in manifest["files"]
        if row["path"] == "src/pipeline_contracts.py"
    )
    if identity(ALPHAGEM_CONTRACTS) != {
        "bytes": expected["bytes"],
        "sha256": expected["sha256"],
    }:
        raise ValueError("Bundled AlphaGEM coordinate validator differs from its source manifest")
    spec = importlib.util.spec_from_file_location("geno2flux_alphagem_contracts", ALPHAGEM_CONTRACTS)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load bundled AlphaGEM coordinate validator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_fasta(path: Path, *, nucleotides: bool) -> dict[str, str]:
    records: dict[str, str] = {}
    record_id: str | None = None
    parts: list[str] = []

    def finish() -> None:
        if record_id is None:
            return
        sequence = "".join(parts).upper()
        if not sequence:
            raise ValueError(f"Empty FASTA record {record_id!r}: {path}")
        if record_id in records:
            raise ValueError(f"Duplicate FASTA ID {record_id!r}: {path}")
        if nucleotides and set(sequence) - NUCLEOTIDES:
            raise ValueError(f"Non-nucleotide genome sequence {record_id!r}")
        if not nucleotides:
            if sequence.endswith("*"):
                sequence = sequence[:-1]
            if not sequence or "*" in sequence or not re.fullmatch(r"[A-Z]+", sequence):
                raise ValueError(f"Invalid protein sequence {record_id!r}")
        records[record_id] = sequence

    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(">"):
                finish()
                record_id = unquote(line[1:].split()[0]) if line[1:].split() else ""
                if not record_id:
                    raise ValueError(f"Blank FASTA ID: {path}")
                parts = []
            elif line.strip():
                if record_id is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                parts.append(line.strip())
    finish()
    if not records:
        raise ValueError(f"Empty FASTA: {path}")
    return records


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in text.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise ValueError(f"Malformed GFF3 attribute {field!r}")
        key, value = field.split("=", 1)
        key = unquote(key)
        if key in result:
            raise ValueError(f"Duplicate GFF3 attribute {key!r}")
        result[key] = value
    return result


def read_gff(path: Path) -> tuple[dict[str, dict], dict[str, list[dict]]]:
    transcripts: dict[str, dict] = {}
    cds_features: list[dict] = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF3 line {line_number} has {len(fields)} fields")
            seqid, _source, kind, start, end, _score, strand, phase, raw_attrs = fields
            if kind not in {"mRNA", "transcript", "CDS"}:
                continue
            attrs = attributes(raw_attrs)
            try:
                left, right = int(start), int(end)
            except ValueError as error:
                raise ValueError(f"Non-integer GFF3 coordinate at line {line_number}") from error
            if left < 1 or right < left or strand not in {"+", "-"}:
                raise ValueError(f"Invalid GFF3 locus at line {line_number}")
            if kind in {"mRNA", "transcript"}:
                identifier = unquote(attrs.get("ID", ""))
                parent_gene = unquote(attrs.get("Parent", ""))
                if not identifier or identifier in transcripts or not parent_gene:
                    raise ValueError(
                        f"Missing/duplicate transcript or gene parent at line {line_number}"
                    )
                transcripts[identifier] = {
                    "seqid": seqid,
                    "strand": strand,
                    "gene": parent_gene,
                }
            else:
                if phase not in {"0", "1", "2"} or not attrs.get("Parent"):
                    raise ValueError(f"Invalid CDS phase/parent at line {line_number}")
                cds_features.append(
                    {
                        "parent": attrs["Parent"],
                        "seqid": seqid,
                        "start": left,
                        "end": right,
                        "strand": strand,
                        "phase": int(phase),
                        "line": line_number,
                    }
                )
    segments: dict[str, list[dict]] = {}
    for cds in cds_features:
        raw_parent = cds.pop("parent")
        decoded_parent = unquote(raw_parent)
        # Preserve predictors such as ANNEVO/Helixer that emit literal commas inside IDs,
        # while still supporting standard comma-separated multi-parent CDS lines.
        if decoded_parent in transcripts:
            parents = [decoded_parent]
        else:
            parents = [unquote(p) for p in raw_parent.split(",") if p]
        if not parents or any(tid not in transcripts for tid in parents):
            raise ValueError(f"Unknown CDS parent at GFF3 line {cds['line']}: {raw_parent!r}")
        for tid in parents:
            segments.setdefault(tid, []).append(cds.copy())
    return transcripts, segments


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def build(
    *,
    genome: Path,
    gff: Path,
    proteins: Path,
    case_id: str,
    output: Path,
    prefix: str = "PROT",
    source_label: str = "annotation",
    keep_original_ids: bool = False,
) -> dict:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", case_id):
        raise ValueError("Unsafe AlphaGEM case ID")
    if not SAFE_PREFIX_RE.fullmatch(prefix):
        raise ValueError(f"Unsafe protein ID prefix: {prefix!r}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Coordinate output already exists: {output}")
    for path in (genome, gff, proteins):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing or symlinked input: {path}")

    observed = {
        name: identity(path)
        for name, path in (("genome", genome), ("gff", gff), ("proteins", proteins))
    }
    vendor = contracts_module()
    genome_records = read_fasta(genome, nucleotides=True)
    protein_records = read_fasta(proteins, nucleotides=False)
    transcripts, cds_by_transcript = read_gff(gff)

    missing_from_gff = sorted(set(protein_records) - set(transcripts))
    if missing_from_gff:
        raise ValueError(
            f"Proteome FASTA contains transcript IDs absent from GFF3 (e.g. {missing_from_gff[:5]})"
        )
    missing_cds = sorted(tid for tid in protein_records if tid not in cds_by_transcript)
    if missing_cds:
        raise ValueError(
            f"Proteome FASTA contains transcripts without CDS features in GFF3 (e.g. {missing_cds[:5]})"
        )

    if keep_original_ids:
        unsafe = [pid for pid in protein_records if not SAFE_PROTEIN_ID_RE.fullmatch(pid)]
        if unsafe:
            raise ValueError(
                f"--keep-original-ids requested, but unsafe FASTA IDs found (e.g. {unsafe[:5]}); "
                "omit --keep-original-ids to normalize them with --prefix"
            )

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=output.name + ".tmp-", dir=output.parent))
    try:
        target_fasta = temporary / "proteome.fasta"
        id_map_path = temporary / "id_map.tsv"
        map_path = temporary / "coordinate-map.tsv"
        segments_path = temporary / "coordinate-segments.tsv"
        target_contract_path = temporary / "target-input-contract.json"
        coordinate_contract_path = temporary / "gapfill-coordinate-contract.json"
        map_rows: list[dict] = []
        segment_rows: list[dict] = []
        id_rows: list[dict] = []
        with target_fasta.open("w", encoding="utf-8") as fasta_stream:
            for index, original_id in enumerate(sorted(protein_records), start=1):
                safe_id = original_id if keep_original_ids else f"{prefix}{index:06d}"
                protein = protein_records[original_id]
                transcript = transcripts[original_id]
                cds_rows = sorted(
                    cds_by_transcript[original_id], key=lambda x: (x["start"], x["end"])
                )
                seqid, strand = transcript["seqid"], transcript["strand"]
                chromosome = genome_records.get(seqid)
                if chromosome is None:
                    raise ValueError(f"GFF3 seqid absent from genome: {seqid}")
                if any(s["seqid"] != seqid or s["strand"] != strand for s in cds_rows):
                    raise ValueError(f"CDS segment disagrees with mRNA: {original_id}")
                for left, right in zip(cds_rows, cds_rows[1:]):
                    if left["end"] >= right["start"]:
                        raise ValueError(f"Overlapping CDS segments: {original_id}")
                if any(s["end"] > len(chromosome) for s in cds_rows):
                    raise ValueError(f"CDS segment beyond genome: {original_id}")
                transcript_order = cds_rows if strand == "+" else list(reversed(cds_rows))
                parts = []
                for segment in transcript_order:
                    part = chromosome[segment["start"] - 1 : segment["end"]]
                    parts.append(part if strand == "+" else vendor.reverse_complement(part))
                coding = "".join(parts)[transcript_order[0]["phase"] :]
                if len(coding) % 3:
                    raise ValueError(f"CDS length is not divisible by three: {original_id}")
                translated = vendor.translate_nucleotide_sequence(
                    coding,
                    translation_table=1,
                    expected_starts_with_methionine=protein.startswith("M"),
                )
                if translated != protein:
                    raise ValueError(f"Genome CDS translation differs from FAA: {original_id}")
                sequence_hash = hashlib.sha256(protein.encode("utf-8")).hexdigest()
                fasta_stream.write(f">{safe_id}\n")
                for offset in range(0, len(protein), 60):
                    fasta_stream.write(protein[offset : offset + 60] + "\n")
                id_rows.append(
                    {
                        "target_protein_id": safe_id,
                        "original_protein_id": original_id,
                        "original_transcript_id": original_id,
                        "predicted_gene_id": transcript["gene"],
                        "genome_seqid": seqid,
                        "strand": strand,
                        "cds_start": cds_rows[0]["start"],
                        "cds_end": cds_rows[-1]["end"],
                        "sequence_sha256": sequence_hash,
                    }
                )
                predicted_cds_id = f"predicted_cds:{original_id}"
                map_rows.append(
                    {
                        "target_protein_id": safe_id,
                        "genome_seqid": seqid,
                        "start": cds_rows[0]["start"],
                        "end": cds_rows[-1]["end"],
                        "strand": strand,
                        "refseq_protein_id": original_id,
                        "refseq_cds_id": predicted_cds_id,
                        "cds_segment_count": len(cds_rows),
                        "cds_segments": ";".join(f"{s['start']}-{s['end']}" for s in cds_rows),
                        "locus_tag": transcript["gene"],
                        "source_database_ids": f"{source_label}_gff3:{original_id}",
                        "mapping_evidence": "same_run_tool_prediction_gff3;strict_genome_translation",
                        "mapping_bucket": "tool_prediction_strict_translation",
                        "local_sequence_sha256": sequence_hash,
                        "reference_sequence_sha256": sequence_hash,
                    }
                )
                for segment_index, segment in enumerate(cds_rows, start=1):
                    segment_rows.append(
                        {
                            "target_protein_id": safe_id,
                            "refseq_protein_id": original_id,
                            "refseq_cds_id": predicted_cds_id,
                            "genome_seqid": seqid,
                            "strand": strand,
                            "segment_index_genomic_order": segment_index,
                            "segment_start": segment["start"],
                            "segment_end": segment["end"],
                            "phase": segment["phase"],
                            "locus_tag": transcript["gene"],
                            "translation_table": 1,
                        }
                    )

        write_tsv(
            id_map_path,
            (
                "target_protein_id",
                "original_protein_id",
                "original_transcript_id",
                "predicted_gene_id",
                "genome_seqid",
                "strand",
                "cds_start",
                "cds_end",
                "sequence_sha256",
            ),
            id_rows,
        )
        write_tsv(map_path, vendor.GAPFILL_COORDINATE_COLUMNS, map_rows)
        write_tsv(segments_path, vendor.GAPFILL_SEGMENT_COLUMNS, segment_rows)
        write_json(
            target_contract_path,
            {
                "schema_version": 1,
                "case_id": case_id,
                "source_fasta": vendor.target_fasta_identity(target_fasta),
                "provenance": {
                    "role": "standardized coding proteome",
                    "source_label": source_label,
                    "source_faa_sha256": observed["proteins"]["sha256"],
                    "id_map_sha256": sha256(id_map_path),
                },
            },
        )
        write_json(
            coordinate_contract_path,
            {
                "schema_version": 1,
                "case_id": case_id,
                "coordinate_semantics": {
                    "coordinate_system": "1-based-inclusive",
                    "segment_order": "genomic_start_ascending",
                    "hit_overlap": "cds_segment_union",
                },
                "inputs": {
                    "target_fasta": identity(target_fasta),
                    "raw_genome": identity(genome),
                    "coordinate_map": identity(map_path),
                    "coordinate_segments": identity(segments_path),
                },
                "provenance": {
                    "audit_script_sha256": sha256(Path(__file__)),
                    "resource_manifest_sha256": sha256(ALPHAGEM_SOURCE_MANIFEST),
                    "audit_command": [
                        "build_coordinate_contract.py",
                        "--case-id",
                        case_id,
                        "--prefix",
                        prefix,
                    ],
                    "official_sources": [
                        f"Coding GFF3 ({source_label}); sha256={observed['gff']['sha256']}",
                        f"Genome assembly; sha256={observed['genome']['sha256']}",
                    ],
                    "source_type": "user_supplied_gff3_strict_translation",
                },
            },
        )
        vendor.validate_target_input_contract(target_contract_path, target_fasta)
        vendor.validate_gapfill_coordinate_contract(
            coordinate_contract_path,
            fasta_path=target_fasta,
            genome_path=genome,
            coordinate_map_path=map_path,
            coordinate_segments_path=segments_path,
            expected_case_id=case_id,
        )
        audited = vendor.validate_gapfill_genome(genome, target_fasta, map_path, segments_path)
        if (
            audited["recovery_translation_mismatches"]
            or audited["target_proteins"] != len(protein_records)
        ):
            raise ValueError(
                "AlphaGEM coordinate validator did not strictly accept every target protein"
            )
        # Keep a copy of the genome alongside the contracts so --inputs directory is self-contained
        shutil.copy2(genome, temporary / "genome.fasta")
        receipt = {
            "status": "strict_contract_validated",
            "source_label": source_label,
            "case_id": case_id,
            "prefix": prefix,
            "keep_original_ids": keep_original_ids,
            "input_identities": observed,
            "gff_total_transcripts": len(transcripts),
            "target_proteins": len(protein_records),
            "non_selected_gff_transcripts": len(set(transcripts) - set(protein_records)),
            "coordinate_records": audited["coordinate_records"],
            "coordinate_segment_records": audited["coordinate_segment_records"],
            "strict_translation_matches": len(protein_records),
            "translation_table": 1,
            "files": {
                p.name: identity(p)
                for p in sorted(temporary.iterdir())
                if p.is_file()
            },
        }
        write_json(temporary / "receipt.json", receipt)
        os.rename(temporary, output)
        return receipt
    except Exception:
        shutil.rmtree(temporary)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, required=True, help="Genome assembly FASTA")
    parser.add_argument("--gff", type=Path, required=True, help="Coding gene/mRNA/CDS GFF3")
    parser.add_argument(
        "--proteins",
        type=Path,
        required=True,
        help="Representative protein FASTA (e.g. longest_isoform.faa)",
    )
    parser.add_argument("--case-id", default="sample", help="Reconstruction/case identifier")
    parser.add_argument("--prefix", default="PROT", help="Safe protein ID prefix (default: PROT)")
    parser.add_argument("--source-label", default="annotation", help="Annotation method/source tag")
    parser.add_argument(
        "--keep-original-ids",
        action="store_true",
        help="Keep original FASTA IDs when they already match safe identifier syntax",
    )
    parser.add_argument(
        "--output",
        "--output-dir",
        dest="output",
        type=Path,
        required=True,
        help="New output directory",
    )
    args = parser.parse_args(argv)
    receipt = build(
        genome=args.genome,
        gff=args.gff,
        proteins=args.proteins,
        case_id=args.case_id,
        output=args.output,
        prefix=args.prefix,
        source_label=args.source_label,
        keep_original_ids=args.keep_original_ids,
    )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "case_id": args.case_id,
                "proteins": receipt["target_proteins"],
                "segments": receipt["coordinate_segment_records"],
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
