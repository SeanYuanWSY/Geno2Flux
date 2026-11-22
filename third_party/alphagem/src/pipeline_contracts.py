"""Small, dependency-free contracts shared by repaired AlphaGEM stages."""

import pickle
import hashlib
import json
import csv
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Dict, Iterable, List, Set


JOB_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")

STANDARD_CODON_TABLE = {
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

GENETIC_CODE_OVERRIDES = {
    1: {},
    4: {"TGA": "W"},
    11: {},
    12: {"CTG": "S"},
}

GENETIC_CODE_STARTS = {
    1: {"ATG", "CTG", "TTG"},
    4: {"ATA", "ATC", "ATG", "ATT", "GTG", "TTG"},
    11: {"ATA", "ATC", "ATG", "ATT", "CTG", "GTG", "TTG"},
    12: {"ATG", "CTG"},
}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(
        str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")
    )[::-1]


def translate_nucleotide_sequence(
    sequence: str,
    *,
    translation_table: int,
    expected_starts_with_methionine: bool,
    amino_acid_overrides=None,
) -> str:
    if translation_table not in GENETIC_CODE_OVERRIDES:
        raise ValueError(
            f"Unsupported audited translation table: {translation_table}"
        )
    if len(sequence) % 3:
        raise ValueError(
            f"CDS length is not divisible by three: {len(sequence)}"
        )
    code = {
        **STANDARD_CODON_TABLE,
        **GENETIC_CODE_OVERRIDES[translation_table],
    }
    amino_acids = []
    for offset in range(0, len(sequence), 3):
        codon = sequence[offset : offset + 3]
        amino_acids.append(code.get(codon, "X"))
    if (
        amino_acids
        and expected_starts_with_methionine
        and sequence[:3] in GENETIC_CODE_STARTS[translation_table]
    ):
        amino_acids[0] = "M"
    for index, amino_acid in dict(amino_acid_overrides or {}).items():
        if (
            not isinstance(index, int)
            or index < 0
            or index >= len(amino_acids)
            or not isinstance(amino_acid, str)
            or not re.fullmatch(r"[A-Z*]", amino_acid)
        ):
            raise ValueError(
                "Invalid audited amino-acid translation override: "
                f"index={index!r}, amino_acid={amino_acid!r}"
            )
        amino_acids[index] = amino_acid
    if amino_acids and amino_acids[-1] == "*":
        amino_acids.pop()
    return "".join(amino_acids)


def translation_matches_target(translated: str, target: str) -> bool:
    """Compare an audited CDS translation while treating genome N as unknown.

    An ``X`` generated from an ambiguous genomic codon is not evidence of a
    contradictory amino acid. A single trailing ``X`` is also accepted when
    an ambiguous terminal codon prevents recognizing the annotated stop.
    """

    def same_length_match(left, right):
        return len(left) == len(right) and all(
            observed == expected or observed == "X"
            for observed, expected in zip(left, right)
        )

    return same_length_match(translated, target) or (
        len(translated) == len(target) + 1
        and translated.endswith("X")
        and same_length_match(translated[:-1], target)
    )


def validate_job_name(value: str) -> str:
    value = str(value)
    if (
        not JOB_NAME_PATTERN.fullmatch(value)
        or value in {".", ".."}
        or ".." in value
    ):
        raise ValueError(
            "Job name must be 1-128 letters, digits, dot, underscore, or "
            "hyphen; path and shell syntax are forbidden"
        )
    return value


def protein_key(protein_id: str) -> str:
    """Return the gene key used by AlphaGEM homology matrices."""

    value = str(protein_id).split()[0]
    parts = value.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return value


def normalize_alphagem_target_id(protein_id: str) -> str:
    """Normalize only the two unambiguous target-ID forms used by AlphaGEM."""

    value = str(protein_id)
    if "|" not in value:
        if not value or "SerialNumber" in value:
            raise ValueError(
                f"Invalid bare AlphaGEM target-protein identifier: {protein_id!r}"
            )
        gene_id = value
    else:
        parts = value.split("|")
        if (
            len(parts) != 3
            or parts[0] != "sp"
            or not parts[1]
            or parts[2] != parts[1]
        ):
            raise ValueError(
                "AlphaGEM target-protein identifiers must be bare IDs or exact "
                f"'sp|ID|ID' values: {protein_id!r}"
            )
        gene_id = parts[1]
    return f"sp|{gene_id}|{gene_id}"


def nonhomolog_query_ids(
    query_ids: Iterable[str],
    homolog_gene_ids: Iterable[str],
) -> Set[str]:
    """Return full query IDs whose normalized gene keys are not homologs."""

    homolog_keys = {protein_key(value) for value in homolog_gene_ids}
    return {
        str(query_id)
        for query_id in query_ids
        if protein_key(query_id) not in homolog_keys
    }


def fasta_record_ids(fasta_path: Path) -> Set[str]:
    ids = set()
    with Path(fasta_path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    if not ids:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    return ids


def fasta_esm_labels(fasta_path: Path) -> Set[str]:
    """Return the exact labels produced by fair-esm's FASTA dataset loader."""

    labels = []
    with Path(fasta_path).open(encoding="utf-8") as handle:
        for line_index, line in enumerate(handle):
            if not line.startswith(">"):
                continue
            label = line[1:].strip()
            labels.append(label or f"seqnum{line_index:09d}")
    if not labels:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    if len(set(labels)) != len(labels):
        raise ValueError(f"Duplicate ESM FASTA labels in {fasta_path}")
    return set(labels)


def assert_partial_embedding_keys(
    embedding: Mapping,
    fasta_path: Path,
) -> int:
    """Reject partial-cache keys that are not exact fair-esm FASTA labels."""

    if not isinstance(embedding, Mapping):
        raise ValueError("Partial embedding cache is not a keyed mapping")
    non_string = [value for value in embedding.keys() if not isinstance(value, str)]
    if non_string:
        raise ValueError(
            "Partial embedding cache contains non-string labels: "
            + ", ".join(repr(value) for value in non_string[:20])
        )
    observed = set(embedding.keys())
    unexpected = observed - fasta_esm_labels(fasta_path)
    if unexpected:
        raise ValueError(
            "Partial embedding contains labels absent from FASTA: "
            + ", ".join(sorted(unexpected)[:20])
        )
    return len(observed)


def fasta_records(fasta_path: Path) -> List[dict]:
    """Read FASTA records without silently collapsing duplicate identifiers."""

    path = Path(fasta_path)
    records = []
    record_id = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    records.append(
                        {"id": record_id, "sequence": "".join(sequence)}
                    )
                record_id = line[1:].strip().split()[0]
                if not record_id:
                    raise ValueError(
                        f"Empty FASTA identifier in {path} at line {line_number}"
                    )
                sequence = []
            else:
                if record_id is None:
                    raise ValueError(
                        f"FASTA sequence precedes the first header in {path}"
                    )
                sequence.append(line)
    if record_id is not None:
        records.append({"id": record_id, "sequence": "".join(sequence)})
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def target_fasta_identity(fasta_path: Path) -> dict:
    """Return exact and canonical identities for one target proteome."""

    path = Path(fasta_path)
    records = fasta_records(path)
    normalized = {}
    for record in records:
        normalized_id = normalize_alphagem_target_id(record["id"])
        key = protein_key(normalized_id)
        if key in normalized:
            raise ValueError(
                "Duplicate normalized target-protein identifier in "
                f"{path}: {key}"
            )
        sequence = str(record["sequence"]).upper()
        if not sequence:
            raise ValueError(f"Empty target-protein sequence in {path}: {key}")
        normalized[key] = sequence
    normalized_ids = "".join(f"{key}\n" for key in sorted(normalized))
    canonical_proteome = "".join(
        f"{key}\t{normalized[key]}\n" for key in sorted(normalized)
    )
    return {
        "sha256": file_sha256(path),
        "bytes": path.stat().st_size,
        "record_count": len(records),
        "normalized_id_sha256": hashlib.sha256(
            normalized_ids.encode("utf-8")
        ).hexdigest(),
        "canonical_proteome_sha256": hashlib.sha256(
            canonical_proteome.encode("utf-8")
        ).hexdigest(),
    }


def validate_target_input_contract(
    contract_path: Path,
    fasta_path: Path,
) -> dict:
    """Pin a formal paper run to one exact, named target-proteome input."""

    contract_path = Path(contract_path)
    if not contract_path.is_file() or contract_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"Target input contract is missing or empty: {contract_path}"
        )
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("schema_version") != 1:
        raise ValueError(
            "Unsupported target input contract schema_version: "
            f"{contract.get('schema_version')!r}"
        )
    case_id = contract.get("case_id")
    if not isinstance(case_id, str) or not case_id.strip():
        raise ValueError("Target input contract requires a non-empty case_id")
    expected = contract.get("source_fasta")
    required_fields = {
        "sha256",
        "bytes",
        "record_count",
        "normalized_id_sha256",
        "canonical_proteome_sha256",
    }
    if not isinstance(expected, dict) or not required_fields.issubset(expected):
        raise ValueError(
            "Target input contract source_fasta is missing required fields: "
            + ", ".join(sorted(required_fields - set(expected or {})))
        )
    provenance = contract.get("provenance")
    if not isinstance(provenance, dict) or not provenance:
        raise ValueError(
            "Target input contract requires a non-empty provenance object"
        )
    observed = target_fasta_identity(fasta_path)
    mismatches = {
        key: {"expected": expected.get(key), "observed": observed[key]}
        for key in sorted(required_fields)
        if expected.get(key) != observed[key]
    }
    if mismatches:
        raise ValueError(
            "Target FASTA does not match its pinned input contract: "
            + json.dumps(mismatches, sort_keys=True)
        )
    return {
        "schema_version": 1,
        "case_id": case_id.strip(),
        "source_fasta": observed,
        "provenance": provenance,
        "contract_sha256": file_sha256(contract_path),
    }


def coordinate_translation_exceptions(contract: Mapping) -> list:
    """Validate and normalize codon reassignments pinned by a coordinate contract."""

    payload = contract.get("translation_exceptions", [])
    if not isinstance(payload, list):
        raise ValueError(
            "Gap-fill coordinate contract translation_exceptions must be a list"
        )
    normalized = []
    seen = set()
    for index, entry in enumerate(payload):
        if not isinstance(entry, dict):
            raise ValueError(
                f"Translation exception {index} must be an object"
            )
        required = {
            "target_protein_id",
            "refseq_cds_id",
            "amino_acid_index_0",
            "genome_seqid",
            "genome_start",
            "genome_end",
            "strand",
            "replacement_amino_acid",
            "source_feature",
        }
        if set(entry) != required:
            raise ValueError(
                "Translation exception fields do not match the formal schema: "
                f"index={index}"
            )
        target_id = protein_key(str(entry["target_protein_id"]).strip())
        refseq_cds_id = str(entry["refseq_cds_id"]).strip()
        genome_seqid = str(entry["genome_seqid"]).strip()
        strand = str(entry["strand"]).strip()
        replacement = str(entry["replacement_amino_acid"]).strip()
        source_feature = str(entry["source_feature"]).strip()
        try:
            amino_acid_index = int(entry["amino_acid_index_0"])
            genome_start = int(entry["genome_start"])
            genome_end = int(entry["genome_end"])
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Translation exception {index} has non-integer coordinates"
            ) from error
        if (
            not target_id
            or not refseq_cds_id
            or not genome_seqid
            or strand not in {"+", "-"}
            or amino_acid_index < 0
            or genome_start < 1
            or genome_end != genome_start + 2
            or not re.fullmatch(r"[A-Z]", replacement)
            or source_feature != "NCBI_GFF_transl_except"
        ):
            raise ValueError(
                f"Invalid formal translation exception at index {index}"
            )
        key = (target_id, amino_acid_index)
        if key in seen:
            raise ValueError(
                f"Duplicate translation exception for {target_id}:{amino_acid_index}"
            )
        seen.add(key)
        normalized.append(
            {
                "target_protein_id": target_id,
                "refseq_cds_id": refseq_cds_id,
                "amino_acid_index_0": amino_acid_index,
                "genome_seqid": genome_seqid,
                "genome_start": genome_start,
                "genome_end": genome_end,
                "strand": strand,
                "replacement_amino_acid": replacement,
                "source_feature": source_feature,
            }
        )
    return normalized


def load_coordinate_translation_exceptions(contract_path: Path) -> list:
    contract = json.loads(Path(contract_path).read_text(encoding="utf-8"))
    return coordinate_translation_exceptions(contract)


def validate_gapfill_coordinate_contract(
    contract_path: Path,
    *,
    fasta_path: Path,
    genome_path: Path,
    coordinate_map_path: Path,
    coordinate_segments_path: Path,
    expected_case_id: str,
) -> dict:
    """Bind coordinate artifacts to one auditable recovery procedure."""

    contract_path = Path(contract_path)
    if not contract_path.is_file() or contract_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"Gap-fill coordinate contract is missing or empty: {contract_path}"
        )
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    schema_version = contract.get("schema_version")
    if schema_version not in {1, 2}:
        raise ValueError(
            "Unsupported gap-fill coordinate contract schema_version: "
            f"{schema_version!r}"
        )
    translation_exceptions = coordinate_translation_exceptions(contract)
    if schema_version == 1 and translation_exceptions:
        raise ValueError(
            "Translation exceptions require coordinate contract schema_version 2"
        )
    if contract.get("case_id") != expected_case_id:
        raise ValueError(
            "Gap-fill coordinate contract case_id does not match the target "
            f"input contract: {contract.get('case_id')!r} != "
            f"{expected_case_id!r}"
        )
    expected_semantics = {
        "coordinate_system": "1-based-inclusive",
        "segment_order": "genomic_start_ascending",
        "hit_overlap": "cds_segment_union",
    }
    if contract.get("coordinate_semantics") != expected_semantics:
        raise ValueError(
            "Gap-fill coordinate contract has unsupported coordinate "
            "semantics"
        )
    supplied_paths = {
        "target_fasta": Path(fasta_path),
        "raw_genome": Path(genome_path),
        "coordinate_map": Path(coordinate_map_path),
        "coordinate_segments": Path(coordinate_segments_path),
    }
    recorded_inputs = contract.get("inputs")
    if not isinstance(recorded_inputs, dict):
        raise ValueError("Gap-fill coordinate contract requires inputs")
    observed_inputs = {}
    for logical_name, path in supplied_paths.items():
        observed = {
            "bytes": path.stat().st_size,
            "sha256": file_sha256(path),
        }
        expected = recorded_inputs.get(logical_name)
        if not isinstance(expected, dict) or any(
            expected.get(key) != value
            for key, value in observed.items()
        ):
            raise ValueError(
                "Gap-fill coordinate input does not match its contract: "
                f"{logical_name}"
            )
        observed_inputs[logical_name] = observed
    provenance = contract.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError(
            "Gap-fill coordinate contract requires provenance"
        )
    for field in (
        "audit_script_sha256",
        "resource_manifest_sha256",
    ):
        value = str(provenance.get(field, ""))
        if not re.fullmatch(r"[0-9a-f]{64}", value):
            raise ValueError(
                f"Invalid coordinate-contract provenance field: {field}"
            )
    audit_command = provenance.get("audit_command")
    official_sources = provenance.get("official_sources")
    if (
        not isinstance(audit_command, list)
        or not audit_command
        or not all(isinstance(value, str) and value for value in audit_command)
    ):
        raise ValueError(
            "Coordinate-contract audit_command must be a non-empty argv list"
        )
    if (
        not isinstance(official_sources, list)
        or not official_sources
        or not all(isinstance(value, str) and value for value in official_sources)
    ):
        raise ValueError(
            "Coordinate-contract official_sources must be a non-empty list"
        )
    if translation_exceptions:
        exception_provenance = provenance.get(
            "translation_exception_provenance"
        )
        if not isinstance(exception_provenance, dict):
            raise ValueError(
                "Translation exceptions require formal provenance"
            )
        for field in (
            "official_gff_sha256",
            "audit_json_sha256",
            "audit_script_sha256",
        ):
            if not re.fullmatch(
                r"[0-9a-f]{64}",
                str(exception_provenance.get(field, "")),
            ):
                raise ValueError(
                    "Invalid translation-exception provenance field: "
                    f"{field}"
                )
    return {
        "schema_version": schema_version,
        "case_id": expected_case_id,
        "coordinate_semantics": expected_semantics,
        "inputs": observed_inputs,
        "provenance": provenance,
        "translation_exceptions": translation_exceptions,
        "contract_sha256": file_sha256(contract_path),
    }


def validate_gapfill_cds(cds_fasta: Path, protein_fasta: Path) -> dict:
    """Validate the explicit nucleotide-CDS contract required by gap filling."""

    cds_records = fasta_records(cds_fasta)
    protein_records = fasta_records(protein_fasta)

    def indexed_keys(records, label):
        keys: Dict[str, str] = {}
        for record in records:
            key = protein_key(record["id"])
            if key in keys:
                raise ValueError(
                    f"Duplicate normalized {label} ID {key!r}: "
                    f"{keys[key]!r} and {record['id']!r}"
                )
            keys[key] = record["id"]
        return keys

    cds_keys = indexed_keys(cds_records, "CDS")
    protein_keys = indexed_keys(protein_records, "protein")
    missing = set(protein_keys) - set(cds_keys)
    extra = set(cds_keys) - set(protein_keys)
    if missing or extra:
        raise ValueError(
            "Gap-fill CDS IDs do not map one-to-one to target protein IDs: "
            f"missing={len(missing)}, extra={len(extra)}"
        )

    nucleotide_alphabet = set("ACGTRYSWKMBDHVN")
    for record in cds_records:
        sequence = record["sequence"].upper()
        invalid = sorted(set(sequence) - nucleotide_alphabet)
        if not sequence:
            raise ValueError(f"Empty CDS sequence for {record['id']!r}")
        if invalid:
            raise ValueError(
                f"Non-nucleotide characters in CDS {record['id']!r}: "
                + "".join(invalid)
            )
        if len(sequence) % 3:
            raise ValueError(
                f"CDS length is not divisible by three for {record['id']!r}"
            )
    return {
        "records": len(cds_records),
        "input_kind": "nucleotide_cds",
        "normalized_ids": sorted(cds_keys),
    }


GAPFILL_COORDINATE_COLUMNS = (
    "target_protein_id",
    "genome_seqid",
    "start",
    "end",
    "strand",
    "refseq_protein_id",
    "refseq_cds_id",
    "cds_segment_count",
    "cds_segments",
    "locus_tag",
    "source_database_ids",
    "mapping_evidence",
    "mapping_bucket",
    "local_sequence_sha256",
    "reference_sequence_sha256",
)

GAPFILL_SEGMENT_COLUMNS = (
    "target_protein_id",
    "refseq_protein_id",
    "refseq_cds_id",
    "genome_seqid",
    "strand",
    "segment_index_genomic_order",
    "segment_start",
    "segment_end",
    "phase",
    "locus_tag",
    "translation_table",
)


def validate_gapfill_genome(
    genome_fasta: Path,
    protein_fasta: Path,
    coordinate_map: Path,
    coordinate_segments: Path,
    translation_exceptions=None,
) -> dict:
    """Validate complete envelope and CDS-segment raw-genome coordinates."""

    genome_records = fasta_records(genome_fasta)
    protein_records = fasta_records(protein_fasta)
    genome_by_id = {}
    nucleotide_alphabet = set("ACGTRYSWKMBDHVN")
    for record in genome_records:
        record_id = record["id"]
        if record_id in genome_by_id:
            raise ValueError(f"Duplicate genome FASTA ID: {record_id}")
        sequence = record["sequence"].upper()
        invalid = sorted(set(sequence) - nucleotide_alphabet)
        if invalid:
            raise ValueError(
                f"Non-nucleotide characters in genome record {record_id!r}: "
                + "".join(invalid)
            )
        genome_by_id[record_id] = sequence

    protein_by_key = {}
    protein_sequence_hash_by_key = {}
    protein_sequence_by_key = {}
    for record in protein_records:
        key = protein_key(record["id"])
        if key in protein_by_key:
            raise ValueError(f"Duplicate normalized target protein ID: {key}")
        protein_by_key[key] = record["id"]
        protein_sequence_hash_by_key[key] = hashlib.sha256(
            str(record["sequence"]).upper().encode("utf-8")
        ).hexdigest()
        protein_sequence_by_key[key] = str(record["sequence"]).upper()

    normalized_translation_exceptions = coordinate_translation_exceptions(
        {"translation_exceptions": list(translation_exceptions or [])}
    )
    exceptions_by_target = {}
    for entry in normalized_translation_exceptions:
        target_key = entry["target_protein_id"]
        if target_key not in protein_by_key:
            raise ValueError(
                "Translation exception refers to an unknown target protein: "
                f"{target_key}"
            )
        exceptions_by_target.setdefault(target_key, []).append(entry)

    mapping_path = Path(coordinate_map)
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GAPFILL_COORDINATE_COLUMNS:
            raise ValueError(
                "Gap-fill coordinate-map header must exactly equal: "
                + "\t".join(GAPFILL_COORDINATE_COLUMNS)
            )
        rows = []
        mapped_keys = {}
        for line_number, row in enumerate(reader, start=2):
            target_id = str(row["target_protein_id"]).strip()
            target_key = protein_key(target_id)
            if not target_id or target_key not in protein_by_key:
                raise ValueError(
                    f"Unknown target protein at {mapping_path}:{line_number}: "
                    f"{target_id!r}"
                )
            if target_key in mapped_keys:
                raise ValueError(
                    f"Duplicate target protein coordinate mapping {target_key!r}: "
                    f"lines {mapped_keys[target_key]} and {line_number}"
                )
            seqid = str(row["genome_seqid"]).strip()
            if seqid not in genome_by_id:
                raise ValueError(
                    f"Unknown genome seqid at {mapping_path}:{line_number}: "
                    f"{seqid!r}"
                )
            try:
                start = int(row["start"])
                end = int(row["end"])
                segment_count = int(row["cds_segment_count"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Non-integer genomic coordinate at "
                    f"{mapping_path}:{line_number}"
                ) from error
            if not (1 <= start <= end <= len(genome_by_id[seqid])):
                raise ValueError(
                    f"Out-of-range genomic coordinate at "
                    f"{mapping_path}:{line_number}: {seqid}:{start}-{end}"
                )
            if segment_count < 1:
                raise ValueError(
                    f"Invalid cds_segment_count at {mapping_path}:{line_number}"
                )
            strand = str(row["strand"]).strip()
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid strand at {mapping_path}:{line_number}: {strand!r}"
                )
            for column in (
                "refseq_protein_id",
                "refseq_cds_id",
                "cds_segments",
                "locus_tag",
                "source_database_ids",
                "mapping_evidence",
                "mapping_bucket",
            ):
                if not str(row[column]).strip():
                    raise ValueError(
                        f"Empty {column} at {mapping_path}:{line_number}"
                    )
            local_hash = str(row["local_sequence_sha256"]).strip()
            reference_hash = str(
                row["reference_sequence_sha256"]
            ).strip()
            for label, value in (
                ("local_sequence_sha256", local_hash),
                ("reference_sequence_sha256", reference_hash),
            ):
                if not re.fullmatch(r"[0-9a-f]{64}", value):
                    raise ValueError(
                        f"Invalid {label} at {mapping_path}:{line_number}"
                    )
            if local_hash != reference_hash:
                raise ValueError(
                    "Target and reference protein sequence hashes differ at "
                    f"{mapping_path}:{line_number}"
                )
            if local_hash != protein_sequence_hash_by_key[target_key]:
                raise ValueError(
                    "Coordinate-map local sequence hash does not match the "
                    f"target FASTA for {target_key}"
                )
            mapped_keys[target_key] = line_number
            rows.append(
                {
                    **row,
                    "target_protein_id": target_key,
                    "genome_seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "cds_segment_count": segment_count,
                    "local_sequence_sha256": local_hash,
                    "reference_sequence_sha256": reference_hash,
                }
            )

    missing = set(protein_by_key) - set(mapped_keys)
    if missing:
        raise ValueError(
            "Gap-fill coordinate map is incomplete for the target proteome: "
            f"missing={len(missing)}; examples={sorted(missing)[:20]}"
        )
    envelope_by_target = {
        row["target_protein_id"]: row
        for row in rows
    }
    segments_path = Path(coordinate_segments)
    with segments_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GAPFILL_SEGMENT_COLUMNS:
            raise ValueError(
                "Gap-fill coordinate-segment header must exactly equal: "
                + "\t".join(GAPFILL_SEGMENT_COLUMNS)
            )
        segment_rows = []
        segments_by_target = {}
        seen_segments = set()
        for line_number, row in enumerate(reader, start=2):
            target_key = protein_key(
                str(row["target_protein_id"]).strip()
            )
            if target_key not in protein_by_key:
                raise ValueError(
                    "Unknown target protein in coordinate segments at "
                    f"{segments_path}:{line_number}: {target_key!r}"
                )
            seqid = str(row["genome_seqid"]).strip()
            strand = str(row["strand"]).strip()
            try:
                segment_index = int(row["segment_index_genomic_order"])
                start = int(row["segment_start"])
                end = int(row["segment_end"])
                phase = int(row["phase"])
                translation_table = int(row["translation_table"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "Non-integer coordinate-segment field at "
                    f"{segments_path}:{line_number}"
                ) from error
            if seqid not in genome_by_id:
                raise ValueError(
                    "Unknown genome seqid in coordinate segments at "
                    f"{segments_path}:{line_number}: {seqid!r}"
                )
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid segment strand at {segments_path}:{line_number}"
                )
            if (
                segment_index < 1
                or phase not in {0, 1, 2}
                or translation_table < 1
            ):
                raise ValueError(
                    "Invalid segment index or phase at "
                    f"{segments_path}:{line_number}"
                )
            if not (1 <= start <= end <= len(genome_by_id[seqid])):
                raise ValueError(
                    "Out-of-range coordinate segment at "
                    f"{segments_path}:{line_number}: {seqid}:{start}-{end}"
                )
            for column in (
                "refseq_protein_id",
                "refseq_cds_id",
                "locus_tag",
            ):
                if not str(row[column]).strip():
                    raise ValueError(
                        f"Empty {column} at {segments_path}:{line_number}"
                    )
            envelope = envelope_by_target[target_key]
            for column, segment_value, envelope_value in (
                ("genome_seqid", seqid, envelope["genome_seqid"]),
                ("strand", strand, envelope["strand"]),
                (
                    "refseq_protein_id",
                    str(row["refseq_protein_id"]).strip(),
                    str(envelope["refseq_protein_id"]).strip(),
                ),
                (
                    "refseq_cds_id",
                    str(row["refseq_cds_id"]).strip(),
                    str(envelope["refseq_cds_id"]).strip(),
                ),
            ):
                if segment_value != envelope_value:
                    raise ValueError(
                        "Coordinate segment disagrees with its envelope for "
                        f"{target_key}: {column}={segment_value!r}, "
                        f"envelope={envelope_value!r}"
                    )
            segment_key = (
                target_key,
                segment_index,
                seqid,
                start,
                end,
                strand,
            )
            if segment_key in seen_segments:
                raise ValueError(
                    "Duplicate coordinate segment at "
                    f"{segments_path}:{line_number}: {segment_key}"
                )
            seen_segments.add(segment_key)
            normalized = {
                **row,
                "target_protein_id": target_key,
                "genome_seqid": seqid,
                "strand": strand,
                "segment_index_genomic_order": segment_index,
                "segment_start": start,
                "segment_end": end,
                "start": start,
                "end": end,
                "phase": phase,
                "translation_table": translation_table,
                "mapping_evidence": envelope["mapping_evidence"],
            }
            segment_rows.append(normalized)
            segments_by_target.setdefault(target_key, []).append(normalized)

    missing_segments = set(protein_by_key) - set(segments_by_target)
    extra_segments = set(segments_by_target) - set(protein_by_key)
    if missing_segments or extra_segments:
        raise ValueError(
            "Gap-fill coordinate segments do not cover the target proteome: "
            f"missing={len(missing_segments)}, extra={len(extra_segments)}"
        )
    recovery_translation_mismatches = []
    applied_translation_exceptions = []
    for target_key, target_segments in segments_by_target.items():
        envelope = envelope_by_target[target_key]
        ordered = sorted(
            target_segments,
            key=lambda row: row["segment_index_genomic_order"],
        )
        expected_indices = list(range(1, len(ordered) + 1))
        observed_indices = [
            row["segment_index_genomic_order"] for row in ordered
        ]
        if observed_indices != expected_indices:
            raise ValueError(
                f"Non-contiguous segment indices for {target_key}: "
                f"{observed_indices}"
            )
        if len(ordered) != envelope["cds_segment_count"]:
            raise ValueError(
                f"Coordinate segment count disagrees with envelope for "
                f"{target_key}"
            )
        starts = [row["segment_start"] for row in ordered]
        if starts != sorted(starts):
            raise ValueError(
                "Coordinate segments are not in genomic-start order for "
                f"{target_key}"
            )
        for left, right in zip(ordered, ordered[1:]):
            if left["segment_end"] >= right["segment_start"]:
                raise ValueError(
                    f"Overlapping coordinate segments for {target_key}"
                )
        serialized_segments = ";".join(
            f"{row['segment_start']}-{row['segment_end']}"
            for row in ordered
        )
        if serialized_segments != str(envelope["cds_segments"]).strip():
            raise ValueError(
                f"Coordinate segment list disagrees with envelope for "
                f"{target_key}"
            )
        locus_tags = {
            str(row["locus_tag"]).strip() for row in ordered
        }
        if locus_tags != {str(envelope["locus_tag"]).strip()}:
            raise ValueError(
                f"Coordinate segment locus_tag disagrees with envelope for "
                f"{target_key}"
            )
        translation_tables = {
            row["translation_table"] for row in ordered
        }
        if len(translation_tables) != 1:
            raise ValueError(
                f"Coordinate segments use multiple translation tables for "
                f"{target_key}"
            )
        if (
            ordered[0]["segment_start"] != envelope["start"]
            or ordered[-1]["segment_end"] != envelope["end"]
        ):
            raise ValueError(
                f"Coordinate segment union does not match envelope for "
                f"{target_key}"
            )
        transcript_order = (
            ordered
            if envelope["strand"] == "+"
            else list(reversed(ordered))
        )
        oriented_parts = []
        oriented_positions = []
        for segment in transcript_order:
            genomic_sequence = genome_by_id[segment["genome_seqid"]][
                segment["segment_start"] - 1 : segment["segment_end"]
            ]
            genomic_positions = list(
                range(segment["segment_start"], segment["segment_end"] + 1)
            )
            if segment["strand"] == "-":
                genomic_sequence = reverse_complement(genomic_sequence)
                genomic_positions.reverse()
            oriented_parts.append(genomic_sequence)
            oriented_positions.extend(genomic_positions)
        coding_sequence = "".join(oriented_parts)
        initial_phase = transcript_order[0]["phase"]
        if initial_phase:
            coding_sequence = coding_sequence[initial_phase:]
            oriented_positions = oriented_positions[initial_phase:]
        envelope_mapping_bucket_len = str(envelope.get("mapping_bucket", "")).strip()
        is_recovery_entry_len = envelope_mapping_bucket_len.startswith("recovery_")
        if not coding_sequence or (len(coding_sequence) % 3 and not is_recovery_entry_len):
            raise ValueError(
                f"Reconstructed CDS length is invalid for {target_key}: "
                f"{len(coding_sequence)}"
            )
        if len(coding_sequence) % 3 and is_recovery_entry_len:
            # Pad recovery CDS to a multiple of 3 for translation
            coding_sequence = coding_sequence + "N" * (3 - len(coding_sequence) % 3)
            oriented_positions.extend([None] * (3 - len(oriented_positions) % 3))
        translation_table = next(iter(translation_tables))
        target_protein = protein_sequence_by_key[target_key]
        amino_acid_overrides = {}
        for exception in exceptions_by_target.get(target_key, []):
            amino_acid_index = exception["amino_acid_index_0"]
            nucleotide_offset = amino_acid_index * 3
            codon = coding_sequence[nucleotide_offset : nucleotide_offset + 3]
            codon_positions = oriented_positions[
                nucleotide_offset : nucleotide_offset + 3
            ]
            if (
                exception["refseq_cds_id"]
                != str(envelope["refseq_cds_id"]).strip()
                or exception["genome_seqid"] != envelope["genome_seqid"]
                or exception["strand"] != envelope["strand"]
                or len(codon) != 3
                or None in codon_positions
                or min(codon_positions) != exception["genome_start"]
                or max(codon_positions) != exception["genome_end"]
            ):
                raise ValueError(
                    "Translation exception does not bind to the audited CDS "
                    f"coordinates for {target_key}:{amino_acid_index}"
                )
            code = {
                **STANDARD_CODON_TABLE,
                **GENETIC_CODE_OVERRIDES[translation_table],
            }
            if code.get(codon, "X") != "*":
                raise ValueError(
                    "Translation exception is not bound to a raw stop codon "
                    f"for {target_key}:{amino_acid_index}: codon={codon}"
                )
            replacement = exception["replacement_amino_acid"]
            if (
                amino_acid_index >= len(target_protein)
                or target_protein[amino_acid_index] != replacement
            ):
                raise ValueError(
                    "Translation exception replacement does not match the "
                    f"target protein for {target_key}:{amino_acid_index}"
                )
            amino_acid_overrides[amino_acid_index] = replacement
            applied_translation_exceptions.append(
                {
                    **exception,
                    "raw_codon": codon,
                    "genome_positions_in_translation_order": codon_positions,
                }
            )
        translated = translate_nucleotide_sequence(
            coding_sequence,
            translation_table=translation_table,
            expected_starts_with_methionine=(
                target_protein.startswith("M")
            ),
            amino_acid_overrides=amino_acid_overrides,
        )
        if not translation_matches_target(translated, target_protein):
            mapping_bucket = str(
                envelope.get("mapping_bucket", "")
            ).strip()
            if not mapping_bucket.startswith("recovery_"):
                raise ValueError(
                    "Genome CDS-segment translation does not match an "
                    "explicit target-protein mapping for "
                    f"{target_key}: mapping_bucket={mapping_bucket!r}, "
                    f"table={translation_table}, "
                    f"translated_length={len(translated)}, "
                    f"target_length={len(target_protein)}"
                )
            mismatch_record = {
                "target_protein_id": target_key,
                "mapping_bucket": mapping_bucket,
                "translation_table": translation_table,
                "translated_length": len(translated),
                "target_length": len(target_protein),
            }
            recovery_translation_mismatches.append(mismatch_record)
            import warnings
            warnings.warn(
                "Recovered genome CDS-segment translation does not match the "
                f"target protein for {target_key}: table={translation_table}, "
                f"translated_length={len(translated)}, "
                f"target_length={len(target_protein)}",
                stacklevel=2,
            )
    return {
        "input_kind": "raw_genome_with_envelope_and_segment_maps",
        "genome_records": len(genome_records),
        "target_proteins": len(protein_by_key),
        "coordinate_records": len(rows),
        "coordinate_segment_records": len(segment_rows),
        "recovery_translation_mismatches": recovery_translation_mismatches,
        "translation_exceptions_applied": applied_translation_exceptions,
        "rows": rows,
        "segment_rows": segment_rows,
    }


def map_tblastn_hits_to_target_proteins(
    rows,
    coordinate_rows,
    *,
    threshold,
    operator,
    tie_break="strict",
    tie_audit=None,
    unmapped_hit_policy="error",
    unmapped_audit=None,
):
    """Map each query's highest qualifying genomic hit to one target protein.

    tie_break="strict" (default) refuses queries whose equal-highest rows map
    to different target proteins. tie_break="deterministic" instead chooses
    the candidate with the largest coordinate overlap in base pairs, then the
    lexicographically smallest target_protein_id on a remaining tie; every
    broken tie is recorded both on the selected record ("tie_break_audit")
    and, when a list is passed as tie_audit, aggregated there.

    unmapped_hit_policy="error" (default, author behavior) raises when an
    equal-top row overlaps no predicted target coordinate on the same strand.
    unmapped_hit_policy="skip_with_audit" skips only those rows instead and
    records each skipped row in unmapped_audit; a query whose equal-top rows
    all lack mapping then carries no gap evidence (it simply does not appear
    in the returned selection, which is the nonfindgenes path). A single HSP
    overlapping multiple target proteins stays fatal in every mode.
    """
    if tie_break not in {"strict", "deterministic"}:
        raise ValueError(f"Unsupported gapfill hit tie_break: {tie_break!r}")
    if unmapped_hit_policy not in {"error", "skip_with_audit"}:
        raise ValueError(
            f"Unsupported unmapped_hit_policy: {unmapped_hit_policy!r}"
        )
    coordinates_by_seqid = {}
    for coordinate in coordinate_rows:
        coordinates_by_seqid.setdefault(
            str(coordinate["genome_seqid"]),
            [],
        ).append(coordinate)

    qualifying_by_query = {}
    for row in rows:
        if passes_threshold(row["bitscore"], threshold, operator):
            qualifying_by_query.setdefault(str(row["query"]), []).append(row)

    selected = {}
    for query, query_rows in qualifying_by_query.items():
        ranked = sorted(
            query_rows,
            key=lambda row: (
                -float(row["bitscore"]),
                float(row["evalue"]),
                str(row["subject"]),
                min(int(row["sstart"]), int(row["send"])),
                max(int(row["sstart"]), int(row["send"])),
            ),
        )
        best_score = float(ranked[0]["bitscore"])
        best_rows = [
            row for row in ranked if float(row["bitscore"]) == best_score
        ]
        mapped_best_rows = []
        for row in best_rows:
            hit_start = min(int(row["sstart"]), int(row["send"]))
            hit_end = max(int(row["sstart"]), int(row["send"]))
            hit_strand = "+" if int(row["sstart"]) <= int(row["send"]) else "-"
            overlap_by_target = {}
            coordinate_by_target = {}
            for coordinate in coordinates_by_seqid.get(str(row["subject"]), []):
                if coordinate["strand"] != hit_strand:
                    continue
                overlap = max(
                    0,
                    min(hit_end, int(coordinate["end"]))
                    - max(hit_start, int(coordinate["start"]))
                    + 1,
                )
                if overlap:
                    target_id = str(coordinate["target_protein_id"])
                    overlap_by_target[target_id] = (
                        overlap_by_target.get(target_id, 0) + overlap
                    )
                    coordinate_by_target.setdefault(target_id, coordinate)
            row_candidates = [
                (
                    overlap,
                    target_id,
                    row,
                    coordinate_by_target[target_id],
                )
                for target_id, overlap in sorted(overlap_by_target.items())
            ]
            if not row_candidates:
                if unmapped_hit_policy == "skip_with_audit":
                    if unmapped_audit is not None:
                        unmapped_audit.append({
                            "query": query,
                            "subject": str(row["subject"]),
                            "sstart": int(hit_start),
                            "send": int(hit_end),
                            "strand": hit_strand,
                            "bitscore": best_score,
                            "evalue": float(row["evalue"]),
                            "reason": "no_overlapping_target_coordinate",
                        })
                    continue
                raise ValueError(
                    f"Highest tblastn hit for {query!r} does not map to a "
                    "target protein coordinate; every equal-top row must map "
                    f"(bitscore={best_score}, subject={row['subject']!r}, "
                    f"start={hit_start}, end={hit_end}, strand={hit_strand})"
                )
            row_target_ids = {value[1] for value in row_candidates}
            if len(row_target_ids) != 1:
                raise ValueError(
                    f"Ambiguous equal-highest tblastn hit for {query!r}: "
                    "one genomic HSP overlaps multiple target proteins: "
                    + ", ".join(sorted(row_target_ids))
                )
            row_max_overlap = max(value[0] for value in row_candidates)
            row_finalists = [
                value
                for value in row_candidates
                if value[0] == row_max_overlap
            ]
            mapped_best_rows.append(
                sorted(
                    row_finalists,
                    key=lambda value: (
                        value[1],
                        float(value[2]["evalue"]),
                        str(value[2]["subject"]),
                        min(int(value[2]["sstart"]), int(value[2]["send"])),
                        max(int(value[2]["sstart"]), int(value[2]["send"])),
                    ),
                )[0]
            )
        if not mapped_best_rows:
            # Only reachable under skip_with_audit: every equal-top row of
            # this query lacked an overlapping target coordinate, so the
            # query simply carries no gap evidence (the nonfindgenes path).
            if unmapped_audit is not None:
                unmapped_audit.append({
                    "query": query,
                    "subject": None,
                    "sstart": None,
                    "send": None,
                    "strand": None,
                    "bitscore": None,
                    "evalue": None,
                    "reason": "all_equal_top_hits_unmapped",
                })
            continue
        target_ids = {value[1] for value in mapped_best_rows}
        if len(target_ids) != 1:
            if tie_break == "strict":
                raise ValueError(
                    f"Ambiguous highest tblastn hit for {query!r}: "
                    + ", ".join(sorted(target_ids))
                )
            # Deterministic tie-break: largest coordinate overlap in bp,
            # then the lexicographically smallest target_protein_id. The
            # final row pick among the chosen target's finalists follows the
            # existing deterministic ordering (evalue, subject, coordinates).
            best_overlap = max(value[0] for value in mapped_best_rows)
            finalists = [
                value
                for value in mapped_best_rows
                if value[0] == best_overlap
            ]
            chosen_target = min(value[1] for value in finalists)
            chosen_rows = [
                value for value in finalists if value[1] == chosen_target
            ]
            selected_overlap, target_id, row, coordinate = sorted(
                chosen_rows,
                key=lambda value: (
                    float(value[2]["evalue"]),
                    str(value[2]["subject"]),
                    min(int(value[2]["sstart"]), int(value[2]["send"])),
                    max(int(value[2]["sstart"]), int(value[2]["send"])),
                ),
            )[0]
            audit_record = {
                "query": query,
                "tied_targets": sorted(target_ids),
                "chosen": chosen_target,
                "tie_break": "deterministic",
                "overlap_bp": int(selected_overlap),
            }
            if tie_audit is not None:
                tie_audit.append(audit_record)
        else:
            selected_overlap, target_id, row, coordinate = sorted(
                mapped_best_rows,
                key=lambda value: (
                    value[1],
                    float(value[2]["evalue"]),
                    str(value[2]["subject"]),
                    min(int(value[2]["sstart"]), int(value[2]["send"])),
                    max(int(value[2]["sstart"]), int(value[2]["send"])),
                ),
            )[0]
            audit_record = None
        selected[query] = {
            **dict(row),
            "genome_subject": str(row["subject"]),
            "target_protein_id": target_id,
            "subject": target_id,
            "coordinate_evidence": str(coordinate["mapping_evidence"]),
            "coordinate_overlap_bp": int(selected_overlap),
        }
        if audit_record is not None:
            selected[query]["tie_break_audit"] = audit_record
    return selected


def passes_threshold(value: float, threshold: float, operator: str) -> bool:
    if operator == "gt":
        return float(value) > float(threshold)
    if operator == "ge":
        return float(value) >= float(threshold)
    raise ValueError(f"Unsupported threshold operator: {operator}")


def best_gapfill_hits(rows, threshold: float, operator: str):
    """Select one deterministic highest-bit-score subject for each query."""

    selected = {}
    for row in rows:
        bitscore = float(row["bitscore"])
        if not passes_threshold(bitscore, threshold, operator):
            continue
        query = str(row["query"])
        candidate = (
            -bitscore,
            float(row["evalue"]),
            str(row["subject"]),
        )
        if query not in selected or candidate < selected[query][0]:
            selected[query] = (candidate, dict(row))
    return {query: value[1] for query, value in selected.items()}


def assert_embedding_matches_fasta(
    embedding_path: Path,
    fasta_path: Path,
) -> int:
    """Reject a cache whose keys differ from fair-esm's exact FASTA labels."""

    with Path(embedding_path).open("rb") as handle:
        embedding = pickle.load(handle)
    if not isinstance(embedding, Mapping):
        raise ValueError(f"Embedding cache is not a keyed mapping: {embedding_path}")
    non_string = [value for value in embedding.keys() if not isinstance(value, str)]
    if non_string:
        raise ValueError(
            "Embedding cache contains non-string labels: "
            + ", ".join(repr(value) for value in non_string[:20])
        )
    embedding_ids = set(embedding.keys())
    fasta_labels = fasta_esm_labels(fasta_path)
    if embedding_ids != fasta_labels:
        missing = fasta_labels - embedding_ids
        extra = embedding_ids - fasta_labels
        raise ValueError(
            "Embedding cache does not match ESM FASTA labels: "
            f"missing={len(missing)}, extra={len(extra)}, "
            f"embedding={embedding_path}, fasta={fasta_path}"
        )
    return len(fasta_labels)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_prott5_xl_uniref50_snapshot(snapshot_path: Path) -> dict:
    """Reject incomplete or wrong local ProtT5-XL-UniRef50 snapshots."""

    snapshot = Path(snapshot_path)
    if snapshot.is_symlink():
        raise ValueError(
            f"ProtTrans snapshot root must not be a symlink: {snapshot}"
        )
    if not snapshot.is_dir():
        raise NotADirectoryError(f"ProtTrans snapshot is not a directory: {snapshot}")
    snapshot_resolved = snapshot.resolve()
    symlink_entries = [
        str(path)
        for path in snapshot.rglob("*")
        if path.is_symlink()
    ]
    if symlink_entries:
        raise ValueError(
            "ProtTrans snapshot must not contain symlinks: "
            + ", ".join(symlink_entries[:20])
        )

    def require_contained_regular_file(path, label):
        if path.is_symlink():
            raise ValueError(
                f"ProtTrans {label} must not be a symlink: {path}"
            )
        try:
            path.resolve().relative_to(snapshot_resolved)
        except ValueError as error:
            raise ValueError(
                f"ProtTrans {label} escapes the snapshot directory: {path}"
            ) from error
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(
                f"ProtTrans snapshot is missing a non-empty {label}: {snapshot}"
            )

    config_path = snapshot / "config.json"
    tokenizer_path = snapshot / "spiece.model"
    for required in (config_path, tokenizer_path):
        require_contained_regular_file(required, required.name)

    config = json.loads(config_path.read_text(encoding="utf-8"))
    expected_config = {
        "model_type": "t5",
        "d_model": 1024,
        "d_ff": 16384,
        "num_layers": 24,
        "num_heads": 32,
    }
    mismatches = {
        key: {"expected": expected, "observed": config.get(key)}
        for key, expected in expected_config.items()
        if config.get(key) != expected
    }
    if mismatches:
        raise ValueError(
            "ProtTrans snapshot is not Rostlab/prot_t5_xl_uniref50: "
            + json.dumps(mismatches, sort_keys=True)
        )

    index_candidates = (
        snapshot / "model.safetensors.index.json",
        snapshot / "pytorch_model.bin.index.json",
    )
    index_path = next((path for path in index_candidates if path.is_file()), None)
    if index_path is not None:
        require_contained_regular_file(index_path, index_path.name)
        index = json.loads(index_path.read_text(encoding="utf-8"))
        weight_map = index.get("weight_map")
        if not isinstance(weight_map, dict) or not weight_map:
            raise ValueError(f"Empty or invalid weight_map in {index_path}")
        weight_names = sorted({str(value) for value in weight_map.values()})
        unsafe = [
            name
            for name in weight_names
            if Path(name).name != name or name in {".", ".."}
        ]
        if unsafe:
            raise ValueError(
                "Unsafe ProtTrans weight shard names: "
                + ", ".join(unsafe[:20])
            )
        weight_files = [snapshot / name for name in weight_names]
    else:
        weight_files = [
            path
            for path in (
                snapshot / "model.safetensors",
                snapshot / "pytorch_model.bin",
            )
            if path.is_file()
        ]
    if not weight_files:
        raise FileNotFoundError(
            f"ProtTrans snapshot contains no PyTorch or safetensors weights: {snapshot}"
        )
    symlinked = [str(path) for path in weight_files if path.is_symlink()]
    if symlinked:
        raise ValueError(
            "ProtTrans weight shards must not be symlinks: "
            + ", ".join(symlinked[:20])
        )
    missing_or_empty = [
        str(path)
        for path in weight_files
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing_or_empty:
        raise FileNotFoundError(
            "ProtTrans snapshot has missing or empty weight shards: "
            + ", ".join(missing_or_empty[:20])
        )
    for path in weight_files:
        require_contained_regular_file(path, path.name)
    return {
        "model_identifier": "Rostlab/prot_t5_xl_uniref50",
        "config": expected_config,
        "tokenizer_file": tokenizer_path.name,
        "weight_files": [path.name for path in weight_files],
    }


def assert_embedding_provenance(
    embedding_path: Path,
    fasta_path: Path,
    metadata_path: Path,
    *,
    model_identifier: str,
    model_path: Path = None,
    pooling: str = "mean_residue_excluding_bos_eos",
) -> dict:
    """Validate the cache bytes, sequences, model and embedding semantics."""
    embedding_path = Path(embedding_path)
    fasta_path = Path(fasta_path)
    metadata_path = Path(metadata_path)
    if not metadata_path.is_file():
        raise ValueError(
            f"Embedding provenance sidecar is missing: {metadata_path}"
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected = {
        "schema_version": 1,
        "complete": True,
        "fasta_sha256": file_sha256(fasta_path),
        "embedding_sha256": file_sha256(embedding_path),
        "model_identifier": Path(model_identifier).stem,
        "representation_layer": 33,
        "pooling": pooling,
        "max_residues": 1022,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"Embedding provenance mismatch for {key}: "
                f"expected={value!r}, observed={metadata.get(key)!r}"
            )
    if model_path is not None and Path(model_path).is_file():
        model_hash = file_sha256(model_path)
        if metadata.get("model_sha256") != model_hash:
            raise ValueError(
                "Embedding provenance does not match the ESM model weights"
            )
    assert_embedding_matches_fasta(embedding_path, fasta_path)
    return metadata


def write_embedding_provenance(
    embedding_path: Path,
    fasta_path: Path,
    metadata_path: Path,
    *,
    model_identifier: str,
    model_path: Path = None,
    token_budget: int,
    generator_path: Path = None,
    pooling: str = "mean_residue_excluding_bos_eos",
) -> dict:
    payload = {
        "schema_version": 1,
        "complete": True,
        "fasta_sha256": file_sha256(fasta_path),
        "embedding_sha256": file_sha256(embedding_path),
        "model_identifier": Path(model_identifier).stem,
        "model_sha256": (
            file_sha256(model_path)
            if model_path is not None and Path(model_path).is_file()
            else None
        ),
        "representation_layer": 33,
        "pooling": pooling,
        "max_residues": 1022,
        "token_budget": int(token_budget),
        "generator_sha256": (
            file_sha256(generator_path)
            if generator_path is not None and Path(generator_path).is_file()
            else None
        ),
    }
    Path(metadata_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return payload


def passes_usalign_filter(
    tm_score_1: float,
    tm_score_2: float,
    coverage: float,
    plddt_1: float,
    plddt_2: float,
    *,
    tm_threshold: float,
    coverage_threshold: float,
    plddt_threshold: float,
) -> bool:
    return (
        (float(tm_score_1) + float(tm_score_2)) / 2 >= tm_threshold
        and float(coverage) >= coverage_threshold
        and float(plddt_1) >= plddt_threshold
        and float(plddt_2) >= plddt_threshold
    )


def filter_blast_alignments(rows, coverage, identity):
    """Apply the inclusive cut-offs stated in the AlphaGEM Methods."""
    return [
        row
        for row in rows
        if float(row["COV"]) >= coverage and float(row["PID"]) >= identity
    ]


def reciprocal_best_hits(forward_rows, reverse_rows):
    """Return the released implementation's PID-first reciprocal best hits.

    After coverage filtering, the author code uses ``PID.idxmax()`` for each
    query. Therefore the highest percent identity wins and an equal-PID tie is
    resolved by the first row in the DIAMOND output. This preserves the
    released selection semantics; bit score and lexical subject order are not
    additional tie-breakers.
    """

    def best_by_query(rows):
        grouped = {}
        for row in rows:
            grouped.setdefault(row["gene"], []).append(row)
        best = {}
        for query, candidates in grouped.items():
            maximum_identity = max(
                float(row["PID"]) for row in candidates
            )
            best[query] = next(
                row
                for row in candidates
                if float(row["PID"]) == maximum_identity
            )
        return best

    forward_best = best_by_query(forward_rows)
    reverse_best = best_by_query(reverse_rows)
    result = []
    for query, forward in forward_best.items():
        row = dict(forward)
        reverse = reverse_best.get(row["subject"])
        if reverse is None:
            continue
        row["BBH"] = "<=>" if reverse["subject"] == query else "->"
        result.append(row)
    return result


def reaction_is_mass_and_charge_balanced(imbalance):
    """COBRA reports both elemental and charge imbalance in this mapping."""
    return not imbalance


def gapfill_media_plan(strategy, growth_medium, supports_full):
    if growth_medium not in {"min", "full"}:
        raise ValueError(f"Unsupported growth medium: {growth_medium}")
    if strategy == "selected_medium":
        return (growth_medium,)
    if strategy == "legacy_dual_medium":
        return ("min", "full") if supports_full else ("min",)
    raise ValueError(f"Unsupported gap-fill strategy: {strategy}")


def select_by_pca_distance(data, dataindex, distance):
    """Apply the published 2-D PCA nearest-distance rule."""
    import numpy as np

    values = np.asarray(data, dtype=float)
    labels = np.asarray(dataindex)
    if len(values) != len(labels):
        raise ValueError("ProtTrans data and labels must have the same length")
    if len(values) <= 1:
        return labels
    centered = values - values.mean(axis=0, keepdims=True)
    _, _, components = np.linalg.svd(centered, full_matrices=False)
    component_count = min(2, len(values), values.shape[1])
    reduced = centered @ components[:component_count].T
    distances = np.linalg.norm(reduced - reduced[0], axis=1)
    nearest = float(np.min(distances[1:]))
    return labels[distances <= distance * nearest]


def resolve_eggnog_diamond_db(
    refname: str,
    data_dir: Path,
    *,
    mode: str = "paper_taxon",
) -> Path:
    """Resolve an explicit paper or legacy eggNOG DIAMOND database."""

    paper_databases = {
        "yeast": "fungi.dmnd",
        "ecoli": "bacteria.dmnd",
        "strco": "bacteria.dmnd",
        "synechocystis": "bacteria.dmnd",
        "human": "mammalia.dmnd",
    }
    packaged_databases = {
        "yeast": "yeast.dmnd",
        "ecoli": "bacteria.dmnd",
        "strco": "bacteria.dmnd",
        "synechocystis": "bacteria.dmnd",
        "human": "mammalia.dmnd",
    }
    if mode == "paper_taxon":
        database_names = paper_databases
    elif mode == "packaged_fallback":
        database_names = packaged_databases
    else:
        raise ValueError(f"Unsupported eggNOG database mode: {mode}")
    try:
        database = Path(data_dir) / database_names[refname]
    except KeyError as exc:
        raise ValueError(f"Unsupported eggNOG reference: {refname}") from exc
    if not database.is_file():
        raise FileNotFoundError(
            f"Required eggNOG DIAMOND database not found: {database}"
        )
    if mode == "paper_taxon" and refname == "yeast":
        packaged_fallback = Path(data_dir) / packaged_databases["yeast"]
        if (
            packaged_fallback.is_file()
            and file_sha256(database) == file_sha256(packaged_fallback)
        ):
            raise ValueError(
                "The required fungi.dmnd is byte-identical to the packaged "
                "yeast.dmnd fallback; this is not an independent paper-taxonomy "
                "database. Use mode='packaged_fallback' and disclose the "
                "substitution, or provide a verified fungi.dmnd."
            )
    return database
