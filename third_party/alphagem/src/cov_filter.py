"""Sequence-alignment filtering used by the PLMSearch route.

The public repository silently continued after failed BLAST commands and reused
bundled database sidecars that are incompatible with recent BLAST releases.
This implementation creates both databases inside the current run directory,
fails on command errors, and applies the inclusive thresholds from the paper.
"""

import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from tqdm import tqdm

from src.pipeline_contracts import (
    filter_blast_alignments,
    reciprocal_best_hits,
)


BLAST_COLUMNS = [
    "gene",
    "subject",
    "PID",
    "alnLength",
    "mismatchCount",
    "gapOpenCount",
    "queryStart",
    "queryEnd",
    "subjectStart",
    "subjectEnd",
    "eVal",
    "bitScore",
]

PAIRWISE_RESULT_COLUMNS = [
    "Gene1",
    "Gene2",
    "E-value",
    "Score",
    "Identity",
    "Alignment Length",
    "COV1",
    "COV2",
]


def _protein_id(value):
    parts = str(value).split("|")
    return parts[1] if len(parts) >= 3 else parts[0]


def _load_fasta(path):
    return {
        _protein_id(record.id): record.seq
        for record in SeqIO.parse(str(path), "fasta")
    }


def _validate_candidate_workbook(
    candidates,
    gene_dict,
    reference_sequences,
    target_sequences,
):
    """Validate the author workbook without treating its Excel index as data."""

    missing_columns = {0, 1}.difference(candidates.columns)
    if missing_columns:
        raise ValueError(
            "PLMSearch candidate workbook must contain numeric columns 0 "
            "(reference locus) and 1 (target entry)"
        )
    if candidates[[0, 1]].isna().any(axis=None):
        raise ValueError("PLMSearch candidate workbook contains empty gene IDs")

    reference_loci = set(candidates[0])
    target_genes = set(candidates[1])
    unknown_loci = reference_loci.difference(gene_dict)
    if unknown_loci:
        raise KeyError(
            "Candidate reference loci are missing from the reference table: "
            f"{sorted(unknown_loci)[:5]}"
        )
    reference_genes = {gene_dict[locus] for locus in reference_loci}
    missing_reference_sequences = reference_genes.difference(reference_sequences)
    if missing_reference_sequences:
        raise KeyError(
            "Candidate reference sequences are missing: "
            f"{sorted(missing_reference_sequences)[:5]}"
        )
    missing_target_sequences = target_genes.difference(target_sequences)
    if missing_target_sequences:
        raise KeyError(
            "Candidate target sequences are missing: "
            f"{sorted(missing_target_sequences)[:5]}"
        )


def _candidate_gene_pair(row):
    """Read the author workbook's numeric gene columns by label, not position."""

    return row.loc[0], row.loc[1]


def result_blast(name, refname, threshold):
    reference = pd.read_excel(f"./data_available/{refname}.xlsx")
    gene_dict = reference.set_index("Gene Names (ordered locus)")["Entry"].to_dict()
    candidates = pd.read_excel(
        f"./working/{name}/{name}_ss_homolog{threshold}.xlsx"
    )
    reference_sequences = _load_fasta(f"./data_available/{refname}.fasta")
    target_sequences = _load_fasta(f"./working/{name}/{name}.fasta")
    _validate_candidate_workbook(
        candidates,
        gene_dict,
        reference_sequences,
        target_sequences,
    )
    workdir = Path("working") / name
    results = []

    for _, row in candidates.iterrows():
        reference_locus, target_gene = _candidate_gene_pair(row)
        reference_gene = gene_dict[reference_locus]
        if reference_gene not in reference_sequences:
            raise KeyError(f"Reference sequence is missing: {reference_gene}")
        if target_gene not in target_sequences:
            raise KeyError(f"Target sequence is missing: {target_gene}")
        reference_temp = workdir / "temp1.fasta"
        target_temp = workdir / "temp2.fasta"
        reference_temp.write_text(
            f">{reference_gene}\n{reference_sequences[reference_gene]}\n",
            encoding="utf-8",
        )
        target_temp.write_text(
            f">{target_gene}\n{target_sequences[target_gene]}\n",
            encoding="utf-8",
        )
        output = workdir / "temp_blast_results.txt"
        blastp = NcbiblastpCommandline(
            query=str(reference_temp),
            subject=str(target_temp),
            outfmt=6,
            out=str(output),
        )
        blastp()
        best = None
        with output.open(encoding="utf-8") as handle:
            for line in tqdm(handle, leave=False):
                fields = line.rstrip().split("\t")
                evalue = float(fields[10])
                if best is None or evalue < best["E-value"]:
                    best = {
                        "Gene1": reference_locus,
                        "Gene2": target_gene,
                        "E-value": evalue,
                        "Score": float(fields[11]),
                        "Identity": float(fields[2]),
                        "Alignment Length": int(fields[3]),
                    }
        if best is not None:
            results.append(best)

    table = pd.DataFrame(results, columns=PAIRWISE_RESULT_COLUMNS)
    target_lengths = get_gene_lens(f"{name}.fasta", f"./working/{name}")
    reference_lengths = get_gene_lens(
        f"{refname}.fasta", "./data_available"
    )
    for index, row in table.iterrows():
        reference_length = reference_lengths[gene_dict[row["Gene1"]]]
        target_length = target_lengths[row["Gene2"]]
        table.loc[index, "COV1"] = row["Alignment Length"] / reference_length
        table.loc[index, "COV2"] = row["Alignment Length"] / target_length
    table.to_csv(
        workdir / f"{name}_ss_{threshold}_blast_results.csv",
        index=False,
    )
    return table


def get_gene_lens(query, in_folder="prots"):
    return {
        _protein_id(record.id): len(record.seq)
        for record in SeqIO.parse(str(Path(in_folder) / query), "fasta")
    }


def get_gene_lens_frame(query, in_folder="prots"):
    return pd.DataFrame(
        [
            {"gene": gene, "gene_length": length}
            for gene, length in get_gene_lens(query, in_folder).items()
        ]
    )


def make_blast_db(fasta, db_prefix, db_type="prot"):
    fasta = Path(fasta)
    db_prefix = Path(db_prefix)
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "makeblastdb",
        "-in",
        str(fasta),
        "-dbtype",
        db_type,
        "-out",
        str(db_prefix),
        "-parse_seqids",
    ]
    print("running:", " ".join(command))
    subprocess.run(command, check=True)
    suffix = ".pin" if db_type == "prot" else ".nin"
    expected = Path(f"{db_prefix}{suffix}")
    if not expected.is_file():
        raise RuntimeError(f"makeblastdb did not create {expected}")
    return db_prefix


def make_diamond_db(fasta, db_prefix):
    fasta = Path(fasta)
    db_prefix = Path(db_prefix)
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "diamond",
        "makedb",
        "--in",
        str(fasta),
        "--db",
        str(db_prefix),
    ]
    print("running:", " ".join(command))
    subprocess.run(command, check=True)
    expected = Path(f"{db_prefix}.dmnd")
    if not expected.is_file():
        raise RuntimeError(f"diamond makedb did not create {expected}")
    return db_prefix


def run_blastp(seq, db, out, outfmt=6, evalue=0.001, threads=1):
    command = [
        "blastp",
        "-db",
        str(db),
        "-query",
        str(seq),
        "-out",
        str(out),
        "-evalue",
        str(evalue),
        "-outfmt",
        str(outfmt),
        "-num_threads",
        str(threads),
    ]
    print("running:", " ".join(command))
    subprocess.run(command, check=True)
    output = Path(out)
    if not output.is_file():
        raise RuntimeError(f"blastp did not create {output}")
    return output


def run_diamond_blastp(seq, db, out, evalue=0.001, threads=1):
    command = [
        "diamond",
        "blastp",
        "--db",
        str(db),
        "--query",
        str(seq),
        "--out",
        str(out),
        "--evalue",
        str(evalue),
        "--outfmt",
        "6",
        "--threads",
        str(threads),
    ]
    print("running:", " ".join(command))
    subprocess.run(command, check=True)
    output = Path(out)
    if not output.is_file():
        raise RuntimeError(f"diamond blastp did not create {output}")
    return output


def _read_blast_table(path, lengths):
    path = Path(path)
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS + ["gene_length", "COV"])
    table = pd.read_csv(path, sep="\t", names=BLAST_COLUMNS)
    table["gene"] = table["gene"].map(_protein_id)
    table["subject"] = table["subject"].map(_protein_id)
    table = pd.merge(
        table,
        lengths,
        on="gene",
        how="left",
        validate="many_to_one",
    )
    if table["gene_length"].isna().any():
        missing = sorted(table.loc[table["gene_length"].isna(), "gene"].unique())
        raise ValueError(f"BLAST query length missing for: {missing[:10]}")
    table["COV"] = table["alnLength"] / table["gene_length"]
    return table


def get_bbh(
    query,
    subject,
    in_fold1,
    in_fold2,
    cov_thre,
    backend="diamond",
    threads=1,
):
    workdir = Path(in_fold1)
    query_fasta = workdir / query
    subject_fasta = Path(in_fold2) / subject
    db_dir = workdir / "blastdb"
    if backend == "diamond":
        make_db = make_diamond_db
        run_search = run_diamond_blastp
    elif backend == "blastp":
        make_db = make_blast_db
        run_search = run_blastp
    else:
        raise ValueError(f"Unsupported BBH backend: {backend}")
    query_db = make_db(query_fasta, db_dir / "target")
    subject_db = make_db(subject_fasta, db_dir / "reference")
    forward_path = run_search(
        query_fasta,
        subject_db,
        workdir / "target_bench.txt",
        threads=threads,
    )
    reverse_path = run_search(
        subject_fasta,
        query_db,
        workdir / "bench_target.txt",
        threads=threads,
    )
    forward = _read_blast_table(
        forward_path, get_gene_lens_frame(query, in_fold1)
    )
    reverse = _read_blast_table(
        reverse_path, get_gene_lens_frame(subject, in_fold2)
    )
    forward.to_csv(workdir / "bbh_cov.csv", index=False)
    reverse.to_csv(workdir / "bbh2_cov.csv", index=False)
    return forward, reverse


def bbh_analysis(query, subject, in_fold1, in_fold2, cov_thre):
    del query, subject, in_fold2
    forward = pd.read_csv(Path(in_fold1) / "bbh_cov.csv")
    reverse = pd.read_csv(Path(in_fold1) / "bbh2_cov.csv")
    forward_rows = filter_blast_alignments(
        forward.to_dict("records"), coverage=cov_thre, identity=0
    )
    reverse_rows = filter_blast_alignments(
        reverse.to_dict("records"), coverage=cov_thre, identity=0
    )
    rows = reciprocal_best_hits(forward_rows, reverse_rows)
    output = pd.DataFrame(rows, columns=list(forward.columns) + ["BBH"])
    output.to_csv(
        Path(in_fold1) / f"bbh_parsed_{cov_thre}.csv", index=False
    )
    return output


def bbh(
    name,
    refname,
    threshold,
    cov_thre=0.25,
    pid_thre=50,
    direction="bidirect",
    backend="diamond",
    threads=1,
):
    reference = pd.read_excel(f"./data_available/{refname}.xlsx")
    entry_to_gene = reference.set_index("Entry")[
        "Gene Names (ordered locus)"
    ].to_dict()
    query = f"{name}.fasta"
    subject = f"{refname}.fasta"
    target_dir = f"./working/{name}"
    reference_dir = "./data_available"
    get_bbh(
        query,
        subject,
        target_dir,
        reference_dir,
        cov_thre,
        backend=backend,
        threads=threads,
    )
    table = bbh_analysis(
        query, subject, target_dir, reference_dir, cov_thre
    )
    if direction == "bidirect":
        table = table[table["BBH"] == "<=>"]
    elif direction != "direct":
        raise ValueError(f"Unsupported BLAST direction: {direction}")
    if not table.empty:
        table = table[table["PID"] >= pid_thre].copy()
        table["subject"] = table["subject"].replace(entry_to_gene)
        table.rename(
            columns={"subject": "Gene1", "gene": "Gene2"},
            inplace=True,
        )
    table.to_excel(
        f"./working/{name}/{name}_ss{threshold}_bbh_pid{pid_thre}"
        f"_cov{cov_thre}_{direction}.xlsx",
        index=False,
    )


def cov_filter(name, refname, threshold, cov_threshold, id_threshold):
    del refname
    path = f"./working/{name}/{name}_ss_{threshold}_blast_results.csv"
    table = pd.read_csv(path)
    table = table[
        (table["COV1"] >= cov_threshold)
        & (table["COV2"] >= cov_threshold)
        & (table["Identity"] >= id_threshold)
    ]
    table.to_csv(
        f"./working/{name}/{name}_ss_{threshold}_blast_cov"
        f"{cov_threshold}_id{id_threshold}.csv",
        index=False,
    )
    return table


def merge(
    name,
    refname,
    threshold,
    cov_threshold,
    id_threshold,
    cov_thre,
    pid_thre,
    direction,
):
    del refname
    bbh_table = pd.read_excel(
        f"./working/{name}/{name}_ss{threshold}_bbh_pid{pid_thre}"
        f"_cov{cov_thre}_{direction}.xlsx"
    )
    ss_table = pd.read_csv(
        f"./working/{name}/{name}_ss_{threshold}_blast_cov"
        f"{cov_threshold}_id{id_threshold}.csv"
    )
    merged = (
        pd.merge(bbh_table, ss_table, on=["Gene1", "Gene2"], how="outer")
        if not bbh_table.empty
        else ss_table
    )
    merged.to_excel(
        f"./working/{name}/{name}_ss{threshold}_cov{cov_threshold}"
        f"_id{id_threshold}_bbh_pid{pid_thre}_cov{cov_thre}"
        f"_{direction}.xlsx",
        index=False,
    )
