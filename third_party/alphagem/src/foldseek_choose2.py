import pandas as pd
import subprocess
import tempfile
from pathlib import Path


TRANSPORTER_IDENTITY_THRESHOLD_PERCENT = 30.0


def sequence_identity_percent(sequence1, sequence2, aligner=None):
    """Return exact-residue identity over the full global-alignment length."""
    seq1 = str(getattr(sequence1, "seq", sequence1))
    seq2 = str(getattr(sequence2, "seq", sequence2))
    if aligner is None:
        from Bio import Align

        aligner = Align.PairwiseAligner()

    alignment = aligner.align(seq1, seq2)[0]
    alignment_length = int(alignment.length)
    if alignment_length == 0:
        raise ValueError("Cannot calculate sequence identity for an empty alignment")

    matches = 0
    for block1, block2 in zip(alignment.aligned[0], alignment.aligned[1]):
        start1, end1 = (int(value) for value in block1)
        start2, end2 = (int(value) for value in block2)
        aligned1 = seq1[start1:end1]
        aligned2 = seq2[start2:end2]
        matches += sum(
            residue1 == residue2
            for residue1, residue2 in zip(aligned1, aligned2)
        )
    return 100.0 * matches / alignment_length


def blastp_identity_percent(sequence1, sequence2, executable="blastp"):
    """Return the best-HSP percent identity reported by paper-version BLASTP."""
    seq1 = str(getattr(sequence1, "seq", sequence1))
    seq2 = str(getattr(sequence2, "seq", sequence2))
    if not seq1 or not seq2:
        raise ValueError("Cannot run BLASTP on an empty sequence")
    with tempfile.TemporaryDirectory(prefix="alphagem-transporter-") as tmp:
        temp_dir = Path(tmp)
        query = temp_dir / "query.fasta"
        subject = temp_dir / "subject.fasta"
        output = temp_dir / "blast.tsv"
        query.write_text(f">query\n{seq1}\n", encoding="utf-8")
        subject.write_text(f">subject\n{seq2}\n", encoding="utf-8")
        subprocess.run(
            [
                executable,
                "-query",
                str(query),
                "-subject",
                str(subject),
                "-out",
                str(output),
                "-outfmt",
                "6 pident bitscore",
                "-max_hsps",
                "1",
            ],
            check=True,
        )
        if not output.is_file() or output.stat().st_size == 0:
            return 0.0
        hits = []
        for line in output.read_text(encoding="utf-8").splitlines():
            identity, bit_score = line.split("\t")
            hits.append((float(bit_score), float(identity)))
        return max(hits)[1]


def passes_transporter_filter(
    coverage_sum,
    reference_plddt,
    target_plddt,
    tm_score,
    identity_percent,
    coverage_threshold,
    plddt_threshold,
    tm_score_threshold,
):
    return (
        coverage_sum >= 2 * coverage_threshold
        and reference_plddt >= plddt_threshold
        and target_plddt >= plddt_threshold
        and tm_score >= tm_score_threshold
        and identity_percent >= TRANSPORTER_IDENTITY_THRESHOLD_PERCENT
    )


def _passes_nontransporter_filter(
    coverage_sum,
    reference_plddt,
    target_plddt,
    tm_score,
    coverage_threshold,
    plddt_threshold,
    tm_score_threshold,
):
    return (
        coverage_sum >= 2 * coverage_threshold
        and reference_plddt >= plddt_threshold
        and target_plddt >= plddt_threshold
        and tm_score >= tm_score_threshold
    )


def _result_row(gx, index, transporter):
    return {
        "gene1": gx.iat[index, 1],
        "gene2": gx.iat[index, 2],
        "coverage": (gx.iat[index, 5] + gx.iat[index, 6]) / 2,
        "tmscore": gx.iat[index, 4],
        "sumscore": (
            (gx.iat[index, 5] + gx.iat[index, 6]) / 2
            + gx.iat[index, 4]
        ),
        "transporter": transporter,
    }


def foldseek_choose(
    refname,
    refmodel,
    name,
    TMscoretrans,
    coveragetrans,
    TMscore,
    coverage,
    pLDDt,
):
    # These author modules keep reference and target FASTA/model data in module
    # globals. Import lazily so identity helpers remain independently testable.
    import align2
    import model_transpoter_found

    gx = pd.read_excel(
        "working/{}/matrix_foldseek_filtered1_{}.xlsx".format(name, name)
    )
    model_transpoter_found.datapre(refmodel,refname)
    align2.predate(name,refname)
    selected_rows = []
    for i in range(len(gx.index)):
        reference_id = gx.iat[i, 1]
        target_id = gx.iat[i, 2]
        coverage_sum = gx.iat[i, 5] + gx.iat[i, 6]
        reference_plddt = gx.iat[i, 7]
        target_plddt = gx.iat[i, 8]
        tm_score = gx.iat[i, 4]

        if model_transpoter_found.transpoter(reference_id) == 1:
            if reference_id not in align2.refgenes:
                raise KeyError(
                    "Reference FASTA sequence not found for transporter {!r}".format(
                        reference_id
                    )
                )
            if target_id not in align2.tgenes:
                raise KeyError(
                    "Target FASTA sequence not found for transporter pair {!r}".format(
                        target_id
                    )
                )
            identity_percent = blastp_identity_percent(
                align2.refgenes[reference_id],
                align2.tgenes[target_id],
            )
            if passes_transporter_filter(
                coverage_sum,
                reference_plddt,
                target_plddt,
                tm_score,
                identity_percent,
                coveragetrans,
                pLDDt,
                TMscoretrans,
            ):
                selected_rows.append(_result_row(gx, i, True))
        elif _passes_nontransporter_filter(
            coverage_sum,
            reference_plddt,
            target_plddt,
            tm_score,
            coverage,
            pLDDt,
            TMscore,
        ):
            selected_rows.append(_result_row(gx, i, False))

    gx2 = pd.DataFrame(
        selected_rows,
        columns=[
            "gene1",
            "gene2",
            "coverage",
            "tmscore",
            "sumscore",
            "transporter",
        ],
    )
    reference_table = pd.read_excel("data_available/{}.xlsx".format(refname))
    entry_to_model_gene = {
        reference_table.iat[i, 0]: reference_table.iat[i, 2]
        for i in range(len(reference_table.index))
    }
    for i in range(len(gx2.index)):
        reference_entry = gx2.iat[i, 0]
        if reference_entry not in entry_to_model_gene:
            raise KeyError(
                "Reference entry {!r} is absent from data_available/{}.xlsx".format(
                    reference_entry,
                    refname,
                )
            )
        gx2.iat[i, 0] = entry_to_model_gene[reference_entry]
    gx2.to_excel(
        "./working/{}/matrix_foldseek_filtered2_{}.xlsx".format(name, name)
    )


def foldseek_choose2(
    refname,
    refmodel,
    name,
    TMscoretrans=0.7,
    coveragetrans=0.8,
    TMscore=0.7,
    coverage=0.8,
    pLDDt=70,
):
    foldseek_choose(
        refname,
        refmodel,
        name,
        TMscoretrans,
        coveragetrans,
        TMscore,
        coverage,
        pLDDt,
    )
