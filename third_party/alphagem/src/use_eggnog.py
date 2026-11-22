import subprocess
import sys
from pathlib import Path

import pandas as pd

from pipeline_contracts import protein_key, resolve_eggnog_diamond_db


def _nonhomolog_genes(name):
    target = pd.read_excel(f"working/{name}/{name}.xlsx")
    homologs = pd.read_excel(f"working/{name}/matrix_homolog{name}.xlsx")
    homolog_gene_ids = {str(value) for value in homologs.iloc[:, 2].dropna()}
    return [
        str(value)
        for value in target.iloc[:, 0].dropna()
        if str(value) not in homolog_gene_ids
    ]


def _parse_annotations(annotation_path, nonhomolog_genes):
    annotations = pd.read_csv(
        annotation_path,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
    )
    by_gene = {}
    for _, row in annotations.iterrows():
        if len(row) <= 14:
            continue
        by_gene[protein_key(row.iloc[0])] = (row.iloc[10], row.iloc[14])

    rows = []
    for gene in nonhomolog_genes:
        ec, reaction = by_gene.get(gene, ("-", "-"))
        rows.append(
            {
                "gene": gene,
                "EC": ec if pd.notna(ec) and ec else "-",
                "Rec": reaction if pd.notna(reaction) and reaction else "-",
            }
        )
    return pd.DataFrame(rows)


def eggnog(name, refname, i=1, database_mode="paper_taxon"):
    """Annotate the non-homolog FASTA using the database shipped by AlphaGEM."""

    root = Path.cwd()
    working = root / "working" / name
    input_fasta = working / f"{name}_homoleft.fasta"
    if not input_fasta.is_file():
        raise FileNotFoundError(
            f"Non-homolog FASTA not found: {input_fasta}. "
            "Run the annotation preparation step first."
        )

    emapper = root / "tools" / "eggnog_mapper" / "emapper.py"
    data_dir = root / "tools" / "eggnog_mapper" / "data"
    database = resolve_eggnog_diamond_db(
        refname,
        data_dir,
        mode=database_mode,
    )
    if not emapper.is_file():
        raise FileNotFoundError(f"eggNOG-mapper entry point not found: {emapper}")

    temp_dir = root / "temp" / name
    temp_dir.mkdir(parents=True, exist_ok=True)
    output_name = f"test{name}"
    tax_scope = {
        "yeast": "Fungi",
        "human": "Mammalia",
        "ecoli": "Bacteria",
        "strco": "Actinobacteria",
        "synechocystis": "Bacteria",
    }[refname]
    command = [
        sys.executable,
        str(emapper),
        "-m",
        "diamond",
        "-i",
        str(input_fasta),
        "-o",
        output_name,
        "--output_dir",
        str(working),
        "--data_dir",
        str(data_dir),
        "--tax_scope",
        tax_scope,
        "--cpu",
        "0",
        "--dmnd_db",
        str(database),
        "--temp_dir",
        str(temp_dir),
        "--override",
    ]
    subprocess.run(command, check=True)

    annotation_path = working / f"{output_name}.emapper.annotations"
    if not annotation_path.is_file():
        raise RuntimeError(
            f"eggNOG-mapper completed without its annotation file: {annotation_path}"
        )
    result = _parse_annotations(annotation_path, _nonhomolog_genes(name))
    result.to_excel(working / f"eggec2{name}.xlsx")
