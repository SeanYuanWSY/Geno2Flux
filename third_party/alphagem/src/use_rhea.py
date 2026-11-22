import os
import pandas as pd
import shutil
import numpy as np
from Bio import SeqIO

from plmsearch import rhea_search
from pipeline_contracts import nonhomolog_query_ids
from ss_predict_choose import ss_predict_choose_rhea

def process_rhea(name,refname,threshold,nonhomolog_queries):
    rhea_homo = pd.read_excel(
        f"./working/{name}/{name}_rhea_homolog{str(threshold)}.xlsx"
    )
    if rhea_homo.empty:
        print("No rhea homologs found")
        rhea_homo = pd.DataFrame({"ID": [], "Query": []})
    else:
        rhea_homo = rhea_homo.rename(
            columns={
                rhea_homo.columns[-2]: "ID",
                rhea_homo.columns[-1]: "Query",
            }
        )[["ID", "Query"]]
    rhea_homo = rhea_homo[rhea_homo["Query"].isin(nonhomolog_queries)]
    gene_rhea = pd.read_csv(
        "./rhea/rhea2uniprot_sprot.tsv",
        sep="\t",
        dtype={"ID": str, "RHEA_ID": str},
    )[["ID", "RHEA_ID"]].dropna().drop_duplicates()
    rhea_homo["ID"] = rhea_homo["ID"].astype(str)
    rhea_homo = rhea_homo.merge(
        gene_rhea,
        on="ID",
        how="inner",
        validate="many_to_many",
    ).rename(columns={"RHEA_ID": "Rhea_id"})
    rhea_homo = rhea_homo.drop_duplicates()
    rhea_homo.to_excel(
        f"./working/{name}/{name}_rheaid_plmsearch.xlsx",
        index=False,
    )
def rhea(name,refname,threshold=0.9):
    homo = pd.read_excel(f'./working/{name}/matrix_homolog{name}.xlsx')
    query_ids = [
        record.id
        for record in SeqIO.parse(f"./working/{name}/{name}.fasta", "fasta")
    ]
    homolog_gene_ids = homo.iloc[:, 2].dropna().astype(str)
    nonhomolog_queries = nonhomolog_query_ids(query_ids, homolog_gene_ids)

    rhea_search(name,refname,threshold)
    ss_predict_choose_rhea(name,threshold)
    process_rhea(name,refname,threshold,nonhomolog_queries)
