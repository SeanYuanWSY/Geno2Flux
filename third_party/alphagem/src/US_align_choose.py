import pandas as pd
import numpy as np
from pipeline_contracts import passes_usalign_filter


def US_align_choose(
    name,
    filterscore=0.5,
    coverage_threshold=0.8,
    plddt_threshold=70,
):
    gx3 = pd.DataFrame()
    gx2 = pd.read_excel(f'./working/{name}/matrix_USalign_filtered{name}.xlsx')
    required = gx2.iloc[:, [3, 4, 5, 6, 7]].apply(
        pd.to_numeric,
        errors="coerce",
    )
    valid_rows = required.notna().all(axis=1)
    rejected = int((~valid_rows).sum())
    if rejected:
        print(f"Rejecting {rejected} US-align rows with missing metrics")
    gx2 = gx2.loc[valid_rows].copy()
    a = 0
    for i in range(len(gx2.index)):
        if passes_usalign_filter(
            gx2.iat[i, 3],
            gx2.iat[i, 4],
            gx2.iat[i, 7],
            gx2.iat[i, 5],
            gx2.iat[i, 6],
            tm_threshold=filterscore,
            coverage_threshold=coverage_threshold,
            plddt_threshold=plddt_threshold,
        ):
            gx3 = pd.concat([gx3, pd.DataFrame({0: [gx2.iat[i, 1]],
                                                1: [gx2.iat[i, 2]],
                                                2: [gx2.iat[i, 3]],
                                                3: [gx2.iat[i, 4]],
                                                4: [gx2.iat[i, 5]],
                                                5: [gx2.iat[i, 6]]
                                                })])
            a += 1
    gx3.to_excel(f'working/{name}/matrix_homo_part1{name}.xlsx')
