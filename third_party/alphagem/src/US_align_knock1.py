import pandas as pd
import hashlib
import json
import re
import os
import tempfile
from plddt_find import get_plddt_from_pdb as gpl
from functools import partial
import multiprocessing
from tqdm import tqdm
import subprocess


def _file_record(path):
    path = os.fspath(path)
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": os.path.getsize(path), "sha256": digest.hexdigest()}


def _write_json_atomic(path, payload):
    path = os.path.abspath(os.fspath(path))
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=os.path.dirname(path),
            prefix=".{}.".format(os.path.basename(path)),
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = handle.name
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


def _pair_frame(row_index, metrics):
    values = [
        gx.iat[row_index, 1],
        gx.iat[row_index, 2],
        *metrics,
    ]
    return pd.DataFrame([values], columns=[1, 2, 3, 4, 5, 6, 7])


def _failure(
    row_index,
    resolved_reference_id,
    target_id,
    reason,
    **details,
):
    return {
        "frame": _pair_frame(row_index, [None, None, None, None, None]),
        "failure": {
            **details,
            "reason": reason,
            "reference_id": str(gx.iat[row_index, 1]),
            "resolved_reference_id": str(resolved_reference_id),
            "row_index": int(row_index),
            "target_id": str(target_id),
        },
    }


def parse_usalign_metrics(output):
    tm_scores = re.findall(r"TM-score=\s*([\d\.]+)", output)
    aligned_match = re.search(r"Aligned length=\s*(\d+)", output)
    length_matches = re.findall(
        r"TM-score=\s*[\d\.]+\s+\(normalized by length of Structure_[12]: L=(\d+),",
        output
    )

    if len(tm_scores) < 2:
        raise ValueError("US-align output does not contain two TM-scores")

    tm_score_1 = float(tm_scores[0])
    tm_score_2 = float(tm_scores[1])

    if aligned_match and len(length_matches) >= 2:
        aligned_length = int(aligned_match.group(1))
        len_1 = int(length_matches[0])
        len_2 = int(length_matches[1])
        avg_coverage = round(((aligned_length / len_1) + (aligned_length / len_2)) / 2, 4)
    else:
        raise ValueError("US-align output does not contain alignment coverage")


    return tm_score_1, tm_score_2, avg_coverage


def worker(i, gx, yeadict, name, refname, path):
    if gx.iat[i, 1]=='b448':
        print(duiying1(gx.iat[i, 1], yeadict))
    return tdblast(
        duiying1(gx.iat[i, 1], yeadict),
        gx.iat[i, 2],
        i,
        name,
        refname,
        path,
    )


def parallel_processing(
    gx,
    yeadict,
    name,
    refname,
    path,
):
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    worker_with_args = partial(
        worker,
        gx=gx,
        yeadict=yeadict,
        name=name,
        refname=refname,
        path=path,
    )
    results = list(tqdm(pool.imap(worker_with_args, range(len(gx.index))), total=len(gx.index)))
    pool.close()
    pool.join()

    return results

def read_gx(name):
    global gx
    gx = pd.read_excel(f'working/{name}/matrix_orthofinder{name}.xlsx')

def tdblast(
    cmd1,
    cmd2,
    i,
    name,
    refname,
    path='',
):
    pathwd=os.getcwd()
    if path=='':
        path=f'{pathwd}/struct_data/taryeast/{name}'
    if cmd1==0 or cmd2==0:
        raise RuntimeError(
            "US-align failed for pair at row {}: "
            "missing_identifier".format(i)
        )

    usalign_exec = os.path.join(pathwd, "tools/USalign/USalign")
    input_file_1 = os.path.join(path, f"AF-{cmd2}-F1-model_v4.pdb")
    input_file_2 = os.path.join(pathwd, f"struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb")

    missing_structures = []
    if not os.path.isfile(input_file_2):
        missing_structures.append("reference")
    if not os.path.isfile(input_file_1):
        missing_structures.append("target")
    if missing_structures:
        return _failure(
            i,
            cmd1,
            cmd2,
            "missing_structure",
            missing_structures=missing_structures,
        )

    cmd = [usalign_exec, input_file_1, input_file_2]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "US-align failed for pair at row {}: "
            "nonzero exit {}".format(i, error.returncode)
        ) from error
    try:
        d, e, avg_cov = parse_usalign_metrics(res.stdout)
    except ValueError as error:
        raise RuntimeError(
            "US-align failed for pair at row {}: {}".format(
                i,
                error,
            )
        ) from error

    # print(f"Score 1: {d}, Score 2: {e}")
    # print('USalign done', cmd1, cmd2)
    try:
        f=gpl(f'{path}/AF-{cmd2}-F1-model_v4.pdb')
        g=gpl(f'{pathwd}/struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb')
    except (OSError, ValueError) as error:
        raise RuntimeError(
            "US-align failed for pair at row {}: pLDDT "
            "read failure".format(i)
        ) from error
    if f is None or g is None:
        raise RuntimeError(
            "US-align failed for pair at row {}: "
            "missing pLDDT".format(i)
        )
    return {
        "frame": _pair_frame(i, [d, e, g, f, avg_cov]),
        "failure": None,
    }

def duiying1(a,yea):
    return yea.get(a, a)
def duiying2(a,tar):
    for i in range(len(tar.index)):
        if a==tar.iat[i,1]:
            return tar.iat[i,0]
    return 0

def US_align_find(name,path,refname,path2, *, write_audit=False):
    yea = pd.read_excel(f'data_available/{refname}.xlsx')
    yeadict={yea.iat[i,2]:yea.iat[i,0] for i in range(len(yea.index))}
    read_gx(name)
    gx2 = pd.DataFrame()
    tar = pd.read_excel(f'working/{name}/{name}.xlsx')
    results=[]
    print('start usalign alignment')
    results = parallel_processing(
        gx,
        yeadict,
        name,
        refname,
        path,
    )
    gx2=pd.concat([result["frame"] for result in results])
    gx2.index=range(len(gx2.index))
    print(gx2)
    matrix_path = f'working/{name}/matrix_USalign_filtered{name}.xlsx'
    gx2.to_excel(matrix_path)
    if write_audit:
        failures = [
            result["failure"]
            for result in results
            if result["failure"] is not None
        ]
        input_path = f"working/{name}/matrix_orthofinder{name}.xlsx"
        _write_json_atomic(
            f"working/{name}/usalign_execution.json",
            {
                "failed_pair_count": len(failures),
                "failures": failures,
                "matrix_output": _file_record(matrix_path),
                "orthofinder_input": _file_record(input_path),
                "successful_pair_count": len(results) - len(failures),
                "total_pair_count": len(results),
            },
        )
