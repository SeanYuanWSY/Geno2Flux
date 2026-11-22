import re
from pathlib import Path
from multiprocessing import Pool

import pandas as pd
from plddt_find import get_plddt_from_pdb as gpl
from tqdm import tqdm


ALPHAFOLD_PDB = re.compile(r"^AF-(.+)-F1-model(?:_v\d+)?\.pdb$", re.IGNORECASE)


def _protein_id(filename):
    match = ALPHAFOLD_PDB.match(filename)
    if not match:
        raise ValueError(
            "Structure filename must use AlphaFold naming "
            "'AF-<Entry>-F1-model[_vN].pdb'; got {!r}".format(filename)
        )
    return match.group(1)


def process_file(path, filename):
    pdb_path = Path(path) / filename
    return _protein_id(filename), gpl(str(pdb_path))


def _pdb_files(path, label):
    directory = Path(path).expanduser().resolve()
    if not directory.is_dir():
        raise NotADirectoryError(
            "{} structure directory not found: {}".format(label, directory)
        )
    files = sorted(
        item.name for item in directory.iterdir() if item.suffix.lower() == ".pdb"
    )
    if not files:
        raise FileNotFoundError(
            "{} structure directory contains no PDB files: {}".format(
                label,
                directory,
            )
        )
    for filename in files:
        _protein_id(filename)
    return directory, files


def _plddt_dict(path, files, description):
    with Pool() as pool:
        results = list(
            tqdm(
                pool.starmap(process_file, [(str(path), filename) for filename in files]),
                desc=description,
            )
        )
    values = dict(results)
    invalid = sorted(key for key, value in values.items() if value is None)
    if invalid:
        raise ValueError(
            "Could not read pLDDT values for: {}".format(", ".join(invalid))
        )
    return values


def generate_dict(path, refname, pathwd):
    target_path, target_files = _pdb_files(path, "Target")
    reference_path, reference_files = _pdb_files(
        Path(pathwd) / "struct_data" / refname,
        "Reference",
    )
    return (
        _plddt_dict(target_path, target_files, "generate target pLDDT"),
        _plddt_dict(reference_path, reference_files, "generate reference pLDDT"),
    )


def read_gx(name):
    return pd.read_excel(
        "./working/{}/matrix_foldseek_{}.xlsx".format(name, name)
    )


def plddt_datacreate(name, refname, path=''):
    pathwd = Path.cwd()
    if not path:
        path = pathwd / "struct_data" / "taryeast" / name
    return generate_dict(path, refname, pathwd)


def tdblast(cmd1, cmd2, i, gx, tardict, refdict):
    if cmd2 not in tardict:
        raise KeyError(
            "Target structure {!r} is absent from the explicit staged "
            "structure directory".format(cmd2)
        )
    if cmd1 not in refdict:
        raise KeyError(
            "Reference structure {!r} is absent from struct_data".format(cmd1)
        )
    return {
        0: gx.iat[i, 1],
        1: gx.iat[i, 2],
        2: gx.iat[i, 3],
        3: gx.iat[i, 4],
        4: gx.iat[i, 5],
        5: gx.iat[i, 6],
        6: refdict[cmd1],
        7: tardict[cmd2],
    }


def foldseek_choose(name, refname, path=''):
    gx = read_gx(name)
    tardict, refdict = plddt_datacreate(name, refname, path)
    rows = []
    for i in tqdm(range(len(gx.index)), "find pLDDT"):
        rows.append(
            tdblast(
                gx.iat[i, 1],
                gx.iat[i, 2],
                i,
                gx,
                tardict,
                refdict,
            )
        )
    gx2 = pd.DataFrame(rows, columns=range(8))
    print(gx2)
    gx2.to_excel(
        "./working/{}/matrix_foldseek_filtered1_{}.xlsx".format(name, name)
    )
