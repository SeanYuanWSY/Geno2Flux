# Third-party software and data

This repository provides orchestration, input-contract validation, compartment
refinement, configurable product-module mounting, and dual-metric metabolic
target prioritization code.

An AlphaGEM source snapshot is included at `third_party/alphagem/`, with its
original GPL-3.0 license, source commit, and SHA-256 file manifest
(`SOURCE_MANIFEST.json`). Draft reconstruction (unified-pipeline steps 1–3)
invokes `AlphaGEM.py` in an isolated runtime directory. Gap-filling (step 5)
imports `src/gapfilling.py` and the repository-authored HiGHS solver shim
`src/optlang_highs_interface.py` in-process. The frozen Yeast9 model in
`third_party/alphagem/models/yeast-GEM.xml` retains its
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) attribution to the
[yeast-GEM project](https://github.com/SysBioChalmers/yeast-GEM).

External tools, pretrained weights, and reference databases remain separate
dependencies obtained from their official sources and are not bundled in this
repository.

| Resource | Official source |
|---|---|
| AlphaGEM (bundled source snapshot; GPL-3.0) | https://github.com/hongzhonglu/AlphaGEM |
| AlphaGEM external resources | https://doi.org/10.6084/m9.figshare.30630524 |
| yeast-GEM / Yeast9 v9.0.1 (bundled template model; CC BY 4.0) | https://github.com/SysBioChalmers/yeast-GEM |
| RAVEN Toolbox (algorithm reference for `published_fseof.py`; GPL-3.0-only) | https://github.com/SysBioChalmers/RAVEN |
| COBRApy / optlang | [COBRApy](https://github.com/opencobra/cobrapy), [optlang](https://optlang.readthedocs.io/en/latest/) |
| HiGHS / highspy (default LP/MILP backend; MIT) | https://highs.dev/ |
| GLPK (optional LP/MILP backend) | https://www.gnu.org/software/glpk/ |
| Biopython / SciPy / NumPy / pandas | [Biopython](https://biopython.org/), [SciPy](https://scipy.org/) |
| DeepLoc 2.0 / 2.1 (external localization predictor; academic/commercial terms) | [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) / [DeepLoc 2.0](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) · [academic download request](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All) |
| DIAMOND (external sequence aligner for paralog clustering; GPL-3.0) | https://github.com/bbuchfink/diamond |
| NCBI BLAST+ (external `tblastn` / `makeblastdb` for genome gap-fill evidence) | [docs](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) · [2.16.0 downloads](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/) |
| gffread (external GFF3/CDS extraction tool) | https://github.com/gpertea/gffread |
| ANNEVO / Helixer / funannotate (optional external genome annotation tools) | [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Helixer](https://github.com/usadellab/Helixer), [funannotate](https://github.com/nextgenusfs/funannotate) |
| RepeatModeler / RepeatMasker (optional external repeat masking tools) | https://www.repeatmasker.org/ |
| ESM-1b / PLMSearch / CLEAN / DeepECtransformer / eggNOG-mapper (external AlphaGEM stages) | [ESM](https://github.com/facebookresearch/esm), [PLMSearch](https://github.com/maovshao/PLMSearch), [CLEAN](https://github.com/tttianhao/CLEAN), [DeepECtransformer](https://github.com/kaistsystemsbiology/DeepProZyme), [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) |
| Rhea / ChEBI / MetaNetX (reference reaction pool cross-references) | [Rhea](https://www.rhea-db.org/), [MetaNetX](https://www.metanetx.org/) |

The Python dependency lockfile `uv.lock` records exact package versions and
distribution hashes. Sample genomes, private strain assemblies, pretrained
weights, and external database archives are excluded from this repository.

The bundled `yeast-GEM.xml` is the unchanged Yeast9 v9.0.1 reference model
provided with AlphaGEM. Credit the Yeast-GEM contributors and
[Yeast9 (2024)](https://doi.org/10.1038/s44320-024-00060-7), retain the
[CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/), and identify
any subsequent changes.

The bundled `reaction_pool_final.json` is the Git-tracked reaction pool from
the GPL-3.0 AlphaGEM source repository. Retain the AlphaGEM license and credit
Rhea, ChEBI, and MetaNetX/MNXref for the reference data.

## External resource distribution

External tool packages, pretrained weights, and annotation databases must be
obtained separately from their official sources:

- **DeepLoc 2.0 / 2.1**: DTU's academic license does not permit redistributing
  the software or weights inside third-party packages. Obtain a license directly
  from DTU Health Tech, or supply a precomputed localization TSV/CSV table.
- **ANNEVO**: carries non-commercial terms and is obtained separately from its
  official repository.
- **AlphaGEM Figshare archive**: linked at its official DOI rather than mirrored
  in this repository.

`bash run.sh fetch-resources` automates fetching and checksum-verifying the
publicly redistributable upstream assets directly from their official hosts
(Figshare, NCBI FTP, GitHub, Meta/FAIR ESM, eggNOG, Zenodo). It never
redistributes or bypasses DTU's DeepLoc license — DeepLoc still requires an
individual academic request and a user-supplied local tarball.

## Included code and models

Copyright (c) 2026 Geno2Flux contributors. The Geno2Flux-owned source is
licensed under `GPL-3.0-only`; see `LICENSE`. This grant does not relicense
third-party software, models, or data. Distribution of the combined program,
including its GPL-covered AlphaGEM modules (`third_party/alphagem/`) and RAVEN
FSEOF port (`src/geno2flux/published_fseof.py`), complies with `GPL-3.0`.

The algorithmic port in `src/geno2flux/published_fseof.py` follows RAVEN
`FSEOF.m` and `solveLP.m` at commit `0e26f3ba558f86115c830d1c4994b10e4d788f37`
and carries `GPL-3.0-only`. The upstream RAVEN license is retained verbatim in
`third_party/raven-license/LICENSE.md`.
