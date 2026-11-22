# AlphaGEM

AlphaGEM reconstructs genome-scale metabolic models from protein sequences using sequence homology, protein structure alignment and functional annotation.

## Usage

Run from the repository root:

```sh
uv run python third_party/alphagem/AlphaGEM.py --help
uv run python scripts/prepare_alphagem_runtime.py --runtime-root /data/alphagem-runtime --prepare
```

The preparation command creates a separate source runtime. Install the required environments, executables, weights and databases using the [installation instructions](UPSTREAM_README.md). These external resources are supplied separately.

Preparing source is not the same as preparing reconstruction data. Fresh reconstruction also needs `data_available/<refname>.fasta` and `data_available/<refname>.xlsx`. When the unified pipeline points `alphagem_runtime_root` at this external runtime, these files must already be present because that runtime is read-only. With `alphagem_runtime_root: null`, the pipeline creates its own runtime copy and copies the configured `reference_fasta` and `reference_xlsx` into it.


| Included resource | Purpose |
| --- | --- |
| `AlphaGEM.py` | Reconstruction command-line interface. |
| `src/` | Reconstruction, reaction-pool, gap-filling and solver modules. |
| `models/yeast-GEM.xml` | Yeast9 reference model. |
| `models/reaction_pool_final.json` | Reference reaction pool. |
| `SOURCE_MANIFEST.json` | Source identity, modification notices and file checksums. |

The bundled source is derived from the [authors' AlphaGEM repository](https://github.com/hongzhonglu/AlphaGEM), frozen from an intermediate reproduction snapshot (identified by `source_commit` in `SOURCE_MANIFEST.json`, not an upstream commit SHA), and includes the Geno2Flux modifications and per-file SHA-256 checksums recorded in `SOURCE_MANIFEST.json`.

AlphaGEM is licensed under [GPL-3.0](LICENSE). The [Yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) reference model retains [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); its SBML identifier specifies version 9.0.1, and the manifest identifies the exact bundled bytes.
