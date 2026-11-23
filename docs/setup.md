# Setup & External Resources

[← Back to README](../README.md) · [简体中文首页](../README.zh-CN.md)

Geno2Flux splits lightweight local analysis from heavy genome-annotation and AlphaGEM reconstruction dependencies, so you only install what your workflow actually needs.

---

## 1. Platform Support & Per-Stage Requirements

| Stage | Included after `bash run.sh setup` | Additional external requirements | Platform |
| --- | --- | --- | --- |
| Core environment & tests (`setup`, `make test`, `check-env`) | Locked Python 3.11–3.12 env via [`uv`](https://docs.astral.sh/uv/) (`uv.lock`), bundled [HiGHS](https://highs.dev/) (`highspy`) for LP **and** MILP | [`uv`](https://docs.astral.sh/uv/) (auto-manages Python 3.11–3.12) | Linux, macOS |
| GFF3 standardization & proteome prep (`predict-proteome --tool gff`, `prepare-inputs`) | Built-in pure-Python CDS splicer fallback when `gffread` is absent | Optional `gffread`; external predictor (`ANNEVO`, `Helixer`, `Tiberius`, or `funannotate`) only for de novo gene prediction | Linux, macOS (predictors vary) |
| Target screening on an existing GEM (`target`, `fseof`, `paralogs`, `knockout`, `filter`) | Full HiGHS LP backend | DIAMOND 2.1.10 — or a precomputed `diamond_self_hits` table | Linux, macOS |
| Minimum-set gap-filling (`gapfill`) | HiGHS MILP backend | NCBI BLAST+ 2.16.0 (`makeblastdb`, `tblastn`; via `fetch-resources --profile gapfill`) | Linux, macOS |
| Subcellular localization refinement (`localize`) | Auto-detection + automatic fallback to default compartments | Individually licensed DeepLoc 2.0/2.1 from DTU — or a precomputed `localization_tsv` | Linux, macOS |
| Fresh AlphaGEM reconstruction (Steps 1–3, `AlphaGEM.py --profile paper`) | GPL-3.0 AlphaGEM snapshot (`third_party/alphagem/`), Yeast9 template, reference reaction pool | Separate AlphaGEM & DeepEC conda environments, ~13 GB Figshare data, ESM-1b weights, and the verified eggNOG `fungi.dmnd` contract | **Linux x86_64 only** |

Check per-stage readiness at any time:

```sh
bash run.sh check-env
```

---

## 2. One-Click External Resource Download (`fetch-resources`)

`bash run.sh fetch-resources` downloads external binaries and databases directly from their official upstream hosts. Pass `--dry-run` first to preview every action.

Modular `--profile` values — `screening`, `gapfill`, `alphagem`, `annotation`, `deeploc`, `all` — compose cleanly (`gapfill` implies `screening`, and `alphagem` implies `gapfill`):

```sh
# Screening only: DIAMOND into tools/bin/
bash run.sh fetch-resources --profile screening

# Gap-filling: DIAMOND + NCBI BLAST+ 2.16.0 (tblastn, makeblastdb, blastp) into tools/bin/
bash run.sh fetch-resources --profile gapfill

# Fresh AlphaGEM reconstruction data: stages third_party/alphagem/ + downloads & verifies
# Figshare data_source.tar.xz (~13 GB) + ESM-1b weights + eggNOG (~41 GB eggnog.db + fungi.dmnd).
# Omit --with-eggnog only if eggNOG is already staged. Add --with-prott5 only for structure_alignment mode.
bash run.sh fetch-resources --profile alphagem --runtime-root /work/alphagem-runtime --with-eggnog

# De novo genome annotation: ANNEVO source + ANNEVO_Fungi.pt (add --with-helixer if needed)
bash run.sh fetch-resources --profile annotation --tool-root /work/annotation-tools
```

### What each profile installs

- **`screening`**: [DIAMOND](https://github.com/bbuchfink/diamond) v2.1.10 into `tools/bin/` on `linux-x86_64`. On macOS and `linux-aarch64`, no official DIAMOND release binary exists, so `fetch-resources` prints manual installation instructions (`conda install -c conda-forge -c bioconda diamond=2.1.10` or `brew install diamond`).
- **`gapfill`**: DIAMOND plus [NCBI BLAST+ 2.16.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/) (`tblastn`, `makeblastdb`, `blastp`) into `tools/bin/` from official archives for `linux-x86_64`, `linux-aarch64`, `darwin-aarch64`, and `darwin-x86_64`. (Note: `gapfill` itself uses BLAST+ only; DIAMOND is included because `--profile gapfill` implies `screening`.)
- **`alphagem`** (`--runtime-root`, optional `--with-eggnog` / `--with-prott5`):
  - Stages the verified `third_party/alphagem/` snapshot via `scripts/prepare_alphagem_runtime.py` and creates a clean local source-only git commit.
  - Downloads and verifies the ~13 GB [AlphaGEM Figshare](https://doi.org/10.6084/m9.figshare.30630524) `data_source.tar.xz`, unpacking `data_available/` (`yeast.fasta`, `yeast.xlsx`), `rhea/rhea_protembedding1~6.pkl`, `rhea/rhea2uniprot_sprot.tsv`, `rhea/rhea.rdf`, `plmsearchtools/`, and `tools/`.
  - Downloads [ESM-1b](https://github.com/facebookresearch/esm) weights (`esm1b_t33_650M_UR50S.pt`, ~7.8 GB, plus `esm1b_t33_650M_UR50S-contact-regression.pt`).
  - `--with-eggnog` downloads and decompresses `eggnog.db` (~41 GB, SHA-256 verified) and builds `fungi.dmnd` (exactly `879,076,733` bytes, SHA-256 verified) under `tools/eggnog_mapper/data/` ([eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)). Intermediate download files (`eggnog.db.gz`, `e5.proteomes.faa`, `e5.taxid_info.tsv`, `e5.proteomes.fungi.faa`) stay isolated in `<runtime>/var/downloads/eggnog/` so `tools/eggnog_mapper/data/` stays compliant with `--profile paper`.
  - `--with-prott5` downloads `Rostlab/prot_t5_xl_uniref50` via `huggingface-cli`/`hf` (only needed when `alphagem_mode` is `"structure_alignment"`).
- **`annotation`** (`--tool-root`, optional `--with-helixer`): fetches [ANNEVO](https://github.com/xjtu-omics/ANNEVO) source + `ANNEVO_Fungi.pt` (requires a separate PyTorch environment).

### Download verification coverage

Verification matches what each upstream publisher provides:

- **Pinned SHA-256 / MD5 digests**: Figshare `data_source.tar.xz` (SHA-256 + MD5 + size), NCBI BLAST+ 2.16.0 archives (MD5), ESM-1b contact-regression file (MD5 + size), and eggNOG assets — decompressed `eggnog.db` (SHA-256 + size), `e5.proteomes.faa` and `e5.taxid_info.tsv` (SHA-256), filtered Fungi FASTA (SHA-256), and built `fungi.dmnd` (SHA-256 + size).
- **Byte-size check only** (no digest published upstream): ESM-1b main weights (`esm1b_t33_650M_UR50S.pt`, ~7.8 GB) and `ANNEVO_Fungi.pt`.
- **Download / extraction check only** (no upstream checksum): DIAMOND v2.1.10 tarball, ANNEVO GitHub `main` source archive, `eggnog.db.gz` container (the decompressed `eggnog.db` is SHA-256 verified), optional Helixer weights, and HuggingFace `prot_t5_xl_uniref50`.

---

## 3. Optional DeepLoc 2.0 / 2.1 Setup & Fallback

DeepLoc 2.0 / 2.1 requires an individual academic license from DTU Health Tech and is never redistributed in this repository. Request access directly from DTU:

- DeepLoc 2.1: [web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) · [download request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All)
- DeepLoc 2.0: [web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) · [download request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.0&packageversion=2.0&platform=All)

Geno2Flux supports both paths automatically:

- **With DeepLoc (auto-detected)**: After downloading the DTU tarball, install it into [`tools/deeploc/`](../tools/deeploc/README.md) with:
  ```sh
  bash run.sh fetch-resources --profile deeploc --deeploc-tarball /path/to/DeepLoc-2.1.tar.gz
  ```
  You can also place or symlink `deeploc2` under `tools/deeploc/` (`tools/deeploc/bin/deeploc2` or `tools/deeploc/deeploc2`), set `G2F_DEEPLOC_BIN=/path/to/deeploc2`, put `deeploc2` on `PATH`, or provide a precomputed `localization_tsv`.
- **Without DeepLoc (graceful fallback)**: If neither `localization_tsv` nor DeepLoc is available, Geno2Flux logs a notice, retains the draft/template model's default compartments (`"deeploc_status": "fallback_default_compartments"` in `localization_refinement_audit.json`), and continues through gap-filling and target screening. To enforce strict failure instead, set `"require_deeploc": true` in your config (or pass `--require-deeploc` to `bash run.sh stage localize`).

---

## 4. Separate AlphaGEM & DeepEC Environments (Fresh Steps 1–3 Only)

`bash run.sh setup` installs only the lightweight Geno2Flux core and HiGHS solver. Running fresh AlphaGEM reconstruction (`AlphaGEM.py --profile paper`) strictly requires **Linux x86_64** and two dedicated Python environments:

1. **AlphaGEM environment** (`"alphagem_python"` in config or `G2F_ALPHAGEM_PYTHON`):
   Create from [`third_party/alphagem/AlphaGEM.yml`](../third_party/alphagem/AlphaGEM.yml) (Linux x86_64; pins Python 3.10.15, COBRApy 0.29.1, BLAST+ 2.16.0, DIAMOND 2.1.10, Foldseek `9.427df8a`), then install PyTorch, `fair-esm==2.0.0`, and [CLEAN](https://github.com/tttianhao/CLEAN) from an isolated copy — see [`third_party/alphagem/setup.sh`](../third_party/alphagem/setup.sh). Do **not** install in-place inside `tools/CLEAN/app` or leave `__pycache__`/`.pyc` files under `tools/`.
2. **DeepEC environment** (`"deepec_python"`):
   Create from [`third_party/alphagem/deepecenvironment.yml`](../third_party/alphagem/deepecenvironment.yml) (Linux x86_64; pins Python 3.7.7, `torch==1.7.0`, `transformers==3.5.1`) for [DeepECtransformer](https://github.com/kaistsystemsbiology/DeepProZyme).
3. **`fungi.dmnd` and `eggnog_static_data` contract gate**:
   The unified pipeline invokes `AlphaGEM.py --profile paper`, which requires `fungi.dmnd` at exactly `879,076,733` bytes with the SHA-256 produced by Linux x86_64 DIAMOND 2.1.10. Building it on another OS or DIAMOND version causes `fetch-resources --with-eggnog` to exit with `contract_mismatch` (exit code 2) and `check-env` to report `BLOCKED`. `check-env --profile alphagem` also verifies `eggnog.db` (size + SHA-256) and the 9-file `eggnog_static_data` bundle under `tools/eggnog_mapper/data/` — it must contain exactly those 9 archive files (including the archive's 472 KB `e5.proteomes.faa`) plus `eggnog.db` and `fungi.dmnd`, with no download intermediates (`e5.taxid_info.tsv`, `e5.proteomes.fungi.faa`) or extra/overwritten files.
4. **Optional `structure_alignment` mode**:
   The default mode is `plmsearch` ([PLMSearch](https://github.com/maovshao/PLMSearch)). Using `"alphagem_mode": "structure_alignment"` additionally requires `prottrans_model_path`, `structure_maplist`, and `structure_dir` (containing `*.pdb`) in config, plus Foldseek `9.427df8a` and the Linux x86_64 `USalign` binary from `data_source.tar.xz`. `structure_maplist` is an `.xlsx` workbook containing `Entry`, `Entry Name`, and a supported structure column (such as `Structure Data`, `AlphaFoldDB`, or `PDB`) mapping each UniProt entry to its structure file under `structure_dir`.
5. **Verify the environment handoff**:
   ```sh
   bash run.sh check-env --profile alphagem \
     --runtime-root /work/alphagem-runtime \
     --alphagem-python /path/to/AlphaGEM/bin/python \
     --deepec-python /path/to/deepectransformer/bin/python
   ```
6. **Runtime ownership rule**:
   When `"alphagem_runtime_root": null` (pipeline-managed runtime), set `reference_fasta` and `reference_xlsx` to your `<refname>.fasta` and `<refname>.xlsx` files; the pipeline copies them into its managed `data_available/`. When pointing to a separately prepared runtime, place both files under `<runtime>/data_available/` beforehand — external runtimes are treated as read-only.

The repository already ships the AlphaGEM source snapshot, the CC BY 4.0 Yeast9 template (`third_party/alphagem/models/yeast-GEM.xml`), and the reference reaction pool (`third_party/alphagem/models/reaction_pool_final.json`). External databases and trained weights are fetched on demand via `fetch-resources`.
