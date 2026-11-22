<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/geno2flux-logo-dark.svg">
    <img alt="Geno2Flux" src="assets/geno2flux-logo-light.svg" width="420">
  </picture>
</p>

<p align="center"><strong>English</strong> · <a href="README.zh-CN.md">简体中文</a></p>

Geno2Flux reconstructs genome-scale metabolic models (GEMs) for eukaryotic/fungal strains from your own genome/proteome data and prioritizes metabolic engineering targets using a dual-metric strategy:

1. **Genome-to-GEM reconstruction & refinement**: genome-to-proteome prediction & coding GFF3 standardization (`ANNEVO` / `Helixer` / `Tiberius` / `funannotate` + RefSeq/UniProt protein evidence, or direct GFF3 import) → AlphaGEM + Yeast9 draft reconstruction → subcellular localization compartment refinement (auto-detected DeepLoc 2.0/2.1, precomputed table, or automatic graceful fallback to default compartments) → genome-backed minimum-set gap-filling (NCBI BLAST+ + HiGHS).
2. **Dual-metric target prioritization**: configurable product pathway or endogenous reaction target → 10-step RAVEN-compatible FSEOF overexpression scanning → DIAMOND paralog-family clustering → single-reaction, single-gene, and homologous-family joint knockout screening → growth-essentiality filtering.

---

## 1. Included vs. External Requirements & Platform Boundaries

| Stage | Included in this repo | Requires external setup |
| --- | --- | --- |
| Core environment & tests (`setup`, `make test`, `check-env`) | Python 3.11–3.12 via [`uv`](https://docs.astral.sh/uv/) (`uv.lock`), bundled [HiGHS](https://highs.dev/) (`highspy`) for LP **and** MILP | [`uv`](https://docs.astral.sh/uv/) must be installed first — it auto-manages Python 3.11–3.12 for `bash run.sh setup`; works on Linux and macOS |
| GFF3 standardization & proteome prep (`predict-proteome --tool gff`, `prepare-inputs`) | Pure-Python CDS splicer fallback when `gffread` is absent | Optional `gffread`; an external predictor (ANNEVO / Helixer / Tiberius / funannotate) only for *fresh* gene prediction |
| Target prioritization on an existing GEM (`target`, `fseof`, `paralogs`, `knockout`, `filter`) | Full HiGHS LP backend | DIAMOND 2.1.10 — or supply precomputed `diamond_self_hits` to bypass DIAMOND |
| Minimum-set gap-filling (`gapfill`) | HiGHS MILP backend | NCBI BLAST+ 2.16.0 (`makeblastdb`, `tblastn`; `fetch-resources --profile gapfill`) — `--profile gapfill` also installs DIAMOND because it implies `screening`, but `gapfill` itself does not use DIAMOND |
| DeepLoc localization refinement (`localize`) | Auto-detection and automatic fallback to default compartments | Individually licensed DeepLoc 2.0/2.1 tarball from DTU — or a precomputed `localization_tsv` |
| Fresh AlphaGEM reconstruction (Steps 1–3, `AlphaGEM.py --profile paper`) | GPL-3.0 AlphaGEM snapshot in `third_party/alphagem/`, Yeast9 template, reaction pool | **Linux x86_64 only**, plus separate AlphaGEM and DeepEC conda environments, ~13 GB Figshare data, ESM-1b weights, and the verified eggNOG `fungi.dmnd` contract |

**Boundary in one sentence:** [`uv`](https://docs.astral.sh/uv/) (which auto-manages Python 3.11–3.12) is the only prerequisite for `bash run.sh setup`; on Linux and macOS, `git clone` + `bash run.sh setup` is then sufficient for tests, GFF3 standardization, and the HiGHS LP stages (`target`, `fseof`, `knockout`, `filter`; `localize` falls back to default compartments); full target prioritization on an existing GEM (`localization_order: "skip"` + `gapfilled_model`) additionally requires DIAMOND 2.1.10 or a precomputed `diamond_self_hits` table for paralog clustering, `gapfill` additionally requires NCBI BLAST+ 2.16.0 (`makeblastdb`, `tblastn`; `fetch-resources --profile gapfill` installs BLAST+ plus DIAMOND because the profile implies `screening`), and **de novo genome-to-GEM reconstruction (Steps 1–3) does not run from a bare clone** — it requires the external resources and Linux x86_64 environments described in Section 5.

---

## 2. Installation & Verification

Sync the locked Python environment with [`uv`](https://docs.astral.sh/uv/) (Python 3.11–3.12) and verify your local installation with the automated test suite:

```sh
bash run.sh setup
make test
```

The synced environment already bundles the open-source [HiGHS](https://highs.dev/) solver (`highspy`), covering both **LP** (FSEOF scanning, knockout screening, growth-essentiality filtering) and **MILP** (minimum-set gap-filling) — no commercial Gurobi license is required.

Inspect per-stage readiness (Python interpreters, external binaries, downloaded resources) at any time:

```sh
bash run.sh check-env
```

See all CLI entry points at any time with:

```sh
bash run.sh help
```

---

## 3. Target Configuration (Existing GEM or Custom Pathway)

Geno2Flux does not assume a product, pathway, or strain. You must configure Step 6 in your JSON config (or `scripts/run_stage.py target`) in one of two ways:

1. **Existing reaction in the model** (no new reactions added):
   ```json
   {
     "target_reaction": "r_2056",
     "target_module": null
   }
   ```

   `r_2056` is the actual succinate-exchange reaction ID in the bundled Yeast9 model. Its BiGG annotation is `EX_succ_e`. IDs are model-specific, and annotation aliases cannot be used as reaction IDs. To inspect candidate exchange reactions in another SBML model:
   ```sh
   uv run --frozen python - <<'PY'
   import cobra
   model = cobra.io.read_sbml_model("path/to/model.xml")
   for reaction in model.exchanges:
       print(reaction.id, reaction.name, reaction.bounds, reaction.annotation)
   PY
   ```
2. **Custom stoichiometric module** (inline object or path to a `.json` file; internal reactions are checked for elemental and charge balance before screening):
   ```json
   {
     "target_module": {
       "mode": "custom",
       "target_reaction": "DM_target_c",
       "expected_metabolites": {
         "chebi30031_c": { "formula": "C4H4O4", "charge": -2, "compartment": "c" }
       },
       "new_metabolites": [
         { "id": "target_product_c", "formula": "C4H4O4", "charge": -2, "compartment": "c" }
       ],
       "reactions": [
         { "id": "TARGET_SYN", "stoichiometry": { "chebi30031_c": -1.0, "target_product_c": 1.0 } },
         { "id": "DM_target_c", "stoichiometry": { "target_product_c": -1.0 } }
       ]
     }
   }
   ```

**Existing-GEM fast path.** If you already have a curated SBML/YAML/JSON GEM and a proteome FASTA, run target prioritization directly (skipping AlphaGEM reconstruction, localization, and gap-filling) by setting:

```json
{
  "localization_order": "skip",
  "gapfilled_model": "path/to/model.xml"
}
```

Only DIAMOND is needed for paralog clustering in this mode — or supply a precomputed `diamond_self_hits` table to bypass DIAMOND entirely.

**SHA-256 reuse gate.** Reusing `proteome_fasta`, `withgaps_model`, `gapfilled_model`, `localization_tsv`, or `template_model` requires the matching `*_sha256` field (`proteome_sha256`, `withgaps_sha256`, `gapfilled_model_sha256`, `localization_sha256`, `template_model_sha256`; likewise `diamond_self_hits_sha256` for `diamond_self_hits`). For exploratory runs you may instead set `"allow_unverified_reuse": true` explicitly.

---

## 4. DeepLoc Plug-and-Play Slot & Automatic Fallback

DeepLoc 2.0 / 2.1 requires an individual academic license from DTU Health Tech and cannot be redistributed inside public repositories. Request access directly from DTU — [DeepLoc 2.1 web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) · [2.1 download request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All), or [DeepLoc 2.0 web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) · [2.0 download request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.0&packageversion=2.0&platform=All). Geno2Flux handles both cases seamlessly:

- **With DeepLoc (zero-config auto-detection)**: Once you obtain the DTU tarball, run `bash run.sh fetch-resources --profile deeploc --deeploc-tarball /path/to/DeepLoc-2.1.tar.gz` for a one-click install into [`tools/deeploc/`](tools/deeploc/README.md) — or manually place/symlink `deeploc2` under `tools/deeploc/` (e.g. `tools/deeploc/bin/deeploc2` or `tools/deeploc/deeploc2`), set `G2F_DEEPLOC_BIN=/path/to/deeploc2`, or place `deeploc2` on `PATH`. You can also supply a precomputed `localization_tsv`. Geno2Flux detects it automatically and refines missing compartments.
- **Without DeepLoc (automatic default fallback)**: If no `localization_tsv` is supplied and DeepLoc is not installed, Geno2Flux **does not crash** — it logs a notice, keeps the draft/template model's default compartments (`"deeploc_status": "fallback_default_compartments"` in `localization_refinement_audit.json`), and runs the rest of the pipeline cleanly. Set `"require_deeploc": true` in config (or `--require-deeploc` on `run_stage.py localize`) if you want strict failure instead.

---

## 5. External Resources & AlphaGEM Environment Setup

### One-click resource download & readiness preflight

`bash run.sh fetch-resources` fetches external tools and databases straight from their official upstream hosts (`--dry-run` previews every action first). Verification strength differs per asset:

- **Pinned SHA-256 / MD5 digests**: Figshare `data_source.tar.xz` (SHA-256 + MD5 + size), NCBI BLAST+ 2.16.0 archives (MD5), the ESM-1b contact-regression file (MD5 + size), and the eggNOG assets — `eggnog.db` (SHA-256 + size after decompression), `e5.proteomes.faa` and `e5.taxid_info.tsv` (SHA-256), the filtered Fungi FASTA (SHA-256 against the pinned contract), and the built `fungi.dmnd` (SHA-256 + size).
- **Byte-size check only** (upstream publishes no digest): the ESM-1b main weights (`esm1b_t33_650M_UR50S.pt`, ~7.8 GB) and `ANNEVO_Fungi.pt`.
- **No published upstream checksum** (verified by successful download/extraction): the DIAMOND v2.1.10 Linux tarball, the ANNEVO GitHub `main` source archive, the `eggnog.db.gz` container (its *decompressed* `eggnog.db` is SHA-256 verified), optional Helixer weights, and HuggingFace `prot_t5_xl_uniref50`.

Modular `--profile` values — `screening`, `gapfill`, `alphagem`, `annotation`, `deeploc`, `all` — compose so that `gapfill` implies `screening` and `alphagem` implies `gapfill`:

```sh
# Screening only: DIAMOND into tools/bin/
bash run.sh fetch-resources --profile screening

# Gap-filling: DIAMOND + NCBI BLAST+ 2.16.0 (tblastn, makeblastdb, blastp) into tools/bin/
bash run.sh fetch-resources --profile gapfill

# AlphaGEM data (fresh reconstruction): stages the verified third_party/alphagem/
# source snapshot + downloads & verifies the ~13 GB Figshare data_source.tar.xz
# + ESM-1b weights + eggNOG (~41 GB eggnog.db + fungi.dmnd, required before draft
# output; omit --with-eggnog only if already staged). Add --with-prott5 only for
# alphagem_mode=structure_alignment.
bash run.sh fetch-resources --profile alphagem --runtime-root /work/alphagem-runtime --with-eggnog

# De novo genome annotation: ANNEVO source + ANNEVO_Fungi.pt (+ --with-helixer)
bash run.sh fetch-resources --profile annotation --tool-root /work/annotation-tools
```

What each profile installs:

- **`screening`**: [DIAMOND](https://github.com/bbuchfink/diamond) v2.1.10 into `tools/bin/` on `linux-x86_64`. On macOS and `linux-aarch64`, no official DIAMOND release binary exists, so `fetch-resources` prints manual install instructions instead (`conda install -c conda-forge -c bioconda diamond=2.1.10` or `brew install diamond`).
- **`gapfill`**: DIAMOND plus [NCBI BLAST+ 2.16.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/) (`tblastn`, `makeblastdb`, `blastp`) into `tools/bin/`, from official archives for `linux-x86_64`, `linux-aarch64`, `darwin-aarch64`, and `darwin-x86_64`.
- **`alphagem`** (`--runtime-root`, optional `--with-eggnog` / `--with-prott5`):
  - Stages the verified `third_party/alphagem/` source snapshot via `prepare_alphagem_runtime.py` and initializes a clean local source-only git commit.
  - Downloads and verifies the ~13 GB [AlphaGEM Figshare](https://doi.org/10.6084/m9.figshare.30630524) `data_source.tar.xz`, unpacking `data_available/` (`yeast.fasta`, `yeast.xlsx`), `rhea/rhea_protembedding1~6.pkl`, `rhea/rhea2uniprot_sprot.tsv`, `rhea/rhea.rdf`, `plmsearchtools/`, and `tools/`.
  - Downloads [ESM-1b](https://github.com/facebookresearch/esm) weights (`esm1b_t33_650M_UR50S.pt`, ~7.8 GB, plus `esm1b_t33_650M_UR50S-contact-regression.pt`).
  - `--with-eggnog` downloads and decompresses `eggnog.db` (~41 GB, SHA-256 verified) and builds `fungi.dmnd` (exactly `879,076,733` bytes, SHA-256 verified) under `tools/eggnog_mapper/data/` ([eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)). Intermediate files (`eggnog.db.gz`, `e5.proteomes.faa`, `e5.taxid_info.tsv`, `e5.proteomes.fungi.faa`) stay isolated in `<runtime>/var/downloads/eggnog/` so `tools/eggnog_mapper/data/` remains compliant with `--profile paper`.
  - `--with-prott5` downloads `Rostlab/prot_t5_xl_uniref50` via `huggingface-cli`/`hf` — only needed for `alphagem_mode: "structure_alignment"`.
- **`annotation`** (`--tool-root`, optional `--with-helixer`): fetches [ANNEVO](https://github.com/xjtu-omics/ANNEVO) source + `ANNEVO_Fungi.pt` (requires a separate PyTorch environment).

### Separate AlphaGEM & DeepEC environments (fresh Steps 1–3 only)

`bash run.sh setup` (`uv`) installs only the lightweight Geno2Flux core + HiGHS solver. Running fresh AlphaGEM Steps 1–3 (`AlphaGEM.py --profile paper`) strictly requires **Linux x86_64** and two separate Python environments:

1. **AlphaGEM environment** (`"alphagem_python"` in config or `G2F_ALPHAGEM_PYTHON`): create from [`third_party/alphagem/AlphaGEM.yml`](third_party/alphagem/AlphaGEM.yml) (Linux x86_64; pins Python 3.10.15, COBRApy 0.29.1, BLAST+ 2.16.0, DIAMOND 2.1.10, Foldseek `9.427df8a`), **and** install PyTorch, `fair-esm==2.0.0`, and [CLEAN](https://github.com/tttianhao/CLEAN) from an isolated copy — see [`third_party/alphagem/setup.sh`](third_party/alphagem/setup.sh); do **not** install in-place inside `tools/CLEAN/app` or leave `__pycache__`/`.pyc` files under `tools/`.
2. **DeepEC environment** (`"deepec_python"`): create from [`third_party/alphagem/deepecenvironment.yml`](third_party/alphagem/deepecenvironment.yml) (Linux x86_64; pins Python 3.7.7, `torch==1.7.0`, `transformers==3.5.1`) for [DeepECtransformer](https://github.com/kaistsystemsbiology/DeepProZyme).
3. **`fungi.dmnd` contract gate**: the unified pipeline invokes `AlphaGEM.py --profile paper`, which requires `fungi.dmnd` at exactly `879,076,733` bytes with the SHA-256 matching Linux x86_64 DIAMOND 2.1.10. Building it on a non-Linux platform or with another DIAMOND version makes `fetch-resources --with-eggnog` exit with `contract_mismatch` (exit code 2), and `check-env` reports `BLOCKED`. `check-env --profile alphagem` also validates `eggnog.db` (size + SHA-256) and the 9-file `eggnog_static_data` bundle contract under `tools/eggnog_mapper/data/` — the bundle must hold exactly those 9 files plus `eggnog.db` and `fungi.dmnd`, with intermediates (`eggnog.db.gz`, `e5.proteomes.faa`, `e5.taxid_info.tsv`, `e5.proteomes.fungi.faa`) rejected.
4. **`alphagem_mode: "structure_alignment"`** additionally requires `prottrans_model_path`, `structure_maplist`, and `structure_dir` (containing `*.pdb`) in the config, plus Foldseek `9.427df8a` and the Linux x86_64 `USalign` binary shipped in `data_source.tar.xz`. `structure_maplist` is an Excel workbook (`.xlsx`) containing `Entry`, `Entry Name`, and a supported structure column (e.g. `Structure Data` / `AlphaFoldDB` / `PDB`) that maps each UniProt entry to its structure file under `structure_dir`. ([PLMSearch](https://github.com/maovshao/PLMSearch) is the default `plmsearch` mode.)
5. Verify the environment handoff with:
   ```sh
   bash run.sh check-env --profile alphagem --runtime-root /work/alphagem-runtime \
     --alphagem-python /path/to/AlphaGEM/bin/python \
     --deepec-python /path/to/deepectransformer/bin/python
   ```
6. **Runtime ownership rule**: for a pipeline-managed runtime (`"alphagem_runtime_root": null`), set `reference_fasta` and `reference_xlsx` to your `<refname>.fasta` / `<refname>.xlsx` files; the pipeline copies them into its own `data_available/`. For a separately prepared runtime, place both files under `<runtime>/data_available/` first — that runtime is treated as read-only.

The repository already includes the AlphaGEM source snapshot, the CC BY 4.0 Yeast9 template (`third_party/alphagem/models/yeast-GEM.xml`), and the reference reaction pool (`third_party/alphagem/models/reaction_pool_final.json`). It does not include reconstruction databases or trained weights — `fetch-resources` fetches them on demand.

---

## 6. Running the Workflow

### Step 1 — Predict proteome from genome (or standardize an existing coding GFF3)

Optionally download NCBI RefSeq or UniProt reference proteomes for homology-guided predictors (`funannotate` / `miniprot`):

```sh
bash run.sh fetch-reference-proteomes \
  --uniprot-proteome YEAST=UP000002311 \
  --output /work/ref_proteomes
```

Then run `predict-proteome` in one command to predict genes (via `annevo`, `helixer`, `tiberius`, `funannotate`, `custom`, or an existing `gff`), validate/extract the longest-isoform proteome, run QC, and emit the ready-to-use `/work/pred/inputs/` AlphaGEM target & coordinate contracts:

```sh
# From an existing coding GFF3 (works even without gffread installed):
bash run.sh predict-proteome \
  --tool gff \
  --genome /work/raw/genome.fasta \
  --gff /work/raw/prediction.gff3 \
  --case-id sample_strain \
  --prefix PROT \
  --output /work/pred

# Or directly invoke ANNEVO / Helixer / Tiberius / funannotate:
bash run.sh predict-proteome \
  --tool annevo \
  --genome /work/raw/genome.fasta \
  --annevo-dir /path/to/ANNEVO \
  --model-weights /path/to/ANNEVO_Fungi.pt \
  --tool-python /path/to/annevo-env/bin/python \
  --case-id sample_strain \
  --output /work/pred
```

De novo predictors from their official sources: [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Helixer](https://github.com/usadellab/Helixer), [Tiberius](https://github.com/Gaius-Augustus/Tiberius), [funannotate](https://github.com/nextgenusfs/funannotate).

### Step 2 — Run the unified pipeline (or individual stages)

Copy `configs/unified_pipeline.example.json` to `/work/run.local.json` and point its proteome/coordinate fields at `/work/pred/inputs/`. For a fresh AlphaGEM run, merge the `config_snippet` printed by `fetch-resources --profile alphagem` into `/work/run.local.json` — it fills in `alphagem_runtime_root`, `reference_fasta`, `reference_xlsx`, and `esm_model_path`, and emits placeholder strings for `alphagem_python` and `deepec_python` that you must replace with the actual Python paths of your AlphaGEM and DeepEC environments (keep `prottrans_model_path: null` for the default `plmsearch` mode). For exploratory runs, either record the SHA-256 hashes or set `"allow_unverified_reuse": true`. Then execute:

```sh
bash run.sh metabolic --config /work/run.local.json --output /work/run-001
```

Key outputs generated under `/work/run-001/<order>/`:

- `localization_refinement_audit.json` — pulled compartment reactions (or default-compartment fallback receipt).
- `target_model_<order>.xml` & `target_module_receipt.json` — mass- and charge-balanced target model.
- `fseof/all_genes_fseof_essential_filtered.tsv` — non-essential overexpression candidates ranked by 10-step FSEOF flux slopes.
- `paralog/paralog_family_rankings_essential_filtered.tsv` — homologous gene families clustered via DIAMOND self-alignment and ranked by aggregate FSEOF scores.
- `knockout/single_gene_ko_essential_filtered.tsv` & `knockout/family_ko_essential_filtered.tsv` — single-gene and homologous-family joint knockout rankings excluding growth-essential genes.

Alternatively, prepare a standalone AlphaGEM runtime (`uv run --frozen python scripts/prepare_alphagem_runtime.py --prepare --runtime-root /work/alphagem-runtime`) and run stages one at a time via `bash run.sh stage`:

```sh
export G2F_RUNTIME=/work/alphagem-runtime
export G2F_STAGE=/work/stages-001

# After AlphaGEM draft reconstruction (--name sample_strain --stop-after reaction_pool):
bash run.sh stage localize --name sample_strain --inputs /work/pred/inputs
bash run.sh stage gapfill  --name sample_strain --inputs /work/pred/inputs
bash run.sh stage target   --inputs /work/pred/inputs --target-reaction r_2056
bash run.sh stage fseof    --inputs /work/pred/inputs --growth-reaction r_4041
bash run.sh stage paralogs --inputs /work/pred/inputs
bash run.sh stage knockout --inputs /work/pred/inputs
bash run.sh stage filter   --inputs /work/pred/inputs --growth-reaction r_4041
```

---

## 7. Repository Layout

- `src/geno2flux/` — unified pipeline, subcellular compartment refinement (with DeepLoc auto-detection & default fallback), configurable target module, RAVEN-compatible FSEOF port, DIAMOND paralog grouping, knockout screening, growth-essentiality filtering, and native HiGHS backend.
- `scripts/` — CLI utilities for genome-to-proteome prediction (`predict_proteome.py`), RefSeq/UniProt reference proteome fetching (`fetch_reference_proteomes.py`), external resource download & readiness preflight (`fetch_external_resources.py`), GFF3 standardization (`standardize_annotation.py`), coordinate contract building (`build_coordinate_contract.py`), AlphaGEM runtime staging (`prepare_alphagem_runtime.py`), and step-by-step execution (`run_stage.py`).
- `tools/deeploc/` — plug-and-play drop-in slot for users who obtain academic DeepLoc 2.0/2.1 access from DTU (gitignored except `README.md`).
- `tools/bin/` — drop-in directory for downloaded external binaries (DIAMOND, BLAST+) that `run.sh` prepends to `PATH` automatically (gitignored except `README.md`).
- `configs/unified_pipeline.example.json` — annotated template configuration for full genome-to-target runs.
- `third_party/alphagem/` — bundled GPL-3.0 AlphaGEM source snapshot with HiGHS solver interface, coordinate contracts, and CC BY 4.0 `yeast-GEM.xml` template.

---

## 8. License & Third-Party Notices

- Geno2Flux is released under the GPL-3.0 license — see [LICENSE](LICENSE).
- `third_party/alphagem/` bundles the GPL-3.0 AlphaGEM source snapshot; the Yeast9 template is licensed CC BY 4.0. DeepLoc 2.0/2.1 remains under its individual DTU academic license and is never redistributed.
- Full attribution for every bundled and externally fetched component: [THIRD_PARTY.md](THIRD_PARTY.md). Python dependencies: [requirements.txt](requirements.txt).
