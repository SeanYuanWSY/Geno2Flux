# Configuration & CLI Reference

[← Back to README](../README.md) · [简体中文首页](../README.zh-CN.md)

Start from the annotated template [`configs/unified_pipeline.example.json`](../configs/unified_pipeline.example.json) and copy it to a local file such as `/work/run.local.json`.

---

## 1. Configuring the Target Product (Step 6)

Geno2Flux makes no assumptions about your strain or target product. Configure either an existing model reaction or a custom heterologous pathway:

### Option 1: Existing reaction in the model

```json
{
  "target_reaction": "r_2056",
  "target_module": null
}
```

`r_2056` is the succinate-exchange reaction ID in the bundled Yeast9 template (its BiGG annotation is `EX_succ_e`). Reaction IDs are model-specific — annotation aliases cannot be passed as `target_reaction`. To list candidate exchange reactions in your own SBML model:

```sh
uv run --frozen python - <<'PY'
import cobra
model = cobra.io.read_sbml_model("path/to/model.xml")
for reaction in model.exchanges:
    print(reaction.id, reaction.name, reaction.bounds, reaction.annotation)
PY
```

### Option 2: Custom stoichiometric pathway (`target_module`)

Pass an inline JSON object or a path to a `.json` file. Internal reactions are checked for elemental and charge balance before screening:

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

---

## 2. Existing-GEM Fast Path & SHA-256 Reuse Gate

### Minimal config for an existing GEM

If you already have a curated SBML/YAML/JSON model and matching proteome FASTA, skip reconstruction, localization, and gap-filling with:

```json
{
  "localization_order": "skip",
  "gapfilled_model": "/work/model.xml",
  "proteome_fasta": "/work/proteome.fasta",
  "id_map": null,
  "biomass": "<growth_reaction_id>",
  "target_reaction": "<target_reaction_id>",
  "compute_growth_essentiality": true,
  "allow_unverified_reuse": true
}
```

Only DIAMOND is needed for paralog clustering in this mode — or supply a precomputed `diamond_self_hits` table to skip DIAMOND altogether.

### SHA-256 verification gate

When reusing pre-existing files (`proteome_fasta`, `withgaps_model`, `gapfilled_model`, `localization_tsv`, `template_model`, or `diamond_self_hits`), Geno2Flux expects the corresponding `*_sha256` checksum field:

- `proteome_sha256`
- `withgaps_sha256`
- `gapfilled_model_sha256`
- `localization_sha256`
- `template_model_sha256`
- `diamond_self_hits_sha256`

For exploratory runs where you want to skip hash pinning, set `"allow_unverified_reuse": true` explicitly.

---

## 3. Full Genome-to-Target Workflow

### Step 1: Predict proteome & build coordinate contracts

Optionally fetch NCBI RefSeq or UniProt reference proteomes for homology-guided predictors (`funannotate` / `miniprot`):

```sh
bash run.sh fetch-reference-proteomes \
  --uniprot-proteome YEAST=UP000002311 \
  --output /work/ref_proteomes
```

Run `predict-proteome` to invoke a predictor (`annevo`, `helixer`, `tiberius`, `funannotate`, `custom`) or standardize an existing coding `gff`, extract the longest-isoform proteome, run QC, and emit `/work/pred/inputs/`:

```sh
# From an existing coding GFF3 (works even without gffread):
bash run.sh predict-proteome \
  --tool gff \
  --genome /work/raw/genome.fasta \
  --gff /work/raw/prediction.gff3 \
  --case-id sample_strain \
  --prefix PROT \
  --output /work/pred

# Or invoke a de novo predictor directly (example: ANNEVO):
bash run.sh predict-proteome \
  --tool annevo \
  --genome /work/raw/genome.fasta \
  --annevo-dir /path/to/ANNEVO \
  --model-weights /path/to/ANNEVO_Fungi.pt \
  --tool-python /path/to/annevo-env/bin/python \
  --case-id sample_strain \
  --output /work/pred
```

Supported upstream predictors: [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Helixer](https://github.com/usadellab/Helixer), [Tiberius](https://github.com/Gaius-Augustus/Tiberius), [funannotate](https://github.com/nextgenusfs/funannotate).

### Step 2: Run the unified pipeline

Point the proteome and coordinate fields in `/work/run.local.json` at `/work/pred/inputs/`. For a fresh AlphaGEM run, merge the `config_snippet` printed by `bash run.sh fetch-resources --profile alphagem`:

- It populates `alphagem_runtime_root`, `reference_fasta`, `reference_xlsx`, and `esm_model_path`.
- Replace the placeholder values for `alphagem_python` and `deepec_python` with your actual environment Python paths.
- Keep `"prottrans_model_path": null` for the default `plmsearch` mode.

Then run:

```sh
bash run.sh metabolic --config /work/run.local.json --output /work/run-001
```

### Output receipts & tables

Generated under `/work/run-001/<order>/` (e.g., `pre-gapfill/`, or `screening/` when `localization_order` is `"skip"`):

- `localization_refinement_audit.json` — compartment reactions pulled from DeepLoc, or the default-compartment fallback receipt.
- `target_model_<order_slug>.xml` & `target_module_receipt.json` — final SBML model with the target reaction or custom module mounted (plus elemental/charge balance checks).
- `fseof/all_genes_fseof_essential_filtered.tsv` — non-essential overexpression candidates ranked by 10-step FSEOF flux slopes.
- `paralog/paralog_family_rankings_essential_filtered.tsv` — homologous gene families clustered via DIAMOND self-alignment and ranked by aggregate FSEOF score.
- `knockout/single_gene_ko_essential_filtered.tsv` & `knockout/family_ko_essential_filtered.tsv` — single-gene and homologous-family joint knockout rankings with growth-essential genes filtered out.

The four `*_essential_filtered.tsv` tables are generated when `"compute_growth_essentiality": true` or `essential_genes_tsv` is set; if both are omitted, only the unfiltered ranking tables are written.

---

## 4. Running Individual Stages

To inspect or rerun stages one by one, prepare a standalone AlphaGEM runtime first (`uv run --frozen python scripts/prepare_alphagem_runtime.py --prepare --runtime-root /work/alphagem-runtime`), then use `bash run.sh stage`:

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

Run `bash run.sh help` to list every subcommand.

---

## 5. Repository Layout

- `src/geno2flux/` — unified pipeline, subcellular compartment refinement (DeepLoc auto-detection & fallback), configurable target module, RAVEN-compatible FSEOF port, DIAMOND paralog grouping, knockout screening, growth-essentiality filtering, and native HiGHS backend.
- `scripts/` — CLI entry points: genome-to-proteome prediction (`predict_proteome.py`), RefSeq/UniProt reference proteome fetching (`fetch_reference_proteomes.py`), external resource download & preflight (`fetch_external_resources.py`), GFF3 standardization (`standardize_annotation.py`), coordinate contract building (`build_coordinate_contract.py`), AlphaGEM runtime staging (`prepare_alphagem_runtime.py`), and step-by-step execution (`run_stage.py`).
- `tools/deeploc/` — drop-in slot for individually licensed DeepLoc 2.0/2.1 from DTU (gitignored except `README.md`).
- `tools/bin/` — drop-in directory for downloaded external binaries (DIAMOND, BLAST+), automatically prepended to `PATH` by `run.sh` (gitignored except `README.md`).
- `configs/unified_pipeline.example.json` — annotated configuration template for full genome-to-target runs.
- `third_party/alphagem/` — bundled GPL-3.0 AlphaGEM source snapshot with HiGHS solver interface, coordinate contracts, CC BY 4.0 `yeast-GEM.xml` template, and `models/reaction_pool_final.json`.
