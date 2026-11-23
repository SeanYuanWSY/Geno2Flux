<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/geno2flux-logo-dark.svg">
    <img alt="Geno2Flux" src="assets/geno2flux-logo-light.svg" width="420">
  </picture>
</p>

<p align="center"><strong>English</strong> · <a href="README.zh-CN.md">简体中文</a></p>

Geno2Flux takes a fungal genome—or a genome-scale metabolic model (GEM) you already have—and produces ranked metabolic engineering targets for both overexpression and knockout. It bundles the open-source [HiGHS](https://highs.dev/) solver (`highspy`) for LP and MILP, so no commercial Gurobi license is needed.

---

## How it works

1. **Genome to proteome** — Predict genes from a genome assembly (`ANNEVO`, `Helixer`, `Tiberius`, or `funannotate`), or import an existing coding GFF3.
2. **Draft GEM & refinement** — Build a draft model with AlphaGEM + Yeast9, optionally refine subcellular compartments with DeepLoc (automatically falls back to default compartments if absent), and run genome-backed gap-filling with BLAST+ and HiGHS.
3. **Target pathway** — Pick an endogenous reaction already in your model, or mount a custom product pathway.
4. **Target prioritization** — Rank overexpression candidates with 10-step FSEOF, cluster paralog families with DIAMOND, screen single-gene and family knockouts, and filter out growth-essential genes.

---

## Quick Start

Requires [`uv`](https://docs.astral.sh/uv/) (which manages Python 3.11–3.12 automatically). Works on Linux and macOS.

```sh
bash run.sh setup
make test
```

### Option A: You already have a GEM

If you already have a model (SBML/YAML/JSON) and matching proteome FASTA, you can skip draft reconstruction and gap-filling and go straight to target screening.

1. Install DIAMOND for paralog clustering (or pass a precomputed `diamond_self_hits` table in your config):
   ```sh
   bash run.sh fetch-resources --profile screening
   ```
2. Copy [`configs/unified_pipeline.example.json`](configs/unified_pipeline.example.json) to `/work/run.local.json` and point it at your model, proteome, growth reaction, and target reaction:
   ```json
   {
     "localization_order": "skip",
     "gapfilled_model": "/work/model.xml",
     "proteome_fasta": "/work/proteome.fasta",
     "id_map": null,
     "biomass": "<your_growth_reaction_id>",
     "target_reaction": "<your_target_reaction_id>",
     "compute_growth_essentiality": true,
     "allow_unverified_reuse": true
   }
   ```
   Need to insert a heterologous pathway instead of an existing reaction? Configure `target_module` as shown in [Configuration & CLI Reference](docs/configuration.md).
3. Run the screening pipeline:
   ```sh
   bash run.sh metabolic --config /work/run.local.json --output /work/run-001
   ```

### Option B: Starting from a raw genome

Fresh AlphaGEM reconstruction runs on **Linux x86_64** and requires external databases plus separate AlphaGEM/DeepEC environments—see [Setup & External Resources](docs/setup.md) before running this path.

1. Prepare the proteome and coordinate files (from a GFF3, or via `annevo` / `helixer` / `tiberius` / `funannotate`):
   ```sh
   bash run.sh predict-proteome \
     --tool gff \
     --genome /work/raw/genome.fasta \
     --gff /work/raw/prediction.gff3 \
     --case-id sample_strain \
     --prefix PROT \
     --output /work/pred
   ```
2. Update `/work/run.local.json` with `/work/pred/inputs/` and your AlphaGEM paths, then run:
   ```sh
   bash run.sh metabolic --config /work/run.local.json --output /work/run-001
   ```

Check environment readiness or available commands at any time with `bash run.sh check-env` and `bash run.sh help`.

---

## Key Outputs

Under `/work/run-001/<order>/`:

| File | What it gives you |
| --- | --- |
| `fseof/all_genes_fseof_essential_filtered.tsv` | Non-essential overexpression candidates ranked by 10-step FSEOF flux slopes |
| `paralog/paralog_family_rankings_essential_filtered.tsv` | Homologous gene families clustered by DIAMOND and ranked by aggregate FSEOF score |
| `knockout/single_gene_ko_essential_filtered.tsv` | Single-gene knockout candidates (growth-essential genes removed) |
| `knockout/family_ko_essential_filtered.tsv` | Paralog-family joint knockout candidates (growth-essential genes removed) |
| `target_model_<order_slug>.xml` | Final SBML model with your target reaction or custom pathway mounted |

> The four `*_essential_filtered.tsv` tables require `"compute_growth_essentiality": true` (enabled in the minimal config and template above) or an external `essential_genes_tsv`; otherwise only the unfiltered ranking tables are written.

---

## Documentation

- **[Setup & External Resources](docs/setup.md)** — Platform support, one-click `fetch-resources` profiles, optional DeepLoc setup, and Linux AlphaGEM/DeepEC environments.
- **[Configuration & CLI Reference](docs/configuration.md)** — Custom product pathways (`target_module`), SHA-256 verification gates, step-by-step `stage` runs, and repository layout.
- **[Third-Party Software & Data](THIRD_PARTY.md)** · **[License (GPL-3.0)](LICENSE)**
