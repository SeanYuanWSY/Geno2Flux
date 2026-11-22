# DeepLoc Plug-and-Play Slot (`tools/deeploc/`)

DeepLoc 2.0 / 2.1 (DTU Health Tech) requires an individual academic license from DTU, so Geno2Flux does **not** bundle DeepLoc binaries or weights inside this repository.

Instead, Geno2Flux reserves this directory as a **zero-configuration drop-in slot** and automatically falls back when DeepLoc is absent:

## Official DTU Links

Request academic access directly from DTU Health Tech — Geno2Flux never redistributes DeepLoc:

- **DeepLoc 2.1**: [web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) · [academic download request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All)
- **DeepLoc 2.0**: [web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) · [academic download request form](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.0&packageversion=2.0&platform=All)

## 1. How Auto-Detection Works

When `localization_tsv` / `--localization-table` is not supplied, Geno2Flux checks for DeepLoc in this exact order:

1. Explicit `deeploc_cmd` in your JSON config (or `--deeploc-bin` on `run_stage.py localize`);
2. `G2F_DEEPLOC_BIN` environment variable (e.g. `export G2F_DEEPLOC_BIN=/opt/deeploc2/bin/deeploc2`);
3. Drop-in binary or symlink in **this directory**:
   - `tools/deeploc/bin/deeploc2` (e.g. a virtual environment or conda prefix installed at `tools/deeploc`)
   - `tools/deeploc/deeploc2` (e.g. a symlink or wrapper script)
4. `deeploc2` on your system `PATH`.

## 2. How to Plug DeepLoc In (Once You Receive Academic Access from DTU)

Pick whichever is easiest for your machine:

- **Option A — One-click installer** (recommended): once you receive the `DeepLoc-2.x.tar.gz` tarball from DTU, run
  ```sh
  bash run.sh fetch-resources --profile deeploc --deeploc-tarball /path/to/DeepLoc-2.1.tar.gz
  ```
  which installs it into this `tools/deeploc/` slot so `deeploc2` is auto-detected.
- **Option B — Install or symlink into this slot manually** (ignored by `.gitignore` so weights are never committed):
  ```bash
  # Example 1: symlink an existing deeploc2 executable here
  ln -s /path/to/your/env/bin/deeploc2 tools/deeploc/deeploc2

  # Example 2: create a dedicated venv right at tools/deeploc
  python3 -m venv tools/deeploc
  tools/deeploc/bin/pip install /path/to/DeepLoc-2.1.tar.gz
  ```
- **Option C — Export an environment variable**:
  ```bash
  export G2F_DEEPLOC_BIN=/path/to/deeploc2
  ```
- **Option D — Pass a precomputed DeepLoc CSV/TSV**:
  Set `"localization_tsv": "path/to/results_xxx.csv"` in your config JSON.

## 3. What Happens When DeepLoc Is Not Installed?

If DeepLoc is not found and no `localization_tsv` is given, Geno2Flux **does not fail**. It logs a clear notice, preserves the draft/template model's default compartments (`"deeploc_status": "fallback_default_compartments"` in `localization_refinement_audit.json` and `pipeline_manifest.json`), and continues seamlessly through gap-filling, FSEOF, paralog clustering, and knockout prioritization.

To enforce strict failure when DeepLoc is missing, set `"require_deeploc": true` in your config JSON (or pass `--require-deeploc` to `bash run.sh stage localize`).
