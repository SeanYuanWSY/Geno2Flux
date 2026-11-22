# External CLI Binary Slot (`tools/bin/`)

This directory is a **gitignored, auto-detected drop-in slot** for external
command-line tools used by Geno2Flux: `diamond` (sequence screening and
`fungi.dmnd` construction) plus `tblastn`, `makeblastdb`, and `blastp` (NCBI
BLAST+ gap-filling). Only this README is committed.

## How it is populated

```bash
# DIAMOND v2.1.10 (Linux x86_64 tarball; conda/brew instructions on macOS)
bash run.sh fetch-resources --profile screening

# DIAMOND + NCBI BLAST+ 2.16.0 platform-matched binaries
bash run.sh fetch-resources --profile gapfill
```

Downloaded archives are verified against official MD5 checksums before any
binary is moved into place.

## How it is detected

`run.sh` prepends this directory to `PATH` automatically, and
`bash run.sh check-env` resolves binaries in this order:

1. `--bin-dir` (explicit override);
2. `G2F_TOOLS_BIN` environment variable;
3. `tools/bin/` (this directory);
4. system `PATH`.

Nothing else in this directory is tracked by git, so downloaded binaries are
never committed.
