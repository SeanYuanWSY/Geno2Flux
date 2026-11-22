#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export UV_CACHE_DIR="${UV_CACHE_DIR:-$ROOT/.cache/uv}"
# Drop-in tool binaries (fetched via 'run.sh fetch-resources') take precedence.
if [ -d "$ROOT/tools/bin" ]; then
  PATH="$ROOT/tools/bin:$PATH"
fi
if [ -n "${G2F_TOOLS_BIN:-}" ] && [ -d "$G2F_TOOLS_BIN" ]; then
  PATH="$G2F_TOOLS_BIN:$PATH"
fi
export PATH

case "${1:-help}" in
  setup|setup-metabolic)
    cd "$ROOT"
    uv sync --frozen --extra test
    ;;
  predict-proteome)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test python "$ROOT/scripts/predict_proteome.py" "$@"
    ;;
  fetch-resources)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test python "$ROOT/scripts/fetch_external_resources.py" "$@"
    ;;
  check-env)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test python "$ROOT/scripts/fetch_external_resources.py" --check "$@"
    ;;
  fetch-reference-proteomes)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test python "$ROOT/scripts/fetch_reference_proteomes.py" "$@"
    ;;
  metabolic|pipeline)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test geno2flux pipeline "$@"
    ;;
  stage)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test python "$ROOT/scripts/run_stage.py" "$@"
    ;;
  prepare-inputs)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test python "$ROOT/scripts/build_coordinate_contract.py" "$@"
    ;;
  metabolic-tool|cli)
    shift
    cd "$ROOT"
    exec uv run --frozen --extra test geno2flux "$@"
    ;;
  help|--help|-h)
    printf '%s\n' \
      'bash run.sh setup                                 Sync locked Python dependencies via uv' \
      'bash run.sh predict-proteome [ARGS...]            Predict/standardize genome -> proteome + AlphaGEM contracts' \
      'bash run.sh fetch-reference-proteomes [ARGS...]   Fetch NCBI RefSeq / UniProt reference proteomes' \
      'bash run.sh fetch-resources [ARGS...]            Download external tools/resources into tools/bin (opt-in: --profile all)' \
      'bash run.sh check-env [ARGS...]                  Preflight-check interpreters, tools, and runtime resources' \
      'bash run.sh prepare-inputs [ARGS...]              Build AlphaGEM target & coordinate contracts from genome + GFF3 + FASTA' \
      'bash run.sh metabolic --config FILE --output DIR  Run unified Geno2Flux pipeline' \
      'bash run.sh stage <subcommand> [ARGS...]          Run an individual stage (localize/gapfill/target/fseof/paralogs/knockout/filter)'
    ;;
  *) printf 'Unknown command: %s\n' "$1" >&2; exit 2 ;;
esac
