#!/usr/bin/env bash
set -euo pipefail
BASE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"
PYTHON="${PYTHON:-python}"
if [ "${DEBUG:-0}" = "1" ]; then set -x; fi
log() { printf '[setup] %s\n' "$*"; }
foldseek_version_matches() {
  printf '%s\n' "$1" |
    grep -Eq '(^|[^[:alnum:].])9[.]427df8a([^[:alnum:].]|$)'
}
diamond_version_matches() {
  printf '%s\n' "$1" |
    grep -Eq '(^|[^[:alnum:].])2[.]1[.]10([^[:alnum:].]|$)'
}
trap 'echo "[setup] FAILED at line $LINENO"; exit 1' ERR
log "Starting setup in $BASE_DIR"
OS_RELEASE_FILE="${ALPHAGEM_OS_RELEASE_FILE:-/etc/os-release}"
if [ -r "$OS_RELEASE_FILE" ]; then
  . "$OS_RELEASE_FILE"
  if [ "${ID:-}" != "ubuntu" ]; then
    echo "Ubuntu only. Detected: ${ID:-unknown}" >&2
    exit 1
  fi
  log "OS check passed: ${ID:-unknown} ${VERSION_ID:-}"
else
  echo "Ubuntu only. OS release file not found: $OS_RELEASE_FILE" >&2
  exit 1
fi
if [ ! -d "./tools" ] || [ ! -d "./data_available" ] || [ ! -d "./plmsearchtools" ] || [ ! -d "./struct_data" ]; then
  if [ -f "data_source.tar.xz" ]; then
    if [ ! -f "data_source.tar.xz.sha256" ]; then
      echo "data_source.tar.xz.sha256 not found; refusing unverified extraction" >&2
      exit 1
    fi
    log "Verifying data_source.tar.xz"
    sha256sum --check data_source.tar.xz.sha256
    log "Extracting data_source.tar.xz"
    tar -xJf data_source.tar.xz
    log "Extraction complete"
  else
    log "data_source.tar.xz not found, skip extraction"
  fi
else
  log "Data directories already present, skip extraction"
fi
INSTALL_FOLDSEEK=0
if command -v foldseek >/dev/null 2>&1; then
  FOLDSEEK_VERSION="$(foldseek version 2>&1)"
  if foldseek_version_matches "$FOLDSEEK_VERSION"; then
    log "Foldseek version verified: $FOLDSEEK_VERSION"
  else
    log "Replacing incompatible Foldseek: $FOLDSEEK_VERSION"
    INSTALL_FOLDSEEK=1
  fi
else
  INSTALL_FOLDSEEK=1
fi
if [ "$INSTALL_FOLDSEEK" = "1" ]; then
  if command -v conda >/dev/null 2>&1; then
    log "Installing Foldseek 9.427df8a via conda"
    conda install -y -c conda-forge -c bioconda foldseek=9.427df8a
    log "foldseek install finished"
  else
    echo "conda not found; cannot install Foldseek 9.427df8a" >&2
    exit 1
  fi
fi
if ! command -v foldseek >/dev/null 2>&1; then
  echo "Foldseek is unavailable after installation" >&2
  exit 1
fi
FOLDSEEK_VERSION="$(foldseek version 2>&1)"
if ! foldseek_version_matches "$FOLDSEEK_VERSION"; then
  echo "Foldseek version mismatch after installation: $FOLDSEEK_VERSION" >&2
  exit 1
fi

INSTALL_DIAMOND=0
if command -v diamond >/dev/null 2>&1; then
  DIAMOND_VERSION="$(diamond version 2>&1)"
  if diamond_version_matches "$DIAMOND_VERSION"; then
    log "DIAMOND version verified: $DIAMOND_VERSION"
  else
    log "Replacing incompatible DIAMOND: $DIAMOND_VERSION"
    INSTALL_DIAMOND=1
  fi
else
  INSTALL_DIAMOND=1
fi
if [ "$INSTALL_DIAMOND" = "1" ]; then
  if command -v conda >/dev/null 2>&1; then
    log "Installing DIAMOND 2.1.10 via conda"
    conda install -y -c conda-forge -c bioconda diamond=2.1.10
    log "diamond install finished"
  else
    echo "conda not found; cannot install DIAMOND 2.1.10" >&2
    exit 1
  fi
fi
if ! command -v diamond >/dev/null 2>&1; then
  echo "DIAMOND is unavailable after installation" >&2
  exit 1
fi
DIAMOND_VERSION="$(diamond version 2>&1)"
if ! diamond_version_matches "$DIAMOND_VERSION"; then
  echo "DIAMOND version mismatch after installation: $DIAMOND_VERSION" >&2
  exit 1
fi
if [ -d "./tools/eggnog_mapper" ]; then
  if [ ! -d "./tools/eggnog_mapper/data" ] || [ ! -f "./tools/eggnog_mapper/data/eggnog.db" ]; then
    log "Downloading eggNOG data"
    pushd ./tools/eggnog_mapper >/dev/null
    yes | "$PYTHON" download_eggnog_data.py
    popd >/dev/null
    log "eggNOG data ready"
  else
    log "eggNOG data already present"
  fi
else
  log "eggnog_mapper directory not found, skip"
fi
if [ -d "./tools/CLEAN/app" ]; then
  log "Installing CLEAN from an isolated copy"
  (
    CLEAN_SETUP_COPY="$(mktemp -d)"
    trap 'rm -rf -- "$CLEAN_SETUP_COPY"' EXIT
    mkdir "$CLEAN_SETUP_COPY/app"
    cp -a ./tools/CLEAN/app/. "$CLEAN_SETUP_COPY/app/"
    cd "$CLEAN_SETUP_COPY/app"
    "$PYTHON" build.py install
  )
  if "$PYTHON" - <<'PY'
import sys, importlib.util
sys.exit(0 if (importlib.util.find_spec('CLEAN') or importlib.util.find_spec('clean')) else 1)
PY
  then
    log "CLEAN import available"
  else
    log "CLEAN import not found after install"
  fi
else
  log "CLEAN directory not found, skip"
fi
if [ -d "./tools/OrthoFinder" ]; then
  log "Installing OrthoFinder from an isolated copy"
  (
    ORTHOFINDER_SETUP_COPY="$(mktemp -d)"
    trap 'rm -rf -- "$ORTHOFINDER_SETUP_COPY"' EXIT
    mkdir "$ORTHOFINDER_SETUP_COPY/OrthoFinder"
    cp -a ./tools/OrthoFinder/. "$ORTHOFINDER_SETUP_COPY/OrthoFinder/"
    cd "$ORTHOFINDER_SETUP_COPY/OrthoFinder"
    "$PYTHON" setup.py install
  )
else
  log "OrthoFinder directory not found, skip"
fi
log "Setup completed successfully"
