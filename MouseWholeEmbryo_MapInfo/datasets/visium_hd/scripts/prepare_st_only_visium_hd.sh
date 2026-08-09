#!/bin/bash
# Export Sections.npy and spatial coordinates for the existing ST-only workflow.
# This does not run max-info or build Map Information results.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REPO_DIR="$(cd "$PIPELINE_DIR/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
PYTHONPATH="${REPO_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "$PYTHON_BIN" \
  "$PIPELINE_DIR/shared/scripts/prepare_st_only_full_mapinfo.py" \
  --input "$PIPELINE_DIR/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad" \
  --output-dir "$PIPELINE_DIR/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo" \
  --section MouseEmbryo_VisiumHD \
  --expected-cells 449700 \
  --overwrite
