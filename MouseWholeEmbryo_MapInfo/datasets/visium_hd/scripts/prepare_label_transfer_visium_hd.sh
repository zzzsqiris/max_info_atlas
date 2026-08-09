#!/bin/bash
# Prepare the existing WholeEmbryo label-transfer workflow for Visium HD.
# This creates normalized shared-gene features and UGE chunks; it does not
# compute kNN neighbors or transfer labels.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REPO_DIR="$(cd "$PIPELINE_DIR/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
PYTHONPATH="${REPO_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "$PYTHON_BIN" \
  "$PIPELINE_DIR/shared/scripts/prepare_label_transfer.py" \
  --reference "$PIPELINE_DIR/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad" \
  --target "$PIPELINE_DIR/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad" \
  --output-dir "$PIPELINE_DIR/datasets/visium_hd/results/label_transfer/work" \
  --expected-reference-cells 215284 \
  --expected-target-cells 449700 \
  --expected-shared-genes 0 \
  --target-cells-per-job 500 \
  --overwrite
