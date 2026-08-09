#!/bin/bash
#$ -cwd
#$ -o MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/logs/prepare_visium_hd.$JOB_ID.log
#$ -j y
#$ -l h_data=64G,h_rt=8:00:00
#$ -pe shared 4
#$ -N visium_hd_prepare_transfer

# Hoffman2 job wrapper. This prepares shared-gene features and chunks only;
# it does not compute kNN neighbors or transfer labels.
set -euo pipefail

. /u/local/Modules/default/init/modules.sh
module load conda
conda activate xenium-mapinfo

export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"

mkdir -p MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/logs
bash MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/prepare_label_transfer_visium_hd.sh

python - <<'PY'
import json
import math
from pathlib import Path

import numpy as np

work = Path("MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/work")
with open(work / "manifest.json") as handle:
    manifest = json.load(handle)
reference = np.load(work / "reference_features.npy", mmap_mode="r")
target = np.load(work / "target_features.npy", mmap_mode="r")
assert manifest["source_reference_cells"] == 215284, manifest
assert manifest["target_cells"] == 449700, manifest
assert manifest["initial_shared_genes"] > 0, manifest
assert manifest["target_cells_for_knn"] + manifest["excluded_zero_variance_target_cells"] == 449700, manifest
assert manifest["jobs"] == math.ceil(manifest["target_cells_for_knn"] / manifest["target_cells_per_job"]), manifest
assert reference.shape[1] == target.shape[1] == manifest["retained_shared_genes"]
print("verified Visium HD prepared transfer:", reference.shape, target.shape)
print("initial shared genes:", manifest["initial_shared_genes"])
print("retained shared genes:", manifest["retained_shared_genes"])
print("full target cells:", manifest["target_cells"])
print("zero-variance target cells:", manifest["excluded_zero_variance_target_cells"])
print("jobs:", manifest["jobs"])
PY
