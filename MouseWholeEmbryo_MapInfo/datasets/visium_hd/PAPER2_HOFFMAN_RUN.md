# Paper 2 Visium HD: Hoffman2 preparation and run commands

This note prepares Paper 2 for the existing WholeEmbryo workflow. It does not
run reference clustering, label transfer, ST-only Map Information, imputation,
or any local heavy preprocessing.

## Local files verified

- `MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad`
  - 215,284 reference cells; read-only AnnData open succeeded.
- `MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy`
  - 215,284 entries; one section: `Qiu_E15.5_embryo_59`.
- `MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad`
  - 449,700 spatial cells, 33,670 genes; sparse CSR AnnData.
- `MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/Sections.npy`
  - 449,700 entries; one section: `MouseEmbryo_VisiumHD`.
- `MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/xy_coordinates/MouseEmbryo_VisiumHD_XY.npy`
  - 449,700 × 2 coordinates; aligned to the Visium AnnData row order.

The large shared-gene feature arrays and kNN chunks are intentionally not
generated locally. They are created by the Hoffman preparation job below.

## Files to upload to the matching Hoffman repository

The required data inputs are:

```text
MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad
MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy
MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad
MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/Sections.npy
MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/xy_coordinates/MouseEmbryo_VisiumHD_XY.npy
```

Upload these changed/required workflow files if the Hoffman checkout does not
already contain the same versions:

```text
MouseWholeEmbryo_MapInfo/reference/whole_embryo/config/e15_5_whole_embryo_reference_clustering.yaml
MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/build_reference_features_job.sh
MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/prepare_reference_clustering_sections.py
MouseWholeEmbryo_MapInfo/shared/scripts/prepare_label_transfer.py
MouseWholeEmbryo_MapInfo/shared/scripts/run_label_transfer_chunk.py
MouseWholeEmbryo_MapInfo/shared/scripts/submit_label_transfer.py
MouseWholeEmbryo_MapInfo/shared/scripts/check_label_transfer.py
MouseWholeEmbryo_MapInfo/shared/scripts/combine_label_transfer.py
MouseWholeEmbryo_MapInfo/shared/scripts/prepare_st_only_full_mapinfo.py
MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/prepare_label_transfer_visium_hd.sh
MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/prepare_label_transfer_visium_hd_job.sh
MouseWholeEmbryo_MapInfo/datasets/visium_hd/config/visium_hd_st_only_map_info.yaml
```

The raw Visium archives are not needed for Label Transfer or ST-only Map
Information once the converted AnnData and spatial arrays have been uploaded.
Do not upload `results/label_transfer/work`, `chunk_results`, or Map
Information results from this local checkout; they do not exist yet.

From the local repository root, an explicit upload is:

```bash
LOCAL=/Users/shiqizhong/Code/max_info_atlas
REMOTE=sqzhong@hoffman2.idre.ucla.edu:~/project-rwollman/max_info_atlas

ssh sqzhong@hoffman2.idre.ucla.edu \
  'mkdir -p ~/project-rwollman/max_info_atlas/MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed \
    ~/project-rwollman/max_info_atlas/MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/xy_coordinates \
    ~/project-rwollman/max_info_atlas/MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts \
    ~/project-rwollman/max_info_atlas/MouseWholeEmbryo_MapInfo/datasets/visium_hd/config'

rsync -av --progress \
  "$LOCAL/MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad" \
  "$LOCAL/MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy" \
  "$REMOTE/MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/"

rsync -av --progress \
  "$LOCAL/MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad" \
  "$REMOTE/MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/"

rsync -av --progress \
  "$LOCAL/MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/Sections.npy" \
  "$LOCAL/MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/xy_coordinates/MouseEmbryo_VisiumHD_XY.npy" \
  "$REMOTE/MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/"

(cd "$LOCAL" && rsync -av --progress --relative \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/config/e15_5_whole_embryo_reference_clustering.yaml \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/build_reference_features_job.sh \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/prepare_reference_clustering_sections.py \
  MouseWholeEmbryo_MapInfo/shared/scripts/prepare_label_transfer.py \
  MouseWholeEmbryo_MapInfo/shared/scripts/run_label_transfer_chunk.py \
  MouseWholeEmbryo_MapInfo/shared/scripts/submit_label_transfer.py \
  MouseWholeEmbryo_MapInfo/shared/scripts/check_label_transfer.py \
  MouseWholeEmbryo_MapInfo/shared/scripts/combine_label_transfer.py \
  MouseWholeEmbryo_MapInfo/shared/scripts/prepare_st_only_full_mapinfo.py \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/prepare_label_transfer_visium_hd.sh \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/prepare_label_transfer_visium_hd_job.sh \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/config/visium_hd_st_only_map_info.yaml \
  "$REMOTE")
```

If the remote checkout is known to be the same code version, upload only the
listed files; do not replace the whole repository.

## Hoffman2 commands

Run these from:

```bash
cd ~/project-rwollman/max_info_atlas
. /u/local/Modules/default/init/modules.sh
module load conda
conda activate xenium-mapinfo
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
```

### 1. Verify the uploaded inputs

```bash
ls -lh \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/Sections.npy \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/xy_coordinates/MouseEmbryo_VisiumHD_XY.npy

python - <<'PY'
import anndata as ad
from pathlib import Path

paths = [
    Path('MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad'),
    Path('MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad'),
]
for path in paths:
    a = ad.read_h5ad(path, backed='r')
    try:
        print(path, a.shape, getattr(a.X, 'format', None))
    finally:
        a.file.close()
PY
```

### 2. Verify or rebuild reference clustering

The reorganization moves the existing Hoffman2 reference-clustering directory
into the path below. Count the label arrays first. If the count is `50`, reuse
that completed reference result and do not resubmit the clustering jobs.

```bash
find MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/clustering/LeidenPCA50Correlation \
  -path '*/Qiu_E15.5_embryo_59.npy' -type f | wc -l
```

Only use the rebuild commands below if the existing reference result is
missing or incomplete and a rebuild is intentionally required.

The reference sections file is already supplied. If it ever needs to be
rebuilt, use the existing validator before `max-info run prepare`:

```bash
python MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/prepare_reference_clustering_sections.py
```

Prepare the reference job lists, then submit the existing checked feature job:

```bash
max-info run prepare \
  --config MouseWholeEmbryo_MapInfo/reference/whole_embryo/config/e15_5_whole_embryo_reference_clustering.yaml

qsub MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/build_reference_features_job.sh
```

After the feature job completes and its validation message is present:

```bash
max-info run submit \
  --config MouseWholeEmbryo_MapInfo/reference/whole_embryo/config/e15_5_whole_embryo_reference_clustering.yaml \
  --step graphs

max-info run submit \
  --config MouseWholeEmbryo_MapInfo/reference/whole_embryo/config/e15_5_whole_embryo_reference_clustering.yaml \
  --step clustering
```

Verify the 50 reference label arrays before starting target transfer:

```bash
find MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/clustering/LeidenPCA50Correlation \
  -path '*/Qiu_E15.5_embryo_59.npy' -type f | wc -l
```

Expected count: `50`.

### 3. Prepare Visium HD Label Transfer inputs

This is the first Paper 2-specific qsub job. It creates normalized shared-gene
feature arrays and target chunks, but does not calculate neighbors:

```bash
mkdir -p MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/logs
qsub MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/prepare_label_transfer_visium_hd_job.sh
```

The Visium wrapper intentionally records the observed reference/target gene
intersection instead of assuming the old 290-gene MERFISH intersection. It
still uses the canonical `normalize_total(10000)`, `log1p`, correlation metric,
`k=15`, and zero-variance source-index handling.

### 4. Submit and validate target kNN label transfer

Preview first; add `--submit` only when ready:

```bash
python MouseWholeEmbryo_MapInfo/shared/scripts/submit_label_transfer.py \
  --chunks-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/work/chunks \
  --work-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/work \
  --results-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/chunk_results \
  --logs-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/logs \
  --temp-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/temp \
  --job-name visium_hd_label_transfer \
  --memory 8G \
  --runtime 12:00:00 \
  --conda-env xenium-mapinfo
```

Submit the same command with `--submit`, then check all chunks:

```bash
python MouseWholeEmbryo_MapInfo/shared/scripts/check_label_transfer.py \
  --work-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/work \
  --results-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/chunk_results \
  --k 15
```

After the checker reports 100% complete, combine the 50 reference resolutions
into a Visium target AnnData:

```bash
python MouseWholeEmbryo_MapInfo/shared/scripts/combine_label_transfer.py \
  --target MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad \
  --labels-dir MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/clustering/LeidenPCA50Correlation \
  --work-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/work \
  --chunks-dir MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/chunk_results \
  --output MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_label_transfer.h5ad \
  --summary MouseWholeEmbryo_MapInfo/datasets/visium_hd/results/label_transfer/label_transfer_summary.csv \
  --expected-source-reference-cells 215284 \
  --expected-target-cells 449700 \
  --reference-section Qiu_E15.5_embryo_59 \
  --k 15 \
  --overwrite
```

This document stops before imputation and Map Information execution.
