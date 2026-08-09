# One-time Hoffman2 reorganization

These commands reorganize the currently verified Hoffman2 directories without
permanently deleting files. Run them from a local terminal. They deliberately
move the existing reference clustering results into the active shared
reference area before archiving the rejected MERFISH workflow.

## 1. Reorganize existing Hoffman2 directories

```bash
ssh sqzhong@hoffman2.idre.ucla.edu '
set -euo pipefail
cd ~/project-rwollman/max_info_atlas

test ! -e MouseWholeEmbryo_MapInfo
test -d LabelTransfer
test -d LabelTransfer_E15.5
test -d LabelTransfer_E15.5_WholeEmbryo
test -d MouseWholeEmbryo_VisiumHD
test -f LabelTransfer_E15.5_WholeEmbryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad
test -f LabelTransfer_E15.5_WholeEmbryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy

mkdir -p \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/results \
  MouseWholeEmbryo_MapInfo/datasets \
  MouseWholeEmbryo_MapInfo/shared \
  MouseWholeEmbryo_MapInfo/archive

mv LabelTransfer_E15.5_WholeEmbryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/
mv LabelTransfer_E15.5_WholeEmbryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/

if [ -d LabelTransfer_E15.5_WholeEmbryo/results/e15_5_whole_embryo_reference_clustering ]; then
  mv LabelTransfer_E15.5_WholeEmbryo/results/e15_5_whole_embryo_reference_clustering \
    MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering
fi

mv MouseWholeEmbryo_VisiumHD MouseWholeEmbryo_MapInfo/datasets/visium_hd
mkdir -p \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/data/raw \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/data/processed \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/scripts \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/config \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/results \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/logs \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/figures \
  MouseWholeEmbryo_MapInfo/datasets/stereo_seq/tmp
mv LabelTransfer MouseWholeEmbryo_MapInfo/archive/starmap
mv LabelTransfer_E15.5 MouseWholeEmbryo_MapInfo/archive/merfish_ec_only
mv LabelTransfer_E15.5_WholeEmbryo MouseWholeEmbryo_MapInfo/archive/merfish_whole_embryo

find MouseWholeEmbryo_MapInfo -maxdepth 3 -type d -print | sort
du -sh MouseWholeEmbryo_MapInfo/reference MouseWholeEmbryo_MapInfo/datasets/* MouseWholeEmbryo_MapInfo/archive/*
'
```

If any preflight `test` fails, the command stops before creating or moving
anything. Do not rerun it blindly after a partial migration; inspect the remote
tree first.

## 2. Upload the reorganized workflow files

The large reference and Visium inputs already exist on Hoffman2 and are moved
in place by step 1. Upload only the reorganized code, configs, and documents:

```bash
cd /Users/shiqizhong/Code/max_info_atlas

rsync -av --progress --relative \
  MouseWholeEmbryo_MapInfo/README.md \
  MouseWholeEmbryo_MapInfo/HOFFMAN_REORGANIZE.md \
  MouseWholeEmbryo_MapInfo/shared/scripts/ \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/config/ \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/ \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/config/ \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/scripts/ \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/PAPER2_HOFFMAN_RUN.md \
  sqzhong@hoffman2.idre.ucla.edu:~/project-rwollman/max_info_atlas/
```

This command does not use `--delete`, so it cannot remove remote results.

## 3. Verify the remote active inputs

```bash
ssh sqzhong@hoffman2.idre.ucla.edu '
set -euo pipefail
cd ~/project-rwollman/max_info_atlas

ls -lh \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad \
  MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_sections.npy \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/Sections.npy \
  MouseWholeEmbryo_MapInfo/datasets/visium_hd/data/processed/visium_hd_st_only_mapinfo/xy_coordinates/MouseEmbryo_VisiumHD_XY.npy

find MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/clustering/LeidenPCA50Correlation \
  -path "*/Qiu_E15.5_embryo_59.npy" -type f | wc -l
'
```

The final count should be `50` if the previously completed reference
clustering was moved successfully. This migration does not submit jobs or run
preprocessing.
