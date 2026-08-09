#!/bin/bash
#$ -cwd
#$ -o MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/logs/features/manual_reference_features.$JOB_ID.log
#$ -j y
#$ -l h_data=48G,h_rt=24:00:00
#$ -pe shared 4
#$ -N qiu_e15_5_features

set -euo pipefail

. /u/local/Modules/default/init/modules.sh
module load conda
conda activate xenium-mapinfo

export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"

chunk_file="MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/jobs/features/chunks/chunk_0001_of_0001.txt"
output_file="MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/features/features_pca50.npy"

if [ ! -f "$chunk_file" ]; then
    echo "Missing feature chunk: $chunk_file"
    exit 1
fi

max-info run-features --chunk-file "$chunk_file"

if [ ! -f "$output_file" ]; then
    echo "Feature worker ended without creating: $output_file"
    exit 1
fi

python -c "import numpy as np; p='$output_file'; x=np.load(p, mmap_mode='r'); assert x.shape == (215284, 50), x.shape; assert x.dtype == np.float32, x.dtype; print('verified PCA features:', x.shape, x.dtype)"
