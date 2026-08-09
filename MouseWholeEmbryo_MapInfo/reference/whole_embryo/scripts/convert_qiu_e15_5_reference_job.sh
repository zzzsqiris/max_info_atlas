#!/bin/bash
#$ -cwd
#$ -o MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/conversion/logs/convert_qiu_e15_5_reference.$JOB_ID.log
#$ -j y
#$ -l h_data=24G,h_rt=24:00:00
#$ -pe shared 4
#$ -N qiu_e15_5_convert

set -euo pipefail

. /u/local/Modules/default/init/modules.sh
module load conda
conda activate xenium-mapinfo

export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"

python MouseWholeEmbryo_MapInfo/reference/whole_embryo/scripts/convert_qiu_e15_5_reference.py
