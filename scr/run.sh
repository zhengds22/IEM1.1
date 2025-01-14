#!/bin/bash
#SBATCH -J test
#SBATCH -N 1
#SBATCH -p cnmix
#SBATCH -o stdout1.%j
#SBATCH -e stderr1.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=56 #根据实际情况设置
#SBATCH -x ibc11b11n[01-38]

source /home/tongdan/anaconda3/etc/profile.d/conda.sh
#conda activate ccus
module load soft/gurobi/config

current_path=$(pwd)
cd "$current_path"

python core_model.py