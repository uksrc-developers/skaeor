#!/bin/bash
#SBATCH --job-name test_pspec
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 10
#SBATCH --time 24:00:00
#SBATCH --mem 128GB
#SBATCH --output ./slurm-out/%J.out

source ~/.bashrc
conda activate h1c-idr3
cd /project/power-spectra/h1c-idr3/subset/

echo "start: $(date)"
pspec_dir=/project/src/skaeor/hera/azimuth/subset
${pspec_dir}/../../pspec_pipe.py ${pspec_dir}/pspec_params_LPXLTK.yaml

echo "end: $(date)"
