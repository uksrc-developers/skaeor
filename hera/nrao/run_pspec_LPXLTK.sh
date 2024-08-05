#!/bin/bash
#SBATCH --job-name test_pspec
#SBATCH --partition hera
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 15
#SBATCH --time 24:00:00
#SBATCH --mem 128GB
#SBATCH --output ./slurm-out/%J.out

source ~/.bashrc
conda activate h1c-idr3
cd /lustre/aoc/projects/hera/jburba/uksrc/hera/power-spectra/h1c-idr3/pspec

echo "start: $(date)"
pspec_dir=/lustre/aoc/projects/hera/jburba/uksrc/skaeor/hera
${pspec_dir}/../pspec_pipe.py ${pspec_dir}/pspec_params_LPXLTK_NRAO.yaml

echo "end: $(date)"
