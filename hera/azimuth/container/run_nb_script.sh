#!/bin/bash
#SBATCH --job-name test_notebook
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 10
#SBATCH --time 59-00:00:00
#SBATCH --mem 152G
#SBATCH --output ./slurm-out/%J.out

source ~/.bashrc
mamba activate h1c-idr3
which python
echo

cd /project/src/skaeor/hera/azimuth/container
pwd
echo

python -u /project/src/skaeor/hera/ps_notebook.py \
    --HERA_data_dir /project/HERA_data/ \
    /project/power-spectra/h1c-idr3/full/pspec.grp1.of1.LPXLTK.h5
