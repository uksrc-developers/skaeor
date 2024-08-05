#!/bin/bash

#source ~/.bashrc
#conda activate h1c-idr3
which python

cd /project/power-spectra/h1c-idr3/subset/

echo "start: $(date)" 
pspec_dir=/project/src/skaeor/hera/azimuth
${pspec_dir}/../pspec_pipe.py ${pspec_dir}/pspec_params_LPXLTK.yaml

echo "end: $(date)"
