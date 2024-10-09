# H1C IDR3 Power Spectrum Pipeline Containerization

**This page is under construction.**

## Building the container

The container can be built from the provided `apptainer` definition file via
```
sudo apptainer build hera-pspec.sif hera-pspec-mambaorg.def | tee build_log.txt
```
The use of the pipe to `tee` in this case saves the stdout from the build process to a text file for verification of installed python packages and versions without having to access the container directly.

## Interacting with the built container

To interact with the built singularity container in a bash shell, use the following command
```
apptainer shell  --shell /bin/bash hera-pspec.sif
```

## Running the container

After building the container, the two other files of interest are

- `pspec_params_LPXLTK.yaml`: yaml file containing file paths and analysis parameters for `hera_pspec`
- `run_pspec_LPXLTK.sh`: slurm sbatch script which calls the power spectrum code, `pspec_pipe.py` (this script is already built into the container at `/opt/hera/pspec_pipe.py`)

The analysis of the full set of H1C IDR3 data using the container can then be run (on Azimuth) via
```
sbatch run_pspec_LPXLTK.sh --bind /project/HERA_data:/data --bind /project/power-spectra/h1c-idr3/full hera-pspec.sif pspec_params_LPXLTK.yaml
```
where the two `--bind` calls in this case do the following:

1. `--bind /project/HERA_data:/data` binds the directory containing the H1C IDR3 data to `/data` inside the container.  On Azimuth, these H1C IDR3 data are stored in `/project/HERA_data`.

2. `--bind /project/power-spectra/h1c-idr3/full` binds the output directory for log and data files to the matching path inside the container.  If you wish to change the output directory, you must change `work_dir` and `out_dir` arguments nested under `io` in `pspec_params_LPXLTK.yaml` and replace `/project/power-spectra/h1c-idr3/full` with the desired output path.

This analysis was designed to run on an Azimuth Slurm platform with a single `memory.medium` compute node with 54 CPUs and 172 GB RAM.  Hence, by default, `run_pspec_LPXLTK.sh` requests 54 CPUs and 160 GB RAM.  These quanities can be adjusted, but decreasing the number of CPUs will increase the run time.

By default, the output files from `hera_pspec` will be stored in the `ska-teal-eor` tenancy on Azimuth at
```
/project/power-spectra/h1c-idr3/full/
```
The files generated within this directory are described in the table below.

| File Name | File Size | Description |
| --------- | --------- | ----------- |
| `pspec.grp1.of1.LPXLTK.h5` | --- | `hera_pspec` HDF5 file containing per-baseline delay power spectra |
| `pspipe_out_LPXLTK.log` | --- | Log file containing input parameter values and execution timing statements |
| `zen.grp1.of1.autos.Tsys.LPXLTK.uvh5` | --- | `pyuvdata` file containing noise spectra (used to calculate power spectrum uncertainties) |

## UKSRC related links and information

### Confluence pages

- [Containerizing the power spectrum pipeline](https://confluence.skatelescope.org/x/ojazEQ)

### Jira tickets

- [TEAL-711](https://jira.skatelescope.org/browse/TEAL-711): containerize HERA step 7 (power spectrum) pipeline
