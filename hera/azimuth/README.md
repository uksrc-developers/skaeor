# Azimuth tests

After successfully running the `hera_pspec` pipeline on NRAO, we transferred the same subset of data used on NRAO to Azimuth and ran a similar test on a general.medium Linux with ssh platform.  So far, we have only tested the `hera_pspec` step, the jupyter notebook step will follow in the future (please see the [H1C IDR3 power spectrum pipeline](https://confluence.skatelescope.org/display/SRCSC/H1C+IDR3+power+spectrum+pipeline) confluence page for more details).  This test is described in the following subsections.

## python environment

Please note that a python environment has already been built in the `/project` space within the `ska-teal-eor` tenancy on Azimuth.  To use this python environment on a fresh Azimuth platform, you must first execute the following command
```
tail -n +2 /project/src/activate_mamba.sh >> ~/.bashrc
```
Close your connection and reconnect to your Azimuth platform and `mamba` will be initialized.  Initializing `mamba` is required before the sbatch script, `run_pspec_LPXLTK.sh`, will work properly.

Please see the section titled "Using the python environment on Azimuth" on the [Generating IDR3 power spectra on Azimuth](https://confluence.skatelescope.org/display/SRCSC/Generating+IDR3+power+spectra+on+Azimuth) confluence page for more details.

## Run `hera_pspec`

In this test, JB ran `hera_pspec` on a subset of the full dataset used by HERA to ensure the python environment can run `hera_pspec` and generate the files that get fed into the notebook.  JB used the same subset of the full HERA data (4 out of 20 files) used in the test on NRAO (see the [README](https://github.com/uksrc-developers/skaeor/blob/main/hera/nrao/README.md) in the `../nrao/` directory and/or the [H1C IDR3 Power Spectra](https://confluence.skatelescope.org/display/SRCSC/H1C+IDR3+Power+Spectra) confluence page for more details).  For this test the three files of interest are

- `pspec_params_LPXLTK.yaml`: yaml file containing file paths and analysis parameters for `hera_pspec`.  The changes made in this file to reduce the frequency and baseline axes are marked with in-line comments.  The parameter values are identical to that used in `../nrao/` with the exception of the file paths which have been updated to reflect the location of the transferred data on Azimuth.
- `../pspec_pipe.py`: python file which reads in the parameters in `pspec_params_LPXLTK.yaml` and runs `hera_pspec`
- `run_pspec_LPXLTK.sh`: slurm sbatch script which calls `../pspec_pipe.py`

These files were used to run `hera_pspec` on Azimuth via the following steps: 

1. Create a directory for the output from slurm via
    ```
    mkdir slurm-out
    ```

2. Submit the sbatch script to slurm via
    ```
    sbatch run_pspec_LPXLTK.sh
    ```
Relevant compute info can be found in the table below.  The max RAM usage was obtained via Slurm's MaxRSS.

| Input data (GB) | CPUs | Max RAM (GB) | Duration (HH:MM:SS:) | Output data |
| --------------- | ---- | ------------ | -------------------- | ----------- |
| 6.4             | 10   | 8.61         | 21:07:35             | 6.7         |

The output files from `hera_pspec` are stored in the `ska-teal-eor` tenancy on Azimuth at
```
/project/power-spectra/h1c-idr3/subset
```
The files generated within this directory are described in the table below.

| File Name | File Size | Description |
| --------- | --------- | ----------- |
| `pspec.grp1.of1.LPXLTK.h5` | 6.6 GB | `hera_pspec` HDF5 file containing per-baseline delay power spectra |
| `pspipe_out_LPXLTK.log` | 33 KB | Log file containing input parameter values and execution timing statements |
| `zen.grp1.of1.autos.Tsys.LPXLTK.uvh5` | 66 MB | `pyuvdata` file containing noise spectra (used to calculate power spectrum uncertainties) |

# UKSRC related links and information

**Confluence pages**

- [H1C IDR3 power spectrum pipeline](https://confluence.skatelescope.org/display/SRCSC/H1C+IDR3+power+spectrum+pipeline)
- [Generating IDR3 power spectra on Azimuth](https://confluence.skatelescope.org/display/SRCSC/Generating+IDR3+power+spectra+on+Azimuth)

**Jira tickets**

- [TEAL-662](https://jira.skatelescope.org/browse/TEAL-662): transfer data subset from NRAO
- [TEAL-671](https://jira.skatelescope.org/browse/TEAL-671): set up and test python environment on Azimuth
- [TEAL-672](https://jira.skatelescope.org/browse/TEAL-672): run power spectrum pipeline (subset of data, Azimuth)
- [TEAL-679](https://jira.skatelescope.org/browse/TEAL-679): run power spectrum pipeline (subset of data, Azimuth) [Continued from TEAL-672]
