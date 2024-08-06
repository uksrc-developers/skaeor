# Azimuth tests

After successfully running the `hera_pspec` pipeline on NRAO, we transferred the same subset of data used on NRAO to Azimuth and ran a similar test on a general.medium Linux with ssh platform.  So far, we have only tested the `hera_pspec` step, the jupyter notebook step will follow in the future.  This test is described in the following subsections.  More information is available on the [H1C IDR3 Power Spectra](https://confluence.skatelescope.org/display/SRCSC/H1C+IDR3+Power+Spectra) confluence page.

## Run `hera_pspec`

In this test, JB ran `hera_pspec` on a subset of the full dataset used by HERA to ensure the python environment can run `hera_pspec` and generate the files that get fed into the notebook.  JB used the same subset of the full HERA data (4 out of 20 files) used in the test on NRAO (see the [README](https://github.com/uksrc-developers/skaeor/blob/main/hera/nrao/README.md) in the `../nrao/` directory and/or the [H1C IDR3 Power Spectra](https://confluence.skatelescope.org/display/SRCSC/H1C+IDR3+Power+Spectra) confluence page for more details).  For this test the two files of interest in this directory are

- `pspec_params_LPXLTK.yaml`: yaml file containing file paths and analysis parameters for `hera_pspec`.  The changes made in this file to reduce the frequency and baseline axes are marked with in-line comments.  The parameter values are identical to that used in `../nrao/` with the exception of the file paths which have been updated to reflect the location of the transferred data on Azimuth.
- `../pspec_pipe.py`: python file which reads in the parameters in `pspec_params_LPXLTK.yaml` and runs `hera_pspec`
- `run_pspec_LPXLTK.sh`: bash script which calls `../pspec_pipe.py`

At the time of this testing, all of the `general.medium` workstations on Azimuth were in use, so JB could not create a slurm platform with a `general.medium` compute node.  The RAM on a `general.small` compute node, at 4 GB, was insufficient for this test.  The identical test as run on NRAO in `../nrao/` required ~5 GB of RAM, for reference.

These files were used to run `hera_pspec` on Azimuth via the following steps: 

1. Start a tmux instance
    ```
    tmux new -s h1c-idr3
    ```
    JB used `tmux` to avoid losing any work as this test was run on a Linux with ssh platform over an ssh connection which could potentially fail while the code was running.

2. Activate the `h1c-idr3` python environment
    ```
    mamba activate h1c-idr3
    ```
    Please note that when using a new Azimuth platform, `mamba` must first be initialized before the above command will work. 
 Please see the section titled "Using the python environment on Azimuth" on the [Generating IDR3 power spectra on Azimuth](https://confluence.skatelescope.org/display/SRCSC/Generating+IDR3+power+spectra+on+Azimuth) confluence page for more details.

3. Navigate to the cloned `uksrc-developers/skaeor` repo
    ```
    cd /project/src/skaeor/hera/azimuth
    ```

4. Run the `hera_pspec` pipeline and redirect the stdout to a file in the created directory `/project/power-spectra/h1c-idr3/subset` via
    ```
    ./run_pspec_LPXLTK.sh > /project/power-spectra/h1c-idr3/subset/stdout.txt
    ```

# UKSRC related links and information

**Confluence pages**

- [H1C IDR3 Power Spectra](https://confluence.skatelescope.org/display/SRCSC/H1C+IDR3+Power+Spectra)
- [Generating IDR3 power spectra on Azimuth](https://confluence.skatelescope.org/display/SRCSC/Generating+IDR3+power+spectra+on+Azimuth)

**Jira tickets**

- [TEAL-662](https://jira.skatelescope.org/browse/TEAL-662): transfer data subset from NRAO
- [TEAL-671](https://jira.skatelescope.org/browse/TEAL-671): set up and test python environment on Azimuth
- [TEAL-672](https://jira.skatelescope.org/browse/TEAL-672): run power spectrum pipeline (subset of data, Azimuth)
