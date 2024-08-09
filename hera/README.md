# HERA Analysis Pipeline

This directory contains information and scripts pertaining to the Hydrogen Epoch of Reionization Array ([HERA](https://reionization.org/)) analysis pipeline.

## HERA data naming convention

HERA labels its seasonal data according the the observing season and internal data release number.  For example, the latest published HERA results come from the H1C IDR3 dataset.  In this example, H1C refers to HERA's first observing season and IDR3 refers to internal data release 3.

## Steps in the full HERA pipeline

The full HERA pipeline contains the following steps

1. Calibration
2. Radio Frequency Interference (RFI) flagging
3. Nigh-to-night LST binning
4. Delay inpainting
5. Systematics subtraction
6. Time averaging
7. Power spectrum estimation
8. Astrophysical parameter estimation

**For now, this demonstrator case is only focusing on the delay power spectrum estimation step.** 

In the future, more analysis steps will be added and documented.

## Power spectrum estimation

For a detailed description of the HERA delay power spectrum pipeline see

- [HERA Collaboration 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...925..221A/abstract): HERA's first power spectrum results which details the full HERA analysis pipeline and results from H1C IDR2 data.
- [HERA Collaboration 2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...945..124H/abstract): HERA's latest power spectrum results using H1C IDR3 data.  This paper only details changes to the analysis that have been made from the 2022 paper using H1C IDR2 data.

**This repository focuses on reproducing the H1C IDR3 power spectrum results from [HERA Collaboration 2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...945..124H/abstract).**

### Power spectrum substeps

There are two substeps which fall under the "power spectrum estimation" umbrella for HERA.

1. Run [`hera_pspec`](https://github.com/HERA-Team/hera_pspec) to calculated per-baseline delay power spectra
2. Run a [jupyter notebook](https://github.com/HERA-Team/H1C_IDR3_Power_Spectra/blob/main/SPOILERS/All_Epochs_Power_Spectra/H1C_IDR3_Power_Spectra.ipynb) to take the per-baseline delay power spectra and form a one-dimensional, spherically-averaged (in Fourier space) power spectrum

The one-dimensional power spectrum is the ultimate output of HERA's power spectrum pipeline.

## python environment

To create the python environment required for the HERA power spectrum pipeline, we have included a `conda` environment yaml file for convenience.  The environment can be created via
```
conda env create -f environment.yaml
```
This environment can then be activated via
```
conda activate h1c-idr3
```

## Initial testing on NRAO

Before transferring data to UKSRC resources, tests were run on the NRAO system to verify the python environment.  The files associated with these tests are stored in the `nrao/` directory.  Please see `nrao/README.md` for more details on these tests.

## Follow up testing on Azimuth

An identical test to that run on NRAO was then run on Azimuth.  The files associated with these tests are stored in the `azimuth/` directory.  Please see `azimuth/README.md` for more details on these tests.

# UKSRC related links and information

**Confluence pages**

All relevant confluence pages are nested under [HERA Pipeline Overview](https://confluence.skatelescope.org/display/SRCSC/HERA+Pipeline+Overview).

**Jira links**

All relevant tickets are linked within the SKAEOR demonstrator case epic, [TEAL-617](https://jira.skatelescope.org/browse/TEAL-617).
