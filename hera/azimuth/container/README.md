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

## UKSRC related links and information

### Confluence pages

- [On The Fly (OTF) Observing Mode](https://confluence.skatelescope.org/x/4ce2E)

### Jira tickets

- [TEAL-711](https://jira.skatelescope.org/browse/TEAL-711): containerize HERA step 7 (power spectrum) pipeline
