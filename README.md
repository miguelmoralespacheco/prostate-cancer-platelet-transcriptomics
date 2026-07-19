# Prostate cancer platelet transcriptomics

This repository contains reproducible scripts for platelet-associated transcriptomic analyses in prostate cancer.

## Modules

- `single_cell_platelet_signature/`: derivation of a platelet-associated transcriptional signature from Human Cell Atlas single-cell bone marrow and blood references.
- `bulk_tumor_validation/`: validation of the platelet-associated transcriptional score in TCGA-PRAD and Friedrich (GSE134051), including cross-cohort analyses.

## Data policy

Large raw data files are not distributed in this repository. The exact H5 inventory and checksums are documented in [`single_cell_platelet_signature/Inputs/DATA_SOURCES.tsv`](single_cell_platelet_signature/Inputs/DATA_SOURCES.tsv); obtain the matrices from the [Human Cell Atlas project](https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79) and place them in the expected `Inputs/*/Raw/` folders before running the scripts.

## Current status

The current module provides the scripts required to reproduce the single-cell platelet-associated transcriptional signature derivation and Reactome pathway annotation.
