# Prostate cancer platelet transcriptomics

This repository contains reproducible scripts for platelet-associated transcriptomic analyses in prostate cancer.

## Modules

- `single_cell_platelet_signature/`: derivation of a platelet-associated transcriptional signature from Human Cell Atlas single-cell bone marrow and blood references.

## Data policy

Large raw data files are not distributed in this repository. Raw 10x Genomics H5 matrices must be downloaded from the original source and placed in the expected `Inputs/*/Raw/` folders before running the scripts.

## Current status

The current module provides the scripts required to reproduce the single-cell platelet-associated transcriptional signature derivation and Reactome pathway annotation.
