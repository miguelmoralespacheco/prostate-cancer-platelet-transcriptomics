# Prostate cancer platelet transcriptomics

This repository contains reproducible scripts for platelet-associated transcriptomic analyses in prostate cancer.

## Modules

- [`single_cell_platelet_signature/`](single_cell_platelet_signature/README.md): derivation of a platelet-associated transcriptional signature from Human Cell Atlas single-cell bone marrow and blood references.
- [`bulk_tumor_validation/tcga_prad/`](bulk_tumor_validation/tcga_prad/README.md): TCGA-PRAD bulk-tumor validation.
- [`bulk_tumor_validation/friedrich_gse134051/`](bulk_tumor_validation/friedrich_gse134051/README.md): Friedrich cohort validation, identified as Friedrich (GSE134051).
- [`bulk_tumor_validation/cross_cohort_tcga_friedrich/`](bulk_tumor_validation/cross_cohort_tcga_friedrich/README.md): TCGA-Friedrich cross-cohort analyses.

See the [bulk-tumor validation overview](bulk_tumor_validation/README.md) for module ownership and execution order.

## Data policy

Large raw data files are not distributed in this repository. The exact H5 inventory and checksums are documented in [`single_cell_platelet_signature/Inputs/DATA_SOURCES.tsv`](single_cell_platelet_signature/Inputs/DATA_SOURCES.tsv); obtain the matrices from the [Human Cell Atlas project](https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79) and place them in the expected `Inputs/*/Raw/` folders before running the scripts.

## Current status

The repository includes the implemented single-cell signature workflow and public bulk-tumor validation modules. Local cohort inputs are not distributed in Git. Generated results and figures are regenerable and ignored by Git.
