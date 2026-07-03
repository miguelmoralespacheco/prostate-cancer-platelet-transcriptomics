# Single-cell platelet-associated transcriptional signature

This module derives a platelet-associated transcriptional signature from Human Cell Atlas single-cell bone marrow and blood references.

## Required input structure

Raw 10x Genomics H5 files are not included in this repository. Place the required files according to the manifest files:

- `Inputs/BoneMarrow/Manifests/manifest_BM_HiSeq9.tsv`
- `Inputs/Blood/Manifests/manifest_BL_HiSeq9.tsv`

Expected raw data folders:

- `Inputs/BoneMarrow/Raw/`
- `Inputs/Blood/Raw/`

## Script order

Run from this module directory:

```bash
export SCORE_CREATION_DIR=$(pwd)

Rscript Scripts/01_build_bonemarrow_platelet_reference.R
Rscript Scripts/02_build_blood_platelet_reference.R
Rscript Scripts/03_define_platelet_associated_signature.R
Rscript Scripts/04_reactome_ORA_platelet_associated_signature.R
```

## Expected outputs

The scripts generate:

- `Results_BoneMarrow/`
- `Results_Blood/`
- `Results_MergeSignature/`
- `Results_Reactome/`

Generated results are intentionally not tracked in Git.
