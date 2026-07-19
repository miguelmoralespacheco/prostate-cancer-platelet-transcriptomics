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
Rscript Scripts/05_export_public_signature_resource.R
```

Scripts 01-04 generate and annotate the internal signature-derivation outputs.
Script 05 validates those text-readable outputs and exports the canonical public
TSV and metadata JSON to the repository `resources/` directory.

`SCORE_CREATION_DIR` identifies the signature-construction project root and
defaults to the current working directory. `PUBLIC_RESOURCE_DIR` optionally
overrides the public output directory; when unset, script 05 uses `resources/`
at the parent of this module. Existing public files are protected unless
`PUBLIC_RESOURCE_OVERWRITE=true` is set explicitly.

## Expected outputs

The scripts generate:

- `Results_BoneMarrow/`
- `Results_Blood/`
- `Results_MergeSignature/`
- `Results_Reactome/`

Generated results are intentionally not tracked in Git.

`Seurat::AddModuleScore` is used only for diagnostic scoring of the single-cell
bone marrow and blood reference objects. It is not a universal bulk-tumor
scoring method; downstream cohort pipelines must document their expression
scale, gene aggregation, missing-gene handling, and within-cohort
standardization.
