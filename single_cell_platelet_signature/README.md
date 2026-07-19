# Single-cell platelet-associated transcriptional signature

This module derives a platelet-associated transcriptional signature from Human Cell Atlas single-cell bone marrow and blood references.

## Human Cell Atlas data source

The single-cell reference data were obtained from the Human Cell Atlas project
**A single cell immune cell atlas of human hematopoietic system**. Project UUID:
`cc95ff89-2e68-4a08-a234-480eca21ce79`. See the [HCA project page][hca-project-page].

The analysis used 16 raw 10x Genomics feature-barcode H5 matrices: eight adult
bone marrow samples and eight adult peripheral blood samples. The selected
matrices were generated using 10x Genomics 3' v2 chemistry and Illumina HiSeq X
sequencing. Bone marrow and peripheral blood references are processed
independently from manifest-defined inputs.

Raw H5 files are not distributed in this repository. Exact filenames, expected
local paths, byte sizes, and SHA-256 values are provided in
[`Inputs/DATA_SOURCES.tsv`](Inputs/DATA_SOURCES.tsv). Obtain the files from the
[HCA project page][hca-project-page] and place them in:

- `Inputs/BoneMarrow/Raw/`
- `Inputs/Blood/Raw/`

Before running scripts 01 and 02, verify candidate files against the published
checksums. Run this example from the module directory on macOS or Linux:

```bash
set -o pipefail

tail -n +2 Inputs/DATA_SOURCES.tsv |
while IFS=$'\t' read -r sample_id tissue local_file_name expected_relative_path \
  source_project source_project_uuid source_project_url library_preparation \
  sequencing_platform expected_sha256 size_bytes; do
  if [[ ! -f "$expected_relative_path" ]]; then
    printf 'MISSING\t%s\t%s\n' "$sample_id" "$expected_relative_path" >&2
    exit 1
  fi

  observed_sha256=$(shasum -a 256 "$expected_relative_path" | awk '{print $1}')
  if [[ "$observed_sha256" != "$expected_sha256" ]]; then
    printf 'MISMATCH\t%s\t%s\n' "$sample_id" "$expected_relative_path" >&2
    exit 1
  fi

  printf 'OK\t%s\t%s\n' "$sample_id" "$expected_relative_path"
done
```

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

[hca-project-page]: https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79
