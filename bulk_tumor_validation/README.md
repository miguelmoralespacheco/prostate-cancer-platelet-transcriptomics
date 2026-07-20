# Bulk-tumor validation

This block applies the repository-level canonical platelet-associated
transcriptional signature through three active modules:

- [`tcga_prad`](tcga_prad/README.md): TCGA-PRAD cohort analysis;
- [`friedrich_gse134051`](friedrich_gse134051/README.md): Friedrich (GSE134051)
  cohort analysis;
- [`cross_cohort_tcga_friedrich`](cross_cohort_tcga_friedrich/README.md):
  cross-cohort comparison of TCGA-PRAD and Friedrich (GSE134051).

The two cohort modules are independent. The CrossCohort module consumes their
tracked, cohort-owned text contracts and does not depend on full cohort RDS
objects.

All modules read the canonical signature from
`../resources/platelet_associated_transcriptional_signature.tsv` through a
repository-root-resolved path. Independent module-local signature copies are
not maintained.

## Inputs and provenance

Local large inputs and serialized analytical objects are not distributed. Each
cohort documents input provenance, evidence status, and expected local artifacts
in:

- `Inputs/README.md`;
- `Inputs/DATA_SOURCES.tsv`;
- `Inputs/EXPECTED_INPUTS.tsv`.

Unresolved historical acquisition details remain explicitly labeled as
unresolved or partially resolved in the source manifests; they are not
reconstructed or inferred.

## Output policy

Regenerable outputs are written below `Results/generated/` and
`Figures/generated/` and are ignored by Git.

## Top-level run sequence

1. Prepare validated local inputs according to the cohort input manifests.
2. Run the TCGA and Friedrich cohort modules using their documented sequences.
3. Refresh CrossCohort contracts only when completed cohort outputs change.
4. Run the two CrossCohort consumers.

See the three module READMEs above for their detailed script sequences.
