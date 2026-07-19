# Cross-cohort TCGA-PRAD and Friedrich (GSE134051)

This module will consume standardized, text-readable exports from the TCGA-PRAD
and Friedrich (GSE134051) cohort modules.

The planned minimum contracts are:

- `sample_scores.tsv`;
- `hallmark_gsea.tsv`;
- `cohort_summary.tsv`.

These contracts are planned and are not yet implemented. The module must not
read full TCGA or Friedrich RDS objects, duplicated signature files, figures, or
manually copied intermediate outputs. It will use the repository-level
canonical signature at
`../../resources/platelet_associated_transcriptional_signature.tsv`, or the
equivalent repository-root-resolved path.
