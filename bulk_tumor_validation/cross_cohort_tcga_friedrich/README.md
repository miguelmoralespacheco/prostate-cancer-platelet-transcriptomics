# Cross-cohort TCGA-PRAD and Friedrich (GSE134051)

This module consumes standardized, text-readable exports from the TCGA-PRAD
and Friedrich (GSE134051) cohort modules rather than cohort RDS objects.

Required cohort-owned contracts:

- `tcga_prad/Contracts/sample_scores.tsv`
- `tcga_prad/Contracts/hallmark_gsea.tsv`
- `friedrich_gse134051/Contracts/sample_scores.tsv`
- `friedrich_gse134051/Contracts/hallmark_gsea.tsv`

Schemas:

- `sample_scores.tsv`: `cohort_id`, `score_z`, `emt_ssgsea`
- `hallmark_gsea.tsv`: `cohort_id`, `analysis_id`, `pathway`, `nes`, `padj`

`cohort_summary.tsv` is not required and is not created. Sample identifiers are
validated upstream and omitted from the public contracts. The current contracts
materialize completed platelet-score, EMT and Q1/Q4 Hallmark results without
recomputing scientific models.

The two cross-cohort analytical consumers are:

- `Scripts/01_Q1Q4_Hallmark_dotplot.R`, which reads the two cohort-owned
  `hallmark_gsea.tsv` contracts and reproduces the selected Q1/Q4 Hallmark
  comparison.
- `Scripts/02_EMT_score_correlations.R`, which reads the two cohort-owned
  `sample_scores.tsv` contracts and reproduces the cohort-specific platelet
  score versus Hallmark EMT correlations.

Both scripts write only to ignored, regenerable locations:

- `Results/generated/`
- `Figures/generated/`

Run order:

1. Refresh the four cohort-owned contracts with the existing contract adapter
   only when completed source outputs change.
2. Run `Rscript Scripts/01_Q1Q4_Hallmark_dotplot.R`.
3. Run `Rscript Scripts/02_EMT_score_correlations.R`.

The consumers do not read cohort RDS objects, expression matrices, clinical
tables, or historical result directories, and they do not rerun differential
expression, enrichment, ssGSEA, or survival analyses.

The existing contract adapter is needed only when completed cohort outputs have
changed and the four tracked contracts must be refreshed.
