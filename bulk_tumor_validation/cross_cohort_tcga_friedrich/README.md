# Cross-cohort TCGA-PRAD and Friedrich (GSE134051)

This module will consume standardized, text-readable exports from the TCGA-PRAD
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
recomputing scientific models. Porting the two CrossCohort analytical scripts
to consume these contracts is a later patch.
