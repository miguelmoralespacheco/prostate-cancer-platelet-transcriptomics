# Friedrich (GSE134051) platelet score validation

- Cohort: Friedrich (GSE134051).
- Final analytic cohort: 164 primary tumors.
- Expression data: normalized log-expression values.
- Signature availability: 38 of 41 genes.
- Unavailable genes: `FTL`, `CCL5`, and `NCOA4`.
- Score: arithmetic mean across the 38 available genes, followed by
  within-cohort z-standardization.
- Planned analytical branches: Q1/Q4 and continuous `Score_z`.
- Median splitting is used only for Kaplan-Meier visualization.
- Clinical survival analyses are exploratory because only 25 overall-survival
  events are available.

This module will use the repository-level canonical signature at
`../../resources/platelet_associated_transcriptional_signature.tsv`.
Data acquisition and portable scripts are pending later controlled migration
phases.
