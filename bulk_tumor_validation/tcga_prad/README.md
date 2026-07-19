# TCGA-PRAD platelet score validation

- Cohort: TCGA prostate adenocarcinoma.
- Final analytic cohort: 497 primary tumors.
- Signature availability: 41 of 41 genes.
- Expression scale: `log2(TPM + 1)`.
- Score: arithmetic mean across the 41 genes, followed by within-cohort
  z-standardization.
- Q1/Q4 is the primary extreme-group transcriptomic comparison.
- Continuous `Score_z` is the complementary robustness model.
- TME-adjusted and residualized analyses are TCGA-specific sensitivity
  analyses.
- Median splitting is used only for Kaplan-Meier visualization.

This module will use the repository-level canonical signature at
`../../resources/platelet_associated_transcriptional_signature.tsv`.
Data-acquisition instructions and portable scripts are pending later controlled
migration phases.
