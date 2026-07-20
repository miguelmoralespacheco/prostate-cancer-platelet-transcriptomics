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

This module uses the repository-level canonical signature at
`../../resources/platelet_associated_transcriptional_signature.tsv`.
Portable cohort scripts are included. Raw and serialized inputs are not
distributed.

## Inputs and provenance

Local input requirements and provenance are documented in:

- [`Inputs/README.md`](Inputs/README.md);
- [`Inputs/DATA_SOURCES.tsv`](Inputs/DATA_SOURCES.tsv);
- [`Inputs/EXPECTED_INPUTS.tsv`](Inputs/EXPECTED_INPUTS.tsv).

Unresolved original GDC acquisition details remain labeled as unresolved or
partially resolved in `DATA_SOURCES.tsv`; no undocumented historical query or
acquisition workflow is asserted.

## Execution order

`Scripts/00_config_paths.R` is sourced by the analytical scripts and is not an
independent scientific analysis.

Initial score:

1. `Scripts/Score/01_compute_platelet_score_41genes.R`

Primary Q1/Q4 branch:

2. `Scripts/Q1Q4/01_deseq2_Q1Q4_41genes.R`
3. `Scripts/Q1Q4/02_gsea_Q1Q4.R`
4. `Scripts/Q1Q4/03_figures_Q1Q4_locked_style.R`

Continuous and sensitivity branch:

2. `Scripts/Continuous/01_deseq2_continuous_score.R`
3. `Scripts/Continuous/02_gsea_continuous_score.R`
4. `Scripts/Continuous/03_qc_tme_adjusted_residual_continuous_score.R`
5. `Scripts/Continuous/04_generate_continuous_EMT_GSEA_figure.R`

Additional analyses after the score and their required local or generated
upstream inputs are available:

- `Scripts/Heatmap_Gradient/01_emt_score_gradient_heatmap.R` uses the score
  master table, count matrix, and gene map.
- `Scripts/Clinical/TCGA_clinical_final.R` uses the score master table,
  expression container, and local clinical sources.

## Outputs

Generated outputs are written below `Results/generated/` and
`Figures/generated/` and are not tracked by Git.
