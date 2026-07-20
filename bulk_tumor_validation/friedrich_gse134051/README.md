# Friedrich (GSE134051) platelet score validation

- Cohort: Friedrich (GSE134051).
- Platform/source annotation: GPL26898 with normalized `gex.logq` expression.
- Final analytic cohort: 164 primary tumors.
- Expression data: normalized log-expression values.
- Signature availability: 38 of 41 genes.
- Unavailable genes: `FTL`, `CCL5`, and `NCOA4`.
- Score: arithmetic mean across the 38 available genes, followed by
  within-cohort z-standardization.
- Active analytical branches: rank-based Q1/Q4 and continuous `Score_z`.
- The Q1/Q4 branch uses the bottom and top `floor(n/4)` samples, yielding 41
  `LOW_Q1` and 41 `HIGH_Q4` primary tumors.
- Median splitting is used only for Kaplan-Meier visualization.
- Clinical survival analyses are exploratory because only 25 overall-survival
  events are available.

This module uses the repository-level canonical signature at
`../../resources/platelet_associated_transcriptional_signature.tsv`.
Portable cohort scripts are included. Raw and serialized inputs are not
distributed.

## Inputs and provenance

Local input requirements and provenance are documented in:

- [`Inputs/README.md`](Inputs/README.md);
- [`Inputs/DATA_SOURCES.tsv`](Inputs/DATA_SOURCES.tsv);
- [`Inputs/EXPECTED_INPUTS.tsv`](Inputs/EXPECTED_INPUTS.tsv).

The exact historical `curatedPCaData` package/resource version remains
unresolved and is labeled accordingly in `DATA_SOURCES.tsv`.

## Execution order

`Scripts/00_config_paths.R` is sourced by the analytical scripts and is not an
independent scientific analysis.

Initial score:

1. `Scripts/Score/01_compute_Friedrich_GSE134051_platelet_score_41genes.R`

Primary Q1/Q4 branch:

2. `Scripts/Q1Q4/01_limma_Q1Q4_41genes.R`
3. `Scripts/Q1Q4/02_gsea_Q1Q4.R`
4. `Scripts/Q1Q4/03_figures_Q1Q4_Friedrich_GSE134051.R`

Continuous branch:

2. `Scripts/Continuous/01_limma_continuous_41genes.R`
3. `Scripts/Continuous/02_gsea_continuous.R`
4. `Scripts/Continuous/03_figures_continuous_Friedrich_GSE134051.R`

Clinical analysis after the score:

- `Scripts/Clinical/01_Friedrich_GSE134051_clinical_final.R`

## Outputs

Generated outputs are written below `Results/generated/` and
`Figures/generated/` and are not tracked by Git.
