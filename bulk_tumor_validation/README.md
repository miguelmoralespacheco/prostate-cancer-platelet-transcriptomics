# Bulk-tumor validation

This block applies the repository-level canonical platelet-associated
transcriptional signature to:

- TCGA-PRAD;
- Friedrich (GSE134051);
- cross-cohort comparison of TCGA-PRAD and Friedrich (GSE134051).

`tcga_prad` and `friedrich_gse134051` are independent cohort modules.
`cross_cohort_tcga_friedrich` will consume standardized, text-readable exports
from both modules and must not depend on full cohort RDS objects.

All modules must read the canonical signature from
`../../resources/platelet_associated_transcriptional_signature.tsv`, relative
to each cohort module, or from the equivalent repository-root-resolved path.
Independent module-local signature copies must not be maintained.

Portable cohort-specific analytical scripts are included for TCGA-PRAD and
Friedrich (GSE134051). Raw and serialized inputs are not distributed, and the
analytical workflows have not yet been executed from this public repository.

TCGA acquisition documentation and Friedrich acquisition/preparation
documentation remain pending. Cross-cohort scripts also remain pending until
standardized text-readable input contracts are defined.
