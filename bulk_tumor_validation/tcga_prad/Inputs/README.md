# TCGA-PRAD inputs

Input data are not distributed in Git. The analytical scripts expect the seven
local files listed in `EXPECTED_INPUTS.tsv`; `DATA_SOURCES.tsv` records the
source and provenance status associated with each input family.

The listed checksums identify the current validated local artifacts. The exact
initial GDC query and release were not historically recovered, so users must
not assume an undocumented acquisition workflow. Clinical inputs and the
Hallmark snapshot require local placement and must remain untracked.

No participant-level data are documented in this repository.
