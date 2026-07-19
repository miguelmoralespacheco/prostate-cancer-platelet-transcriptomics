#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD | Official 41-gene platelet-associated transcriptional score
#
# Initial cohort:
#   - 501 primary tumor samples.
#   - 497 unique patients.
#
# Duplicate-patient resolution:
#   - For patients with more than one primary tumor sample, retain the sample
#     with the largest library size from the "unstranded" count assay.
#   - Use the complete barcode as a deterministic tie-breaker.
#
# Score:
#   Score_raw = mean(log2(TPM + 1)) across available signature genes.
#   Score_z   = scale(Score_raw) across the final 497-patient cohort.
#
# This script does not run differential expression, GSEA, clinical models,
# or plotting.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 0. Frozen analytical constants
################################################################################

EXPECTED_INITIAL_SAMPLES <- 501L
EXPECTED_UNIQUE_PATIENTS <- 497L
EXPECTED_DUPLICATE_PATIENTS <- 4L
EXPECTED_EXCLUDED_SAMPLES <- 4L
EXPECTED_FINAL_SAMPLES <- 497L

EXPECTED_SIGNATURE_GENES <- 41L
MIN_SIGNATURE_GENES_PRESENT <- 35L

PRIMARY_SAMPLE_TYPE <- "Primary Tumor"
TPM_ASSAY <- "tpm_unstrand"
COUNT_ASSAY <- "unstranded"

################################################################################
# 1. Project paths
################################################################################

.repo_override <- Sys.getenv("PLATELET_REPO_ROOT", unset = "")
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (nzchar(.repo_override)) {
  .config_path <- file.path(
    normalizePath(.repo_override, winslash = "/", mustWork = TRUE),
    "bulk_tumor_validation", "tcga_prad", "Scripts", "00_config_paths.R"
  )
} else {
  if (length(.script_args) != 1L) {
    stop("Set PLATELET_REPO_ROOT when sourcing this script.", call. = FALSE)
  }
  .script_file <- normalizePath(
    sub("^--file=", "", .script_args[[1]]), winslash = "/", mustWork = TRUE
  )
  .config_path <- file.path(dirname(dirname(.script_file)), "00_config_paths.R")
}
source(.config_path, local = FALSE)
project_dir <- MODULE_DIR
rm(.repo_override, .script_args, .config_path)
if (exists(".script_file")) rm(.script_file)

script_relative_path <- file.path(
  "Scripts", "Score", "01_compute_platelet_score_41genes.R"
)

input_paths <- list(
  signature = file.path(
    RESOURCE_DIR,
    "platelet_associated_transcriptional_signature.tsv"
  ),
  se_tumor = file.path(
    INPUT_DIR, "LocalLarge",
    "TCGA_PRAD_SE_tumorOnly.rds"
  ),
  metadata = file.path(
    INPUT_DIR, "Metadata",
    "TCGA_PRAD_colData_tumorOnly_clean.csv"
  ),
  gene_map = file.path(
    INPUT_DIR, "Metadata",
    "MAP_ENSGversion_to_geneName.csv"
  )
)

output_dir <- file.path(GENERATED_RESULTS_DIR, "Score")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_paths <- list(
  master_csv = file.path(
    output_dir, "TCGA_primary_master_table_41genes.csv"
  ),
  master_rds = file.path(
    output_dir, "TCGA_primary_master_table_41genes.rds"
  ),
  duplicate_resolution = file.path(
    output_dir, "TCGA_duplicate_patient_resolution.csv"
  ),
  signature_availability = file.path(
    output_dir, "TCGA_signature_gene_availability.csv"
  ),
  metadata_json = file.path(
    output_dir, "TCGA_score_metadata.json"
  ),
  qc = file.path(
    output_dir, "TCGA_score_QC.txt"
  )
)

################################################################################
# 2. Helpers
################################################################################

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

warnings_vec <- character(0)

add_warning <- function(x) {
  warnings_vec <<- unique(c(warnings_vec, x))
  warning(x, call. = FALSE)
  invisible(x)
}

fail_with_qc <- function(msg_text) {
  failure_lines <- c(
    "TCGA-PRAD official 41-gene platelet score QC",
    "================================================",
    "",
    paste0("Date/time: ", as.character(Sys.time())),
    paste0("Project dir: ", project_dir),
    paste0("Failure reason: ", msg_text),
    "",
    "Final status: FAIL"
  )

  writeLines(failure_lines, output_paths$qc, useBytes = TRUE)
  stop(msg_text, call. = FALSE)
}

fail <- function(...) {
  fail_with_qc(paste0(...))
}

require_input <- function(path, label) {
  if (!file.exists(path)) {
    fail("Missing required input: ", label, " at ", path)
  }
}

read_required_csv <- function(path, label) {
  require_input(path, label)

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_csv <- function(x, path) {
  utils::write.csv(
    as.data.frame(x),
    path,
    row.names = FALSE,
    quote = TRUE
  )
}

detect_gene_column <- function(df) {
  normalized_names <- tolower(gsub("[^a-z0-9]+", "", names(df)))

  candidates <- c(
    "gene",
    "gene_symbol",
    "genesymbol",
    "symbol",
    "hgncsymbol",
    "hugo",
    "genename",
    "gene_name"
  )

  normalized_candidates <- tolower(
    gsub("[^a-z0-9]+", "", candidates)
  )

  hit <- which(normalized_names %in% normalized_candidates)

  if (length(hit) == 0) {
    fail(
      "Could not detect gene symbol column in signature. Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }

  names(df)[hit[1]]
}

################################################################################
# 3. Required packages
################################################################################

required_packages <- c(
  "SummarizedExperiment",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  fail(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". The script does not install packages."
  )
}

################################################################################
# 4. Validate input files
################################################################################

msg("Project dir:", project_dir)

for (nm in names(input_paths)) {
  require_input(input_paths[[nm]], nm)
}

################################################################################
# 5. Load and validate the official signature
################################################################################

msg("[1/8] Loading official 41-gene signature")

signature_df <- read_canonical_platelet_signature()

signature_gene_col <- detect_gene_column(signature_df)

signature_genes <- trimws(
  as.character(signature_df[[signature_gene_col]])
)

signature_genes <- signature_genes[
  !is.na(signature_genes) & signature_genes != ""
]

signature_keys <- toupper(signature_genes)

if (anyDuplicated(signature_keys)) {
  duplicated_signature_genes <- unique(
    signature_genes[duplicated(signature_keys)]
  )

  fail(
    "Official signature contains duplicated genes after case-insensitive ",
    "normalization: ",
    paste(duplicated_signature_genes, collapse = ", ")
  )
}

if (length(signature_keys) != EXPECTED_SIGNATURE_GENES) {
  fail(
    "Official signature must contain exactly ",
    EXPECTED_SIGNATURE_GENES,
    " unique genes. Observed: ",
    length(signature_keys)
  )
}

msg("Signature genes:", length(signature_keys))

################################################################################
# 6. Load and validate the tumor-only SummarizedExperiment
################################################################################

msg("[2/8] Loading tumor-only SummarizedExperiment")

se_tumor <- readRDS(input_paths$se_tumor)

if (!inherits(se_tumor, "SummarizedExperiment")) {
  fail(
    "TCGA_PRAD_SE_tumorOnly.rds must be a SummarizedExperiment or ",
    "RangedSummarizedExperiment. Class: ",
    paste(class(se_tumor), collapse = ", ")
  )
}

assays_available <- SummarizedExperiment::assayNames(se_tumor)

msg(
  "Assays available:",
  paste(assays_available, collapse = ", ")
)

for (required_assay in c(TPM_ASSAY, COUNT_ASSAY)) {
  if (!required_assay %in% assays_available) {
    fail(
      "Required assay '", required_assay,
      "' was not found. Available assays: ",
      paste(assays_available, collapse = ", ")
    )
  }
}

tpm_all <- SummarizedExperiment::assay(se_tumor, TPM_ASSAY)
counts_all <- SummarizedExperiment::assay(se_tumor, COUNT_ASSAY)

if (is.null(dim(tpm_all)) || length(dim(tpm_all)) != 2) {
  fail("The TPM assay is not a two-dimensional matrix-like object.")
}

if (is.null(dim(counts_all)) || length(dim(counts_all)) != 2) {
  fail("The unstranded count assay is not a two-dimensional matrix-like object.")
}

if (!identical(dim(tpm_all), dim(counts_all))) {
  fail("TPM and unstranded count assays do not have identical dimensions.")
}

if (!identical(colnames(tpm_all), colnames(counts_all))) {
  fail("TPM and unstranded count assays do not have identical sample order.")
}

sample_ids_all <- colnames(tpm_all)

if (
  is.null(sample_ids_all) ||
  any(is.na(sample_ids_all)) ||
  any(sample_ids_all == "")
) {
  fail("Expression assays must contain non-empty sample barcodes.")
}

if (anyDuplicated(sample_ids_all)) {
  fail("Expression assays contain duplicated sample barcodes.")
}

if (length(sample_ids_all) != EXPECTED_INITIAL_SAMPLES) {
  fail(
    "Expected ", EXPECTED_INITIAL_SAMPLES,
    " initial primary tumor samples. Observed: ",
    length(sample_ids_all)
  )
}

################################################################################
# 7. Load and align tumor metadata
################################################################################

msg("[3/8] Loading and aligning tumor metadata")

metadata_df <- read_required_csv(
  input_paths$metadata,
  "clean tumor metadata"
)

required_metadata_columns <- c(
  "barcode",
  "sample_type"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(metadata_df)
)

if (length(missing_metadata_columns) > 0) {
  fail(
    "Metadata is missing required columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

metadata_df$barcode <- trimws(as.character(metadata_df$barcode))
metadata_df$sample_type <- trimws(as.character(metadata_df$sample_type))

if (
  any(is.na(metadata_df$barcode)) ||
  any(metadata_df$barcode == "")
) {
  fail("Metadata contains missing or empty barcode values.")
}

if (anyDuplicated(metadata_df$barcode)) {
  fail("Metadata contains duplicated barcode values.")
}

if (!setequal(sample_ids_all, metadata_df$barcode)) {
  missing_in_metadata <- setdiff(sample_ids_all, metadata_df$barcode)
  missing_in_se <- setdiff(metadata_df$barcode, sample_ids_all)

  fail(
    "Cannot align metadata to expression samples. Missing in metadata: ",
    length(missing_in_metadata),
    "; missing in SummarizedExperiment: ",
    length(missing_in_se),
    "."
  )
}

metadata_aligned_all <- metadata_df[
  match(sample_ids_all, metadata_df$barcode),
  ,
  drop = FALSE
]

if (!identical(metadata_aligned_all$barcode, sample_ids_all)) {
  fail("Metadata alignment failed after matching by barcode.")
}

sample_type_counts <- table(
  metadata_aligned_all$sample_type,
  useNA = "ifany"
)

if (
  length(sample_type_counts) != 1 ||
  names(sample_type_counts)[1] != PRIMARY_SAMPLE_TYPE ||
  as.integer(sample_type_counts[1]) != EXPECTED_INITIAL_SAMPLES
) {
  fail(
    "The frozen input cohort must contain exactly ",
    EXPECTED_INITIAL_SAMPLES,
    " samples labeled '", PRIMARY_SAMPLE_TYPE,
    "'. Observed: ",
    paste(
      names(sample_type_counts),
      as.integer(sample_type_counts),
      sep = "=",
      collapse = "; "
    )
  )
}

################################################################################
# 8. Resolve duplicated patients deterministically
################################################################################

msg("[4/8] Resolving duplicated primary samples per patient")

if ("patient" %in% names(metadata_aligned_all)) {
  patient_id_all <- trimws(
    as.character(metadata_aligned_all$patient)
  )
} else {
  patient_id_all <- substr(sample_ids_all, 1, 12)

  add_warning(
    "Metadata does not contain a patient column; patient IDs were derived ",
    "from the first 12 characters of each TCGA barcode."
  )
}

if (
  any(is.na(patient_id_all)) ||
  any(patient_id_all == "")
) {
  fail("Patient identifiers contain missing or empty values.")
}

n_unique_patients <- length(unique(patient_id_all))

if (n_unique_patients != EXPECTED_UNIQUE_PATIENTS) {
  fail(
    "Expected ", EXPECTED_UNIQUE_PATIENTS,
    " unique patients. Observed: ",
    n_unique_patients
  )
}

patient_sample_counts <- table(patient_id_all)

duplicate_patient_ids <- names(
  patient_sample_counts[patient_sample_counts > 1]
)

n_duplicate_patients <- length(duplicate_patient_ids)

if (n_duplicate_patients != EXPECTED_DUPLICATE_PATIENTS) {
  fail(
    "Expected ", EXPECTED_DUPLICATE_PATIENTS,
    " patients with multiple primary tumor samples. Observed: ",
    n_duplicate_patients
  )
}

library_size_all <- colSums(counts_all, na.rm = TRUE)
names(library_size_all) <- colnames(counts_all)

if (any(!is.finite(library_size_all))) {
  fail("Non-finite unstranded library sizes were detected.")
}

sample_selection <- data.frame(
  patient_id = patient_id_all,
  sample_id = sample_ids_all,
  barcode = sample_ids_all,
  sample_type = metadata_aligned_all$sample_type,
  unstranded_library_size = as.numeric(
    library_size_all[sample_ids_all]
  ),
  patient_sample_count = as.integer(
    patient_sample_counts[patient_id_all]
  ),
  original_assay_order = seq_along(sample_ids_all),
  stringsAsFactors = FALSE
)

sample_selection$max_library_size_within_patient <- ave(
  sample_selection$unstranded_library_size,
  sample_selection$patient_id,
  FUN = max
)

sample_selection$n_samples_at_max_library_size <- ave(
  sample_selection$unstranded_library_size,
  sample_selection$patient_id,
  FUN = function(x) sum(x == max(x))
)

sample_selection <- sample_selection[
  order(
    sample_selection$patient_id,
    -sample_selection$unstranded_library_size,
    sample_selection$barcode
  ),
  ,
  drop = FALSE
]

sample_selection$selection_rank_within_patient <- ave(
  seq_len(nrow(sample_selection)),
  sample_selection$patient_id,
  FUN = seq_along
)

sample_selection$retained <- (
  sample_selection$selection_rank_within_patient == 1L
)

sample_selection$selection_reason <- ifelse(
  sample_selection$patient_sample_count == 1L,
  "only_primary_sample",
  ifelse(
    sample_selection$retained &
      sample_selection$n_samples_at_max_library_size == 1L,
    "retained_max_unstranded_library_size",
    ifelse(
      sample_selection$retained &
        sample_selection$n_samples_at_max_library_size > 1L,
      "retained_barcode_tiebreak",
      ifelse(
        sample_selection$unstranded_library_size <
          sample_selection$max_library_size_within_patient,
        "excluded_lower_unstranded_library_size",
        "excluded_barcode_tiebreak"
      )
    )
  )
)

retained_sample_ids <- sample_selection$sample_id[
  sample_selection$retained
]

excluded_sample_ids <- sample_selection$sample_id[
  !sample_selection$retained
]

if (length(excluded_sample_ids) != EXPECTED_EXCLUDED_SAMPLES) {
  fail(
    "Expected ", EXPECTED_EXCLUDED_SAMPLES,
    " excluded duplicate samples. Observed: ",
    length(excluded_sample_ids)
  )
}

if (length(retained_sample_ids) != EXPECTED_FINAL_SAMPLES) {
  fail(
    "Expected ", EXPECTED_FINAL_SAMPLES,
    " retained unique-patient samples. Observed: ",
    length(retained_sample_ids)
  )
}

keep_idx <- which(sample_ids_all %in% retained_sample_ids)

if (length(keep_idx) != EXPECTED_FINAL_SAMPLES) {
  fail(
    "Final sample selection could not be matched back to the ",
    "expression assays."
  )
}

# Preserve the original assay order after duplicate resolution.
tpm <- tpm_all[, keep_idx, drop = FALSE]

metadata_aligned <- metadata_aligned_all[
  keep_idx,
  ,
  drop = FALSE
]

sample_ids <- colnames(tpm)
patient_id <- patient_id_all[keep_idx]
library_size <- library_size_all[sample_ids]

if (!identical(metadata_aligned$barcode, sample_ids)) {
  fail(
    "Metadata and expression assay became misaligned after ",
    "duplicate resolution."
  )
}

if (anyDuplicated(patient_id)) {
  fail("Duplicated patients remain after deterministic sample selection.")
}

if (length(sample_ids) != EXPECTED_FINAL_SAMPLES) {
  fail(
    "Final analytical cohort must contain exactly ",
    EXPECTED_FINAL_SAMPLES,
    " samples. Observed: ",
    length(sample_ids)
  )
}

write_csv(
  sample_selection,
  output_paths$duplicate_resolution
)

msg("Initial primary samples:", length(sample_ids_all))
msg("Unique patients:", n_unique_patients)
msg("Duplicate patients resolved:", n_duplicate_patients)
msg("Excluded duplicate samples:", length(excluded_sample_ids))
msg("Final analytical samples:", length(sample_ids))

################################################################################
# 9. Load and validate the Ensembl-to-symbol map
################################################################################

msg("[5/8] Mapping Ensembl rows to official gene symbols")

map_df <- read_required_csv(
  input_paths$gene_map,
  "Ensembl-version to gene-symbol map"
)

required_map_columns <- c(
  "ensg_version",
  "gene_name"
)

missing_map_columns <- setdiff(
  required_map_columns,
  names(map_df)
)

if (length(missing_map_columns) > 0) {
  fail(
    "Gene map is missing required columns: ",
    paste(missing_map_columns, collapse = ", ")
  )
}

map_df$ensg_version <- trimws(as.character(map_df$ensg_version))
map_df$gene_name <- trimws(as.character(map_df$gene_name))

map_df <- map_df[
  !is.na(map_df$ensg_version) &
    map_df$ensg_version != "",
  ,
  drop = FALSE
]

if (anyDuplicated(map_df$ensg_version)) {
  duplicated_ensg <- unique(
    map_df$ensg_version[duplicated(map_df$ensg_version)]
  )

  conflicting_ensg <- duplicated_ensg[
    vapply(
      duplicated_ensg,
      function(id) {
        symbols <- unique(
          map_df$gene_name[
            map_df$ensg_version == id &
              !is.na(map_df$gene_name) &
              map_df$gene_name != ""
          ]
        )

        length(symbols) > 1
      },
      logical(1)
    )
  ]

  if (length(conflicting_ensg) > 0) {
    fail(
      "Gene map contains Ensembl-version identifiers mapped to multiple ",
      "gene symbols. Examples: ",
      paste(utils::head(conflicting_ensg, 20), collapse = ", ")
    )
  }

  map_df <- map_df[
    !duplicated(map_df$ensg_version),
    ,
    drop = FALSE
  ]
}

tpm_row_ids <- rownames(tpm)

if (
  is.null(tpm_row_ids) ||
  any(is.na(tpm_row_ids)) ||
  any(tpm_row_ids == "")
) {
  fail("TPM assay must contain non-empty Ensembl row identifiers.")
}

map_aligned <- map_df[
  match(tpm_row_ids, map_df$ensg_version),
  ,
  drop = FALSE
]

mapped_symbols <- trimws(as.character(map_aligned$gene_name))
mapped_symbol_keys <- toupper(mapped_symbols)

signature_row_mask <- (
  !is.na(mapped_symbol_keys) &
    mapped_symbol_keys != "" &
    mapped_symbol_keys %in% signature_keys
)

if (sum(signature_row_mask) == 0) {
  fail(
    "No expression rows could be mapped to genes in the official signature."
  )
}

signature_row_keys <- mapped_symbol_keys[signature_row_mask]

signature_expression_rows <- as.matrix(
  tpm[signature_row_mask, , drop = FALSE]
)

if (any(!is.finite(signature_expression_rows))) {
  fail(
    "Non-finite TPM values were detected in signature-mapped ",
    "expression rows."
  )
}

if (any(signature_expression_rows < 0)) {
  fail(
    "Negative TPM values were detected in signature-mapped ",
    "expression rows."
  )
}

################################################################################
# 10. Collapse duplicate gene symbols and determine signature availability
################################################################################

signature_row_counts <- table(signature_row_keys)

signature_tpm_sum <- rowsum(
  signature_expression_rows,
  group = signature_row_keys,
  reorder = FALSE
)

rows_per_symbol <- as.numeric(
  signature_row_counts[rownames(signature_tpm_sum)]
)

signature_tpm_by_symbol <- sweep(
  signature_tpm_sum,
  1,
  rows_per_symbol,
  "/"
)

n_rows_per_signature_gene <- as.integer(
  signature_row_counts[
    match(signature_keys, names(signature_row_counts))
  ]
)

n_rows_per_signature_gene[
  is.na(n_rows_per_signature_gene)
] <- 0L

signature_availability <- data.frame(
  signature_order = seq_along(signature_genes),
  gene_symbol = signature_genes,
  gene_key = signature_keys,
  present_in_TCGA = n_rows_per_signature_gene > 0L,
  n_ensembl_rows_mapped = n_rows_per_signature_gene,
  stringsAsFactors = FALSE
)

present_mask <- signature_availability$present_in_TCGA

present_genes <- signature_genes[present_mask]
present_keys <- signature_keys[present_mask]
missing_genes <- signature_genes[!present_mask]

n_present <- length(present_genes)
n_missing <- length(missing_genes)

if (n_present < MIN_SIGNATURE_GENES_PRESENT) {
  fail(
    "Fewer than ", MIN_SIGNATURE_GENES_PRESENT,
    " official signature genes are available in TCGA: ",
    n_present, "/", EXPECTED_SIGNATURE_GENES
  )
}

if (n_missing > 0) {
  add_warning(
    paste0(
      "The following official signature genes were unavailable in TCGA: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

n_duplicate_signature_symbols <- sum(
  signature_availability$n_ensembl_rows_mapped > 1L
)

write_csv(
  signature_availability,
  output_paths$signature_availability
)

msg("Signature genes present:", n_present, "/", EXPECTED_SIGNATURE_GENES)

if (n_missing > 0) {
  msg("Missing signature genes:", paste(missing_genes, collapse = ", "))
}

msg(
  "Signature symbols represented by multiple Ensembl rows:",
  n_duplicate_signature_symbols
)

################################################################################
# 11. Calculate the official platelet-associated transcriptional score
################################################################################

msg("[6/8] Calculating the official platelet-associated score")

tpm_signature <- signature_tpm_by_symbol[
  present_keys,
  ,
  drop = FALSE
]

if (!identical(colnames(tpm_signature), sample_ids)) {
  fail("Signature TPM matrix is not aligned to the final sample cohort.")
}

log2_tpm_signature <- log2(tpm_signature + 1)

score_raw <- colMeans(
  log2_tpm_signature,
  na.rm = FALSE
)

if (any(!is.finite(score_raw))) {
  fail("Score_raw contains non-finite values.")
}

if (!is.finite(stats::sd(score_raw)) || stats::sd(score_raw) == 0) {
  fail("Score_raw has zero or non-finite variance.")
}

score_z <- as.numeric(scale(score_raw))
names(score_z) <- names(score_raw)

if (any(!is.finite(score_z))) {
  fail("Score_z contains non-finite values.")
}

score_z_mean <- mean(score_z)
score_z_sd <- stats::sd(score_z)

if (abs(score_z_mean) > 1e-10) {
  add_warning(
    paste0(
      "Score_z mean differs from zero beyond numerical tolerance: ",
      format(score_z_mean, digits = 8)
    )
  )
}

if (abs(score_z_sd - 1) > 1e-10) {
  add_warning(
    paste0(
      "Score_z standard deviation differs from one beyond numerical ",
      "tolerance: ",
      format(score_z_sd, digits = 8)
    )
  )
}

quartile_cutoffs <- stats::quantile(
  score_raw,
  probs = c(0.25, 0.75),
  na.rm = TRUE,
  names = FALSE,
  type = 7
)

names(quartile_cutoffs) <- c(
  "Q1_25pct",
  "Q3_75pct"
)

platelet_score_group <- ifelse(
  score_raw <= quartile_cutoffs["Q1_25pct"],
  "LOW_Q1",
  ifelse(
    score_raw >= quartile_cutoffs["Q3_75pct"],
    "HIGH_Q4",
    "MID"
  )
)

score_quartile <- ifelse(
  platelet_score_group == "LOW_Q1",
  "Q1",
  ifelse(
    platelet_score_group == "HIGH_Q4",
    "Q4",
    "Q2_Q3"
  )
)

group_counts <- table(
  factor(
    platelet_score_group,
    levels = c("LOW_Q1", "MID", "HIGH_Q4")
  )
)

if (
  as.integer(group_counts["LOW_Q1"]) != 125L ||
  as.integer(group_counts["HIGH_Q4"]) != 125L
) {
  add_warning(
    paste0(
      "Q1/Q4 group sizes differ from the expected 125/125, likely because ",
      "of tied score values. Observed LOW_Q1=",
      as.integer(group_counts["LOW_Q1"]),
      ", HIGH_Q4=",
      as.integer(group_counts["HIGH_Q4"]),
      "."
    )
  )
}

msg(
  "Score groups:",
  paste(
    names(group_counts),
    as.integer(group_counts),
    sep = "=",
    collapse = "; "
  )
)

################################################################################
# 12. Build the final master table
################################################################################

msg("[7/8] Building the final score master table")

score_table <- data.frame(
  sample_id = sample_ids,
  barcode = sample_ids,
  patient_id = patient_id,
  sample_type = metadata_aligned$sample_type,
  unstranded_library_size = as.numeric(library_size),
  Score_raw = as.numeric(score_raw[sample_ids]),
  Score_z = as.numeric(score_z[sample_ids]),
  score_quartile = as.character(score_quartile[sample_ids]),
  platelet_score_group = as.character(
    platelet_score_group[sample_ids]
  ),
  n_signature_genes_available = n_present,
  n_signature_genes_missing = n_missing,
  stringsAsFactors = FALSE
)

if ("case_id" %in% names(metadata_aligned)) {
  score_table$case_id <- as.character(metadata_aligned$case_id)
}

if (anyDuplicated(score_table$sample_id)) {
  fail("Final score table contains duplicated sample IDs.")
}

if (anyDuplicated(score_table$patient_id)) {
  fail("Final score table contains duplicated patient IDs.")
}

if (nrow(score_table) != EXPECTED_FINAL_SAMPLES) {
  fail(
    "Final score table must contain exactly ",
    EXPECTED_FINAL_SAMPLES,
    " rows. Observed: ",
    nrow(score_table)
  )
}

if (
  any(!is.finite(score_table$Score_raw)) ||
  any(!is.finite(score_table$Score_z))
) {
  fail("Final score table contains non-finite score values.")
}

write_csv(score_table, output_paths$master_csv)
saveRDS(score_table, output_paths$master_rds)

################################################################################
# 13. Metadata
################################################################################

duplicate_resolution_subset <- sample_selection[
  sample_selection$patient_sample_count > 1,
  c(
    "patient_id",
    "sample_id",
    "unstranded_library_size",
    "selection_rank_within_patient",
    "retained",
    "selection_reason"
  ),
  drop = FALSE
]

retained_duplicate_samples <- duplicate_resolution_subset$sample_id[
  duplicate_resolution_subset$retained
]

excluded_duplicate_samples <- duplicate_resolution_subset$sample_id[
  !duplicate_resolution_subset$retained
]

group_counts_list <- as.list(
  stats::setNames(
    as.integer(group_counts),
    names(group_counts)
  )
)

metadata_json <- list(
  date_time = as.character(Sys.time()),
  project_dir = project_dir,
  script = script_relative_path,
  project = "tcga_prad",
  cohort = "TCGA-PRAD primary tumors, one sample per patient",
  analysis = "official 41-gene platelet-associated transcriptional score",
  score_formula = list(
    Score_raw = paste0(
      "mean(log2(TPM + 1)) across available official ",
      "signature genes"
    ),
    Score_z = paste0(
      "standardized Score_raw across the final ",
      EXPECTED_FINAL_SAMPLES,
      "-patient cohort"
    )
  ),
  assays = list(
    score_assay = TPM_ASSAY,
    duplicate_resolution_assay = COUNT_ASSAY,
    assays_available = as.character(assays_available)
  ),
  inputs = input_paths,
  cohort_selection = list(
    initial_primary_samples = length(sample_ids_all),
    unique_patients = n_unique_patients,
    duplicate_patients = n_duplicate_patients,
    duplicate_patient_ids = duplicate_patient_ids,
    selection_rule = paste0(
      "For patients with multiple primary tumor samples, retain the ",
      "sample with the largest unstranded library size; use barcode ",
      "as a deterministic tie-breaker."
    ),
    retained_duplicate_samples = retained_duplicate_samples,
    excluded_duplicate_samples = excluded_duplicate_samples,
    excluded_samples = length(excluded_sample_ids),
    final_samples = length(sample_ids)
  ),
  signature = list(
    input_genes = length(signature_keys),
    genes_present = n_present,
    genes_missing = n_missing,
    present_genes = present_genes,
    missing_genes = missing_genes,
    duplicate_symbol_strategy = paste0(
      "Mean TPM across Ensembl-version rows mapping to the same ",
      "case-insensitive gene symbol."
    ),
    symbols_with_multiple_ensembl_rows = n_duplicate_signature_symbols
  ),
  score_distribution = list(
    Score_raw_min = min(score_table$Score_raw),
    Score_raw_median = stats::median(score_table$Score_raw),
    Score_raw_max = max(score_table$Score_raw),
    Score_z_mean = score_z_mean,
    Score_z_sd = score_z_sd
  ),
  quartiles = list(
    Q1_25pct = unname(quartile_cutoffs["Q1_25pct"]),
    Q3_75pct = unname(quartile_cutoffs["Q3_75pct"]),
    group_counts = group_counts_list
  ),
  outputs = output_paths,
  reproducibility = list(
    R_version = R.version.string,
    SummarizedExperiment_version = as.character(
      utils::packageVersion("SummarizedExperiment")
    ),
    jsonlite_version = as.character(
      utils::packageVersion("jsonlite")
    )
  ),
  warnings = warnings_vec,
  qc_status = if (length(warnings_vec) == 0) {
    "PASS"
  } else {
    "PASS_WITH_WARNINGS"
  }
)

jsonlite::write_json(
  metadata_json,
  output_paths$metadata_json,
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null",
  digits = NA
)

################################################################################
# 14. Final QC report
################################################################################

msg("[8/8] Writing final QC report")

final_status <- if (length(warnings_vec) == 0) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

duplicate_qc_lines <- capture.output(
  print(
    duplicate_resolution_subset,
    row.names = FALSE
  )
)

qc_lines <- c(
  "TCGA-PRAD official 41-gene platelet score QC",
  "================================================",
  "",
  paste0("Date/time: ", as.character(Sys.time())),
  paste0("Project dir: ", project_dir),
  paste0("Script: ", script_relative_path),
  "",
  "1. Inputs",
  paste0("  Signature: ", input_paths$signature),
  paste0("  Tumor-only SummarizedExperiment: ", input_paths$se_tumor),
  paste0("  Tumor metadata: ", input_paths$metadata),
  paste0("  Gene map: ", input_paths$gene_map),
  "",
  "2. Expression assays",
  paste0("  Available assays: ", paste(assays_available, collapse = ", ")),
  paste0("  TPM assay used for score: ", TPM_ASSAY),
  paste0("  Count assay used for duplicate resolution: ", COUNT_ASSAY),
  "",
  "3. Cohort",
  paste0("  Initial primary tumor samples: ", length(sample_ids_all)),
  paste0("  Unique patients: ", n_unique_patients),
  paste0("  Patients with duplicated primary samples: ", n_duplicate_patients),
  paste0("  Excluded duplicate samples: ", length(excluded_sample_ids)),
  paste0("  Final unique-patient samples: ", length(sample_ids)),
  "",
  "4. Duplicate-patient resolution",
  paste0(
    "  Rule: retain the sample with the largest unstranded library size; ",
    "use barcode as deterministic tie-breaker."
  ),
  duplicate_qc_lines,
  "",
  "5. Official signature",
  paste0("  Signature genes in input: ", length(signature_keys)),
  paste0("  Signature genes present: ", n_present, "/", EXPECTED_SIGNATURE_GENES),
  paste0(
    "  Missing genes: ",
    if (n_missing == 0) "none" else paste(missing_genes, collapse = ", ")
  ),
  paste0(
    "  Signature symbols represented by multiple Ensembl rows: ",
    n_duplicate_signature_symbols
  ),
  paste0(
    "  Duplicate-symbol handling: mean TPM across Ensembl rows mapping ",
    "to the same case-insensitive gene symbol."
  ),
  "",
  "6. Score",
  "  Score_raw = mean(log2(TPM + 1)) across available signature genes",
  paste0(
    "  Score_z = scale(Score_raw) across the final ",
    EXPECTED_FINAL_SAMPLES,
    "-patient cohort"
  ),
  paste0("  Score_raw minimum: ", format(min(score_table$Score_raw), digits = 10)),
  paste0(
    "  Score_raw median: ",
    format(stats::median(score_table$Score_raw), digits = 10)
  ),
  paste0("  Score_raw maximum: ", format(max(score_table$Score_raw), digits = 10)),
  paste0("  Score_z mean: ", format(score_z_mean, digits = 10)),
  paste0("  Score_z SD: ", format(score_z_sd, digits = 10)),
  "",
  "7. Quartiles",
  paste0(
    "  Q1 25% cutoff: ",
    format(quartile_cutoffs["Q1_25pct"], digits = 10)
  ),
  paste0(
    "  Q3 75% cutoff: ",
    format(quartile_cutoffs["Q3_75pct"], digits = 10)
  ),
  paste0(
    "  LOW_Q1: ",
    as.integer(group_counts["LOW_Q1"])
  ),
  paste0(
    "  MID: ",
    as.integer(group_counts["MID"])
  ),
  paste0(
    "  HIGH_Q4: ",
    as.integer(group_counts["HIGH_Q4"])
  ),
  "",
  "8. Outputs",
  paste0("  Master CSV: ", output_paths$master_csv),
  paste0("  Master RDS: ", output_paths$master_rds),
  paste0(
    "  Duplicate resolution table: ",
    output_paths$duplicate_resolution
  ),
  paste0(
    "  Signature availability table: ",
    output_paths$signature_availability
  ),
  paste0("  Metadata JSON: ", output_paths$metadata_json),
  paste0("  QC report: ", output_paths$qc),
  "",
  "9. Warnings",
  if (length(warnings_vec) == 0) {
    "  None"
  } else {
    paste0("  - ", warnings_vec)
  },
  "",
  paste0("Final status: ", final_status)
)

writeLines(
  qc_lines,
  output_paths$qc,
  useBytes = TRUE
)

required_outputs <- unlist(
  output_paths,
  use.names = FALSE
)

missing_outputs <- required_outputs[
  !file.exists(required_outputs)
]

if (length(missing_outputs) > 0) {
  fail(
    "Required outputs were not generated: ",
    paste(missing_outputs, collapse = ", ")
  )
}

msg("Master score table:", output_paths$master_csv)
msg("Duplicate-resolution audit:", output_paths$duplicate_resolution)
msg("Signature availability:", output_paths$signature_availability)
msg("Metadata:", output_paths$metadata_json)
msg("QC:", output_paths$qc)
msg("Final status:", final_status)
