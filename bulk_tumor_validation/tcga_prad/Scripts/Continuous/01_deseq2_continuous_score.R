#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD | Continuous platelet-associated transcriptional score model
#
# Cohort:
#   497 primary tumors, one sample per patient.
#
# Model:
#   DESeq2 design: ~ Score_z
#
# Interpretation:
#   log2FoldChange represents the expression change associated with a
#   one-standard-deviation increase in the platelet-associated score.
#
# Score-associated gene definition:
#   padj < 0.05 and |apeglm-shrunken log2FoldChange| > 0.25.
#
# GSEA:
#   The downstream GSEA analysis must use the complete no-signature table and
#   rank genes by the unshrunken DESeq2 Wald statistic.
#
# This is an unadjusted molecular association model. It is complementary to the
# primary Q1/Q4 transcriptomic analysis and must not be interpreted as causal or
# independent of tumor-microenvironment composition.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 0. Frozen analytical constants
################################################################################

EXPECTED_SCORE_SAMPLES <- 497L
EXPECTED_COUNT_SAMPLES <- 501L
EXPECTED_EXCLUDED_COUNT_SAMPLES <- 4L
EXPECTED_SIGNATURE_GENES <- 41L

PADJ_THRESHOLD <- 0.05
ABS_LOG2FC_THRESHOLD <- 0.25

PREFILTER_MIN_COUNT <- 10L
PREFILTER_MIN_SAMPLES <- 10L

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
  "Scripts", "Continuous", "01_deseq2_continuous_score.R"
)

input_paths <- list(
  score = file.path(
    GENERATED_RESULTS_DIR, "Score",
    "TCGA_primary_master_table_41genes.csv"
  ),
  counts = file.path(
    INPUT_DIR, "LocalLarge",
    "TCGA_PRAD_COUNTS_tumorOnly.rds"
  ),
  gene_map = file.path(
    INPUT_DIR, "Metadata",
    "MAP_ENSGversion_to_geneName.csv"
  ),
  signature = file.path(
    RESOURCE_DIR,
    "platelet_associated_transcriptional_signature.tsv"
  )
)

tables_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous", "Tables", "DESeq2"
)

objects_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous", "Objects"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous", "Logs"
)

invisible(lapply(
  c(tables_dir, objects_dir, logs_dir),
  dir.create,
  showWarnings = FALSE,
  recursive = TRUE
))

output_paths <- list(
  full_results = file.path(
    tables_dir,
    "DESeq2_continuous_score_full_results.csv"
  ),
  no_signature = file.path(
    tables_dir,
    "DESeq2_continuous_score_no_signature_genes.csv"
  ),
  associated_no_signature = file.path(
    tables_dir,
    paste0(
      "DESeq2_continuous_score_no_signature_associated_",
      "FDR005_abslog2FC025.csv"
    )
  ),
  sample_metadata = file.path(
    tables_dir,
    "DESeq2_continuous_score_sample_metadata.csv"
  ),
  summary = file.path(
    tables_dir,
    "DESeq2_continuous_score_summary.csv"
  ),
  metadata = file.path(
    tables_dir,
    "DESeq2_continuous_score_metadata.json"
  ),
  dds = file.path(
    objects_dir,
    "DESeq2_continuous_score_dds.rds"
  ),
  qc = file.path(
    logs_dir,
    "DESeq2_continuous_score_QC.txt"
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

fail_with_qc <- function(message_text) {
  failure_lines <- c(
    "TCGA-PRAD continuous platelet score DESeq2 QC",
    "================================================",
    "",
    paste0("Date/time: ", as.character(Sys.time())),
    paste0("Project dir: ", project_dir),
    paste0("Script: ", script_relative_path),
    paste0("Failure reason: ", message_text),
    "",
    "Final status: FAIL"
  )

  writeLines(failure_lines, output_paths$qc, useBytes = TRUE)
  stop(message_text, call. = FALSE)
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
      "Could not detect the gene-symbol column in the signature. ",
      "Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }

  names(df)[hit[1]]
}

summarize_results <- function(
    df,
    result_set,
    n_samples,
    genes_before,
    genes_after,
    score_summary
) {
  data.frame(
    result_set = result_set,
    n_genes_tested = nrow(df),
    n_padj_lt_0.05 = sum(
      !is.na(df$padj) & df$padj < PADJ_THRESHOLD
    ),
    n_score_associated_FDR005_abslog2FC025 = sum(
      df$score_association %in% c(
        "POSITIVE_ASSOCIATION",
        "NEGATIVE_ASSOCIATION"
      ),
      na.rm = TRUE
    ),
    n_positive_association = sum(
      df$score_association == "POSITIVE_ASSOCIATION",
      na.rm = TRUE
    ),
    n_negative_association = sum(
      df$score_association == "NEGATIVE_ASSOCIATION",
      na.rm = TRUE
    ),
    n_signature_gene_rows = sum(
      df$is_signature_gene,
      na.rm = TRUE
    ),
    n_samples = n_samples,
    n_genes_before_prefilter = genes_before,
    n_genes_after_prefilter = genes_after,
    Score_z_min = unname(score_summary["min"]),
    Score_z_mean = unname(score_summary["mean"]),
    Score_z_sd = unname(score_summary["sd"]),
    Score_z_max = unname(score_summary["max"]),
    stringsAsFactors = FALSE
  )
}

################################################################################
# 3. Required packages
################################################################################

required_packages <- c(
  "DESeq2",
  "apeglm",
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
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      ". The script does not install packages."
    ),
    call. = FALSE
  )
}

################################################################################
# 4. Validate required inputs
################################################################################

msg("Project dir:", project_dir)

for (nm in names(input_paths)) {
  require_input(input_paths[[nm]], nm)
}

################################################################################
# 5. Load and validate the frozen score cohort
################################################################################

msg("[1/8] Loading the canonical 497-patient score cohort")

score_df <- read_required_csv(
  input_paths$score,
  "canonical 41-gene score master table"
)

required_score_columns <- c(
  "sample_id",
  "patient_id",
  "sample_type",
  "Score_raw",
  "Score_z"
)

missing_score_columns <- setdiff(
  required_score_columns,
  names(score_df)
)

if (length(missing_score_columns) > 0) {
  fail(
    "Score table is missing required columns: ",
    paste(missing_score_columns, collapse = ", ")
  )
}

score_df$sample_id <- trimws(as.character(score_df$sample_id))
score_df$patient_id <- trimws(as.character(score_df$patient_id))
score_df$sample_type <- trimws(as.character(score_df$sample_type))

score_df$Score_raw <- suppressWarnings(
  as.numeric(score_df$Score_raw)
)

score_df$Score_z <- suppressWarnings(
  as.numeric(score_df$Score_z)
)

if (nrow(score_df) != EXPECTED_SCORE_SAMPLES) {
  fail(
    "Canonical score table must contain exactly ",
    EXPECTED_SCORE_SAMPLES,
    " samples. Observed: ",
    nrow(score_df)
  )
}

if (
  any(is.na(score_df$sample_id)) ||
  any(score_df$sample_id == "")
) {
  fail("Score table contains missing or empty sample IDs.")
}

if (
  any(is.na(score_df$patient_id)) ||
  any(score_df$patient_id == "")
) {
  fail("Score table contains missing or empty patient IDs.")
}

if (anyDuplicated(score_df$sample_id)) {
  fail("Score table contains duplicated sample IDs.")
}

if (anyDuplicated(score_df$patient_id)) {
  fail("Score table contains duplicated patient IDs.")
}

if (!all(score_df$sample_type == "Primary Tumor")) {
  fail("The continuous cohort contains non-primary tumor samples.")
}

if (
  any(!is.finite(score_df$Score_raw)) ||
  any(!is.finite(score_df$Score_z))
) {
  fail("Score table contains non-finite Score_raw or Score_z values.")
}

score_summary <- c(
  min = min(score_df$Score_z),
  mean = mean(score_df$Score_z),
  sd = stats::sd(score_df$Score_z),
  max = max(score_df$Score_z)
)

if (abs(score_summary["mean"]) > 1e-8) {
  fail(
    "Score_z mean differs from zero beyond tolerance: ",
    format(score_summary["mean"], digits = 10)
  )
}

if (abs(score_summary["sd"] - 1) > 1e-8) {
  fail(
    "Score_z standard deviation differs from one beyond tolerance: ",
    format(score_summary["sd"], digits = 10)
  )
}

msg(
  "Score_z summary:",
  paste(
    names(score_summary),
    format(score_summary, digits = 7),
    sep = "=",
    collapse = "; "
  )
)

################################################################################
# 6. Load and align the raw count matrix
################################################################################

msg("[2/8] Loading and aligning the raw count matrix")

counts <- readRDS(input_paths$counts)

if (is.data.frame(counts)) {
  counts <- as.matrix(counts)
}

if (is.null(dim(counts)) || length(dim(counts)) != 2) {
  fail("Counts input is not a two-dimensional matrix-like object.")
}

counts_dimensions_before_alignment <- dim(counts)

if (ncol(counts) != EXPECTED_COUNT_SAMPLES) {
  fail(
    "Expected ", EXPECTED_COUNT_SAMPLES,
    " columns in the tumor-only count matrix. Observed: ",
    ncol(counts)
  )
}

if (
  is.null(colnames(counts)) ||
  any(is.na(colnames(counts))) ||
  any(colnames(counts) == "")
) {
  fail("Counts matrix contains missing or empty sample identifiers.")
}

if (
  is.null(rownames(counts)) ||
  any(is.na(rownames(counts))) ||
  any(rownames(counts) == "")
) {
  fail("Counts matrix contains missing or empty Ensembl identifiers.")
}

if (anyDuplicated(colnames(counts))) {
  fail("Counts matrix contains duplicated sample columns.")
}

if (anyDuplicated(rownames(counts))) {
  fail("Counts matrix contains duplicated Ensembl rows.")
}

missing_in_counts <- setdiff(
  score_df$sample_id,
  colnames(counts)
)

if (length(missing_in_counts) > 0) {
  fail(
    "Canonical score samples are missing from the counts matrix: ",
    paste(missing_in_counts, collapse = ", ")
  )
}

excluded_count_samples <- setdiff(
  colnames(counts),
  score_df$sample_id
)

if (
  length(excluded_count_samples) !=
  EXPECTED_EXCLUDED_COUNT_SAMPLES
) {
  fail(
    "Expected ", EXPECTED_EXCLUDED_COUNT_SAMPLES,
    " count-matrix samples outside the frozen 497-patient cohort. Observed: ",
    length(excluded_count_samples)
  )
}

counts_aligned <- counts[
  ,
  match(score_df$sample_id, colnames(counts)),
  drop = FALSE
]

if (!identical(colnames(counts_aligned), score_df$sample_id)) {
  fail("Counts alignment failed after matching by sample ID.")
}

counts_dimensions_after_alignment <- dim(counts_aligned)

rm(counts)
gc(verbose = FALSE)

if (any(!is.finite(counts_aligned))) {
  fail("Aligned count matrix contains non-finite values.")
}

if (any(counts_aligned < 0)) {
  fail("Aligned count matrix contains negative values.")
}

non_integer_values <- sum(
  abs(counts_aligned - round(counts_aligned)) > 1e-8
)

if (non_integer_values > 0) {
  fail(
    "Count matrix contains ",
    non_integer_values,
    " non-integer values."
  )
}

counts_aligned <- round(counts_aligned)
storage.mode(counts_aligned) <- "integer"

msg(
  "Counts dimensions before alignment:",
  paste(counts_dimensions_before_alignment, collapse = " x ")
)

msg(
  "Counts dimensions after alignment:",
  paste(counts_dimensions_after_alignment, collapse = " x ")
)

msg(
  "Excluded duplicate-patient count samples:",
  paste(excluded_count_samples, collapse = ", ")
)

################################################################################
# 7. Load the signature and gene map
################################################################################

msg("[3/8] Loading the official signature and gene map")

signature_df <- read_canonical_platelet_signature()

signature_gene_column <- detect_gene_column(signature_df)

signature_genes <- trimws(
  as.character(signature_df[[signature_gene_column]])
)

signature_genes <- signature_genes[
  !is.na(signature_genes) &
    signature_genes != ""
]

signature_keys <- toupper(signature_genes)

if (anyDuplicated(signature_keys)) {
  fail(
    "Official signature contains duplicated genes after ",
    "case-insensitive normalization."
  )
}

if (length(signature_keys) != EXPECTED_SIGNATURE_GENES) {
  fail(
    "Official signature must contain exactly ",
    EXPECTED_SIGNATURE_GENES,
    " genes. Observed: ",
    length(signature_keys)
  )
}

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

map_df$ensg_version <- trimws(
  as.character(map_df$ensg_version)
)

map_df$gene_name <- trimws(
  as.character(map_df$gene_name)
)

map_df <- map_df[
  !is.na(map_df$ensg_version) &
    map_df$ensg_version != "",
  ,
  drop = FALSE
]

if (anyDuplicated(map_df$ensg_version)) {
  duplicated_ids <- unique(
    map_df$ensg_version[
      duplicated(map_df$ensg_version)
    ]
  )

  conflicting_ids <- duplicated_ids[
    vapply(
      duplicated_ids,
      function(id) {
        mapped_symbols <- unique(
          map_df$gene_name[
            map_df$ensg_version == id &
              !is.na(map_df$gene_name) &
              map_df$gene_name != ""
          ]
        )

        length(mapped_symbols) > 1
      },
      logical(1)
    )
  ]

  if (length(conflicting_ids) > 0) {
    fail(
      "Gene map contains Ensembl IDs mapped to multiple symbols. Examples: ",
      paste(utils::head(conflicting_ids, 20), collapse = ", ")
    )
  }

  map_df <- map_df[
    !duplicated(map_df$ensg_version),
    ,
    drop = FALSE
  ]
}

################################################################################
# 8. Prepare colData and apply the expression prefilter
################################################################################

msg("[4/8] Preparing continuous-model sample metadata")

coldata <- data.frame(
  sample_id = score_df$sample_id,
  patient_id = score_df$patient_id,
  sample_type = score_df$sample_type,
  Score_raw = score_df$Score_raw,
  Score_z = score_df$Score_z,
  row.names = score_df$sample_id,
  stringsAsFactors = FALSE
)

if ("platelet_score_group" %in% names(score_df)) {
  coldata$platelet_score_group <- as.character(
    score_df$platelet_score_group
  )
}

if (!identical(rownames(coldata), colnames(counts_aligned))) {
  fail("DESeq2 colData and aligned counts are not identically ordered.")
}

write_csv(
  coldata,
  output_paths$sample_metadata
)

msg(
  "[5/8] Applying prefilter: counts >= ",
  PREFILTER_MIN_COUNT,
  " in at least ",
  PREFILTER_MIN_SAMPLES,
  " samples"
)

n_genes_before_prefilter <- nrow(counts_aligned)

keep <- rowSums(
  counts_aligned >= PREFILTER_MIN_COUNT
) >= PREFILTER_MIN_SAMPLES

counts_prefiltered <- counts_aligned[
  keep,
  ,
  drop = FALSE
]

n_genes_after_prefilter <- nrow(counts_prefiltered)
n_genes_removed_prefilter <- n_genes_before_prefilter -
  n_genes_after_prefilter

if (n_genes_after_prefilter < 1000) {
  fail(
    "Too few genes remained after prefiltering: ",
    n_genes_after_prefilter
  )
}

msg("Genes before prefilter:", n_genes_before_prefilter)
msg("Genes after prefilter:", n_genes_after_prefilter)
msg("Genes removed by prefilter:", n_genes_removed_prefilter)

rm(counts_aligned)
gc(verbose = FALSE)

################################################################################
# 9. Run the continuous DESeq2 model
################################################################################

msg("[6/8] Running DESeq2 continuous model: ~ Score_z")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_prefiltered,
  colData = coldata,
  design = ~ Score_z
)

dds <- DESeq2::DESeq(
  dds,
  quiet = FALSE
)

results_names <- DESeq2::resultsNames(dds)
expected_coefficient <- "Score_z"

if (!expected_coefficient %in% results_names) {
  fail(
    "Expected DESeq2 coefficient was not found: ",
    expected_coefficient,
    ". Available coefficients: ",
    paste(results_names, collapse = ", ")
  )
}

coefficient_name <- expected_coefficient

res <- DESeq2::results(
  dds,
  name = coefficient_name,
  alpha = PADJ_THRESHOLD
)

msg(
  "Running apeglm shrinkage for coefficient:",
  coefficient_name
)

res_shrunken <- tryCatch(
  DESeq2::lfcShrink(
    dds,
    coef = coefficient_name,
    res = res,
    type = "apeglm"
  ),
  error = function(e) e
)

if (inherits(res_shrunken, "error")) {
  fail(
    "apeglm shrinkage failed: ",
    conditionMessage(res_shrunken)
  )
}

res_df <- as.data.frame(res)
res_shrunken_df <- as.data.frame(res_shrunken)

if (!identical(rownames(res_df), rownames(res_shrunken_df))) {
  fail(
    "Shrunken and unshrunken DESeq2 results are not identically aligned."
  )
}

res_df$ensembl_id <- rownames(res_df)

res_df$log2FoldChange_shrunken <- as.numeric(
  res_shrunken_df$log2FoldChange
)

if (any(!is.finite(res_df$stat))) {
  fail(
    "Continuous DESeq2 Wald statistics contain non-finite values."
  )
}

if (any(!is.finite(res_df$log2FoldChange_shrunken))) {
  fail(
    "Continuous apeglm-shrunken log2 fold changes contain non-finite values."
  )
}

saveRDS(
  dds,
  output_paths$dds
)

################################################################################
# 10. Annotate and classify score-associated genes
################################################################################

msg("[7/8] Annotating continuous-model results")

map_aligned <- map_df[
  match(res_df$ensembl_id, map_df$ensg_version),
  ,
  drop = FALSE
]

gene_name_raw <- trimws(
  as.character(map_aligned$gene_name)
)

missing_gene_name <- (
  is.na(gene_name_raw) |
    gene_name_raw == ""
)

gene_name <- gene_name_raw
gene_name[missing_gene_name] <-
  res_df$ensembl_id[missing_gene_name]

is_signature_gene <- (
  !missing_gene_name &
    toupper(gene_name) %in% signature_keys
)

score_association <- ifelse(
  !is.na(res_df$padj) &
    res_df$padj < PADJ_THRESHOLD &
    res_df$log2FoldChange_shrunken > ABS_LOG2FC_THRESHOLD,
  "POSITIVE_ASSOCIATION",
  ifelse(
    !is.na(res_df$padj) &
      res_df$padj < PADJ_THRESHOLD &
      res_df$log2FoldChange_shrunken < -ABS_LOG2FC_THRESHOLD,
    "NEGATIVE_ASSOCIATION",
    "NS"
  )
)

full_results <- data.frame(
  ensembl_id = res_df$ensembl_id,
  gene_name = gene_name,
  gene_name_is_ensembl_fallback = missing_gene_name,
  baseMean = as.numeric(res_df$baseMean),
  log2FoldChange = as.numeric(res_df$log2FoldChange),
  lfcSE = as.numeric(res_df$lfcSE),
  stat = as.numeric(res_df$stat),
  pvalue = as.numeric(res_df$pvalue),
  padj = as.numeric(res_df$padj),
  log2FoldChange_shrunken = as.numeric(
    res_df$log2FoldChange_shrunken
  ),
  shrinkage_method = "apeglm",
  is_signature_gene = is_signature_gene,
  score_association = score_association,
  stringsAsFactors = FALSE
)

full_results <- full_results[
  order(
    full_results$padj,
    -abs(full_results$log2FoldChange_shrunken),
    full_results$gene_name,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

no_signature_results <- full_results[
  !full_results$is_signature_gene,
  ,
  drop = FALSE
]

associated_no_signature <- no_signature_results[
  no_signature_results$score_association %in% c(
    "POSITIVE_ASSOCIATION",
    "NEGATIVE_ASSOCIATION"
  ),
  ,
  drop = FALSE
]

n_signature_rows_in_results <- sum(
  full_results$is_signature_gene,
  na.rm = TRUE
)

n_unique_signature_genes_tested <- length(
  unique(
    toupper(
      full_results$gene_name[
        full_results$is_signature_gene
      ]
    )
  )
)

n_missing_gene_names <- sum(
  full_results$gene_name_is_ensembl_fallback,
  na.rm = TRUE
)

if (
  n_unique_signature_genes_tested !=
  EXPECTED_SIGNATURE_GENES
) {
  fail(
    "Expected all ", EXPECTED_SIGNATURE_GENES,
    " signature genes to be represented in the continuous results. Observed: ",
    n_unique_signature_genes_tested
  )
}

n_positive <- sum(
  associated_no_signature$score_association ==
    "POSITIVE_ASSOCIATION"
)

n_negative <- sum(
  associated_no_signature$score_association ==
    "NEGATIVE_ASSOCIATION"
)

n_fdr_full <- sum(
  !is.na(full_results$padj) &
    full_results$padj < PADJ_THRESHOLD
)

n_fdr_no_signature <- sum(
  !is.na(no_signature_results$padj) &
    no_signature_results$padj < PADJ_THRESHOLD
)

summary_df <- rbind(
  summarize_results(
    full_results,
    "full",
    nrow(score_df),
    n_genes_before_prefilter,
    n_genes_after_prefilter,
    score_summary
  ),
  summarize_results(
    no_signature_results,
    "no_signature",
    nrow(score_df),
    n_genes_before_prefilter,
    n_genes_after_prefilter,
    score_summary
  ),
  summarize_results(
    associated_no_signature,
    "associated_no_signature",
    nrow(score_df),
    n_genes_before_prefilter,
    n_genes_after_prefilter,
    score_summary
  )
)

write_csv(
  full_results,
  output_paths$full_results
)

write_csv(
  no_signature_results,
  output_paths$no_signature
)

write_csv(
  associated_no_signature,
  output_paths$associated_no_signature
)

write_csv(
  summary_df,
  output_paths$summary
)

################################################################################
# 11. Metadata
################################################################################

msg("[8/8] Writing metadata and QC")

final_status <- if (length(warnings_vec) == 0) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

metadata <- list(
  date_time = as.character(Sys.time()),
  project_dir = project_dir,
  script = script_relative_path,
  project = "tcga_prad",
  analysis = "continuous platelet-associated transcriptional score model",
  analytic_hierarchy = "complementary robustness analysis",
  interpretation = paste0(
    "Unadjusted molecular association model. The coefficient represents ",
    "expression change per one-standard-deviation increase in Score_z. ",
    "Results are not causal and are not independent of TME composition."
  ),
  model = list(
    method = "DESeq2",
    design = "~ Score_z",
    coefficient = coefficient_name,
    statistical_test = "Wald test",
    coefficient_unit = "one standard deviation increase in Score_z"
  ),
  effect_size = list(
    unshrunken_column = "log2FoldChange",
    shrunken_column = "log2FoldChange_shrunken",
    shrinkage_method = "apeglm",
    interpretation_column = "log2FoldChange_shrunken"
  ),
  thresholds = list(
    padj = PADJ_THRESHOLD,
    absolute_shrunken_log2FC = ABS_LOG2FC_THRESHOLD,
    terminology = "score-associated genes, not DEGs"
  ),
  GSEA = list(
    input_table = output_paths$no_signature,
    ranking_metric = "stat",
    ranking_source = "unshrunken DESeq2 Wald statistic",
    thresholded_gene_table_used = FALSE
  ),
  prefilter = list(
    rule = paste0(
      "rowSums(counts >= ",
      PREFILTER_MIN_COUNT,
      ") >= ",
      PREFILTER_MIN_SAMPLES
    ),
    minimum_count = PREFILTER_MIN_COUNT,
    minimum_samples = PREFILTER_MIN_SAMPLES,
    genes_before = n_genes_before_prefilter,
    genes_after = n_genes_after_prefilter,
    genes_removed = n_genes_removed_prefilter
  ),
  inputs = input_paths,
  cohort = list(
    primary_tumors = nrow(score_df),
    unique_patients = length(unique(score_df$patient_id)),
    count_matrix_samples = counts_dimensions_before_alignment[2],
    excluded_duplicate_samples = excluded_count_samples
  ),
  Score_z = as.list(score_summary),
  count_dimensions = list(
    before_alignment = as.integer(
      counts_dimensions_before_alignment
    ),
    after_alignment = as.integer(
      counts_dimensions_after_alignment
    )
  ),
  results = list(
    genes_tested = nrow(full_results),
    genes_without_signature = nrow(no_signature_results),
    genes_FDR005_full = n_fdr_full,
    genes_FDR005_no_signature = n_fdr_no_signature,
    associated_no_signature_FDR005_abslog2FC025 =
      nrow(associated_no_signature),
    positive_score_association = n_positive,
    negative_score_association = n_negative,
    signature_gene_rows_in_results =
      n_signature_rows_in_results,
    unique_signature_genes_tested =
      n_unique_signature_genes_tested,
    rows_using_Ensembl_fallback =
      n_missing_gene_names
  ),
  outputs = output_paths,
  reproducibility = list(
    R_version = R.version.string,
    DESeq2_version = as.character(
      utils::packageVersion("DESeq2")
    ),
    apeglm_version = as.character(
      utils::packageVersion("apeglm")
    ),
    jsonlite_version = as.character(
      utils::packageVersion("jsonlite")
    )
  ),
  warnings = warnings_vec,
  qc_status = final_status
)

jsonlite::write_json(
  metadata,
  output_paths$metadata,
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null",
  digits = NA
)

################################################################################
# 12. Final QC report
################################################################################

qc_lines <- c(
  "TCGA-PRAD continuous platelet score DESeq2 QC",
  "================================================",
  "",
  paste0("Date/time: ", as.character(Sys.time())),
  paste0("Project dir: ", project_dir),
  paste0("Script: ", script_relative_path),
  "",
  "1. Inputs",
  paste0("  Canonical score table: ", input_paths$score),
  paste0("  Count matrix: ", input_paths$counts),
  paste0("  Gene map: ", input_paths$gene_map),
  paste0("  Official signature: ", input_paths$signature),
  "",
  "2. Cohort",
  paste0("  Primary tumors: ", nrow(score_df)),
  paste0(
    "  Unique patients: ",
    length(unique(score_df$patient_id))
  ),
  paste0(
    "  Count-matrix samples before alignment: ",
    counts_dimensions_before_alignment[2]
  ),
  paste0(
    "  Samples excluded during alignment: ",
    length(excluded_count_samples)
  ),
  paste0(
    "  Excluded sample IDs: ",
    paste(excluded_count_samples, collapse = ", ")
  ),
  paste0(
    "  Count-matrix samples after alignment: ",
    counts_dimensions_after_alignment[2]
  ),
  "",
  "3. Score_z",
  paste0(
    "  Minimum: ",
    format(score_summary["min"], digits = 10)
  ),
  paste0(
    "  Mean: ",
    format(score_summary["mean"], digits = 10)
  ),
  paste0(
    "  SD: ",
    format(score_summary["sd"], digits = 10)
  ),
  paste0(
    "  Maximum: ",
    format(score_summary["max"], digits = 10)
  ),
  "",
  "4. Prefilter",
  paste0(
    "  Rule: rowSums(counts >= ",
    PREFILTER_MIN_COUNT,
    ") >= ",
    PREFILTER_MIN_SAMPLES
  ),
  paste0("  Genes before prefilter: ", n_genes_before_prefilter),
  paste0("  Genes after prefilter: ", n_genes_after_prefilter),
  paste0("  Genes removed: ", n_genes_removed_prefilter),
  "",
  "5. Continuous DESeq2 model",
  "  Design: ~ Score_z",
  paste0("  Coefficient: ", coefficient_name),
  "  Statistical test: Wald test",
  "  Coefficient unit: one-standard-deviation increase in Score_z",
  "  Shrinkage method: apeglm",
  "",
  "6. Interpretation",
  "  This is a complementary unadjusted molecular association model.",
  "  It is not the primary Q1/Q4 transcriptomic contrast.",
  "  Results are termed score-associated genes, not DEGs.",
  "  Results are not causal and are not independent of TME composition.",
  "",
  "7. Thresholds",
  paste0("  padj < ", PADJ_THRESHOLD),
  paste0(
    "  |log2FoldChange_shrunken| > ",
    ABS_LOG2FC_THRESHOLD
  ),
  "",
  "8. Signature exclusion",
  paste0("  Official signature genes: ", length(signature_keys)),
  paste0(
    "  Signature-gene rows in full results: ",
    n_signature_rows_in_results
  ),
  paste0(
    "  Unique signature genes tested: ",
    n_unique_signature_genes_tested
  ),
  paste0(
    "  Rows using Ensembl ID as gene-name fallback: ",
    n_missing_gene_names
  ),
  "",
  "9. Results",
  paste0("  Genes tested: ", nrow(full_results)),
  paste0(
    "  Genes after signature exclusion: ",
    nrow(no_signature_results)
  ),
  paste0(
    "  Genes FDR < 0.05, full: ",
    n_fdr_full
  ),
  paste0(
    "  Genes FDR < 0.05, no signature: ",
    n_fdr_no_signature
  ),
  paste0(
    "  Score-associated genes, no signature, FDR < 0.05 and ",
    "|shrunken log2FC| > 0.25: ",
    nrow(associated_no_signature)
  ),
  paste0(
    "  Positive association with Score_z: ",
    n_positive
  ),
  paste0(
    "  Negative association with Score_z: ",
    n_negative
  ),
  "",
  "10. GSEA contract",
  paste0(
    "  Input table: ",
    output_paths$no_signature
  ),
  "  Ranking metric: unshrunken DESeq2 Wald statistic (stat)",
  "  Thresholded associated-gene table must not be used for GSEA.",
  "",
  "11. Summary table",
  paste(
    capture.output(
      print(summary_df, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "12. Outputs",
  paste0("  Full results: ", output_paths$full_results),
  paste0("  No-signature results: ", output_paths$no_signature),
  paste0(
    "  Associated no-signature genes: ",
    output_paths$associated_no_signature
  ),
  paste0("  Sample metadata: ", output_paths$sample_metadata),
  paste0("  Summary: ", output_paths$summary),
  paste0("  Metadata: ", output_paths$metadata),
  paste0("  DESeq2 object: ", output_paths$dds),
  paste0("  QC: ", output_paths$qc),
  "",
  "13. Warnings",
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

msg("Full results:", output_paths$full_results)
msg("No-signature results:", output_paths$no_signature)
msg(
  "Associated no-signature genes:",
  output_paths$associated_no_signature
)
msg("Sample metadata:", output_paths$sample_metadata)
msg("DESeq2 object:", output_paths$dds)
msg("Metadata:", output_paths$metadata)
msg("QC:", output_paths$qc)
msg("Final status:", final_status)
