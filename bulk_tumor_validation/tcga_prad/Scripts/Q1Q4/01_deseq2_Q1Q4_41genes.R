#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD | Q1 vs Q4 differential expression analysis
#
# Score:
#   Official 41-gene platelet-associated transcriptional score.
#
# Cohort:
#   Final frozen cohort of 497 primary tumors, one sample per patient.
#   Q1: 125 samples.
#   Q4: 125 samples.
#
# Differential-expression model:
#   DESeq2 design: ~ platelet_score_group
#   Contrast: HIGH_Q4 vs LOW_Q1
#
# Interpretation:
#   Statistical testing and GSEA ranking use the unshrunken DESeq2 Wald test.
#   DEG effect-size classification uses apeglm-shrunken log2 fold change.
#
# DEG definition:
#   padj < 0.05 and |shrunken log2FC| > 1.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 0. Frozen analytical constants
################################################################################

EXPECTED_SCORE_SAMPLES <- 497L
EXPECTED_LOW_Q1 <- 125L
EXPECTED_HIGH_Q4 <- 125L
EXPECTED_Q1Q4_TOTAL <- 250L
EXPECTED_SIGNATURE_GENES <- 41L

PADJ_THRESHOLD <- 0.05
ABS_LOG2FC_THRESHOLD <- 1

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
  "Scripts", "Q1Q4", "01_deseq2_Q1Q4_41genes.R"
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
  GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "DESeq2"
)

objects_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Objects"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Logs"
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
    "DESeq2_Q1Q4_41genes_full_results.csv"
  ),
  no_signature = file.path(
    tables_dir,
    "DESeq2_Q1Q4_41genes_no_signature_genes.csv"
  ),
  significant_no_signature = file.path(
    tables_dir,
    "DESeq2_Q1Q4_41genes_no_signature_DEG_FDR005_abslog2FC1.csv"
  ),
  sample_metadata = file.path(
    tables_dir,
    "DESeq2_Q1Q4_sample_metadata.csv"
  ),
  summary = file.path(
    tables_dir,
    "DESeq2_Q1Q4_41genes_summary.csv"
  ),
  metadata = file.path(
    tables_dir,
    "DESeq2_Q1Q4_41genes_metadata.json"
  ),
  dds = file.path(
    objects_dir,
    "DESeq2_Q1Q4_dds.rds"
  ),
  vst = file.path(
    objects_dir,
    "DESeq2_Q1Q4_vst.rds"
  ),
  qc = file.path(
    logs_dir,
    "DESeq2_Q1Q4_41genes_QC.txt"
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
    "TCGA-PRAD DESeq2 Q1 vs Q4 QC",
    "================================",
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
    group_counts,
    genes_before,
    genes_after
) {
  data.frame(
    result_set = result_set,
    n_genes_tested = nrow(df),
    n_padj_lt_0.05 = sum(
      !is.na(df$padj) & df$padj < PADJ_THRESHOLD
    ),
    n_DEG_FDR005_abslog2FC1 = sum(
      df$regulation %in% c("UP", "DOWN"),
      na.rm = TRUE
    ),
    n_up_FDR005_abslog2FC1 = sum(
      df$regulation == "UP",
      na.rm = TRUE
    ),
    n_down_FDR005_abslog2FC1 = sum(
      df$regulation == "DOWN",
      na.rm = TRUE
    ),
    n_signature_gene_rows = sum(
      df$is_signature_gene,
      na.rm = TRUE
    ),
    n_LOW_Q1 = as.integer(group_counts["LOW_Q1"]),
    n_HIGH_Q4 = as.integer(group_counts["HIGH_Q4"]),
    n_samples_total = sum(group_counts),
    n_genes_before_prefilter = genes_before,
    n_genes_after_prefilter = genes_after,
    stringsAsFactors = FALSE
  )
}

################################################################################
# 3. Required packages
################################################################################

required_packages <- c(
  "DESeq2",
  "apeglm",
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
# 5. Load and validate the frozen score cohort
################################################################################

msg("[1/8] Loading the official score cohort")

score_df <- read_required_csv(
  input_paths$score,
  "official 41-gene score master table"
)

required_score_columns <- c(
  "sample_id",
  "patient_id",
  "Score_raw",
  "Score_z",
  "platelet_score_group"
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
score_df$Score_raw <- suppressWarnings(as.numeric(score_df$Score_raw))
score_df$Score_z <- suppressWarnings(as.numeric(score_df$Score_z))
score_df$platelet_score_group <- as.character(
  score_df$platelet_score_group
)

if (nrow(score_df) != EXPECTED_SCORE_SAMPLES) {
  fail(
    "Official score table must contain exactly ",
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

if (
  any(!is.finite(score_df$Score_raw)) ||
  any(!is.finite(score_df$Score_z))
) {
  fail("Score table contains non-finite Score_raw or Score_z values.")
}

allowed_groups <- c("LOW_Q1", "MID", "HIGH_Q4")

unexpected_groups <- setdiff(
  unique(score_df$platelet_score_group),
  allowed_groups
)

if (length(unexpected_groups) > 0) {
  fail(
    "Unexpected platelet score groups were found: ",
    paste(unexpected_groups, collapse = ", ")
  )
}

score_q <- score_df[
  score_df$platelet_score_group %in% c("LOW_Q1", "HIGH_Q4"),
  ,
  drop = FALSE
]

score_q$platelet_score_group <- factor(
  score_q$platelet_score_group,
  levels = c("LOW_Q1", "HIGH_Q4")
)

group_counts <- table(score_q$platelet_score_group)

if (
  as.integer(group_counts["LOW_Q1"]) != EXPECTED_LOW_Q1 ||
  as.integer(group_counts["HIGH_Q4"]) != EXPECTED_HIGH_Q4 ||
  nrow(score_q) != EXPECTED_Q1Q4_TOTAL
) {
  fail(
    "Unexpected Q1/Q4 sample counts. Expected LOW_Q1=",
    EXPECTED_LOW_Q1,
    ", HIGH_Q4=",
    EXPECTED_HIGH_Q4,
    ", total=",
    EXPECTED_Q1Q4_TOTAL,
    ". Observed: ",
    paste(
      names(group_counts),
      as.integer(group_counts),
      sep = "=",
      collapse = "; "
    ),
    "; total=",
    nrow(score_q)
  )
}

if (anyDuplicated(score_q$patient_id)) {
  fail("The Q1/Q4 cohort contains duplicated patients.")
}

msg(
  "Q1/Q4 samples:",
  paste(
    names(group_counts),
    as.integer(group_counts),
    sep = "=",
    collapse = "; "
  )
)

################################################################################
# 6. Load and align the count matrix
################################################################################

msg("[2/8] Loading and aligning the unstranded count matrix")

counts <- readRDS(input_paths$counts)

if (is.data.frame(counts)) {
  counts <- as.matrix(counts)
}

if (is.null(dim(counts)) || length(dim(counts)) != 2) {
  fail("Counts input is not a two-dimensional matrix-like object.")
}

counts_dimensions_before_alignment <- dim(counts)

if (
  is.null(colnames(counts)) ||
  any(is.na(colnames(counts))) ||
  any(colnames(counts) == "")
) {
  fail("Counts matrix must contain non-empty sample column names.")
}

if (
  is.null(rownames(counts)) ||
  any(is.na(rownames(counts))) ||
  any(rownames(counts) == "")
) {
  fail("Counts matrix must contain non-empty Ensembl row identifiers.")
}

if (anyDuplicated(colnames(counts))) {
  fail("Counts matrix contains duplicated sample columns.")
}

if (anyDuplicated(rownames(counts))) {
  fail("Counts matrix contains duplicated Ensembl row identifiers.")
}

missing_in_counts <- setdiff(
  score_q$sample_id,
  colnames(counts)
)

if (length(missing_in_counts) > 0) {
  fail(
    "Q1/Q4 samples are missing from the counts matrix: ",
    paste(missing_in_counts, collapse = ", ")
  )
}

counts_q <- counts[
  ,
  match(score_q$sample_id, colnames(counts)),
  drop = FALSE
]

if (!identical(colnames(counts_q), score_q$sample_id)) {
  fail("Counts alignment failed after matching by sample ID.")
}

counts_dimensions_after_alignment <- dim(counts_q)

rm(counts)
gc(verbose = FALSE)

if (any(!is.finite(counts_q))) {
  fail("The aligned count matrix contains non-finite values.")
}

if (any(counts_q < 0)) {
  fail("The aligned count matrix contains negative values.")
}

non_integer_count_values <- sum(
  abs(counts_q - round(counts_q)) > 1e-8
)

if (non_integer_count_values > 0) {
  fail(
    "The count matrix contains ",
    non_integer_count_values,
    " non-integer values."
  )
}

counts_q <- round(counts_q)
storage.mode(counts_q) <- "integer"

################################################################################
# 7. Load the signature and gene map
################################################################################

msg("[3/8] Loading the signature and gene map")

signature_df <- read_canonical_platelet_signature()

signature_gene_column <- detect_gene_column(signature_df)

signature_genes <- trimws(
  as.character(signature_df[[signature_gene_column]])
)

signature_genes <- signature_genes[
  !is.na(signature_genes) & signature_genes != ""
]

signature_keys <- toupper(signature_genes)

if (anyDuplicated(signature_keys)) {
  fail(
    "The official signature contains duplicated genes after ",
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

required_map_columns <- c("ensg_version", "gene_name")

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
  duplicated_ids <- unique(
    map_df$ensg_version[duplicated(map_df$ensg_version)]
  )

  conflicting_ids <- duplicated_ids[
    vapply(
      duplicated_ids,
      function(id) {
        mapped_names <- unique(
          map_df$gene_name[
            map_df$ensg_version == id &
              !is.na(map_df$gene_name) &
              map_df$gene_name != ""
          ]
        )

        length(mapped_names) > 1
      },
      logical(1)
    )
  ]

  if (length(conflicting_ids) > 0) {
    fail(
      "Gene map contains Ensembl IDs mapped to multiple gene symbols. ",
      "Examples: ",
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
# 8. Prepare sample metadata and apply the expression prefilter
################################################################################

msg("[4/8] Preparing DESeq2 sample metadata")

coldata <- data.frame(
  sample_id = score_q$sample_id,
  patient_id = score_q$patient_id,
  platelet_score_group = factor(
    as.character(score_q$platelet_score_group),
    levels = c("LOW_Q1", "HIGH_Q4")
  ),
  Score_raw = score_q$Score_raw,
  Score_z = score_q$Score_z,
  row.names = score_q$sample_id,
  stringsAsFactors = FALSE
)

if (!identical(rownames(coldata), colnames(counts_q))) {
  fail("DESeq2 colData and count matrix are not identically aligned.")
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

n_genes_before_prefilter <- nrow(counts_q)

keep <- rowSums(
  counts_q >= PREFILTER_MIN_COUNT
) >= PREFILTER_MIN_SAMPLES

counts_prefiltered <- counts_q[
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

rm(counts_q)
gc(verbose = FALSE)

################################################################################
# 9. Run DESeq2 and apeglm shrinkage
################################################################################

msg("[6/8] Running DESeq2: HIGH_Q4 vs LOW_Q1")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_prefiltered,
  colData = coldata,
  design = ~ platelet_score_group
)

dds <- DESeq2::DESeq(
  dds,
  quiet = FALSE
)

results_names <- DESeq2::resultsNames(dds)

expected_coefficient <- "platelet_score_group_HIGH_Q4_vs_LOW_Q1"

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

msg("Running apeglm shrinkage with coefficient:", coefficient_name)

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
  fail("Shrunken and unshrunken DESeq2 results are not identically aligned.")
}

res_df$ensembl_id <- rownames(res_df)

res_df$log2FoldChange_shrunken <- as.numeric(
  res_shrunken_df$log2FoldChange
)

if (all(!is.finite(res_df$log2FoldChange_shrunken))) {
  fail("All apeglm-shrunken log2 fold changes are non-finite.")
}

################################################################################
# 10. Generate VST object for downstream figures
################################################################################

msg("Generating variance-stabilized expression object")

vst_object <- DESeq2::vst(
  dds,
  blind = FALSE
)

if (!identical(colnames(vst_object), rownames(coldata))) {
  fail("VST object and Q1/Q4 sample metadata are not aligned.")
}

################################################################################
# 11. Annotate and classify DESeq2 results
################################################################################

msg("[7/8] Annotating DESeq2 results")

map_aligned <- map_df[
  match(res_df$ensembl_id, map_df$ensg_version),
  ,
  drop = FALSE
]

gene_name_raw <- trimws(as.character(map_aligned$gene_name))

missing_gene_name <- (
  is.na(gene_name_raw) |
    gene_name_raw == ""
)

gene_name <- gene_name_raw
gene_name[missing_gene_name] <- res_df$ensembl_id[missing_gene_name]

is_signature_gene <- (
  !missing_gene_name &
    toupper(gene_name) %in% signature_keys
)

regulation <- ifelse(
  !is.na(res_df$padj) &
    res_df$padj < PADJ_THRESHOLD &
    res_df$log2FoldChange_shrunken > ABS_LOG2FC_THRESHOLD,
  "UP",
  ifelse(
    !is.na(res_df$padj) &
      res_df$padj < PADJ_THRESHOLD &
      res_df$log2FoldChange_shrunken < -ABS_LOG2FC_THRESHOLD,
    "DOWN",
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
  regulation = regulation,
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

significant_no_signature <- no_signature_results[
  no_signature_results$regulation %in% c("UP", "DOWN"),
  ,
  drop = FALSE
]

n_signature_rows_in_results <- sum(
  full_results$is_signature_gene,
  na.rm = TRUE
)

n_missing_gene_names <- sum(
  full_results$gene_name_is_ensembl_fallback,
  na.rm = TRUE
)

signature_genes_tested <- unique(
  full_results$gene_name[
    full_results$is_signature_gene
  ]
)

n_signature_genes_tested <- length(
  unique(toupper(signature_genes_tested))
)

summary_df <- rbind(
  summarize_results(
    full_results,
    "full",
    group_counts,
    n_genes_before_prefilter,
    n_genes_after_prefilter
  ),
  summarize_results(
    no_signature_results,
    "no_signature",
    group_counts,
    n_genes_before_prefilter,
    n_genes_after_prefilter
  ),
  summarize_results(
    significant_no_signature,
    "significant_no_signature",
    group_counts,
    n_genes_before_prefilter,
    n_genes_after_prefilter
  )
)

################################################################################
# 12. Write tables and analytical objects
################################################################################

msg("[8/8] Writing DESeq2 outputs")

write_csv(
  full_results,
  output_paths$full_results
)

write_csv(
  no_signature_results,
  output_paths$no_signature
)

write_csv(
  significant_no_signature,
  output_paths$significant_no_signature
)

write_csv(
  summary_df,
  output_paths$summary
)

saveRDS(
  dds,
  output_paths$dds
)

saveRDS(
  vst_object,
  output_paths$vst
)

################################################################################
# 13. Metadata
################################################################################

n_significant_fdr_full <- sum(
  !is.na(full_results$padj) &
    full_results$padj < PADJ_THRESHOLD
)

n_significant_fdr_no_signature <- sum(
  !is.na(no_signature_results$padj) &
    no_signature_results$padj < PADJ_THRESHOLD
)

n_up <- sum(
  significant_no_signature$regulation == "UP"
)

n_down <- sum(
  significant_no_signature$regulation == "DOWN"
)

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
  cohort = paste0(
    "TCGA-PRAD primary tumors, one sample per patient; ",
    "Q1 versus Q4 of the official 41-gene score"
  ),
  analytic_hierarchy = "primary transcriptomic analysis",
  model = list(
    method = "DESeq2",
    design = "~ platelet_score_group",
    reference_group = "LOW_Q1",
    contrast = "HIGH_Q4 vs LOW_Q1",
    coefficient = coefficient_name,
    statistical_test = "Wald test"
  ),
  effect_size = list(
    unshrunken_column = "log2FoldChange",
    shrunken_column = "log2FoldChange_shrunken",
    shrinkage_method = "apeglm",
    DEG_interpretation_column = "log2FoldChange_shrunken"
  ),
  thresholds = list(
    padj = PADJ_THRESHOLD,
    absolute_shrunken_log2FC = ABS_LOG2FC_THRESHOLD
  ),
  GSEA = list(
    ranking_metric = "stat",
    ranking_source = "unshrunken DESeq2 Wald statistic",
    input_table = output_paths$no_signature
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
  samples = list(
    frozen_score_cohort = nrow(score_df),
    LOW_Q1 = as.integer(group_counts["LOW_Q1"]),
    HIGH_Q4 = as.integer(group_counts["HIGH_Q4"]),
    total_Q1Q4 = nrow(score_q)
  ),
  count_dimensions = list(
    before_alignment = as.integer(
      counts_dimensions_before_alignment
    ),
    after_Q1Q4_alignment = as.integer(
      counts_dimensions_after_alignment
    )
  ),
  results = list(
    genes_tested = nrow(full_results),
    genes_without_signature = nrow(no_signature_results),
    genes_FDR005_full = n_significant_fdr_full,
    genes_FDR005_no_signature = n_significant_fdr_no_signature,
    DEG_no_signature_FDR005_abslog2FC1 =
      nrow(significant_no_signature),
    DEG_up_HIGH_Q4 = n_up,
    DEG_down_HIGH_Q4 = n_down,
    signature_gene_rows_in_results =
      n_signature_rows_in_results,
    unique_signature_genes_tested =
      n_signature_genes_tested,
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
    SummarizedExperiment_version = as.character(
      utils::packageVersion("SummarizedExperiment")
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
# 14. QC report
################################################################################

qc_lines <- c(
  "TCGA-PRAD DESeq2 Q1 vs Q4 QC",
  "================================",
  "",
  paste0("Date/time: ", as.character(Sys.time())),
  paste0("Project dir: ", project_dir),
  paste0("Script: ", script_relative_path),
  "",
  "1. Inputs",
  paste0("  Score table: ", input_paths$score),
  paste0("  Counts matrix: ", input_paths$counts),
  paste0("  Gene map: ", input_paths$gene_map),
  paste0("  Official signature: ", input_paths$signature),
  "",
  "2. Frozen cohort",
  paste0("  Score cohort: ", nrow(score_df)),
  paste0("  LOW_Q1: ", as.integer(group_counts["LOW_Q1"])),
  paste0("  HIGH_Q4: ", as.integer(group_counts["HIGH_Q4"])),
  paste0("  Total Q1/Q4: ", nrow(score_q)),
  paste0(
    "  Unique patients in Q1/Q4: ",
    length(unique(score_q$patient_id))
  ),
  "",
  "3. Count alignment",
  paste0(
    "  Dimensions before alignment: ",
    paste(counts_dimensions_before_alignment, collapse = " x ")
  ),
  paste0(
    "  Dimensions after Q1/Q4 alignment: ",
    paste(counts_dimensions_after_alignment, collapse = " x ")
  ),
  "  Sample order identical between counts and colData: TRUE",
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
  "5. DESeq2",
  "  Design: ~ platelet_score_group",
  "  Reference group: LOW_Q1",
  "  Contrast: HIGH_Q4 vs LOW_Q1",
  paste0("  Coefficient: ", coefficient_name),
  "  Statistical test: Wald test",
  "  GSEA ranking metric: unshrunken stat",
  "  DEG effect size: apeglm-shrunken log2 fold change",
  "",
  "6. Thresholds",
  paste0("  padj < ", PADJ_THRESHOLD),
  paste0(
    "  |log2FoldChange_shrunken| > ",
    ABS_LOG2FC_THRESHOLD
  ),
  "",
  "7. Signature exclusion",
  paste0("  Official signature genes: ", length(signature_keys)),
  paste0(
    "  Signature-gene rows in tested results: ",
    n_signature_rows_in_results
  ),
  paste0(
    "  Unique signature genes tested: ",
    n_signature_genes_tested
  ),
  paste0(
    "  Signature rows excluded from no-signature table: ",
    n_signature_rows_in_results
  ),
  paste0(
    "  Rows using Ensembl ID as gene-name fallback: ",
    n_missing_gene_names
  ),
  "",
  "8. Results",
  paste0("  Genes tested: ", nrow(full_results)),
  paste0(
    "  Genes tested after signature exclusion: ",
    nrow(no_signature_results)
  ),
  paste0(
    "  Genes FDR < 0.05, full: ",
    n_significant_fdr_full
  ),
  paste0(
    "  Genes FDR < 0.05, no signature: ",
    n_significant_fdr_no_signature
  ),
  paste0(
    "  DEG no signature, FDR < 0.05 and |shrunken log2FC| > 1: ",
    nrow(significant_no_signature)
  ),
  paste0("  UP in HIGH_Q4: ", n_up),
  paste0("  DOWN in HIGH_Q4: ", n_down),
  "",
  "9. Summary table",
  paste(
    capture.output(
      print(summary_df, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "10. Outputs",
  paste0("  Full results: ", output_paths$full_results),
  paste0("  No-signature results: ", output_paths$no_signature),
  paste0(
    "  Significant no-signature DEGs: ",
    output_paths$significant_no_signature
  ),
  paste0("  Sample metadata: ", output_paths$sample_metadata),
  paste0("  Summary: ", output_paths$summary),
  paste0("  Metadata: ", output_paths$metadata),
  paste0("  DESeq2 object: ", output_paths$dds),
  paste0("  VST object: ", output_paths$vst),
  paste0("  QC: ", output_paths$qc),
  "",
  "11. Warnings",
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

required_outputs <- unlist(output_paths, use.names = FALSE)

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
msg("Significant DEGs:", output_paths$significant_no_signature)
msg("DESeq2 object:", output_paths$dds)
msg("VST object:", output_paths$vst)
msg("QC:", output_paths$qc)
msg("Final status:", final_status)
