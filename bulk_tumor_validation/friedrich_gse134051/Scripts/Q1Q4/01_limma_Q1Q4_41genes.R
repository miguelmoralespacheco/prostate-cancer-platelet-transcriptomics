#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Q1/Q4 limma differential expression
#
# This script compares gene expression between HIGH_Q4 and LOW_Q1 Friedrich/GSE134051 primary
# tumors using the official 41-gene platelet-associated transcriptional score.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

COHORT_ID <- "Friedrich_GSE134051"
COHORT_LABEL <- "Friedrich/GSE134051"
COHORT_SHORT <- "Friedrich"
SOURCE_DATASET <- "GSE134051"
SOURCE_STUDY <- "Friedrich et al."
SOURCE_ASSAY <- "gex.logq"
EXPECTED_TOTAL_N <- 255L
EXPECTED_PRIMARY_N <- 164L
EXPECTED_EXPRESSION_GENES <- 23097L

message("[1/6] Loading inputs")

.repo_override <- Sys.getenv("PLATELET_REPO_ROOT", unset = "")
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (nzchar(.repo_override)) {
  .config_path <- file.path(
    normalizePath(.repo_override, winslash = "/", mustWork = TRUE),
    "bulk_tumor_validation", "friedrich_gse134051", "Scripts", "00_config_paths.R"
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

input_expression <- file.path(
  INPUT_DIR, "LocalLarge",
  "Friedrich_GSE134051_expression_logq.rds"
)

input_q1q4 <- file.path(
  GENERATED_RESULTS_DIR, "Score",
  "Friedrich_GSE134051_primary_Q1Q4_groups_41genes.csv"
)

input_signature <- file.path(
  RESOURCE_DIR,
  "platelet_associated_transcriptional_signature.tsv"
)

tables_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4",
  "Tables"
)

limma_dir <- file.path(
  tables_dir,
  "LIMMA"
)

rds_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4",
  "Objects"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4",
  "Logs"
)

output_full <- file.path(
  limma_dir,
  "LIMMA_Q1Q4_41genes_full_results.csv"
)

output_no_signature <- file.path(
  limma_dir,
  "LIMMA_Q1Q4_41genes_no_signature_genes.csv"
)

output_significant_no_signature <- file.path(
  limma_dir,
  "LIMMA_Q1Q4_41genes_no_signature_DEG_FDR005_abslogFC1.csv"
)

output_summary <- file.path(
  limma_dir,
  "LIMMA_Q1Q4_41genes_summary.csv"
)

output_metadata <- file.path(
  limma_dir,
  "LIMMA_Q1Q4_41genes_metadata.json"
)

output_fit <- file.path(
  rds_dir,
  "LIMMA_Q1Q4_41genes_fit.rds"
)

output_qc <- file.path(
  logs_dir,
  "LIMMA_Q1Q4_41genes_QC.txt"
)

dir.create(
  tables_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  limma_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  rds_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  logs_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

warnings_vec <- character(0)

fail_clear <- function(msg) {
  stop(msg, call. = FALSE)
}

add_warning <- function(msg) {
  warnings_vec <<- c(warnings_vec, msg)
  warning(msg, call. = FALSE)
}

write_csv_base <- function(x, path) {
  utils::write.csv(
    x,
    path,
    row.names = FALSE,
    quote = TRUE
  )
}

read_csv_base <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_metadata_json <- function(metadata, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    fail_clear(
      "Package 'jsonlite' is required to write metadata JSON."
    )
  }

  jsonlite::write_json(
    metadata,
    path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
}

detect_gene_column <- function(signature_df) {
  candidates <- c(
    "gene_name",
    "gene_symbol",
    "symbol",
    "gene",
    "external_gene_name"
  )

  lower_names <- tolower(names(signature_df))
  idx <- match(candidates, lower_names)
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0) {
    fail_clear(
      paste0(
        "Could not detect gene symbol column in signature. Accepted names: ",
        paste(candidates, collapse = ", ")
      )
    )
  }

  names(signature_df)[idx[1]]
}

extract_signature_genes <- function(signature_df) {
  gene_col <- detect_gene_column(signature_df)

  genes <- trimws(
    as.character(signature_df[[gene_col]])
  )

  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]

  duplicated_genes <- unique(
    genes[duplicated(genes)]
  )

  if (length(duplicated_genes) > 0) {
    fail_clear(
      paste0(
        "Duplicate genes detected in official signature: ",
        paste(duplicated_genes, collapse = ", ")
      )
    )
  }

  if (length(genes) != 41) {
    fail_clear(
      paste0(
        "Official signature must contain exactly 41 unique genes. Observed: ",
        length(genes)
      )
    )
  }

  genes
}

is_default_rownames <- function(x) {
  rn <- rownames(x)

  is.null(rn) ||
    identical(
      rn,
      as.character(seq_len(nrow(x)))
    )
}

coerce_expression_matrix <- function(expr_obj) {
  expr_class <- paste(
    class(expr_obj),
    collapse = ","
  )

  if (is.matrix(expr_obj)) {
    expr_mat <- expr_obj

  } else if (is.data.frame(expr_obj)) {
    expr_df <- expr_obj

    if (
      is_default_rownames(expr_df) &&
      ncol(expr_df) >= 2
    ) {
      first_col <- expr_df[[1]]

      first_col_chr <- trimws(
        as.character(first_col)
      )

      first_col_nonempty <- !any(
        is.na(first_col_chr) |
          first_col_chr == ""
      )

      first_col_numeric <- suppressWarnings(
        all(
          !is.na(
            as.numeric(first_col_chr)
          )
        )
      )

      first_col_unique <- !anyDuplicated(
        first_col_chr
      )

      if (
        first_col_nonempty &&
        !first_col_numeric &&
        first_col_unique
      ) {
        rownames(expr_df) <- first_col_chr

        expr_df <- expr_df[
          ,
          -1,
          drop = FALSE
        ]

        add_warning(
          paste0(
            "Expression data.frame had default rownames; ",
            "first non-numeric unique column was used as gene rownames."
          )
        )
      }
    }

    expr_mat <- as.matrix(expr_df)

  } else {
    fail_clear(
      paste0(
        "Expression object must be a matrix or data.frame. Observed class: ",
        expr_class
      )
    )
  }

  storage.mode(expr_mat) <- "numeric"
  attr(expr_mat, "input_class") <- expr_class

  expr_mat
}

validate_friedrich_expression_identity <- function(expr_mat) {
  expected_dim <- c(EXPECTED_EXPRESSION_GENES, EXPECTED_TOTAL_N)
  forbidden_id_pattern <- paste0("^ICGC", "_PCA")

  if (!identical(dim(expr_mat), expected_dim)) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: expected ",
        EXPECTED_EXPRESSION_GENES,
        " genes x ",
        EXPECTED_TOTAL_N,
        " samples for ",
        COHORT_LABEL,
        "; observed ",
        nrow(expr_mat),
        " x ",
        ncol(expr_mat),
        "."
      )
    )
  }

  sample_ids <- colnames(expr_mat)
  if (
    any(!grepl("^GSM[0-9]+$", sample_ids)) ||
    any(grepl(forbidden_id_pattern, sample_ids, ignore.case = TRUE)) ||
    anyDuplicated(sample_ids)
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: ",
        COHORT_LABEL,
        " expression must contain unique GSM sample identifiers only."
      )
    )
  }

  if (!is.numeric(expr_mat)) {
    fail_clear(
      paste0(
        "Possible incorrect input: ",
        COHORT_LABEL,
        " expression must be numeric continuous log-expression."
      )
    )
  }

  finite_values <- expr_mat[is.finite(expr_mat)]
  if (
    length(finite_values) == 0 ||
    all(abs(finite_values - round(finite_values)) < 1e-8)
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input: ",
        SOURCE_ASSAY,
        " must contain continuous log-expression values, not integer counts."
      )
    )
  }

  invisible(TRUE)
}

validate_expression_orientation <- function(
    expr_mat,
    sample_names,
    signature_genes
) {
  if (!is.matrix(expr_mat)) {
    fail_clear(
      "Expression object could not be converted to matrix."
    )
  }

  if (
    nrow(expr_mat) == 0 ||
    ncol(expr_mat) == 0
  ) {
    fail_clear(
      "Expression matrix is empty."
    )
  }

  if (
    is.null(rownames(expr_mat)) ||
    any(is.na(rownames(expr_mat))) ||
    any(rownames(expr_mat) == "")
  ) {
    fail_clear(
      "Expression matrix must have non-empty gene rownames."
    )
  }

  if (
    is.null(colnames(expr_mat)) ||
    any(is.na(colnames(expr_mat))) ||
    any(colnames(expr_mat) == "")
  ) {
    fail_clear(
      "Expression matrix must have non-empty sample colnames."
    )
  }

  row_gene_hits <- sum(
    signature_genes %in% rownames(expr_mat)
  )

  col_gene_hits <- sum(
    signature_genes %in% colnames(expr_mat)
  )

  col_sample_hits <- sum(
    sample_names %in% colnames(expr_mat)
  )

  row_sample_hits <- sum(
    sample_names %in% rownames(expr_mat)
  )

  if (
    row_sample_hits > col_sample_hits ||
    col_gene_hits > row_gene_hits
  ) {
    fail_clear(
      paste0(
        "Expression matrix appears transposed or malformed. ",
        "sample hits row/col = ",
        row_sample_hits,
        "/",
        col_sample_hits,
        "; signature gene hits row/col = ",
        row_gene_hits,
        "/",
        col_gene_hits,
        ". This script does not transpose automatically."
      )
    )
  }

  if (row_gene_hits == 0) {
    fail_clear(
      paste0(
        "No official signature genes found in expression rownames; ",
        "matrix orientation or gene symbols may be wrong."
      )
    )
  }

  invisible(TRUE)
}

write_qc_report <- function(path, qc) {
  sink(path)

  cat("Friedrich/GSE134051 Q1/Q4 limma differential expression QC\n")
  cat("==================================================\n\n")

  cat("1. Inputs\n")
  cat(
    "  Expression input: ",
    qc$inputs$expression,
    "\n",
    sep = ""
  )
  cat(
    "  Q1/Q4 input: ",
    qc$inputs$q1q4,
    "\n",
    sep = ""
  )
  cat(
    "  Signature input: ",
    qc$inputs$signature,
    "\n\n",
    sep = ""
  )

  cat("2. Sample universe\n")
  cat(
    "  Total Q1/Q4 samples: ",
    qc$samples$n_total,
    "\n",
    sep = ""
  )
  cat(
    "  LOW_Q1: ",
    qc$samples$n_low,
    "\n",
    sep = ""
  )
  cat(
    "  HIGH_Q4: ",
    qc$samples$n_high,
    "\n",
    sep = ""
  )
  cat(
    "  MID_Q2Q3 rows excluded: ",
    qc$samples$n_mid_excluded,
    "\n",
    sep = ""
  )
  cat(
    "  Other rows excluded: ",
    qc$samples$n_other_excluded,
    "\n\n",
    sep = ""
  )

  cat("3. Expression matrix\n")
  cat(
    "  Input class: ",
    qc$expression$input_class,
    "\n",
    sep = ""
  )
  cat(
    "  Original dimensions: ",
    qc$expression$n_genes_original,
    " genes x ",
    qc$expression$n_samples_original,
    " samples\n",
    sep = ""
  )
  cat(
    "  Q1/Q4 dimensions before filtering: ",
    qc$expression$n_genes_q1q4,
    " genes x ",
    qc$expression$n_samples_q1q4,
    " samples\n\n",
    sep = ""
  )

  cat("4. Filtering\n")
  cat(
    "  Genes removed all/non-finite: ",
    qc$filter$n_removed_nonfinite,
    "\n",
    sep = ""
  )
  cat(
    "  Genes removed zero/non-finite variance: ",
    qc$filter$n_removed_zero_variance,
    "\n",
    sep = ""
  )
  cat(
    "  Genes tested: ",
    qc$filter$n_genes_tested,
    "\n\n",
    sep = ""
  )

  cat("5. Limma model\n")
  cat("  Method: limma\n")
  cat("  Design: ~ 0 + platelet_score_group\n")
  cat("  Contrast: HIGH_Q4 - LOW_Q1\n")
  cat("  Direction positive: higher in HIGH_Q4\n")
  cat("  Direction negative: higher in LOW_Q1\n")
  cat(
    "  Fit completed: ",
    qc$limma$fit_completed,
    "\n\n",
    sep = ""
  )

  cat("6. Signature gene handling\n")
  cat(
    "  Signature genes total: ",
    qc$signature$n_total,
    "\n",
    sep = ""
  )
  cat(
    "  Signature genes present in results: ",
    qc$signature$n_present_in_results,
    "\n",
    sep = ""
  )
  cat(
    "  Signature genes removed from no-signature table: ",
    qc$signature$n_removed_from_no_signature,
    "\n",
    sep = ""
  )
  cat(
    "  No-signature table contains signature genes: ",
    qc$signature$no_signature_contains_signature,
    "\n\n",
    sep = ""
  )

  cat("7. Differential expression summary\n")
  cat(
    "  Full FDR<0.05: ",
    qc$de$n_full_fdr005,
    "\n",
    sep = ""
  )
  cat(
    "  Full FDR<0.05 & abs(logFC)>1: ",
    qc$de$n_full_fdr005_abslogfc1,
    "\n",
    sep = ""
  )
  cat(
    "  Full up FDR<0.05 & abs(logFC)>1: ",
    qc$de$n_full_up_fdr005_abslogfc1,
    "\n",
    sep = ""
  )
  cat(
    "  Full down FDR<0.05 & abs(logFC)>1: ",
    qc$de$n_full_down_fdr005_abslogfc1,
    "\n",
    sep = ""
  )
  cat(
    "  No-signature genes tested: ",
    qc$de$n_no_signature_genes_tested,
    "\n",
    sep = ""
  )
  cat(
    "  No-signature FDR<0.05: ",
    qc$de$n_no_signature_fdr005,
    "\n",
    sep = ""
  )
  cat(
    "  No-signature FDR<0.05 & abs(logFC)>1: ",
    qc$de$n_no_signature_fdr005_abslogfc1,
    "\n",
    sep = ""
  )
  cat(
    "  No-signature up FDR<0.05 & abs(logFC)>1: ",
    qc$de$n_no_signature_up_fdr005_abslogfc1,
    "\n",
    sep = ""
  )
  cat(
    "  No-signature down FDR<0.05 & abs(logFC)>1: ",
    qc$de$n_no_signature_down_fdr005_abslogfc1,
    "\n\n",
    sep = ""
  )

  cat("8. Warnings\n")

  if (length(qc$warnings) == 0) {
    cat("  None\n\n")
  } else {
    for (w in qc$warnings) {
      cat(
        "  - ",
        w,
        "\n",
        sep = ""
      )
    }

    cat("\n")
  }

  cat(
    "9. Final status: ",
    qc$qc_status,
    "\n",
    sep = ""
  )

  sink()
}

make_empty_qc <- function() {
  list(
    inputs = list(
      expression = input_expression,
      q1q4 = input_q1q4,
      signature = input_signature
    ),
    samples = list(
      n_total = NA_integer_,
      n_low = NA_integer_,
      n_high = NA_integer_,
      n_mid_excluded = NA_integer_,
      n_other_excluded = NA_integer_
    ),
    expression = list(
      input_class = NA_character_,
      n_genes_original = NA_integer_,
      n_samples_original = NA_integer_,
      n_genes_q1q4 = NA_integer_,
      n_samples_q1q4 = NA_integer_
    ),
    filter = list(
      n_removed_nonfinite = NA_integer_,
      n_removed_zero_variance = NA_integer_,
      n_genes_tested = NA_integer_
    ),
    limma = list(
      fit_completed = FALSE
    ),
    signature = list(
      n_total = NA_integer_,
      n_present_in_results = NA_integer_,
      n_removed_from_no_signature = NA_integer_,
      no_signature_contains_signature = NA
    ),
    de = list(
      n_full_fdr005 = NA_integer_,
      n_full_fdr005_abslogfc1 = NA_integer_,
      n_full_up_fdr005_abslogfc1 = NA_integer_,
      n_full_down_fdr005_abslogfc1 = NA_integer_,
      n_no_signature_genes_tested = NA_integer_,
      n_no_signature_fdr005 = NA_integer_,
      n_no_signature_fdr005_abslogfc1 = NA_integer_,
      n_no_signature_up_fdr005_abslogfc1 = NA_integer_,
      n_no_signature_down_fdr005_abslogfc1 = NA_integer_
    ),
    warnings = warnings_vec,
    qc_status = "FAIL"
  )
}

fail_with_qc <- function(
    msg,
    qc = make_empty_qc()
) {
  qc$warnings <- warnings_vec
  qc$qc_status <- "FAIL"

  write_qc_report(
    output_qc,
    qc
  )

  fail_clear(msg)
}

derive_qc_status <- function(
    n_low,
    n_high,
    n_genes_tested,
    warnings_vec,
    no_signature_contains_signature,
    fit_completed,
    full_written,
    no_signature_written
) {
  if (
    n_low < 30 ||
    n_high < 30
  ) {
    return("FAIL")
  }

  if (n_genes_tested < 5000) {
    return("FAIL")
  }

  if (isTRUE(no_signature_contains_signature)) {
    return("FAIL")
  }

  if (
    !isTRUE(fit_completed) ||
    !isTRUE(full_written) ||
    !isTRUE(no_signature_written)
  ) {
    return("FAIL")
  }

  if (
    n_low != 41 ||
    n_high != 41
  ) {
    return("WARNING")
  }

  if (n_genes_tested <= 10000) {
    return("WARNING")
  }

  if (length(warnings_vec) > 0) {
    return("WARNING")
  }

  "PASS"
}

summarize_de <- function(res) {
  sig <- res$significant_FDR005
  sig_abs <- res$significant_FDR005_abslogFC1

  list(
    n_fdr005 = sum(
      sig,
      na.rm = TRUE
    ),
    n_fdr005_abslogfc1 = sum(
      sig_abs,
      na.rm = TRUE
    ),
    n_up_fdr005_abslogfc1 = sum(
      sig_abs &
        res$logFC > 0,
      na.rm = TRUE
    ),
    n_down_fdr005_abslogfc1 = sum(
      sig_abs &
        res$logFC < 0,
      na.rm = TRUE
    )
  )
}

required_packages <- c(
  "limma",
  "jsonlite"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    fail_with_qc(
      paste0(
        "Required package not installed: ",
        pkg
      )
    )
  }
}

for (
  fp in c(
    input_expression,
    input_q1q4,
    input_signature
  )
) {
  if (!file.exists(fp)) {
    fail_with_qc(
      paste0(
        "Missing required input: ",
        fp
      )
    )
  }
}

q1q4_raw <- read_csv_base(input_q1q4)
signature_df <- read_canonical_platelet_signature()
signature_genes <- extract_signature_genes(signature_df)
expr_obj <- readRDS(input_expression)
expr_mat <- coerce_expression_matrix(expr_obj)
validate_friedrich_expression_identity(expr_mat)

if (!"sample_name" %in% names(q1q4_raw)) {
  fail_with_qc(
    "Q1/Q4 input must contain column: sample_name"
  )
}

if (!"platelet_score_group" %in% names(q1q4_raw)) {
  fail_with_qc(
    "Q1/Q4 input must contain column: platelet_score_group"
  )
}

if (anyDuplicated(q1q4_raw$sample_name)) {
  dup_samples <- unique(
    q1q4_raw$sample_name[
      duplicated(q1q4_raw$sample_name)
    ]
  )

  fail_with_qc(
    paste0(
      "Duplicate sample_name values in Q1/Q4 input: ",
      paste(dup_samples, collapse = ", ")
    )
  )
}

message("[2/6] Validating Q1/Q4 samples")

q1q4_raw$sample_name <- as.character(
  q1q4_raw$sample_name
)

forbidden_id_pattern <- paste0("^ICGC", "_PCA")
if (
  any(!grepl("^GSM[0-9]+$", q1q4_raw$sample_name)) ||
  any(grepl(forbidden_id_pattern, q1q4_raw$sample_name, ignore.case = TRUE))
) {
  fail_with_qc(
    paste0(
      "Possible incorrect input or cohort contamination: Q1/Q4 sample_name ",
      "values must be GSM identifiers from ",
      COHORT_LABEL,
      "."
    )
  )
}

q1q4_raw$platelet_score_group <- as.character(
  q1q4_raw$platelet_score_group
)

n_mid_excluded <- sum(
  q1q4_raw$platelet_score_group == "MID_Q2Q3",
  na.rm = TRUE
)

if (n_mid_excluded > 0) {
  add_warning(
    paste0(
      "Excluding MID_Q2Q3 rows from Q1/Q4 limma analysis: ",
      n_mid_excluded
    )
  )
}

known_groups <- c(
  "LOW_Q1",
  "HIGH_Q4",
  "MID_Q2Q3"
)

n_other_excluded <- sum(
  !q1q4_raw$platelet_score_group %in% known_groups
)

if (n_other_excluded > 0) {
  add_warning(
    paste0(
      "Excluding rows with non-Q1/Q4 platelet_score_group values: ",
      n_other_excluded
    )
  )
}

q1q4 <- q1q4_raw[
  q1q4_raw$platelet_score_group %in%
    c(
      "LOW_Q1",
      "HIGH_Q4"
    ),
  ,
  drop = FALSE
]

q1q4$platelet_score_group <- factor(
  q1q4$platelet_score_group,
  levels = c(
    "LOW_Q1",
    "HIGH_Q4"
  )
)

n_low <- sum(
  q1q4$platelet_score_group == "LOW_Q1"
)

n_high <- sum(
  q1q4$platelet_score_group == "HIGH_Q4"
)

n_total <- nrow(q1q4)

sample_qc <- list(
  n_total = n_total,
  n_low = n_low,
  n_high = n_high,
  n_mid_excluded = n_mid_excluded,
  n_other_excluded = n_other_excluded
)

if (n_total == 0) {
  qc <- make_empty_qc()
  qc$samples <- sample_qc

  fail_with_qc(
    "No LOW_Q1/HIGH_Q4 samples available after filtering.",
    qc
  )
}

if (
  n_low < 30 ||
  n_high < 30
) {
  qc <- make_empty_qc()
  qc$samples <- sample_qc

  fail_with_qc(
    paste0(
      "Insufficient group size for limma: LOW_Q1=",
      n_low,
      ", HIGH_Q4=",
      n_high
    ),
    qc
  )
}

if (
  n_low != 41 ||
  n_high != 41 ||
  n_total != 82
) {
  add_warning(
    paste0(
      "Expected LOW_Q1=41, HIGH_Q4=41, total=82; observed LOW_Q1=",
      n_low,
      ", HIGH_Q4=",
      n_high,
      ", total=",
      n_total
    )
  )
}

message("[3/6] Preparing expression matrix")

tryCatch(
  validate_expression_orientation(
    expr_mat,
    q1q4$sample_name,
    signature_genes
  ),
  error = function(e) {
    qc <- make_empty_qc()
    qc$samples <- sample_qc
    qc$signature$n_total <- length(signature_genes)
    qc$expression$input_class <- attr(expr_mat, "input_class")
    qc$expression$n_genes_original <- nrow(expr_mat)
    qc$expression$n_samples_original <- ncol(expr_mat)

    fail_with_qc(
      conditionMessage(e),
      qc
    )
  }
)

if (anyDuplicated(rownames(expr_mat))) {
  dup_genes <- unique(
    rownames(expr_mat)[
      duplicated(rownames(expr_mat))
    ]
  )

  qc <- make_empty_qc()
  qc$samples <- sample_qc
  qc$signature$n_total <- length(signature_genes)

  fail_with_qc(
    paste0(
      "Duplicate gene rownames in expression matrix: ",
      paste(
        utils::head(dup_genes, 25),
        collapse = ", "
      ),
      if (length(dup_genes) > 25) {
        " ..."
      } else {
        ""
      }
    ),
    qc
  )
}

missing_samples <- setdiff(
  q1q4$sample_name,
  colnames(expr_mat)
)

if (length(missing_samples) > 0) {
  qc <- make_empty_qc()
  qc$samples <- sample_qc
  qc$signature$n_total <- length(signature_genes)
  qc$expression$input_class <- attr(expr_mat, "input_class")
  qc$expression$n_genes_original <- nrow(expr_mat)
  qc$expression$n_samples_original <- ncol(expr_mat)

  fail_with_qc(
    paste0(
      "Q1/Q4 samples missing from expression matrix: ",
      paste(missing_samples, collapse = ", ")
    ),
    qc
  )
}

expr_q <- expr_mat[
  ,
  q1q4$sample_name,
  drop = FALSE
]

if (
  !identical(
    colnames(expr_q),
    q1q4$sample_name
  )
) {
  fail_with_qc(
    "Expression columns could not be aligned with Q1/Q4 metadata."
  )
}

finite_counts <- rowSums(
  is.finite(expr_q)
)

keep_has_finite <- finite_counts > 0

n_removed_nonfinite <- sum(
  !keep_has_finite
)

row_variance <- vapply(
  seq_len(nrow(expr_q)),
  function(i) {
    vals <- expr_q[i, ]
    vals <- vals[is.finite(vals)]

    if (length(vals) < 2) {
      return(NA_real_)
    }

    stats::var(vals)
  },
  numeric(1)
)

keep_finite_variance <- is.finite(row_variance) &
  row_variance > 0

n_removed_zero_variance <- sum(
  keep_has_finite &
    !keep_finite_variance
)

keep_genes <- keep_has_finite &
  keep_finite_variance

expr_test <- expr_q[
  keep_genes,
  ,
  drop = FALSE
]

if (nrow(expr_test) < 5000) {
  qc <- make_empty_qc()
  qc$samples <- sample_qc

  qc$expression <- list(
    input_class = attr(expr_mat, "input_class"),
    n_genes_original = nrow(expr_mat),
    n_samples_original = ncol(expr_mat),
    n_genes_q1q4 = nrow(expr_q),
    n_samples_q1q4 = ncol(expr_q)
  )

  qc$filter <- list(
    n_removed_nonfinite = n_removed_nonfinite,
    n_removed_zero_variance = n_removed_zero_variance,
    n_genes_tested = nrow(expr_test)
  )

  qc$signature$n_total <- length(signature_genes)

  fail_with_qc(
    paste0(
      "Too few genes remain after filtering: ",
      nrow(expr_test)
    ),
    qc
  )
}

if (nrow(expr_test) <= 10000) {
  add_warning(
    paste0(
      "Genes tested between 5000 and 10000: ",
      nrow(expr_test)
    )
  )
}

expr_test[
  !is.finite(expr_test)
] <- NA_real_

message("[4/6] Running limma")

design <- stats::model.matrix(
  ~ 0 + platelet_score_group,
  data = q1q4
)

colnames(design) <- sub(
  "^platelet_score_group",
  "",
  colnames(design)
)

if (
  !identical(
    colnames(design),
    c(
      "LOW_Q1",
      "HIGH_Q4"
    )
  )
) {
  fail_with_qc(
    paste0(
      "Unexpected design columns: ",
      paste(
        colnames(design),
        collapse = ", "
      )
    )
  )
}

fit_completed <- FALSE

fit <- tryCatch(
  {
    fit0 <- limma::lmFit(
      expr_test,
      design
    )

    contrast_matrix <- limma::makeContrasts(
      HIGH_Q4 - LOW_Q1,
      levels = design
    )

    fit1 <- limma::contrasts.fit(
      fit0,
      contrast_matrix
    )

    fit2 <- limma::eBayes(fit1)

    fit_completed <<- TRUE

    fit2
  },
  error = function(e) {
    fail_with_qc(
      paste0(
        "limma failed: ",
        conditionMessage(e)
      )
    )
  }
)

tt <- limma::topTable(
  fit,
  number = Inf,
  sort.by = "none"
)

tt$gene_name <- rownames(tt)

required_limma_cols <- c(
  "gene_name",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)

missing_limma_cols <- setdiff(
  required_limma_cols,
  names(tt)
)

if (length(missing_limma_cols) > 0) {
  fail_with_qc(
    paste0(
      "limma output missing required columns: ",
      paste(
        missing_limma_cols,
        collapse = ", "
      )
    )
  )
}

full_results <- data.frame(
  gene_name = tt$gene_name,
  logFC = tt$logFC,
  AveExpr = tt$AveExpr,
  t_stat = tt$t,
  P.Value = tt$P.Value,
  adj.P.Val = tt$adj.P.Val,
  B = tt$B,
  stringsAsFactors = FALSE
)

full_results$direction <- ifelse(
  full_results$logFC > 0,
  "HIGH_Q4",
  ifelse(
    full_results$logFC < 0,
    "LOW_Q1",
    "ZERO"
  )
)

full_results$significant_FDR005 <-
  full_results$adj.P.Val < 0.05

full_results$significant_FDR005_abslogFC1 <-
  full_results$adj.P.Val < 0.05 &
  abs(full_results$logFC) > 1

full_results$is_signature_gene <-
  full_results$gene_name %in% signature_genes

full_results <- full_results[
  order(full_results$gene_name),
  ,
  drop = FALSE
]

no_signature_results <- full_results[
  !full_results$is_signature_gene,
  ,
  drop = FALSE
]

significant_no_signature_results <- no_signature_results[
  !is.na(
    no_signature_results$significant_FDR005_abslogFC1
  ) &
    no_signature_results$significant_FDR005_abslogFC1,
  ,
  drop = FALSE
]

significant_no_signature_results <-
  significant_no_signature_results[
    order(
      significant_no_signature_results$adj.P.Val,
      -abs(significant_no_signature_results$logFC),
      significant_no_signature_results$gene_name
    ),
    ,
    drop = FALSE
  ]

no_signature_contains_signature <- any(
  no_signature_results$gene_name %in%
    signature_genes
)

significant_no_signature_contains_signature <- any(
  significant_no_signature_results$gene_name %in%
    signature_genes
)

if (no_signature_contains_signature) {
  fail_with_qc(
    "No-signature results still contain official signature genes."
  )
}

if (significant_no_signature_contains_signature) {
  fail_with_qc(
    paste0(
      "Significant no-signature DEG table still contains ",
      "official signature genes."
    )
  )
}

n_signature_present <- sum(
  full_results$is_signature_gene
)

n_signature_removed <- nrow(full_results) -
  nrow(no_signature_results)

full_de <- summarize_de(full_results)

no_sig_de <- summarize_de(
  no_signature_results
)

prewrite_full_ready <- nrow(full_results) > 0

prewrite_no_signature_ready <-
  nrow(no_signature_results) > 0

prewrite_significant_no_signature_ready <-
  is.data.frame(significant_no_signature_results)

qc_status <- derive_qc_status(
  n_low = n_low,
  n_high = n_high,
  n_genes_tested = nrow(full_results),
  warnings_vec = warnings_vec,
  no_signature_contains_signature =
    no_signature_contains_signature,
  fit_completed = fit_completed,
  full_written = prewrite_full_ready,
  no_signature_written =
    prewrite_no_signature_ready
)

summary_df <- data.frame(
  metric = c(
    "cohort",
    "analysis",
    "analytic_universe",
    "expression_input",
    "q1q4_input",
    "signature_input",
    "n_samples_total",
    "n_LOW_Q1",
    "n_HIGH_Q4",
    "n_expression_genes_original",
    "n_expression_genes_q1q4",
    "n_genes_removed_nonfinite",
    "n_genes_removed_zero_variance",
    "n_genes_tested",
    "n_signature_genes_total",
    "n_signature_genes_present_in_results",
    "n_signature_genes_removed_from_no_signature",
    "n_full_FDR005",
    "n_full_FDR005_abslogFC1",
    "n_full_up_FDR005_abslogFC1",
    "n_full_down_FDR005_abslogFC1",
    "n_no_signature_genes_tested",
    "n_no_signature_FDR005",
    "n_no_signature_FDR005_abslogFC1",
    "n_no_signature_up_FDR005_abslogFC1",
    "n_no_signature_down_FDR005_abslogFC1",
    "contrast",
    "direction_positive",
    "direction_negative",
    "qc_status"
  ),
  value = as.character(
    c(
      "Friedrich/GSE134051",
      "Q1Q4_limma_41genes",
      "primary tumors only",
      input_expression,
      input_q1q4,
      input_signature,
      n_total,
      n_low,
      n_high,
      nrow(expr_mat),
      nrow(expr_q),
      n_removed_nonfinite,
      n_removed_zero_variance,
      nrow(full_results),
      length(signature_genes),
      n_signature_present,
      n_signature_removed,
      full_de$n_fdr005,
      full_de$n_fdr005_abslogfc1,
      full_de$n_up_fdr005_abslogfc1,
      full_de$n_down_fdr005_abslogfc1,
      nrow(no_signature_results),
      no_sig_de$n_fdr005,
      no_sig_de$n_fdr005_abslogfc1,
      no_sig_de$n_up_fdr005_abslogfc1,
      no_sig_de$n_down_fdr005_abslogfc1,
      "HIGH_Q4 - LOW_Q1",
      "higher in HIGH_Q4",
      "higher in LOW_Q1",
      qc_status
    )
  ),
  stringsAsFactors = FALSE
)

metadata <- list(
  date_time = as.character(Sys.time()),
  project_dir = project_dir,
  inputs = list(
    expression = input_expression,
    q1q4_groups = input_q1q4,
    signature = input_signature
  ),
  outputs = list(
    full_results = output_full,
    no_signature_results = output_no_signature,
    significant_no_signature_results =
      output_significant_no_signature,
    summary = output_summary,
    metadata = output_metadata,
    fit = output_fit,
    qc = output_qc
  ),
  method = "limma",
  contrast = "HIGH_Q4 - LOW_Q1",
  direction_positive = "higher in HIGH_Q4",
  direction_negative = "higher in LOW_Q1",
  deg_thresholds = list(
    adjusted_p_value = 0.05,
    absolute_log2_fold_change = 1,
    signature_genes_excluded = TRUE
  ),
  sample_counts = list(
    n_samples_total = n_total,
    n_LOW_Q1 = n_low,
    n_HIGH_Q4 = n_high
  ),
  expression_dimensions = list(
    original_genes = nrow(expr_mat),
    original_samples = ncol(expr_mat),
    q1q4_genes_before_filter = nrow(expr_q),
    q1q4_samples = ncol(expr_q),
    genes_tested = nrow(full_results)
  ),
  differential_expression = list(
    n_no_signature_FDR005_abslogFC1 =
      nrow(significant_no_signature_results),
    n_no_signature_up_FDR005_abslogFC1 =
      no_sig_de$n_up_fdr005_abslogfc1,
    n_no_signature_down_FDR005_abslogFC1 =
      no_sig_de$n_down_fdr005_abslogfc1
  ),
  filter_details = list(
    n_genes_removed_nonfinite =
      n_removed_nonfinite,
    n_genes_removed_zero_variance =
      n_removed_zero_variance
  ),
  signature_genes = signature_genes,
  qc_status = qc_status
)

qc <- list(
  inputs = list(
    expression = input_expression,
    q1q4 = input_q1q4,
    signature = input_signature
  ),
  samples = sample_qc,
  expression = list(
    input_class = attr(expr_mat, "input_class"),
    n_genes_original = nrow(expr_mat),
    n_samples_original = ncol(expr_mat),
    n_genes_q1q4 = nrow(expr_q),
    n_samples_q1q4 = ncol(expr_q)
  ),
  filter = list(
    n_removed_nonfinite =
      n_removed_nonfinite,
    n_removed_zero_variance =
      n_removed_zero_variance,
    n_genes_tested =
      nrow(full_results)
  ),
  limma = list(
    fit_completed = fit_completed,
    full_results_ready =
      prewrite_full_ready,
    no_signature_results_ready =
      prewrite_no_signature_ready,
    significant_no_signature_results_ready =
      prewrite_significant_no_signature_ready
  ),
  signature = list(
    n_total = length(signature_genes),
    n_present_in_results =
      n_signature_present,
    n_removed_from_no_signature =
      n_signature_removed,
    no_signature_contains_signature =
      no_signature_contains_signature
  ),
  de = list(
    n_full_fdr005 =
      full_de$n_fdr005,
    n_full_fdr005_abslogfc1 =
      full_de$n_fdr005_abslogfc1,
    n_full_up_fdr005_abslogfc1 =
      full_de$n_up_fdr005_abslogfc1,
    n_full_down_fdr005_abslogfc1 =
      full_de$n_down_fdr005_abslogfc1,
    n_no_signature_genes_tested =
      nrow(no_signature_results),
    n_no_signature_fdr005 =
      no_sig_de$n_fdr005,
    n_no_signature_fdr005_abslogfc1 =
      no_sig_de$n_fdr005_abslogfc1,
    n_no_signature_up_fdr005_abslogfc1 =
      no_sig_de$n_up_fdr005_abslogfc1,
    n_no_signature_down_fdr005_abslogfc1 =
      no_sig_de$n_down_fdr005_abslogfc1
  ),
  warnings = warnings_vec,
  qc_status = qc_status
)

message("[5/6] Writing outputs")

write_csv_base(
  full_results,
  output_full
)

write_csv_base(
  no_signature_results,
  output_no_signature
)

write_csv_base(
  significant_no_signature_results,
  output_significant_no_signature
)

write_csv_base(
  summary_df,
  output_summary
)

write_metadata_json(
  metadata,
  output_metadata
)

saveRDS(
  fit,
  output_fit
)

write_qc_report(
  output_qc,
  qc
)

if (qc_status == "FAIL") {
  fail_clear(
    "Final status: FAIL"
  )
}

message(
  "[6/6] Final status: ",
  qc_status
)
