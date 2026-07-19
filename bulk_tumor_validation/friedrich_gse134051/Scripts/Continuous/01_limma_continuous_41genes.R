#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Continuous limma association analysis
#
# Model: expression ~ Score_z
# Score: official 41-gene platelet-associated transcriptional score
#
# Friedrich/GSE134051 expression is already log-transformed. This script fits limma directly
# on the log-expression matrix and does not use DESeq2, edgeR, voom, or count
# transformations.
#
# Scientific interpretation:
# - logFC is the expression slope per 1-SD increase in Score_z.
# - Positive coefficients indicate higher expression at higher Score_z.
# - Negative coefficients indicate higher expression at lower Score_z.
# - This is a continuous gene-score association analysis, not a group-based
#   differential-expression comparison.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

message("[1/6] Loading inputs")

###############################################################################
# 1. Analysis constants
###############################################################################

COHORT_ID <- "Friedrich_GSE134051"
COHORT_LABEL <- "Friedrich/GSE134051"
COHORT_SHORT <- "Friedrich"
SOURCE_DATASET <- "GSE134051"
SOURCE_STUDY <- "Friedrich et al."
SOURCE_ASSAY <- "gex.logq"
EXPECTED_TOTAL_N <- 255L
EXPECTED_PRIMARY_N <- 164L
EXPECTED_EXPRESSION_GENES <- 23097L
FDR_THRESHOLD <- 0.05
ABS_LOGFC_THRESHOLD <- 0.25

###############################################################################
# 2. Paths
###############################################################################

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

input_score <- file.path(
  GENERATED_RESULTS_DIR, "Score",
  "Friedrich_GSE134051_primary_master_table_41genes.csv"
)

input_signature <- file.path(
  RESOURCE_DIR,
  "platelet_associated_transcriptional_signature.tsv"
)

tables_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables"
)

limma_dir <- file.path(
  tables_dir,
  "LIMMA"
)

rds_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Objects"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Logs"
)

output_full <- file.path(
  limma_dir,
  "LIMMA_continuous_41genes_full_results.csv"
)

output_no_signature <- file.path(
  limma_dir,
  "LIMMA_continuous_41genes_no_signature_genes.csv"
)

output_significant_no_signature <- file.path(
  limma_dir,
  paste0(
    "LIMMA_continuous_41genes_no_signature_",
    "significant_FDR005_abslogFC025.csv"
  )
)

output_summary <- file.path(
  limma_dir,
  "LIMMA_continuous_41genes_summary.csv"
)

output_metadata <- file.path(
  limma_dir,
  "LIMMA_continuous_41genes_metadata.json"
)

output_fit <- file.path(
  rds_dir,
  "LIMMA_continuous_41genes_fit.rds"
)

output_qc <- file.path(
  logs_dir,
  "LIMMA_continuous_41genes_QC.txt"
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

###############################################################################
# 3. Helpers
###############################################################################

warnings_vec <- character(0)

fail_with_qc <- function(msg) {
  writeLines(
    c(
      "Friedrich/GSE134051 continuous limma association QC",
      "==========================================",
      "",
      paste0("Date: ", as.character(Sys.time())),
      paste0("Project dir: ", project_dir),
      paste0("Failure reason: ", msg),
      "",
      "Final status: FAIL"
    ),
    con = output_qc,
    useBytes = TRUE
  )

  stop(msg, call. = FALSE)
}

add_warning <- function(msg) {
  warnings_vec <<- c(
    warnings_vec,
    msg
  )

  warning(
    msg,
    call. = FALSE
  )
}

require_packages <- function(pkgs) {
  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) > 0) {
    fail_with_qc(
      paste0(
        "Missing required package(s): ",
        paste(missing, collapse = ", ")
      )
    )
  }
}

read_csv_base <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_csv_base <- function(x, path) {
  utils::write.csv(
    x,
    path,
    row.names = FALSE,
    quote = TRUE
  )
}

write_json <- function(x, path) {
  jsonlite::write_json(
    x,
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
    fail_with_qc(
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
    fail_with_qc(
      paste0(
        "Duplicate genes detected in official signature: ",
        paste(duplicated_genes, collapse = ", ")
      )
    )
  }

  if (length(genes) != 41) {
    fail_with_qc(
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
      first_col_chr <- trimws(
        as.character(expr_df[[1]])
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
    fail_with_qc(
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
    fail_with_qc(
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
    fail_with_qc(
      paste0(
        "Possible incorrect input or cohort contamination: ",
        COHORT_LABEL,
        " expression must contain unique GSM sample identifiers only."
      )
    )
  }

  if (!is.numeric(expr_mat)) {
    fail_with_qc(
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
    fail_with_qc(
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
    fail_with_qc(
      "Expression object could not be converted to matrix."
    )
  }

  if (
    nrow(expr_mat) == 0 ||
    ncol(expr_mat) == 0
  ) {
    fail_with_qc(
      "Expression matrix is empty."
    )
  }

  if (
    is.null(rownames(expr_mat)) ||
    any(is.na(rownames(expr_mat))) ||
    any(rownames(expr_mat) == "")
  ) {
    fail_with_qc(
      "Expression matrix must have non-empty gene rownames."
    )
  }

  if (
    is.null(colnames(expr_mat)) ||
    any(is.na(colnames(expr_mat))) ||
    any(colnames(expr_mat) == "")
  ) {
    fail_with_qc(
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
    fail_with_qc(
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
    fail_with_qc(
      paste0(
        "No official signature genes found in expression rownames; ",
        "matrix orientation or gene symbols may be wrong."
      )
    )
  }

  invisible(TRUE)
}

has_non_primary_rows <- function(score_df) {
  sample_type_cols <- intersect(
    c(
      "sample_type",
      "sample_category",
      "sample_class",
      "tissue_type"
    ),
    names(score_df)
  )

  if (length(sample_type_cols) == 0) {
    return(FALSE)
  }

  values <- unlist(
    score_df[sample_type_cols],
    use.names = FALSE
  )

  values <- as.character(values)

  any(
    grepl(
      "BPH|normal|benign",
      values,
      ignore.case = TRUE
    ),
    na.rm = TRUE
  )
}

summarize_results <- function(res) {
  significant_fdr <- (
    !is.na(res$adj.P.Val) &
      res$adj.P.Val < FDR_THRESHOLD
  )

  significant_effect <- (
    significant_fdr &
      abs(res$logFC) > ABS_LOGFC_THRESHOLD
  )

  list(
    n_genes = nrow(res),
    n_fdr005 = sum(
      significant_fdr,
      na.rm = TRUE
    ),
    n_fdr005_abslogfc025 = sum(
      significant_effect,
      na.rm = TRUE
    ),
    n_positive_fdr005 = sum(
      significant_fdr &
        res$logFC > 0,
      na.rm = TRUE
    ),
    n_negative_fdr005 = sum(
      significant_fdr &
        res$logFC < 0,
      na.rm = TRUE
    ),
    n_positive_fdr005_abslogfc025 = sum(
      significant_effect &
        res$logFC > 0,
      na.rm = TRUE
    ),
    n_negative_fdr005_abslogfc025 = sum(
      significant_effect &
        res$logFC < 0,
      na.rm = TRUE
    )
  )
}

write_qc_report <- function(qc) {
  sink(output_qc)

  cat("Friedrich/GSE134051 continuous limma association QC\n")
  cat("==========================================\n\n")

  cat(
    "Date: ",
    as.character(Sys.time()),
    "\n",
    sep = ""
  )

  cat(
    "Project dir: ",
    project_dir,
    "\n\n",
    sep = ""
  )

  cat("1. Inputs\n")

  cat(
    "  Expression input: ",
    qc$inputs$expression,
    "\n",
    sep = ""
  )

  cat(
    "  Score input: ",
    qc$inputs$score,
    "\n",
    sep = ""
  )

  cat(
    "  Signature input: ",
    qc$inputs$signature,
    "\n\n",
    sep = ""
  )

  cat("2. Expression matrix\n")

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
    "  Aligned primary dimensions: ",
    qc$expression$n_genes_aligned,
    " genes x ",
    qc$expression$n_samples_aligned,
    " samples\n\n",
    sep = ""
  )

  cat("3. Samples and score\n")

  cat(
    "  Expected primary samples: ",
    EXPECTED_PRIMARY_N,
    "\n",
    sep = ""
  )

  cat(
    "  Primary samples used: ",
    qc$samples$n_primary_used,
    "\n",
    sep = ""
  )

  cat(
    "  Score_z finite n: ",
    qc$samples$n_score_z_finite,
    "\n",
    sep = ""
  )

  cat(
    "  Score_z min: ",
    qc$samples$score_z_min,
    "\n",
    sep = ""
  )

  cat(
    "  Score_z mean: ",
    qc$samples$score_z_mean,
    "\n",
    sep = ""
  )

  cat(
    "  Score_z sd: ",
    qc$samples$score_z_sd,
    "\n",
    sep = ""
  )

  cat(
    "  Score_z max: ",
    qc$samples$score_z_max,
    "\n\n",
    sep = ""
  )

  cat("4. Filtering\n")

  cat(
    "  Genes removed nonfinite: ",
    qc$filter$n_removed_nonfinite,
    "\n",
    sep = ""
  )

  cat(
    "  Genes removed zero variance: ",
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

  cat("5. Signature handling\n")

  cat(
    "  Signature genes total: ",
    qc$signature$n_total,
    "\n",
    sep = ""
  )

  cat(
    "  Signature genes detected in expression: ",
    qc$signature$n_detected_expression,
    "\n",
    sep = ""
  )

  cat(
    "  Signature genes excluded from no-signature results: ",
    qc$signature$n_excluded_no_signature,
    "\n",
    sep = ""
  )

  cat(
    "  No-signature results contain signature genes: ",
    qc$signature$no_signature_contains_signature,
    "\n",
    sep = ""
  )

  cat(
    "  Significant no-signature table contains signature genes: ",
    qc$signature$significant_contains_signature,
    "\n\n",
    sep = ""
  )

  cat("6. Limma model\n")
  cat("  Method: limma\n")
  cat("  Design: ~ Score_z\n")
  cat("  Coefficient: Score_z\n")
  cat("  Interpretation: logFC is slope per 1 SD increase in Score_z\n")
  cat("  Positive logFC/t_stat: higher expression with higher Score_z\n")
  cat("  Negative logFC/t_stat: higher expression with lower Score_z\n\n")

  cat("7. Full results\n")

  cat(
    "  Genes tested: ",
    qc$results$full$n_genes,
    "\n",
    sep = ""
  )

  cat(
    "  FDR < 0.05: ",
    qc$results$full$n_fdr005,
    "\n",
    sep = ""
  )

  cat(
    "  FDR < 0.05 and abs(logFC) > 0.25: ",
    qc$results$full$n_fdr005_abslogfc025,
    "\n",
    sep = ""
  )

  cat(
    "  Positive FDR < 0.05: ",
    qc$results$full$n_positive_fdr005,
    "\n",
    sep = ""
  )

  cat(
    "  Negative FDR < 0.05: ",
    qc$results$full$n_negative_fdr005,
    "\n\n",
    sep = ""
  )

  cat("8. No-signature results\n")

  cat(
    "  Genes tested: ",
    qc$results$no_signature$n_genes,
    "\n",
    sep = ""
  )

  cat(
    "  FDR < 0.05: ",
    qc$results$no_signature$n_fdr005,
    "\n",
    sep = ""
  )

  cat(
    "  FDR < 0.05 and abs(logFC) > 0.25: ",
    qc$results$no_signature$n_fdr005_abslogfc025,
    "\n",
    sep = ""
  )

  cat(
    "  Positive FDR < 0.05 and abs(logFC) > 0.25: ",
    qc$results$no_signature$n_positive_fdr005_abslogfc025,
    "\n",
    sep = ""
  )

  cat(
    "  Negative FDR < 0.05 and abs(logFC) > 0.25: ",
    qc$results$no_signature$n_negative_fdr005_abslogfc025,
    "\n\n",
    sep = ""
  )

  cat("9. Outputs\n")

  cat(
    "  Full results: ",
    output_full,
    "\n",
    sep = ""
  )

  cat(
    "  No-signature results: ",
    output_no_signature,
    "\n",
    sep = ""
  )

  cat(
    "  Significant no-signature results: ",
    output_significant_no_signature,
    "\n",
    sep = ""
  )

  cat(
    "  Summary: ",
    output_summary,
    "\n",
    sep = ""
  )

  cat(
    "  Metadata: ",
    output_metadata,
    "\n",
    sep = ""
  )

  cat(
    "  Fit RDS: ",
    output_fit,
    "\n\n",
    sep = ""
  )

  cat("10. Warnings\n")

  if (length(qc$warnings) == 0) {
    cat("  None\n")
  } else {
    for (w in qc$warnings) {
      cat(
        "  - ",
        w,
        "\n",
        sep = ""
      )
    }
  }

  cat(
    "\nFinal status: ",
    qc$status,
    "\n",
    sep = ""
  )

  sink()
}

###############################################################################
# 4. Load and validate inputs
###############################################################################

require_packages(
  c(
    "limma",
    "jsonlite"
  )
)

for (
  fp in c(
    input_expression,
    input_score,
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

score_df <- read_csv_base(input_score)
signature_df <- read_canonical_platelet_signature()
signature_genes <- extract_signature_genes(signature_df)

expr_obj <- readRDS(input_expression)
expr_mat <- coerce_expression_matrix(expr_obj)
validate_friedrich_expression_identity(expr_mat)

required_score_cols <- c(
  "sample_name",
  "Score_z"
)

missing_score_cols <- setdiff(
  required_score_cols,
  names(score_df)
)

if (length(missing_score_cols) > 0) {
  fail_with_qc(
    paste0(
      "Score input missing required columns: ",
      paste(missing_score_cols, collapse = ", ")
    )
  )
}

if (has_non_primary_rows(score_df)) {
  fail_with_qc(
    paste0(
      "Score master table appears to contain BPH, normal, or benign samples; ",
      "expected primary prostate tumors only."
    )
  )
}

score_df$sample_name <- trimws(
  as.character(score_df$sample_name)
)

forbidden_id_pattern <- paste0("^ICGC", "_PCA")
if (
  any(!grepl("^GSM[0-9]+$", score_df$sample_name)) ||
  any(grepl(forbidden_id_pattern, score_df$sample_name, ignore.case = TRUE))
) {
  fail_with_qc(
    paste0(
      "Possible incorrect input or cohort contamination: score sample_name ",
      "values must be GSM identifiers from ",
      COHORT_LABEL,
      "."
    )
  )
}

score_df$Score_z <- suppressWarnings(
  as.numeric(score_df$Score_z)
)

if (anyDuplicated(score_df$sample_name)) {
  dup_samples <- unique(
    score_df$sample_name[
      duplicated(score_df$sample_name)
    ]
  )

  fail_with_qc(
    paste0(
      "Duplicate sample_name values in score input: ",
      paste(dup_samples, collapse = ", ")
    )
  )
}

if (
  any(
    is.na(score_df$sample_name) |
    score_df$sample_name == ""
  )
) {
  fail_with_qc(
    "Score input contains missing or empty sample_name values."
  )
}

if (nrow(score_df) != EXPECTED_PRIMARY_N) {
  fail_with_qc(
    paste0(
      "Expected exactly ",
      EXPECTED_PRIMARY_N,
      " Friedrich/GSE134051 primary samples in score master table; observed: ",
      nrow(score_df)
    )
  )
}

n_score_z_finite <- sum(
  is.finite(score_df$Score_z)
)

if (n_score_z_finite != EXPECTED_PRIMARY_N) {
  fail_with_qc(
    paste0(
      "Expected finite Score_z for all ",
      EXPECTED_PRIMARY_N,
      " primary samples; observed finite values: ",
      n_score_z_finite
    )
  )
}

score_metadata <- score_df

score_z_mean <- mean(score_metadata$Score_z)
score_z_sd <- stats::sd(score_metadata$Score_z)

if (
  !is.finite(score_z_sd) ||
  score_z_sd == 0
) {
  fail_with_qc(
    "Score_z has zero or non-finite variance; continuous limma cannot be fitted."
  )
}

if (abs(score_z_mean) > 0.05) {
  add_warning(
    paste0(
      "Score_z mean differs from zero by more than 0.05: ",
      score_z_mean
    )
  )
}

if (abs(score_z_sd - 1) > 0.05) {
  add_warning(
    paste0(
      "Score_z standard deviation differs from one by more than 0.05: ",
      score_z_sd
    )
  )
}

###############################################################################
# 5. Align and filter expression
###############################################################################

message("[2/6] Validating expression matrix")

validate_expression_orientation(
  expr_mat,
  score_metadata$sample_name,
  signature_genes
)

if (anyDuplicated(rownames(expr_mat))) {
  dup_genes <- unique(
    rownames(expr_mat)[
      duplicated(rownames(expr_mat))
    ]
  )

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
    )
  )
}

missing_samples <- setdiff(
  score_metadata$sample_name,
  colnames(expr_mat)
)

if (length(missing_samples) > 0) {
  fail_with_qc(
    paste0(
      "Primary score samples missing from expression matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  )
}

expr_aligned <- expr_mat[
  ,
  score_metadata$sample_name,
  drop = FALSE
]

if (
  !identical(
    colnames(expr_aligned),
    score_metadata$sample_name
  )
) {
  fail_with_qc(
    "Expression columns could not be aligned to score sample_name."
  )
}

message("[3/6] Filtering genes")

keep_finite <- rowSums(
  !is.finite(expr_aligned)
) == 0

n_removed_nonfinite <- sum(
  !keep_finite
)

expr_finite <- expr_aligned[
  keep_finite,
  ,
  drop = FALSE
]

row_variance <- apply(
  expr_finite,
  1,
  stats::var
)

keep_variance <- is.finite(row_variance) &
  row_variance > 0

n_removed_zero_variance <- sum(
  !keep_variance
)

expr_test <- expr_finite[
  keep_variance,
  ,
  drop = FALSE
]

if (nrow(expr_test) < 5000) {
  fail_with_qc(
    paste0(
      "Too few genes available for limma after filtering: ",
      nrow(expr_test)
    )
  )
}

###############################################################################
# 6. Fit continuous limma model
###############################################################################

message("[4/6] Fitting limma continuous model")

design <- stats::model.matrix(
  ~ Score_z,
  data = score_metadata
)

rownames(design) <- score_metadata$sample_name

if (!"Score_z" %in% colnames(design)) {
  fail_with_qc(
    "Design matrix does not contain coefficient Score_z."
  )
}

fit_completed <- FALSE

fit <- tryCatch(
  {
    fit0 <- limma::lmFit(
      expr_test,
      design
    )

    fit1 <- limma::eBayes(fit0)

    fit_completed <- TRUE

    fit1
  },
  error = function(e) {
    fail_with_qc(
      paste0(
        "Continuous limma model failed: ",
        conditionMessage(e)
      )
    )
  }
)

if (!fit_completed) {
  fail_with_qc(
    "Continuous limma model did not complete."
  )
}

tt <- limma::topTable(
  fit,
  coef = "Score_z",
  number = Inf,
  sort.by = "none"
)

required_limma_cols <- c(
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
      "limma topTable output lacks required columns: ",
      paste(missing_limma_cols, collapse = ", ")
    )
  )
}

results <- data.frame(
  gene_name = rownames(tt),
  logFC = as.numeric(tt$logFC),
  AveExpr = as.numeric(tt$AveExpr),
  t_stat = as.numeric(tt$t),
  P.Value = as.numeric(tt$P.Value),
  adj.P.Val = as.numeric(tt$adj.P.Val),
  B = as.numeric(tt$B),
  stringsAsFactors = FALSE
)

results$direction <- ifelse(
  results$logFC > 0,
  "Higher_Score_z",
  ifelse(
    results$logFC < 0,
    "Lower_Score_z",
    "No_direction"
  )
)

results$significant_FDR005 <- (
  !is.na(results$adj.P.Val) &
    results$adj.P.Val < FDR_THRESHOLD
)

results$significant_FDR005_abslogFC025 <- (
  results$significant_FDR005 &
    abs(results$logFC) > ABS_LOGFC_THRESHOLD
)

results$is_signature_gene <- (
  results$gene_name %in% signature_genes
)

results <- results[
  order(
    results$adj.P.Val,
    -abs(results$t_stat),
    results$gene_name
  ),
  ,
  drop = FALSE
]

no_signature_results <- results[
  !results$is_signature_gene,
  ,
  drop = FALSE
]

significant_no_signature_results <- no_signature_results[
  no_signature_results$significant_FDR005_abslogFC025,
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

significant_contains_signature <- any(
  significant_no_signature_results$gene_name %in%
    signature_genes
)

if (no_signature_contains_signature) {
  fail_with_qc(
    "No-signature continuous results still contain official signature genes."
  )
}

if (significant_contains_signature) {
  fail_with_qc(
    paste0(
      "Significant no-signature continuous table still contains ",
      "official signature genes."
    )
  )
}

n_signature_detected <- sum(
  signature_genes %in% rownames(expr_aligned)
)

n_signature_excluded <- sum(
  results$is_signature_gene
)

if (n_signature_detected != n_signature_excluded) {
  fail_with_qc(
    paste0(
      "Signature exclusion mismatch: detected in expression = ",
      n_signature_detected,
      "; excluded from results = ",
      n_signature_excluded,
      "."
    )
  )
}

full_summary <- summarize_results(results)

no_signature_summary <- summarize_results(
  no_signature_results
)

###############################################################################
# 7. Summary, metadata and QC
###############################################################################

summary_df <- data.frame(
  result_set = c(
    "full",
    "no_signature"
  ),
  n_genes = c(
    full_summary$n_genes,
    no_signature_summary$n_genes
  ),
  n_FDR005 = c(
    full_summary$n_fdr005,
    no_signature_summary$n_fdr005
  ),
  n_FDR005_abslogFC025 = c(
    full_summary$n_fdr005_abslogfc025,
    no_signature_summary$n_fdr005_abslogfc025
  ),
  n_positive_FDR005 = c(
    full_summary$n_positive_fdr005,
    no_signature_summary$n_positive_fdr005
  ),
  n_negative_FDR005 = c(
    full_summary$n_negative_fdr005,
    no_signature_summary$n_negative_fdr005
  ),
  n_positive_FDR005_abslogFC025 = c(
    full_summary$n_positive_fdr005_abslogfc025,
    no_signature_summary$n_positive_fdr005_abslogfc025
  ),
  n_negative_FDR005_abslogFC025 = c(
    full_summary$n_negative_fdr005_abslogfc025,
    no_signature_summary$n_negative_fdr005_abslogfc025
  ),
  stringsAsFactors = FALSE
)

qc_status <- if (
  length(warnings_vec) > 0
) {
  "WARNING"
} else {
  "PASS"
}

metadata_json <- list(
  date_time = as.character(Sys.time()),
  project_dir = project_dir,
  project = "friedrich_gse134051",
  cohort = "Friedrich/GSE134051 primary prostate tumors",
  analysis = "continuous limma gene-score association",
  analytic_hierarchy = "complementary robustness analysis",
  score = paste0(
    "official 41-gene platelet-associated ",
    "transcriptional score"
  ),
  model = "expression ~ Score_z",
  coefficient = "Score_z",
  inputs = list(
    expression = input_expression,
    score = input_score,
    signature = input_signature
  ),
  outputs = list(
    full_results = output_full,
    no_signature_results = output_no_signature,
    significant_no_signature_results =
      output_significant_no_signature,
    summary = output_summary,
    metadata = output_metadata,
    fit_rds = output_fit,
    qc = output_qc
  ),
  thresholds = list(
    adjusted_p_value = FDR_THRESHOLD,
    absolute_log2_expression_slope =
      ABS_LOGFC_THRESHOLD,
    signature_genes_excluded_from_interpretation =
      TRUE
  ),
  interpretation = list(
    logFC = "expression slope per 1 SD increase in Score_z",
    positive = paste0(
      "higher expression with higher platelet-associated ",
      "transcriptional score"
    ),
    negative = paste0(
      "higher expression with lower platelet-associated ",
      "transcriptional score"
    )
  ),
  sample_counts = list(
    expected_primary = EXPECTED_PRIMARY_N,
    primary_used = nrow(score_metadata)
  ),
  expression_dimensions = list(
    original_genes = nrow(expr_mat),
    original_samples = ncol(expr_mat),
    aligned_genes = nrow(expr_aligned),
    aligned_samples = ncol(expr_aligned),
    genes_tested = nrow(results)
  ),
  signature_handling = list(
    signature_genes_total = length(signature_genes),
    signature_genes_detected = n_signature_detected,
    signature_genes_excluded = n_signature_excluded
  ),
  association_results = list(
    no_signature_FDR005 =
      no_signature_summary$n_fdr005,
    no_signature_FDR005_abslogFC025 =
      no_signature_summary$n_fdr005_abslogfc025,
    positive_no_signature_FDR005_abslogFC025 =
      no_signature_summary$n_positive_fdr005_abslogfc025,
    negative_no_signature_FDR005_abslogFC025 =
      no_signature_summary$n_negative_fdr005_abslogfc025
  ),
  qc_status = qc_status
)

qc <- list(
  inputs = list(
    expression = input_expression,
    score = input_score,
    signature = input_signature
  ),
  expression = list(
    input_class = attr(expr_mat, "input_class"),
    n_genes_original = nrow(expr_mat),
    n_samples_original = ncol(expr_mat),
    n_genes_aligned = nrow(expr_aligned),
    n_samples_aligned = ncol(expr_aligned)
  ),
  samples = list(
    n_primary_used = nrow(score_metadata),
    n_score_z_finite = n_score_z_finite,
    score_z_min = min(score_metadata$Score_z),
    score_z_mean = score_z_mean,
    score_z_sd = score_z_sd,
    score_z_max = max(score_metadata$Score_z)
  ),
  filter = list(
    n_removed_nonfinite = n_removed_nonfinite,
    n_removed_zero_variance = n_removed_zero_variance,
    n_genes_tested = nrow(results)
  ),
  signature = list(
    n_total = length(signature_genes),
    n_detected_expression = n_signature_detected,
    n_excluded_no_signature = n_signature_excluded,
    no_signature_contains_signature =
      no_signature_contains_signature,
    significant_contains_signature =
      significant_contains_signature
  ),
  results = list(
    full = full_summary,
    no_signature = no_signature_summary
  ),
  warnings = warnings_vec,
  status = qc_status
)

###############################################################################
# 8. Write outputs
###############################################################################

message("[5/6] Writing outputs")

write_csv_base(
  results,
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

write_json(
  metadata_json,
  output_metadata
)

saveRDS(
  fit,
  output_fit
)

write_qc_report(qc)

message(
  "[6/6] Final status: ",
  qc_status
)
