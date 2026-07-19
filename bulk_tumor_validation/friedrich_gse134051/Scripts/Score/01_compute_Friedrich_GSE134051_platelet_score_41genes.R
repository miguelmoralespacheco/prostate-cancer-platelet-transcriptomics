#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Official 41-gene platelet-associated transcriptional score
#
# This script calculates the official 41-gene platelet-associated transcriptional
# score in Friedrich/GSE134051 using primary prostate cancer samples only.

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

###############################################################################
# 1. Paths
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

input_expression <- file.path(INPUT_DIR, "LocalLarge", "Friedrich_GSE134051_expression_logq.rds")
input_clinical   <- file.path(INPUT_DIR, "LocalLarge", "Friedrich_GSE134051_clinical.rds")
input_signature  <- file.path(RESOURCE_DIR, "platelet_associated_transcriptional_signature.tsv")

output_dir <- file.path(GENERATED_RESULTS_DIR, "Score")

output_score <- file.path(output_dir, "Friedrich_GSE134051_primary_platelet_score_41genes.csv")
output_master <- file.path(output_dir, "Friedrich_GSE134051_primary_master_table_41genes.csv")
output_q1q4 <- file.path(output_dir, "Friedrich_GSE134051_primary_Q1Q4_groups_41genes.csv")
output_mapping <- file.path(output_dir, "Friedrich_GSE134051_primary_platelet_score_41genes_gene_mapping.csv")
output_summary <- file.path(output_dir, "Friedrich_GSE134051_primary_platelet_score_41genes_summary.csv")
output_metadata <- file.path(output_dir, "Friedrich_GSE134051_primary_platelet_score_41genes_metadata.json")
output_qc <- file.path(output_dir, "Friedrich_GSE134051_primary_platelet_score_41genes_QC.txt")

###############################################################################
# 2. Helpers
###############################################################################

warnings_vec <- character(0)

fail_clear <- function(msg) {
  stop(msg, call. = FALSE)
}

add_warning <- function(msg) {
  warnings_vec <<- c(warnings_vec, msg)
  warning(msg, call. = FALSE)
}

write_csv_base <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE)
}

read_signature_table <- function(path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    return(readr::read_csv(path, show_col_types = FALSE))
  }

  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_metadata_json <- function(metadata, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    fail_clear("Package 'jsonlite' is required to write metadata JSON.")
  }

  jsonlite::write_json(
    metadata,
    path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
}

is_default_rownames <- function(x) {
  rn <- rownames(x)
  is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))
}

coerce_expression_matrix <- function(expr_obj) {
  expr_class <- paste(class(expr_obj), collapse = ",")

  if (is.matrix(expr_obj)) {
    expr_mat <- expr_obj
  } else if (is.data.frame(expr_obj)) {
    expr_df <- expr_obj

    if (is_default_rownames(expr_df) && ncol(expr_df) >= 2) {
      first_col <- expr_df[[1]]
      first_col_chr <- trimws(as.character(first_col))

      first_col_nonempty <- !any(is.na(first_col_chr) | first_col_chr == "")
      first_col_numeric <- suppressWarnings(all(!is.na(as.numeric(first_col_chr))))
      first_col_unique <- !anyDuplicated(first_col_chr)

      if (first_col_nonempty && !first_col_numeric && first_col_unique) {
        rownames(expr_df) <- first_col_chr
        expr_df <- expr_df[, -1, drop = FALSE]

        add_warning(
          "Expression data.frame had default rownames; first non-numeric unique column was used as gene rownames."
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

validate_expression_matrix <- function(expr_mat) {
  if (!is.matrix(expr_mat)) {
    fail_clear("Expression object could not be converted to matrix.")
  }

  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    fail_clear("Expression matrix is empty.")
  }

  if (
    is.null(rownames(expr_mat)) ||
    any(is.na(rownames(expr_mat))) ||
    any(rownames(expr_mat) == "")
  ) {
    fail_clear("Expression matrix must have non-empty gene rownames.")
  }

  if (
    is.null(colnames(expr_mat)) ||
    any(is.na(colnames(expr_mat))) ||
    any(colnames(expr_mat) == "")
  ) {
    fail_clear("Expression matrix must have non-empty sample colnames.")
  }

  finite_fraction <- mean(is.finite(expr_mat))

  if (!is.finite(finite_fraction) || finite_fraction < 0.80) {
    fail_clear(sprintf("Too few finite expression values: %.3f", finite_fraction))
  }

  invisible(TRUE)
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
    length(sample_ids) != EXPECTED_TOTAL_N ||
    any(!grepl("^GSM[0-9]+$", sample_ids)) ||
    any(grepl(forbidden_id_pattern, sample_ids, ignore.case = TRUE))
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: ",
        COHORT_LABEL,
        " expression must contain only GSM sample identifiers and no foreign cohort identifiers."
      )
    )
  }

  if (!is.numeric(expr_mat)) {
    fail_clear(
      paste0(
        "Possible incorrect input: ",
        COHORT_LABEL,
        " expression must be a numeric continuous log-expression matrix."
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
        COHORT_LABEL,
        " assay ",
        SOURCE_ASSAY,
        " must be continuous log-expression, not integer counts."
      )
    )
  }

  invisible(TRUE)
}

validate_friedrich_clinical_identity <- function(clinical_df) {
  forbidden_id_pattern <- paste0("ICGC", "_PCA")
  required_identity_columns <- c(
    "sample_name",
    "sample_type",
    "study_name"
  )
  missing_identity_columns <- setdiff(
    required_identity_columns,
    names(clinical_df)
  )

  if (length(missing_identity_columns) > 0) {
    fail_clear(
      paste0(
        "Possible incorrect clinical input: missing Friedrich identity columns: ",
        paste(missing_identity_columns, collapse = ", "),
        "."
      )
    )
  }

  if (nrow(clinical_df) != EXPECTED_TOTAL_N) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: expected ",
        EXPECTED_TOTAL_N,
        " clinical records for ",
        COHORT_LABEL,
        "; observed ",
        nrow(clinical_df),
        "."
      )
    )
  }

  sample_ids <- trimws(as.character(clinical_df$sample_name))
  if (
    any(!grepl("^GSM[0-9]+$", sample_ids)) ||
    any(grepl(forbidden_id_pattern, sample_ids, ignore.case = TRUE)) ||
    anyDuplicated(sample_ids)
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: ",
        COHORT_LABEL,
        " clinical sample_name values must be unique GSM identifiers."
      )
    )
  }

  if ("alt_sample_name" %in% names(clinical_df)) {
    alt_ids <- trimws(as.character(clinical_df$alt_sample_name))
    if (
      any(is.na(alt_ids) | !grepl("^RIB", alt_ids)) ||
      any(grepl(forbidden_id_pattern, alt_ids, ignore.case = TRUE))
    ) {
      fail_clear(
        paste0(
          "Possible incorrect input or cohort contamination: ",
          COHORT_LABEL,
          " alt_sample_name values must be RIB identifiers."
        )
      )
    }
  }

  study_name <- trimws(as.character(clinical_df$study_name))
  if (
    any(is.na(study_name) | study_name == "") ||
    any(!grepl("Friedrich", study_name, ignore.case = TRUE))
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: study_name is not compatible with ",
        SOURCE_STUDY,
        "."
      )
    )
  }

  sample_type <- trimws(as.character(clinical_df$sample_type))
  observed_counts <- table(
    factor(sample_type, levels = c("primary", "normal", "BPH"))
  )
  expected_counts <- c(primary = 164L, normal = 52L, BPH = 39L)

  if (
    any(!sample_type %in% names(expected_counts)) ||
    !identical(as.integer(observed_counts), as.integer(expected_counts))
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: expected sample_type distribution ",
        "primary=164, normal=52, BPH=39 for ",
        COHORT_LABEL,
        "."
      )
    )
  }

  invisible(TRUE)
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

extract_official_genes <- function(signature_df) {
  gene_col <- detect_gene_column(signature_df)

  genes_raw <- trimws(as.character(signature_df[[gene_col]]))
  genes <- genes_raw[!is.na(genes_raw) & genes_raw != ""]

  duplicated_genes <- unique(genes[duplicated(genes)])

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

orient_expression_if_needed <- function(expr_mat, official_genes) {
  input_class <- attr(expr_mat, "input_class")

  row_hits <- sum(official_genes %in% rownames(expr_mat))
  col_hits <- sum(official_genes %in% colnames(expr_mat))

  if (row_hits == 0 && col_hits >= 10) {
    fail_clear(
      paste0(
        "Expression matrix appears transposed: 0 official genes in rownames and ",
        col_hits,
        " in colnames. Expected genes as rows and samples as columns. ",
        "Please inspect the expression object before scoring."
      )
    )
  }

  if (col_hits > row_hits && col_hits >= 10) {
    add_warning(
      paste0(
        "Expression matrix has more official gene matches in colnames (",
        col_hits,
        ") than rownames (",
        row_hits,
        "). Keeping original orientation, but this should be manually reviewed."
      )
    )
  }

  if (row_hits == 0) {
    fail_clear(
      paste0(
        "No official platelet-associated signature genes were found in expression rownames. ",
        "Cannot compute score."
      )
    )
  }

  attr(expr_mat, "input_class") <- input_class

  expr_mat
}

build_gene_mapping <- function(expr_genes, official_genes) {
  duplicated_expr_genes <- unique(expr_genes[duplicated(expr_genes)])

  if (length(duplicated_expr_genes) > 0) {
    fail_clear(
      paste0(
        "Duplicate gene symbols detected in expression rownames: ",
        paste(utils::head(duplicated_expr_genes, 25), collapse = ", "),
        if (length(duplicated_expr_genes) > 25) " ..." else ""
      )
    )
  }

  expr_lower <- tolower(expr_genes)
  official_lower <- tolower(official_genes)

  ci_match <- vapply(
    seq_along(official_genes),
    function(i) {
      hits <- expr_genes[expr_lower == official_lower[i]]

      if (length(hits) == 0) {
        return(NA_character_)
      }

      paste(hits, collapse = ";")
    },
    character(1)
  )

  matched_exact <- ifelse(
    official_genes %in% expr_genes,
    official_genes,
    NA_character_
  )

  data.frame(
    official_gene = official_genes,
    present_exact = official_genes %in% expr_genes,
    matched_expression_gene = matched_exact,
    present_case_insensitive = !is.na(ci_match),
    case_insensitive_match = ci_match,
    used_for_score = official_genes %in% expr_genes,
    stringsAsFactors = FALSE
  )
}

derive_qc_status <- function(n_used) {
  if (n_used >= 35) {
    return("PASS")
  }

  if (n_used >= 30) {
    return("WARNING")
  }

  "FAIL"
}

derive_binary_event <- function(x, variable_name) {
  if (is.numeric(x) || is.integer(x)) {
    return(ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L)))
  }

  x_chr <- tolower(trimws(as.character(x)))

  event <- rep(NA_integer_, length(x_chr))
  event[x_chr %in% c("1", "yes", "true", "dead", "death", "deceased", "event")] <- 1L
  event[x_chr %in% c("0", "no", "false", "alive", "living", "censored")] <- 0L

  if (all(is.na(event))) {
    add_warning(
      paste0(
        "Could not confidently derive binary event variable from ",
        variable_name,
        "."
      )
    )
  }

  event
}

assign_primary_q1q4 <- function(df) {
  df <- df[order(df$Score_z, df$sample_name), , drop = FALSE]

  n_primary <- nrow(df)
  q_n <- floor(n_primary / 4)

  if (q_n < 1) {
    fail_clear("Too few primary samples to define Q1/Q4 groups.")
  }

  df$score_quartile <- "Q2Q3"
  df$platelet_score_group <- "MID_Q2Q3"

  low_idx <- seq_len(q_n)
  high_idx <- seq.int(n_primary - q_n + 1, n_primary)

  df$score_quartile[low_idx] <- "Q1"
  df$platelet_score_group[low_idx] <- "LOW_Q1"

  df$score_quartile[high_idx] <- "Q4"
  df$platelet_score_group[high_idx] <- "HIGH_Q4"

  df <- df[order(df$sample_name), , drop = FALSE]

  df
}

write_qc_report <- function(path, qc) {
  sink(path)

  cat("Friedrich/GSE134051 official 41-gene platelet-associated transcriptional score QC\n")
  cat("=====================================================================\n\n")

  cat("Input paths:\n")
  cat("  Expression: ", qc$input_expression, "\n", sep = "")
  cat("  Clinical:   ", qc$input_clinical, "\n", sep = "")
  cat("  Signature:  ", qc$input_signature, "\n\n", sep = "")

  cat("Expression:\n")
  cat("  Input class: ", qc$expression_class, "\n", sep = "")
  cat(
    "  Original dimensions: ",
    qc$n_expression_genes_original,
    " genes x ",
    qc$n_expression_samples_original,
    " samples\n",
    sep = ""
  )
  cat(
    "  Primary-only dimensions: ",
    qc$n_expression_genes_primary,
    " genes x ",
    qc$n_primary_samples,
    " samples\n",
    sep = ""
  )
  cat("  Finite fraction, primary expression: ", qc$finite_fraction_primary, "\n\n", sep = "")

  cat("Clinical filtering:\n")
  cat("  Clinical input class: ", qc$clinical_class, "\n", sep = "")
  cat("  Clinical samples total: ", qc$n_clinical_samples_total, "\n", sep = "")
  cat("  Score-clinical expression intersection: ", qc$n_expression_clinical_intersection, "\n", sep = "")
  cat("  BPH samples excluded: ", qc$n_bph_excluded, "\n", sep = "")
  cat("  Normal samples excluded: ", qc$n_normal_excluded, "\n", sep = "")
  cat("  Primary samples retained: ", qc$n_primary_samples, "\n\n", sep = "")

  cat("Signature mapping:\n")
  cat("  Signature genes total: ", qc$n_official_signature_genes, "\n", sep = "")
  cat("  Genes used: ", qc$n_genes_used_for_score, "\n", sep = "")
  cat("  Genes missing: ", qc$n_missing_signature_genes, "\n", sep = "")
  cat("  Coverage %: ", qc$coverage_pct, "\n\n", sep = "")

  cat("Genes missing list:\n")
  if (length(qc$genes_missing) == 0) {
    cat("  None\n\n")
  } else {
    cat("  ", paste(qc$genes_missing, collapse = ", "), "\n\n", sep = "")
  }

  cat("Score summary, primary samples only:\n")
  if (!is.null(qc$score_summary)) {
    cat("  Score_raw min: ", qc$score_summary$score_raw_min, "\n", sep = "")
    cat("  Score_raw median: ", qc$score_summary$score_raw_median, "\n", sep = "")
    cat("  Score_raw max: ", qc$score_summary$score_raw_max, "\n", sep = "")
    cat("  Score_z mean: ", qc$score_summary$score_z_mean, "\n", sep = "")
    cat("  Score_z sd: ", qc$score_summary$score_z_sd, "\n\n", sep = "")
  } else {
    cat("  Not computed because QC status is FAIL.\n\n")
  }

  cat("Primary Q1/Q4 groups:\n")
  if (!is.null(qc$q1q4_summary)) {
    cat("  LOW_Q1: ", qc$q1q4_summary$n_LOW_Q1, "\n", sep = "")
    cat("  MID_Q2Q3: ", qc$q1q4_summary$n_MID_Q2Q3, "\n", sep = "")
    cat("  HIGH_Q4: ", qc$q1q4_summary$n_HIGH_Q4, "\n", sep = "")
    cat("  Q1Q4 total: ", qc$q1q4_summary$n_Q1Q4_total, "\n\n", sep = "")
  } else {
    cat("  Not available.\n\n")
  }

  cat("Overall survival, primary samples:\n")
  if (!is.null(qc$os_summary)) {
    cat("  OS usable: ", qc$os_summary$n_OS_usable, "\n", sep = "")
    cat("  OS events: ", qc$os_summary$n_OS_events, "\n", sep = "")
    cat("  OS censored: ", qc$os_summary$n_OS_censored, "\n\n", sep = "")
  } else {
    cat("  Not available.\n\n")
  }

  cat("Warnings:\n")
  if (length(qc$warnings) == 0) {
    cat("  None\n\n")
  } else {
    for (w in qc$warnings) {
      cat("  - ", w, "\n", sep = "")
    }
    cat("\n")
  }

  cat("Final status: ", qc$qc_status, "\n", sep = "")

  sink()
}

###############################################################################
# 3. Input checks
###############################################################################

if (!file.exists(input_expression)) {
  fail_clear(paste0("Missing expression input: ", input_expression))
}

if (!file.exists(input_clinical)) {
  fail_clear(paste0("Missing clinical input: ", input_clinical))
}

if (!file.exists(input_signature)) {
  fail_clear(paste0("Missing signature input: ", input_signature))
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  fail_clear("Package 'jsonlite' is required to write metadata JSON.")
}

###############################################################################
# 4. Load expression, clinical metadata and signature
###############################################################################

expr_obj <- readRDS(input_expression)
expr_mat <- coerce_expression_matrix(expr_obj)
validate_expression_matrix(expr_mat)

n_expression_genes_original <- nrow(expr_mat)
n_expression_samples_original <- ncol(expr_mat)

clinical_obj <- readRDS(input_clinical)

clinical_class <- paste(class(clinical_obj), collapse = ",")

clinical_df <- tryCatch(
  {
    as.data.frame(clinical_obj, stringsAsFactors = FALSE)
  },
  error = function(e) {
    fail_clear(
      paste0(
        "Clinical object could not be converted to data.frame. Observed class: ",
        clinical_class,
        ". Error: ",
        conditionMessage(e)
      )
    )
  }
)

attr(clinical_df, "input_class") <- clinical_class

if (!is.data.frame(clinical_df)) {
  fail_clear(
    paste0(
      "Clinical object conversion did not return a data.frame. Observed class: ",
      clinical_class
    )
  )
}

validate_friedrich_clinical_identity(clinical_df)

if (!"sample_name" %in% names(clinical_df)) {
  fail_clear("Clinical table must contain column: sample_name")
}

if (!"sample_type" %in% names(clinical_df)) {
  fail_clear("Clinical table must contain column: sample_type")
}

if (anyDuplicated(clinical_df$sample_name) > 0) {
  fail_clear("Duplicated sample_name values detected in clinical table.")
}

if (anyDuplicated(colnames(expr_mat)) > 0) {
  fail_clear("Duplicated sample names detected in expression matrix colnames.")
}

signature_df <- read_canonical_platelet_signature()
official_genes <- extract_official_genes(signature_df)

expr_mat <- orient_expression_if_needed(expr_mat, official_genes)
validate_expression_matrix(expr_mat)
validate_friedrich_expression_identity(expr_mat)

###############################################################################
# 5. Keep primary samples only
###############################################################################

message("[2/6] Filtering primary samples")

clinical_df$sample_type <- trimws(as.character(clinical_df$sample_type))

n_clinical_samples_total <- nrow(clinical_df)

samples_in_both <- intersect(colnames(expr_mat), clinical_df$sample_name)
n_expression_clinical_intersection <- length(samples_in_both)

expression_only_samples <- setdiff(
  colnames(expr_mat),
  clinical_df$sample_name
)
clinical_only_samples <- setdiff(
  clinical_df$sample_name,
  colnames(expr_mat)
)

if (
  n_expression_clinical_intersection != EXPECTED_TOTAL_N ||
  length(expression_only_samples) > 0 ||
  length(clinical_only_samples) > 0
) {
  fail_clear(
    paste0(
      "Possible incorrect input or cohort contamination: expression and clinical ",
      "records must match completely for all ",
      EXPECTED_TOTAL_N,
      " ",
      COHORT_LABEL,
      " samples."
    )
  )
}

clinical_overlap <- clinical_df[clinical_df$sample_name %in% samples_in_both, , drop = FALSE]

n_bph_excluded <- sum(clinical_overlap$sample_type == "BPH", na.rm = TRUE)
n_normal_excluded <- sum(clinical_overlap$sample_type == "normal", na.rm = TRUE)

primary_clinical <- clinical_overlap[
  clinical_overlap$sample_type == "primary",
  ,
  drop = FALSE
]

if (nrow(primary_clinical) == 0) {
  fail_clear("No primary samples detected using sample_type == 'primary'.")
}

if (nrow(primary_clinical) != EXPECTED_PRIMARY_N) {
  fail_clear(
    paste0(
      "Possible incorrect input or cohort contamination: expected ",
      EXPECTED_PRIMARY_N,
      " primary samples for ",
      COHORT_LABEL,
      "; observed ",
      nrow(primary_clinical)
    )
  )
}

primary_samples <- primary_clinical$sample_name

missing_primary_in_expression <- setdiff(primary_samples, colnames(expr_mat))
if (length(missing_primary_in_expression) > 0) {
  fail_clear(
    paste0(
      "Primary clinical samples missing from expression matrix: ",
      paste(missing_primary_in_expression, collapse = ", ")
    )
  )
}

expr_primary <- expr_mat[, primary_samples, drop = FALSE]
validate_expression_matrix(expr_primary)

if (
  !identical(colnames(expr_primary), primary_samples) ||
  anyDuplicated(primary_samples)
) {
  fail_clear(
    paste0(
      "Possible incorrect input or cohort contamination: primary expression and ",
      "clinical samples are not uniquely aligned in the same order."
    )
  )
}

finite_fraction_primary <- round(mean(is.finite(expr_primary)), 6)

###############################################################################
# 6. Gene mapping
###############################################################################

message("[3/6] Mapping genes")

expr_genes <- rownames(expr_primary)
gene_mapping <- build_gene_mapping(expr_genes, official_genes)

n_official <- length(official_genes)
n_present_exact <- sum(gene_mapping$present_exact)

genes_used <- gene_mapping$official_gene[gene_mapping$used_for_score]
genes_missing <- gene_mapping$official_gene[!gene_mapping$used_for_score]

n_used <- length(genes_used)
n_missing <- length(genes_missing)
coverage_pct <- round(100 * n_used / n_official, 2)

qc_status <- derive_qc_status(n_used)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write_csv_base(gene_mapping, output_mapping)

expression_class <- attr(expr_mat, "input_class")
if (is.null(expression_class) || is.na(expression_class) || expression_class == "") {
  expression_class <- paste(class(expr_mat), collapse = ",")
}

qc_base <- list(
  input_expression = input_expression,
  input_clinical = input_clinical,
  input_signature = input_signature,
  expression_class = expression_class,
  clinical_class = attr(clinical_df, "input_class"),
  n_expression_genes_original = n_expression_genes_original,
  n_expression_samples_original = n_expression_samples_original,
  n_expression_genes_primary = nrow(expr_primary),
  n_primary_samples = ncol(expr_primary),
  finite_fraction_primary = finite_fraction_primary,
  n_clinical_samples_total = n_clinical_samples_total,
  n_expression_clinical_intersection = n_expression_clinical_intersection,
  n_bph_excluded = n_bph_excluded,
  n_normal_excluded = n_normal_excluded,
  n_official_signature_genes = n_official,
  n_present_exact = n_present_exact,
  n_genes_used_for_score = n_used,
  n_missing_signature_genes = n_missing,
  coverage_pct = coverage_pct,
  genes_missing = genes_missing,
  score_summary = NULL,
  q1q4_summary = NULL,
  os_summary = NULL,
  warnings = warnings_vec,
  qc_status = qc_status
)

if (qc_status == "FAIL") {
  write_qc_report(output_qc, qc_base)

  fail_clear(
    paste0(
      "Final status: FAIL. Only gene mapping and QC were written. Genes used: ",
      n_used
    )
  )
}

if (qc_status == "WARNING") {
  add_warning(
    "QC status is WARNING because only 30-34 official genes were available for scoring."
  )

  qc_base$warnings <- warnings_vec
}

###############################################################################
# 7. Compute primary-only score
###############################################################################

message("[4/6] Computing primary-only score")

score_raw <- colMeans(expr_primary[genes_used, , drop = FALSE], na.rm = TRUE)

if (any(!is.finite(score_raw)) || any(is.na(score_raw))) {
  fail_clear("Score_raw contains NA or non-finite values.")
}

score_z <- as.numeric(scale(score_raw))

if (any(!is.finite(score_z)) || any(is.na(score_z))) {
  fail_clear("Score_z contains NA or non-finite values.")
}

score_df <- data.frame(
  sample_name = names(score_raw),
  Score_raw = as.numeric(score_raw),
  Score_z = score_z,
  stringsAsFactors = FALSE
)

score_df <- score_df[order(score_df$sample_name), , drop = FALSE]

primary_master <- merge(
  primary_clinical,
  score_df,
  by = "sample_name",
  all = FALSE,
  sort = FALSE
)

if (nrow(primary_master) != nrow(primary_clinical)) {
  fail_clear("Primary master table lost samples during clinical-score merge.")
}

primary_master <- assign_primary_q1q4(primary_master)

q1q4_df <- primary_master[
  primary_master$platelet_score_group %in% c("LOW_Q1", "HIGH_Q4"),
  ,
  drop = FALSE
]

group_counts <- table(primary_master$platelet_score_group)

required_groups <- c("LOW_Q1", "MID_Q2Q3", "HIGH_Q4")
if (!all(required_groups %in% names(group_counts))) {
  fail_clear("Q1/Q4 group assignment failed.")
}

expected_qn <- floor(nrow(primary_master) / 4)

if (
  as.integer(group_counts[["LOW_Q1"]]) != expected_qn ||
  as.integer(group_counts[["HIGH_Q4"]]) != expected_qn
) {
  fail_clear("Q1/Q4 group sizes are inconsistent with expected quartile size.")
}

###############################################################################
# 8. Overall survival variables
###############################################################################

if ("overall_survival_status" %in% names(primary_master)) {
  primary_master$OS_event <- derive_binary_event(
    primary_master$overall_survival_status,
    "overall_survival_status"
  )
} else {
  add_warning("Column overall_survival_status not found; OS_event was not created.")
  primary_master$OS_event <- NA_integer_
}

if ("days_to_overall_survival" %in% names(primary_master)) {
  primary_master$OS_time_days <- suppressWarnings(
    as.numeric(primary_master$days_to_overall_survival)
  )
  primary_master$OS_time_months <- primary_master$OS_time_days / 30.4375
} else {
  add_warning("Column days_to_overall_survival not found; OS time variables were not created.")
  primary_master$OS_time_days <- NA_real_
  primary_master$OS_time_months <- NA_real_
}

q1q4_df <- primary_master[
  primary_master$platelet_score_group %in% c("LOW_Q1", "HIGH_Q4"),
  ,
  drop = FALSE
]

n_os_usable <- sum(
  is.finite(primary_master$OS_time_months) &
    primary_master$OS_time_months > 0 &
    !is.na(primary_master$OS_event)
)

n_os_events <- sum(primary_master$OS_event == 1L, na.rm = TRUE)
n_os_censored <- sum(primary_master$OS_event == 0L, na.rm = TRUE)

score_summary <- list(
  score_raw_min = min(score_df$Score_raw),
  score_raw_median = stats::median(score_df$Score_raw),
  score_raw_max = max(score_df$Score_raw),
  score_z_mean = mean(score_df$Score_z),
  score_z_sd = stats::sd(score_df$Score_z)
)

q1q4_summary <- list(
  n_LOW_Q1 = as.integer(group_counts[["LOW_Q1"]]),
  n_MID_Q2Q3 = as.integer(group_counts[["MID_Q2Q3"]]),
  n_HIGH_Q4 = as.integer(group_counts[["HIGH_Q4"]]),
  n_Q1Q4_total = nrow(q1q4_df)
)

os_summary <- list(
  n_OS_usable = n_os_usable,
  n_OS_events = n_os_events,
  n_OS_censored = n_os_censored
)

summary_df <- data.frame(
  metric = c(
    "cohort",
    "score_name",
    "analytic_universe",
    "expression_input",
    "clinical_input",
    "signature_input",
    "n_expression_genes_original",
    "n_expression_samples_original",
    "n_clinical_samples_total",
    "n_expression_clinical_intersection",
    "n_bph_excluded",
    "n_normal_excluded",
    "n_primary_samples",
    "n_expression_genes_primary",
    "n_official_signature_genes",
    "n_present_exact",
    "n_genes_used_for_score",
    "n_missing_signature_genes",
    "coverage_pct",
    "score_raw_min",
    "score_raw_median",
    "score_raw_max",
    "score_z_mean",
    "score_z_sd",
    "n_LOW_Q1",
    "n_MID_Q2Q3",
    "n_HIGH_Q4",
    "n_Q1Q4_total",
    "n_OS_usable",
    "n_OS_events",
    "n_OS_censored",
    "qc_status"
  ),
  value = as.character(
    c(
      "Friedrich/GSE134051",
      "official_41_gene_platelet_associated_transcriptional_score",
      "primary prostate cancer samples only; sample_type == 'primary'",
      input_expression,
      input_clinical,
      input_signature,
      n_expression_genes_original,
      n_expression_samples_original,
      n_clinical_samples_total,
      n_expression_clinical_intersection,
      n_bph_excluded,
      n_normal_excluded,
      ncol(expr_primary),
      nrow(expr_primary),
      n_official,
      n_present_exact,
      n_used,
      n_missing,
      coverage_pct,
      score_summary$score_raw_min,
      score_summary$score_raw_median,
      score_summary$score_raw_max,
      score_summary$score_z_mean,
      score_summary$score_z_sd,
      q1q4_summary$n_LOW_Q1,
      q1q4_summary$n_MID_Q2Q3,
      q1q4_summary$n_HIGH_Q4,
      q1q4_summary$n_Q1Q4_total,
      os_summary$n_OS_usable,
      os_summary$n_OS_events,
      os_summary$n_OS_censored,
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
    clinical = input_clinical,
    signature = input_signature
  ),
  outputs = list(
    score = output_score,
    master = output_master,
    q1q4 = output_q1q4,
    gene_mapping = output_mapping,
    summary = output_summary,
    metadata = output_metadata,
    qc = output_qc
  ),
  analytic_universe = list(
    cohort = "Friedrich/GSE134051",
    filter = "sample_type == 'primary'",
    excluded_sample_types = c("BPH", "normal"),
    n_primary_samples = ncol(expr_primary)
  ),
  expression_dimensions = list(
    original_genes = n_expression_genes_original,
    original_samples = n_expression_samples_original,
    primary_genes = nrow(expr_primary),
    primary_samples = ncol(expr_primary)
  ),
  score_formula = paste(
    "Score_raw = mean log-transformed expression across official genes present exactly;",
    "Score_z = z-score of Score_raw across primary Friedrich/GSE134051 samples only"
  ),
  q1q4_rule = paste(
    "Samples ordered by primary-only Score_z;",
    "bottom floor(n/4) assigned LOW_Q1 and top floor(n/4) assigned HIGH_Q4"
  ),
  official_genes = official_genes,
  genes_used = genes_used,
  genes_missing = genes_missing,
  q1q4_summary = q1q4_summary,
  os_summary = os_summary,
  qc_status = qc_status
)

qc_full <- qc_base
qc_full$warnings <- warnings_vec
qc_full$qc_status <- qc_status
qc_full$score_summary <- score_summary
qc_full$q1q4_summary <- q1q4_summary
qc_full$os_summary <- os_summary

###############################################################################
# 9. Write outputs
###############################################################################

message("[5/6] Writing outputs")

write_csv_base(score_df, output_score)
write_csv_base(primary_master, output_master)
write_csv_base(q1q4_df, output_q1q4)
write_csv_base(summary_df, output_summary)
write_metadata_json(metadata, output_metadata)
write_qc_report(output_qc, qc_full)

message("[6/6] Final status: ", qc_status)
