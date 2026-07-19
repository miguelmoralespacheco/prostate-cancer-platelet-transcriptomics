#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Q1/Q4 Hallmark GSEA
#
# This script runs Hallmark GSEA for Friedrich/GSE134051 Q1/Q4 using limma results from the
# official 41-gene platelet-associated transcriptional score workflow.
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

input_limma <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4",
  "Tables",
  "LIMMA",
  "LIMMA_Q1Q4_41genes_no_signature_genes.csv"
)
input_signature <- file.path(RESOURCE_DIR, "platelet_associated_transcriptional_signature.tsv")

gsea_dir <- file.path(GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "GSEA")
logs_dir <- file.path(GENERATED_RESULTS_DIR, "Q1Q4", "Logs")

output_full <- file.path(gsea_dir, "GSEA_Q1Q4_no_signature_Hallmark_full.csv")
output_fdr005 <- file.path(gsea_dir, "GSEA_Q1Q4_no_signature_Hallmark_FDR005.csv")
output_emt <- file.path(gsea_dir, "GSEA_Q1Q4_no_signature_Hallmark_EMT_row.csv")
output_metadata <- file.path(gsea_dir, "GSEA_Q1Q4_metadata.json")
output_qc <- file.path(logs_dir, "GSEA_Q1Q4_QC.txt")

dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

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

read_csv_base <- function(path) {
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
empty_gsea_table <- function() {
  data.frame(
    pathway = character(0),
    pval = numeric(0),
    padj = numeric(0),
    log2err = numeric(0),
    ES = numeric(0),
    NES = numeric(0),
    nMoreExtreme = numeric(0),
    size = integer(0),
    leadingEdge = character(0),
    direction = character(0),
    interpretation = character(0),
    stringsAsFactors = FALSE
  )
}

detect_gene_column <- function(signature_df) {
  candidates <- c("gene_name", "gene_symbol", "symbol", "gene", "external_gene_name")
  lower_names <- tolower(names(signature_df))
  idx <- match(candidates, lower_names)
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0) {
    fail_clear(paste0(
      "Could not detect gene symbol column in signature. Accepted names: ",
      paste(candidates, collapse = ", ")
    ))
  }

  names(signature_df)[idx[1]]
}

extract_signature_genes <- function(signature_df) {
  gene_col <- detect_gene_column(signature_df)
  genes <- trimws(as.character(signature_df[[gene_col]]))
  genes <- genes[!is.na(genes) & genes != ""]

  duplicated_genes <- unique(genes[duplicated(genes)])
  if (length(duplicated_genes) > 0) {
    fail_clear(paste0(
      "Duplicate genes detected in official signature: ",
      paste(duplicated_genes, collapse = ", ")
    ))
  }

  if (length(genes) != 41) {
    fail_clear(paste0("Official signature must contain exactly 41 unique genes. Observed: ", length(genes)))
  }

  genes
}

coerce_signature_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  x_chr <- tolower(trimws(as.character(x)))
  ifelse(x_chr %in% c("true", "t", "1", "yes"), TRUE,
         ifelse(x_chr %in% c("false", "f", "0", "no"), FALSE, NA))
}

load_hallmark_sets <- function() {
  hallmark <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
    error = function(e1) {
      tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
        error = function(e2) {
          fail_clear(paste0(
            "Could not load MSigDB Hallmark gene sets. category error: ",
            conditionMessage(e1), "; collection error: ", conditionMessage(e2)
          ))
        }
      )
    }
  )

  if (!all(c("gs_name", "gene_symbol") %in% names(hallmark))) {
    fail_clear("msigdbr Hallmark table must contain columns gs_name and gene_symbol.")
  }

  hallmark <- hallmark[!is.na(hallmark$gs_name) & !is.na(hallmark$gene_symbol), , drop = FALSE]
  hallmark$gs_name <- as.character(hallmark$gs_name)
  hallmark$gene_symbol <- as.character(hallmark$gene_symbol)

  lapply(split(hallmark$gene_symbol, hallmark$gs_name), unique)
}

format_gsea_result <- function(fg) {
  if (nrow(fg) == 0) {
    return(empty_gsea_table())
  }

  fg <- as.data.frame(fg, stringsAsFactors = FALSE)

  if (!"log2err" %in% names(fg)) {
    fg$log2err <- NA_real_
  }
  if (!"nMoreExtreme" %in% names(fg)) {
    fg$nMoreExtreme <- NA_real_
  }

  if ("leadingEdge" %in% names(fg)) {
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) {
      paste(as.character(x), collapse = ";")
    }, character(1))
  } else {
    fg$leadingEdge <- NA_character_
  }

  fg$direction <- ifelse(fg$NES > 0, "HIGH_Q4", ifelse(fg$NES < 0, "LOW_Q1", "ZERO"))
  fg$interpretation <- ifelse(
    fg$NES > 0,
    "enriched in HIGH_Q4",
    ifelse(fg$NES < 0, "enriched in LOW_Q1", "no direction")
  )

  keep_cols <- c(
    "pathway", "pval", "padj", "log2err", "ES", "NES",
    "nMoreExtreme", "size", "leadingEdge", "direction", "interpretation"
  )
  missing_cols <- setdiff(keep_cols, names(fg))
  for (col in missing_cols) {
    fg[[col]] <- NA
  }

  fg <- fg[, keep_cols, drop = FALSE]
  fg[order(-fg$NES, fg$padj), , drop = FALSE]
}

top_pathway <- function(gsea_df, positive = TRUE) {
  if (nrow(gsea_df) == 0) {
    return(list(pathway = NA_character_, NES = NA_real_, padj = NA_real_))
  }

  if (positive) {
    sub <- gsea_df[gsea_df$NES > 0, , drop = FALSE]
    sub <- sub[order(-sub$NES, sub$padj), , drop = FALSE]
  } else {
    sub <- gsea_df[gsea_df$NES < 0, , drop = FALSE]
    sub <- sub[order(sub$NES, sub$padj), , drop = FALSE]
  }

  if (nrow(sub) == 0) {
    return(list(pathway = NA_character_, NES = NA_real_, padj = NA_real_))
  }

  list(pathway = sub$pathway[1], NES = sub$NES[1], padj = sub$padj[1])
}

write_qc_report <- function(path, qc) {
  sink(path)
  cat("Friedrich/GSE134051 Q1/Q4 Hallmark GSEA QC\n")
  cat("=================================\n\n")

  cat("1. Inputs\n")
  cat("  Limma no-signature input: ", qc$inputs$limma, "\n", sep = "")
  cat("  Signature input: ", qc$inputs$signature, "\n\n", sep = "")

  cat("2. Ranking\n")
  cat("  Ranking variable: t_stat\n")
  cat("  Ranking fallback: none\n")
  cat("  n_input_rows: ", qc$ranking$n_input_rows, "\n", sep = "")
  cat("  n_ranked_genes: ", qc$ranking$n_ranked_genes, "\n", sep = "")
  cat("  n_removed_nonfinite_tstat: ", qc$ranking$n_removed_nonfinite_tstat, "\n", sep = "")
  cat("  n_duplicate_tstat_values: ", qc$ranking$n_duplicate_tstat_values, "\n", sep = "")
  cat("  t_stat_min: ", qc$ranking$t_stat_min, "\n", sep = "")
  cat("  t_stat_median: ", qc$ranking$t_stat_median, "\n", sep = "")
  cat("  t_stat_max: ", qc$ranking$t_stat_max, "\n\n", sep = "")

  cat("3. Signature exclusion\n")
  cat("  Signature genes total: ", qc$signature$n_signature_genes, "\n", sep = "")
  cat("  Input is_signature_gene TRUE rows: ", qc$signature$n_input_signature_flag_true, "\n", sep = "")
  cat("  Official signature genes found in no-signature input: ",
      qc$signature$n_official_signature_genes_in_input, "\n\n", sep = "")

  cat("4. Hallmark collection\n")
  cat("  Collection: MSigDB Hallmark\n")
  cat("  n_pathways_tested: ", qc$hallmark$n_pathways_tested, "\n\n", sep = "")

  cat("5. fgsea parameters\n")
  cat("  minSize: 15\n")
  cat("  maxSize: 500\n")
  cat("  eps: 0\n\n")
  cat("  seed: ", GSEA_SEED, "\n", sep = "")
  cat("  nproc: 1\n")
  cat("  R version: ", R.version.string, "\n", sep = "")
  cat("  fgsea version: ",
      as.character(utils::packageVersion("fgsea")), "\n", sep = "")
  cat("  msigdbr version: ",
      as.character(utils::packageVersion("msigdbr")), "\n\n", sep = "")

  cat("6. GSEA results\n")
  cat("  fgsea completed: ", qc$gsea$completed, "\n", sep = "")
  cat("  n_pathways_FDR005: ", qc$gsea$n_pathways_fdr005, "\n", sep = "")
  cat("  n_positive_NES: ", qc$gsea$n_positive_nes, "\n", sep = "")
  cat("  n_negative_NES: ", qc$gsea$n_negative_nes, "\n", sep = "")
  cat("  top_positive_pathway: ", qc$gsea$top_positive_pathway, "\n", sep = "")
  cat("  top_positive_NES: ", qc$gsea$top_positive_nes, "\n", sep = "")
  cat("  top_positive_padj: ", qc$gsea$top_positive_padj, "\n", sep = "")
  cat("  top_negative_pathway: ", qc$gsea$top_negative_pathway, "\n", sep = "")
  cat("  top_negative_NES: ", qc$gsea$top_negative_nes, "\n", sep = "")
  cat("  top_negative_padj: ", qc$gsea$top_negative_padj, "\n\n", sep = "")

  cat("7. EMT pathway\n")
  cat("  EMT found: ", qc$emt$found, "\n", sep = "")
  cat("  EMT_NES: ", qc$emt$NES, "\n", sep = "")
  cat("  EMT_padj: ", qc$emt$padj, "\n", sep = "")
  cat("  EMT_direction: ", qc$emt$direction, "\n\n", sep = "")

  cat("8. Warnings\n")
  if (length(qc$warnings) == 0) {
    cat("  None\n\n")
  } else {
    for (w in qc$warnings) {
      cat("  - ", w, "\n", sep = "")
    }
    cat("\n")
  }

  cat("9. Final status: ", qc$qc_status, "\n", sep = "")
  sink()
}

make_empty_qc <- function() {
  list(
    inputs = list(limma = input_limma, signature = input_signature),
    ranking = list(
      n_input_rows = NA_integer_,
      n_ranked_genes = NA_integer_,
      n_removed_nonfinite_tstat = NA_integer_,
      n_duplicate_tstat_values = NA_integer_,
      t_stat_min = NA_real_,
      t_stat_median = NA_real_,
      t_stat_max = NA_real_
    ),
    signature = list(
      n_signature_genes = NA_integer_,
      n_input_signature_flag_true = NA_integer_,
      n_official_signature_genes_in_input = NA_integer_
    ),
    hallmark = list(n_pathways_tested = NA_integer_),
    gsea = list(
      completed = FALSE,
      n_pathways_fdr005 = NA_integer_,
      n_positive_nes = NA_integer_,
      n_negative_nes = NA_integer_,
      top_positive_pathway = NA_character_,
      top_positive_nes = NA_real_,
      top_positive_padj = NA_real_,
      top_negative_pathway = NA_character_,
      top_negative_nes = NA_real_,
      top_negative_padj = NA_real_
    ),
    emt = list(found = FALSE, NES = NA_real_, padj = NA_real_, direction = NA_character_),
    warnings = warnings_vec,
    qc_status = "FAIL"
  )
}

fail_with_qc <- function(msg, qc = make_empty_qc()) {
  qc$warnings <- warnings_vec
  qc$qc_status <- "FAIL"
  write_qc_report(output_qc, qc)
  fail_clear(msg)
}

derive_qc_status <- function(qc, full_written, metadata_written, qc_written) {
  if (is.na(qc$ranking$n_ranked_genes) || qc$ranking$n_ranked_genes < 5000) {
    return("FAIL")
  }
  if (is.na(qc$hallmark$n_pathways_tested) || qc$hallmark$n_pathways_tested < 40) {
    return("FAIL")
  }
  if (!isTRUE(qc$gsea$completed)) {
    return("FAIL")
  }
  if (!isTRUE(full_written) || !isTRUE(metadata_written) || !isTRUE(qc_written)) {
    return("FAIL")
  }
  if (qc$ranking$n_ranked_genes <= 10000) {
    return("WARNING")
  }
  if (qc$hallmark$n_pathways_tested <= 44) {
    return("WARNING")
  }
  if (!isTRUE(qc$emt$found)) {
    return("WARNING")
  }
  if (isTRUE(qc$gsea$n_pathways_fdr005 == 0)) {
    return("WARNING")
  }
  if (length(qc$warnings) > 0) {
    return("WARNING")
  }
  "PASS"
}

required_packages <- c("fgsea", "msigdbr", "jsonlite")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    fail_with_qc(paste0("Required package not installed: ", pkg))
  }
}

for (fp in c(input_limma, input_signature)) {
  if (!file.exists(fp)) {
    fail_with_qc(paste0("Missing required input: ", fp))
  }
}

limma_df <- read_csv_base(input_limma)
signature_df <- read_canonical_platelet_signature()

signature_genes <- tryCatch(
  extract_signature_genes(signature_df),
  error = function(e) {
    fail_with_qc(conditionMessage(e))
  }
)

required_cols <- c("gene_name", "t_stat", "logFC", "adj.P.Val", "is_signature_gene")
missing_cols <- setdiff(required_cols, names(limma_df))
if (length(missing_cols) > 0) {
  fail_with_qc(paste0("Limma no-signature input missing required columns: ", paste(missing_cols, collapse = ", ")))
}

limma_df$gene_name <- trimws(as.character(limma_df$gene_name))
limma_df$t_stat <- suppressWarnings(as.numeric(limma_df$t_stat))
limma_df$logFC <- suppressWarnings(as.numeric(limma_df$logFC))
limma_df$adj.P.Val <- suppressWarnings(as.numeric(limma_df$adj.P.Val))
limma_df$is_signature_gene <- coerce_signature_flag(limma_df$is_signature_gene)

if (any(is.na(limma_df$gene_name) | limma_df$gene_name == "")) {
  fail_with_qc("Limma no-signature input contains missing or empty gene_name values.")
}
if (anyDuplicated(limma_df$gene_name)) {
  duplicated_genes <- unique(limma_df$gene_name[duplicated(limma_df$gene_name)])
  fail_with_qc(paste0(
    "Duplicate gene_name values detected: ",
    paste(utils::head(duplicated_genes, 25), collapse = ", "),
    if (length(duplicated_genes) > 25) " ..." else ""
  ))
}
if (all(!is.finite(limma_df$t_stat))) {
  fail_with_qc("t_stat is present but all values are NA or non-finite.")
}
if (any(limma_df$is_signature_gene %in% TRUE, na.rm = TRUE)) {
  fail_with_qc("Limma no-signature input contains rows flagged as is_signature_gene TRUE.")
}

official_in_input <- intersect(limma_df$gene_name, signature_genes)
if (length(official_in_input) > 0) {
  qc <- make_empty_qc()
  qc$signature$n_signature_genes <- length(signature_genes)
  qc$signature$n_input_signature_flag_true <- sum(limma_df$is_signature_gene %in% TRUE, na.rm = TRUE)
  qc$signature$n_official_signature_genes_in_input <- length(official_in_input)
  fail_with_qc(paste0(
    "Limma no-signature input contains official signature genes: ",
    paste(official_in_input, collapse = ", ")
  ), qc)
}

message("[2/6] Building ranked statistic")

n_input_rows <- nrow(limma_df)
rank_df <- limma_df[is.finite(limma_df$t_stat), c("gene_name", "t_stat"), drop = FALSE]
n_removed_nonfinite <- n_input_rows - nrow(rank_df)
if (n_removed_nonfinite > 0) {
  add_warning(paste0("Removed rows with non-finite t_stat: ", n_removed_nonfinite))
}
if (nrow(rank_df) < 5000) {
  qc <- make_empty_qc()
  qc$ranking$n_input_rows <- n_input_rows
  qc$ranking$n_ranked_genes <- nrow(rank_df)
  qc$ranking$n_removed_nonfinite_tstat <- n_removed_nonfinite
  qc$signature$n_signature_genes <- length(signature_genes)
  qc$signature$n_input_signature_flag_true <- sum(limma_df$is_signature_gene %in% TRUE, na.rm = TRUE)
  qc$signature$n_official_signature_genes_in_input <- length(official_in_input)
  fail_with_qc(paste0("Too few ranked genes after t_stat filtering: ", nrow(rank_df)), qc)
}
if (nrow(rank_df) <= 10000) {
  add_warning(paste0("Ranked gene count between 5000 and 10000: ", nrow(rank_df)))
}

n_duplicate_tstat_values <- sum(duplicated(rank_df$t_stat))
ranks <- stats::setNames(rank_df$t_stat, rank_df$gene_name)
ranks <- sort(ranks, decreasing = TRUE)

ranking_qc <- list(
  n_input_rows = n_input_rows,
  n_ranked_genes = length(ranks),
  n_removed_nonfinite_tstat = n_removed_nonfinite,
  n_duplicate_tstat_values = n_duplicate_tstat_values,
  t_stat_min = min(ranks, na.rm = TRUE),
  t_stat_median = stats::median(ranks, na.rm = TRUE),
  t_stat_max = max(ranks, na.rm = TRUE)
)

message("[3/6] Loading Hallmark gene sets")

hallmark_list_raw <- tryCatch(
  load_hallmark_sets(),
  error = function(e) {
    qc <- make_empty_qc()
    qc$ranking <- ranking_qc
    qc$signature$n_signature_genes <- length(signature_genes)
    qc$signature$n_input_signature_flag_true <- sum(limma_df$is_signature_gene %in% TRUE, na.rm = TRUE)
    qc$signature$n_official_signature_genes_in_input <- length(official_in_input)
    fail_with_qc(conditionMessage(e), qc)
  }
)

hallmark_list <- lapply(hallmark_list_raw, function(gs) {
  unique(gs[gs %in% names(ranks)])
})
hallmark_list <- hallmark_list[vapply(hallmark_list, length, integer(1)) >= 15]
hallmark_list <- hallmark_list[vapply(hallmark_list, length, integer(1)) <= 500]
n_pathways_tested <- length(hallmark_list)

if (n_pathways_tested < 40) {
  qc <- make_empty_qc()
  qc$ranking <- ranking_qc
  qc$signature$n_signature_genes <- length(signature_genes)
  qc$signature$n_input_signature_flag_true <- sum(limma_df$is_signature_gene %in% TRUE, na.rm = TRUE)
  qc$signature$n_official_signature_genes_in_input <- length(official_in_input)
  qc$hallmark$n_pathways_tested <- n_pathways_tested
  fail_with_qc(paste0("Too few Hallmark pathways available after ranking overlap and size filters: ", n_pathways_tested), qc)
}
if (n_pathways_tested <= 44) {
  add_warning(paste0("Hallmark pathways tested between 40 and 44: ", n_pathways_tested))
}


message("[4/6] Running fgsea")

GSEA_SEED <- 123
set.seed(GSEA_SEED)

fg <- tryCatch(
  fgsea::fgsea(
    pathways = hallmark_list,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    nproc = 1
  ),
  error = function(e) {
    qc <- make_empty_qc()
    qc$ranking <- ranking_qc
    qc$signature$n_signature_genes <- length(signature_genes)
    qc$signature$n_input_signature_flag_true <- sum(
      limma_df$is_signature_gene %in% TRUE,
      na.rm = TRUE
    )
    qc$signature$n_official_signature_genes_in_input <- length(official_in_input)
    qc$hallmark$n_pathways_tested <- n_pathways_tested
    fail_with_qc(paste0("fgsea failed: ", conditionMessage(e)), qc)
  }
)

gsea_full <- format_gsea_result(fg)
gsea_fdr005 <- gsea_full[!is.na(gsea_full$padj) & gsea_full$padj < 0.05, , drop = FALSE]
if (nrow(gsea_fdr005) == 0) {
  add_warning("No Hallmark pathways were significant at FDR < 0.05.")
}

emt_pathway <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
emt_row <- gsea_full[gsea_full$pathway == emt_pathway, , drop = FALSE]
emt_found <- nrow(emt_row) > 0
if (!emt_found) {
  add_warning("EMT pathway was not found in Hallmark GSEA results.")
  emt_row <- empty_gsea_table()
}

top_pos <- top_pathway(gsea_full, positive = TRUE)
top_neg <- top_pathway(gsea_full, positive = FALSE)

emt_nes <- if (emt_found) emt_row$NES[1] else NA_real_
emt_padj <- if (emt_found) emt_row$padj[1] else NA_real_
emt_direction <- if (emt_found) emt_row$direction[1] else NA_character_

qc <- list(
  inputs = list(limma = input_limma, signature = input_signature),
  ranking = ranking_qc,
  signature = list(
    n_signature_genes = length(signature_genes),
    n_input_signature_flag_true = sum(limma_df$is_signature_gene %in% TRUE, na.rm = TRUE),
    n_official_signature_genes_in_input = length(official_in_input)
  ),
  hallmark = list(n_pathways_tested = n_pathways_tested),
  gsea = list(
    completed = TRUE,
    n_pathways_fdr005 = nrow(gsea_fdr005),
    n_positive_nes = sum(gsea_full$NES > 0, na.rm = TRUE),
    n_negative_nes = sum(gsea_full$NES < 0, na.rm = TRUE),
    top_positive_pathway = top_pos$pathway,
    top_positive_nes = top_pos$NES,
    top_positive_padj = top_pos$padj,
    top_negative_pathway = top_neg$pathway,
    top_negative_nes = top_neg$NES,
    top_negative_padj = top_neg$padj
  ),
  emt = list(
    found = emt_found,
    NES = emt_nes,
    padj = emt_padj,
    direction = emt_direction
  ),
  warnings = warnings_vec,
  qc_status = "PENDING"
)

metadata <- list(
  date_time = as.character(Sys.time()),
  project_dir = project_dir,
  inputs = list(
    limma_no_signature = input_limma,
    signature = input_signature
  ),
  outputs = list(
    full = output_full,
    fdr005 = output_fdr005,
    emt = output_emt,
    metadata = output_metadata,
    qc = output_qc
  ),
  method = "fgsea preranked",
  collection = "MSigDB Hallmark",
  ranking_variable = "t_stat",
  ranking_fallback = "none",
  contrast = "HIGH_Q4 - LOW_Q1",
  direction_positive = "enriched in HIGH_Q4",
  direction_negative = "enriched in LOW_Q1",
  n_input_rows = n_input_rows,
  n_ranked_genes = length(ranks),
  n_removed_nonfinite_tstat = n_removed_nonfinite,
  n_duplicate_tstat_values = n_duplicate_tstat_values,
  n_pathways_tested = n_pathways_tested,
  n_pathways_FDR005 = nrow(gsea_fdr005),
  EMT_NES = emt_nes,
  EMT_padj = emt_padj,
  EMT_direction = emt_direction,
  qc_status = "PENDING",
  reproducibility = list(
    seed = GSEA_SEED,
    nproc = 1L,
    R_version = R.version.string,
    fgsea_version = as.character(utils::packageVersion("fgsea")),
    msigdbr_version = as.character(utils::packageVersion("msigdbr"))
  )
)

message("[5/6] Writing outputs")

full_written <- FALSE
metadata_written <- FALSE
qc_written <- FALSE

tryCatch({
  write_csv_base(gsea_full, output_full)
  write_csv_base(gsea_fdr005, output_fdr005)
  write_csv_base(emt_row, output_emt)

  full_written <- all(file.exists(c(
    output_full,
    output_fdr005,
    output_emt
  )))
}, error = function(e) {
  fail_with_qc(
    paste0("Could not write GSEA CSV outputs: ", conditionMessage(e)),
    qc
  )
})

if (!full_written) {
  fail_with_qc("One or more required GSEA CSV outputs were not created.", qc)
}

# Los errores reales de escritura se controlan inmediatamente abajo.
qc$qc_status <- derive_qc_status(
  qc = qc,
  full_written = TRUE,
  metadata_written = TRUE,
  qc_written = TRUE
)
metadata$qc_status <- qc$qc_status

tryCatch({
  write_metadata_json(metadata, output_metadata)
  metadata_written <- file.exists(output_metadata)
}, error = function(e) {
  fail_clear(paste0(
    "Could not write metadata JSON: ",
    conditionMessage(e)
  ))
})

tryCatch({
  write_qc_report(output_qc, qc)
  qc_written <- file.exists(output_qc)
}, error = function(e) {
  fail_clear(paste0(
    "Could not write QC report: ",
    conditionMessage(e)
  ))
})

if (!full_written || !metadata_written || !qc_written) {
  fail_clear("One or more required GSEA outputs could not be written.")
}

if (qc$qc_status == "FAIL") {
  fail_clear("Final status: FAIL")
}

message("[6/6] Final status: ", qc$qc_status)
