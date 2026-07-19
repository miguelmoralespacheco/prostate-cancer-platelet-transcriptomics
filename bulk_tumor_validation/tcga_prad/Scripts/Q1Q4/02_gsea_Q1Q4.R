#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD | Q1 vs Q4 Hallmark GSEA
#
# Input:
#   Complete no-signature DESeq2 table from the primary Q1/Q4 analysis.
#
# Ranking:
#   DESeq2 Wald statistic from HIGH_Q4 vs LOW_Q1.
#   Multiple Ensembl rows mapping to the same gene symbol are collapsed by
#   retaining the row with the largest absolute Wald statistic.
#
# Interpretation:
#   NES > 0: enrichment in HIGH_Q4.
#   NES < 0: enrichment in LOW_Q1.
#
# GSEA:
#   MSigDB Hallmark collection.
#   fgsea, seed = 123, nproc = 1, minSize = 15, maxSize = 500, eps = 0.
#
# This script does not generate figures.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 0. Frozen analytical constants
################################################################################

EXPECTED_INPUT_ROWS <- 26966L
EXPECTED_RANKED_GENES <- 26931L
EXPECTED_HALLMARK_PATHWAYS <- 50L

GSEA_SEED <- 123L
GSEA_NPROC <- 1L
GSEA_MIN_SIZE <- 15L
GSEA_MAX_SIZE <- 500L
GSEA_EPS <- 0

PADJ_THRESHOLD <- 0.05

EMT_PATHWAY <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"

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
  "Scripts", "Q1Q4", "02_gsea_Q1Q4.R"
)

input_paths <- list(
  de_table = file.path(
    GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "DESeq2",
    "DESeq2_Q1Q4_41genes_no_signature_genes.csv"
  )
)

gene_set_dir <- file.path(
  INPUT_DIR, "GeneSets"
)

gsea_tables_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "GSEA"
)

gsea_objects_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Objects", "GSEA"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Logs"
)

invisible(lapply(
  c(gsea_tables_dir, gsea_objects_dir, logs_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

hallmark_snapshot_path <- file.path(
  gene_set_dir,
  "MSigDB_Hallmark_Homo_sapiens_frozen.csv"
)

output_paths <- list(
  rank_table = file.path(
    gsea_tables_dir,
    "GSEA_Q1Q4_DESeq2_Wald_stat_gene_symbol_rank.csv"
  ),
  full_results = file.path(
    gsea_tables_dir,
    "GSEA_Q1Q4_no_signature_Hallmark_full.csv"
  ),
  significant_results = file.path(
    gsea_tables_dir,
    "GSEA_Q1Q4_no_signature_Hallmark_FDR005.csv"
  ),
  metadata = file.path(
    gsea_tables_dir,
    "GSEA_Q1Q4_metadata.json"
  ),
  rank_rds = file.path(
    gsea_objects_dir,
    "GSEA_Q1Q4_DESeq2_Wald_stat_rank.rds"
  ),
  pathways_rds = file.path(
    gsea_objects_dir,
    "GSEA_Q1Q4_Hallmark_pathways.rds"
  ),
  results_rds = file.path(
    gsea_objects_dir,
    "GSEA_Q1Q4_Hallmark_results.rds"
  ),
  qc = file.path(
    logs_dir,
    "GSEA_Q1Q4_QC.txt"
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
    "TCGA-PRAD Q1 vs Q4 Hallmark GSEA QC",
    "======================================",
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

parse_logical_strict <- function(x, label) {
  if (is.logical(x)) {
    return(x)
  }

  value <- tolower(trimws(as.character(x)))

  output <- rep(NA, length(value))
  output[value %in% c("true", "t", "1")] <- TRUE
  output[value %in% c("false", "f", "0")] <- FALSE

  invalid <- is.na(output) & !is.na(value) & value != ""

  if (any(invalid)) {
    fail(
      "Column '", label,
      "' contains values that cannot be interpreted as logical: ",
      paste(unique(value[invalid]), collapse = ", ")
    )
  }

  output
}

clean_pathway_name <- function(pathway) {
  output <- gsub("^HALLMARK_", "", pathway)
  output <- gsub("_", " ", output)
  tools::toTitleCase(tolower(output))
}

collapse_leading_edge <- function(x) {
  vapply(
    x,
    function(genes) paste(genes, collapse = ";"),
    character(1)
  )
}

format_named_values <- function(x, n = 10L, decreasing = TRUE) {
  if (length(x) == 0) {
    return("none")
  }

  selected <- sort(x, decreasing = decreasing)
  selected <- utils::head(selected, n)

  paste(
    paste(
      names(selected),
      signif(as.numeric(selected), 6),
      sep = "\t"
    ),
    collapse = "\n"
  )
}

format_pathways <- function(df, n = 10L, decreasing = TRUE) {
  if (nrow(df) == 0) {
    return("none")
  }

  pathway_order <- order(
    df$NES,
    decreasing = decreasing,
    na.last = NA
  )

  selected <- utils::head(
    df[pathway_order, , drop = FALSE],
    n
  )

  paste(
    paste(
      selected$pathway,
      signif(selected$NES, 6),
      signif(selected$padj, 6),
      sep = "\t"
    ),
    collapse = "\n"
  )
}

################################################################################
# 3. Required packages
################################################################################

required_packages <- c(
  "fgsea",
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
# 4. Validate the DESeq2 input
################################################################################

msg("Project dir:", project_dir)
msg("[1/6] Loading the complete no-signature DESeq2 table")

de_table <- read_required_csv(
  input_paths$de_table,
  "complete no-signature DESeq2 result table"
)

required_columns <- c(
  "ensembl_id",
  "gene_name",
  "gene_name_is_ensembl_fallback",
  "stat",
  "padj",
  "is_signature_gene"
)

missing_columns <- setdiff(
  required_columns,
  names(de_table)
)

if (length(missing_columns) > 0) {
  fail(
    "DESeq2 input table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(de_table) != EXPECTED_INPUT_ROWS) {
  fail(
    "Expected ", EXPECTED_INPUT_ROWS,
    " rows in the no-signature DESeq2 table. Observed: ",
    nrow(de_table)
  )
}

de_table$ensembl_id <- trimws(as.character(de_table$ensembl_id))
de_table$gene_name <- trimws(as.character(de_table$gene_name))
de_table$stat <- suppressWarnings(as.numeric(de_table$stat))

de_table$is_signature_gene <- parse_logical_strict(
  de_table$is_signature_gene,
  "is_signature_gene"
)

de_table$gene_name_is_ensembl_fallback <- parse_logical_strict(
  de_table$gene_name_is_ensembl_fallback,
  "gene_name_is_ensembl_fallback"
)

if (
  any(is.na(de_table$ensembl_id)) ||
  any(de_table$ensembl_id == "")
) {
  fail("DESeq2 input contains missing or empty Ensembl identifiers.")
}

if (anyDuplicated(de_table$ensembl_id)) {
  fail("DESeq2 input contains duplicated Ensembl identifiers.")
}

if (
  any(is.na(de_table$gene_name)) ||
  any(de_table$gene_name == "")
) {
  fail("DESeq2 input contains missing or empty gene names.")
}

if (any(de_table$is_signature_gene %in% TRUE, na.rm = TRUE)) {
  fail(
    "The GSEA input still contains official signature genes. ",
    "Only the no-signature DESeq2 table is valid."
  )
}

if (
  any(
    de_table$gene_name_is_ensembl_fallback %in% TRUE,
    na.rm = TRUE
  )
) {
  fail(
    "The GSEA input contains Ensembl identifiers used as gene-symbol ",
    "fallbacks. Hallmark GSEA requires valid gene symbols."
  )
}

if (any(!is.finite(de_table$stat))) {
  fail(
    "The DESeq2 Wald statistic contains missing or non-finite values."
  )
}

################################################################################
# 5. Prepare the gene-symbol-level ranking
################################################################################

msg("[2/6] Preparing the DESeq2 Wald-statistic ranking")

de_table$gene_key <- toupper(de_table$gene_name)

gene_row_counts <- table(de_table$gene_key)

ranking_order <- order(
  -abs(de_table$stat),
  de_table$ensembl_id
)

rank_table <- de_table[
  ranking_order,
  c("ensembl_id", "gene_name", "gene_key", "stat"),
  drop = FALSE
]

rank_table <- rank_table[
  !duplicated(rank_table$gene_key),
  ,
  drop = FALSE
]

rank_table$n_ensembl_rows_for_gene_symbol <- as.integer(
  gene_row_counts[rank_table$gene_key]
)

duplicates_removed <- nrow(de_table) - nrow(rank_table)

ranks <- stats::setNames(
  rank_table$stat,
  rank_table$gene_key
)

ranks <- sort(
  ranks[is.finite(ranks)],
  decreasing = TRUE
)

rank_table <- rank_table[
  match(names(ranks), rank_table$gene_key),
  ,
  drop = FALSE
]

rank_table$rank_position <- seq_len(nrow(rank_table))

rank_table <- rank_table[
  ,
  c(
    "rank_position",
    "gene_key",
    "gene_name",
    "ensembl_id",
    "stat",
    "n_ensembl_rows_for_gene_symbol"
  ),
  drop = FALSE
]

if (length(ranks) != EXPECTED_RANKED_GENES) {
  fail(
    "Expected ", EXPECTED_RANKED_GENES,
    " unique ranked gene symbols. Observed: ",
    length(ranks)
  )
}

if (anyDuplicated(names(ranks))) {
  fail("The final GSEA ranking contains duplicated gene symbols.")
}

if (!all(diff(ranks) <= 0)) {
  fail("The final GSEA ranking is not sorted in decreasing order.")
}

n_tied_stat_values <- sum(duplicated(unname(ranks)))

if (n_tied_stat_values > 0) {
  add_warning(
    paste0(
      "The final ranking contains ",
      n_tied_stat_values,
      " tied Wald-statistic values."
    )
  )
}

write_csv(rank_table, output_paths$rank_table)
saveRDS(ranks, output_paths$rank_rds)

msg("Input DESeq2 rows:", nrow(de_table))
msg("Unique ranked gene symbols:", length(ranks))
msg("Duplicate gene-symbol rows removed:", duplicates_removed)

################################################################################
# 6. Load or create the frozen Hallmark collection
################################################################################

if (!file.exists(hallmark_snapshot_path)) {
  fail(
    "Required frozen MSigDB Hallmark snapshot is absent: ",
    hallmark_snapshot_path,
    ". This workflow does not replace it with a current msigdbr download."
  )
}

hallmark_table <- read_required_csv(
  hallmark_snapshot_path,
  "frozen MSigDB Hallmark snapshot"
)

hallmark_source_mode <- "existing_frozen_snapshot"

required_hallmark_columns <- c(
  "gs_name",
  "gene_symbol"
)

missing_hallmark_columns <- setdiff(
  required_hallmark_columns,
  names(hallmark_table)
)

if (length(missing_hallmark_columns) > 0) {
  fail(
    "Frozen Hallmark snapshot is missing required columns: ",
    paste(missing_hallmark_columns, collapse = ", ")
  )
}

hallmark_table$gs_name <- trimws(
  as.character(hallmark_table$gs_name)
)

hallmark_table$gene_symbol <- toupper(
  trimws(as.character(hallmark_table$gene_symbol))
)

invalid_hallmark_rows <- (
  is.na(hallmark_table$gs_name) |
    hallmark_table$gs_name == "" |
    is.na(hallmark_table$gene_symbol) |
    hallmark_table$gene_symbol == ""
)

if (any(invalid_hallmark_rows)) {
  fail(
    "Frozen Hallmark snapshot contains missing pathway or gene-symbol values."
  )
}

hallmark_table <- unique(hallmark_table)

pathways <- split(
  hallmark_table$gene_symbol,
  hallmark_table$gs_name
)

pathways <- lapply(
  pathways,
  function(x) sort(unique(x))
)

pathways <- pathways[
  sort(names(pathways))
]

if (length(pathways) != EXPECTED_HALLMARK_PATHWAYS) {
  fail(
    "Expected ", EXPECTED_HALLMARK_PATHWAYS,
    " Hallmark pathways. Observed: ",
    length(pathways)
  )
}

if (!EMT_PATHWAY %in% names(pathways)) {
  fail(
    "The frozen Hallmark collection does not contain ",
    EMT_PATHWAY,
    "."
  )
}

hallmark_snapshot_md5 <- unname(
  tools::md5sum(hallmark_snapshot_path)
)

hallmark_db_versions <- if (
  "db_version" %in% names(hallmark_table)
) {
  sort(unique(as.character(hallmark_table$db_version)))
} else {
  character(0)
}

saveRDS(pathways, output_paths$pathways_rds)

msg("Hallmark source mode:", hallmark_source_mode)
msg("Hallmark pathways:", length(pathways))
msg("Hallmark snapshot MD5:", hallmark_snapshot_md5)

################################################################################
# 7. Run deterministic Hallmark GSEA
################################################################################

msg("[4/6] Running deterministic Hallmark GSEA")

set.seed(GSEA_SEED)

fgsea_result <- fgsea::fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = GSEA_MIN_SIZE,
  maxSize = GSEA_MAX_SIZE,
  eps = GSEA_EPS,
  nproc = GSEA_NPROC
)

fgsea_result <- as.data.frame(fgsea_result)

if (nrow(fgsea_result) != EXPECTED_HALLMARK_PATHWAYS) {
  fail(
    "Expected results for ", EXPECTED_HALLMARK_PATHWAYS,
    " Hallmark pathways. Observed: ",
    nrow(fgsea_result)
  )
}

if (anyDuplicated(fgsea_result$pathway)) {
  fail("GSEA output contains duplicated pathway names.")
}

if (
  any(!is.finite(fgsea_result$NES)) ||
  any(!is.finite(fgsea_result$ES))
) {
  fail("GSEA output contains non-finite ES or NES values.")
}

fgsea_result$direction <- ifelse(
  fgsea_result$NES > 0,
  "enriched_in_HIGH_Q4",
  "enriched_in_LOW_Q1"
)

fgsea_result$pathway_clean <- clean_pathway_name(
  fgsea_result$pathway
)

fgsea_result <- fgsea_result[
  order(
    fgsea_result$padj,
    -abs(fgsea_result$NES),
    fgsea_result$pathway,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

significant_result <- fgsea_result[
  !is.na(fgsea_result$padj) &
    fgsea_result$padj < PADJ_THRESHOLD,
  ,
  drop = FALSE
]

emt_result <- fgsea_result[
  fgsea_result$pathway == EMT_PATHWAY,
  ,
  drop = FALSE
]

if (nrow(emt_result) != 1L) {
  fail(
    "Expected exactly one EMT Hallmark result. Observed: ",
    nrow(emt_result)
  )
}

saveRDS(
  fgsea_result,
  output_paths$results_rds
)

fgsea_csv <- fgsea_result
fgsea_csv$leadingEdge <- collapse_leading_edge(
  fgsea_csv$leadingEdge
)

significant_csv <- significant_result
significant_csv$leadingEdge <- collapse_leading_edge(
  significant_csv$leadingEdge
)

write_csv(
  fgsea_csv,
  output_paths$full_results
)

write_csv(
  significant_csv,
  output_paths$significant_results
)

msg("Pathways tested:", nrow(fgsea_result))
msg("Significant Hallmarks:", nrow(significant_result))
msg(
  "EMT:",
  paste0(
    "NES=", signif(emt_result$NES, 7),
    "; padj=", signif(emt_result$padj, 7),
    "; direction=", emt_result$direction
  )
)

################################################################################
# 8. Metadata
################################################################################

msg("[5/6] Writing metadata")

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
  analysis = "Q1 vs Q4 Hallmark GSEA",
  analytic_hierarchy = "primary transcriptomic analysis",
  input = list(
    DESeq2_table = input_paths$de_table,
    no_signature_table = TRUE,
    input_rows = nrow(de_table),
    signature_genes_present = FALSE,
    Ensembl_fallback_rows = sum(
      de_table$gene_name_is_ensembl_fallback,
      na.rm = TRUE
    )
  ),
  ranking = list(
    metric = "DESeq2 Wald statistic",
    contrast = "HIGH_Q4 vs LOW_Q1",
    column = "stat",
    direction = paste0(
      "positive statistic and NES indicate HIGH_Q4; ",
      "negative statistic and NES indicate LOW_Q1"
    ),
    input_rows = nrow(de_table),
    ranked_gene_symbols = length(ranks),
    duplicate_gene_symbol_rows_removed = duplicates_removed,
    duplicate_rule = paste0(
      "retain the Ensembl row with the largest absolute Wald statistic; ",
      "use Ensembl ID as deterministic tie-breaker"
    ),
    tied_stat_values = n_tied_stat_values
  ),
  gene_sets = list(
    database = "MSigDB Hallmark collection",
    species = "Homo sapiens",
    snapshot_file = hallmark_snapshot_path,
    snapshot_source_mode = hallmark_source_mode,
    snapshot_md5 = hallmark_snapshot_md5,
    db_versions = hallmark_db_versions,
    pathways = length(pathways)
  ),
  GSEA = list(
    function_name = "fgsea::fgsea",
    seed = GSEA_SEED,
    nproc = GSEA_NPROC,
    minSize = GSEA_MIN_SIZE,
    maxSize = GSEA_MAX_SIZE,
    eps = GSEA_EPS,
    pathways_tested = nrow(fgsea_result),
    pathways_FDR005 = nrow(significant_result)
  ),
  EMT = list(
    pathway = EMT_PATHWAY,
    ES = unname(emt_result$ES),
    NES = unname(emt_result$NES),
    pvalue = unname(emt_result$pval),
    padj = unname(emt_result$padj),
    size = unname(emt_result$size),
    direction = unname(emt_result$direction)
  ),
  outputs = output_paths,
  reproducibility = list(
    R_version = R.version.string,
    fgsea_version = as.character(
      utils::packageVersion("fgsea")
    ),
    msigdbr_version = "not_used_frozen_snapshot",
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
# 9. Final QC report
################################################################################

msg("[6/6] Writing final QC report")

qc_lines <- c(
  "TCGA-PRAD Q1 vs Q4 Hallmark GSEA QC",
  "======================================",
  "",
  paste0("Date/time: ", as.character(Sys.time())),
  paste0("Project dir: ", project_dir),
  paste0("Script: ", script_relative_path),
  "",
  "1. Input",
  paste0("  DESeq2 table: ", input_paths$de_table),
  paste0("  Input rows: ", nrow(de_table)),
  "  Official signature genes present: FALSE",
  paste0(
    "  Ensembl fallback rows: ",
    sum(de_table$gene_name_is_ensembl_fallback, na.rm = TRUE)
  ),
  "",
  "2. Ranking",
  "  Metric: DESeq2 Wald statistic",
  "  Contrast: HIGH_Q4 vs LOW_Q1",
  paste0("  Input Ensembl rows: ", nrow(de_table)),
  paste0("  Ranked gene symbols: ", length(ranks)),
  paste0(
    "  Duplicate gene-symbol rows removed: ",
    duplicates_removed
  ),
  paste0("  Tied statistic values: ", n_tied_stat_values),
  paste0(
    "  Duplicate rule: largest absolute Wald statistic; ",
    "Ensembl ID as deterministic tie-breaker"
  ),
  "",
  "3. Hallmark gene sets",
  paste0("  Snapshot: ", hallmark_snapshot_path),
  paste0("  Snapshot source mode: ", hallmark_source_mode),
  paste0("  Snapshot MD5: ", hallmark_snapshot_md5),
  paste0(
    "  MSigDB database version: ",
    if (length(hallmark_db_versions) == 0) {
      "not reported in snapshot"
    } else {
      paste(hallmark_db_versions, collapse = ", ")
    }
  ),
  paste0("  Hallmark pathways: ", length(pathways)),
  "",
  "4. GSEA",
  "  Function: fgsea::fgsea",
  paste0("  Seed: ", GSEA_SEED),
  paste0("  nproc: ", GSEA_NPROC),
  paste0("  minSize: ", GSEA_MIN_SIZE),
  paste0("  maxSize: ", GSEA_MAX_SIZE),
  paste0("  eps: ", GSEA_EPS),
  paste0("  Pathways tested: ", nrow(fgsea_result)),
  paste0(
    "  Significant pathways, FDR < ",
    PADJ_THRESHOLD,
    ": ",
    nrow(significant_result)
  ),
  "  NES > 0: enriched in HIGH_Q4",
  "  NES < 0: enriched in LOW_Q1",
  "",
  "5. EMT",
  paste0("  Pathway: ", EMT_PATHWAY),
  paste0("  ES: ", signif(emt_result$ES, 8)),
  paste0("  NES: ", signif(emt_result$NES, 8)),
  paste0("  p-value: ", signif(emt_result$pval, 8)),
  paste0("  FDR: ", signif(emt_result$padj, 8)),
  paste0("  Size: ", emt_result$size),
  paste0("  Direction: ", emt_result$direction),
  "",
  "6. Top positive ranked genes",
  format_named_values(
    ranks,
    n = 10L,
    decreasing = TRUE
  ),
  "",
  "7. Top negative ranked genes",
  format_named_values(
    ranks,
    n = 10L,
    decreasing = FALSE
  ),
  "",
  "8. Top positive NES pathways",
  format_pathways(
    fgsea_result,
    n = 10L,
    decreasing = TRUE
  ),
  "",
  "9. Top negative NES pathways",
  format_pathways(
    fgsea_result,
    n = 10L,
    decreasing = FALSE
  ),
  "",
  "10. Outputs",
  paste0("  Gene-symbol rank table: ", output_paths$rank_table),
  paste0("  Full Hallmark results: ", output_paths$full_results),
  paste0(
    "  Significant Hallmark results: ",
    output_paths$significant_results
  ),
  paste0("  Rank RDS: ", output_paths$rank_rds),
  paste0("  Pathway RDS: ", output_paths$pathways_rds),
  paste0("  GSEA result RDS: ", output_paths$results_rds),
  paste0("  Metadata: ", output_paths$metadata),
  paste0("  QC: ", output_paths$qc),
  "",
  "11. Figures",
  "  None",
  "",
  "12. Warnings",
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

required_outputs <- c(
  hallmark_snapshot_path,
  unlist(output_paths, use.names = FALSE)
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

msg("Frozen Hallmark snapshot:", hallmark_snapshot_path)
msg("Rank table:", output_paths$rank_table)
msg("Full GSEA results:", output_paths$full_results)
msg("Significant GSEA results:", output_paths$significant_results)
msg("Metadata:", output_paths$metadata)
msg("QC:", output_paths$qc)
msg("Final status:", final_status)
