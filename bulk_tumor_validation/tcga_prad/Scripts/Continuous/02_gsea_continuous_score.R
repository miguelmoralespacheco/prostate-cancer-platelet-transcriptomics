#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD | Continuous platelet-associated score Hallmark GSEA
#
# Input:
#   Complete no-signature DESeq2 continuous-model table.
#
# Ranking:
#   Unshrunken DESeq2 Wald statistic (stat).
#
# Duplicate gene-symbol rule:
#   Retain the Ensembl row with the largest absolute Wald statistic.
#   Use Ensembl ID as deterministic tie-breaker.
#
# Gene sets:
#   Exact frozen MSigDB Hallmark snapshot already used by Q1/Q4.
#
# Reproducibility:
#   seed = 123
#   nproc = 1
#   minSize = 15
#   maxSize = 500
#   eps = 0
#
# This script does not:
#   - rerun DESeq2;
#   - use the thresholded score-associated-gene table;
#   - query msigdbr;
#   - cap, transform, or jitter the Wald statistics;
#   - compare against legacy analyses;
#   - generate figures.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 0. Frozen analytical expectations
################################################################################

EXPECTED_INPUT_ROWS <- 28460L
EXPECTED_RANKED_GENES <- 28406L
EXPECTED_DUPLICATE_ROWS_REMOVED <- 54L
EXPECTED_HALLMARK_PATHWAYS <- 50L

EXPECTED_HALLMARK_MD5 <- "56ef50e187bafc696128554ec7406702"

GSEA_SEED <- 123L
GSEA_NPROC <- 1L
GSEA_MIN_SIZE <- 15L
GSEA_MAX_SIZE <- 500L
GSEA_EPS <- 0

FDR_THRESHOLD <- 0.05

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
  "Scripts",
  "Continuous",
  "02_gsea_continuous_score.R"
)

input_paths <- list(
  continuous_no_signature = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Tables",
    "DESeq2",
    "DESeq2_continuous_score_no_signature_genes.csv"
  ),
  hallmark_snapshot = file.path(
    INPUT_DIR, "GeneSets",
    "MSigDB_Hallmark_Homo_sapiens_frozen.csv"
  )
)

tables_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables",
  "GSEA"
)

objects_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Objects",
  "GSEA"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Logs"
)

invisible(lapply(
  c(tables_dir, objects_dir, logs_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

output_paths <- list(
  rank_table = file.path(
    tables_dir,
    "GSEA_continuous_DESeq2_Wald_stat_gene_symbol_rank.csv"
  ),
  full_results = file.path(
    tables_dir,
    "GSEA_continuous_no_signature_Hallmark_full.csv"
  ),
  significant_results = file.path(
    tables_dir,
    "GSEA_continuous_no_signature_Hallmark_FDR005.csv"
  ),
  metadata = file.path(
    tables_dir,
    "GSEA_continuous_metadata.json"
  ),
  rank_rds = file.path(
    objects_dir,
    "GSEA_continuous_DESeq2_Wald_stat_rank.rds"
  ),
  pathways_rds = file.path(
    objects_dir,
    "GSEA_continuous_Hallmark_pathways.rds"
  ),
  results_rds = file.path(
    objects_dir,
    "GSEA_continuous_Hallmark_results.rds"
  ),
  qc = file.path(
    logs_dir,
    "GSEA_continuous_QC.txt"
  )
)

################################################################################
# 2. Helpers
################################################################################

msg <- function(...) {
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "|",
    ...,
    "\n"
  )
}

fail_with_qc <- function(message_text) {
  failure_lines <- c(
    "TCGA-PRAD continuous Hallmark GSEA QC",
    "=======================================",
    "",
    paste0("Date/time: ", as.character(Sys.time())),
    paste0("Project dir: ", project_dir),
    paste0("Script: ", script_relative_path),
    paste0("Failure reason: ", message_text),
    "",
    "Final status: FAIL"
  )

  writeLines(
    failure_lines,
    output_paths$qc,
    useBytes = TRUE
  )

  stop(message_text, call. = FALSE)
}

fail <- function(...) {
  fail_with_qc(paste0(...))
}

require_input <- function(path, label) {
  if (!file.exists(path)) {
    fail(
      "Missing required input: ",
      label,
      " at ",
      path
    )
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
    output <- x
  } else {
    value <- tolower(trimws(as.character(x)))
    output <- rep(NA, length(value))

    output[value %in% c("true", "t", "1")] <- TRUE
    output[value %in% c("false", "f", "0")] <- FALSE

    invalid <- is.na(output) &
      !is.na(value) &
      value != ""

    if (any(invalid)) {
      fail(
        "Column '",
        label,
        "' contains invalid logical values: ",
        paste(unique(value[invalid]), collapse = ", ")
      )
    }
  }

  if (any(is.na(output))) {
    fail(
      "Column '",
      label,
      "' contains missing logical values."
    )
  }

  output
}

clean_hallmark_label <- function(pathway) {
  pathway_clean <- gsub(
    "^HALLMARK_",
    "",
    pathway
  )

  pathway_clean <- gsub(
    "_",
    " ",
    pathway_clean
  )

  pathway_clean <- tools::toTitleCase(
    tolower(pathway_clean)
  )

  pathway_clean <- gsub(
    "Tnfa",
    "TNF-alpha",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Nfkb",
    "NFkB",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Il6 Jak Stat3",
    "IL-6/JAK/STAT3",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Il2 Stat5",
    "IL2/STAT5",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Kras",
    "KRAS",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Tgf Beta",
    "TGF-beta",
    pathway_clean
  )

  pathway_clean <- gsub(
    "P53",
    "p53",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Myc",
    "MYC",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Signaling",
    "signaling",
    pathway_clean
  )

  pathway_clean <- gsub(
    "Via",
    "via",
    pathway_clean
  )

  pathway_clean
}

format_scientific <- function(x, digits = 4) {
  if (length(x) == 0 || is.na(x)) {
    return("NA")
  }

  format(
    x,
    scientific = TRUE,
    digits = digits,
    trim = TRUE
  )
}

format_top_pathways <- function(
    df,
    decreasing = TRUE,
    n = 10L
) {
  if (nrow(df) == 0) {
    return("none")
  }

  ord <- order(
    df$NES,
    decreasing = decreasing,
    na.last = NA
  )

  selected <- utils::head(
    df[ord, , drop = FALSE],
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

validate_output <- function(path, label) {
  if (!file.exists(path)) {
    fail(
      label,
      " was not generated: ",
      path
    )
  }

  file_size <- file.info(path)$size

  if (is.na(file_size) || file_size <= 0) {
    fail(
      label,
      " is empty or invalid: ",
      path
    )
  }
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
# 4. Validate inputs
################################################################################

msg("Project dir:", project_dir)

for (nm in names(input_paths)) {
  require_input(input_paths[[nm]], nm)
}

################################################################################
# 5. Load and validate the continuous DESeq2 table
################################################################################

msg("[1/7] Loading continuous no-signature DESeq2 table")

de <- read_required_csv(
  input_paths$continuous_no_signature,
  "continuous no-signature DESeq2 table"
)

required_de_columns <- c(
  "ensembl_id",
  "gene_name",
  "gene_name_is_ensembl_fallback",
  "stat",
  "is_signature_gene"
)

missing_de_columns <- setdiff(
  required_de_columns,
  names(de)
)

if (length(missing_de_columns) > 0) {
  fail(
    "Continuous DESeq2 table is missing required columns: ",
    paste(missing_de_columns, collapse = ", ")
  )
}

if (nrow(de) != EXPECTED_INPUT_ROWS) {
  fail(
    "Expected ",
    EXPECTED_INPUT_ROWS,
    " rows in the continuous no-signature input. Observed: ",
    nrow(de)
  )
}

de$ensembl_id <- trimws(
  as.character(de$ensembl_id)
)

de$gene_name <- trimws(
  as.character(de$gene_name)
)

de$stat <- suppressWarnings(
  as.numeric(de$stat)
)

de$is_signature_gene <- parse_logical_strict(
  de$is_signature_gene,
  "is_signature_gene"
)

de$gene_name_is_ensembl_fallback <- parse_logical_strict(
  de$gene_name_is_ensembl_fallback,
  "gene_name_is_ensembl_fallback"
)

if (
  any(is.na(de$ensembl_id)) ||
  any(de$ensembl_id == "")
) {
  fail(
    "Continuous DESeq2 input contains missing or empty Ensembl IDs."
  )
}

if (
  any(is.na(de$gene_name)) ||
  any(de$gene_name == "")
) {
  fail(
    "Continuous DESeq2 input contains missing or empty gene names."
  )
}

if (anyDuplicated(de$ensembl_id)) {
  fail(
    "Continuous DESeq2 input contains duplicated Ensembl IDs."
  )
}

if (any(de$is_signature_gene)) {
  fail(
    "Continuous no-signature input still contains signature genes."
  )
}

if (any(de$gene_name_is_ensembl_fallback)) {
  fail(
    "Continuous input contains Ensembl fallback gene names."
  )
}

if (any(!is.finite(de$stat))) {
  fail(
    "Continuous DESeq2 input contains non-finite Wald statistics."
  )
}

msg("Input rows:", nrow(de))
msg("Signature rows detected:", sum(de$is_signature_gene))
msg(
  "Ensembl fallback rows:",
  sum(de$gene_name_is_ensembl_fallback)
)

################################################################################
# 6. Build the deterministic gene-symbol ranking
################################################################################

msg("[2/7] Building deterministic Wald-statistic ranking")

rank_input <- data.frame(
  ensembl_id = de$ensembl_id,
  gene_symbol = toupper(de$gene_name),
  stat = de$stat,
  stringsAsFactors = FALSE
)

rank_input <- rank_input[
  order(
    rank_input$gene_symbol,
    -abs(rank_input$stat),
    rank_input$ensembl_id
  ),
  ,
  drop = FALSE
]

rank_df <- rank_input[
  !duplicated(rank_input$gene_symbol),
  ,
  drop = FALSE
]

duplicate_rows_removed <- nrow(rank_input) - nrow(rank_df)

if (
  duplicate_rows_removed !=
  EXPECTED_DUPLICATE_ROWS_REMOVED
) {
  fail(
    "Expected ",
    EXPECTED_DUPLICATE_ROWS_REMOVED,
    " duplicate gene-symbol rows to be removed. Observed: ",
    duplicate_rows_removed
  )
}

if (nrow(rank_df) != EXPECTED_RANKED_GENES) {
  fail(
    "Expected ",
    EXPECTED_RANKED_GENES,
    " unique ranked genes. Observed: ",
    nrow(rank_df)
  )
}

rank_df <- rank_df[
  order(
    -rank_df$stat,
    rank_df$gene_symbol
  ),
  ,
  drop = FALSE
]

rank_df$rank <- seq_len(nrow(rank_df))

rank_df <- rank_df[
  ,
  c(
    "rank",
    "gene_symbol",
    "ensembl_id",
    "stat"
  ),
  drop = FALSE
]

ranks <- stats::setNames(
  rank_df$stat,
  rank_df$gene_symbol
)

if (
  anyDuplicated(names(ranks)) ||
  any(!is.finite(ranks))
) {
  fail(
    "Final gene-level ranking contains duplicated names or non-finite values."
  )
}

tied_stat_rows <- sum(
  duplicated(rank_df$stat)
)

write_csv(
  rank_df,
  output_paths$rank_table
)

saveRDS(
  ranks,
  output_paths$rank_rds
)

msg("Valid input rows:", nrow(rank_input))
msg("Ranked unique genes:", length(ranks))
msg("Duplicate rows removed:", duplicate_rows_removed)
msg("Repeated statistic values beyond first occurrence:", tied_stat_rows)

################################################################################
# 7. Load and validate the frozen Hallmark snapshot
################################################################################

msg("[3/7] Loading frozen Hallmark snapshot")

hallmark_md5 <- unname(
  tools::md5sum(
    input_paths$hallmark_snapshot
  )
)

if (!identical(hallmark_md5, EXPECTED_HALLMARK_MD5)) {
  fail(
    "Frozen Hallmark snapshot MD5 differs from the Q1/Q4 contract. ",
    "Expected: ",
    EXPECTED_HALLMARK_MD5,
    "; observed: ",
    hallmark_md5
  )
}

hallmark <- read_required_csv(
  input_paths$hallmark_snapshot,
  "frozen Hallmark snapshot"
)

required_hallmark_columns <- c(
  "gs_name",
  "gene_symbol"
)

missing_hallmark_columns <- setdiff(
  required_hallmark_columns,
  names(hallmark)
)

if (length(missing_hallmark_columns) > 0) {
  fail(
    "Frozen Hallmark snapshot is missing required columns: ",
    paste(missing_hallmark_columns, collapse = ", ")
  )
}

hallmark$gs_name <- trimws(
  as.character(hallmark$gs_name)
)

hallmark$gene_symbol <- toupper(
  trimws(
    as.character(hallmark$gene_symbol)
  )
)

hallmark <- hallmark[
  !is.na(hallmark$gs_name) &
    hallmark$gs_name != "" &
    !is.na(hallmark$gene_symbol) &
    hallmark$gene_symbol != "",
  ,
  drop = FALSE
]

hallmark_pairs <- unique(
  hallmark[
    ,
    c("gs_name", "gene_symbol"),
    drop = FALSE
  ]
)

pathways <- split(
  hallmark_pairs$gene_symbol,
  hallmark_pairs$gs_name
)

pathways <- lapply(
  pathways,
  unique
)

pathways <- pathways[
  order(names(pathways))
]

if (length(pathways) != EXPECTED_HALLMARK_PATHWAYS) {
  fail(
    "Expected ",
    EXPECTED_HALLMARK_PATHWAYS,
    " Hallmark pathways. Observed: ",
    length(pathways)
  )
}

if (!EMT_PATHWAY %in% names(pathways)) {
  fail(
    "Frozen Hallmark snapshot does not contain ",
    EMT_PATHWAY,
    "."
  )
}

hallmark_db_versions <- character(0)

if ("db_version" %in% names(hallmark)) {
  hallmark_db_versions <- unique(
    trimws(
      as.character(hallmark$db_version)
    )
  )

  hallmark_db_versions <- hallmark_db_versions[
    !is.na(hallmark_db_versions) &
      hallmark_db_versions != ""
  ]
}

saveRDS(
  pathways,
  output_paths$pathways_rds
)

msg("Hallmark snapshot MD5:", hallmark_md5)
msg("Hallmark pathways:", length(pathways))

################################################################################
# 8. Run deterministic Hallmark GSEA
################################################################################

msg("[4/7] Running deterministic Hallmark GSEA")

set.seed(GSEA_SEED)

warning_messages <- character(0)

set.seed(GSEA_SEED)

fgsea_raw <- withCallingHandlers(
  fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = GSEA_MIN_SIZE,
    maxSize = GSEA_MAX_SIZE,
    eps = GSEA_EPS,
    nproc = GSEA_NPROC
  ),
  warning = function(w) {
    warning_messages <<- unique(
      c(
        warning_messages,
        conditionMessage(w)
      )
    )

    invokeRestart("muffleWarning")
  }
)

fgsea_results <- as.data.frame(
  fgsea_raw
)

required_result_columns <- c(
  "pathway",
  "pval",
  "padj",
  "log2err",
  "ES",
  "NES",
  "size",
  "leadingEdge"
)

missing_result_columns <- setdiff(
  required_result_columns,
  names(fgsea_results)
)

if (length(missing_result_columns) > 0) {
  fail(
    "fgsea output is missing required columns: ",
    paste(missing_result_columns, collapse = ", ")
  )
}

fgsea_results$direction <- ifelse(
  fgsea_results$NES > 0,
  "enriched_with_higher_Score_z",
  "enriched_with_lower_Score_z"
)

fgsea_results$pathway_clean <- clean_hallmark_label(
  fgsea_results$pathway
)

fgsea_results <- fgsea_results[
  order(
    fgsea_results$padj,
    -abs(fgsea_results$NES),
    fgsea_results$pathway,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

if (nrow(fgsea_results) != EXPECTED_HALLMARK_PATHWAYS) {
  fail(
    "Expected ",
    EXPECTED_HALLMARK_PATHWAYS,
    " Hallmark pathways to be tested. Observed: ",
    nrow(fgsea_results)
  )
}

fgsea_significant <- fgsea_results[
  !is.na(fgsea_results$padj) &
    fgsea_results$padj < FDR_THRESHOLD,
  ,
  drop = FALSE
]

emt_result <- fgsea_results[
  fgsea_results$pathway == EMT_PATHWAY,
  ,
  drop = FALSE
]

if (nrow(emt_result) != 1L) {
  fail(
    "Expected exactly one EMT result. Observed: ",
    nrow(emt_result)
  )
}

saveRDS(
  fgsea_results,
  output_paths$results_rds
)

################################################################################
# 9. Write GSEA tables
################################################################################

msg("[5/7] Writing Hallmark GSEA tables")

fgsea_full_csv <- fgsea_results

fgsea_full_csv$leadingEdge <- vapply(
  fgsea_full_csv$leadingEdge,
  function(x) {
    paste(
      as.character(x),
      collapse = ";"
    )
  },
  character(1)
)

fgsea_significant_csv <- fgsea_significant

fgsea_significant_csv$leadingEdge <- vapply(
  fgsea_significant_csv$leadingEdge,
  function(x) {
    paste(
      as.character(x),
      collapse = ";"
    )
  },
  character(1)
)

write_csv(
  fgsea_full_csv,
  output_paths$full_results
)

write_csv(
  fgsea_significant_csv,
  output_paths$significant_results
)

msg(
  "Significant Hallmarks at FDR < 0.05:",
  nrow(fgsea_significant)
)

msg(
  "EMT NES:",
  signif(emt_result$NES, 8),
  "| EMT FDR:",
  format_scientific(emt_result$padj, digits = 6)
)

################################################################################
# 10. Metadata
################################################################################

msg("[6/7] Writing metadata")

final_status <- if (length(warning_messages) == 0) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

metadata <- list(
  date_time = as.character(Sys.time()),
  project_dir = project_dir,
  script = script_relative_path,
  project = "tcga_prad",
  analysis = "continuous platelet-associated score Hallmark GSEA",
  analytic_hierarchy = "complementary robustness analysis",
  input = list(
    DESeq2_table = input_paths$continuous_no_signature,
    no_signature_table = TRUE,
    input_rows = nrow(de),
    signature_genes_present = any(de$is_signature_gene),
    Ensembl_fallback_rows = sum(
      de$gene_name_is_ensembl_fallback
    )
  ),
  ranking = list(
    metric = "DESeq2 Wald statistic",
    model = "~ Score_z",
    coefficient = "Score_z",
    coefficient_unit = "one standard deviation increase in Score_z",
    column = "stat",
    direction = paste0(
      "positive statistic and NES indicate association with higher Score_z; ",
      "negative statistic and NES indicate association with lower Score_z"
    ),
    input_rows = nrow(rank_input),
    ranked_gene_symbols = length(ranks),
    duplicate_gene_symbol_rows_removed = duplicate_rows_removed,
    duplicate_rule = paste0(
      "retain the Ensembl row with the largest absolute Wald statistic; ",
      "use Ensembl ID as deterministic tie-breaker"
    ),
    repeated_statistic_values_beyond_first = tied_stat_rows,
    ranking_transformed = FALSE,
    ranking_capped = FALSE,
    ranking_jittered = FALSE
  ),
  gene_sets = list(
    database = "MSigDB Hallmark collection",
    species = "Homo sapiens",
    snapshot_file = input_paths$hallmark_snapshot,
    snapshot_source_mode = "existing_frozen_snapshot",
    snapshot_md5 = hallmark_md5,
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
    pathways_tested = nrow(fgsea_results),
    pathways_FDR005 = nrow(fgsea_significant)
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
  interpretation = paste0(
    "This is an unadjusted continuous molecular association analysis. ",
    "Positive NES indicates enrichment among genes associated with higher ",
    "Score_z; negative NES indicates enrichment among genes associated with ",
    "lower Score_z. Results are not causal and are not independent of TME ",
    "composition."
  ),
  outputs = output_paths,
  reproducibility = list(
    R_version = R.version.string,
    fgsea_version = as.character(
      utils::packageVersion("fgsea")
    ),
    jsonlite_version = as.character(
      utils::packageVersion("jsonlite")
    )
  ),
  warnings = warning_messages,
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
# 11. Final QC report
################################################################################

msg("[7/7] Writing final QC report")

qc_lines <- c(
  "TCGA-PRAD continuous Hallmark GSEA QC",
  "=======================================",
  "",
  paste0("Date/time: ", as.character(Sys.time())),
  paste0("Project dir: ", project_dir),
  paste0("Script: ", script_relative_path),
  "",
  "1. Analytical hierarchy",
  "  Continuous model status: complementary robustness analysis",
  "  Primary transcriptomic analysis: Q1 vs Q4",
  "",
  "2. Input",
  paste0(
    "  Continuous no-signature table: ",
    input_paths$continuous_no_signature
  ),
  paste0("  Input rows: ", nrow(de)),
  paste0(
    "  Signature genes detected: ",
    sum(de$is_signature_gene)
  ),
  paste0(
    "  Ensembl fallback rows: ",
    sum(de$gene_name_is_ensembl_fallback)
  ),
  "",
  "3. Ranking",
  "  Metric: unshrunken DESeq2 Wald statistic",
  "  Model: ~ Score_z",
  "  Coefficient unit: one standard deviation increase in Score_z",
  paste0("  Valid input rows: ", nrow(rank_input)),
  paste0("  Unique ranked genes: ", length(ranks)),
  paste0(
    "  Duplicate symbol rows removed: ",
    duplicate_rows_removed
  ),
  paste0(
    "  Repeated statistic values beyond first occurrence: ",
    tied_stat_rows
  ),
  "  Duplicate rule: largest absolute Wald statistic; Ensembl ID tie-break",
  "  Ranking capped: FALSE",
  "  Ranking jittered: FALSE",
  "  Ranking otherwise transformed: FALSE",
  "",
  "4. Hallmark gene sets",
  paste0(
    "  Frozen snapshot: ",
    input_paths$hallmark_snapshot
  ),
  paste0("  Snapshot MD5: ", hallmark_md5),
  paste0(
    "  Snapshot database version: ",
    if (length(hallmark_db_versions) == 0) {
      "not recorded in snapshot"
    } else {
      paste(hallmark_db_versions, collapse = "; ")
    }
  ),
  paste0("  Hallmark pathways: ", length(pathways)),
  "  msigdbr queried by this script: FALSE",
  "",
  "5. GSEA",
  "  Function: fgsea::fgsea",
  paste0("  seed: ", GSEA_SEED),
  paste0("  nproc: ", GSEA_NPROC),
  paste0("  minSize: ", GSEA_MIN_SIZE),
  paste0("  maxSize: ", GSEA_MAX_SIZE),
  paste0("  eps: ", GSEA_EPS),
  paste0(
    "  Pathways tested: ",
    nrow(fgsea_results)
  ),
  paste0(
    "  Pathways with FDR < 0.05: ",
    nrow(fgsea_significant)
  ),
  "",
  "6. EMT",
  paste0("  Pathway: ", EMT_PATHWAY),
  paste0("  ES: ", signif(emt_result$ES, 10)),
  paste0("  NES: ", signif(emt_result$NES, 10)),
  paste0(
    "  p-value: ",
    format_scientific(emt_result$pval, digits = 8)
  ),
  paste0(
    "  FDR: ",
    format_scientific(emt_result$padj, digits = 8)
  ),
  paste0("  Size: ", emt_result$size),
  paste0("  Direction: ", emt_result$direction),
  "",
  "7. Top positive NES pathways",
  format_top_pathways(
    fgsea_results,
    decreasing = TRUE,
    n = 10L
  ),
  "",
  "8. Top negative NES pathways",
  format_top_pathways(
    fgsea_results,
    decreasing = FALSE,
    n = 10L
  ),
  "",
  "9. Outputs",
  paste0("  Rank table: ", output_paths$rank_table),
  paste0("  Full GSEA results: ", output_paths$full_results),
  paste0(
    "  Significant GSEA results: ",
    output_paths$significant_results
  ),
  paste0("  Metadata: ", output_paths$metadata),
  paste0("  Rank RDS: ", output_paths$rank_rds),
  paste0("  Pathways RDS: ", output_paths$pathways_rds),
  paste0("  Results RDS: ", output_paths$results_rds),
  paste0("  QC: ", output_paths$qc),
  "",
  "10. Excluded operations",
  "  DESeq2 rerun: FALSE",
  "  Thresholded 6,038-gene table used: FALSE",
  "  Legacy comparison: FALSE",
  "  Rank capping: FALSE",
  "  Rank jitter: FALSE",
  "  Figures generated: none",
  "",
  "11. Warnings",
  if (length(warning_messages) == 0) {
    "  None"
  } else {
    paste0("  - ", warning_messages)
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

for (nm in names(output_paths)) {
  validate_output(
    output_paths[[nm]],
    nm
  )
}

msg("Ranked genes:", length(ranks))
msg("Hallmark pathways tested:", nrow(fgsea_results))
msg(
  "Significant Hallmarks:",
  nrow(fgsea_significant)
)
msg(
  "EMT NES:",
  signif(emt_result$NES, 8)
)
msg(
  "EMT FDR:",
  format_scientific(emt_result$padj, digits = 6)
)
msg("Metadata:", output_paths$metadata)
msg("QC:", output_paths$qc)
msg("Final status:", final_status)
