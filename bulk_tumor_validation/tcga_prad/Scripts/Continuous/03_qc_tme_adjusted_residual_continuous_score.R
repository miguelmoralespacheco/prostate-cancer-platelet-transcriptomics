#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD platelet-associated transcriptional score module
# TCGA-only sensitivity analysis:
# Continuous Score_z TME-adjusted and TME-residualized models
#
# Script:
# Scripts/Continuous/03_qc_tme_adjusted_residual_continuous_score.R
#
# Important:
# TME scores are transcriptomic proxy scores derived from bulk RNA-seq.
# They are not measured cellular fractions.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 1. Packages
################################################################################

required_packages <- c(
  "DESeq2",
  "apeglm",
  "fgsea",
  "GSVA",
  "ggplot2",
  "patchwork",
  "SummarizedExperiment"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them outside this script and rerun.",
    call. = FALSE
  )
}

################################################################################
# 2. Constants
################################################################################

EXPECTED_SCORE_SAMPLES <- 497L
EXPECTED_SIGNATURE_GENES <- 41L
EXPECTED_HALLMARK_PATHWAYS <- 50L

HALLMARK_EXPECTED_MD5 <- "56ef50e187bafc696128554ec7406702"

GSEA_SEED <- 123L
GSEA_NPROC <- 1L
GSEA_MIN_SIZE <- 15L
GSEA_MAX_SIZE <- 500L
GSEA_EPS <- 0

PADJ_THRESHOLD <- 0.05
ABS_SHRUNKEN_LFC_THRESHOLD <- 0.25

COL_REG_LINE <- "#2166AC"
COL_CONF_FILL <- "#92C5DE"

################################################################################
# 3. Helpers
################################################################################

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

warnings_log <- character()

add_warning <- function(...) {
  text <- paste0(...)
  warnings_log <<- c(warnings_log, text)
  warning(text, call. = FALSE)
}

read_csv_required <- function(path, label) {
  if (!file.exists(path)) fail("Missing ", label, ": ", path)

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

require_columns <- function(df, columns, label) {
  missing_columns <- setdiff(columns, colnames(df))

  if (length(missing_columns) > 0L) {
    fail(
      label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(x))
}

scale_numeric <- function(x, label) {
  x <- as_numeric_clean(x)

  if (anyNA(x) || any(!is.finite(x))) {
    fail(label, " contains missing or non-finite values.")
  }

  if (!is.finite(stats::sd(x)) || stats::sd(x) <= 0) {
    fail(label, " has zero or non-finite standard deviation.")
  }

  as.numeric(scale(x))
}

csv_safe <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  list_columns <- vapply(df, is.list, logical(1))

  if (any(list_columns)) {
    df[list_columns] <- lapply(
      df[list_columns],
      function(x) {
        vapply(
          x,
          function(xx) paste(as.character(xx), collapse = ";"),
          character(1)
        )
      }
    )
  }

  df
}

clean_pathway_name <- function(x) {
  x <- gsub("^HALLMARK_", "", as.character(x))
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

save_pdf <- function(plot, path, width_cm, height_cm) {
  grDevices::pdf(
    file = path,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    family = "Helvetica",
    useDingbats = FALSE
  )

  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

################################################################################
# 4. Paths
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

msg("Project directory:", project_dir)

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
  ),

  hallmark_snapshot = file.path(
    INPUT_DIR, "GeneSets",
    "MSigDB_Hallmark_Homo_sapiens_frozen.csv"
  ),

  unadjusted_gsea = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Tables",
    "GSEA",
    "GSEA_continuous_no_signature_Hallmark_full.csv"
  )
)

output_dirs <- list(
  tables = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Tables",
    "TME_Adjusted"
  ),

  figures = file.path(
    GENERATED_FIGURES_DIR, "Continuous",
    "TME_Adjusted"
  ),

  objects = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Objects",
    "TME_Adjusted"
  ),

  logs = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Logs"
  )
)

invisible(
  lapply(
    output_dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

output_paths <- list(
  tme_coverage = file.path(
    output_dirs$tables,
    "TME_geneSet_coverage.csv"
  ),

  tme_scores = file.path(
    output_dirs$tables,
    "TME_ssGSEA_scores_perSample.csv"
  ),

  tme_correlations = file.path(
    output_dirs$tables,
    "COR_Scorez_vs_TME_ssGSEA.csv"
  ),

  tme_covariate_correlations = file.path(
    output_dirs$tables,
    "TME_covariate_correlation.csv"
  ),

  deseq_tme_with_signature = file.path(
    output_dirs$tables,
    "DESeq2_TMEadjusted_WITH_SIGNATURE.csv"
  ),

  deseq_tme_no_signature = file.path(
    output_dirs$tables,
    "DESeq2_TMEadjusted_NO_SIGNATURE.csv"
  ),

  gsea_tme_full = file.path(
    output_dirs$tables,
    "GSEA_HALLMARK_TMEadjusted_full.csv"
  ),

  gsea_tme_significant = file.path(
    output_dirs$tables,
    "GSEA_HALLMARK_TMEadjusted_significant.csv"
  ),

  residual_summary = file.path(
    output_dirs$tables,
    "Score_residualized_TMEonly_model_summary.txt"
  ),

  residual_scores = file.path(
    output_dirs$tables,
    "Score_residualized_TMEonly_perSample.csv"
  ),

  deseq_residual_with_signature = file.path(
    output_dirs$tables,
    "DESeq2_RESIDUAL_TMEonly_WITH_SIGNATURE.csv"
  ),

  deseq_residual_no_signature = file.path(
    output_dirs$tables,
    "DESeq2_RESIDUAL_TMEonly_NO_SIGNATURE.csv"
  ),

  gsea_residual_full = file.path(
    output_dirs$tables,
    "GSEA_HALLMARK_RESIDUAL_TMEonly_full.csv"
  ),

  gsea_residual_significant = file.path(
    output_dirs$tables,
    "GSEA_HALLMARK_RESIDUAL_TMEonly_significant.csv"
  ),

  emt_summary = file.path(
    output_dirs$tables,
    "STEP3_HALLMARK_EMT_GSEA_platelet_score_models_summary.csv"
  ),

  dds_tme = file.path(
    output_dirs$objects,
    "DESeq2_dds_TMEadjusted.rds"
  ),

  dds_residual = file.path(
    output_dirs$objects,
    "DESeq2_dds_RESIDUAL_TMEonly.rds"
  ),

  figure_tme_stromal = file.path(
    output_dirs$figures,
    "STEP3_TME_StromalScore_correlation.pdf"
  ),

  figure_tme_immune = file.path(
    output_dirs$figures,
    "STEP3_TME_ImmuneScore_correlation.pdf"
  ),

  figure_tme_endothelial = file.path(
    output_dirs$figures,
    "STEP3_TME_EndothelialScore_correlation.pdf"
  ),

  figure_tme_epithelial = file.path(
    output_dirs$figures,
    "STEP3_TME_EpithelialScore_correlation.pdf"
  ),

  figure_emt_comparison = file.path(
    output_dirs$figures,
    "STEP3_HALLMARK_EMT_GSEA_NES_comparison.pdf"
  ),

  qc_log = file.path(
    output_dirs$logs,
    "03_qc_tme_adjusted_residual_continuous_score_QC.txt"
  )
)

for (nm in names(input_paths)) {
  if (!file.exists(input_paths[[nm]])) {
    fail("Missing required input '", nm, "': ", input_paths[[nm]])
  }
}

################################################################################
# 5. Input handling
################################################################################

align_exact_samples <- function(reference_ids, target_ids, label) {
  reference_ids <- trimws(as.character(reference_ids))
  target_ids <- trimws(as.character(target_ids))

  if (anyNA(reference_ids) || any(reference_ids == "")) {
    fail("Reference sample IDs contain missing or empty values.")
  }

  if (anyNA(target_ids) || any(target_ids == "")) {
    fail(label, " sample IDs contain missing or empty values.")
  }

  if (anyDuplicated(reference_ids)) {
    fail("Reference sample IDs are duplicated.")
  }

  if (anyDuplicated(target_ids)) {
    fail(label, " sample IDs are duplicated.")
  }

  index <- match(reference_ids, target_ids)

  if (anyNA(index)) {
    fail(
      "Exact sample alignment failed for ",
      label,
      ". Missing examples: ",
      paste(utils::head(reference_ids[is.na(index)], 8L), collapse = ", "),
      ". Barcode truncation is not allowed."
    )
  }

  index
}

load_counts_matrix <- function(path) {
  counts <- readRDS(path)

  if (inherits(counts, "SummarizedExperiment")) {
    assay_names <- SummarizedExperiment::assayNames(counts)
    preferred <- c("unstranded", "counts", "raw_counts")
    assay_name <- preferred[preferred %in% assay_names][1L]

    if (length(assay_name) == 0L || is.na(assay_name)) {
      assay_name <- assay_names[1L]
    }

    msg("Using count assay:", assay_name)

    counts <- SummarizedExperiment::assay(
      counts,
      assay_name
    )
  }

  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }

  if (is.null(dim(counts)) || length(dim(counts)) != 2L) {
    fail("Counts input is not a two-dimensional matrix.")
  }

  if (is.null(rownames(counts)) || anyNA(rownames(counts))) {
    fail("Counts matrix lacks valid gene row names.")
  }

  if (is.null(colnames(counts)) || anyNA(colnames(counts))) {
    fail("Counts matrix lacks valid sample column names.")
  }

  if (anyDuplicated(rownames(counts))) {
    fail("Counts matrix contains duplicated gene identifiers.")
  }

  if (anyDuplicated(colnames(counts))) {
    fail("Counts matrix contains duplicated sample identifiers.")
  }

  if (anyNA(counts) || any(counts < 0)) {
    fail("Counts matrix contains missing or negative values.")
  }

  if (!is.integer(counts)) {
    if (any(abs(counts - round(counts)) > 1e-6)) {
      fail("Counts matrix is not integer-like.")
    }

    storage.mode(counts) <- "integer"
  }

  counts
}

detect_signature_column <- function(df) {
  candidates <- c(
    "gene",
    "gene_name",
    "gene_symbol",
    "symbol",
    "hgnc_symbol",
    "hugo"
  )

  hit <- candidates[candidates %in% colnames(df)]

  if (length(hit) > 0L) {
    return(hit[1L])
  }

  normalized_names <- tolower(
    gsub("[^a-z0-9]+", "", colnames(df))
  )

  normalized_candidates <- tolower(
    gsub("[^a-z0-9]+", "", candidates)
  )

  hit_index <- which(
    normalized_names %in% normalized_candidates
  )

  if (length(hit_index) == 0L) {
    fail(
      "Could not detect signature gene column. Columns: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  colnames(df)[hit_index[1L]]
}

prepare_gene_map <- function(df) {
  require_columns(
    df,
    c("ensg_version", "gene_name"),
    "gene map"
  )

  df$ensg_version <- trimws(
    as.character(df$ensg_version)
  )

  df$gene_name <- toupper(
    trimws(as.character(df$gene_name))
  )

  df <- df[
    !is.na(df$ensg_version) &
      df$ensg_version != "",
    ,
    drop = FALSE
  ]

  duplicated_ids <- unique(
    df$ensg_version[duplicated(df$ensg_version)]
  )

  if (length(duplicated_ids) > 0L) {
    conflicting_ids <- duplicated_ids[
      vapply(
        duplicated_ids,
        function(id) {
          symbols <- unique(
            df$gene_name[
              df$ensg_version == id &
                !is.na(df$gene_name) &
                df$gene_name != ""
            ]
          )

          length(symbols) > 1L
        },
        logical(1)
      )
    ]

    if (length(conflicting_ids) > 0L) {
      fail(
        "Conflicting duplicate mappings for: ",
        paste(utils::head(conflicting_ids, 8L), collapse = ", ")
      )
    }

    df <- df[
      !duplicated(df$ensg_version),
      ,
      drop = FALSE
    ]
  }

  df
}

################################################################################
# 6. Hallmark snapshot
################################################################################

load_hallmark_snapshot <- function(path) {
  observed_md5 <- unname(
    tools::md5sum(path)
  )

  if (!identical(observed_md5, HALLMARK_EXPECTED_MD5)) {
    fail(
      "Hallmark snapshot MD5 mismatch. Expected ",
      HALLMARK_EXPECTED_MD5,
      "; observed ",
      observed_md5,
      "."
    )
  }

  hallmark <- read_csv_required(
    path,
    "Hallmark snapshot"
  )

  pathway_candidates <- c(
    "gs_name",
    "pathway",
    "gene_set",
    "geneset"
  )

  gene_candidates <- c(
    "gene_symbol",
    "gene_name",
    "symbol",
    "gene"
  )

  pathway_column <- pathway_candidates[
    pathway_candidates %in% colnames(hallmark)
  ][1L]

  gene_column <- gene_candidates[
    gene_candidates %in% colnames(hallmark)
  ][1L]

  if (
    length(pathway_column) == 0L ||
    length(gene_column) == 0L ||
    is.na(pathway_column) ||
    is.na(gene_column)
  ) {
    fail(
      "Could not identify Hallmark pathway/gene columns. Columns: ",
      paste(colnames(hallmark), collapse = ", ")
    )
  }

  hallmark$pathway_internal <- trimws(
    as.character(hallmark[[pathway_column]])
  )

  hallmark$gene_internal <- toupper(
    trimws(as.character(hallmark[[gene_column]]))
  )

  hallmark <- hallmark[
    !is.na(hallmark$pathway_internal) &
      hallmark$pathway_internal != "" &
      !is.na(hallmark$gene_internal) &
      hallmark$gene_internal != "",
    ,
    drop = FALSE
  ]

  pathways <- split(
    hallmark$gene_internal,
    hallmark$pathway_internal
  )

  pathways <- lapply(pathways, unique)
  pathways <- pathways[order(names(pathways))]

  if (length(pathways) != EXPECTED_HALLMARK_PATHWAYS) {
    fail(
      "Expected ",
      EXPECTED_HALLMARK_PATHWAYS,
      " Hallmark pathways; observed ",
      length(pathways),
      "."
    )
  }

  list(
    pathways = pathways,
    md5 = observed_md5,
    source_rows = nrow(hallmark)
  )
}

################################################################################
# 7. DESeq2 annotation
################################################################################

annotate_deseq_results <- function(
    mle_result,
    shrunken_result,
    gene_map,
    signature_symbols,
    signature_ensg
) {
  mle_df <- as.data.frame(mle_result)
  shrunk_df <- as.data.frame(shrunken_result)

  if (!identical(rownames(mle_df), rownames(shrunk_df))) {
    fail("MLE and shrunken DESeq2 results are not aligned.")
  }

  ensembl_id <- rownames(mle_df)

  gene_name <- gene_map$gene_name[
    match(ensembl_id, gene_map$ensg_version)
  ]

  fallback <- is.na(gene_name) | gene_name == ""
  gene_name[fallback] <- ensembl_id[fallback]

  result_df <- data.frame(
    ensembl_id = ensembl_id,
    gene_name = gene_name,
    gene_name_is_ensembl_fallback = fallback,
    baseMean = as_numeric_clean(mle_df$baseMean),
    log2FoldChange = as_numeric_clean(mle_df$log2FoldChange),
    lfcSE = as_numeric_clean(mle_df$lfcSE),
    stat = as_numeric_clean(mle_df$stat),
    pvalue = as_numeric_clean(mle_df$pvalue),
    padj = as_numeric_clean(mle_df$padj),
    log2FoldChange_shrunken = as_numeric_clean(
      shrunk_df$log2FoldChange
    ),
    shrinkage_method = "apeglm",
    stringsAsFactors = FALSE
  )

  result_df$is_signature_gene <-
    toupper(result_df$gene_name) %in% signature_symbols |
    result_df$ensembl_id %in% signature_ensg

  result_df$score_association <- ifelse(
    !is.na(result_df$padj) &
      result_df$padj < PADJ_THRESHOLD &
      result_df$log2FoldChange_shrunken >
      ABS_SHRUNKEN_LFC_THRESHOLD,
    "POSITIVE",
    ifelse(
      !is.na(result_df$padj) &
        result_df$padj < PADJ_THRESHOLD &
        result_df$log2FoldChange_shrunken <
        -ABS_SHRUNKEN_LFC_THRESHOLD,
      "NEGATIVE",
      "NS"
    )
  )

  result_df <- result_df[
    order(
      result_df$padj,
      -abs(result_df$stat),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(result_df) <- NULL
  result_df
}

validate_signature_representation <- function(
    result_df,
    signature_symbols,
    label
) {
  represented <- unique(
    toupper(
      result_df$gene_name[
        result_df$is_signature_gene &
          !result_df$gene_name_is_ensembl_fallback
      ]
    )
  )

  missing_symbols <- setdiff(
    signature_symbols,
    represented
  )

  if (length(missing_symbols) > 0L) {
    fail(
      label,
      " does not contain all signature genes. Missing: ",
      paste(missing_symbols, collapse = ", ")
    )
  }
}

################################################################################
# 8. GSEA
################################################################################

prepare_rank <- function(df) {
  require_columns(
    df,
    c(
      "ensembl_id",
      "gene_name",
      "gene_name_is_ensembl_fallback",
      "stat"
    ),
    "DESeq2 GSEA input"
  )

  rank_df <- data.frame(
    ensembl_id = trimws(
      as.character(df$ensembl_id)
    ),
    gene_name = toupper(
      trimws(as.character(df$gene_name))
    ),
    fallback = as.logical(
      df$gene_name_is_ensembl_fallback
    ),
    stat = as_numeric_clean(df$stat),
    stringsAsFactors = FALSE
  )

  input_rows <- nrow(rank_df)

  rank_df <- rank_df[
    !rank_df$fallback &
      !is.na(rank_df$gene_name) &
      rank_df$gene_name != "" &
      !is.na(rank_df$ensembl_id) &
      rank_df$ensembl_id != "" &
      is.finite(rank_df$stat),
    ,
    drop = FALSE
  ]

  usable_rows <- nrow(rank_df)

  rank_df <- rank_df[
    order(
      rank_df$gene_name,
      -abs(rank_df$stat),
      rank_df$ensembl_id
    ),
    ,
    drop = FALSE
  ]

  duplicates_removed <- sum(
    duplicated(rank_df$gene_name)
  )

  rank_df <- rank_df[
    !duplicated(rank_df$gene_name),
    ,
    drop = FALSE
  ]

  rank_df <- rank_df[
    order(
      -rank_df$stat,
      rank_df$gene_name,
      rank_df$ensembl_id
    ),
    ,
    drop = FALSE
  ]

  ranks <- stats::setNames(
    rank_df$stat,
    rank_df$gene_name
  )

  if (
    length(ranks) == 0L ||
    anyDuplicated(names(ranks)) ||
    any(!is.finite(ranks))
  ) {
    fail("Invalid final GSEA ranking.")
  }

  list(
    ranks = ranks,
    input_rows = input_rows,
    usable_rows = usable_rows,
    ranked_symbols = length(ranks),
    duplicate_rows_removed = duplicates_removed
  )
}

run_fgsea <- function(pathways, ranks, model_label, positive_label) {
  fgsea_warnings <- character()

  set.seed(GSEA_SEED)

  result <- withCallingHandlers(
    fgsea::fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = GSEA_MIN_SIZE,
      maxSize = GSEA_MAX_SIZE,
      eps = GSEA_EPS,
      nproc = GSEA_NPROC
    ),
    warning = function(w) {
      fgsea_warnings <<- c(
        fgsea_warnings,
        conditionMessage(w)
      )

      invokeRestart("muffleWarning")
    }
  )

  if (length(fgsea_warnings) > 0L) {
    for (warning_text in unique(fgsea_warnings)) {
      add_warning(
        model_label,
        " fgsea warning: ",
        warning_text
      )
    }
  }

  result <- as.data.frame(
    result,
    stringsAsFactors = FALSE
  )

  require_columns(
    result,
    c(
      "pathway",
      "pval",
      "padj",
      "ES",
      "NES",
      "size",
      "leadingEdge"
    ),
    paste0(model_label, " fgsea result")
  )

  if (nrow(result) != EXPECTED_HALLMARK_PATHWAYS) {
    fail(
      model_label,
      " returned ",
      nrow(result),
      " pathways; expected ",
      EXPECTED_HALLMARK_PATHWAYS,
      "."
    )
  }

  result <- result[
    order(
      result$padj,
      -abs(result$NES),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  result$direction <- ifelse(
    result$NES > 0,
    positive_label,
    paste0("Opposite to ", positive_label)
  )

  result$pathway_clean <- clean_pathway_name(
    result$pathway
  )

  result$model <- model_label

  csv_safe(result)
}

write_gsea_results <- function(
    result,
    full_path,
    significant_path
) {
  significant <- result[
    !is.na(result$padj) &
      result$padj < PADJ_THRESHOLD,
    ,
    drop = FALSE
  ]

  emt <- result[
    result$pathway ==
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ,
    drop = FALSE
  ]

  if (nrow(emt) != 1L) {
    fail("Could not identify a unique Hallmark EMT result.")
  }

  utils::write.csv(
    csv_safe(result),
    full_path,
    row.names = FALSE
  )

  utils::write.csv(
    csv_safe(significant),
    significant_path,
    row.names = FALSE
  )

  list(
    n_tested = nrow(result),
    n_significant = nrow(significant),
    emt = emt
  )
}

################################################################################
# 9. TME proxy scoring
################################################################################

collapse_expression_by_symbol <- function(expression_matrix, gene_map) {
  symbols <- gene_map$gene_name[
    match(
      rownames(expression_matrix),
      gene_map$ensg_version
    )
  ]

  keep <- !is.na(symbols) & symbols != ""

  expression_matrix <- expression_matrix[
    keep,
    ,
    drop = FALSE
  ]

  symbols <- symbols[keep]

  row_variance <- apply(
    expression_matrix,
    1L,
    stats::var
  )

  split_indices <- split(
    seq_along(symbols),
    symbols
  )

  selected <- vapply(
    split_indices,
    function(indices) {
      indices[which.max(row_variance[indices])]
    },
    integer(1)
  )

  collapsed <- expression_matrix[
    selected,
    ,
    drop = FALSE
  ]

  rownames(collapsed) <- names(selected)
  collapsed
}

calculate_gene_set_coverage <- function(gene_sets, available_symbols) {
  do.call(
    rbind,
    lapply(
      names(gene_sets),
      function(set_name) {
        defined <- unique(
          toupper(gene_sets[[set_name]])
        )

        present <- intersect(
          defined,
          available_symbols
        )

        missing <- setdiff(
          defined,
          available_symbols
        )

        data.frame(
          gene_set = set_name,
          n_genes_defined = length(defined),
          n_genes_present = length(present),
          n_genes_missing = length(missing),
          genes_present = paste(present, collapse = ";"),
          genes_missing = paste(missing, collapse = ";"),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

run_ssgsea <- function(expression_matrix, gene_sets) {
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    parameter <- tryCatch(
      GSVA::ssgseaParam(
        exprData = expression_matrix,
        geneSets = gene_sets,
        normalize = TRUE
      ),
      error = function(e) {
        GSVA::ssgseaParam(
          expr = expression_matrix,
          geneSets = gene_sets
        )
      }
    )

    return(
      tryCatch(
        GSVA::gsva(parameter, verbose = FALSE),
        error = function(e) GSVA::gsva(parameter)
      )
    )
  }

  tryCatch(
    GSVA::gsva(
      expression_matrix,
      gene_sets,
      method = "ssgsea",
      ssgsea.norm = TRUE,
      verbose = FALSE
    ),
    error = function(e) {
      GSVA::gsva(
        expression_matrix,
        gene_sets,
        method = "ssgsea",
        ssgsea.norm = TRUE
      )
    }
  )
}

calculate_correlations <- function(
    data,
    score_column,
    variable_columns
) {
  results <- list()

  for (method in c("pearson", "spearman")) {
    for (variable in variable_columns) {
      complete <- is.finite(data[[score_column]]) &
        is.finite(data[[variable]])

      test <- suppressWarnings(
        stats::cor.test(
          data[[score_column]][complete],
          data[[variable]][complete],
          method = method,
          exact = FALSE
        )
      )

      results[[length(results) + 1L]] <- data.frame(
        variable = variable,
        method = method,
        n_complete = sum(complete),
        correlation = unname(test$estimate),
        pvalue = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}

################################################################################
# 10. EMT summary
################################################################################

extract_emt <- function(
    result,
    model,
    statistical_model,
    direction,
    interpretation
) {
  require_columns(
    result,
    c("pathway", "NES", "padj"),
    paste0(model, " GSEA result")
  )

  emt <- result[
    result$pathway ==
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ,
    drop = FALSE
  ]

  if (nrow(emt) != 1L) {
    fail(
      "Could not identify a unique Hallmark EMT result for ",
      model,
      "."
    )
  }

  data.frame(
    model = model,
    statistical_model = statistical_model,
    pathway = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    NES = as_numeric_clean(emt$NES),
    ES = if ("ES" %in% colnames(emt)) {
      as_numeric_clean(emt$ES)
    } else {
      NA_real_
    },
    pval = if ("pval" %in% colnames(emt)) {
      as_numeric_clean(emt$pval)
    } else {
      NA_real_
    },
    padj = as_numeric_clean(emt$padj),
    log2err = if ("log2err" %in% colnames(emt)) {
      as_numeric_clean(emt$log2err)
    } else {
      NA_real_
    },
    direction = direction,
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

################################################################################
# 11. Figures
################################################################################

plot_tme_correlations <- function(
    tme_scores,
    correlations,
    output_paths
) {
  get_spearman_stats <- function(variable) {
    hit <- correlations[
      correlations$variable == variable &
        correlations$method == "spearman",
      ,
      drop = FALSE
    ]

    if (nrow(hit) != 1L) {
      fail(
        "Could not identify a unique Spearman result for ",
        variable,
        "."
      )
    }

    list(
      rho = hit$correlation[1L],
      padj = hit$padj_BH[1L]
    )
  }

  theme_small <- ggplot2::theme_classic(
    base_family = "Helvetica",
    base_size = 5
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = 4.3,
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        size = 4.6,
        color = "black"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.20,
        color = "black"
      ),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.20
      ),
      panel.grid = ggplot2::element_blank(),
      plot.margin = grid::unit(
        c(1.2, 1.2, 1.2, 1.2),
        "mm"
      )
    )

  make_panel <- function(variable, x_label) {
    panel_data <- tme_scores[
      is.finite(tme_scores[[variable]]) &
        is.finite(tme_scores$Score_z),
      ,
      drop = FALSE
    ]

    if (nrow(panel_data) == 0L) {
      fail(
        "No complete observations available for ",
        variable,
        "."
      )
    }

    x_range <- range(
      panel_data[[variable]],
      na.rm = TRUE
    )

    y_range <- range(
      panel_data$Score_z,
      na.rm = TRUE
    )

    if (
      !all(is.finite(x_range)) ||
      diff(x_range) <= 0 ||
      !all(is.finite(y_range)) ||
      diff(y_range) <= 0
    ) {
      fail(
        "Invalid plotting range for ",
        variable,
        "."
      )
    }

    spearman_stats <- get_spearman_stats(
      variable
    )

    x_position <- x_range[1L] +
      0.04 * diff(x_range)

    y_position_rho <- y_range[1L] +
      0.105 * diff(y_range)

    y_position_p <- y_range[1L] +
      0.025 * diff(y_range)

    ggplot2::ggplot(
      panel_data,
      ggplot2::aes(
        x = .data[[variable]],
        y = Score_z
      )
    ) +
      ggplot2::geom_point(
        size = 0.55,
        shape = 16,
        color = "black",
        alpha = 0.55,
        stroke = 0
      ) +
      ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.25,
        color = COL_REG_LINE,
        fill = COL_CONF_FILL
      ) +
      ggplot2::labs(
        x = x_label,
        y = "Platelet score (z)"
      ) +
      ggplot2::annotate(
        "text",
        x = x_position,
        y = y_position_rho,
        label = paste0(
          "Spearman rho = ",
          sprintf("%.2f", spearman_stats$rho)
        ),
        hjust = 0,
        vjust = 0,
        size = 4.7 / ggplot2::.pt,
        color = "black",
        family = "Helvetica"
      ) +
      ggplot2::annotate(
        "text",
        x = x_position,
        y = y_position_p,
        label = paste0(
          "'BH-adjusted'~italic(P)~'='~'",
          formatC(
            spearman_stats$padj,
            format = "e",
            digits = 2
          ),
          "'"
        ),
        parse = TRUE,
        hjust = 0,
        vjust = 0,
        size = 4.7 / ggplot2::.pt,
        color = "black",
        family = "Helvetica"
      ) +
      theme_small
  }

  stromal_plot <- make_panel(
    "StromalScore_ESTlike",
    "Stromal ssGSEA score"
  )

  immune_plot <- make_panel(
    "ImmuneScore_ESTlike",
    "Immune ssGSEA score"
  )

  endothelial_plot <- make_panel(
    "EndothelialScore",
    "Endothelial ssGSEA score"
  )

  epithelial_plot <- make_panel(
    "EpithelialScore",
    "Epithelial ssGSEA score"
  )

  save_pdf(
    stromal_plot,
    output_paths$stromal,
    width_cm = 4.6,
    height_cm = 4.2
  )

  save_pdf(
    immune_plot,
    output_paths$immune,
    width_cm = 4.6,
    height_cm = 4.2
  )

  save_pdf(
    endothelial_plot,
    output_paths$endothelial,
    width_cm = 4.6,
    height_cm = 4.2
  )

  save_pdf(
    epithelial_plot,
    output_paths$epithelial,
    width_cm = 4.6,
    height_cm = 4.2
  )
}

plot_emt_comparison <- function(summary_df, output_path) {
  summary_df$model_label <- summary_df$model

  summary_df$model_label[
    summary_df$model == "Continuous unadjusted"
  ] <- "Platelet score"

  summary_df$model_label[
    summary_df$model == "Continuous TME-adjusted"
  ] <- "Platelet score + TME\nadjustment"

  summary_df$model_label[
    summary_df$model == "Continuous TME-residualized"
  ] <- "TME-residualized\nplatelet score"

  summary_df$model_label <- factor(
    summary_df$model_label,
    levels = c(
      "TME-residualized\nplatelet score",
      "Platelet score + TME\nadjustment",
      "Platelet score"
    )
  )

  plot <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = NES,
      y = model_label
    )
  ) +
    ggplot2::geom_col(
      width = 0.68,
      fill = "#A9C4F0",
      color = "grey55",
      linewidth = 0.22
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-0.14, 4.1),
      breaks = 0:4,
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    ggplot2::labs(
      x = "Hallmark EMT normalized\nenrichment score",
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_family = "Helvetica",
      base_size = 5
    ) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(
        linewidth = 0.20,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.20,
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 4.5,
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 5,
        color = "black"
      ),
      axis.title.x = ggplot2::element_text(
        size = 5,
        color = "black"
      ),
      panel.grid = ggplot2::element_blank(),
      plot.margin = grid::unit(
        c(2, 2, 2, 2),
        "mm"
      )
    )

  save_pdf(
    plot,
    output_path,
    width_cm = 4.6,
    height_cm = 2.6
  )
}

################################################################################
# 12. Load canonical score cohort
################################################################################

msg("Loading canonical score cohort")

score_df <- read_csv_required(
  input_paths$score,
  "canonical score table"
)

require_columns(
  score_df,
  c(
    "sample_id",
    "patient_id",
    "Score_raw",
    "Score_z"
  ),
  "canonical score table"
)

score_df$sample_id <- trimws(
  as.character(score_df$sample_id)
)

score_df$patient_id <- trimws(
  as.character(score_df$patient_id)
)

score_df$Score_raw <- as_numeric_clean(
  score_df$Score_raw
)

score_df$Score_z <- as_numeric_clean(
  score_df$Score_z
)

if (nrow(score_df) != EXPECTED_SCORE_SAMPLES) {
  fail(
    "Expected ",
    EXPECTED_SCORE_SAMPLES,
    " score rows; observed ",
    nrow(score_df),
    "."
  )
}

if (
  anyNA(score_df$sample_id) ||
  any(score_df$sample_id == "") ||
  anyDuplicated(score_df$sample_id)
) {
  fail("Canonical score table must contain unique sample_id values.")
}

if (
  anyNA(score_df$patient_id) ||
  any(score_df$patient_id == "") ||
  anyDuplicated(score_df$patient_id)
) {
  fail(
    "Canonical score table must contain ",
    EXPECTED_SCORE_SAMPLES,
    " unique patient_id values."
  )
}

if (
  anyNA(score_df$Score_raw) ||
  anyNA(score_df$Score_z) ||
  any(!is.finite(score_df$Score_raw)) ||
  any(!is.finite(score_df$Score_z))
) {
  fail("Score_raw and Score_z must be finite.")
}

if (abs(mean(score_df$Score_z)) > 1e-8) {
  fail("Score_z mean differs from zero beyond tolerance.")
}

if (abs(stats::sd(score_df$Score_z) - 1) > 1e-8) {
  fail("Score_z standard deviation differs from one beyond tolerance.")
}

################################################################################
# 13. Counts, map, signature and Hallmark
################################################################################

msg("Loading and aligning count matrix")

counts <- load_counts_matrix(
  input_paths$counts
)

counts_dimensions_before_alignment <- dim(counts)

count_index <- align_exact_samples(
  score_df$sample_id,
  colnames(counts),
  "count matrix"
)

counts <- counts[
  ,
  count_index,
  drop = FALSE
]

if (!identical(colnames(counts), score_df$sample_id)) {
  fail("Count matrix alignment failed.")
}

raw_gene_count <- nrow(counts)

msg("Loading gene map and signature")

gene_map <- prepare_gene_map(
  read_csv_required(
    input_paths$gene_map,
    "gene map"
  )
)

signature_df <- read_canonical_platelet_signature()

signature_column <- detect_signature_column(
  signature_df
)

signature_symbols <- unique(
  toupper(
    trimws(
      as.character(signature_df[[signature_column]])
    )
  )
)

signature_symbols <- signature_symbols[
  !is.na(signature_symbols) &
    signature_symbols != ""
]

if (length(signature_symbols) != EXPECTED_SIGNATURE_GENES) {
  fail(
    "Expected ",
    EXPECTED_SIGNATURE_GENES,
    " unique signature genes; observed ",
    length(signature_symbols),
    "."
  )
}

mapped_signature_symbols <- intersect(
  signature_symbols,
  unique(gene_map$gene_name)
)

if (length(mapped_signature_symbols) != EXPECTED_SIGNATURE_GENES) {
  fail(
    "Signature genes missing from gene map: ",
    paste(
      setdiff(
        signature_symbols,
        mapped_signature_symbols
      ),
      collapse = ", "
    )
  )
}

signature_ensg <- unique(
  gene_map$ensg_version[
    gene_map$gene_name %in% signature_symbols
  ]
)

msg("Loading frozen Hallmark snapshot")

hallmark_snapshot <- load_hallmark_snapshot(
  input_paths$hallmark_snapshot
)

hallmark_pathways <- hallmark_snapshot$pathways

msg("Loading unadjusted continuous GSEA")

unadjusted_gsea <- read_csv_required(
  input_paths$unadjusted_gsea,
  "unadjusted continuous GSEA"
)

require_columns(
  unadjusted_gsea,
  c(
    "pathway",
    "NES",
    "padj"
  ),
  "unadjusted continuous GSEA"
)

if (nrow(unadjusted_gsea) != EXPECTED_HALLMARK_PATHWAYS) {
  fail(
    "Unadjusted GSEA table must contain ",
    EXPECTED_HALLMARK_PATHWAYS,
    " pathways."
  )
}

################################################################################
# 14. Prefilter and VST
################################################################################

msg("Applying prefilter: rowSums(counts >= 10) >= 10")

genes_before_prefilter <- nrow(counts)

keep <- rowSums(counts >= 10) >= 10

counts_prefiltered <- counts[
  keep,
  ,
  drop = FALSE
]

genes_after_prefilter <- nrow(
  counts_prefiltered
)

rm(counts)
invisible(gc(verbose = FALSE))

coldata_base <- data.frame(
  Score_raw = score_df$Score_raw,
  Score_z = score_df$Score_z,
  row.names = score_df$sample_id,
  stringsAsFactors = FALSE
)

dds_vst <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_prefiltered,
  colData = coldata_base,
  design = ~ 1
)

dds_vst <- DESeq2::estimateSizeFactors(
  dds_vst
)

vst_object <- DESeq2::vst(
  dds_vst,
  blind = TRUE
)

vst_matrix <- SummarizedExperiment::assay(
  vst_object
)

expression_by_symbol <- collapse_expression_by_symbol(
  vst_matrix,
  gene_map
)

################################################################################
# 15. TME transcriptomic proxies
################################################################################

msg("Calculating TME transcriptomic proxy scores")

tme_gene_sets <- list(
  StromalScore_ESTlike = c(
    "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1",
    "FAP", "PDGFRB", "TAGLN", "SPARC", "VCAN"
  ),

  ImmuneScore_ESTlike = c(
    "PTPRC", "LST1", "TYROBP", "LYZ", "HLA-DRA",
    "HLA-DRB1", "CD74", "CSF1R", "FCER1G", "CTSS"
  ),

  EndothelialScore = c(
    "PECAM1", "VWF", "KDR", "EMCN", "ESAM",
    "RAMP2", "PLVAP", "ENG", "CD34", "KLF2"
  ),

  EpithelialScore = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "MUC1",
    "TACSTD2", "EHF", "CLDN3", "CLDN4", "KRT7"
  )
)

tme_coverage <- calculate_gene_set_coverage(
  tme_gene_sets,
  rownames(expression_by_symbol)
)

utils::write.csv(
  tme_coverage,
  output_paths$tme_coverage,
  row.names = FALSE
)

incomplete_sets <- tme_coverage$gene_set[
  tme_coverage$n_genes_missing > 0
]

if (length(incomplete_sets) > 0L) {
  add_warning(
    "Incomplete TME proxy coverage: ",
    paste(incomplete_sets, collapse = ", ")
  )
}

tme_matrix <- run_ssgsea(
  expression_by_symbol,
  tme_gene_sets
)

tme_raw <- as.data.frame(
  t(tme_matrix),
  stringsAsFactors = FALSE
)

tme_raw$sample_id <- rownames(tme_raw)

tme_index <- align_exact_samples(
  score_df$sample_id,
  tme_raw$sample_id,
  "TME proxy matrix"
)

tme_raw <- tme_raw[
  tme_index,
  ,
  drop = FALSE
]

tme_columns <- c(
  "StromalScore_ESTlike",
  "ImmuneScore_ESTlike",
  "EndothelialScore",
  "EpithelialScore"
)

require_columns(
  tme_raw,
  tme_columns,
  "TME proxy matrix"
)

tme_scores <- data.frame(
  sample_id = score_df$sample_id,
  patient_id = score_df$patient_id,
  Score_raw = score_df$Score_raw,
  Score_z = score_df$Score_z,
  stringsAsFactors = FALSE
)

for (column in tme_columns) {
  tme_scores[[column]] <- as_numeric_clean(
    tme_raw[[column]]
  )

  if (
    anyNA(tme_scores[[column]]) ||
    any(!is.finite(tme_scores[[column]]))
  ) {
    fail(column, " contains non-finite values.")
  }
}

utils::write.csv(
  tme_scores,
  output_paths$tme_scores,
  row.names = FALSE
)

tme_correlations <- calculate_correlations(
  tme_scores,
  "Score_z",
  tme_columns
)
tme_correlations$padj_BH <- NA_real_

for (method_name in unique(tme_correlations$method)) {
  idx <- tme_correlations$method == method_name

  tme_correlations$padj_BH[idx] <- stats::p.adjust(
    tme_correlations$pvalue[idx],
    method = "BH"
  )
}

utils::write.csv(
  tme_correlations,
  output_paths$tme_correlations,
  row.names = FALSE
)

################################################################################
# 16. TME-adjusted DESeq2
################################################################################

msg("Running TME-adjusted continuous DESeq2 model")

tme_coldata <- data.frame(
  StromalScore_ESTlike = scale_numeric(
    tme_scores$StromalScore_ESTlike,
    "StromalScore_ESTlike"
  ),
  ImmuneScore_ESTlike = scale_numeric(
    tme_scores$ImmuneScore_ESTlike,
    "ImmuneScore_ESTlike"
  ),
  EndothelialScore = scale_numeric(
    tme_scores$EndothelialScore,
    "EndothelialScore"
  ),
  Score_z = tme_scores$Score_z,
  row.names = tme_scores$sample_id,
  stringsAsFactors = FALSE
)

tme_coldata <- tme_coldata[
  colnames(counts_prefiltered),
  ,
  drop = FALSE
]

tme_design <- ~
  StromalScore_ESTlike +
  ImmuneScore_ESTlike +
  EndothelialScore +
  Score_z

tme_model_matrix <- stats::model.matrix(
  tme_design,
  data = tme_coldata
)

if (qr(tme_model_matrix)$rank < ncol(tme_model_matrix)) {
  fail("TME-adjusted model matrix is rank deficient.")
}

tme_covariate_correlations <- stats::cor(
  tme_coldata,
  use = "pairwise.complete.obs"
)

utils::write.csv(
  data.frame(
    variable = rownames(tme_covariate_correlations),
    tme_covariate_correlations,
    row.names = NULL,
    check.names = FALSE
  ),
  output_paths$tme_covariate_correlations,
  row.names = FALSE
)

dds_tme <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_prefiltered,
  colData = tme_coldata,
  design = tme_design
)

dds_tme <- DESeq2::DESeq(
  dds_tme,
  quiet = FALSE
)

saveRDS(
  dds_tme,
  output_paths$dds_tme
)

tme_coefficient <- DESeq2::resultsNames(
  dds_tme
)

tme_coefficient <- tme_coefficient[
  grepl("Score_z", tme_coefficient, fixed = TRUE)
]

if (length(tme_coefficient) != 1L) {
  fail(
    "Could not identify unique Score_z coefficient. resultsNames: ",
    paste(DESeq2::resultsNames(dds_tme), collapse = ", ")
  )
}

tme_mle <- DESeq2::results(
  dds_tme,
  name = tme_coefficient,
  alpha = PADJ_THRESHOLD
)

tme_shrunken <- DESeq2::lfcShrink(
  dds_tme,
  coef = tme_coefficient,
  type = "apeglm"
)

tme_with_signature <- annotate_deseq_results(
  tme_mle,
  tme_shrunken,
  gene_map,
  signature_symbols,
  signature_ensg
)

validate_signature_representation(
  tme_with_signature,
  signature_symbols,
  "TME-adjusted model"
)

tme_no_signature <- tme_with_signature[
  !tme_with_signature$is_signature_gene,
  ,
  drop = FALSE
]

utils::write.csv(
  tme_with_signature,
  output_paths$deseq_tme_with_signature,
  row.names = FALSE
)

utils::write.csv(
  tme_no_signature,
  output_paths$deseq_tme_no_signature,
  row.names = FALSE
)

################################################################################
# 17. TME-adjusted GSEA
################################################################################

msg("Running TME-adjusted Hallmark GSEA")

tme_rank <- prepare_rank(
  tme_no_signature
)

gsea_tme <- run_fgsea(
  hallmark_pathways,
  tme_rank$ranks,
  "Continuous TME-adjusted",
  "Higher Score_z in TME-adjusted model"
)

tme_gsea_info <- write_gsea_results(
  gsea_tme,
  output_paths$gsea_tme_full,
  output_paths$gsea_tme_significant
)

################################################################################
# 18. Residualized score
################################################################################

msg("Residualizing Score_z against TME proxies")

residual_data <- data.frame(
  sample_id = tme_scores$sample_id,
  Score_z = tme_scores$Score_z,
  StromalScore_ESTlike = scale_numeric(
    tme_scores$StromalScore_ESTlike,
    "StromalScore_ESTlike"
  ),
  ImmuneScore_ESTlike = scale_numeric(
    tme_scores$ImmuneScore_ESTlike,
    "ImmuneScore_ESTlike"
  ),
  EndothelialScore = scale_numeric(
    tme_scores$EndothelialScore,
    "EndothelialScore"
  ),
  stringsAsFactors = FALSE
)

residual_fit <- stats::lm(
  Score_z ~
    StromalScore_ESTlike +
    ImmuneScore_ESTlike +
    EndothelialScore,
  data = residual_data
)

residual_data$Score_resid <- as.numeric(
  stats::residuals(residual_fit)
)

residual_data$Score_resid_z <- scale_numeric(
  residual_data$Score_resid,
  "Score_resid"
)

residual_r2 <- summary(residual_fit)$r.squared
residual_adjusted_r2 <- summary(residual_fit)$adj.r.squared

writeLines(
  c(
    "TCGA-only platelet-score TME residualization",
    "============================================",
    "",
    paste("Date/time:", as.character(Sys.time())),
    "",
    paste0(
      "Model: Score_z ~ StromalScore_ESTlike + ",
      "ImmuneScore_ESTlike + EndothelialScore"
    ),
    "",
    "TME variables are transcriptomic proxies, not measured cell fractions.",
    "",
    capture.output(summary(residual_fit))
  ),
  con = output_paths$residual_summary,
  useBytes = TRUE
)

utils::write.csv(
  residual_data,
  output_paths$residual_scores,
  row.names = FALSE
)

################################################################################
# 19. Residualized DESeq2
################################################################################

msg("Running TME-residualized continuous DESeq2 model")

residual_coldata <- data.frame(
  Score_resid_z = residual_data$Score_resid_z,
  row.names = residual_data$sample_id,
  stringsAsFactors = FALSE
)

residual_coldata <- residual_coldata[
  colnames(counts_prefiltered),
  ,
  drop = FALSE
]

dds_residual <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_prefiltered,
  colData = residual_coldata,
  design = ~ Score_resid_z
)

dds_residual <- DESeq2::DESeq(
  dds_residual,
  quiet = FALSE
)

saveRDS(
  dds_residual,
  output_paths$dds_residual
)

residual_coefficient <- DESeq2::resultsNames(
  dds_residual
)

residual_coefficient <- residual_coefficient[
  grepl(
    "Score_resid_z",
    residual_coefficient,
    fixed = TRUE
  )
]

if (length(residual_coefficient) != 1L) {
  fail(
    "Could not identify unique Score_resid_z coefficient. resultsNames: ",
    paste(DESeq2::resultsNames(dds_residual), collapse = ", ")
  )
}

residual_mle <- DESeq2::results(
  dds_residual,
  name = residual_coefficient,
  alpha = PADJ_THRESHOLD
)

residual_shrunken <- DESeq2::lfcShrink(
  dds_residual,
  coef = residual_coefficient,
  type = "apeglm"
)

residual_with_signature <- annotate_deseq_results(
  residual_mle,
  residual_shrunken,
  gene_map,
  signature_symbols,
  signature_ensg
)

validate_signature_representation(
  residual_with_signature,
  signature_symbols,
  "TME-residualized model"
)

residual_no_signature <- residual_with_signature[
  !residual_with_signature$is_signature_gene,
  ,
  drop = FALSE
]

utils::write.csv(
  residual_with_signature,
  output_paths$deseq_residual_with_signature,
  row.names = FALSE
)

utils::write.csv(
  residual_no_signature,
  output_paths$deseq_residual_no_signature,
  row.names = FALSE
)

################################################################################
# 20. Residualized GSEA
################################################################################

msg("Running TME-residualized Hallmark GSEA")

residual_rank <- prepare_rank(
  residual_no_signature
)

gsea_residual <- run_fgsea(
  hallmark_pathways,
  residual_rank$ranks,
  "Continuous TME-residualized",
  "Higher TME-residualized Score_z"
)

residual_gsea_info <- write_gsea_results(
  gsea_residual,
  output_paths$gsea_residual_full,
  output_paths$gsea_residual_significant
)

################################################################################
# 21. EMT comparison
################################################################################

msg("Building three-model Hallmark EMT comparison")

emt_summary <- rbind(
  extract_emt(
    unadjusted_gsea,
    model = "Continuous unadjusted",
    statistical_model = "expression ~ Score_z",
    direction = "Higher Score_z",
    interpretation = "Unadjusted continuous platelet-score association."
  ),

  extract_emt(
    gsea_tme,
    model = "Continuous TME-adjusted",
    statistical_model = paste0(
      "expression ~ StromalScore_ESTlike + ImmuneScore_ESTlike + ",
      "EndothelialScore + Score_z"
    ),
    direction = "Higher Score_z after TME proxy adjustment",
    interpretation = paste0(
      "Association after adjustment for stromal, immune and ",
      "endothelial transcriptomic proxy scores."
    )
  ),

  extract_emt(
    gsea_residual,
    model = "Continuous TME-residualized",
    statistical_model = "expression ~ Score_resid_z",
    direction = "Higher TME-residualized Score_z",
    interpretation = paste0(
      "Association with the score component residualized against ",
      "stromal, immune and endothelial transcriptomic proxies."
    )
  )
)

utils::write.csv(
  emt_summary,
  output_paths$emt_summary,
  row.names = FALSE
)

################################################################################
# 22. Figures
################################################################################

msg("Generating TME correlation figure")

plot_tme_correlations(
  tme_scores,
  tme_correlations,
  list(
    stromal = output_paths$figure_tme_stromal,
    immune = output_paths$figure_tme_immune,
    endothelial = output_paths$figure_tme_endothelial,
    epithelial = output_paths$figure_tme_epithelial
  )
)

msg("Generating EMT comparison figure")

plot_emt_comparison(
  emt_summary,
  output_paths$figure_emt_comparison
)

################################################################################
# 23. Final QC
################################################################################

count_fdr <- function(df) {
  sum(!is.na(df$padj) & df$padj < PADJ_THRESHOLD)
}

count_associated <- function(df, direction = NULL) {
  if (is.null(direction)) {
    return(
      sum(
        df$score_association %in% c(
          "POSITIVE",
          "NEGATIVE"
        )
      )
    )
  }

  sum(df$score_association == direction)
}

emt_text <- function(df) {
  paste0(
    "ES=", signif(df$ES[1L], 6),
    "; NES=", signif(df$NES[1L], 6),
    "; padj=", signif(df$padj[1L], 6)
  )
}

required_outputs <- unname(
  unlist(output_paths[
    names(output_paths) != "qc_log"
  ])
)

missing_outputs <- required_outputs[
  !file.exists(required_outputs)
]

if (length(missing_outputs) > 0L) {
  fail(
    "Missing generated outputs: ",
    paste(missing_outputs, collapse = " | ")
  )
}

final_status <- if (length(warnings_log) == 0L) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

qc_lines <- c(
  "TCGA-PRAD continuous platelet-associated transcriptional score",
  "TME-adjusted and TME-residualized sensitivity analysis",
  "============================================================",
  "",
  paste("Date/time:", as.character(Sys.time())),
  paste("Project directory:", project_dir),
  paste("Final status:", final_status),
  "",
  "ANALYTICAL SCOPE",
  "----------------",
  "TCGA-only sensitivity analysis.",
  "TME variables are transcriptomic proxies, not measured cellular fractions.",
  "",
  "COHORT",
  "------",
  paste("Canonical samples:", nrow(score_df)),
  paste("Unique patients:", length(unique(score_df$patient_id))),
  paste(
    "Count dimensions before alignment:",
    paste(counts_dimensions_before_alignment, collapse = " x ")
  ),
  paste("Aligned count samples:", ncol(counts_prefiltered)),
  "Alignment method: exact full sample_id",
  "Barcode truncation: not used",
  "",
  "GENES",
  "-----",
  paste("Raw genes:", raw_gene_count),
  paste("Genes before prefilter:", genes_before_prefilter),
  paste("Genes after prefilter:", genes_after_prefilter),
  "Prefilter: rowSums(counts >= 10) >= 10",
  paste("Signature genes:", length(signature_symbols)),
  "",
  "HALLMARK",
  "--------",
  paste("Snapshot:", input_paths$hallmark_snapshot),
  paste("MD5:", hallmark_snapshot$md5),
  paste("Pathways:", length(hallmark_pathways)),
  paste("Seed:", GSEA_SEED),
  paste("nproc:", GSEA_NPROC),
  paste("minSize:", GSEA_MIN_SIZE),
  paste("maxSize:", GSEA_MAX_SIZE),
  paste("eps:", GSEA_EPS),
  "Dynamic msigdbr query: not used",
  "Ranking cap: not used",
  "Ranking jitter: not used",
  "Ranking transformation: not used",
  "",
  "TME PROXY COVERAGE",
  "------------------",
  paste(
    capture.output(
      print(
        tme_coverage[
          ,
          c(
            "gene_set",
            "n_genes_defined",
            "n_genes_present",
            "n_genes_missing"
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "SCORE_Z VERSUS TME PROXIES",
  "--------------------------",
  paste(
    capture.output(
      print(
        tme_correlations,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "TME-ADJUSTED MODEL",
  "------------------",
  paste0(
    "Design: ~ StromalScore_ESTlike + ImmuneScore_ESTlike + ",
    "EndothelialScore + Score_z"
  ),
  paste("Coefficient:", tme_coefficient),
  paste("Rows WITH signature:", nrow(tme_with_signature)),
  paste("Rows NO signature:", nrow(tme_no_signature)),
  paste("Fallback rows:", sum(tme_with_signature$gene_name_is_ensembl_fallback)),
  paste("FDR-significant NO signature:", count_fdr(tme_no_signature)),
  paste("Associated NO signature:", count_associated(tme_no_signature)),
  paste("Positive:", count_associated(tme_no_signature, "POSITIVE")),
  paste("Negative:", count_associated(tme_no_signature, "NEGATIVE")),
  paste("GSEA ranked symbols:", tme_rank$ranked_symbols),
  paste("Duplicate symbol rows removed:", tme_rank$duplicate_rows_removed),
  paste("Significant Hallmarks:", tme_gsea_info$n_significant),
  paste("EMT:", emt_text(tme_gsea_info$emt)),
  "",
  "TME-RESIDUALIZED MODEL",
  "----------------------",
  paste(
    "Residualization R-squared:",
    signif(residual_r2, 6)
  ),
  paste(
    "Residualization adjusted R-squared:",
    signif(residual_adjusted_r2, 6)
  ),
  paste("Coefficient:", residual_coefficient),
  paste("Rows WITH signature:", nrow(residual_with_signature)),
  paste("Rows NO signature:", nrow(residual_no_signature)),
  paste(
    "Fallback rows:",
    sum(residual_with_signature$gene_name_is_ensembl_fallback)
  ),
  paste(
    "FDR-significant NO signature:",
    count_fdr(residual_no_signature)
  ),
  paste(
    "Associated NO signature:",
    count_associated(residual_no_signature)
  ),
  paste(
    "Positive:",
    count_associated(
      residual_no_signature,
      "POSITIVE"
    )
  ),
  paste(
    "Negative:",
    count_associated(
      residual_no_signature,
      "NEGATIVE"
    )
  ),
  paste(
    "GSEA ranked symbols:",
    residual_rank$ranked_symbols
  ),
  paste(
    "Duplicate symbol rows removed:",
    residual_rank$duplicate_rows_removed
  ),
  paste(
    "Significant Hallmarks:",
    residual_gsea_info$n_significant
  ),
  paste(
    "EMT:",
    emt_text(residual_gsea_info$emt)
  ),
  "",
  "OUTPUTS",
  "-------",
  paste(
    "TME stromal correlation figure:",
    output_paths$figure_tme_stromal
  ),
  paste(
    "TME immune correlation figure:",
    output_paths$figure_tme_immune
  ),
  paste(
    "TME endothelial correlation figure:",
    output_paths$figure_tme_endothelial
  ),
  paste(
    "TME epithelial correlation figure:",
    output_paths$figure_tme_epithelial
  ),
  paste("EMT comparison figure:", output_paths$figure_emt_comparison),
  paste("TME-adjusted DDS:", output_paths$dds_tme),
  paste("Residualized DDS:", output_paths$dds_residual),
  "",
  "WARNINGS",
  "--------",
  if (length(warnings_log) == 0L) {
    "None"
  } else {
    paste(" -", unique(warnings_log))
  },
  "",
  "SESSION INFORMATION",
  "-------------------",
  paste(
    capture.output(sessionInfo()),
    collapse = "\n"
  ),
  "",
  paste("FINAL STATUS:", final_status)
)

writeLines(
  qc_lines,
  con = output_paths$qc_log,
  useBytes = TRUE
)

msg("QC log written:", output_paths$qc_log)
msg("Final status:", final_status)
