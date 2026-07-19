#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Continuous Hallmark GSEA
#
# Input rank: limma t_stat from the complete no-signature continuous results.
# Score: official 41-gene platelet-associated transcriptional score.
#
# This is a complementary continuous robustness analysis.
# All non-signature genes with finite t statistics are used for ranking.
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

message("[1/5] Loading inputs")

###############################################################################
# 1. Constants
###############################################################################

GSEA_SEED <- 123L
GSEA_NPROC <- 1L
MIN_SIZE <- 15L
MAX_SIZE <- 500L
FDR_THRESHOLD <- 0.05

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

input_limma <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables",
  "LIMMA",
  "LIMMA_continuous_41genes_no_signature_genes.csv"
)

input_signature <- file.path(
  RESOURCE_DIR,
  "platelet_associated_transcriptional_signature.tsv"
)

gsea_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables",
  "GSEA"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Logs"
)

output_full <- file.path(
  gsea_dir,
  "GSEA_continuous_no_signature_Hallmark_full.csv"
)

output_fdr <- file.path(
  gsea_dir,
  "GSEA_continuous_no_signature_Hallmark_FDR005.csv"
)

output_emt <- file.path(
  gsea_dir,
  "GSEA_continuous_no_signature_Hallmark_EMT_row.csv"
)

output_metadata <- file.path(
  gsea_dir,
  "GSEA_continuous_metadata.json"
)

output_qc <- file.path(
  logs_dir,
  "GSEA_continuous_QC.txt"
)

dir.create(
  gsea_dir,
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
      "Friedrich/GSE134051 continuous Hallmark GSEA QC",
      "=====================================",
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

  stop(
    msg,
    call. = FALSE
  )
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

  lower_names <- tolower(
    names(signature_df)
  )

  idx <- match(
    candidates,
    lower_names
  )

  idx <- idx[
    !is.na(idx)
  ]

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
  gene_col <- detect_gene_column(
    signature_df
  )

  genes <- trimws(
    as.character(
      signature_df[[gene_col]]
    )
  )

  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]

  duplicated_genes <- unique(
    genes[
      duplicated(genes)
    ]
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

coerce_signature_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  x_chr <- tolower(
    trimws(
      as.character(x)
    )
  )

  ifelse(
    x_chr %in% c(
      "true",
      "t",
      "1",
      "yes"
    ),
    TRUE,
    ifelse(
      x_chr %in% c(
        "false",
        "f",
        "0",
        "no"
      ),
      FALSE,
      NA
    )
  )
}

load_hallmark_pathways <- function() {
  hallmark <- tryCatch(
    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "H"
    ),
    error = function(e1) {
      tryCatch(
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = "H"
        ),
        error = function(e2) {
          fail_with_qc(
            paste0(
              "Could not load MSigDB Hallmark gene sets. ",
              "collection error: ",
              conditionMessage(e1),
              "; category error: ",
              conditionMessage(e2)
            )
          )
        }
      )
    }
  )

  required_cols <- c(
    "gs_name",
    "gene_symbol"
  )

  missing_cols <- setdiff(
    required_cols,
    names(hallmark)
  )

  if (length(missing_cols) > 0) {
    fail_with_qc(
      paste0(
        "msigdbr Hallmark table missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  hallmark <- hallmark[
    !is.na(hallmark$gs_name) &
      !is.na(hallmark$gene_symbol),
    ,
    drop = FALSE
  ]

  hallmark$gs_name <- as.character(
    hallmark$gs_name
  )

  hallmark$gene_symbol <- trimws(
    as.character(
      hallmark$gene_symbol
    )
  )

  hallmark <- hallmark[
    hallmark$gene_symbol != "",
    ,
    drop = FALSE
  ]

  lapply(
    split(
      hallmark$gene_symbol,
      hallmark$gs_name
    ),
    unique
  )
}

format_gsea_result <- function(fg) {
  fg <- as.data.frame(
    fg,
    stringsAsFactors = FALSE
  )

  if (nrow(fg) == 0) {
    fail_with_qc(
      "fgsea returned an empty result table."
    )
  }

  if (!"log2err" %in% names(fg)) {
    fg$log2err <- NA_real_
  }

  if (!"nMoreExtreme" %in% names(fg)) {
    fg$nMoreExtreme <- NA_real_
  }

  if ("leadingEdge" %in% names(fg)) {
    fg$leadingEdge <- vapply(
      fg$leadingEdge,
      function(x) {
        paste(
          as.character(x),
          collapse = ";"
        )
      },
      character(1)
    )
  } else {
    fg$leadingEdge <- NA_character_
  }

  fg$direction <- ifelse(
    fg$NES > 0,
    "Higher_Score_z",
    ifelse(
      fg$NES < 0,
      "Lower_Score_z",
      "Zero_NES"
    )
  )

  fg$interpretation <- ifelse(
    fg$NES > 0,
    paste0(
      "enriched with higher platelet-associated ",
      "transcriptional score"
    ),
    ifelse(
      fg$NES < 0,
      paste0(
        "enriched with lower platelet-associated ",
        "transcriptional score"
      ),
      "no enrichment direction"
    )
  )

  keep_cols <- c(
    "pathway",
    "pval",
    "padj",
    "log2err",
    "ES",
    "NES",
    "nMoreExtreme",
    "size",
    "leadingEdge",
    "direction",
    "interpretation"
  )

  missing_cols <- setdiff(
    keep_cols,
    names(fg)
  )

  for (col in missing_cols) {
    fg[[col]] <- NA
  }

  fg <- fg[
    ,
    keep_cols,
    drop = FALSE
  ]

  fg[
    order(
      -fg$NES,
      fg$padj,
      fg$pathway
    ),
    ,
    drop = FALSE
  ]
}

write_qc_report <- function(qc) {
  sink(output_qc)

  cat("Friedrich/GSE134051 continuous Hallmark GSEA QC\n")
  cat("=====================================\n\n")

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
    "  LIMMA no-signature input: ",
    qc$inputs$limma,
    "\n",
    sep = ""
  )

  cat(
    "  Signature input: ",
    qc$inputs$signature,
    "\n\n",
    sep = ""
  )

  cat("2. Ranking\n")
  cat("  Ranking variable: t_stat\n")
  cat("  Ranking fallback: none\n")

  cat(
    "  Input rows: ",
    qc$ranking$n_input_rows,
    "\n",
    sep = ""
  )

  cat(
    "  Ranked genes: ",
    qc$ranking$n_ranked_genes,
    "\n",
    sep = ""
  )

  cat(
    "  Non-finite rows removed: ",
    qc$ranking$n_nonfinite_removed,
    "\n",
    sep = ""
  )

  cat(
    "  Duplicate gene names: ",
    qc$ranking$n_duplicate_genes,
    "\n",
    sep = ""
  )

  cat(
    "  Duplicate t-statistic values: ",
    qc$ranking$n_duplicate_tstat_values,
    "\n",
    sep = ""
  )

  cat(
    "  Signature genes in input: ",
    qc$ranking$n_signature_in_input,
    "\n",
    sep = ""
  )

  cat(
    "  t_stat minimum: ",
    qc$ranking$t_stat_min,
    "\n",
    sep = ""
  )

  cat(
    "  t_stat median: ",
    qc$ranking$t_stat_median,
    "\n",
    sep = ""
  )

  cat(
    "  t_stat maximum: ",
    qc$ranking$t_stat_max,
    "\n\n",
    sep = ""
  )

  cat("3. GSEA\n")
  cat("  Method: fgsea preranked\n")
  cat("  Collection: MSigDB Hallmark\n")

  cat(
    "  minSize: ",
    MIN_SIZE,
    "\n",
    sep = ""
  )

  cat(
    "  maxSize: ",
    MAX_SIZE,
    "\n",
    sep = ""
  )

  cat("  eps: 0\n")

  cat(
    "  Seed: ",
    GSEA_SEED,
    "\n",
    sep = ""
  )

  cat(
    "  nproc: ",
    GSEA_NPROC,
    "\n",
    sep = ""
  )

  cat(
    "  Pathways tested: ",
    qc$gsea$n_pathways_tested,
    "\n",
    sep = ""
  )

  cat(
    "  Pathways FDR < 0.05: ",
    qc$gsea$n_pathways_fdr005,
    "\n",
    sep = ""
  )

  cat(
    "  Positive NES: ",
    qc$gsea$n_positive_nes,
    "\n",
    sep = ""
  )

  cat(
    "  Negative NES: ",
    qc$gsea$n_negative_nes,
    "\n\n",
    sep = ""
  )

  cat("4. EMT\n")

  cat(
    "  EMT NES: ",
    qc$emt$NES,
    "\n",
    sep = ""
  )

  cat(
    "  EMT padj: ",
    qc$emt$padj,
    "\n",
    sep = ""
  )

  cat(
    "  EMT direction: ",
    qc$emt$direction,
    "\n\n",
    sep = ""
  )

  cat("5. Reproducibility\n")

  cat(
    "  R version: ",
    R.version.string,
    "\n",
    sep = ""
  )

  cat(
    "  fgsea version: ",
    as.character(
      utils::packageVersion("fgsea")
    ),
    "\n",
    sep = ""
  )

  cat(
    "  msigdbr version: ",
    as.character(
      utils::packageVersion("msigdbr")
    ),
    "\n\n",
    sep = ""
  )

  cat("6. Outputs\n")

  cat(
    "  Full table: ",
    output_full,
    "\n",
    sep = ""
  )

  cat(
    "  FDR005 table: ",
    output_fdr,
    "\n",
    sep = ""
  )

  cat(
    "  EMT row: ",
    output_emt,
    "\n",
    sep = ""
  )

  cat(
    "  Metadata: ",
    output_metadata,
    "\n\n",
    sep = ""
  )

  cat("7. Warnings\n")

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
    "fgsea",
    "msigdbr",
    "jsonlite"
  )
)

for (
  fp in c(
    input_limma,
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

limma_df <- read_csv_base(
  input_limma
)

signature_df <- read_canonical_platelet_signature()

signature_genes <- extract_signature_genes(
  signature_df
)

required_cols <- c(
  "gene_name",
  "t_stat",
  "is_signature_gene"
)

missing_cols <- setdiff(
  required_cols,
  names(limma_df)
)

if (length(missing_cols) > 0) {
  fail_with_qc(
    paste0(
      "LIMMA no-signature input missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

###############################################################################
# 5. Prepare ranking
###############################################################################

message("[2/5] Preparing t_stat ranking")

limma_df$gene_name <- trimws(
  as.character(
    limma_df$gene_name
  )
)

limma_df$t_stat <- suppressWarnings(
  as.numeric(
    limma_df$t_stat
  )
)

limma_df$is_signature_gene <- coerce_signature_flag(
  limma_df$is_signature_gene
)

if (
  any(
    is.na(limma_df$gene_name) |
    limma_df$gene_name == ""
  )
) {
  fail_with_qc(
    "LIMMA no-signature input contains missing or empty gene_name values."
  )
}

n_duplicate_genes <- sum(
  duplicated(limma_df$gene_name)
)

if (n_duplicate_genes > 0) {
  duplicated_genes <- unique(
    limma_df$gene_name[
      duplicated(limma_df$gene_name)
    ]
  )

  fail_with_qc(
    paste0(
      "Duplicate gene_name values detected in LIMMA input: ",
      paste(
        utils::head(duplicated_genes, 25),
        collapse = ", "
      ),
      if (length(duplicated_genes) > 25) {
        " ..."
      } else {
        ""
      }
    )
  )
}

if (
  any(
    limma_df$is_signature_gene %in% TRUE,
    na.rm = TRUE
  )
) {
  fail_with_qc(
    paste0(
      "No-signature input contains rows flagged as signature genes: ",
      sum(
        limma_df$is_signature_gene %in% TRUE,
        na.rm = TRUE
      )
    )
  )
}

n_signature_in_input <- sum(
  limma_df$gene_name %in%
    signature_genes
)

if (n_signature_in_input > 0) {
  fail_with_qc(
    paste0(
      "No-signature input contains official signature genes: ",
      n_signature_in_input
    )
  )
}

n_input_rows <- nrow(
  limma_df
)

rank_df <- limma_df[
  is.finite(limma_df$t_stat),
  c(
    "gene_name",
    "t_stat"
  ),
  drop = FALSE
]

n_nonfinite_removed <- (
  n_input_rows -
    nrow(rank_df)
)

if (n_nonfinite_removed > 0) {
  add_warning(
    paste0(
      "Removed rows with non-finite t_stat: ",
      n_nonfinite_removed
    )
  )
}

if (nrow(rank_df) < 5000) {
  fail_with_qc(
    paste0(
      "Too few finite ranked genes for GSEA: ",
      nrow(rank_df)
    )
  )
}

n_duplicate_tstat_values <- sum(
  duplicated(rank_df$t_stat)
)

ranks <- stats::setNames(
  rank_df$t_stat,
  rank_df$gene_name
)

ranks <- sort(
  ranks,
  decreasing = TRUE
)

ranking_qc <- list(
  n_input_rows = n_input_rows,
  n_ranked_genes = length(ranks),
  n_nonfinite_removed = n_nonfinite_removed,
  n_duplicate_genes = n_duplicate_genes,
  n_duplicate_tstat_values =
    n_duplicate_tstat_values,
  n_signature_in_input =
    n_signature_in_input,
  t_stat_min = min(
    ranks,
    na.rm = TRUE
  ),
  t_stat_median = stats::median(
    ranks,
    na.rm = TRUE
  ),
  t_stat_max = max(
    ranks,
    na.rm = TRUE
  )
)

###############################################################################
# 6. Load Hallmark collection
###############################################################################

message("[3/5] Running fgsea preranked Hallmark analysis")

hallmark_raw <- load_hallmark_pathways()

pathways <- lapply(
  hallmark_raw,
  function(gs) {
    unique(
      gs[
        gs %in% names(ranks)
      ]
    )
  }
)

pathway_sizes <- vapply(
  pathways,
  length,
  integer(1)
)

pathways <- pathways[
  pathway_sizes >= MIN_SIZE &
    pathway_sizes <= MAX_SIZE
]

n_pathways_tested <- length(
  pathways
)

if (n_pathways_tested < 40) {
  fail_with_qc(
    paste0(
      "Too few Hallmark pathways available after overlap and size filtering: ",
      n_pathways_tested
    )
  )
}

###############################################################################
# 7. Run deterministic GSEA
###############################################################################

set.seed(
  GSEA_SEED
)

fg <- tryCatch(
  fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    eps = 0,
    nproc = GSEA_NPROC
  ),
  error = function(e) {
    fail_with_qc(
      paste0(
        "fgsea failed: ",
        conditionMessage(e)
      )
    )
  }
)

gsea_full <- format_gsea_result(
  fg
)

gsea_fdr <- gsea_full[
  !is.na(gsea_full$padj) &
    gsea_full$padj < FDR_THRESHOLD,
  ,
  drop = FALSE
]

if (nrow(gsea_fdr) == 0) {
  add_warning(
    "No Hallmark pathways were significant at FDR < 0.05."
  )
}

emt_pathway <- paste0(
  "HALLMARK_EPITHELIAL_",
  "MESENCHYMAL_TRANSITION"
)

emt_row <- gsea_full[
  gsea_full$pathway == emt_pathway,
  ,
  drop = FALSE
]

if (nrow(emt_row) != 1) {
  fail_with_qc(
    paste0(
      emt_pathway,
      " was not found exactly once in GSEA output."
    )
  )
}

emt_direction <- as.character(
  emt_row$direction[1]
)

###############################################################################
# 8. Metadata and QC
###############################################################################

qc_status <- if (
  length(warnings_vec) > 0
) {
  "WARNING"
} else {
  "PASS"
}

metadata <- list(
  date_time = as.character(
    Sys.time()
  ),
  project_dir = project_dir,
  project = "friedrich_gse134051",
  cohort = "Friedrich/GSE134051 primary prostate tumors",
  analysis = "continuous Hallmark GSEA",
  analytic_hierarchy =
    "complementary robustness analysis",
  score = paste0(
    "official 41-gene platelet-associated ",
    "transcriptional score"
  ),
  method = "fgsea preranked",
  collection = "MSigDB Hallmark",
  model = "expression ~ Score_z",
  ranking_variable = "t_stat",
  ranking_fallback = "none",
  inputs = list(
    limma_no_signature = input_limma,
    signature = input_signature
  ),
  outputs = list(
    full = output_full,
    fdr005 = output_fdr,
    emt = output_emt,
    metadata = output_metadata,
    qc = output_qc
  ),
  parameters = list(
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    eps = 0,
    FDR_threshold = FDR_THRESHOLD
  ),
  ranking = list(
    n_input_rows = n_input_rows,
    n_ranked_genes = length(ranks),
    n_nonfinite_removed =
      n_nonfinite_removed,
    n_duplicate_tstat_values =
      n_duplicate_tstat_values
  ),
  results = list(
    n_pathways_tested =
      nrow(gsea_full),
    n_pathways_FDR005 =
      nrow(gsea_fdr),
    EMT_NES =
      as.numeric(emt_row$NES[1]),
    EMT_padj =
      as.numeric(emt_row$padj[1]),
    EMT_direction =
      emt_direction
  ),
  interpretation = list(
    positive_NES = paste0(
      "enriched with higher platelet-associated ",
      "transcriptional score"
    ),
    negative_NES = paste0(
      "enriched with lower platelet-associated ",
      "transcriptional score"
    )
  ),
  reproducibility = list(
    seed = GSEA_SEED,
    nproc = GSEA_NPROC,
    R_version = R.version.string,
    fgsea_version = as.character(
      utils::packageVersion("fgsea")
    ),
    msigdbr_version = as.character(
      utils::packageVersion("msigdbr")
    )
  ),
  qc_status = qc_status
)

qc <- list(
  inputs = list(
    limma = input_limma,
    signature = input_signature
  ),
  ranking = ranking_qc,
  gsea = list(
    n_pathways_tested =
      nrow(gsea_full),
    n_pathways_fdr005 =
      nrow(gsea_fdr),
    n_positive_nes = sum(
      gsea_full$NES > 0,
      na.rm = TRUE
    ),
    n_negative_nes = sum(
      gsea_full$NES < 0,
      na.rm = TRUE
    )
  ),
  emt = list(
    NES = as.numeric(
      emt_row$NES[1]
    ),
    padj = as.numeric(
      emt_row$padj[1]
    ),
    direction = emt_direction
  ),
  warnings = warnings_vec,
  status = qc_status
)

###############################################################################
# 9. Write outputs
###############################################################################

message("[4/5] Writing outputs")

write_csv_base(
  gsea_full,
  output_full
)

write_csv_base(
  gsea_fdr,
  output_fdr
)

write_csv_base(
  emt_row,
  output_emt
)

write_json(
  metadata,
  output_metadata
)

write_qc_report(
  qc
)

message(
  "[5/5] Final status: ",
  qc_status
)
