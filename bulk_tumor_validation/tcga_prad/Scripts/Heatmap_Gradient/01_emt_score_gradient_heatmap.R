#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD platelet-associated transcriptional score module
# Curated EMT score-gradient heatmap
#
# Script:
# Scripts/Heatmap_Gradient/01_emt_score_gradient_heatmap.R
#
# Analytical scope:
#   1. Use the canonical 497-patient platelet-score cohort.
#   2. Order samples by the official continuous Score_z.
#   3. Divide the ordered cohort into five balanced descriptive bins.
#   4. Calculate mean VST expression per bin for a curated EMT marker panel.
#   5. Row-standardize the five-bin expression means.
#   6. Generate the main epithelial/mesenchymal score-gradient heatmap.
#
# Important:
#   This is a descriptive five-bin visualization of the Score_z gradient.
#   It is not a continuous regression model, does not perform differential-
#   expression testing across bins and is not a Q1-versus-Q4 analysis.
#
#   The marker panel is a curated epithelial/mesenchymal set and is not the
#   MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION gene set.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

################################################################################
# 1. Required packages
################################################################################

required_packages <- c(
  "DESeq2",
  "matrixStats",
  "ComplexHeatmap",
  "circlize",
  "SummarizedExperiment"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
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
# 2. Constants and locked style
################################################################################

EXPECTED_SCORE_SAMPLES <- 497L
NUMBER_OF_BINS <- 5L

BASE_FAMILY <- "Helvetica"
TEXT_SIZE <- 5

COL_UP <- "#E22B27"
COL_DOWN <- "#3B76B7"

COL_EPI <- "#3BB247"
COL_MES <- "#D7DF23"

################################################################################
# 3. Helpers
################################################################################

msg <- function(...) {
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "|",
    ...,
    "\n"
  )
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

warnings_log <- character()
generated_figures <- character()
generated_tables <- character()

add_warning <- function(...) {
  warning_text <- paste0(...)

  warnings_log <<- unique(
    c(warnings_log, warning_text)
  )

  warning(
    warning_text,
    call. = FALSE
  )
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    fail(
      "Missing ",
      label,
      ": ",
      path
    )
  }
}

read_csv_required <- function(path, label) {
  require_file(path, label)

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

require_columns <- function(df, columns, label) {
  missing_columns <- setdiff(
    columns,
    colnames(df)
  )

  if (length(missing_columns) > 0L) {
    fail(
      label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

write_table <- function(
    x,
    path,
    row_names = FALSE
) {
  utils::write.csv(
    x,
    path,
    row.names = row_names
  )

  generated_tables <<- unique(
    c(generated_tables, path)
  )

  invisible(path)
}

load_counts_matrix <- function(path) {
  counts <- readRDS(path)

  if (inherits(counts, "SummarizedExperiment")) {
    assay_names <- SummarizedExperiment::assayNames(
      counts
    )

    if (length(assay_names) == 0L) {
      fail(
        "Counts SummarizedExperiment contains no assays."
      )
    }

    preferred_assays <- c(
      "unstranded",
      "counts",
      "raw_counts"
    )

    assay_name <- preferred_assays[
      preferred_assays %in% assay_names
    ][1L]

    if (
      length(assay_name) == 0L ||
      is.na(assay_name)
    ) {
      assay_name <- assay_names[1L]
    }

    msg(
      "Counts object is a SummarizedExperiment; using assay:",
      assay_name
    )

    counts <- SummarizedExperiment::assay(
      counts,
      assay_name
    )
  }

  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }

  if (
    is.null(dim(counts)) ||
    length(dim(counts)) != 2L
  ) {
    fail(
      "Counts input is not a two-dimensional matrix-like object."
    )
  }

  if (
    is.null(rownames(counts)) ||
    anyNA(rownames(counts)) ||
    any(rownames(counts) == "")
  ) {
    fail(
      "Counts matrix must have complete gene row identifiers."
    )
  }

  if (
    is.null(colnames(counts)) ||
    anyNA(colnames(counts)) ||
    any(colnames(counts) == "")
  ) {
    fail(
      "Counts matrix must have complete sample identifiers."
    )
  }

  if (anyDuplicated(rownames(counts))) {
    fail(
      "Counts matrix contains duplicated gene row identifiers."
    )
  }

  if (anyDuplicated(colnames(counts))) {
    fail(
      "Counts matrix contains duplicated sample identifiers."
    )
  }

  if (anyNA(counts)) {
    fail(
      "Counts matrix contains missing values."
    )
  }

  if (any(counts < 0, na.rm = TRUE)) {
    fail(
      "Counts matrix contains negative values."
    )
  }

  if (!is.integer(counts)) {
    non_integer_values <- abs(
      counts - round(counts)
    ) > 1e-6

    if (any(non_integer_values, na.rm = TRUE)) {
      fail(
        "Counts matrix does not contain integer-like raw counts."
      )
    }

    storage.mode(counts) <- "integer"
  }

  counts
}

prepare_gene_map <- function(map_df) {
  require_columns(
    map_df,
    c(
      "ensg_version",
      "gene_name"
    ),
    "gene map"
  )

  map_df$ensg_version <- trimws(
    as.character(map_df$ensg_version)
  )

  map_df$gene_name <- toupper(
    trimws(as.character(map_df$gene_name))
  )

  map_df <- map_df[
    !is.na(map_df$ensg_version) &
      map_df$ensg_version != "",
    ,
    drop = FALSE
  ]

  duplicated_ensg <- unique(
    map_df$ensg_version[
      duplicated(map_df$ensg_version)
    ]
  )

  if (length(duplicated_ensg) > 0L) {
    conflicting_ids <- duplicated_ensg[
      vapply(
        duplicated_ensg,
        function(ensembl_id) {
          symbols <- unique(
            map_df$gene_name[
              map_df$ensg_version == ensembl_id &
                !is.na(map_df$gene_name) &
                map_df$gene_name != ""
            ]
          )

          length(symbols) > 1L
        },
        logical(1)
      )
    ]

    if (length(conflicting_ids) > 0L) {
      fail(
        "Gene map contains ENSG identifiers mapped to conflicting symbols. ",
        "Examples: ",
        paste(
          utils::head(
            conflicting_ids,
            8L
          ),
          collapse = ", "
        )
      )
    }

    map_df <- map_df[
      !duplicated(map_df$ensg_version),
      ,
      drop = FALSE
    ]
  }

  map_df
}

align_exact_samples <- function(
    reference_ids,
    target_ids,
    target_label
) {
  reference_ids <- trimws(
    as.character(reference_ids)
  )

  target_ids <- trimws(
    as.character(target_ids)
  )

  if (
    anyNA(reference_ids) ||
    any(reference_ids == "")
  ) {
    fail(
      "Canonical sample_id vector contains missing or empty values."
    )
  }

  if (
    anyNA(target_ids) ||
    any(target_ids == "")
  ) {
    fail(
      target_label,
      " contains missing or empty sample identifiers."
    )
  }

  if (anyDuplicated(reference_ids)) {
    fail(
      "Canonical sample_id vector contains duplicated values."
    )
  }

  if (anyDuplicated(target_ids)) {
    fail(
      target_label,
      " contains duplicated sample identifiers."
    )
  }

  alignment_index <- match(
    reference_ids,
    target_ids
  )

  if (anyNA(alignment_index)) {
    missing_ids <- reference_ids[
      is.na(alignment_index)
    ]

    fail(
      "Exact sample alignment failed for ",
      target_label,
      ". Missing examples: ",
      paste(
        utils::head(
          missing_ids,
          8L
        ),
        collapse = ", "
      ),
      ". Barcode normalization or truncation is not allowed."
    )
  }

  alignment_index
}

collapse_by_max_variance <- function(
    expression_matrix
) {
  row_variances <- matrixStats::rowVars(
    expression_matrix
  )

  symbol_indices <- split(
    seq_len(nrow(expression_matrix)),
    rownames(expression_matrix)
  )

  selected_indices <- vapply(
    symbol_indices,
    function(indices) {
      indices[
        which.max(row_variances[indices])
      ]
    },
    integer(1)
  )

  collapsed_matrix <- expression_matrix[
    selected_indices,
    ,
    drop = FALSE
  ]

  rownames(collapsed_matrix) <- names(
    selected_indices
  )

  collapsed_matrix
}

calculate_row_zscores <- function(matrix_input) {
  row_standard_deviations <- matrixStats::rowSds(
    matrix_input
  )

  invalid_rows <- !is.finite(
    row_standard_deviations
  ) |
    row_standard_deviations <= 0

  if (any(invalid_rows)) {
    fail(
      "Cannot calculate row Z-scores for zero-variance genes: ",
      paste(
        rownames(matrix_input)[invalid_rows],
        collapse = ", "
      )
    )
  }

  centered_matrix <- sweep(
    matrix_input,
    1L,
    matrixStats::rowMeans2(matrix_input),
    FUN = "-"
  )

  sweep(
    centered_matrix,
    1L,
    row_standard_deviations,
    FUN = "/"
  )
}

make_balanced_bins <- function(
    n,
    number_of_bins
) {
  if (n < number_of_bins) {
    fail(
      "Cannot divide ",
      n,
      " samples into ",
      number_of_bins,
      " bins."
    )
  }

  base_size <- n %/% number_of_bins
  remainder <- n %% number_of_bins

  bin_sizes <- rep(
    base_size,
    number_of_bins
  )

  if (remainder > 0L) {
    bin_sizes[
      seq_len(remainder)
    ] <- bin_sizes[
      seq_len(remainder)
    ] + 1L
  }

  rep(
    seq_len(number_of_bins),
    times = bin_sizes
  )
}

################################################################################
# 4. Project paths
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

msg(
  "Project directory:",
  project_dir
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
  )
)

output_dirs <- list(
  figures = file.path(
    GENERATED_FIGURES_DIR, "Heatmap_Gradient"
  ),

  tables = file.path(
    GENERATED_RESULTS_DIR, "Heatmap_Gradient",
    "Tables"
  ),

  logs = file.path(
    GENERATED_RESULTS_DIR, "Heatmap_Gradient",
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
  figure = file.path(
    output_dirs$figures,
    "EMT_score_gradient_5bins.pdf"
  ),

  marker_coverage = file.path(
    output_dirs$tables,
    "EMT_marker_coverage.csv"
  ),

  direction_validation = file.path(
    output_dirs$tables,
    "EMT_score_direction_validation.csv"
  ),

  bin_summary = file.path(
    output_dirs$tables,
    "EMT_score_gradient_bin_summary.csv"
  ),

  bin_means = file.path(
    output_dirs$tables,
    "EMT_score_gradient_bin_means.csv"
  ),

  row_zscores = file.path(
    output_dirs$tables,
    "EMT_score_gradient_row_zscores.csv"
  ),

  qc_log = file.path(
    output_dirs$logs,
    "01_emt_score_gradient_heatmap_QC.txt"
  )
)

for (input_name in names(input_paths)) {
  require_file(
    input_paths[[input_name]],
    input_name
  )
}

################################################################################
# 5. Load canonical score cohort
################################################################################

msg(
  "Loading canonical platelet-score cohort"
)

score_df <- read_csv_required(
  input_paths$score,
  "canonical platelet-score table"
)

require_columns(
  score_df,
  c(
    "sample_id",
    "patient_id",
    "Score_raw",
    "Score_z"
  ),
  "canonical platelet-score table"
)

score_df$sample_id <- trimws(
  as.character(score_df$sample_id)
)

score_df$patient_id <- trimws(
  as.character(score_df$patient_id)
)

score_df$Score_raw <- suppressWarnings(
  as.numeric(score_df$Score_raw)
)

score_df$Score_z <- suppressWarnings(
  as.numeric(score_df$Score_z)
)

if (nrow(score_df) != EXPECTED_SCORE_SAMPLES) {
  fail(
    "Expected ",
    EXPECTED_SCORE_SAMPLES,
    " canonical score rows; observed ",
    nrow(score_df),
    "."
  )
}

if (
  anyNA(score_df$sample_id) ||
  any(score_df$sample_id == "") ||
  anyDuplicated(score_df$sample_id)
) {
  fail(
    "Canonical score table must contain unique and complete sample_id values."
  )
}

if (
  anyNA(score_df$patient_id) ||
  any(score_df$patient_id == "") ||
  anyDuplicated(score_df$patient_id)
) {
  fail(
    "Canonical score table must contain one unique patient_id per row."
  )
}

if (
  anyNA(score_df$Score_z) ||
  any(!is.finite(score_df$Score_z))
) {
  fail(
    "Score_z contains missing or non-finite values."
  )
}

################################################################################
# 6. Load and align raw counts
################################################################################

msg(
  "Loading raw tumor count matrix"
)

counts <- load_counts_matrix(
  input_paths$counts
)

counts_dimensions_before_alignment <- dim(
  counts
)

raw_gene_count <- nrow(
  counts
)

count_alignment_index <- align_exact_samples(
  reference_ids = score_df$sample_id,
  target_ids = colnames(counts),
  target_label = "raw count matrix"
)

counts_aligned <- counts[
  ,
  count_alignment_index,
  drop = FALSE
]

if (
  !identical(
    colnames(counts_aligned),
    score_df$sample_id
  )
) {
  fail(
    "Aligned count matrix does not preserve canonical sample order."
  )
}

rm(counts)
invisible(gc(verbose = FALSE))

################################################################################
# 7. Gene map and VST expression
################################################################################

msg(
  "Loading canonical gene map"
)

gene_map <- prepare_gene_map(
  read_csv_required(
    input_paths$gene_map,
    "gene map"
  )
)

msg(
  "Applying prefilter: rowSums(counts >= 10) >= 10"
)

genes_before_prefilter <- nrow(
  counts_aligned
)

prefilter_keep <- rowSums(
  counts_aligned >= 10
) >= 10

counts_prefiltered <- counts_aligned[
  prefilter_keep,
  ,
  drop = FALSE
]

genes_after_prefilter <- nrow(
  counts_prefiltered
)

rm(counts_aligned)
invisible(gc(verbose = FALSE))

coldata <- data.frame(
  Score_z = score_df$Score_z,
  row.names = score_df$sample_id,
  stringsAsFactors = FALSE
)

if (
  !identical(
    rownames(coldata),
    colnames(counts_prefiltered)
  )
) {
  fail(
    "DESeq2 colData and count matrix are not identically aligned."
  )
}

msg(
  "Calculating blind VST expression matrix"
)

dds_vst <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_prefiltered,
  colData = coldata,
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

################################################################################
# 8. Map expression to unique gene symbols
################################################################################

gene_symbols <- gene_map$gene_name[
  match(
    rownames(vst_matrix),
    gene_map$ensg_version
  )
]

mapped_rows <- !is.na(gene_symbols) &
  gene_symbols != ""

expression_by_symbol <- vst_matrix[
  mapped_rows,
  ,
  drop = FALSE
]

rownames(expression_by_symbol) <- gene_symbols[
  mapped_rows
]

expression_by_symbol <- collapse_by_max_variance(
  expression_by_symbol
)

################################################################################
# 9. Curated EMT marker panel
################################################################################

EPITHELIAL_GENES <- c(
  "EPCAM",
  "CDH1",
  "OCLN",
  "DSP",
  "CLDN3",
  "CLDN4",
  "KRT8",
  "KRT18"
)

MESENCHYMAL_GENES <- c(
  "SNAI1",
  "SNAI2",
  "ZEB1",
  "ZEB2",
  "FOXC2",
  "PRRX1",
  "SMAD3",
  "VIM",
  "CDH2",
  "FN1",
  "ITGA5",
  "LGALS3",
  "PLAUR",
  "MMP2",
  "MMP9",
  "SERPINE1"
)

EMT_MARKER_GENES <- unique(
  c(
    EPITHELIAL_GENES,
    MESENCHYMAL_GENES
  )
)

epithelial_present <- intersect(
  EPITHELIAL_GENES,
  rownames(expression_by_symbol)
)

mesenchymal_present <- intersect(
  MESENCHYMAL_GENES,
  rownames(expression_by_symbol)
)

genes_present <- c(
  epithelial_present,
  mesenchymal_present
)

if (length(genes_present) < 8L) {
  fail(
    "Too few curated EMT marker genes are present: ",
    length(genes_present),
    "."
  )
}

marker_coverage <- data.frame(
  marker_class = c(
    "Epithelial",
    "Mesenchymal",
    "All"
  ),

  n_requested = c(
    length(EPITHELIAL_GENES),
    length(MESENCHYMAL_GENES),
    length(EMT_MARKER_GENES)
  ),

  n_present = c(
    length(epithelial_present),
    length(mesenchymal_present),
    length(genes_present)
  ),

  present_genes = c(
    paste(
      epithelial_present,
      collapse = ";"
    ),

    paste(
      mesenchymal_present,
      collapse = ";"
    ),

    paste(
      genes_present,
      collapse = ";"
    )
  ),

  missing_genes = c(
    paste(
      setdiff(
        EPITHELIAL_GENES,
        epithelial_present
      ),
      collapse = ";"
    ),

    paste(
      setdiff(
        MESENCHYMAL_GENES,
        mesenchymal_present
      ),
      collapse = ";"
    ),

    paste(
      setdiff(
        EMT_MARKER_GENES,
        genes_present
      ),
      collapse = ";"
    )
  ),

  stringsAsFactors = FALSE
)

write_table(
  marker_coverage,
  output_paths$marker_coverage
)

expression_emt <- expression_by_symbol[
  genes_present,
  ,
  drop = FALSE
]

################################################################################
# 10. Descriptive direction validation
################################################################################

msg(
  "Calculating descriptive marker-score correlations"
)

score_vector <- score_df$Score_z
names(score_vector) <- score_df$sample_id

calculate_marker_correlations <- function(
    genes,
    marker_class
) {
  correlations <- vapply(
    genes,
    function(gene_name) {
      suppressWarnings(
        stats::cor(
          score_vector,
          expression_emt[
            gene_name,
            names(score_vector)
          ],
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
    },
    numeric(1)
  )

  data.frame(
    gene = genes,
    marker_class = marker_class,
    spearman_rho = as.numeric(correlations),
    stringsAsFactors = FALSE
  )
}

marker_correlations <- rbind(
  calculate_marker_correlations(
    epithelial_present,
    "Epithelial"
  ),

  calculate_marker_correlations(
    mesenchymal_present,
    "Mesenchymal"
  )
)

direction_validation <- do.call(
  rbind,
  lapply(
    c(
      "Epithelial",
      "Mesenchymal"
    ),
    function(marker_class) {
      class_values <- marker_correlations$spearman_rho[
        marker_correlations$marker_class == marker_class
      ]

      data.frame(
        score = "Score_z",
        marker_class = marker_class,
        n_genes = length(class_values),
        mean_spearman_rho = mean(
          class_values,
          na.rm = TRUE
        ),
        median_spearman_rho = stats::median(
          class_values,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_table(
  direction_validation,
  output_paths$direction_validation
)

################################################################################
# 11. Build five-bin score-gradient matrices
################################################################################

msg(
  "Building five balanced Score_z bins"
)

sample_order <- order(
  score_df$Score_z,
  score_df$sample_id
)

gradient_df <- data.frame(
  sample_id = score_df$sample_id[
    sample_order
  ],

  patient_id = score_df$patient_id[
    sample_order
  ],

  Score_z = score_df$Score_z[
    sample_order
  ],

  stringsAsFactors = FALSE
)

gradient_df$bin_number <- make_balanced_bins(
  n = nrow(gradient_df),
  number_of_bins = NUMBER_OF_BINS
)

gradient_df$bin <- factor(
  paste0(
    "B",
    gradient_df$bin_number
  ),
  levels = paste0(
    "B",
    seq_len(NUMBER_OF_BINS)
  )
)

bin_summary <- do.call(
  rbind,
  lapply(
    levels(gradient_df$bin),
    function(bin_label) {
      bin_rows <- gradient_df[
        gradient_df$bin == bin_label,
        ,
        drop = FALSE
      ]

      data.frame(
        bin = bin_label,
        n = nrow(bin_rows),
        score_min = min(bin_rows$Score_z),
        score_median = stats::median(
          bin_rows$Score_z
        ),
        score_max = max(bin_rows$Score_z),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_table(
  bin_summary,
  output_paths$bin_summary
)

bin_mean_matrix <- sapply(
  levels(gradient_df$bin),
  function(bin_label) {
    bin_samples <- gradient_df$sample_id[
      gradient_df$bin == bin_label
    ]

    matrixStats::rowMeans2(
      expression_emt[
        ,
        bin_samples,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
)

colnames(bin_mean_matrix) <- levels(
  gradient_df$bin
)

write_table(
  data.frame(
    gene = rownames(bin_mean_matrix),
    bin_mean_matrix,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  output_paths$bin_means
)

row_zscore_matrix <- calculate_row_zscores(
  bin_mean_matrix
)

gene_order <- c(
  epithelial_present,
  mesenchymal_present
)

row_zscore_matrix <- row_zscore_matrix[
  gene_order,
  ,
  drop = FALSE
]

write_table(
  data.frame(
    gene = rownames(row_zscore_matrix),
    row_zscore_matrix,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  output_paths$row_zscores
)

marker_class <- factor(
  c(
    rep(
      "Epithelial",
      length(epithelial_present)
    ),

    rep(
      "Mesenchymal",
      length(mesenchymal_present)
    )
  ),
  levels = c(
    "Epithelial",
    "Mesenchymal"
  )
)

names(marker_class) <- gene_order

################################################################################
# 12. Heatmap
################################################################################

msg(
  "Generating curated EMT score-gradient heatmap"
)

heatmap_limit <- as.numeric(
  stats::quantile(
    abs(
      as.numeric(row_zscore_matrix)
    ),
    probs = 0.98,
    na.rm = TRUE
  )
)

if (
  !is.finite(heatmap_limit) ||
  heatmap_limit <= 0
) {
  heatmap_limit <- 2
}

heatmap_color_function <- circlize::colorRamp2(
  c(
    -heatmap_limit,
    0,
    heatmap_limit
  ),
  c(
    COL_DOWN,
    "white",
    COL_UP
  )
)

marker_class_colors <- c(
  Epithelial = COL_EPI,
  Mesenchymal = COL_MES
)

ComplexHeatmap::ht_opt(
  ROW_ANNO_PADDING = grid::unit(
    0.4,
    "mm"
  )
)

ComplexHeatmap::ht_opt(
  HEATMAP_LEGEND_PADDING = grid::unit(
    0.75,
    "mm"
  )
)

ComplexHeatmap::ht_opt(
  ANNOTATION_LEGEND_PADDING = grid::unit(
    0.75,
    "mm"
  )
)

ComplexHeatmap::ht_opt(
  TITLE_PADDING = grid::unit(
    c(
      0.7,
      0.7
    ),
    "mm"
  )
)

row_annotation <- ComplexHeatmap::rowAnnotation(
  Class = marker_class,
  col = list(
    Class = marker_class_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = grid::unit(
    1.1,
    "mm"
  ),
  annotation_width = grid::unit(
    1.1,
    "mm"
  )
)

heatmap_object <- ComplexHeatmap::Heatmap(
  row_zscore_matrix,
  name = "Z-score",
  col = heatmap_color_function,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,

  row_split = marker_class,
  row_title = c("", ""),
  row_title_gp = grid::gpar(
    fontsize = 0
  ),
  row_gap = grid::unit(
    1.5,
    "mm"
  ),

  row_names_side = "left",
  row_names_gp = grid::gpar(
    fontfamily = BASE_FAMILY,
    fontsize = 5,
    col = "black"
  ),

  column_title = "Platelet score (low to high)",
  column_title_side = "bottom",
  column_title_gp = grid::gpar(
    fontfamily = BASE_FAMILY,
    fontsize = TEXT_SIZE,
    col = "black",
    fontface = "plain"
  ),

  rect_gp = grid::gpar(
    col = NA
  ),

  border = FALSE,

  width = grid::unit(
    28,
    "mm"
  ),

  height = grid::unit(
    65,
    "mm"
  ),

  right_annotation = row_annotation
)

################################################################################
# 13. Export PDF with locked manual legends
################################################################################

save_heatmap_pdf <- function(
    heatmap_object,
    output_path,
    width_cm,
    height_cm,
    class_colors,
    heatmap_limit,
    color_function
) {
  grDevices::pdf(
    file = output_path,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    family = BASE_FAMILY,
    useDingbats = FALSE
  )

  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  grid::grid.newpage()

  grid::pushViewport(
    grid::viewport(
      x = grid::unit(
        0.35,
        "npc"
      ),
      y = grid::unit(
        0.50,
        "npc"
      ),
      width = grid::unit(
        1,
        "npc"
      ),
      height = grid::unit(
        1,
        "npc"
      ),
      just = c(
        "center",
        "center"
      )
    )
  )

  ComplexHeatmap::draw(
    heatmap_object,
    newpage = FALSE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE,
    padding = grid::unit(
      c(
        2,
        2,
        2,
        4
      ),
      "mm"
    )
  )

  grid::upViewport(0)

  legend_left <- grid::unit(
    0.705,
    "npc"
  )

  square_width <- grid::unit(
    1.5,
    "mm"
  )

  square_height <- grid::unit(
    1.5,
    "mm"
  )

  text_gap <- grid::unit(
    2.2,
    "mm"
  )

  class_x <- legend_left
  class_text_x <- legend_left + text_gap

  grid::grid.rect(
    x = class_x,
    y = grid::unit(
      0.953,
      "npc"
    ),
    width = square_width,
    height = square_height,
    just = "left",
    gp = grid::gpar(
      fill = class_colors["Epithelial"],
      col = NA
    )
  )

  grid::grid.text(
    "Epithelial",
    x = class_text_x,
    y = grid::unit(
      0.953,
      "npc"
    ),
    just = "left",
    gp = grid::gpar(
      fontfamily = BASE_FAMILY,
      fontsize = 4.5,
      col = "black"
    )
  )

  grid::grid.rect(
    x = class_x,
    y = grid::unit(
      0.917,
      "npc"
    ),
    width = square_width,
    height = square_height,
    just = "left",
    gp = grid::gpar(
      fill = class_colors["Mesenchymal"],
      col = NA
    )
  )

  grid::grid.text(
    "Mesenchymal",
    x = class_text_x,
    y = grid::unit(
      0.917,
      "npc"
    ),
    just = "left",
    gp = grid::gpar(
      fontfamily = BASE_FAMILY,
      fontsize = 4.5,
      col = "black"
    )
  )

  z_title_x <- legend_left
  z_bar_x <- legend_left
  z_text_x <- z_bar_x + grid::unit(
    2.7,
    "mm"
  )

  z_title_y <- grid::unit(
    0.184,
    "npc"
  )

  z_bar_y <- grid::unit(
    0.114,
    "npc"
  )

  grid::grid.text(
    "Z-score",
    x = z_title_x,
    y = z_title_y,
    just = "left",
    gp = grid::gpar(
      fontfamily = BASE_FAMILY,
      fontsize = 4.5,
      col = "black"
    )
  )

  z_colors <- color_function(
    seq(
      heatmap_limit,
      -heatmap_limit,
      length.out = 100
    )
  )

  z_raster <- as.raster(
    matrix(
      z_colors,
      ncol = 1
    )
  )

  grid::grid.raster(
    z_raster,
    x = z_bar_x,
    y = z_bar_y,
    width = grid::unit(
      2.0,
      "mm"
    ),
    height = grid::unit(
      6.0,
      "mm"
    ),
    just = "left",
    interpolate = TRUE
  )

  grid::grid.text(
    sprintf(
      "%.1f",
      heatmap_limit
    ),
    x = z_text_x,
    y = z_bar_y + grid::unit(
      3.0,
      "mm"
    ),
    just = "left",
    gp = grid::gpar(
      fontfamily = BASE_FAMILY,
      fontsize = 4,
      col = "black"
    )
  )

  grid::grid.text(
    sprintf(
      "%.1f",
      -heatmap_limit
    ),
    x = z_text_x,
    y = z_bar_y - grid::unit(
      3.0,
      "mm"
    ),
    just = "left",
    gp = grid::gpar(
      fontfamily = BASE_FAMILY,
      fontsize = 4,
      col = "black"
    )
  )

  generated_figures <<- unique(
    c(
      generated_figures,
      output_path
    )
  )

  invisible(output_path)
}

save_heatmap_pdf(
  heatmap_object = heatmap_object,
  output_path = output_paths$figure,
  width_cm = 6,
  height_cm = 7.3,
  class_colors = marker_class_colors,
  heatmap_limit = heatmap_limit,
  color_function = heatmap_color_function
)

################################################################################
# 14. Final validation and QC log
################################################################################

required_outputs <- c(
  output_paths$figure,
  output_paths$marker_coverage,
  output_paths$direction_validation,
  output_paths$bin_summary,
  output_paths$bin_means,
  output_paths$row_zscores
)

missing_outputs <- required_outputs[
  !file.exists(required_outputs)
]

if (length(missing_outputs) > 0L) {
  fail(
    "Required outputs were not generated: ",
    paste(
      missing_outputs,
      collapse = " | "
    )
  )
}

final_status <- if (
  length(warnings_log) == 0L
) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

qc_lines <- c(
  "TCGA-PRAD curated EMT score-gradient heatmap",
  "============================================",
  "",
  paste(
    "Date/time:",
    as.character(Sys.time())
  ),
  paste(
    "Project directory:",
    project_dir
  ),
  paste(
    "Final status:",
    final_status
  ),
  "",
  "ANALYTICAL SCOPE",
  "----------------",
  "Descriptive five-bin visualization of the continuous Score_z gradient.",
  "Continuous regression model: not fitted",
  "Differential-expression testing across bins: not performed",
  "Q1-versus-Q4 comparison: not performed",
  "MSigDB Hallmark EMT gene set: not used",
  "Marker panel: curated epithelial/mesenchymal EMT markers",
  "",
  "CANONICAL INPUTS",
  "----------------",
  paste(
    "Score table:",
    input_paths$score
  ),
  paste(
    "Counts:",
    input_paths$counts
  ),
  paste(
    "Gene map:",
    input_paths$gene_map
  ),
  "",
  "COHORT AND ALIGNMENT",
  "--------------------",
  paste(
    "Canonical score rows:",
    nrow(score_df)
  ),
  paste(
    "Unique patients:",
    length(unique(score_df$patient_id))
  ),
  paste(
    "Counts dimensions before alignment:",
    paste(
      counts_dimensions_before_alignment,
      collapse = " x "
    )
  ),
  paste(
    "Aligned count samples:",
    ncol(counts_prefiltered)
  ),
  "Alignment method: exact full sample_id",
  "Barcode normalization: not used",
  "Barcode truncation: not used",
  "",
  "GENES AND TRANSFORMATION",
  "------------------------",
  paste(
    "Raw genes:",
    raw_gene_count
  ),
  paste(
    "Genes before prefilter:",
    genes_before_prefilter
  ),
  paste(
    "Genes after prefilter:",
    genes_after_prefilter
  ),
  "Prefilter: rowSums(counts >= 10) >= 10",
  "Expression transformation: DESeq2 blind VST",
  "Duplicate gene-symbol handling: maximum VST variance",
  "",
  "CURATED EMT MARKERS",
  "-------------------",
  paste(
    "Markers requested:",
    length(EMT_MARKER_GENES)
  ),
  paste(
    "Markers present:",
    length(genes_present)
  ),
  paste(
    "Epithelial markers present:",
    paste(
      epithelial_present,
      collapse = ", "
    )
  ),
  paste(
    "Mesenchymal markers present:",
    paste(
      mesenchymal_present,
      collapse = ", "
    )
  ),
  "",
  "FIVE-BIN GRADIENT",
  "-----------------",
  paste(
    "Number of bins:",
    NUMBER_OF_BINS
  ),
  "Bin labels: B1, B2, B3, B4, B5",
  "Direction: B1 = lowest Score_z; B5 = highest Score_z",
  "Tie handling: Score_z followed by exact sample_id",
  paste(
    capture.output(
      print(
        bin_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "DIRECTION VALIDATION",
  "--------------------",
  paste(
    capture.output(
      print(
        direction_validation,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "HEATMAP",
  "-------",
  "Rows: curated EMT markers ordered epithelial then mesenchymal",
  "Columns: mean VST expression in five ordered score bins",
  "Displayed values: row Z-scores across the five bin means",
  paste(
    "Color limit:",
    signif(
      heatmap_limit,
      6
    )
  ),
  paste(
    "Figure:",
    output_paths$figure
  ),
  "",
  "GENERATED TABLES",
  "----------------",
  if (
    length(generated_tables) == 0L
  ) {
    "None"
  } else {
    paste(
      " -",
      generated_tables
    )
  },
  "",
  "GENERATED FIGURES",
  "-----------------",
  if (
    length(generated_figures) == 0L
  ) {
    "None"
  } else {
    paste(
      " -",
      generated_figures
    )
  },
  "",
  "WARNINGS",
  "--------",
  if (
    length(warnings_log) == 0L
  ) {
    "None"
  } else {
    paste(
      " -",
      warnings_log
    )
  },
  "",
  "SESSION INFORMATION",
  "-------------------",
  paste(
    capture.output(
      sessionInfo()
    ),
    collapse = "\n"
  ),
  "",
  paste(
    "FINAL STATUS:",
    final_status
  )
)

writeLines(
  qc_lines,
  con = output_paths$qc_log,
  useBytes = TRUE
)

msg(
  "Figure written:",
  output_paths$figure
)

msg(
  "QC log written:",
  output_paths$qc_log
)

msg(
  "Final status:",
  final_status
)
