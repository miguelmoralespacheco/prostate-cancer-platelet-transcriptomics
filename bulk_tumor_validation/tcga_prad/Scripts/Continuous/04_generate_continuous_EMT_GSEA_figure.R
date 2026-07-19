#!/usr/bin/env Rscript
###############################################################################
# TCGA-PRAD | Continuous EMT GSEA figure
#
# Generates only the EMT running-enrichment figure for the complementary
# continuous platelet-associated transcriptional score analysis.
#
# No DESeq2 or fgsea model is run here. The figure is reconstructed from the
# canonical no-signature Wald-statistic rank, frozen Hallmark pathways, and
# canonical continuous GSEA results.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

message("[1/5] Resolving canonical inputs")

###############################################################################
# 1. Paths and locked constants
###############################################################################

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

script_path <- file.path(
  project_dir,
  "Scripts",
  "Continuous",
  "04_generate_continuous_EMT_GSEA_figure.R"
)

input_paths <- list(
  rank_rds = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Objects",
    "GSEA",
    "GSEA_continuous_DESeq2_Wald_stat_rank.rds"
  ),
  hallmark_pathways_rds = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Objects",
    "GSEA",
    "GSEA_continuous_Hallmark_pathways.rds"
  ),
  gsea_full = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Tables",
    "GSEA",
    "GSEA_continuous_no_signature_Hallmark_full.csv"
  ),
  gsea_metadata = file.path(
    GENERATED_RESULTS_DIR, "Continuous",
    "Tables",
    "GSEA",
    "GSEA_continuous_metadata.json"
  ),
  signature = file.path(
    RESOURCE_DIR,
    "platelet_associated_transcriptional_signature.tsv"
  )
)

figures_dir <- file.path(
  GENERATED_FIGURES_DIR, "Continuous",
  "GSEA"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Logs"
)

output_pdf <- file.path(
  figures_dir,
  "TCGA_continuous_EMT_GSEA_plot.pdf"
)

output_qc <- file.path(
  logs_dir,
  "TCGA_continuous_EMT_GSEA_figure_QC.txt"
)

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

EMT_PATHWAY <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
ES_TOLERANCE <- 1e-6

BASE_FAMILY <- "Helvetica"
GSEA_FONT_SIZE <- 5
GSEA_TITLE_SIZE <- 5
GSEA_AXIS_SIZE <- 4.35
GSEA_LINE_LW <- 0.16
GSEA_BORDER_LW <- 0.20

COL_ES <- "#95B300"
COL_NEG <- "#1C5DA3"
COL_POS <- "#DB0C07"

PDF_WIDTH_CM <- 4.2
PDF_HEIGHT_CM <- 3.9

warnings_vec <- character(0)

###############################################################################
# 2. Helpers
###############################################################################

fail_with_qc <- function(msg) {
  input_lines <- unlist(
    lapply(
      names(input_paths),
      function(nm) {
        paste0("  ", nm, ": ", input_paths[[nm]])
      }
    ),
    use.names = FALSE
  )

  writeLines(
    c(
      "TCGA-PRAD continuous EMT GSEA figure QC",
      "========================================",
      "",
      paste0("Date: ", as.character(Sys.time())),
      paste0("Project dir: ", project_dir),
      paste0("Script: ", script_path),
      "",
      "Canonical inputs:",
      input_lines,
      "",
      paste0("Output PDF: ", output_pdf),
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
  warnings_vec <<- unique(c(warnings_vec, as.character(msg)))
}

read_csv_base <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

format_fdr <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }

  if (x < 0.001) {
    return(format(x, scientific = TRUE, digits = 2))
  }

  sprintf("%.3f", x)
}

require_columns <- function(tab, required, label) {
  missing_columns <- setdiff(required, names(tab))

  if (length(missing_columns) > 0) {
    fail_with_qc(
      paste0(
        label,
        " is missing required columns: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
}

###############################################################################
# 3. Input validation
###############################################################################

for (path in unname(unlist(input_paths))) {
  if (!file.exists(path)) {
    fail_with_qc(paste0("Missing required input: ", path))
  }
}

required_packages <- c("ggplot2", "jsonlite", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  fail_with_qc(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

gsea_metadata <- jsonlite::read_json(
  input_paths$gsea_metadata,
  simplifyVector = TRUE
)

metadata_model <- as.character(gsea_metadata$ranking$model)
metadata_rank_column <- as.character(gsea_metadata$ranking$column)
metadata_rank_metric <- as.character(gsea_metadata$ranking$metric)
metadata_duplicate_rule <- as.character(gsea_metadata$ranking$duplicate_rule)
metadata_duplicate_rows <- suppressWarnings(
  as.integer(gsea_metadata$ranking$duplicate_gene_symbol_rows_removed)
)

if (!identical(metadata_model, "~ Score_z")) {
  fail_with_qc(
    paste0("Canonical GSEA metadata model is not ~ Score_z: ", metadata_model)
  )
}

if (!identical(metadata_rank_column, "stat")) {
  fail_with_qc(
    paste0("Canonical GSEA ranking column is not stat: ", metadata_rank_column)
  )
}

if (
  length(metadata_rank_metric) != 1L ||
  is.na(metadata_rank_metric) ||
  !grepl("Wald statistic", metadata_rank_metric, fixed = TRUE)
) {
  fail_with_qc("Canonical GSEA metadata does not identify a Wald-statistic rank.")
}

if (
  !isTRUE(gsea_metadata$input$no_signature_table) ||
  !identical(gsea_metadata$input$signature_genes_present, FALSE)
) {
  fail_with_qc(
    "Canonical GSEA metadata does not confirm a no-signature input table."
  )
}

if (
  length(metadata_duplicate_rule) != 1L ||
  is.na(metadata_duplicate_rule) ||
  !nzchar(metadata_duplicate_rule)
) {
  fail_with_qc("Canonical ranking duplicate rule is undocumented.")
}

ranks <- readRDS(input_paths$rank_rds)

if (!is.numeric(ranks)) {
  fail_with_qc("Canonical GSEA rank object is not numeric.")
}

if (is.null(names(ranks)) || length(names(ranks)) != length(ranks)) {
  fail_with_qc("Canonical GSEA rank object lacks one gene name per value.")
}

if (any(is.na(names(ranks)) | trimws(names(ranks)) == "")) {
  fail_with_qc("Canonical GSEA rank object contains empty gene names.")
}

if (any(!is.finite(ranks))) {
  fail_with_qc("Canonical GSEA rank object contains non-finite values.")
}

if (anyDuplicated(names(ranks))) {
  fail_with_qc(
    paste0(
      "Canonical GSEA rank object still contains duplicated gene names ",
      "despite its documented upstream rule: ",
      metadata_duplicate_rule
    )
  )
}

if (length(ranks) < 2L || any(diff(as.numeric(ranks)) > 0)) {
  fail_with_qc(
    "Canonical GSEA rank is not ordered from highest to lowest Wald statistic."
  )
}

metadata_ranked_genes <- suppressWarnings(
  as.integer(gsea_metadata$ranking$ranked_gene_symbols)
)

if (
  is.finite(metadata_ranked_genes) &&
  length(ranks) != metadata_ranked_genes
) {
  fail_with_qc(
    paste0(
      "Canonical rank length differs from metadata: rank = ",
      length(ranks),
      "; metadata = ",
      metadata_ranked_genes,
      "."
    )
  )
}

signature <- read_canonical_platelet_signature()
signature_candidates <- c(
  "gene_name",
  "gene_symbol",
  "symbol",
  "gene",
  "external_gene_name"
)
signature_column <- signature_candidates[
  signature_candidates %in% names(signature)
][1]

if (is.na(signature_column)) {
  fail_with_qc(
    paste0(
      "Official signature lacks a recognized gene column: ",
      paste(signature_candidates, collapse = ", ")
    )
  )
}

signature_genes <- unique(
  trimws(as.character(signature[[signature_column]]))
)
signature_genes <- signature_genes[
  !is.na(signature_genes) & nzchar(signature_genes)
]

if (length(signature_genes) != 41L) {
  fail_with_qc(
    paste0(
      "Expected 41 unique official signature genes. Observed: ",
      length(signature_genes)
    )
  )
}

signature_overlap <- intersect(
  toupper(signature_genes),
  toupper(names(ranks))
)

if (length(signature_overlap) > 0) {
  fail_with_qc(
    paste0(
      "Official signature genes remain in the canonical rank: ",
      paste(signature_overlap, collapse = ", ")
    )
  )
}

hallmark_pathways <- readRDS(input_paths$hallmark_pathways_rds)

if (!is.list(hallmark_pathways) || is.null(names(hallmark_pathways))) {
  fail_with_qc("Frozen Hallmark pathway object is not a named list.")
}

if (sum(names(hallmark_pathways) == EMT_PATHWAY) != 1L) {
  fail_with_qc(
    paste0(
      "Frozen Hallmark object must contain EMT exactly once. Observed: ",
      sum(names(hallmark_pathways) == EMT_PATHWAY)
    )
  )
}

emt_genes <- unique(
  trimws(as.character(hallmark_pathways[[EMT_PATHWAY]]))
)
emt_genes <- emt_genes[!is.na(emt_genes) & nzchar(emt_genes)]

gsea_full <- read_csv_base(input_paths$gsea_full)
require_columns(
  gsea_full,
  c("pathway", "ES", "NES", "padj", "direction"),
  "Canonical full GSEA table"
)

gsea_full$pathway <- as.character(gsea_full$pathway)
emt_result <- gsea_full[
  gsea_full$pathway == EMT_PATHWAY,
  ,
  drop = FALSE
]

if (nrow(emt_result) != 1L) {
  fail_with_qc(
    paste0(
      "Canonical full GSEA table must contain EMT exactly once. Observed: ",
      nrow(emt_result)
    )
  )
}

canonical_es <- suppressWarnings(as.numeric(emt_result$ES[1]))
canonical_nes <- suppressWarnings(as.numeric(emt_result$NES[1]))
canonical_fdr <- suppressWarnings(as.numeric(emt_result$padj[1]))
canonical_direction <- as.character(emt_result$direction[1])

if (!all(is.finite(c(canonical_es, canonical_nes, canonical_fdr)))) {
  fail_with_qc("Canonical EMT ES, NES, or FDR is non-finite.")
}

if (
  length(canonical_direction) != 1L ||
  is.na(canonical_direction) ||
  !nzchar(canonical_direction)
) {
  fail_with_qc("Canonical EMT enrichment direction is missing.")
}

###############################################################################
# 4. Manual running enrichment score reconstruction
###############################################################################

message("[2/5] Reconstructing the canonical EMT running score")

emt_genes_in_rank <- intersect(emt_genes, names(ranks))
number_ranked_genes <- length(ranks)
number_emt_genes <- length(emt_genes_in_rank)

if (number_emt_genes < 10L) {
  fail_with_qc(
    paste0(
      "Fewer than 10 EMT genes are represented in the canonical rank: ",
      number_emt_genes
    )
  )
}

if (number_ranked_genes <= number_emt_genes) {
  fail_with_qc(
    paste0(
      "Invalid ranking dimensions: ranked genes = ",
      number_ranked_genes,
      "; EMT genes = ",
      number_emt_genes,
      "."
    )
  )
}

if ("size" %in% names(emt_result)) {
  canonical_size <- suppressWarnings(as.integer(emt_result$size[1]))

  if (is.finite(canonical_size) && canonical_size != number_emt_genes) {
    fail_with_qc(
      paste0(
        "Canonical EMT size differs from pathway overlap: table = ",
        canonical_size,
        "; reconstructed = ",
        number_emt_genes,
        "."
      )
    )
  }
}

hit_positions <- which(names(ranks) %in% emt_genes_in_rank)
hit_indicator <- logical(number_ranked_genes)
hit_indicator[hit_positions] <- TRUE

hit_weights <- abs(ranks[hit_positions])
hit_weight_sum <- sum(hit_weights)

if (!is.finite(hit_weight_sum) || hit_weight_sum == 0) {
  fail_with_qc("EMT hit weights are zero or non-finite.")
}

hit_weights <- hit_weights / hit_weight_sum
miss_penalty <- 1 / (number_ranked_genes - number_emt_genes)

running_score <- numeric(number_ranked_genes)
current_score <- 0
hit_index <- 1L

for (i in seq_len(number_ranked_genes)) {
  if (hit_indicator[i]) {
    current_score <- current_score + hit_weights[hit_index]
    hit_index <- hit_index + 1L
  } else {
    current_score <- current_score - miss_penalty
  }

  running_score[i] <- current_score
}

positive_es <- max(running_score, na.rm = TRUE)
negative_es <- min(running_score, na.rm = TRUE)

manual_es <- if (abs(positive_es) >= abs(negative_es)) {
  positive_es
} else {
  negative_es
}

manual_es_difference <- abs(manual_es - canonical_es)

if (!is.finite(manual_es_difference)) {
  fail_with_qc("Manual and canonical EMT ES values could not be compared.")
}

if (manual_es_difference > ES_TOLERANCE) {
  fail_with_qc(
    paste0(
      "Manual EMT ES does not match the canonical fgsea ES within ",
      format(ES_TOLERANCE, scientific = TRUE),
      ". Manual ES = ",
      format(manual_es, digits = 16),
      "; canonical ES = ",
      format(canonical_es, digits = 16),
      "; absolute difference = ",
      format(manual_es_difference, scientific = TRUE, digits = 16),
      "."
    )
  )
}

###############################################################################
# 5. Locked three-panel figure
###############################################################################

message("[3/5] Building the locked three-panel figure")

metrics_block <- paste(
  c(
    paste0("NES = ", round(canonical_nes, 2)),
    paste0("FDR = ", format_fdr(canonical_fdr))
  ),
  collapse = "\n"
)

df_es <- data.frame(
  rank = seq_len(number_ranked_genes),
  ES = running_score
)

df_hits <- data.frame(rank = hit_positions)

df_metric <- data.frame(
  rank = seq_len(number_ranked_genes),
  metric = as.numeric(ranks)
)

lower_quantile <- as.numeric(
  stats::quantile(df_metric$metric, 0.01, na.rm = TRUE)
)
upper_quantile <- as.numeric(
  stats::quantile(df_metric$metric, 0.99, na.rm = TRUE)
)

df_metric$metric_clip <- pmax(
  lower_quantile,
  pmin(upper_quantile, df_metric$metric)
)

df_strip <- data.frame(
  rank = df_metric$rank,
  y = 0,
  z = df_metric$metric_clip
)

theme_nature <- ggplot2::theme_classic(
  base_size = GSEA_FONT_SIZE,
  base_family = BASE_FAMILY
) +
  ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      linewidth = GSEA_LINE_LW,
      colour = "grey90"
    ),
    axis.ticks = ggplot2::element_line(linewidth = GSEA_LINE_LW),
    axis.ticks.length = grid::unit(0.5, "mm"),
    panel.border = ggplot2::element_rect(
      linewidth = GSEA_BORDER_LW,
      colour = "black",
      fill = NA
    ),
    plot.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.title = ggplot2::element_text(
      face = "plain",
      size = GSEA_TITLE_SIZE,
      hjust = 0.5,
      margin = ggplot2::margin(b = 0)
    ),
    axis.title = ggplot2::element_text(size = GSEA_AXIS_SIZE),
    axis.text = ggplot2::element_text(size = GSEA_AXIS_SIZE),
    plot.margin = grid::unit(c(1, 1, 1, 1), "mm")
  )

p_es <- ggplot2::ggplot(
  df_es,
  ggplot2::aes(x = rank, y = ES)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = GSEA_LINE_LW,
    color = "black",
    alpha = 0.6
  ) +
  ggplot2::geom_line(
    linewidth = 0.75,
    color = COL_ES,
    lineend = "round"
  ) +
  ggplot2::labs(
    title = "Hallmark epithelial-mesenchymal transition",
    x = NULL,
    y = "Running enrichment\nscore"
  ) +
  ggplot2::annotate(
    "text",
    x = Inf,
    y = Inf,
    label = metrics_block,
    hjust = 1.1,
    vjust = 1.26,
    lineheight = 0.84,
    family = BASE_FAMILY,
    size = GSEA_FONT_SIZE / ggplot2::.pt
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_nature +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      face = "plain",
      size = GSEA_TITLE_SIZE,
      hjust = 0.5,
      margin = ggplot2::margin(b = 0.40, unit = "mm")
    ),
    plot.margin = grid::unit(c(0.8, 1, 0.5, 1), "mm")
  )

p_hits <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = df_hits,
    ggplot2::aes(x = rank, xend = rank, y = 1, yend = 0),
    linewidth = 0.10,
    color = "black",
    alpha = 0.60
  ) +
  ggplot2::geom_tile(
    data = df_strip,
    ggplot2::aes(x = rank, y = -0.33, fill = z),
    height = 0.55,
    width = 1
  ) +
  ggplot2::scale_fill_gradient2(
    low = COL_NEG,
    mid = "white",
    high = COL_POS,
    midpoint = 0,
    guide = "none"
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = -0.92,
    label = "Platelet score high",
    hjust = 0,
    vjust = 0.5,
    family = BASE_FAMILY,
    fontface = "plain",
    size = 4 / ggplot2::.pt
  ) +
  ggplot2::annotate(
    "text",
    x = number_ranked_genes,
    y = -0.92,
    label = "Platelet score low",
    hjust = 1,
    vjust = 0.5,
    family = BASE_FAMILY,
    fontface = "plain",
    size = 4 / ggplot2::.pt
  ) +
  ggplot2::coord_cartesian(
    xlim = c(1, number_ranked_genes),
    ylim = c(-0.95, 1.02),
    clip = "off"
  ) +
  ggplot2::theme_void(
    base_family = BASE_FAMILY,
    base_size = GSEA_FONT_SIZE
  ) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.margin = grid::unit(c(0.2, 1, 0.2, 1), "mm")
  )

p_metric <- ggplot2::ggplot(
  df_metric,
  ggplot2::aes(x = rank, y = metric)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = GSEA_LINE_LW,
    color = "black",
    alpha = 0.75
  ) +
  ggplot2::geom_area(fill = "grey85", alpha = 1) +
  ggplot2::geom_line(linewidth = 0.20, color = "grey55") +
  ggplot2::labs(
    x = "Rank in ordered dataset",
    y = "Ranking metric\n(Wald statistic)"
  ) +
  theme_nature +
  ggplot2::theme(
    plot.margin = grid::unit(c(0.5, 1, 1, 1), "mm")
  )

emt_gsea_plot <- patchwork::wrap_plots(
  p_es,
  p_hits,
  p_metric,
  ncol = 1,
  heights = c(2.0, 0.55, 0.85)
) &
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.box.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.key = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.margin = ggplot2::margin(0.6, 0.6, 0.6, 0.6, "pt")
  )

message("[4/5] Writing the requested PDF")

withCallingHandlers(
  ggplot2::ggsave(
    filename = output_pdf,
    plot = emt_gsea_plot,
    width = PDF_WIDTH_CM,
    height = PDF_HEIGHT_CM,
    units = "cm",
    device = grDevices::pdf,
    family = BASE_FAMILY,
    useDingbats = FALSE,
    bg = "transparent"
  ),
  warning = function(w) {
    add_warning(conditionMessage(w))
  }
)

pdf_info <- file.info(output_pdf)

if (
  !file.exists(output_pdf) ||
  nrow(pdf_info) != 1L ||
  is.na(pdf_info$size) ||
  pdf_info$size <= 0
) {
  fail_with_qc(paste0("Output PDF is missing or empty: ", output_pdf))
}

###############################################################################
# 6. QC report
###############################################################################

message("[5/5] Writing QC report")

final_status <- if (length(warnings_vec) > 0) "WARNING" else "PASS"

warning_lines <- if (length(warnings_vec) == 0) {
  "  none"
} else {
  paste0("  - ", warnings_vec)
}

qc_lines <- c(
  "TCGA-PRAD continuous EMT GSEA figure QC",
  "========================================",
  "",
  paste0("Date: ", as.character(Sys.time())),
  paste0("Project dir: ", project_dir),
  paste0("Script: ", script_path),
  "",
  "Canonical inputs used:",
  paste0("  Rank RDS: ", input_paths$rank_rds),
  paste0("  Frozen Hallmark pathways RDS: ", input_paths$hallmark_pathways_rds),
  paste0("  Full continuous GSEA table: ", input_paths$gsea_full),
  paste0("  Continuous GSEA metadata: ", input_paths$gsea_metadata),
  paste0("  Official 41-gene signature: ", input_paths$signature),
  "",
  "Analytical identity:",
  paste0("  Continuous model: ", metadata_model),
  "  Cohort: complete TCGA-PRAD tumor cohort used by the canonical continuous model",
  paste0("  Ranking variable: ", metadata_rank_metric, " (", metadata_rank_column, ")"),
  "  Ranking order: descending Wald statistic",
  "  DESeq2 rerun: no",
  "  fgsea rerun: no",
  "  Platelet score recalculation: no",
  "",
  "Ranking and signature checks:",
  paste0("  Total ranked genes: ", number_ranked_genes),
  paste0("  Empty gene names: ", sum(trimws(names(ranks)) == "")),
  paste0("  Duplicated ranked gene names: ", anyDuplicated(names(ranks))),
  paste0("  Canonical duplicate rows removed upstream: ", metadata_duplicate_rows),
  paste0("  Canonical duplicate rule: ", metadata_duplicate_rule),
  paste0("  Official signature genes: ", length(signature_genes)),
  paste0("  Official signature genes remaining in rank: ", length(signature_overlap)),
  "  Signature exclusion confirmed: TRUE",
  "",
  "EMT reconstruction:",
  paste0("  Pathway: ", EMT_PATHWAY),
  paste0("  EMT genes in frozen pathway: ", length(emt_genes)),
  paste0("  EMT genes represented in rank: ", number_emt_genes),
  paste0("  Manual ES: ", format(manual_es, digits = 16)),
  paste0("  Canonical ES: ", format(canonical_es, digits = 16)),
  paste0(
    "  Absolute ES difference: ",
    format(manual_es_difference, scientific = TRUE, digits = 16)
  ),
  paste0("  Required ES tolerance: ", format(ES_TOLERANCE, scientific = TRUE)),
  paste0("  NES: ", format(canonical_nes, digits = 16)),
  paste0("  FDR: ", format(canonical_fdr, scientific = TRUE, digits = 16)),
  paste0("  Enrichment direction: ", canonical_direction),
  "",
  "Locked figure specification:",
  paste0("  PDF dimensions: ", PDF_WIDTH_CM, " x ", PDF_HEIGHT_CM, " cm"),
  "  Background: transparent",
  paste0("  Font family: ", BASE_FAMILY),
  paste0("  Title size: ", GSEA_TITLE_SIZE, " pt"),
  paste0("  Axis number size: ", GSEA_AXIS_SIZE, " pt"),
  paste0("  Axis title size: ", GSEA_AXIS_SIZE, " pt"),
  paste0("  Panel border linewidth: ", GSEA_BORDER_LW),
  paste0("  Structural linewidth: ", GSEA_LINE_LW),
  paste0("  EMT curve: ", COL_ES, "; linewidth 0.75"),
  "  EMT hits: linewidth 0.10; alpha 0.60",
  paste0("  Rank strip: ", COL_NEG, " / white / ", COL_POS),
  "  Ranking area: grey85; ranking line: grey55; linewidth 0.20",
  "  Patchwork heights: 2.0, 0.55, 0.85",
  "  External plot margin: 0.6 pt per side",
  "",
  "Output:",
  paste0("  PDF: ", output_pdf),
  paste0("  PDF size bytes: ", pdf_info$size),
  paste0("  QC: ", output_qc),
  "",
  "Warnings:",
  warning_lines,
  "",
  paste0("Final status: ", final_status)
)

writeLines(qc_lines, con = output_qc, useBytes = TRUE)

message("Output PDF: ", output_pdf)
message("QC report: ", output_qc)
message("Final status: ", final_status)

if (!identical(final_status, "PASS")) {
  quit(status = 2L, save = "no")
}
