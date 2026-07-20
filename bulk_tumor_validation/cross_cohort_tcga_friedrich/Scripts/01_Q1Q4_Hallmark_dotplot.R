#!/usr/bin/env Rscript
###############################################################################
# Combined TCGA-PRAD + Friedrich/GSE134051 Hallmark dot plot
#
# Cohort-wise pathway-level visualization for selected Hallmark pathways.
#
# Rows: selected Hallmark pathways
# Columns: TCGA-PRAD and Friedrich/GSE134051
# Dot color: NES
# Dot size: -log10(FDR)
#
# Positive NES = enriched in platelet score-high tumors
# Negative NES = enriched in platelet score-low tumors
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

###############################################################################
# 1) Portable paths
###############################################################################

get_script_path <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_arg) != 1L) {
    stop("Could not determine script path from Rscript.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

derive_repo_root <- function(script_path) {
  override <- Sys.getenv("PLATELET_REPO_ROOT", unset = "")
  root <- if (nzchar(override)) {
    normalizePath(override, mustWork = TRUE)
  } else {
    normalizePath(
      dirname(dirname(dirname(dirname(script_path)))),
      mustWork = TRUE
    )
  }

  expected <- file.path(
    root,
    "bulk_tumor_validation",
    "cross_cohort_tcga_friedrich",
    "Scripts",
    basename(script_path)
  )
  if (
    !file.exists(expected) ||
    !identical(normalizePath(expected, mustWork = TRUE), script_path)
  ) {
    stop("Repository root is incompatible with the script path.", call. = FALSE)
  }
  root
}

script_path <- get_script_path()
repo_root <- derive_repo_root(script_path)

COHORT_TCGA_ID <- "TCGA_PRAD"
COHORT_TCGA_LABEL <- "TCGA-PRAD"
COHORT_FRIEDRICH_ID <- "Friedrich_GSE134051"
COHORT_FRIEDRICH_LABEL <- "Friedrich/GSE134051"
COHORT_FRIEDRICH_DISPLAY_LABEL <- "Friedrich (GSE134051)"

input_tcga_gsea <- file.path(
  repo_root,
  "bulk_tumor_validation", "tcga_prad", "Contracts", "hallmark_gsea.tsv"
)
input_friedrich_gsea <- file.path(
  repo_root,
  "bulk_tumor_validation", "friedrich_gse134051", "Contracts",
  "hallmark_gsea.tsv"
)

results_dir <- file.path(
  repo_root,
  "bulk_tumor_validation", "cross_cohort_tcga_friedrich",
  "Results", "generated", "Q1Q4_Hallmark"
)
figures_dir <- file.path(
  repo_root,
  "bulk_tumor_validation", "cross_cohort_tcga_friedrich",
  "Figures", "generated", "Q1Q4_Hallmark"
)

output_pdf <- file.path(
  figures_dir,
  "Combined_TCGA_Friedrich_Hallmark_SELECTED_EXTENDED_dotplot_MAIN.pdf"
)
output_table <- file.path(
  results_dir,
  "Combined_TCGA_Friedrich_Hallmark_SELECTED_EXTENDED_dotplot_table.csv"
)
output_qc <- file.path(
  results_dir,
  "Combined_TCGA_Friedrich_Hallmark_SELECTED_EXTENDED_dotplot_QC.txt"
)

###############################################################################
# 2) Validation and helpers
###############################################################################

fail_clear <- function(...) {
  stop(paste0(...), call. = FALSE)
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    fail_clear("Missing required package(s): ", paste(missing, collapse = ", "))
  }
}

read_contract <- function(path, cohort_id) {
  if (!file.exists(path)) {
    fail_clear("Missing Hallmark contract for ", cohort_id, ".")
  }

  tab <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expected_names <- c("cohort_id", "analysis_id", "pathway", "nes", "padj")
  if (!identical(names(tab), expected_names)) {
    fail_clear(cohort_id, " Hallmark contract has an unexpected schema.")
  }
  if (
    nrow(tab) != 50L ||
    anyNA(tab$cohort_id) ||
    !identical(unique(as.character(tab$cohort_id)), cohort_id) ||
    anyNA(tab$analysis_id) ||
    !identical(unique(as.character(tab$analysis_id)), "Q1Q4")
  ) {
    fail_clear(cohort_id, " Hallmark contract has an unexpected cohort, analysis or row count.")
  }

  tab$pathway <- trimws(as.character(tab$pathway))
  tab$nes <- suppressWarnings(as.numeric(tab$nes))
  tab$padj <- suppressWarnings(as.numeric(tab$padj))

  if (
    anyNA(tab$pathway) ||
    any(tab$pathway == "") ||
    any(!grepl("^HALLMARK_[A-Z0-9_]+$", tab$pathway)) ||
    anyDuplicated(tab$pathway) ||
    anyDuplicated(tab[c("cohort_id", "analysis_id", "pathway")]) ||
    anyNA(tab$nes) ||
    any(!is.finite(tab$nes)) ||
    anyNA(tab$padj) ||
    any(!is.finite(tab$padj)) ||
    any(tab$padj < 0 | tab$padj > 1)
  ) {
    fail_clear(cohort_id, " Hallmark contract contains malformed or invalid values.")
  }

  tab
}

clean_hallmark_label <- function(pathway) {
  x <- gsub("^HALLMARK_", "", pathway)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("Tnfa", "TNF-alpha", x)
  x <- gsub("Nfkb", "NFkB", x)
  x <- gsub("Il6 Jak Stat3", "IL-6/JAK/STAT3", x)
  x <- gsub("Il2 Stat5", "IL2/STAT5", x)
  x <- gsub("Kras", "KRAS", x)
  x <- gsub("Tgf Beta", "TGF-beta", x)
  x <- gsub("P53", "p53", x)
  x <- gsub("Myc", "MYC", x)
  x <- gsub("Signaling", "signaling", x)
  x <- gsub("Via", "via", x)
  x
}

display_hallmark_label <- function(x) {
  x <- gsub(
    "Epithelial Mesenchymal Transition",
    "Epithelial–Mesenchymal Transition",
    x
  )
  x <- gsub("Inflammatory Response", "Inflammatory Response", x)
  x <- gsub("KRAS Signaling Up", "KRAS signaling Up", x)
  x <- gsub("KRAS Signaling Dn", "KRAS signaling Dn", x)
  x <- gsub("TGF-beta Signaling", "TGF-beta signaling", x)
  x <- gsub("Androgen Response", "Androgen Response", x)
  x <- gsub("MYC Targets V1", "MYC Targets V1", x)
  x
}

prep_gsea <- function(df, cohort_label) {
  df$pathway <- as.character(df$pathway)
  df$NES <- suppressWarnings(as.numeric(df$nes))
  df$padj <- suppressWarnings(as.numeric(df$padj))
  df <- df[
    !is.na(df$pathway) &
      df$pathway != "" &
      is.finite(df$NES) &
      !is.na(df$padj),
    ,
    drop = FALSE
  ]
  df$cohort <- cohort_label
  df$pathway_clean <- clean_hallmark_label(df$pathway)
  df$direction <- ifelse(df$NES > 0, "Score-high", "Score-low")
  df$neg_log10_FDR_raw <- -log10(pmax(df$padj, 1e-50))
  df$neg_log10_FDR <- pmin(df$neg_log10_FDR_raw, 50)
  df
}

###############################################################################
# 3) Style and selected pathways
###############################################################################

BASE_FAMILY <- "Helvetica"

COL_LOW_NES  <- "#1C5DA3"
COL_MID_NES  <- "white"
COL_HIGH_NES <- "#DB0C07"

TEXT_SIZE <- 5
AXIS_TEXT_SIZE <- 5
LEGEND_TEXT_SIZE <- 4
LINE_SIZE <- 0.20

COLOR_LIMIT <- 4

keep_main_ids <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_COAGULATION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_KRAS_SIGNALING_DN"
)

###############################################################################
# 4) Load, validate and prepare contracts
###############################################################################

require_packages(c("ggplot2", "scales", "grid"))

tcga_raw <- read_contract(input_tcga_gsea, COHORT_TCGA_ID)
friedrich_raw <- read_contract(input_friedrich_gsea, COHORT_FRIEDRICH_ID)

if (!setequal(tcga_raw$pathway, friedrich_raw$pathway)) {
  fail_clear("TCGA and Friedrich contracts do not contain identical Hallmark sets.")
}
if (
  !all(keep_main_ids %in% tcga_raw$pathway) ||
  !all(keep_main_ids %in% friedrich_raw$pathway)
) {
  fail_clear("One or more historically selected pathways are missing.")
}

tcga <- prep_gsea(tcga_raw, COHORT_TCGA_LABEL)
friedrich <- prep_gsea(friedrich_raw, COHORT_FRIEDRICH_LABEL)

tcga_sel <- tcga[
  tcga$pathway %in% keep_main_ids & tcga$padj < 0.05,
  ,
  drop = FALSE
]
friedrich_sel <- friedrich[
  friedrich$pathway %in% keep_main_ids & friedrich$padj < 0.05,
  ,
  drop = FALSE
]

n_tcga_selected_significant <- nrow(tcga_sel)
n_friedrich_selected_significant <- nrow(friedrich_sel)
common_selected <- intersect(tcga_sel$pathway, friedrich_sel$pathway)

if (length(common_selected) == 0L) {
  fail_clear("No selected Hallmark pathways were significant in both cohorts.")
}
if (
  n_tcga_selected_significant != 14L ||
  n_friedrich_selected_significant != 15L ||
  length(common_selected) != 14L
) {
  fail_clear("Selected-pathway counts differ from the validated current state.")
}

tcga_sel <- tcga_sel[tcga_sel$pathway %in% common_selected, , drop = FALSE]
friedrich_sel <- friedrich_sel[
  friedrich_sel$pathway %in% common_selected,
  ,
  drop = FALSE
]

plot_cols <- c(
  "pathway",
  "pathway_clean",
  "cohort",
  "NES",
  "padj",
  "direction",
  "neg_log10_FDR_raw",
  "neg_log10_FDR"
)
plot_df <- rbind(
  tcga_sel[, plot_cols, drop = FALSE],
  friedrich_sel[, plot_cols, drop = FALSE]
)

summary_df <- merge(
  tcga_sel[, c("pathway", "pathway_clean", "NES", "padj")],
  friedrich_sel[, c("pathway", "NES", "padj")],
  by = "pathway",
  suffixes = c("_TCGA", "_Friedrich")
)

summary_df$pathway_clean <- clean_hallmark_label(summary_df$pathway)
summary_df$same_direction <- sign(summary_df$NES_TCGA) ==
  sign(summary_df$NES_Friedrich)
summary_df$mean_NES <- rowMeans(
  cbind(summary_df$NES_TCGA, summary_df$NES_Friedrich),
  na.rm = TRUE
)
summary_df$mean_abs_NES <- rowMeans(
  abs(cbind(summary_df$NES_TCGA, summary_df$NES_Friedrich)),
  na.rm = TRUE
)
summary_df$order_block <- ifelse(
  summary_df$NES_TCGA > 0 & summary_df$NES_Friedrich > 0,
  "Concordant_positive",
  "Discordant_or_negative"
)
summary_df <- summary_df[
  order(
    summary_df$order_block != "Concordant_positive",
    -summary_df$mean_NES,
    -summary_df$mean_abs_NES
  ),
  ,
  drop = FALSE
]

pathway_order <- summary_df$pathway
pathway_label_map <- setNames(
  display_hallmark_label(summary_df$pathway_clean),
  summary_df$pathway
)
plot_df$pathway_factor <- factor(
  plot_df$pathway,
  levels = rev(pathway_order),
  labels = rev(pathway_label_map[pathway_order])
)
plot_df$cohort <- factor(
  plot_df$cohort,
  levels = c(COHORT_TCGA_LABEL, COHORT_FRIEDRICH_LABEL),
  labels = c(COHORT_TCGA_LABEL, COHORT_FRIEDRICH_LABEL)
)

if (!setequal(as.character(unique(plot_df$cohort)), c(
  COHORT_TCGA_LABEL,
  COHORT_FRIEDRICH_LABEL
))) {
  fail_clear("Cross-cohort labels are incomplete or unexpected.")
}

###############################################################################
# 5) Historical plot
###############################################################################

plot_df$x_pos <- ifelse(
  plot_df$cohort == COHORT_TCGA_LABEL,
  1.00,
  1.22
)

p_dot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = x_pos,
    y = pathway_factor,
    size = neg_log10_FDR
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(color = NES),
    shape = 16,
    alpha = 1
  ) +
  ggplot2::scale_color_gradient2(
    low = COL_LOW_NES,
    mid = COL_MID_NES,
    high = COL_HIGH_NES,
    midpoint = 0,
    limits = c(-COLOR_LIMIT, COLOR_LIMIT),
    breaks = c(-COLOR_LIMIT, 0, COLOR_LIMIT),
    labels = c(
      paste0("-", COLOR_LIMIT),
      "0",
      paste0(COLOR_LIMIT)
    ),
    oob = scales::squish,
    name = "Normalized\nenrichment\nscore (NES)"
  ) +
  ggplot2::scale_size_continuous(
    range = c(1.18, 3.68),
    limits = c(0, 50),
    breaks = c(10, 25, 50),
    name = expression(-log[10]~"FDR")
  )  +
  ggplot2::scale_x_continuous(
    breaks = c(1.00, 1.22),
    labels = c(COHORT_TCGA_LABEL, COHORT_FRIEDRICH_DISPLAY_LABEL),
    limits = c(0.90, 1.32),
    expand = c(0, 0),
    position = "bottom"
  ) +
  ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(add = 0.80)
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_family = BASE_FAMILY,
    base_size = TEXT_SIZE
  ) +
  ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = AXIS_TEXT_SIZE,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1,
      margin = ggplot2::margin(t = 1.5)
    ),
    axis.text.y = ggplot2::element_text(
      size = AXIS_TEXT_SIZE,
      color = "black",
      lineheight = 0.88,
      margin = ggplot2::margin(r = 1.2)
    ),
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = LINE_SIZE
    ),
    plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(
      size = LEGEND_TEXT_SIZE,
      hjust = 0.5
    ),
    legend.text = ggplot2::element_text(
      size = LEGEND_TEXT_SIZE
    ),
    legend.ticks = ggplot2::element_blank(),
    legend.ticks.length = grid::unit(0, "mm"),
    legend.key.height = grid::unit(1.8, "mm"),
    legend.key.width = grid::unit(1.8, "mm"),
    legend.spacing.y = grid::unit(5.5, "mm"),
    legend.box.spacing = grid::unit(2.5, "mm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(1.0, 0.8, 1.0, 1.0, "mm")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(9, "mm"),
      barwidth = grid::unit(2.0, "mm"),
      frame.colour = NA,
      ticks = FALSE,
      order = 1,
      theme = ggplot2::theme(
        legend.title = ggplot2::element_text(
          size = LEGEND_TEXT_SIZE,
          hjust = 0.5,
          margin = ggplot2::margin(b = 1.7, unit = "mm")
        )
      )
    ),
    size = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      order = 2,
      override.aes = list(
        shape = 21,
        fill = "grey85",
        color = "black",
        stroke = 0.10,
        alpha = 1
      )
    )
  )

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

ggplot2::ggsave(
  filename = output_pdf,
  plot = p_dot,
  width = 5.3,
  height = 6.2,
  units = "cm",
  device = grDevices::pdf,
  family = BASE_FAMILY,
  useDingbats = FALSE,
  bg = "transparent"
)

###############################################################################
# 6) Output table and QC
###############################################################################

summary_out <- summary_df[, c(
  "pathway",
  "pathway_clean",
  "NES_TCGA",
  "padj_TCGA",
  "NES_Friedrich",
  "padj_Friedrich",
  "same_direction",
  "mean_NES",
  "mean_abs_NES",
  "order_block"
)]
utils::write.csv(summary_out, output_table, row.names = FALSE)

qc_lines <- c(
  "Combined TCGA-PRAD + Friedrich/GSE134051 Hallmark dot plot QC",
  "==============================================================",
  "",
  "Inputs:",
  "  tcga_prad/Contracts/hallmark_gsea.tsv",
  "  friedrich_gse134051/Contracts/hallmark_gsea.tsv",
  "",
  paste0("TCGA pathways: ", nrow(tcga)),
  paste0("Friedrich pathways: ", nrow(friedrich)),
  paste0("Selected Hallmarks requested: ", length(keep_main_ids)),
  paste0("Selected significant in TCGA: ", n_tcga_selected_significant),
  paste0("Selected significant in Friedrich: ", n_friedrich_selected_significant),
  paste0("Selected significant in both: ", length(common_selected)),
  "Visible labels: TCGA-PRAD; Friedrich (GSE134051)",
  "",
  paste0("Table: ", basename(output_table)),
  paste0("PDF: ", basename(output_pdf)),
  "",
  "Final status: PASS"
)
writeLines(qc_lines, output_qc, useBytes = TRUE)

message("Validated Hallmarks: TCGA=50; Friedrich=50; shared displayed=14")
message("Final status: PASS")
