#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Q1/Q4 figures
#
# This script generates Q1/Q4 validation figures for Friedrich/GSE134051 using closed
# molecular results from the official 41-gene platelet-associated transcriptional
# score workflow.
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

input_score <- file.path(
  GENERATED_RESULTS_DIR, "Score", "Friedrich_GSE134051_primary_Q1Q4_groups_41genes.csv"
)

input_expression <- file.path(
  INPUT_DIR, "LocalLarge", "Friedrich_GSE134051_expression_logq.rds"
)

input_signature <- file.path(
  RESOURCE_DIR, "platelet_associated_transcriptional_signature.tsv"
)

input_limma <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4",
  "Tables",
  "LIMMA",
  "LIMMA_Q1Q4_41genes_no_signature_genes.csv"
)

input_gsea_full <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "GSEA",
  "GSEA_Q1Q4_no_signature_Hallmark_full.csv"
)

input_gsea_emt <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "GSEA",
  "GSEA_Q1Q4_no_signature_Hallmark_EMT_row.csv"
)

input_gsea_metadata <- file.path(
  GENERATED_RESULTS_DIR, "Q1Q4", "Tables", "GSEA", "GSEA_Q1Q4_metadata.json"
)

figures_dir <- file.path(GENERATED_FIGURES_DIR, "Q1Q4")
pca_dir <- file.path(figures_dir, "PCA")
volcano_dir <- file.path(figures_dir, "Volcano")
gsea_dir <- file.path(figures_dir, "GSEA")
logs_dir <- file.path(GENERATED_RESULTS_DIR, "Q1Q4", "Logs")

output_pca <- file.path(
  pca_dir,
  "PCA_Minimal.pdf"
)

output_volcano <- file.path(
  volcano_dir,
  "Friedrich_GSE134051_Q1Q4_limma_volcano.pdf"
)

output_emt <- file.path(
  gsea_dir,
  "Friedrich_GSE134051_Q1Q4_EMT_GSEA_plot.pdf"
)

output_nes <- file.path(
  gsea_dir,
  "Friedrich_GSE134051_Q1Q4_Hallmark_SELECTED_EXTENDED_NES_barplot_MAIN.pdf"
)

output_qc <- file.path(
  logs_dir,
  "Friedrich_GSE134051_Q1Q4_figures_QC.txt"
)

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 2. Helpers
###############################################################################

warnings_vec <- character(0)
figures_written <- character(0)
qc_failure_written <- FALSE

fail_clear <- function(msg) {
  if (!qc_failure_written) {
    qc_failure_written <<- TRUE

    sink(output_qc)
    on.exit(sink(), add = TRUE)

    cat("Friedrich/GSE134051 Q1/Q4 figures QC\n")
    cat("==========================\n\n")

    cat("Input paths:\n")
    cat("  score: ", input_score, "\n", sep = "")
    cat("  expression: ", input_expression, "\n", sep = "")
    cat("  signature: ", input_signature, "\n", sep = "")
    cat("  limma_no_signature: ", input_limma, "\n", sep = "")
    cat("  gsea_full: ", input_gsea_full, "\n", sep = "")
    cat("  gsea_emt: ", input_gsea_emt, "\n", sep = "")
    cat("  gsea_metadata: ", input_gsea_metadata, "\n", sep = "")

    cat("\nInputs found:\n")
    cat("  score: ", file.exists(input_score), "\n", sep = "")
    cat("  expression: ", file.exists(input_expression), "\n", sep = "")
    cat("  signature: ", file.exists(input_signature), "\n", sep = "")
    cat("  limma_no_signature: ", file.exists(input_limma), "\n", sep = "")
    cat("  gsea_full: ", file.exists(input_gsea_full), "\n", sep = "")
    cat("  gsea_emt: ", file.exists(input_gsea_emt), "\n", sep = "")
    cat("  gsea_metadata: ", file.exists(input_gsea_metadata), "\n", sep = "")

    cat("\nFigures written:\n")
    if (length(figures_written) == 0) {
      cat("  None\n")
    } else {
      for (fp in figures_written) {
        cat("  - ", fp, "\n", sep = "")
      }
    }

    cat("\nWarnings:\n")
    if (length(warnings_vec) == 0) {
      cat("  None\n")
    } else {
      for (w in warnings_vec) {
        cat("  - ", w, "\n", sep = "")
      }
    }

    cat("\nFailure reason: ", msg, "\n", sep = "")
    cat("Final status: FAIL\n")
  }

  stop(msg, call. = FALSE)
}

add_warning <- function(msg) {
  warnings_vec <<- c(warnings_vec, msg)
  warning(msg, call. = FALSE)
}

read_csv_base <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
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
    fail_clear(
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

  genes <- trimws(as.character(signature_df[[gene_col]]))
  genes <- genes[!is.na(genes) & genes != ""]

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

is_default_rownames <- function(x) {
  rn <- rownames(x)

  is.null(rn) ||
    identical(rn, as.character(seq_len(nrow(x))))
}

coerce_expression_matrix <- function(expr_obj) {
  expr_class <- paste(class(expr_obj), collapse = ",")

  if (is.matrix(expr_obj)) {
    expr_mat <- expr_obj

  } else if (is.data.frame(expr_obj)) {
    expr_df <- expr_obj

    if (is_default_rownames(expr_df) && ncol(expr_df) >= 2) {
      first_col_chr <- trimws(as.character(expr_df[[1]]))

      first_col_nonempty <- !any(
        is.na(first_col_chr) |
          first_col_chr == ""
      )

      first_col_numeric <- suppressWarnings(
        all(!is.na(as.numeric(first_col_chr)))
      )

      first_col_unique <- !anyDuplicated(first_col_chr)

      if (
        first_col_nonempty &&
        !first_col_numeric &&
        first_col_unique
      ) {
        rownames(expr_df) <- first_col_chr
        expr_df <- expr_df[, -1, drop = FALSE]

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
    any(!grepl("^GSM[0-9]+$", sample_ids)) ||
    any(grepl(forbidden_id_pattern, sample_ids, ignore.case = TRUE)) ||
    anyDuplicated(sample_ids)
  ) {
    fail_clear(
      paste0(
        "Possible incorrect input or cohort contamination: ",
        COHORT_LABEL,
        " expression must contain unique GSM sample identifiers only."
      )
    )
  }

  if (!is.numeric(expr_mat)) {
    fail_clear(
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
    fail_clear(
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
    fail_clear(
      "Expression object could not be converted to matrix."
    )
  }

  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    fail_clear(
      "Expression matrix is empty."
    )
  }

  if (
    is.null(rownames(expr_mat)) ||
    any(is.na(rownames(expr_mat))) ||
    any(rownames(expr_mat) == "")
  ) {
    fail_clear(
      "Expression matrix must have non-empty gene rownames."
    )
  }

  if (
    is.null(colnames(expr_mat)) ||
    any(is.na(colnames(expr_mat))) ||
    any(colnames(expr_mat) == "")
  ) {
    fail_clear(
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
    fail_clear(
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
    fail_clear(
      paste0(
        "No official signature genes found in expression rownames; ",
        "matrix orientation or gene symbols may be wrong."
      )
    )
  }

  invisible(TRUE)
}

require_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      fail_clear(
        paste0(
          "Required package not installed: ",
          pkg
        )
      )
    }
  }
}

format_fdr <- function(x) {
  if (
    is.null(x) ||
    length(x) == 0 ||
    is.na(x)
  ) {
    return("NA")
  }

  if (x < 0.001) {
    return(
      format(
        x,
        scientific = TRUE,
        digits = 2
      )
    )
  }

  sprintf("%.3f", x)
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

extract_constant <- function(
    lines,
    name,
    default,
    is_string = FALSE
) {
  pattern <- paste0(
    "^\\s*",
    name,
    "\\s*<-\\s*"
  )

  hit <- grep(
    pattern,
    lines,
    value = TRUE
  )

  if (length(hit) == 0) {
    return(default)
  }

  rhs <- sub(
    pattern,
    "",
    hit[1]
  )

  rhs <- sub(
    "\\s*#.*$",
    "",
    rhs
  )

  rhs <- trimws(rhs)

  if (is_string) {
    rhs <- gsub(
      "^['\"]|['\"]$",
      "",
      rhs
    )

    return(rhs)
  }

  val <- suppressWarnings(
    as.numeric(rhs)
  )

  ifelse(
    is.na(val),
    default,
    val
  )
}

extract_pdf_dim <- function(
    lines,
    output_name,
    default_width,
    default_height
) {
  idx <- grep(
    output_name,
    lines,
    fixed = TRUE
  )

  if (length(idx) == 0) {
    return(
      list(
        width = default_width,
        height = default_height,
        detected = FALSE
      )
    )
  }

  window <- lines[
    idx[1]:min(length(lines), idx[1] + 12)
  ]

  width_line <- grep(
    "width\\s*=",
    window,
    value = TRUE
  )

  height_line <- grep(
    "height\\s*=",
    window,
    value = TRUE
  )

  width <- default_width
  height <- default_height
  detected <- TRUE

  if (length(width_line) > 0) {
    width <- suppressWarnings(
      as.numeric(
        sub(
          ".*width\\s*=\\s*([0-9.]+).*",
          "\\1",
          width_line[1]
        )
      )
    )
  }

  if (length(height_line) > 0) {
    height <- suppressWarnings(
      as.numeric(
        sub(
          ".*height\\s*=\\s*([0-9.]+).*",
          "\\1",
          height_line[1]
        )
      )
    )
  }

  if (!is.finite(width)) {
    width <- default_width
    detected <- FALSE
  }

  if (!is.finite(height)) {
    height <- default_height
    detected <- FALSE
  }

  list(
    width = width,
    height = height,
    detected = detected
  )
}

is_valid_colour <- function(x) {
  if (
    is.null(x) ||
    length(x) != 1 ||
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(FALSE)
  }

  tryCatch(
    {
      grDevices::col2rgb(x)
      TRUE
    },
    error = function(e) FALSE
  )
}

validate_style_colours <- function(style) {
  colour_fields <- c(
    "col_low",
    "col_high",
    "col_ns",
    "col_down_gene",
    "col_up_gene",
    "col_es"
  )

  for (field in colour_fields) {
    if (!is_valid_colour(style[[field]])) {
      fail_clear(
        paste0(
          "Invalid colour in style$",
          field,
          ": '",
          as.character(style[[field]]),
          "'."
        )
      )
    }
  }

  invisible(TRUE)
}

###############################################################################
# 3. TCGA locked style
###############################################################################

detect_tcga_style <- function() {
  # Exact values frozen from the authoritative TCGA Q1/Q4 figure script.
  style <- list(
    reference_script = "TCGA Q1/Q4 locked-style values frozen in this script",
    reference_detected = FALSE,

    base_family = "Helvetica",

    title_size = 6,
    text_size = 5,
    font_size = 5,

    line_size = 0.20,
    grid_lw = 0.20,
    border_lw = 0.20,
    point_size = 0.55,

    col_low = "#044C87",
    col_high = "#EF711D",

    col_ns = "gray80",
    col_down_gene = "#1C5DA3",
    col_up_gene = "#DB0C07",

    col_es = "#95B300",

    pca_width = 4.2,
    pca_height = 3.8,

    volcano_width = 4.1,
    volcano_height = 3.7,

    nes_width = 5.7,
    nes_height_min = 5.5,

    emt_width = 4.2,
    emt_height = 3.9,

    fallbacks = character(0)
  )

  validate_style_colours(style)
  style
}

###############################################################################
# 4. Themes and save helpers
###############################################################################

standard_theme_bw <- function(style) {
  ggplot2::theme_bw(
    base_size = style$text_size,
    base_family = style$base_family
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = style$base_family,
        size = style$text_size,
        color = "black"
      ),
      legend.position = "none",
      axis.title = ggplot2::element_text(
        size = 4.6,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = 4.3,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.20,
        color = "black"
      ),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.20
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        linewidth = 0.20,
        colour = "black",
        fill = NA
      ),
      plot.margin = grid::unit(
        c(2, 2, 2, 2),
        "mm"
      )
    )
}

standard_theme_classic <- function(
    style,
    margin_mm = c(2, 2, 2, 2)
) {
  ggplot2::theme_classic(
    base_size = style$text_size
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = style$base_family,
        size = style$text_size
      ),
      legend.position = "none",
      axis.line = ggplot2::element_line(
        linewidth = 0.16,
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        size = style$text_size
      ),
      axis.text = ggplot2::element_text(
        size = style$text_size,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = style$line_size
      ),
      axis.ticks.length = grid::unit(
        0.4,
        "mm"
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = grid::unit(
        margin_mm,
        "mm"
      )
    )
}

theme_nature_emt <- function(style) {
  ggplot2::theme_classic(
    base_size = 5,
    base_family = style$base_family
  ) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.16,
        colour = "grey90"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.16
      ),
      axis.ticks.length = grid::unit(
        0.5,
        "mm"
      ),
      panel.border = ggplot2::element_rect(
        linewidth = 0.20,
        colour = "black",
        fill = NA
      ),
      plot.title = ggplot2::element_text(
        face = "plain",
        size = 5,
        hjust = 0.5,
        margin = ggplot2::margin(b = 0)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 5,
        hjust = 0
      ),
      axis.title = ggplot2::element_text(
        size = 4.35
      ),
      axis.text = ggplot2::element_text(
        size = 4.35
      ),
      plot.margin = grid::unit(
        c(1, 1, 1, 1),
        "mm"
      )
    )
}

save_pdf <- function(
    plot,
    path,
    width,
    height
) {
  ggplot2::ggsave(
    path,
    plot,
    width = width,
    height = height,
    units = "cm",
    device = "pdf"
  )

  figures_written <<- c(
    figures_written,
    path
  )
}

###############################################################################
# 5. EMT genes and plot builders
###############################################################################

load_emt_genes <- function() {
  hallmark <- tryCatch(
    msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "H"
    ),
    error = function(e1) {
      tryCatch(
        msigdbr::msigdbr(
          species = "Homo sapiens",
          collection = "H"
        ),
        error = function(e2) {
          fail_clear(
            paste0(
              "Could not load MSigDB Hallmark gene sets. category error: ",
              conditionMessage(e1),
              "; collection error: ",
              conditionMessage(e2)
            )
          )
        }
      )
    }
  )

  if (
    !all(
      c("gs_name", "gene_symbol") %in%
      names(hallmark)
    )
  ) {
    fail_clear(
      paste0(
        "msigdbr Hallmark table must contain columns ",
        "gs_name and gene_symbol."
      )
    )
  }

  unique(
    as.character(
      hallmark$gene_symbol[
        hallmark$gs_name ==
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
      ]
    )
  )
}

build_volcano <- function(
    limma_df,
    style
) {
  padj_thr <- 0.05
  lfc_thr <- 1

  plot_df <- limma_df

  plot_df$neg_log10_padj <- -log10(
    pmax(
      plot_df$adj.P.Val,
      .Machine$double.xmin
    )
  )

  plot_df$Significance <- "Not significant"

  plot_df$Significance[
    plot_df$adj.P.Val < padj_thr &
      plot_df$logFC > lfc_thr
  ] <- "Upregulated"

  plot_df$Significance[
    plot_df$adj.P.Val < padj_thr &
      plot_df$logFC < -lfc_thr
  ] <- "Downregulated"

  plot_df$Significance <- factor(
    plot_df$Significance,
    levels = c(
      "Downregulated",
      "Not significant",
      "Upregulated"
    )
  )

  xmax <- max(
    abs(plot_df$logFC),
    na.rm = TRUE
  )

  xlim_val <- max(
    3,
    min(
      8,
      1.1 * xmax
    )
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = logFC,
      y = neg_log10_padj
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = Significance
      ),
      size = 0.42,
      alpha = 0.6,
      shape = 16
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -lfc_thr,
        lfc_thr
      ),
      linetype = "dashed",
      linewidth = style$line_size,
      alpha = 0.4,
      color = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_thr),
      linetype = "dashed",
      linewidth = style$line_size,
      alpha = 0.4,
      color = "black"
    ) +
    ggplot2::scale_color_manual(
      breaks = c(
        "Downregulated",
        "Not significant",
        "Upregulated"
      ),
      values = c(
        "Downregulated" = style$col_down_gene,
        "Not significant" = style$col_ns,
        "Upregulated" = style$col_up_gene
      ),
      name = NULL,
      labels = c(
        "Downregulated",
        "NS",
        "Upregulated"
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        -xlim_val,
        xlim_val
      )
    ) +
    ggplot2::labs(
      x = expression(log[2]~"fold change"),
      y = expression(-log[10]~"FDR")
    ) +
    standard_theme_bw(style)
}

build_nes_barplot <- function(
    gsea_full,
    style
) {
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

  fg <- gsea_full[
    !is.na(gsea_full$padj) &
      gsea_full$padj < 0.05 &
      gsea_full$pathway %in% keep_main_ids,
    ,
    drop = FALSE
  ]

  if (nrow(fg) == 0) {
    fail_clear(
      paste0(
        "No selected extended Hallmark pathways ",
        "passed FDR < 0.05."
      )
    )
  }

  fg$direction <- ifelse(
    fg$NES > 0,
    "Platelet Score High",
    "Platelet Score Low"
  )

  fg$pathway_clean <- clean_hallmark_label(
    fg$pathway
  )

  # Visual order:
  # positives first, ordered by decreasing NES;
  # negatives last, ordered by decreasing NES so they remain at the bottom.
  fg <- fg[
    order(
      fg$NES < 0,
      -fg$NES
    ),
    ,
    drop = FALSE
  ]

  fg$pathway_clean <- factor(
    fg$pathway_clean,
    levels = rev(fg$pathway_clean)
  )

  x_min_data <- min(
    c(
      fg$NES,
      0
    ),
    na.rm = TRUE
  )

  x_max_data <- max(
    c(
      fg$NES,
      0
    ),
    na.rm = TRUE
  )

  x_min <- floor(x_min_data)
  x_max <- ceiling(x_max_data)

  x_breaks <- sort(
    unique(
      c(
        seq(
          0,
          x_min,
          by = -2
        ),
        seq(
          0,
          x_max,
          by = 2
        )
      )
    )
  )

  ggplot2::ggplot(
    fg,
    ggplot2::aes(
      x = NES,
      y = pathway_clean,
      fill = direction
    )
  ) +
    ggplot2::geom_col(
      width = 0.62
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.18,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Platelet Score High" = style$col_up_gene,
        "Platelet Score Low" = style$col_down_gene
      ),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        x_min,
        x_max
      ),
      breaks = x_breaks,
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    ggplot2::labs(
      x = "Normalized Enrichment\nScore (NES)",
      y = NULL,
      title = NULL
    ) +
    ggplot2::theme_classic(
      base_family = style$base_family,
      base_size = style$text_size
    ) +
    ggplot2::theme(
      legend.position = "none",

      axis.line.x = ggplot2::element_line(
        linewidth = 0.18,
        color = "black"
      ),

      axis.line.y = ggplot2::element_line(
        linewidth = 0.18,
        color = "black"
      ),

      axis.title.x = ggplot2::element_text(
        size = 4.2
      ),

      axis.title.y = ggplot2::element_blank(),

      axis.text.x = ggplot2::element_text(
        size = 4,
        color = "black"
      ),

      axis.text.y = ggplot2::element_text(
        size = 4.5,
        color = "black"
      ),

      axis.ticks.x = ggplot2::element_line(
        linewidth = 0.13,
        color = "black"
      ),

      axis.ticks.y = ggplot2::element_line(
        linewidth = 0.13,
        color = "black"
      ),

      axis.ticks.length = grid::unit(
        0.35,
        "mm"
      ),

      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      plot.margin = ggplot2::margin(
        0.7,
        1.0,
        0.9,
        1.6,
        "mm"
      )
    )
}

build_pca_minimal <- function(
    expr_mat,
    score_df,
    signature_genes,
    style
) {
  score_q <- score_df[
    score_df$platelet_score_group %in%
      c(
        "LOW_Q1",
        "HIGH_Q4"
      ),
    ,
    drop = FALSE
  ]

  score_q$sample_name <- as.character(
    score_q$sample_name
  )

  score_q$platelet_score_group <- as.character(
    score_q$platelet_score_group
  )

  if (anyDuplicated(score_q$sample_name)) {
    dup_samples <- unique(
      score_q$sample_name[
        duplicated(score_q$sample_name)
      ]
    )

    fail_clear(
      paste0(
        "Duplicate sample_name values in Q1/Q4 score input: ",
        paste(
          dup_samples,
          collapse = ", "
        )
      )
    )
  }

  score_q$Group <- ifelse(
    score_q$platelet_score_group == "LOW_Q1",
    "Platelet Score Low",
    "Platelet Score High"
  )

  score_q$Group <- factor(
    score_q$Group,
    levels = c(
      "Platelet Score Low",
      "Platelet Score High"
    )
  )

  validate_expression_orientation(
    expr_mat,
    score_q$sample_name,
    signature_genes
  )

  if (anyDuplicated(rownames(expr_mat))) {
    dup_genes <- unique(
      rownames(expr_mat)[
        duplicated(rownames(expr_mat))
      ]
    )

    fail_clear(
      paste0(
        "Duplicate gene rownames in expression matrix: ",
        paste(
          utils::head(
            dup_genes,
            25
          ),
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
    score_q$sample_name,
    colnames(expr_mat)
  )

  if (length(missing_samples) > 0) {
    fail_clear(
      paste0(
        "Q1/Q4 samples missing from expression matrix: ",
        paste(
          missing_samples,
          collapse = ", "
        )
      )
    )
  }

  expr_q <- expr_mat[
    ,
    score_q$sample_name,
    drop = FALSE
  ]

  if (
    !identical(
      colnames(expr_q),
      score_q$sample_name
    )
  ) {
    fail_clear(
      paste0(
        "Expression columns could not be aligned ",
        "with Q1/Q4 sample_name."
      )
    )
  }

  n_before_signature <- nrow(expr_q)

  signature_hit <- rownames(expr_q) %in%
    signature_genes

  n_signature_excluded <- sum(signature_hit)

  expr_no_signature <- expr_q[
    !signature_hit,
    ,
    drop = FALSE
  ]

  n_after_signature <- nrow(expr_no_signature)

  keep_finite <- rowSums(
    !is.finite(expr_no_signature)
  ) == 0

  expr_finite <- expr_no_signature[
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

  expr_pca <- expr_finite[
    keep_variance,
    ,
    drop = FALSE
  ]

  if (nrow(expr_pca) < 2) {
    fail_clear(
      paste0(
        "Too few non-signature finite variable genes ",
        "available for PCA."
      )
    )
  }

  if (ncol(expr_pca) < 3) {
    fail_clear(
      "Too few Q1/Q4 samples available for PCA."
    )
  }

  pca <- stats::prcomp(
    t(expr_pca),
    center = TRUE,
    scale. = FALSE
  )

  var_exp <- round(
    (
      pca$sdev^2 /
        sum(pca$sdev^2)
    ) * 100,
    2
  )

  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    sample_name = rownames(pca$x),
    stringsAsFactors = FALSE
  )

  pca_df$Group <- score_q$Group[
    match(
      pca_df$sample_name,
      score_q$sample_name
    )
  ]

  if (any(is.na(pca_df$Group))) {
    fail_clear(
      "PCA metadata join introduced NA group labels."
    )
  }

  x_lab <- max(
    pca_df$PC1,
    na.rm = TRUE
  ) -
    0.005 *
    diff(
      range(
        pca_df$PC1,
        na.rm = TRUE
      )
    )

  y_lab1 <- min(
    pca_df$PC2,
    na.rm = TRUE
  ) -
    0.20 *
    diff(
      range(
        pca_df$PC2,
        na.rm = TRUE
      )
    )

  y_lab2 <- min(
    pca_df$PC2,
    na.rm = TRUE
  ) -
    0.11 *
    diff(
      range(
        pca_df$PC2,
        na.rm = TRUE
      )
    )

  plot <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = Group
    )
  ) +
    ggplot2::geom_point(
      size = 0.9,
      alpha = 0.6,
      shape = 16
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Platelet Score Low" = style$col_low,
        "Platelet Score High" = style$col_high
      ),
      name = NULL
    ) +
    ggplot2::annotate(
      "text",
      x = x_lab,
      y = y_lab2,
      label = "Platelet score low",
      hjust = 1,
      vjust = 0.5,
      color = style$col_low,
      family = style$base_family,
      fontface = "plain",
      size = 4.5 / ggplot2::.pt
    ) +
    ggplot2::annotate(
      "text",
      x = x_lab,
      y = y_lab1,
      label = "Platelet score high",
      hjust = 1,
      vjust = 0.5,
      color = style$col_high,
      family = style$base_family,
      fontface = "plain",
      size = 4.5 / ggplot2::.pt
    ) +
    ggplot2::labs(
      x = paste0(
        "PC1 (",
        var_exp[1],
        "%)"
      ),
      y = paste0(
        "PC2 (",
        var_exp[2],
        "%)"
      )
    ) +
    ggplot2::theme_bw(
      base_size = style$text_size,
      base_family = style$base_family
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = style$base_family,
        size = style$text_size,
        color = "black"
      ),
      legend.position = "none",
      axis.title = ggplot2::element_text(
        size = 4.5,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = 4.3,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.20,
        color = "black"
      ),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.14,
        color = "grey95"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        linewidth = 0.20,
        colour = "black",
        fill = NA
      ),
      plot.margin = grid::unit(
        c(2, 2, 2, 2),
        "mm"
      )
    )

  list(
    plot = plot,
    qc = list(
      n_samples = nrow(score_q),
      n_low = sum(
        score_q$platelet_score_group ==
          "LOW_Q1"
      ),
      n_high = sum(
        score_q$platelet_score_group ==
          "HIGH_Q4"
      ),
      n_genes_before_signature_exclusion =
        n_before_signature,
      n_signature_genes_excluded =
        n_signature_excluded,
      n_genes_after_signature_exclusion =
        n_after_signature,
      n_genes_removed_nonfinite =
        sum(!keep_finite),
      n_genes_removed_zero_variance =
        sum(!keep_variance),
      n_genes_used = nrow(expr_pca),
      pc1_variance = var_exp[1],
      pc2_variance = var_exp[2],
      output_path = output_pca
    )
  )
}

build_emt_gsea_plot <- function(
    limma_df,
    emt_row,
    style
) {
  rnk <- limma_df[
    !is.na(limma_df$gene_name) &
      limma_df$gene_name != "" &
      is.finite(limma_df$t_stat),
    c(
      "gene_name",
      "t_stat"
    ),
    drop = FALSE
  ]

  if (anyDuplicated(rnk$gene_name)) {
    fail_clear(
      paste0(
        "Cannot reconstruct EMT GSEA curve because ",
        "ranking contains duplicated gene_name values."
      )
    )
  }

  stats <- stats::setNames(
    rnk$t_stat,
    rnk$gene_name
  )

  stats <- sort(
    stats[
      is.finite(stats)
    ],
    decreasing = TRUE
  )

  emt_genes <- load_emt_genes()

  emt_genes_in <- intersect(
    emt_genes,
    names(stats)
  )

  Nh <- length(emt_genes_in)
  N <- length(stats)

  if (Nh < 10) {
    fail_clear(
      paste0(
        "Cannot reconstruct EMT GSEA curve; ",
        "too few EMT genes in ranking: ",
        Nh
      )
    )
  }

  hits <- which(
    names(stats) %in%
      emt_genes_in
  )

  hit_set <- logical(N)
  hit_set[hits] <- TRUE

  w_hit <- abs(stats[hits])
  w_hit <- w_hit / sum(w_hit)

  miss_penalty <- 1 / (N - Nh)

  running <- numeric(N)
  rs <- 0
  hit_idx <- 1

  for (i in seq_len(N)) {
    if (hit_set[i]) {
      rs <- rs + w_hit[hit_idx]
      hit_idx <- hit_idx + 1
    } else {
      rs <- rs - miss_penalty
    }

    running[i] <- rs
  }

  ES_pos <- max(
    running,
    na.rm = TRUE
  )

  ES_neg <- min(
    running,
    na.rm = TRUE
  )

  ES <- if (
    abs(ES_pos) >= abs(ES_neg)
  ) {
    ES_pos
  } else {
    ES_neg
  }

  NES_txt <- suppressWarnings(
    as.numeric(
      emt_row$NES[1]
    )
  )

  FDR_txt <- suppressWarnings(
    as.numeric(
      emt_row$padj[1]
    )
  )

  metrics_block <- paste(
    stats::na.omit(
      c(
        if (is.finite(NES_txt)) {
          paste0(
            "NES = ",
            round(
              NES_txt,
              2
            )
          )
        } else {
          paste0(
            "ES = ",
            round(
              ES,
              3
            )
          )
        },
        if (is.finite(FDR_txt)) {
          paste0(
            "FDR = ",
            format_fdr(FDR_txt)
          )
        } else {
          NA_character_
        }
      )
    ),
    collapse = "\n"
  )

  df_es <- data.frame(
    rank = seq_len(N),
    ES = running
  )

  df_hits <- data.frame(
    rank = hits
  )

  df_metric <- data.frame(
    rank = seq_len(N),
    metric = as.numeric(stats)
  )

  q_lo <- as.numeric(
    stats::quantile(
      df_metric$metric,
      0.01,
      na.rm = TRUE
    )
  )

  q_hi <- as.numeric(
    stats::quantile(
      df_metric$metric,
      0.99,
      na.rm = TRUE
    )
  )

  df_metric$metric_clip <- pmax(
    q_lo,
    pmin(
      q_hi,
      df_metric$metric
    )
  )

  df_strip <- data.frame(
    rank = df_metric$rank,
    y = 0,
    z = df_metric$metric_clip
  )

  p_es <- ggplot2::ggplot(
    df_es,
    ggplot2::aes(
      rank,
      ES
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.16,
      color = "black",
      alpha = 0.6
    ) +
    ggplot2::geom_line(
      linewidth = 0.75,
      color = style$col_es,
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
      family = style$base_family,
      size = style$font_size / ggplot2::.pt
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    theme_nature_emt(style) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),

      plot.title = ggplot2::element_text(
        face = "plain",
        size = 5,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 0.40,
          unit = "mm"
        )
      ),

      plot.margin = grid::unit(
        c(
          0.8,
          1,
          0.5,
          1
        ),
        "mm"
      )
    )

  p_hits <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df_hits,
      ggplot2::aes(
        x = rank,
        xend = rank,
        y = 1,
        yend = 0
      ),
      linewidth = 0.10,
      color = "black",
      alpha = 0.60
    ) +
    ggplot2::geom_tile(
      data = df_strip,
      ggplot2::aes(
        x = rank,
        y = -0.33,
        fill = z
      ),
      height = 0.55,
      width = 1
    ) +
    ggplot2::scale_fill_gradient2(
      low = style$col_down_gene,
      mid = "white",
      high = style$col_up_gene,
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
      family = style$base_family,
      fontface = "plain",
      size = 4 / ggplot2::.pt
    ) +
    ggplot2::annotate(
      "text",
      x = N,
      y = -0.92,
      label = "Platelet score low",
      hjust = 1,
      vjust = 0.5,
      family = style$base_family,
      fontface = "plain",
      size = 4 / ggplot2::.pt
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        1,
        N
      ),
      ylim = c(
        -0.95,
        1.02
      ),
      clip = "off"
    ) +
    ggplot2::theme_void(
      base_family = style$base_family,
      base_size = style$font_size
    ) +
    ggplot2::theme(
      plot.margin = grid::unit(
        c(
          0.2,
          1,
          0.2,
          1
        ),
        "mm"
      )
    )

  p_metric <- ggplot2::ggplot(
    df_metric,
    ggplot2::aes(
      rank,
      metric
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = style$line_size,
      color = "black",
      alpha = 0.75
    ) +
    ggplot2::geom_area(
      fill = "grey85",
      alpha = 1
    ) +
    ggplot2::geom_line(
      linewidth = 0.20,
      color = "grey55"
    ) +
    ggplot2::labs(
      x = "Rank in ordered dataset",
      y = "Ranking metric\n(limma t statistic)"
    ) +
    theme_nature_emt(style) +
    ggplot2::theme(
      plot.margin = grid::unit(
        c(
          0.5,
          1,
          1,
          1
        ),
        "mm"
      )
    )

  list(
    p_es = p_es,
    p_hits = p_hits,
    p_metric = p_metric,
    n_emt_genes_in_ranking = Nh
  )
}

###############################################################################
# 6. QC
###############################################################################

write_qc_report <- function(
    path,
    qc
) {
  sink(path)

  cat("Friedrich/GSE134051 Q1/Q4 figures QC\n")
  cat("==========================\n\n")

  cat("Inputs found:\n")

  for (nm in names(qc$inputs_found)) {
    cat(
      "  ",
      nm,
      ": ",
      qc$inputs_found[[nm]],
      "\n",
      sep = ""
    )
  }

  cat("\nVolcano:\n")
  cat(
    "  n_genes: ",
    qc$volcano$n_genes,
    "\n",
    sep = ""
  )
  cat(
    "  n_significant: ",
    qc$volcano$n_significant,
    "\n",
    sep = ""
  )
  cat(
    "  n_up_HIGH_Q4: ",
    qc$volcano$n_up_high,
    "\n",
    sep = ""
  )
  cat(
    "  n_down_LOW_Q1: ",
    qc$volcano$n_down_low,
    "\n",
    sep = ""
  )

  cat("\nPCA:\n")
  cat(
    "  n_pca_samples: ",
    qc$pca$n_pca_samples,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_LOW_Q1: ",
    qc$pca$n_pca_LOW_Q1,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_HIGH_Q4: ",
    qc$pca$n_pca_HIGH_Q4,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_genes_before_signature_exclusion: ",
    qc$pca$n_pca_genes_before_signature_exclusion,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_signature_genes_excluded: ",
    qc$pca$n_pca_signature_genes_excluded,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_genes_after_signature_exclusion: ",
    qc$pca$n_pca_genes_after_signature_exclusion,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_genes_used: ",
    qc$pca$n_pca_genes_used,
    "\n",
    sep = ""
  )
  cat(
    "  PC1 variance percent: ",
    qc$pca$pc1_variance,
    "\n",
    sep = ""
  )
  cat(
    "  PC2 variance percent: ",
    qc$pca$pc2_variance,
    "\n",
    sep = ""
  )
  cat(
    "  PCA output path: ",
    qc$pca$pca_output_path,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_genes_removed_nonfinite: ",
    qc$pca$n_pca_genes_removed_nonfinite,
    "\n",
    sep = ""
  )
  cat(
    "  n_pca_genes_removed_zero_variance: ",
    qc$pca$n_pca_genes_removed_zero_variance,
    "\n",
    sep = ""
  )

  cat("\nHallmark:\n")
  cat(
    "  n_pathways: ",
    qc$gsea$n_pathways,
    "\n",
    sep = ""
  )
  cat(
    "  n_FDR005: ",
    qc$gsea$n_fdr005,
    "\n",
    sep = ""
  )

  cat("\nEMT:\n")
  cat(
    "  EMT_NES: ",
    qc$emt$NES,
    "\n",
    sep = ""
  )
  cat(
    "  EMT_padj: ",
    qc$emt$padj,
    "\n",
    sep = ""
  )
  cat(
    "  EMT_direction: ",
    qc$emt$direction,
    "\n",
    sep = ""
  )
  cat(
    "  EMT genes in ranking: ",
    qc$emt$n_genes_in_ranking,
    "\n",
    sep = ""
  )

  cat("\nFigures written:\n")

  if (length(qc$figures_written) == 0) {
    cat("  None\n")
  } else {
    for (fp in qc$figures_written) {
      cat(
        "  - ",
        fp,
        "\n",
        sep = ""
      )
    }
  }

  cat("\nTCGA reference:\n")
  cat(
    "  Reference script: ",
    qc$tcga_reference$script,
    "\n",
    sep = ""
  )
  cat(
    "  Reference detected: ",
    qc$tcga_reference$detected,
    "\n",
    sep = ""
  )
  cat(
    "  PCA dimensions cm: ",
    qc$tcga_reference$pca_width,
    " x ",
    qc$tcga_reference$pca_height,
    "\n",
    sep = ""
  )
  cat(
    "  Volcano dimensions cm: ",
    qc$tcga_reference$volcano_width,
    " x ",
    qc$tcga_reference$volcano_height,
    "\n",
    sep = ""
  )
  cat(
    "  Hallmark dimensions cm: ",
    qc$tcga_reference$nes_width,
    " x dynamic height\n",
    sep = ""
  )
  cat(
    "  EMT dimensions cm: ",
    qc$tcga_reference$emt_width,
    " x ",
    qc$tcga_reference$emt_height,
    "\n",
    sep = ""
  )

  cat("\nLocked TCGA colours used:\n")
  cat(
    "  COL_LOW: ",
    qc$tcga_colours$col_low,
    "\n",
    sep = ""
  )
  cat(
    "  COL_HIGH: ",
    qc$tcga_colours$col_high,
    "\n",
    sep = ""
  )
  cat(
    "  COL_NS: ",
    qc$tcga_colours$col_ns,
    "\n",
    sep = ""
  )
  cat(
    "  COL_DOWN_GENE: ",
    qc$tcga_colours$col_down_gene,
    "\n",
    sep = ""
  )
  cat(
    "  COL_UP_GENE: ",
    qc$tcga_colours$col_up_gene,
    "\n",
    sep = ""
  )
  cat(
    "  COL_ES: ",
    qc$tcga_colours$col_es,
    "\n",
    sep = ""
  )

  cat("\nPackages used:\n")
  cat(
    "  ggplot2, jsonlite, msigdbr, patchwork, grid/base\n"
  )

  cat("\nWarnings:\n")

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
    qc$qc_status,
    "\n",
    sep = ""
  )

  sink()
}

###############################################################################
# 7. Main
###############################################################################

require_packages(
  c(
    "ggplot2",
    "jsonlite",
    "msigdbr",
    "patchwork"
  )
)

critical_inputs <- c(
  input_score,
  input_expression,
  input_signature,
  input_limma,
  input_gsea_full,
  input_gsea_emt,
  input_gsea_metadata
)

missing_inputs <- critical_inputs[
  !file.exists(critical_inputs)
]

if (length(missing_inputs) > 0) {
  fail_clear(
    paste0(
      "Missing critical input(s): ",
      paste(
        missing_inputs,
        collapse = ", "
      )
    )
  )
}

score_df <- read_csv_base(
  input_score
)

signature_df <- read_canonical_platelet_signature()

signature_genes <- extract_signature_genes(
  signature_df
)

expr_obj <- readRDS(
  input_expression
)

expr_mat <- coerce_expression_matrix(
  expr_obj
)
validate_friedrich_expression_identity(expr_mat)

limma_df <- read_csv_base(
  input_limma
)

gsea_full <- read_csv_base(
  input_gsea_full
)

gsea_emt <- read_csv_base(
  input_gsea_emt
)

gsea_metadata <- jsonlite::read_json(
  input_gsea_metadata,
  simplifyVector = TRUE
)

message("[2/6] Loading frozen graphical style")

style <- detect_tcga_style()

if (length(style$fallbacks) > 0) {
  for (w in style$fallbacks) {
    add_warning(w)
  }
}

message("[3/6] Preparing plot data")

required_limma <- c(
  "gene_name",
  "logFC",
  "adj.P.Val",
  "t_stat",
  "is_signature_gene"
)

missing_limma <- setdiff(
  required_limma,
  names(limma_df)
)

if (length(missing_limma) > 0) {
  fail_clear(
    paste0(
      "Volcano input missing required columns: ",
      paste(
        missing_limma,
        collapse = ", "
      )
    )
  )
}

required_gsea <- c(
  "pathway",
  "NES",
  "padj"
)

missing_gsea <- setdiff(
  required_gsea,
  names(gsea_full)
)

if (length(missing_gsea) > 0) {
  fail_clear(
    paste0(
      "GSEA full input missing required columns: ",
      paste(
        missing_gsea,
        collapse = ", "
      )
    )
  )
}

missing_emt <- setdiff(
  required_gsea,
  names(gsea_emt)
)

if (length(missing_emt) > 0) {
  fail_clear(
    paste0(
      "EMT row input missing required columns: ",
      paste(
        missing_emt,
        collapse = ", "
      )
    )
  )
}

if (
  !"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" %in%
  gsea_emt$pathway
) {
  fail_clear(
    paste0(
      "EMT row input does not contain ",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION."
    )
  )
}

if (
  !all(
    c(
      "sample_name",
      "Score_z",
      "platelet_score_group"
    ) %in% names(score_df)
  )
) {
  fail_clear(
    paste0(
      "Score input must contain sample_name, ",
      "Score_z, and platelet_score_group."
    )
  )
}

score_df$sample_name <- trimws(
  as.character(score_df$sample_name)
)

forbidden_id_pattern <- paste0("^ICGC", "_PCA")
if (
  any(!grepl("^GSM[0-9]+$", score_df$sample_name)) ||
  any(grepl(forbidden_id_pattern, score_df$sample_name, ignore.case = TRUE)) ||
  anyDuplicated(score_df$sample_name)
) {
  fail_clear(
    paste0(
      "Possible incorrect input or cohort contamination: score sample_name ",
      "values must be unique GSM identifiers from ",
      COHORT_LABEL,
      "."
    )
  )
}

limma_df$gene_name <- as.character(
  limma_df$gene_name
)

limma_df$logFC <- suppressWarnings(
  as.numeric(limma_df$logFC)
)

limma_df$adj.P.Val <- suppressWarnings(
  as.numeric(limma_df$adj.P.Val)
)

limma_df$t_stat <- suppressWarnings(
  as.numeric(limma_df$t_stat)
)

gsea_full$pathway <- as.character(
  gsea_full$pathway
)

gsea_full$NES <- suppressWarnings(
  as.numeric(gsea_full$NES)
)

gsea_full$padj <- suppressWarnings(
  as.numeric(gsea_full$padj)
)

gsea_emt$pathway <- as.character(
  gsea_emt$pathway
)

gsea_emt$NES <- suppressWarnings(
  as.numeric(gsea_emt$NES)
)

gsea_emt$padj <- suppressWarnings(
  as.numeric(gsea_emt$padj)
)

score_df$Score_z <- suppressWarnings(
  as.numeric(score_df$Score_z)
)

if (all(!is.finite(limma_df$logFC))) {
  fail_clear(
    "All logFC values are non-finite."
  )
}

if (all(!is.finite(limma_df$adj.P.Val))) {
  fail_clear(
    "All adj.P.Val values are non-finite."
  )
}

if (all(!is.finite(limma_df$t_stat))) {
  fail_clear(
    "All t_stat values are non-finite."
  )
}

volcano_sig <- limma_df$adj.P.Val < 0.05 &
  abs(limma_df$logFC) > 1

n_volcano_sig <- sum(
  volcano_sig,
  na.rm = TRUE
)

n_volcano_up <- sum(
  volcano_sig &
    limma_df$logFC > 0,
  na.rm = TRUE
)

n_volcano_down <- sum(
  volcano_sig &
    limma_df$logFC < 0,
  na.rm = TRUE
)

n_gsea_fdr <- sum(
  !is.na(gsea_full$padj) &
    gsea_full$padj < 0.05
)

emt_row <- gsea_emt[
  gsea_emt$pathway ==
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  ,
  drop = FALSE
]

emt_direction <- if (
  "direction" %in% names(emt_row)
) {
  as.character(
    emt_row$direction[1]
  )
} else if (emt_row$NES[1] > 0) {
  "HIGH_Q4"
} else if (emt_row$NES[1] < 0) {
  "LOW_Q1"
} else {
  "ZERO"
}

message("[4/6] Building figures")

pca_result <- build_pca_minimal(
  expr_mat,
  score_df,
  signature_genes,
  style
)

p_volcano <- build_volcano(
  limma_df,
  style
)

p_nes <- build_nes_barplot(
  gsea_full,
  style
)

emt_plots <- build_emt_gsea_plot(
  limma_df,
  emt_row,
  style
)

message("[5/6] Writing outputs")

save_pdf(
  pca_result$plot,
  output_pca,
  style$pca_width,
  style$pca_height
)

save_pdf(
  p_volcano,
  output_volcano,
  style$volcano_width,
  style$volcano_height
)

keep_main_ids_for_height <- c(
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

n_selected_nes <- sum(
  !is.na(gsea_full$padj) &
    gsea_full$padj < 0.05 &
    gsea_full$pathway %in%
    keep_main_ids_for_height
)

nes_height <- max(
  3.95,
  0.22 *
    max(
      1,
      n_selected_nes
    ) +
    1
)

save_pdf(
  p_nes,
  output_nes,
  4.9,
  nes_height
)

p_emt_gsea <- (
  (
    emt_plots$p_es /
      emt_plots$p_hits /
      emt_plots$p_metric
  ) +
    patchwork::plot_layout(
      heights = c(
        2.0,
        0.55,
        0.85
      )
    )
) &
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      0.6,
      0.6,
      0.6,
      0.6,
      "pt"
    )
  )

ggplot2::ggsave(
  output_emt,
  p_emt_gsea,
  width = style$emt_width,
  height = style$emt_height,
  units = "cm",
  device = "pdf"
)

figures_written <- c(
  figures_written,
  output_emt
)

qc_status <- "PASS"

if (length(warnings_vec) > 0) {
  qc_status <- "WARNING"
}

if (
  nrow(limma_df) == 0 ||
  nrow(gsea_full) == 0 ||
  nrow(emt_row) == 0
) {
  qc_status <- "FAIL"
}

qc <- list(
  inputs_found = list(
    score = file.exists(input_score),
    expression = file.exists(input_expression),
    signature = file.exists(input_signature),
    limma_no_signature = file.exists(input_limma),
    gsea_full = file.exists(input_gsea_full),
    gsea_emt = file.exists(input_gsea_emt),
    gsea_metadata = file.exists(input_gsea_metadata)
  ),
  volcano = list(
    n_genes = nrow(limma_df),
    n_significant = n_volcano_sig,
    n_up_high = n_volcano_up,
    n_down_low = n_volcano_down
  ),
  pca = list(
    n_pca_samples =
      pca_result$qc$n_samples,
    n_pca_LOW_Q1 =
      pca_result$qc$n_low,
    n_pca_HIGH_Q4 =
      pca_result$qc$n_high,
    n_pca_genes_before_signature_exclusion =
      pca_result$qc$n_genes_before_signature_exclusion,
    n_pca_signature_genes_excluded =
      pca_result$qc$n_signature_genes_excluded,
    n_pca_genes_after_signature_exclusion =
      pca_result$qc$n_genes_after_signature_exclusion,
    n_pca_genes_removed_nonfinite =
      pca_result$qc$n_genes_removed_nonfinite,
    n_pca_genes_removed_zero_variance =
      pca_result$qc$n_genes_removed_zero_variance,
    n_pca_genes_used =
      pca_result$qc$n_genes_used,
    pc1_variance =
      pca_result$qc$pc1_variance,
    pc2_variance =
      pca_result$qc$pc2_variance,
    pca_output_path =
      pca_result$qc$output_path
  ),
  gsea = list(
    n_pathways = nrow(gsea_full),
    n_fdr005 = n_gsea_fdr
  ),
  emt = list(
    NES = emt_row$NES[1],
    padj = emt_row$padj[1],
    direction = emt_direction,
    n_genes_in_ranking =
      emt_plots$n_emt_genes_in_ranking
  ),
  figures_written = figures_written,
  tcga_reference = list(
    script = style$reference_script,
    detected = style$reference_detected,
    pca_width = style$pca_width,
    pca_height = style$pca_height,
    volcano_width = style$volcano_width,
    volcano_height = style$volcano_height,
    nes_width = style$nes_width,
    emt_width = style$emt_width,
    emt_height = style$emt_height
  ),
  tcga_colours = list(
    col_low = style$col_low,
    col_high = style$col_high,
    col_ns = style$col_ns,
    col_down_gene = style$col_down_gene,
    col_up_gene = style$col_up_gene,
    col_es = style$col_es
  ),
  gsea_metadata = gsea_metadata,
  warnings = warnings_vec,
  qc_status = qc_status
)

write_qc_report(
  output_qc,
  qc
)

if (qc_status == "FAIL") {
  stop(
    "Final status: FAIL",
    call. = FALSE
  )
}

message(
  "[6/6] Final status: ",
  qc_status
)
