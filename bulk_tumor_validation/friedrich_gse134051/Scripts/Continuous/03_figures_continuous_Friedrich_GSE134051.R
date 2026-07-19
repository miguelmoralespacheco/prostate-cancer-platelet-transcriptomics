#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 | Continuous EMT GSEA figure
#
# Generates only the EMT GSEA enrichment plot for the complementary continuous
# platelet-associated transcriptional score analysis.
#
# Input rank:
#   limma t_stat from the complete no-signature continuous results.
#
# The graphical construction, colours, typography, dimensions, line widths,
# margins and panel proportions are preserved from the previously approved
# Friedrich/GSE134051 continuous GSEA figure.
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

input_limma <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables",
  "LIMMA",
  "LIMMA_continuous_41genes_no_signature_genes.csv"
)

input_gsea_emt <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables",
  "GSEA",
  "GSEA_continuous_no_signature_Hallmark_EMT_row.csv"
)

input_gsea_metadata <- file.path(
  GENERATED_RESULTS_DIR, "Continuous",
  "Tables",
  "GSEA",
  "GSEA_continuous_metadata.json"
)

figures_dir <- file.path(
  GENERATED_FIGURES_DIR, "Continuous"
)

gsea_figures_dir <- file.path(
  figures_dir,
  "GSEA"
)

logs_dir <- file.path(
  GENERATED_RESULTS_DIR, "Continuous", "Logs"
)

output_emt <- file.path(
  gsea_figures_dir,
  "Friedrich_GSE134051_continuous_EMT_GSEA_plot.pdf"
)

output_qc <- file.path(
  logs_dir,
  "Friedrich_GSE134051_continuous_EMT_GSEA_figure_QC.txt"
)

dir.create(
  gsea_figures_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  logs_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

warnings_vec <- character(0)
figures_written <- character(0)

###############################################################################
# 2. General helpers
###############################################################################

fail_with_qc <- function(msg) {
  writeLines(
    c(
      "Friedrich/GSE134051 continuous EMT GSEA figure QC",
      "========================================",
      "",
      paste0(
        "Date: ",
        as.character(Sys.time())
      ),
      paste0(
        "Project dir: ",
        project_dir
      ),
      paste0(
        "Failure reason: ",
        msg
      ),
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
  warnings_vec <<- unique(
    c(
      warnings_vec,
      msg
    )
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
        paste(
          missing,
          collapse = ", "
        )
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

  sprintf(
    "%.3f",
    x
  )
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
    idx[1]:min(
      length(lines),
      idx[1] + 12
    )
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

valid_colour <- function(x) {
  tryCatch(
    {
      is.character(x) &&
        length(x) == 1 &&
        !is.na(x) &&
        nzchar(x) &&
        !is.na(
          grDevices::col2rgb(x)[1]
        )
    },
    error = function(e) FALSE
  )
}

save_pdf <- function(
    plot,
    path,
    width,
    height,
    family = "Helvetica"
) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    device = grDevices::pdf,
    family = family,
    useDingbats = FALSE,
    bg = "transparent"
  )

  figures_written <<- unique(
    c(
      figures_written,
      path
    )
  )

  invisible(path)
}

###############################################################################
# 3. Locked graphical style
###############################################################################

detect_style <- function() {
  # Exact values frozen from the authoritative TCGA locked-style source.
  list(
    reference_script = "TCGA locked-style values frozen in this script",
    reference_type = "embedded exact values",
    reference_detected = FALSE,

    base_family = "Helvetica",
    title_size = 6,
    text_size = 5,
    font_size = 5,

    line_size = 0.20,
    grid_lw = 0.20,
    border_lw = 0.20,

    col_ns_light = "grey85",
    col_neg = "#1C5DA3",
    col_pos = "#DB0C07",
    col_emt = "#95B300",

    emt_width = 4.2,
    emt_height = 3.9,

    fallbacks = character(0)
  )
}

theme_emt <- function(style) {
  ggplot2::theme_classic(
    base_size = style$font_size,
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

      plot.title = ggplot2::element_text(
        face = "plain",
        size = 5,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 0
        )
      ),

      axis.title = ggplot2::element_text(
        size = 4.35
      ),

      axis.text = ggplot2::element_text(
        size = 4.35
      ),

      plot.margin = grid::unit(
        c(
          1,
          1,
          1,
          1
        ),
        "mm"
      )
    )
}

###############################################################################
# 4. Ranking and Hallmark helpers
###############################################################################

prepare_rank <- function(limma_df) {
  rank_df <- limma_df[
    !is.na(limma_df$gene_name) &
      limma_df$gene_name != "" &
      is.finite(limma_df$t_stat),
    c(
      "gene_name",
      "t_stat"
    ),
    drop = FALSE
  ]

  rank_df <- rank_df[
    order(
      -abs(rank_df$t_stat),
      rank_df$gene_name
    ),
    ,
    drop = FALSE
  ]

  rank_df <- rank_df[
    !duplicated(rank_df$gene_name),
    ,
    drop = FALSE
  ]

  ranked_stats <- stats::setNames(
    rank_df$t_stat,
    rank_df$gene_name
  )

  sort(
    ranked_stats[
      is.finite(ranked_stats)
    ],
    decreasing = TRUE
  )
}

load_emt_genes <- function() {
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

  if (
    !all(
      c(
        "gs_name",
        "gene_symbol"
      ) %in% names(hallmark)
    )
  ) {
    fail_with_qc(
      paste0(
        "msigdbr Hallmark table must contain ",
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

###############################################################################
# 5. EMT GSEA figure builder
###############################################################################

build_emt_gsea_plot <- function(
    limma_df,
    emt_row,
    style
) {
  ranked_stats <- prepare_rank(
    limma_df
  )

  emt_genes <- load_emt_genes()

  emt_genes_in <- intersect(
    emt_genes,
    names(ranked_stats)
  )

  Nh <- length(
    emt_genes_in
  )

  N <- length(
    ranked_stats
  )

  if (Nh < 10) {
    fail_with_qc(
      paste0(
        "Too few EMT genes in ranking for EMT figure: ",
        Nh
      )
    )
  }

  if (N <= Nh) {
    fail_with_qc(
      paste0(
        "Invalid ranking dimensions: N = ",
        N,
        "; EMT genes = ",
        Nh,
        "."
      )
    )
  }

  hits <- which(
    names(ranked_stats) %in%
      emt_genes_in
  )

  hit_set <- logical(N)
  hit_set[hits] <- TRUE

  w_hit <- abs(
    ranked_stats[hits]
  )

  if (
    !is.finite(sum(w_hit)) ||
    sum(w_hit) == 0
  ) {
    fail_with_qc(
      paste0(
        "EMT hit weights are zero or non-finite; ",
        "running enrichment score cannot be reconstructed."
      )
    )
  }

  w_hit <- w_hit /
    sum(w_hit)

  miss_penalty <- 1 /
    (N - Nh)

  running <- numeric(N)
  rs <- 0
  hit_idx <- 1

  for (i in seq_len(N)) {
    if (hit_set[i]) {
      rs <- rs +
        w_hit[hit_idx]

      hit_idx <- hit_idx + 1
    } else {
      rs <- rs -
        miss_penalty
    }

    running[i] <- rs
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
          NA_character_
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
    metric = as.numeric(ranked_stats)
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

  ###########################################################################
  # Upper panel: running enrichment score
  ###########################################################################

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
      color = style$col_emt,
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
      size = style$font_size /
        ggplot2::.pt
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    theme_emt(style) +
    ggplot2::theme(
      panel.grid.major =
        ggplot2::element_blank(),

      panel.grid.minor =
        ggplot2::element_blank(),

      axis.text.x =
        ggplot2::element_blank(),

      axis.ticks.x =
        ggplot2::element_blank(),

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

  ###########################################################################
  # Middle panel: hit positions and ranked-expression strip
  ###########################################################################

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
      low = style$col_neg,
      mid = "white",
      high = style$col_pos,
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
      size = 4 /
        ggplot2::.pt
    ) +
    ggplot2::annotate(
      "text",
      x = N,
      y = -0.92,
      label = "Platelet score low",
      hjust = 1,
      vjust = 0.5,
      family = style$base_family,
      size = 4 /
        ggplot2::.pt
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
      plot.background = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

      panel.background = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

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

  ###########################################################################
  # Lower panel: ranking metric
  ###########################################################################

  p_metric <- ggplot2::ggplot(
    df_metric,
    ggplot2::aes(
      rank,
      metric
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.16,
      color = "black",
      alpha = 0.75
    ) +
    ggplot2::geom_area(
      fill = style$col_ns_light,
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
    theme_emt(style) +
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

  ###########################################################################
  # Locked three-panel construction
  ###########################################################################

  combined <- (
    p_es /
      p_hits /
      p_metric
  ) +
    patchwork::plot_layout(
      heights = c(
        2.0,
        0.55,
        0.85
      )
    ) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(
          fill = NA,
          colour = NA
        )
      )
    ) &
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

      panel.background = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

      legend.background = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

      legend.box.background = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

      legend.key = ggplot2::element_rect(
        fill = NA,
        colour = NA
      ),

      plot.margin = ggplot2::margin(
        0.6,
        0.6,
        0.6,
        0.6,
        "pt"
      )
    )

  list(
    plot = combined,
    n_emt_genes_in_ranking = Nh,
    n_ranked_genes = N
  )
}

###############################################################################
# 6. QC report
###############################################################################

write_qc_report <- function(qc) {
  sink(output_qc)

  cat("Friedrich/GSE134051 continuous EMT GSEA figure QC\n")
  cat("========================================\n\n")

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

  cat("\n2. Ranking\n")

  cat(
    "  Ranking variable: t_stat\n"
  )

  cat(
    "  Ranking fallback: none\n"
  )

  cat(
    "  Ranked genes: ",
    qc$ranking$n_ranked_genes,
    "\n",
    sep = ""
  )

  cat(
    "  EMT genes in ranking: ",
    qc$ranking$n_emt_genes_in_ranking,
    "\n",
    sep = ""
  )

  cat("\n3. EMT result\n")

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

  cat("\n4. Locked graphical style\n")

  cat(
    "  Reference script: ",
    qc$style$reference_script,
    "\n",
    sep = ""
  )

  cat(
    "  Reference type: ",
    qc$style$reference_type,
    "\n",
    sep = ""
  )

  cat(
    "  Reference detected: ",
    qc$style$reference_detected,
    "\n",
    sep = ""
  )

  cat(
    "  Font family: ",
    qc$style$base_family,
    "\n",
    sep = ""
  )

  cat(
    "  EMT curve colour: ",
    qc$style$col_emt,
    "\n",
    sep = ""
  )

  cat(
    "  Negative ranking colour: ",
    qc$style$col_neg,
    "\n",
    sep = ""
  )

  cat(
    "  Positive ranking colour: ",
    qc$style$col_pos,
    "\n",
    sep = ""
  )

  cat(
    "  EMT dimensions cm: ",
    qc$style$emt_width,
    " x ",
    qc$style$emt_height,
    "\n",
    sep = ""
  )

  cat(
    "  Patchwork heights: 2.0, 0.55, 0.85\n"
  )

  cat("\n5. Output\n")

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

  cat("\n6. Warnings\n")

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
# 7. Load and validate inputs
###############################################################################

require_packages(
  c(
    "ggplot2",
    "jsonlite",
    "msigdbr",
    "patchwork"
  )
)

for (
  fp in c(
    input_limma,
    input_gsea_emt,
    input_gsea_metadata
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

gsea_emt <- read_csv_base(
  input_gsea_emt
)

gsea_metadata <- jsonlite::read_json(
  input_gsea_metadata,
  simplifyVector = TRUE
)

if (
  !is.null(gsea_metadata$ranking_variable) &&
  !identical(
    as.character(
      gsea_metadata$ranking_variable
    ),
    "t_stat"
  )
) {
  fail_with_qc(
    paste0(
      "GSEA metadata ranking_variable is not t_stat: ",
      as.character(
        gsea_metadata$ranking_variable
      )
    )
  )
}

if (
  !is.null(gsea_metadata$ranking_fallback) &&
  !identical(
    as.character(
      gsea_metadata$ranking_fallback
    ),
    "none"
  )
) {
  fail_with_qc(
    paste0(
      "GSEA metadata ranking_fallback is not none: ",
      as.character(
        gsea_metadata$ranking_fallback
      )
    )
  )
}

if (
  !is.null(gsea_metadata$model) &&
  !identical(
    as.character(
      gsea_metadata$model
    ),
    "expression ~ Score_z"
  )
) {
  fail_with_qc(
    paste0(
      "GSEA metadata model is not expression ~ Score_z: ",
      as.character(
        gsea_metadata$model
      )
    )
  )
}

message("[2/5] Validating input columns")

required_limma <- c(
  "gene_name",
  "t_stat"
)

missing_limma <- setdiff(
  required_limma,
  names(limma_df)
)

if (length(missing_limma) > 0) {
  fail_with_qc(
    paste0(
      "LIMMA input missing required columns: ",
      paste(
        missing_limma,
        collapse = ", "
      )
    )
  )
}

required_emt <- c(
  "pathway",
  "NES",
  "padj"
)

missing_emt <- setdiff(
  required_emt,
  names(gsea_emt)
)

if (length(missing_emt) > 0) {
  fail_with_qc(
    paste0(
      "EMT row input missing required columns: ",
      paste(
        missing_emt,
        collapse = ", "
      )
    )
  )
}

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

gsea_emt$pathway <- as.character(
  gsea_emt$pathway
)

gsea_emt$NES <- suppressWarnings(
  as.numeric(
    gsea_emt$NES
  )
)

gsea_emt$padj <- suppressWarnings(
  as.numeric(
    gsea_emt$padj
  )
)

emt_pathway <- paste0(
  "HALLMARK_EPITHELIAL_",
  "MESENCHYMAL_TRANSITION"
)

emt_row <- gsea_emt[
  gsea_emt$pathway == emt_pathway,
  ,
  drop = FALSE
]

if (nrow(emt_row) != 1) {
  fail_with_qc(
    paste0(
      "EMT row input must contain exactly one ",
      emt_pathway,
      " row. Observed: ",
      nrow(emt_row)
    )
  )
}

if (
  !is.finite(
    emt_row$NES[1]
  )
) {
  fail_with_qc(
    "EMT NES is non-finite."
  )
}

if (
  !is.finite(
    emt_row$padj[1]
  )
) {
  fail_with_qc(
    "EMT adjusted P value is non-finite."
  )
}

if (
  all(
    !is.finite(
      limma_df$t_stat
    )
  )
) {
  fail_with_qc(
    "All t_stat values are non-finite."
  )
}

if (
  any(
    is.na(limma_df$gene_name) |
    limma_df$gene_name == ""
  )
) {
  fail_with_qc(
    "LIMMA input contains missing or empty gene_name values."
  )
}

if (anyDuplicated(limma_df$gene_name)) {
  duplicated_genes <- unique(
    limma_df$gene_name[
      duplicated(
        limma_df$gene_name
      )
    ]
  )

  fail_with_qc(
    paste0(
      "LIMMA input contains duplicated gene_name values: ",
      paste(
        utils::head(
          duplicated_genes,
          25
        ),
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

style <- detect_style()

if (length(style$fallbacks) > 0) {
  for (w in style$fallbacks) {
    add_warning(w)
  }
}

for (
  colour_name in c(
    "col_ns_light",
    "col_neg",
    "col_pos",
    "col_emt"
  )
) {
  colour_value <- style[[colour_name]]

  if (!valid_colour(colour_value)) {
    fail_with_qc(
      paste0(
        "Invalid colour detected in style$",
        colour_name,
        ": ",
        colour_value
      )
    )
  }
}

emt_direction <- if (
  emt_row$NES[1] > 0
) {
  "Higher Score_z"
} else if (
  emt_row$NES[1] < 0
) {
  "Lower Score_z"
} else {
  "Zero NES"
}

###############################################################################
# 8. Build and write EMT figure
###############################################################################

message("[3/5] Building continuous EMT GSEA figure")

emt_plot <- build_emt_gsea_plot(
  limma_df,
  emt_row,
  style
)

message("[4/5] Writing continuous EMT GSEA figure")

# Exact locked output dimensions: 4.2 x 3.9 cm.
save_pdf(
  emt_plot$plot,
  output_emt,
  4.2,
  3.9,
  family = style$base_family
)

if (!file.exists(output_emt)) {
  fail_with_qc(
    paste0(
      "EMT GSEA PDF was not created: ",
      output_emt
    )
  )
}

###############################################################################
# 9. Final QC
###############################################################################

qc_status <- if (
  length(warnings_vec) > 0
) {
  "WARNING"
} else {
  "PASS"
}

qc <- list(
  inputs_found = list(
    limma_no_signature =
      file.exists(input_limma),

    gsea_emt =
      file.exists(input_gsea_emt),

    gsea_metadata =
      file.exists(input_gsea_metadata)
  ),

  ranking = list(
    n_ranked_genes =
      emt_plot$n_ranked_genes,

    n_emt_genes_in_ranking =
      emt_plot$n_emt_genes_in_ranking
  ),

  emt = list(
    NES =
      as.numeric(
        emt_row$NES[1]
      ),

    padj =
      as.numeric(
        emt_row$padj[1]
      ),

    direction =
      emt_direction
  ),

  style = list(
    reference_script =
      style$reference_script,

    reference_type =
      style$reference_type,

    reference_detected =
      style$reference_detected,

    base_family =
      style$base_family,

    col_emt =
      style$col_emt,

    col_neg =
      style$col_neg,

    col_pos =
      style$col_pos,

    emt_width = 4.2,

    emt_height = 3.9
  ),

  gsea_metadata =
    gsea_metadata,

  figures_written =
    figures_written,

  warnings =
    warnings_vec,

  status =
    qc_status
)

write_qc_report(
  qc
)

message(
  "[5/5] Final status: ",
  qc_status
)
