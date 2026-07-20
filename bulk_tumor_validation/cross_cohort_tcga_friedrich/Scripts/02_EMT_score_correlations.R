#!/usr/bin/env Rscript

###############################################################################
# Friedrich/GSE134051 + TCGA-PRAD | Platelet score vs Hallmark EMT ssGSEA
#
# Generate two cohort-specific correlation figures from the approved public
# sample-score contracts. EMT activity and platelet scores were computed and
# aligned upstream before contract materialization.
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

###############################################################################
# 1. Portable paths
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

input_tcga <- file.path(
  repo_root,
  "bulk_tumor_validation", "tcga_prad", "Contracts", "sample_scores.tsv"
)
input_friedrich <- file.path(
  repo_root,
  "bulk_tumor_validation", "friedrich_gse134051", "Contracts",
  "sample_scores.tsv"
)

results_dir <- file.path(
  repo_root,
  "bulk_tumor_validation", "cross_cohort_tcga_friedrich",
  "Results", "generated", "EMT_correlations"
)
figures_dir <- file.path(
  repo_root,
  "bulk_tumor_validation", "cross_cohort_tcga_friedrich",
  "Figures", "generated", "EMT_correlations"
)

output_pairs <- file.path(
  results_dir,
  "Platelet_score_vs_Hallmark_EMT_ssGSEA_TCGA_Friedrich_standardized_pairs.csv"
)
output_stats <- file.path(
  results_dir,
  "Platelet_score_vs_Hallmark_EMT_ssGSEA_TCGA_Friedrich_stats.csv"
)
output_qc <- file.path(
  results_dir,
  "Platelet_score_vs_Hallmark_EMT_ssGSEA_TCGA_Friedrich_QC.txt"
)
output_pdf_tcga <- file.path(
  figures_dir,
  "Platelet_score_vs_Hallmark_EMT_ssGSEA_TCGA_PRAD.pdf"
)
output_pdf_friedrich <- file.path(
  figures_dir,
  "Platelet_score_vs_Hallmark_EMT_ssGSEA_Friedrich_GSE134051.pdf"
)

###############################################################################
# 2. Contract validation and helpers
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

strict_numeric <- function(x, label) {
  if (length(x) == 0L || anyNA(x) || any(trimws(as.character(x)) == "")) {
    fail_clear(label, " contains missing or empty values.")
  }
  out <- suppressWarnings(as.numeric(x))
  if (anyNA(out) || any(!is.finite(out))) {
    fail_clear(label, " must contain finite numeric values only.")
  }
  out
}

read_contract <- function(path, cohort_id, expected_n) {
  if (!file.exists(path)) {
    fail_clear("Missing sample-score contract for ", cohort_id, ".")
  }
  tab <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expected_names <- c("cohort_id", "score_z", "emt_ssgsea")
  if (!identical(names(tab), expected_names)) {
    fail_clear(cohort_id, " sample-score contract has an unexpected schema.")
  }
  if (
    nrow(tab) != expected_n ||
    anyNA(tab$cohort_id) ||
    !identical(unique(as.character(tab$cohort_id)), cohort_id)
  ) {
    fail_clear(cohort_id, " sample-score contract has an unexpected cohort or row count.")
  }

  tab$score_z <- strict_numeric(tab$score_z, paste0(cohort_id, " score_z"))
  tab$emt_ssgsea <- strict_numeric(
    tab$emt_ssgsea,
    paste0(cohort_id, " emt_ssgsea")
  )

  if (
    abs(mean(tab$score_z)) > 1e-8 ||
    abs(stats::sd(tab$score_z) - 1) > 1e-8
  ) {
    fail_clear(cohort_id, " score_z is not standardized within tolerance.")
  }
  if (
    !is.finite(stats::sd(tab$emt_ssgsea)) ||
    stats::sd(tab$emt_ssgsea) <= 0
  ) {
    fail_clear(cohort_id, " emt_ssgsea must have positive finite variation.")
  }

  forbidden <- grepl(
    "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}|GSM[0-9]+|RIB[[:alnum:]_-]*",
    as.character(unlist(tab, use.names = FALSE)),
    ignore.case = TRUE
  )
  if (any(forbidden)) {
    fail_clear(cohort_id, " contract contains an individual technical identifier.")
  }

  tab
}

zscore <- function(x) {
  x <- as.numeric(x)
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(x_sd) || x_sd == 0) {
    fail_clear("Cannot z-score a vector with zero or non-finite standard deviation.")
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / x_sd)
}

format_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return(formatC(p, format = "e", digits = 2))
  }
  sprintf("%.3f", p)
}

format_rho <- function(x) {
  sprintf("%.2f", x)
}

###############################################################################
# 3. Correlation statistics
###############################################################################

compute_stats <- function(tab) {
  cohorts <- unique(as.character(tab$cohort))
  stats_list <- lapply(
    cohorts,
    function(coh) {
      sub <- tab[as.character(tab$cohort) == coh, , drop = FALSE]
      ok <- is.finite(sub$score_z) & is.finite(sub$emt_ssgsea)
      sub <- sub[ok, , drop = FALSE]
      if (nrow(sub) < 3L) {
        fail_clear(coh, " has fewer than 3 observations for Spearman correlation.")
      }
      ct <- suppressWarnings(
        stats::cor.test(
          sub$score_z,
          sub$emt_ssgsea,
          method = "spearman",
          exact = FALSE
        )
      )
      data.frame(
        cohort = coh,
        n = nrow(sub),
        spearman_rho = unname(ct$estimate),
        p_value = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
  do.call(rbind, stats_list)
}

###############################################################################
# 4. Historical plotting
###############################################################################

make_plot <- function(
    tab,
    stats_tab,
    cohort_subset,
    display_title = NULL
) {
  tab <- tab[as.character(tab$cohort) %in% cohort_subset, , drop = FALSE]
  stats_tab <- stats_tab[
    as.character(stats_tab$cohort) %in% cohort_subset,
    ,
    drop = FALSE
  ]

  if (nrow(tab) == 0L || nrow(stats_tab) == 0L) {
    fail_clear("No data available for requested plot cohort subset.")
  }

  cohort_title <- if (is.null(display_title)) {
    as.character(stats_tab$cohort[1L])
  } else {
    display_title
  }

  xr <- range(tab$score_z, na.rm = TRUE)
  yr <- range(tab$emt_ssgsea_z, na.rm = TRUE)
  annotation_x <- xr[1L] + 0.010 * diff(xr)
  annotation_y_rho <- yr[2L] - 0.010 * diff(yr)
  annotation_y_p <- annotation_y_rho - 0.085 * diff(yr)

  annotation_df <- data.frame(
    x = annotation_x,
    rho_label = paste0(
      "Spearman rho = ",
      format_rho(stats_tab$spearman_rho[1L])
    ),
    p_label = paste0(
      "italic(P) == '",
      format_p(stats_tab$p_value[1L]),
      "'"
    ),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(
    tab,
    ggplot2::aes(
      x = score_z,
      y = emt_ssgsea_z
    )
  ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      fill = "#B5E18B",
      linewidth = 0.22,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      color = "#0F2854",
      size = 0.5,
      alpha = 0.70,
      shape = 16
    ) +
    ggplot2::geom_text(
      data = annotation_df,
      ggplot2::aes(
        x = x,
        y = Inf,
        label = rho_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 2.35,
      family = "Helvetica",
      size = 4.5 / ggplot2::.pt
    ) +
    ggplot2::geom_text(
      data = annotation_df,
      ggplot2::aes(
        x = x,
        y = Inf,
        label = p_label
      ),
      inherit.aes = FALSE,
      parse = TRUE,
      hjust = 0,
      vjust = 3.95,
      family = "Helvetica",
      size = 4.5 / ggplot2::.pt
    ) +
    ggplot2::labs(
      title = cohort_title,
      x = "Platelet score (z-score)",
      y = "Hallmark EMT score (z-score)"
    ) +
    ggplot2::theme_classic(
      base_size = 5,
      base_family = "Helvetica"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        size = 5
      ),
      axis.title = ggplot2::element_text(
        size = 4.5
      ),
      axis.text = ggplot2::element_text(
        size = 4,
        color = "black"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.20,
        color = "black"
      ),
      axis.ticks.length = grid::unit(
        0.55,
        "mm"
      ),
      plot.title = ggplot2::element_text(
        size = 5,
        face = "plain",
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 1.0
        )
      ),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.20
      ),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      plot.margin = ggplot2::margin(
        1.0,
        1.5,
        1.5,
        1.5,
        "mm"
      )
    )

  p
}

###############################################################################
# 5. Main
###############################################################################

main <- function() {
  require_packages(c("ggplot2", "grid"))

  tcga <- read_contract(input_tcga, COHORT_TCGA_ID, 497L)
  friedrich <- read_contract(input_friedrich, COHORT_FRIEDRICH_ID, 164L)

  tcga$cohort <- COHORT_TCGA_LABEL
  friedrich$cohort <- COHORT_FRIEDRICH_LABEL
  tcga$emt_ssgsea_z <- zscore(tcga$emt_ssgsea)
  friedrich$emt_ssgsea_z <- zscore(friedrich$emt_ssgsea)

  out_tab <- rbind(tcga, friedrich)
  out_tab$cohort <- factor(
    as.character(out_tab$cohort),
    levels = c(COHORT_TCGA_LABEL, COHORT_FRIEDRICH_LABEL)
  )
  if (anyNA(out_tab$cohort)) {
    fail_clear("Combined cohort field contains an unexpected label.")
  }

  stats_tab <- compute_stats(out_tab)
  pair_out <- out_tab[, c(
    "cohort_id",
    "score_z",
    "emt_ssgsea",
    "emt_ssgsea_z"
  )]

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

  utils::write.csv(pair_out, output_pairs, row.names = FALSE)
  utils::write.csv(stats_tab, output_stats, row.names = FALSE)

  p_tcga <- make_plot(
    out_tab,
    stats_tab,
    cohort_subset = COHORT_TCGA_LABEL
  )
  ggplot2::ggsave(
    filename = output_pdf_tcga,
    plot = p_tcga,
    width = 4.2,
    height = 4.2,
    units = "cm",
    device = grDevices::pdf,
    family = "Helvetica",
    useDingbats = FALSE,
    bg = "transparent"
  )

  p_friedrich <- make_plot(
    out_tab,
    stats_tab,
    cohort_subset = COHORT_FRIEDRICH_LABEL,
    display_title = COHORT_FRIEDRICH_DISPLAY_LABEL
  )
  ggplot2::ggsave(
    filename = output_pdf_friedrich,
    plot = p_friedrich,
    width = 4.2,
    height = 4.2,
    units = "cm",
    device = grDevices::pdf,
    family = "Helvetica",
    useDingbats = FALSE,
    bg = "transparent"
  )

  qc_lines <- c(
    "Platelet score versus Hallmark EMT contract-consumer QC",
    "=======================================================",
    "",
    "Inputs:",
    "  tcga_prad/Contracts/sample_scores.tsv",
    "  friedrich_gse134051/Contracts/sample_scores.tsv",
    "",
    paste0("TCGA paired observations: ", nrow(tcga)),
    paste0("Friedrich paired observations: ", nrow(friedrich)),
    "EMT z-score: derived within cohort using sample mean and sample SD",
    "Correlation: Spearman; exact = FALSE",
    "Visible labels: TCGA-PRAD; Friedrich (GSE134051)",
    "",
    paste0("Pair table: ", basename(output_pairs)),
    paste0("Statistics table: ", basename(output_stats)),
    paste0("TCGA PDF: ", basename(output_pdf_tcga)),
    paste0("Friedrich PDF: ", basename(output_pdf_friedrich)),
    "",
    "Final status: PASS"
  )
  writeLines(qc_lines, output_qc, useBytes = TRUE)

  msg("Validated observations: TCGA=497; Friedrich=164")
  msg("Final status: PASS")
}

main()
