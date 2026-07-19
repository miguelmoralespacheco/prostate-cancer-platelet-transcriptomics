#!/usr/bin/env Rscript
###############################################################################
# Friedrich/GSE134051 clinical final analysis using official 41-gene platelet score.
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

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(broom)
  library(grid)
  library(scales)
})

###############################################################################
# 0) EDITORIAL STYLE
#
# LOCKED: no graphical values in this section should be modified.
###############################################################################

font_use <- "Helvetica"
font_family <- font_use

pt_axis <- 5.0
pt_title <- 5.0
pt_small <- 4.0

lw_axis <- 0.30
lw_tick <- 0.30
lw_line <- 0.45
lw_box <- 0.35
lw_border <- 0.35

pt_dot <- 0.65
alpha_dot <- 0.40

col_low <- "#0F558F"
col_high <- "#E47513"
col_neutral <- "grey85"
col_black <- "black"

pt_to_annotate <- function(pt) pt / 2.8

theme_nature <- function() {
  theme_classic(base_family = font_use, base_size = pt_axis) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = 3.5, color = "black"),
      axis.title = element_text(size = pt_title, color = "black"),
      plot.title = element_text(
        size = pt_title,
        face = "plain",
        hjust = 0.5,
        color = "black"
      ),
      plot.subtitle = element_text(size = pt_small, color = "black"),
      axis.line = element_line(linewidth = lw_axis, color = "black"),
      axis.ticks = element_line(linewidth = lw_tick, color = "black"),
      axis.ticks.length = unit(0.7, "mm"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pt_small, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = pt_axis, color = "black"),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2, "mm")
    )
}

theme_nature_box <- function() {
  theme_nature() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = lw_border
      ),
      axis.line = element_blank()
    )
}

###############################################################################
# 1) PATHS
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

script_fp <- file.path(
  project_dir,
  "Scripts",
  "Clinical",
  "01_Friedrich_GSE134051_clinical_final.R"
)

score_fp <- file.path(
  GENERATED_RESULTS_DIR, "Score",
  "Friedrich_GSE134051_primary_master_table_41genes.csv"
)

clinical_fp <- file.path(
  INPUT_DIR, "LocalLarge",
  "Friedrich_GSE134051_clinical.rds"
)

results_dir <- file.path(
  GENERATED_RESULTS_DIR, "Clinical"
)

tab_dir <- file.path(
  results_dir,
  "Tables"
)

pdf_dir <- file.path(
  GENERATED_FIGURES_DIR, "Clinical"
)

rds_dir <- file.path(
  results_dir,
  "Objects"
)

report_dir <- file.path(
  results_dir,
  "Logs"
)

invisible(
  lapply(
    c(
      results_dir,
      tab_dir,
      pdf_dir,
      rds_dir,
      report_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

###############################################################################
# 2) HELPERS
###############################################################################

warnings_vec <- character()
tables_generated <- character()
figures_generated <- character()
rds_generated <- character()

add_warning <- function(x) {
  message(
    "WARNING: ",
    x
  )

  warnings_vec <<- unique(
    c(
      warnings_vec,
      x
    )
  )

  invisible(x)
}

write_table <- function(x, fp) {
  utils::write.csv(
    as.data.frame(x),
    fp,
    row.names = FALSE
  )

  tables_generated <<- unique(
    c(
      tables_generated,
      fp
    )
  )

  invisible(fp)
}

write_rds_tracked <- function(x, fp) {
  saveRDS(
    x,
    fp
  )

  rds_generated <<- unique(
    c(
      rds_generated,
      fp
    )
  )

  invisible(fp)
}

save_pdf <- function(
    plot,
    fp,
    width_cm,
    height_cm
) {
  ggplot2::ggsave(
    filename = fp,
    plot = plot,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = grDevices::pdf,
    family = font_use,
    useDingbats = FALSE
  )

  figures_generated <<- unique(
    c(
      figures_generated,
      fp
    )
  )

  invisible(fp)
}

save_surv_pdf <- function(
    gp,
    fp,
    width_cm,
    height_cm
) {
  ggplot2::ggsave(
    filename = fp,
    plot = gp,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = grDevices::pdf,
    family = font_use,
    useDingbats = FALSE
  )

  figures_generated <<- unique(
    c(
      figures_generated,
      fp
    )
  )

  invisible(fp)
}

as_num_safe <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 1e-4,
      format(
        p,
        scientific = TRUE,
        digits = 2
      ),
      sprintf(
        "%.4f",
        p
      )
    )
  )
}

fmt_p_fig <- function(p) {
  if (is.na(p)) {
    return("p = NA")
  }

  if (p < 1e-3) {
    return(
      paste0(
        "p = ",
        format(
          p,
          scientific = TRUE,
          digits = 2
        )
      )
    )
  }

  paste0(
    "p = ",
    sprintf(
      "%.3f",
      p
    )
  )
}

plotmath_p_label <- function(p) {
  if (is.na(p)) {
    return("italic(p) == 'NA'")
  }

  p_txt <- if (p < 0.001) {
    formatC(
      p,
      format = "e",
      digits = 2
    )
  } else {
    formatC(
      p,
      format = "f",
      digits = 4
    )
  }

  sprintf(
    "italic(p) == '%s'",
    p_txt
  )
}

scale_numeric <- function(x) {
  x <- as_num_safe(x)

  if (
    sum(
      is.finite(x) &
      !is.na(x)
    ) < 2
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  as.numeric(
    scale(x)
  )
}

count_nonmissing <- function(x) {
  sum(
    !is.na(x) &
      as.character(x) != "",
    na.rm = TRUE
  )
}

count_finite <- function(x) {
  sum(
    is.finite(
      as_num_safe(x)
    ),
    na.rm = TRUE
  )
}

safe_block <- function(label, expr) {
  expr <- substitute(expr)

  tryCatch(
    {
      eval.parent(expr)
      invisible(TRUE)
    },
    error = function(e) {
      add_warning(
        paste0(
          label,
          " failed: ",
          conditionMessage(e)
        )
      )

      invisible(FALSE)
    }
  )
}

require_input <- function(fp, label) {
  if (!file.exists(fp)) {
    stop(
      label,
      " was not found: ",
      fp,
      call. = FALSE
    )
  }
}

pick_first <- function(nms, candidates) {
  hit <- candidates[
    candidates %in% nms
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

validate_friedrich_clinical_identity <- function(clinical_df) {
  forbidden_id_pattern <- paste0("ICGC", "_PCA")
  required_identity_columns <- c(
    "sample_name",
    "sample_type",
    "study_name"
  )
  missing_identity_columns <- setdiff(
    required_identity_columns,
    names(clinical_df)
  )

  if (length(missing_identity_columns) > 0) {
    stop(
      paste0(
        "Possible incorrect clinical input: missing Friedrich identity columns: ",
        paste(missing_identity_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (nrow(clinical_df) != EXPECTED_TOTAL_N) {
    stop(
      paste0(
        "Possible incorrect input or cohort contamination: expected ",
        EXPECTED_TOTAL_N,
        " clinical records for ",
        COHORT_LABEL,
        "; observed ",
        nrow(clinical_df),
        "."
      ),
      call. = FALSE
    )
  }

  sample_ids <- trimws(as.character(clinical_df$sample_name))
  if (
    any(!grepl("^GSM[0-9]+$", sample_ids)) ||
    any(grepl(forbidden_id_pattern, sample_ids, ignore.case = TRUE)) ||
    anyDuplicated(sample_ids)
  ) {
    stop(
      paste0(
        "Possible incorrect input or cohort contamination: ",
        COHORT_LABEL,
        " clinical sample_name values must be unique GSM identifiers."
      ),
      call. = FALSE
    )
  }

  if ("alt_sample_name" %in% names(clinical_df)) {
    alt_ids <- trimws(as.character(clinical_df$alt_sample_name))
    if (
      any(is.na(alt_ids) | !grepl("^RIB", alt_ids)) ||
      any(grepl(forbidden_id_pattern, alt_ids, ignore.case = TRUE))
    ) {
      stop(
        paste0(
          "Possible incorrect input or cohort contamination: ",
          COHORT_LABEL,
          " alt_sample_name values must be RIB identifiers."
        ),
        call. = FALSE
      )
    }
  }

  study_name <- trimws(as.character(clinical_df$study_name))
  if (
    any(is.na(study_name) | study_name == "") ||
    any(!grepl("Friedrich", study_name, ignore.case = TRUE))
  ) {
    stop(
      paste0(
        "Possible incorrect input or cohort contamination: study_name is not compatible with ",
        SOURCE_STUDY,
        "."
      ),
      call. = FALSE
    )
  }

  sample_type <- trimws(as.character(clinical_df$sample_type))
  observed_counts <- table(
    factor(sample_type, levels = c("primary", "normal", "BPH"))
  )
  expected_counts <- c(primary = 164L, normal = 52L, BPH = 39L)

  if (
    any(!sample_type %in% names(expected_counts)) ||
    !identical(as.integer(observed_counts), as.integer(expected_counts))
  ) {
    stop(
      paste0(
        "Possible incorrect input or cohort contamination: expected sample_type distribution ",
        "primary=164, normal=52, BPH=39 for ",
        COHORT_LABEL,
        "."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

require_input(
  score_fp,
  "Official 41-gene platelet score"
)

require_input(
  clinical_fp,
  "Friedrich/GSE134051 clinical table"
)

###############################################################################
# 3) LOAD OFFICIAL Friedrich/GSE134051 41-GENE PLATELET SCORE
###############################################################################

score_raw <- read.csv(
  score_fp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_col <- pick_first(
  names(score_raw),
  c(
    "sample_name",
    "sample_id",
    "Sample",
    "sample",
    "barcode"
  )
)

score_z_col <- pick_first(
  names(score_raw),
  c(
    "Score_z",
    "score_z",
    "platelet_score_z",
    "PlateletScore_z"
  )
)

score_raw_col <- pick_first(
  names(score_raw),
  c(
    "Score_raw",
    "score_raw",
    "platelet_score_raw",
    "PlateletScore_raw"
  )
)

if (is.na(sample_col)) {
  stop(
    "No sample column found in Friedrich/GSE134051 score file.",
    call. = FALSE
  )
}

if (is.na(score_z_col)) {
  stop(
    "No Score_z column found in Friedrich/GSE134051 score file.",
    call. = FALSE
  )
}

if (is.na(score_raw_col)) {
  add_warning(
    paste0(
      "No Score_raw column found in Friedrich/GSE134051 score file; ",
      "Score_raw will be NA."
    )
  )
}

score_df <- score_raw %>%
  transmute(
    sample_name = trimws(
      as.character(
        .data[[sample_col]]
      )
    ),

    Score_raw = if (!is.na(score_raw_col)) {
      as_num_safe(
        .data[[score_raw_col]]
      )
    } else {
      NA_real_
    },

    Score_z = as_num_safe(
      .data[[score_z_col]]
    )
  )

forbidden_id_pattern <- paste0("^ICGC", "_PCA")
if (
  any(!grepl("^GSM[0-9]+$", score_df$sample_name)) ||
  any(grepl(forbidden_id_pattern, score_df$sample_name, ignore.case = TRUE))
) {
  stop(
    paste0(
      "Possible incorrect input or cohort contamination: score sample_name ",
      "values must be GSM identifiers from ",
      COHORT_LABEL,
      "."
    ),
    call. = FALSE
  )
}

if (
  any(
    is.na(score_df$sample_name) |
    score_df$sample_name == ""
  )
) {
  stop(
    "Friedrich/GSE134051 score table contains missing or empty sample_name values.",
    call. = FALSE
  )
}

if (anyDuplicated(score_df$sample_name)) {
  stop(
    "Friedrich/GSE134051 score table contains duplicated sample_name values.",
    call. = FALSE
  )
}

if (!all(is.finite(score_df$Score_z))) {
  stop(
    "Score_z contains non-finite values.",
    call. = FALSE
  )
}

if (nrow(score_df) != EXPECTED_PRIMARY_N) {
  stop(
    paste0(
      "Official score table must contain exactly ",
      EXPECTED_PRIMARY_N,
      " primary samples. Observed: ",
      nrow(score_df),
      "."
    ),
    call. = FALSE
  )
}

###############################################################################
# 4) LOAD Friedrich/GSE134051 CLINICAL TABLE
###############################################################################

clin_raw <- readRDS(
  clinical_fp
)

clin <- as.data.frame(
  clin_raw,
  stringsAsFactors = FALSE
)

validate_friedrich_clinical_identity(clin)

required_clinical_columns <- c(
  "sample_name",
  "sample_type",
  "days_to_overall_survival",
  "overall_survival_status",
  "gleason_grade",
  "days_to_disease_specific_recurrence",
  "disease_specific_recurrence_status",
  "days_to_metastatic_occurrence",
  "metastasis_occurrence_status"
)

missing_clinical_columns <- setdiff(
  required_clinical_columns,
  names(clin)
)

if (length(missing_clinical_columns) > 0) {
  stop(
    paste0(
      "Friedrich/GSE134051 clinical table is missing required columns: ",
      paste(
        missing_clinical_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

clin_primary <- clin %>%
  mutate(
    sample_name = trimws(
      as.character(sample_name)
    ),
    sample_type = as.character(sample_type)
  ) %>%
  filter(
    sample_type == "primary"
  )

if (nrow(clin_primary) == 0) {
  stop(
    "No primary samples found in Friedrich/GSE134051 clinical table.",
    call. = FALSE
  )
}

if (
  any(
    is.na(clin_primary$sample_name) |
    clin_primary$sample_name == ""
  )
) {
  stop(
    "Friedrich/GSE134051 primary clinical table contains missing sample_name values.",
    call. = FALSE
  )
}

if (anyDuplicated(clin_primary$sample_name)) {
  stop(
    paste0(
      "Friedrich/GSE134051 primary clinical table contains duplicated ",
      "sample_name values."
    ),
    call. = FALSE
  )
}

if (nrow(clin_primary) != EXPECTED_PRIMARY_N) {
  stop(
    paste0(
      "Friedrich/GSE134051 primary clinical table must contain exactly ",
      EXPECTED_PRIMARY_N,
      " samples. Observed: ",
      nrow(clin_primary),
      "."
    ),
    call. = FALSE
  )
}

score_only_samples <- setdiff(
  score_df$sample_name,
  clin_primary$sample_name
)

clinical_only_samples <- setdiff(
  clin_primary$sample_name,
  score_df$sample_name
)

if (
  length(score_only_samples) > 0 ||
  length(clinical_only_samples) > 0
) {
  stop(
    paste0(
      "Score and clinical primary sample sets do not match. ",
      "Score-only samples: ",
      if (length(score_only_samples) == 0) {
        "none"
      } else {
        paste(
          score_only_samples,
          collapse = ", "
        )
      },
      ". Clinical-only samples: ",
      if (length(clinical_only_samples) == 0) {
        "none"
      } else {
        paste(
          clinical_only_samples,
          collapse = ", "
        )
      },
      "."
    ),
    call. = FALSE
  )
}

###############################################################################
# 5) MASTER CLINICAL TABLE
###############################################################################

df_clinical_final <- score_df %>%
  inner_join(
    clin_primary,
    by = "sample_name"
  ) %>%
  mutate(
    OS_days = as_num_safe(
      days_to_overall_survival
    ),

    OS_months = OS_days / 30.4375,

    OS_event = as.integer(
      as_num_safe(
        overall_survival_status
      ) > 0
    ),

    gleason_grade_numeric = as_num_safe(
      gleason_grade
    ),

    gleason_score_group = dplyr::case_when(
      is.finite(gleason_grade_numeric) &
        gleason_grade_numeric <= 6 ~ "<=6",

      is.finite(gleason_grade_numeric) &
        gleason_grade_numeric == 7 ~ "7",

      is.finite(gleason_grade_numeric) &
        gleason_grade_numeric == 8 ~ "8",

      is.finite(gleason_grade_numeric) &
        gleason_grade_numeric >= 9 ~ ">=9",

      TRUE ~ NA_character_
    ),

    gleason_score_group = factor(
      gleason_score_group,
      levels = c(
        "<=6",
        "7",
        "8",
        ">=9"
      ),
      ordered = TRUE
    ),

    grade_group_clean = as.character(
      grade_group
    ),

    grade_group_clean = factor(
      grade_group_clean,
      levels = c(
        "<=6",
        "7",
        ">=8"
      ),
      ordered = TRUE
    ),

    platelet_score_median_group = ifelse(
      Score_z <= stats::median(
        Score_z,
        na.rm = TRUE
      ),
      "Low",
      "High"
    ),

    platelet_score_median_group = factor(
      platelet_score_median_group,
      levels = c(
        "Low",
        "High"
      )
    )
  )

if (nrow(df_clinical_final) != EXPECTED_PRIMARY_N) {
  stop(
    paste0(
      "Final Friedrich/GSE134051 clinical table must contain exactly ",
      EXPECTED_PRIMARY_N,
      " samples. Observed after join: ",
      nrow(df_clinical_final),
      "."
    ),
    call. = FALSE
  )
}

if (anyDuplicated(df_clinical_final$sample_name)) {
  stop(
    "Final Friedrich/GSE134051 clinical table contains duplicated sample_name values.",
    call. = FALSE
  )
}

if (!identical(df_clinical_final$sample_name, score_df$sample_name)) {
  stop(
    paste0(
      "Possible incorrect input or cohort contamination: final clinical and score ",
      "samples are not aligned in the same order."
    ),
    call. = FALSE
  )
}

write_table(
  df_clinical_final,
  file.path(
    tab_dir,
    "Friedrich_GSE134051_clinical_final_master_table.csv"
  )
)

write_rds_tracked(
  df_clinical_final,
  file.path(
    rds_dir,
    "Friedrich_GSE134051_clinical_final_master_table.rds"
  )
)

###############################################################################
# 6) CLINICAL COMPLETENESS AND ENDPOINT COUNTS
###############################################################################

clinical_completeness <- tibble(
  variable = c(
    "age_at_initial_diagnosis",
    "gleason_grade",
    "gleason_major",
    "gleason_minor",
    "grade_group",
    "T_pathological",
    "T_clinical",
    "N_stage",
    "M_stage",
    "psa",
    "extraprostatic_extension",
    "perineural_invasion",
    "seminal_vesicle_invasion",
    "angiolymphatic_invasion",
    "tumor_margins_positive",
    "AR_activity",
    "prolaris",
    "decipher",
    "oncotypedx",
    "genome_altered"
  )
) %>%
  mutate(
    n_nonmissing = vapply(
      variable,
      function(v) {
        if (!v %in% names(clin_primary)) {
          return(NA_integer_)
        }

        count_nonmissing(
          clin_primary[[v]]
        )
      },
      integer(1)
    )
  )

write_table(
  clinical_completeness,
  file.path(
    tab_dir,
    "Clinical_variable_completeness.csv"
  )
)

endpoint_event_counts <- tibble(
  endpoint = c(
    "OS",
    "Disease-specific recurrence",
    "Metastatic occurrence"
  ),

  time_col = c(
    "days_to_overall_survival",
    "days_to_disease_specific_recurrence",
    "days_to_metastatic_occurrence"
  ),

  status_col = c(
    "overall_survival_status",
    "disease_specific_recurrence_status",
    "metastasis_occurrence_status"
  ),

  n_with_time = c(
    count_finite(
      clin_primary$days_to_overall_survival
    ),

    count_finite(
      clin_primary$days_to_disease_specific_recurrence
    ),

    count_finite(
      clin_primary$days_to_metastatic_occurrence
    )
  ),

  events = c(
    sum(
      as_num_safe(
        clin_primary$overall_survival_status
      ) == 1,
      na.rm = TRUE
    ),

    sum(
      as_num_safe(
        clin_primary$disease_specific_recurrence_status
      ) == 1,
      na.rm = TRUE
    ),

    sum(
      as_num_safe(
        clin_primary$metastasis_occurrence_status
      ) == 1,
      na.rm = TRUE
    )
  )
)

write_table(
  endpoint_event_counts,
  file.path(
    tab_dir,
    "Endpoint_event_counts.csv"
  )
)

if (
  endpoint_event_counts$events[
    endpoint_event_counts$endpoint == "OS"
  ] < 30
) {
  add_warning(
    paste0(
      "OS has fewer than 30 events; ",
      "all OS analyses are exploratory and unstable."
    )
  )
}

if (
  endpoint_event_counts$n_with_time[
    endpoint_event_counts$endpoint ==
    "Disease-specific recurrence"
  ] == 0
) {
  add_warning(
    paste0(
      "Disease-specific recurrence is unavailable ",
      "in Friedrich/GSE134051 primary clinical data."
    )
  )
}

if (
  endpoint_event_counts$n_with_time[
    endpoint_event_counts$endpoint ==
    "Metastatic occurrence"
  ] == 0
) {
  add_warning(
    paste0(
      "Metastatic occurrence is unavailable ",
      "in Friedrich/GSE134051 primary clinical data."
    )
  )
}

###############################################################################
# 7) MINIMAL CLINICAL ASSOCIATIONS
###############################################################################

spearman_test <- function(
    df,
    variable,
    label
) {
  sub <- df %>%
    select(
      Score_z,
      value = all_of(variable)
    ) %>%
    filter(
      is.finite(Score_z),
      is.finite(value)
    )

  if (nrow(sub) < 10) {
    return(
      tibble(
        endpoint = label,
        test = "Spearman",
        n = nrow(sub),
        estimate = NA_real_,
        p_value = NA_real_
      )
    )
  }

  tt <- suppressWarnings(
    stats::cor.test(
      sub$Score_z,
      sub$value,
      method = "spearman"
    )
  )

  tibble(
    endpoint = label,
    test = "Spearman",
    n = nrow(sub),
    estimate = unname(tt$estimate),
    p_value = tt$p.value
  )
}

kw_test <- function(
    df,
    variable,
    label
) {
  sub <- df %>%
    select(
      Score_z,
      group = all_of(variable)
    ) %>%
    filter(
      is.finite(Score_z),
      !is.na(group)
    )

  if (
    nrow(sub) < 10 ||
    length(
      unique(sub$group)
    ) < 2
  ) {
    return(
      tibble(
        endpoint = label,
        test = "Kruskal-Wallis",
        n = nrow(sub),
        estimate = NA_real_,
        p_value = NA_real_
      )
    )
  }

  tt <- stats::kruskal.test(
    Score_z ~ group,
    data = sub
  )

  tibble(
    endpoint = label,
    test = "Kruskal-Wallis",
    n = nrow(sub),
    estimate = unname(tt$statistic),
    p_value = tt$p.value
  )
}

assoc_rows <- bind_rows(
  spearman_test(
    df_clinical_final,
    "gleason_grade_numeric",
    "Gleason score"
  ),

  kw_test(
    df_clinical_final,
    "gleason_score_group",
    "Grouped Gleason score"
  )
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    ),

    p_label = fmt_p(
      p_value
    ),

    FDR_label = fmt_p(
      FDR
    )
  )

write_table(
  assoc_rows,
  file.path(
    tab_dir,
    "Clinical_assoc_minimal_table.csv"
  )
)

###############################################################################
# 9) FIGURE: SCORE_Z VS GROUPED GLEASON SCORE
#
# LOCKED GRAPHICAL CONSTRUCTION.
###############################################################################

safe_block("Score_z vs grouped Gleason score", {
  df_gl <- df_clinical_final %>%
    filter(
      is.finite(Score_z),
      !is.na(gleason_score_group)
    ) %>%
    mutate(
      gleason_score_group = factor(
        as.character(gleason_score_group),
        levels = c(
          "<=6",
          "7",
          "8",
          ">=9"
        ),
        ordered = TRUE
      ),

      g_rank = as.numeric(
        gleason_score_group
      )
    )

  if (
    nrow(df_gl) < 10 ||
    length(
      unique(df_gl$gleason_score_group)
    ) < 2
  ) {
    stop(
      "insufficient grouped Gleason score levels"
    )
  }

  gleason_test <- suppressWarnings(
    stats::cor.test(
      df_gl$Score_z,
      df_gl$g_rank,
      method = "spearman",
      exact = FALSE
    )
  )

  rho <- unname(
    gleason_test$estimate
  )

  pval <- gleason_test$p.value

  df_sum <- df_gl %>%
    group_by(
      gleason_score_group,
      g_rank
    ) %>%
    summarise(
      n = n(),

      med = median(
        Score_z,
        na.rm = TRUE
      ),

      q1 = quantile(
        Score_z,
        0.25,
        na.rm = TRUE
      ),

      q3 = quantile(
        Score_z,
        0.75,
        na.rm = TRUE
      ),

      .groups = "drop"
    )

  write_table(
    df_sum,
    file.path(
      tab_dir,
      "PlateletScore_grouped_Gleason_score_summary.csv"
    )
  )

  write_table(
    tibble(
      n = nrow(df_gl),
      rho = rho,
      p_value = pval,
      p_label = fmt_p(pval)
    ),
    file.path(
      tab_dir,
      "PlateletScore_grouped_Gleason_score_Spearman.csv"
    )
  )

  gleason_x_labels <- setNames(
    paste0(
      as.character(
        df_sum$gleason_score_group
      ),
      "\n(n=",
      df_sum$n,
      ")"
    ),
    df_sum$g_rank
  )

  gleason_y_quantiles <- stats::quantile(
    df_gl$Score_z,
    c(
      0.01,
      0.99
    ),
    na.rm = TRUE
  )

  gleason_y_padding <- diff(
    gleason_y_quantiles
  ) * 0.35

  gleason_y_limits <- c(
    gleason_y_quantiles[1L] -
      gleason_y_padding,

    gleason_y_quantiles[2L] +
      gleason_y_padding
  )

  gleason_x_min <- min(
    df_sum$g_rank
  )

  gleason_x_max <- max(
    df_sum$g_rank
  )

  gleason_p_text <- if (
    pval < 0.001
  ) {
    formatC(
      pval,
      format = "e",
      digits = 2
    )
  } else {
    sprintf(
      "%.3f",
      pval
    )
  }

  gleason_p_label <- paste0(
    "italic(P) == '",
    gleason_p_text,
    "'"
  )

  p_gleason <- ggplot(
    df_gl,
    aes(
      x = g_rank,
      y = Score_z
    )
  ) +
    geom_point(
      position = position_jitter(
        width = 0.12,
        height = 0
      ),
      size = 0.45,
      alpha = 0.6,
      shape = 16,
      color = "#3E6A8E"
    ) +
    geom_linerange(
      data = df_sum,
      aes(
        x = g_rank,
        ymin = q1,
        ymax = q3
      ),
      inherit.aes = FALSE,
      linewidth = 0.15,
      color = "black"
    ) +
    geom_point(
      data = df_sum,
      aes(
        x = g_rank,
        y = med
      ),
      inherit.aes = FALSE,
      size = 0.60,
      shape = 16,
      color = "black"
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.20,
      color = "#243B6B"
    ) +
    scale_x_continuous(
      breaks = df_sum$g_rank,
      labels = gleason_x_labels,
      limits = c(
        gleason_x_min - 0.35,
        gleason_x_max + 0.35
      ),
      expand = c(
        0,
        0
      )
    ) +
    coord_cartesian(
      ylim = gleason_y_limits
    ) +
    annotate(
      "text",
      x = gleason_x_min - 0.22,
      y = gleason_y_limits[2L] -
        0.018 * diff(gleason_y_limits),
      label = paste0(
        "Spearman rho = ",
        sprintf(
          "%.3f",
          rho
        )
      ),
      hjust = 0,
      vjust = 1,
      size = 4 / 2.845,
      family = font_family,
      color = "black"
    ) +
    annotate(
      "text",
      x = gleason_x_min - 0.22,
      y = gleason_y_limits[2L] -
        0.105 * diff(gleason_y_limits),
      label = gleason_p_label,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 4 / 2.845,
      family = font_family,
      color = "black"
    ) +
    labs(
      x = "Gleason score group",
      y = "Platelet score (z)"
    ) +
    theme_bw(
      base_family = font_family,
      base_size = 5
    ) +
    theme(
      panel.grid = element_blank(),

      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.22
      ),

      axis.text.x = element_text(
        family = font_family,
        color = "black",
        size = 5,
        lineheight = 0.75,
        margin = margin(
          t = 1
        )
      ),

      axis.text.y = element_text(
        family = font_family,
        color = "black",
        size = 5
      ),

      axis.title = element_text(
        family = font_family,
        color = "black",
        size = 5
      ),

      axis.ticks = element_line(
        linewidth = 0.20,
        color = "black"
      ),

      plot.margin = margin(
        1,
        1,
        2,
        1,
        unit = "mm"
      )
    )

  save_pdf(
    p_gleason,
    file.path(
      pdf_dir,
      "PlateletScore_vs_grouped_Gleason_score.pdf"
    ),
    width_cm = 4.2,
    height_cm = 3.5
  )
})

###############################################################################
# 10) SURVIVAL HELPERS
###############################################################################

make_os_df <- function(df) {
  df %>%
    transmute(
      sample_name,

      Score_z = as_num_safe(
        Score_z
      ),

      gleason_grade_numeric = as_num_safe(
        gleason_grade_numeric
      ),

      time = as_num_safe(
        OS_months
      ),

      event = as.integer(
        as_num_safe(
          OS_event
        ) > 0
      )
    ) %>%
    filter(
      is.finite(Score_z),
      is.finite(time),
      time > 0,
      !is.na(event)
    ) %>%
    mutate(
      median_cutpoint = stats::median(
        Score_z,
        na.rm = TRUE
      ),

      score_group = factor(
        ifelse(
          Score_z <= median_cutpoint,
          "Low",
          "High"
        ),
        levels = c(
          "Low",
          "High"
        )
      ),

      Score_z_s = scale_numeric(
        Score_z
      ),

      gleason_s = scale_numeric(
        gleason_grade_numeric
      )
    )
}

logrank_p <- function(fit, df) {
  sd <- survival::survdiff(
    survival::Surv(
      time,
      event
    ) ~ score_group,
    data = df
  )

  stats::pchisq(
    sd$chisq,
    df = length(sd$n) - 1,
    lower.tail = FALSE
  )
}

tidy_cox_or_empty <- function(
    fit,
    endpoint,
    model,
    n,
    events
) {
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      endpoint = endpoint,
      model = model,
      n = n,
      events = events
    ) %>%
    select(
      endpoint,
      model,
      n,
      events,
      everything()
    )
}

###############################################################################
# 11) OS SURVIVAL
#
# LOCKED GRAPHICAL CONSTRUCTION.
###############################################################################

survival_summary_all <- tibble()

safe_block("OS survival", {
  df_os <- make_os_df(
    df_clinical_final
  )

  n_os <- nrow(
    df_os
  )

  events_os <- sum(
    df_os$event == 1L,
    na.rm = TRUE
  )

  if (
    n_os < 10 ||
    events_os < 1 ||
    length(
      unique(df_os$score_group)
    ) < 2
  ) {
    stop(
      "insufficient OS data"
    )
  }

  fit_km <- survival::survfit(
    survival::Surv(
      time,
      event
    ) ~ score_group,
    data = df_os
  )

  p_lr <- logrank_p(
    fit_km,
    df_os
  )

  km_summary <- tibble(
    endpoint = "OS",
    analysis = "KM_median_split",
    model = "logrank",
    n = n_os,
    events = events_os,
    term = "High_vs_Low",
    estimate = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    p_value = p_lr,
    p_label = fmt_p(p_lr),
    median_cutpoint = unique(
      df_os$median_cutpoint
    )[1],
    n_low = sum(
      df_os$score_group == "Low",
      na.rm = TRUE
    ),
    n_high = sum(
      df_os$score_group == "High",
      na.rm = TRUE
    )
  )

  cox_bin <- survival::coxph(
    survival::Surv(
      time,
      event
    ) ~ score_group,
    data = df_os
  )

  cox_bin_tab <- tidy_cox_or_empty(
    cox_bin,
    endpoint = "OS",
    model = "binary_median",
    n = n_os,
    events = events_os
  )

  cox_cont <- survival::coxph(
    survival::Surv(
      time,
      event
    ) ~ Score_z_s,
    data = df_os
  )

  cox_cont_tab <- tidy_cox_or_empty(
    cox_cont,
    endpoint = "OS",
    model = "continuous_Score_z",
    n = n_os,
    events = events_os
  )

  cox_gleason_tab <- NULL

  complete_gleason <- df_os %>%
    filter(
      is.finite(Score_z_s),
      is.finite(gleason_s)
    )

  if (
    nrow(complete_gleason) >= 30 &&
    sum(
      complete_gleason$event == 1L,
      na.rm = TRUE
    ) >= 20
  ) {
    cox_gleason <- survival::coxph(
      survival::Surv(
        time,
        event
      ) ~ Score_z_s + gleason_s,
      data = complete_gleason
    )

    cox_gleason_tab <- tidy_cox_or_empty(
      cox_gleason,
      endpoint = "OS",
      model = "gleason_adjusted",
      n = nrow(complete_gleason),
      events = sum(
        complete_gleason$event == 1L,
        na.rm = TRUE
      )
    )

  } else {
    add_warning(
      paste0(
        "OS Gleason-adjusted Cox skipped due to ",
        "event or complete-case limitations."
      )
    )
  }

  cox_all <- bind_rows(
    cox_bin_tab,
    cox_cont_tab,
    cox_gleason_tab
  )

  survival_summary_all <<- bind_rows(
    km_summary,

    cox_all %>%
      transmute(
        endpoint,

        analysis = case_when(
          model == "binary_median" ~
            "Cox_median_split",

          model == "continuous_Score_z" ~
            "Cox_continuous_score",

          model == "gleason_adjusted" ~
            "Cox_gleason_adjusted",

          TRUE ~ model
        ),

        model,
        n,
        events,
        term,
        estimate,
        conf.low,
        conf.high,

        p_value = p.value,

        p_label = fmt_p(
          p.value
        ),

        median_cutpoint = NA_real_,
        n_low = NA_integer_,
        n_high = NA_integer_
      )
  )

  write_table(
    survival_summary_all,
    file.path(
      tab_dir,
      "Survival_summary_OS.csv"
    )
  )

  cox_bin_plot_tbl <- broom::tidy(
    cox_bin,
    exponentiate = TRUE,
    conf.int = TRUE
  )

  hr_txt <- NA_character_

  if (nrow(cox_bin_plot_tbl) == 1) {
    hr_txt <- paste0(
      "HR = ",
      sprintf(
        "%.2f",
        cox_bin_plot_tbl$estimate[1]
      ),
      " (95% CI ",
      sprintf(
        "%.2f",
        cox_bin_plot_tbl$conf.low[1]
      ),
      "-",
      sprintf(
        "%.2f",
        cox_bin_plot_tbl$conf.high[1]
      ),
      ")"
    )
  }

  p_txt <- if (p_lr < 0.001) {
    "'Log-rank'~italic(P) < 0.001"
  } else {
    paste0(
      "'Log-rank'~italic(P) == '",
      sprintf("%.3f", p_lr),
      "'"
    )
  }

  df_os_plot <- df_os %>%
    mutate(
      score_group_plot = factor(
        ifelse(
          score_group == "Low",
          "Low",
          "High"
        ),
        levels = c(
          "Low",
          "High"
        )
      )
    )

  fit_km_plot <- survival::survfit(
    survival::Surv(
      time,
      event
    ) ~ score_group_plot,
    data = df_os_plot
  )

  cols <- c(
    "Low" = col_low,
    "High" = col_high
  )

  gp <- survminer::ggsurvplot(
    fit_km_plot,
    data = df_os_plot,
    risk.table = FALSE,
    conf.int = FALSE,
    palette = unname(cols),
    legend = "none",
    pval = FALSE,
    size = 0.36,
    censor.size = 0.75,
    censor.shape = 124,
    break.time.by = 50,
    xlab = "Months",
    ylab = "Overall survival (%)",
    ggtheme = theme_classic(
      base_family = font_use,
      base_size = pt_axis
    )
  )

  xmax <- 200

  risk_times <- c(
    0,
    50,
    100,
    150,
    200
  )

  risk_sum <- summary(
    fit_km_plot,
    times = risk_times,
    extend = TRUE
  )

  risk_tbl <- data.frame(
    time = risk_sum$time,
    n_risk = risk_sum$n.risk,
    strata = as.character(
      risk_sum$strata
    ),
    stringsAsFactors = FALSE
  )

  risk_tbl$strata <- sub(
    "^score_group_plot=",
    "",
    risk_tbl$strata
  )

  risk_tbl$y <- ifelse(
    risk_tbl$strata == "Low",
    -0.235,
    -0.305
  )

  risk_tbl$col <- ifelse(
    risk_tbl$strata == "Low",
    cols["Low"],
    cols["High"]
  )

  x_annot <- xmax * 0.020
  y_top <- 0.11
  dy <- 0.070

  label_low <- paste0(
    "Platelet score <= median (n=",
    sum(
      df_os$score_group == "Low",
      na.rm = TRUE
    ),
    ")"
  )

  label_high <- paste0(
    "Platelet score > median (n=",
    sum(
      df_os$score_group == "High",
      na.rm = TRUE
    ),
    ")"
  )

  risk_title_y <- -0.165

  risk_tbl_low <- risk_tbl %>%
    filter(
      strata == "Low"
    )

  risk_tbl_high <- risk_tbl %>%
    filter(
      strata == "High"
    )

  km_plot <- gp$plot +
    coord_cartesian(
      xlim = c(
        0,
        xmax
      ),
      ylim = c(
        0,
        1
      ),
      expand = FALSE,
      clip = "off"
    ) +
    scale_x_continuous(
      breaks = risk_times,
      limits = c(
        0,
        xmax
      ),
      expand = c(
        0,
        0
      ),
      oob = scales::squish
    ) +
    scale_y_continuous(
      breaks = c(
        0,
        0.25,
        0.5,
        0.75,
        1
      ),
      labels = c(
        "0",
        "25",
        "50",
        "75",
        "100"
      )
    ) +
    theme(
      text = element_text(
        family = font_use,
        color = "black"
      ),
      panel.grid = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.20
      ),
      axis.line = element_blank(),
      axis.text = element_text(
        size = 4.0,
        color = "black"
      ),
      axis.title = element_text(
        size = 4.5,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.20,
        color = "black"
      ),
      axis.ticks.length = unit(
        0.5,
        "mm"
      ),
      plot.margin = margin(
        2.2,
        1.8,
        8.0,
        1.8,
        "mm"
      )
    ) +
    annotate(
      "text",
      x = x_annot,
      y = y_top + 2 * dy,
      label = label_low,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      color = col_low
    ) +
    annotate(
      "text",
      x = x_annot,
      y = y_top + 1 * dy,
      label = label_high,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      color = col_high
    ) +
    annotate(
      "text",
      x = x_annot,
      y = y_top + 0 * dy,
      label = hr_txt,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      color = "black"
    ) +
    annotate(
      "text",
      x = x_annot,
      y = y_top - 1.20 * dy,
      label = p_txt,
      parse = TRUE,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      color = "black"
    ) +
    annotate(
      "text",
      x = xmax * 0.5,
      y = risk_title_y,
      label = "Number at risk",
      hjust = 0.5,
      vjust = 0.5,
      size = 4.4 / 2.845,
      color = "black"
    ) +
    geom_text(
      data = risk_tbl_low,
      aes(
        x = time,
        y = y,
        label = n_risk
      ),
      inherit.aes = FALSE,
      size = 4.5 / 2.845,
      family = font_use,
      color = col_low,
      vjust = 0.5
    ) +
    geom_text(
      data = risk_tbl_high,
      aes(
        x = time,
        y = y,
        label = n_risk
      ),
      inherit.aes = FALSE,
      size = 4.5 / 2.845,
      family = font_use,
      color = col_high,
      vjust = 0.5
    )

  save_pdf(
    km_plot,
    file.path(
      pdf_dir,
      "KM_OS_medianStrat.pdf"
    ),
    width_cm = 4.2,
    height_cm = 4.6
  )
})

###############################################################################
# 12) REQUIRED OUTPUT VALIDATION
###############################################################################

required_output_paths <- c(
  file.path(
    tab_dir,
    "Friedrich_GSE134051_clinical_final_master_table.csv"
  ),

  file.path(
    tab_dir,
    "Clinical_variable_completeness.csv"
  ),

  file.path(
    tab_dir,
    "Endpoint_event_counts.csv"
  ),

  file.path(
    tab_dir,
    "Clinical_assoc_minimal_table.csv"
  ),

  file.path(
    tab_dir,
    "PlateletScore_grouped_Gleason_score_summary.csv"
  ),

  file.path(
    tab_dir,
    "PlateletScore_grouped_Gleason_score_Spearman.csv"
  ),

  file.path(
    tab_dir,
    "Survival_summary_OS.csv"
  ),

  file.path(
    pdf_dir,
    "PlateletScore_vs_grouped_Gleason_score.pdf"
  ),

  file.path(
    pdf_dir,
    "KM_OS_medianStrat.pdf"
  ),

  file.path(
    rds_dir,
    "Friedrich_GSE134051_clinical_final_master_table.rds"
  )
)

missing_required_outputs <- required_output_paths[
  !file.exists(
    required_output_paths
  )
]

###############################################################################
# 13) FINAL REPORT
###############################################################################

status <- if (
  length(missing_required_outputs) > 0
) {
  "FAIL"
} else if (
  length(warnings_vec) == 0
) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

report_lines <- c(
  "# Friedrich/GSE134051 Clinical Final Report",
  "",
  paste0(
    "Date: ",
    as.character(Sys.Date())
  ),
  paste0(
    "Status: ",
    status
  ),
  "",
  "## Script",
  script_fp,
  "",
  "## Inputs",
  paste0(
    "- Official 41-gene platelet score: ",
    score_fp
  ),
  paste0(
    "- Friedrich/GSE134051 clinical table: ",
    clinical_fp
  ),
  "",
  "## Official Score Confirmation",
  paste0(
    "- The active score source is the official ",
    "41-gene platelet-associated transcriptional score file."
  ),
  "- Score_z is read directly from the official score table.",
  "- Legacy platelet score files are not used as active inputs.",
  "",
  "## Counts",
  paste0(
    "- Expected primary samples: ",
    EXPECTED_PRIMARY_N
  ),
  paste0(
    "- Final master table samples: ",
    nrow(df_clinical_final)
  ),
  paste0(
    "- Primary clinical samples: ",
    nrow(clin_primary)
  ),
  "",
  "## Endpoint Events",
  paste(
    apply(
      endpoint_event_counts,
      1,
      function(x) {
        paste0(
          "- ",
          x[["endpoint"]],
          ": ",
          x[["events"]],
          " events / ",
          x[["n_with_time"]],
          " with time"
        )
      }
    ),
    collapse = "\n"
  ),
  "",
  "## Clinical Availability",
  paste0(
    "- Available for analysis: Gleason score, ",
    "grouped Gleason score, OS."
  ),
  paste0(
    "- Not available in Friedrich/GSE134051 primary samples: age, PSA, ",
    "T stage, N stage, M stage, disease-specific recurrence, ",
    "metastatic occurrence."
  ),
  "",
  "## Main Method Note",
  paste0(
    "- Clinical analyses in Friedrich/GSE134051 were performed ",
    "as exploratory association analyses."
  ),
  paste0(
    "- The platelet-associated transcriptional score was ",
    "analyzed as a continuous standardized variable for Cox models."
  ),
  paste0(
    "- Median stratification was used only for ",
    "Kaplan-Meier visualization."
  ),
  paste0(
    "- OS analyses are exploratory due to the low event count."
  ),
  paste0(
    "- Clinical outputs should not be interpreted as ",
    "independent clinical biomarker validation."
  ),
  "",
  "## Required Outputs",
  if (length(missing_required_outputs) == 0) {
    "- All required outputs were generated."
  } else {
    paste0(
      "- Missing: ",
      missing_required_outputs
    )
  },
  "",
  "## Tables Generated",
  if (length(tables_generated) == 0) {
    "- None"
  } else {
    paste0(
      "- ",
      tables_generated
    )
  },
  "",
  "## Figures Generated",
  if (length(figures_generated) == 0) {
    "- None"
  } else {
    paste0(
      "- ",
      figures_generated
    )
  },
  "",
  "## RDS Generated",
  if (length(rds_generated) == 0) {
    "- None"
  } else {
    paste0(
      "- ",
      rds_generated
    )
  },
  "",
  "## Warnings",
  if (length(warnings_vec) == 0) {
    "- None"
  } else {
    paste0(
      "- ",
      warnings_vec
    )
  }
)

writeLines(
  report_lines,
  file.path(
    report_dir,
    "REPORT_Friedrich_GSE134051_clinical_final.md"
  ),
  useBytes = TRUE
)

qc_lines <- c(
  "Friedrich/GSE134051 Clinical Final QC",
  "==========================",
  "",
  paste0(
    "Date/time: ",
    as.character(Sys.time())
  ),
  paste0(
    "Status: ",
    status
  ),
  "",
  "Inputs:",
  paste0(
    " - Official 41-gene platelet score: ",
    score_fp
  ),
  paste0(
    " - Friedrich/GSE134051 clinical table: ",
    clinical_fp
  ),
  "",
  "Counts:",
  paste0(
    " - Expected primary samples: ",
    EXPECTED_PRIMARY_N
  ),
  paste0(
    " - Primary clinical samples: ",
    nrow(clin_primary)
  ),
  paste0(
    " - Final master table samples: ",
    nrow(df_clinical_final)
  ),
  "",
  "Endpoint events:",
  paste(
    apply(
      endpoint_event_counts,
      1,
      function(x) {
        paste0(
          " - ",
          x[["endpoint"]],
          ": ",
          x[["events"]],
          " events / ",
          x[["n_with_time"]],
          " with time"
        )
      }
    ),
    collapse = "\n"
  ),
  "",
  "Required outputs:",
  if (length(missing_required_outputs) == 0) {
    " - All required outputs present"
  } else {
    paste0(
      " - Missing: ",
      missing_required_outputs
    )
  },
  "",
  "Tables generated:",
  if (length(tables_generated) == 0) {
    " - None"
  } else {
    paste0(
      " - ",
      tables_generated
    )
  },
  "",
  "Figures generated:",
  if (length(figures_generated) == 0) {
    " - None"
  } else {
    paste0(
      " - ",
      figures_generated
    )
  },
  "",
  "RDS generated:",
  if (length(rds_generated) == 0) {
    " - None"
  } else {
    paste0(
      " - ",
      rds_generated
    )
  },
  "",
  "Warnings:",
  if (length(warnings_vec) == 0) {
    " - None"
  } else {
    paste0(
      " - ",
      warnings_vec
    )
  },
  "",
  paste0(
    "Final status: ",
    status
  )
)

writeLines(
  qc_lines,
  file.path(
    report_dir,
    "Friedrich_GSE134051_clinical_final_QC.txt"
  ),
  useBytes = TRUE
)

if (status == "FAIL") {
  stop(
    paste0(
      "Friedrich/GSE134051 clinical final analysis failed because required outputs ",
      "were not generated: ",
      paste(
        missing_required_outputs,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

message(
  "Friedrich/GSE134051 clinical final script completed with status: ",
  status
)
