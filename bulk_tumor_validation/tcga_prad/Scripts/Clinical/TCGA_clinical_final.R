#!/usr/bin/env Rscript

################################################################################
# TCGA-PRAD clinical analysis using the canonical 41-gene
# platelet-associated transcriptional score
#
# Script:
# Scripts/Clinical/TCGA_clinical_final.R
#
# Analytical contract:
#   - Canonical cohort: 497 primary tumors, one sample per patient.
#   - Median stratification is used only for Kaplan–Meier visualization.
#   - Canonical Score_z is used directly in continuous Cox models.
#   - Score_z is not re-standardized within endpoint-specific subsets.
#   - BCR is the principal TCGA survival endpoint.
#   - PFI and RFS are complementary endpoints.
#   - OS is exploratory because of the low event count.
#   - Clinical results represent association analyses, not biomarker validation.
################################################################################

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  required_packages <- c(
    "SummarizedExperiment",
    "dplyr",
    "tibble",
    "stringr",
    "ggplot2",
    "readxl",
    "survival",
    "survminer",
    "broom"
  )
})

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
# 1. Constants and frozen visual style
################################################################################

EXPECTED_CANONICAL_SAMPLES <- 497L

FONT_FAMILY <- "Helvetica"

COL_LOW <- "#0F558F"
COL_HIGH <- "#E47513"
COL_GLEASON <- "#3E6A8E"
COL_GLEASON_LINE <- "#243B6B"

AXIS_TEXT_PT <- 4.0
AXIS_TITLE_PT <- 4.5
MAIN_TEXT_PT <- 5.0
STAT_TEXT_PT <- 4.0

PANEL_BORDER_LW <- 0.20
AXIS_TICK_LW <- 0.20
KM_LINE_LW <- 0.36
KM_CENSOR_SIZE <- 0.75
KM_CENSOR_SHAPE <- 124

KM_XMAX <- 160
KM_TIMES <- c(0, 40, 80, 120, 160)
KM_BREAK_BY <- 40

MIN_MULTIVARIABLE_EVENTS <- 15L
MIN_MULTIVARIABLE_COMPLETE_CASES <- 30L

################################################################################
# 2. Paths
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

script_fp <- file.path(
  project_dir,
    "Scripts",
    "Clinical",
    "TCGA_clinical_final.R"
)

score_fp <- file.path(
  GENERATED_RESULTS_DIR, "Score",
  "TCGA_primary_master_table_41genes.csv"
)

se_fp <- file.path(
  INPUT_DIR, "LocalLarge",
  "TCGA_PRAD_SE_tumorOnly.rds"
)

clin_txt_fp <- file.path(
  INPUT_DIR, "Clinical",
  "PRAD_clin_merged.txt"
)

cdr_fp <- file.path(
  INPUT_DIR, "Clinical",
  "TCGA-CDR-SupplementalTableS1.xlsx"
)

results_dir <- file.path(GENERATED_RESULTS_DIR, "Clinical")
table_dir <- file.path(results_dir, "Tables")
figure_dir <- file.path(GENERATED_FIGURES_DIR, "Clinical")
object_dir <- file.path(results_dir, "Objects")
log_dir <- file.path(results_dir, "Logs")

invisible(
  lapply(
    c(results_dir, table_dir, figure_dir, object_dir, log_dir),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

output_paths <- list(
  master_csv = file.path(
    table_dir,
    "TCGA_PRAD_clinical_master_table.csv"
  ),

  master_rds = file.path(
    object_dir,
    "TCGA_PRAD_clinical_master_table.rds"
  ),

  clinical_associations = file.path(
    table_dir,
    "Clinical_associations.csv"
  ),

  clinical_availability = file.path(
    table_dir,
    "Clinical_variable_availability.csv"
  ),

  gleason_counts = file.path(
    table_dir,
    "PlateletScore_Gleason_group_counts.csv"
  ),

  gleason_summary = file.path(
    table_dir,
    "PlateletScore_Gleason_summary.csv"
  ),

  gleason_spearman = file.path(
    table_dir,
    "PlateletScore_Gleason_Spearman.csv"
  ),

  t_stage_wilcoxon = file.path(
    table_dir,
    "PlateletScore_T2_vs_T3plus_Wilcoxon.csv"
  ),

  endpoint_counts = file.path(
    table_dir,
    "Endpoint_analysis_counts.csv"
  ),

  survival_summary = file.path(
    table_dir,
    "Survival_summary_all_endpoints.csv"
  ),

  cox_ph_diagnostics = file.path(
    table_dir,
    "Cox_proportional_hazards_diagnostics.csv"
  ),

  figure_gleason = file.path(
    figure_dir,
    "PlateletScore_vs_Gleason.pdf"
  ),

  figure_t_stage = file.path(
    figure_dir,
    "PlateletScore_T2_vs_T3plus.pdf"
  ),

  figure_bcr_km = file.path(
    figure_dir,
    "KM_BCR_medianStrat.pdf"
  ),

  qc_log = file.path(
    log_dir,
    "01_TCGA_clinical_final_QC.txt"
  )
)

################################################################################
# 3. General helpers
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
tables_generated <- character()
figures_generated <- character()
objects_generated <- character()

add_warning <- function(...) {
  warning_text <- paste0(...)
  warnings_log <<- unique(c(warnings_log, warning_text))
  warning(warning_text, call. = FALSE)
  invisible(warning_text)
}

require_input <- function(path, label) {
  if (!file.exists(path)) {
    fail("Missing ", label, ": ", path)
  }
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

write_table <- function(df, path) {
  utils::write.csv(
    as.data.frame(df),
    path,
    row.names = FALSE
  )

  tables_generated <<- unique(
    c(tables_generated, path)
  )

  invisible(path)
}

write_rds <- function(object, path) {
  saveRDS(object, path)

  objects_generated <<- unique(
    c(objects_generated, path)
  )

  invisible(path)
}

save_pdf <- function(plot, path, width_cm, height_cm) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = grDevices::pdf,
    family = FONT_FAMILY,
    useDingbats = FALSE
  )

  figures_generated <<- unique(
    c(figures_generated, path)
  )

  invisible(path)
}

as_numeric_clean <- function(x) {
  suppressWarnings(
    as.numeric(as.character(x))
  )
}

normalize_barcode <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("\\.", "-", x)
}

normalize_patient <- function(x) {
  x <- normalize_barcode(x)

  stringr::str_extract(
    x,
    "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"
  )
}

normalize_t_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- stringr::str_replace_all(x, "\\s+", "")

  stringr::str_extract(
    x,
    "T[0-9]"
  )
}

normalize_n_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- stringr::str_replace_all(x, "\\s+", "")

  stringr::str_extract(
    x,
    "N[0-9]"
  )
}

pick_first <- function(names_vector, candidates) {
  hit <- candidates[candidates %in% names_vector]

  if (length(hit) == 0L) {
    return(NA_character_)
  }

  hit[1L]
}

pick_first_grep <- function(names_vector, patterns) {
  for (pattern in patterns) {
    hit <- grep(
      pattern,
      names_vector,
      value = TRUE,
      ignore.case = TRUE
    )

    if (length(hit) > 0L) {
      return(hit[1L])
    }
  }

  NA_character_
}

format_p <- function(p_value, digits = 2L) {
  if (!is.finite(p_value)) {
    return("NA")
  }

  if (p_value < 0.001) {
    return(
      formatC(
        p_value,
        format = "e",
        digits = digits
      )
    )
  }

  sprintf("%.3f", p_value)
}

format_p_label <- function(p_value) {
  paste0("P = ", format_p(p_value))
}

positive_status <- function(x) {
  x <- tolower(trimws(as.character(x)))

  x %in% c(
    "yes",
    "y",
    "true",
    "1",
    "dead",
    "deceased",
    "recurred",
    "recurrence",
    "progressed",
    "progression",
    "event"
  )
}

negative_status <- function(x) {
  x <- tolower(trimws(as.character(x)))

  x %in% c(
    "no",
    "n",
    "false",
    "0",
    "alive",
    "disease free",
    "disease-free",
    "no recurrence",
    "no progression"
  )
}

minimum_finite <- function(x) {
  x <- as_numeric_clean(x)
  x <- x[is.finite(x) & x >= 0]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  min(x)
}

maximum_finite <- function(x) {
  x <- as_numeric_clean(x)
  x <- x[is.finite(x) & x >= 0]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  max(x)
}

first_finite <- function(x) {
  x <- as_numeric_clean(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  x[1L]
}

first_nonmissing_integer <- function(x) {
  x <- as.integer(x)
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return(NA_integer_)
  }

  x[1L]
}

collapse_binary_status <- function(status, event_date = NULL) {
  if (
    !is.null(event_date) &&
    any(is.finite(as_numeric_clean(event_date)))
  ) {
    return(1L)
  }

  if (any(positive_status(status), na.rm = TRUE)) {
    return(1L)
  }

  if (any(negative_status(status), na.rm = TRUE)) {
    return(0L)
  }

  NA_integer_
}

coerce_binary_event <- function(x) {
  output <- rep(NA_integer_, length(x))
  numeric_x <- as_numeric_clean(x)

  finite_numeric <- is.finite(numeric_x)

  output[
    finite_numeric & numeric_x == 0
  ] <- 0L

  output[
    finite_numeric & numeric_x > 0
  ] <- 1L

  unresolved <- is.na(output)

  output[
    unresolved & positive_status(x)
  ] <- 1L

  output[
    unresolved & negative_status(x)
  ] <- 0L

  output
}

get_numeric_date_column <- function(df, pattern, reducer) {
  columns <- grep(
    pattern,
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(columns) == 0L) {
    return(rep(NA_real_, nrow(df)))
  }

  value_matrix <- sapply(
    df[, columns, drop = FALSE],
    as_numeric_clean
  )

  if (is.null(dim(value_matrix))) {
    value_matrix <- matrix(
      value_matrix,
      ncol = 1L
    )
  }

  apply(
    value_matrix,
    1L,
    reducer
  )
}

extract_grade_pattern <- function(x) {
  as_numeric_clean(
    stringr::str_extract(
      as.character(x),
      "[0-9]+"
    )
  )
}

convert_age_to_years <- function(x) {
  x <- as_numeric_clean(x)

  finite_values <- abs(x[is.finite(x)])

  if (length(finite_values) == 0L) {
    return(rep(NA_real_, length(x)))
  }

  if (stats::median(finite_values) > 150) {
    return(abs(x) / 365.25)
  }

  x
}

make_gleason_group <- function(primary, secondary, total) {
  primary <- as_numeric_clean(primary)
  secondary <- as_numeric_clean(secondary)
  total <- as_numeric_clean(total)

  dplyr::case_when(
    is.finite(primary) &
      is.finite(secondary) &
      primary == 3 &
      secondary == 3 ~ "6",

    is.finite(primary) &
      is.finite(secondary) &
      primary == 3 &
      secondary == 4 ~ "3+4",

    is.finite(primary) &
      is.finite(secondary) &
      primary == 4 &
      secondary == 3 ~ "4+3",

    is.finite(primary) &
      is.finite(secondary) &
      primary + secondary == 8 ~ "8",

    is.finite(primary) &
      is.finite(secondary) &
      primary + secondary >= 9 ~ ">=9",

    (!is.finite(primary) | !is.finite(secondary)) &
      total <= 6 ~ "6",

    (!is.finite(primary) | !is.finite(secondary)) &
      total == 7 ~ "7",

    (!is.finite(primary) | !is.finite(secondary)) &
      total == 8 ~ "8",

    (!is.finite(primary) | !is.finite(secondary)) &
      total >= 9 ~ ">=9",

    TRUE ~ NA_character_
  )
}

################################################################################
# 4. Input checks
################################################################################

require_input(
  score_fp,
  "canonical 41-gene platelet-score table"
)

require_input(
  se_fp,
  "TCGA-PRAD tumor-only SummarizedExperiment"
)

require_input(
  clin_txt_fp,
  "ClinPlus clinical text file"
)

require_input(
  cdr_fp,
  "TCGA-CDR supplemental table"
)

msg("Project directory:", project_dir)

################################################################################
# 5. Load canonical score cohort
################################################################################

msg("Loading canonical 497-patient score cohort")

score_df <- utils::read.csv(
  score_fp,
  stringsAsFactors = FALSE,
  check.names = FALSE
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

score_df <- score_df |>
  dplyr::transmute(
    sample_id = normalize_barcode(sample_id),
    patient_id = normalize_patient(patient_id),
    Score_raw = as_numeric_clean(Score_raw),
    Score_z = as_numeric_clean(Score_z),
    platelet_score_group = if (
      "platelet_score_group" %in% colnames(score_df)
    ) {
      as.character(platelet_score_group)
    } else {
      NA_character_
    }
  )

if (nrow(score_df) != EXPECTED_CANONICAL_SAMPLES) {
  fail(
    "Expected ",
    EXPECTED_CANONICAL_SAMPLES,
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
    "Canonical score table must contain ",
    EXPECTED_CANONICAL_SAMPLES,
    " unique sample_id values."
  )
}

if (
  anyNA(score_df$patient_id) ||
  any(score_df$patient_id == "") ||
  anyDuplicated(score_df$patient_id)
) {
  fail(
    "Canonical score table must contain ",
    EXPECTED_CANONICAL_SAMPLES,
    " unique patient_id values."
  )
}

patient_from_sample <- normalize_patient(
  score_df$sample_id
)

if (
  anyNA(patient_from_sample) ||
  !identical(
    patient_from_sample,
    score_df$patient_id
  )
) {
  fail(
    "Canonical patient_id values are not identical to the patient ",
    "barcodes derived from sample_id."
  )
}

if (
  anyNA(score_df$Score_raw) ||
  anyNA(score_df$Score_z) ||
  any(!is.finite(score_df$Score_raw)) ||
  any(!is.finite(score_df$Score_z))
) {
  fail("Score_raw and Score_z must be finite for all canonical patients.")
}

if (abs(mean(score_df$Score_z)) > 1e-8) {
  fail("Canonical Score_z mean differs from zero beyond tolerance.")
}

if (abs(stats::sd(score_df$Score_z) - 1) > 1e-8) {
  fail("Canonical Score_z standard deviation differs from one.")
}

canonical_score_median <- stats::median(
  score_df$Score_z
)

################################################################################
# 6. Load and align canonical SE clinical covariates
################################################################################

msg("Loading and aligning tumor-only SummarizedExperiment")

se_tumor <- readRDS(se_fp)

if (!inherits(se_tumor, "SummarizedExperiment")) {
  fail(
    "Tumor-only SE input is not a SummarizedExperiment."
  )
}

se_coldata <- as.data.frame(
  SummarizedExperiment::colData(se_tumor),
  stringsAsFactors = FALSE
)

se_sample_ids <- if (
  "barcode" %in% colnames(se_coldata)
) {
  normalize_barcode(se_coldata$barcode)
} else {
  normalize_barcode(colnames(se_tumor))
}

if (
  length(se_sample_ids) != nrow(se_coldata) ||
  anyNA(se_sample_ids) ||
  any(se_sample_ids == "") ||
  anyDuplicated(se_sample_ids)
) {
  fail(
    "SE clinical metadata does not contain unique valid sample IDs."
  )
}

se_index <- match(
  score_df$sample_id,
  se_sample_ids
)

if (anyNA(se_index)) {
  fail(
    "Exact canonical sample alignment failed for SE metadata. ",
    "Missing examples: ",
    paste(
      utils::head(
        score_df$sample_id[is.na(se_index)],
        8L
      ),
      collapse = ", "
    ),
    ". Barcode truncation is not allowed."
  )
}

se_coldata <- se_coldata[
  se_index,
  ,
  drop = FALSE
]

se_coldata$sample_id <- se_sample_ids[se_index]

if (!identical(se_coldata$sample_id, score_df$sample_id)) {
  fail("Aligned SE metadata does not preserve canonical sample order.")
}

se_names <- colnames(se_coldata)

age_col <- pick_first_grep(
  se_names,
  c(
    "^age_at_diagnosis$",
    "^age_at_index$",
    "^paper_Age$",
    "days_to_birth",
    "age_at_initial_pathologic",
    "^age$"
  )
)

psa_col <- pick_first_grep(
  se_names,
  c(
    "preoperative.*psa",
    "prostate_specific_antigen",
    "^psa$",
    "psa"
  )
)

gleason_primary_col <- pick_first_grep(
  se_names,
  c(
    "^primary_gleason_grade$",
    "primary.*gleason",
    "gleason.*primary",
    "primary.*pattern"
  )
)

gleason_secondary_col <- pick_first_grep(
  se_names,
  c(
    "^secondary_gleason_grade$",
    "secondary.*gleason",
    "gleason.*secondary",
    "secondary.*pattern"
  )
)

gleason_sum_col <- pick_first_grep(
  se_names,
  c(
    "^gleason_score$",
    "gleason.*sum",
    "gleason.*combined",
    "gleason.*total"
  )
)

t_col <- pick_first_grep(
  se_names,
  c(
    "^ajcc_pathologic_t$",
    "^paper_ajcc_pathologic_t$",
    "pathologic.*t",
    "ajcc.*pathologic.*t"
  )
)

n_col <- pick_first_grep(
  se_names,
  c(
    "^ajcc_pathologic_n$",
    "^paper_ajcc_pathologic_n$",
    "pathologic.*n",
    "ajcc.*pathologic.*n"
  )
)

ln_col <- pick_first_grep(
  se_names,
  c(
    "positive.*lymph",
    "lymph.*positive",
    "nodes.*positive",
    "lymph_node.*count"
  )
)

detected_columns <- tibble::tibble(
  variable = c(
    "age",
    "PSA",
    "gleason_primary",
    "gleason_secondary",
    "gleason_sum",
    "pathologic_T",
    "pathologic_N",
    "positive_lymph_nodes"
  ),
  source_column = c(
    age_col,
    psa_col,
    gleason_primary_col,
    gleason_secondary_col,
    gleason_sum_col,
    t_col,
    n_col,
    ln_col
  )
)

missing_detected <- detected_columns$variable[
  is.na(detected_columns$source_column)
]

if (length(missing_detected) > 0L) {
  add_warning(
    "SE clinical columns not detected for: ",
    paste(missing_detected, collapse = ", "),
    ". Corresponding analyses will use missing values."
  )
}

age <- if (!is.na(age_col)) {
  convert_age_to_years(se_coldata[[age_col]])
} else {
  rep(NA_real_, nrow(se_coldata))
}

psa <- if (!is.na(psa_col)) {
  as_numeric_clean(se_coldata[[psa_col]])
} else {
  rep(NA_real_, nrow(se_coldata))
}

gleason_primary <- if (!is.na(gleason_primary_col)) {
  extract_grade_pattern(
    se_coldata[[gleason_primary_col]]
  )
} else {
  rep(NA_real_, nrow(se_coldata))
}

gleason_secondary <- if (!is.na(gleason_secondary_col)) {
  extract_grade_pattern(
    se_coldata[[gleason_secondary_col]]
  )
} else {
  rep(NA_real_, nrow(se_coldata))
}

gleason_sum <- if (!is.na(gleason_sum_col)) {
  as_numeric_clean(
    se_coldata[[gleason_sum_col]]
  )
} else {
  rep(NA_real_, nrow(se_coldata))
}

pattern_sum <- gleason_primary + gleason_secondary

gleason_sum[
  !is.finite(gleason_sum) &
    is.finite(pattern_sum)
] <- pattern_sum[
  !is.finite(gleason_sum) &
    is.finite(pattern_sum)
]

gleason_group <- make_gleason_group(
  gleason_primary,
  gleason_secondary,
  gleason_sum
)

if (any(gleason_group == "7", na.rm = TRUE)) {
  add_warning(
    "Some Gleason sum-7 tumors could not be separated into 3+4 versus 4+3 ",
    "because primary/secondary patterns were unavailable."
  )
}

pathologic_t <- if (!is.na(t_col)) {
  normalize_t_stage(
    se_coldata[[t_col]]
  )
} else {
  rep(NA_character_, nrow(se_coldata))
}

pathologic_n <- if (!is.na(n_col)) {
  normalize_n_stage(
    se_coldata[[n_col]]
  )
} else {
  rep(NA_character_, nrow(se_coldata))
}

positive_lymph_nodes <- if (!is.na(ln_col)) {
  as_numeric_clean(
    se_coldata[[ln_col]]
  )
} else {
  rep(NA_real_, nrow(se_coldata))
}

se_covariates <- tibble::tibble(
  sample_id = score_df$sample_id,
  patient_id = score_df$patient_id,
  age = age,
  PSA = psa,
  gleason_primary = gleason_primary,
  gleason_secondary = gleason_secondary,
  gleason_sum = gleason_sum,
  gleason_group = gleason_group,
  T = pathologic_t,
  N = pathologic_n,
  LN_pos = positive_lymph_nodes
) |>
  dplyr::mutate(
    gleason_group = factor(
      gleason_group,
      levels = c(
        "6",
        "3+4",
        "4+3",
        "7",
        "8",
        ">=9"
      ),
      ordered = TRUE
    ),

    T_bin = dplyr::case_when(
      T == "T2" ~ "T2",
      T %in% c("T3", "T4") ~ "T3+",
      TRUE ~ NA_character_
    ),

    N_bin = dplyr::case_when(
      N == "N0" ~ "N0",
      N %in% c("N1", "N2", "N3") ~ "N+",
      TRUE ~ NA_character_
    )
  )

if (
  anyDuplicated(se_covariates$sample_id) ||
  anyDuplicated(se_covariates$patient_id)
) {
  fail(
    "Aligned SE clinical table contains duplicated canonical samples or patients."
  )
}

################################################################################
# 7. ClinPlus OS, RFS and BCR endpoints
################################################################################

msg("Loading ClinPlus OS, RFS and BCR endpoints")

clin_raw <- utils::read.delim(
  clin_txt_fp,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clin_t <- as.data.frame(
  t(clin_raw),
  stringsAsFactors = FALSE
)

colnames(clin_t) <- make.names(
  as.character(
    unlist(clin_t[1L, ])
  ),
  unique = TRUE
)

clin_t <- clin_t[-1L, , drop = FALSE]
clin_t$sample_id_raw <- rownames(clin_t)

clin_t <- tibble::as_tibble(clin_t)

clin_t[] <- lapply(
  clin_t,
  function(x) {
    x <- trimws(as.character(x))

    missing_tokens <- c(
      "",
      "NA",
      "N/A",
      "[Not Available]",
      "[Not Applicable]",
      "[Unknown]",
      "unknown",
      "not reported"
    )

    x[x %in% missing_tokens] <- NA_character_
    x
  }
)

patient_col <- pick_first_grep(
  colnames(clin_t),
  c(
    "^patient\\.bcr_patient_barcode$",
    "bcr_patient_barcode",
    "patient_barcode"
  )
)

vital_col <- pick_first_grep(
  colnames(clin_t),
  c(
    "^patient\\.vital_status$",
    "vital_status"
  )
)

bcr_status_col <- pick_first_grep(
  colnames(clin_t),
  c(
    "biochemical_recurrence"
  )
)

if (is.na(patient_col)) {
  fail(
    "No patient barcode column was found in ClinPlus."
  )
}

days_to_death <- get_numeric_date_column(
  clin_t,
  "days_to_death",
  minimum_finite
)

days_to_last_followup <- get_numeric_date_column(
  clin_t,
  "days_to_last_followup|days_to_last_follow_up",
  maximum_finite
)

days_to_new_tumor <- get_numeric_date_column(
  clin_t,
  "days_to_new_tumor_event",
  minimum_finite
)

days_to_bcr <- get_numeric_date_column(
  clin_t,
  paste0(
    "days_to_first_biochemical_recurrence|",
    "days_to_biochemical_recurrence"
  ),
  minimum_finite
)

vital_status <- if (!is.na(vital_col)) {
  clin_t[[vital_col]]
} else {
  rep(NA_character_, nrow(clin_t))
}

bcr_status <- if (!is.na(bcr_status_col)) {
  clin_t[[bcr_status_col]]
} else {
  rep(NA_character_, nrow(clin_t))
}

clin_rows <- tibble::tibble(
  patient_id = normalize_patient(
    clin_t[[patient_col]]
  ),

  sample_id_raw = clin_t$sample_id_raw,
  vital_status = vital_status,
  bcr_status = bcr_status,
  days_to_death = days_to_death,
  days_to_last_followup = days_to_last_followup,
  days_to_new_tumor = days_to_new_tumor,
  days_to_bcr = days_to_bcr
) |>
  dplyr::filter(
    !is.na(patient_id)
  )

clin_endpoints <- clin_rows |>
  dplyr::group_by(patient_id) |>
  dplyr::summarise(
    days_to_death = minimum_finite(days_to_death),
    days_to_last_followup = maximum_finite(days_to_last_followup),
    days_to_new_tumor = minimum_finite(days_to_new_tumor),
    days_to_bcr = minimum_finite(days_to_bcr),

    death_flag = collapse_binary_status(
      vital_status,
      days_to_death
    ),

    bcr_flag = collapse_binary_status(
      bcr_status,
      days_to_bcr
    ),

    .groups = "drop"
  ) |>
  dplyr::mutate(
    recurrence_flag = dplyr::case_when(
      is.finite(days_to_new_tumor) ~ 1L,
      is.finite(days_to_last_followup) ~ 0L,
      TRUE ~ NA_integer_
    ),

    OS_days = dplyr::case_when(
      death_flag == 1L &
        is.finite(days_to_death) ~ days_to_death,

      death_flag == 0L &
        is.finite(days_to_last_followup) ~ days_to_last_followup,

      TRUE ~ NA_real_
    ),

    RFS_days = dplyr::case_when(
      recurrence_flag == 1L &
        is.finite(days_to_new_tumor) ~ days_to_new_tumor,

      recurrence_flag == 0L &
        is.finite(days_to_last_followup) ~ days_to_last_followup,

      TRUE ~ NA_real_
    ),

    BCR_days = dplyr::case_when(
      bcr_flag == 1L &
        is.finite(days_to_bcr) ~ days_to_bcr,

      bcr_flag == 0L &
        is.finite(days_to_last_followup) ~ days_to_last_followup,

      TRUE ~ NA_real_
    ),

    OS_months = OS_days / 30.4375,
    RFS_months = RFS_days / 30.4375,
    BCR_months = BCR_days / 30.4375
  )

ambiguous_bcr_patients <- sum(
  is.na(clin_endpoints$bcr_flag),
  na.rm = TRUE
)

bcr_event_without_time <- sum(
  clin_endpoints$bcr_flag == 1L &
    !is.finite(clin_endpoints$BCR_months),
  na.rm = TRUE
)

################################################################################
# 8. TCGA-CDR PFI endpoint
################################################################################

msg("Loading TCGA-CDR PFI endpoint")

cdr_raw <- readxl::read_xlsx(cdr_fp)

cdr <- as.data.frame(
  cdr_raw,
  stringsAsFactors = FALSE
)

colnames(cdr) <- make.names(
  colnames(cdr),
  unique = TRUE
)

cdr_patient_col <- pick_first(
  colnames(cdr),
  c(
    "bcr_patient_barcode",
    "patient",
    "patient_id",
    "PatientID"
  )
)

if (is.na(cdr_patient_col)) {
  cdr_patient_col <- pick_first_grep(
    colnames(cdr),
    c(
      "bcr.*patient.*barcode",
      "patient"
    )
  )
}

pfi_time_col <- pick_first(
  colnames(cdr),
  c(
    "PFI.time",
    "PFI_time",
    "PFI.time.days",
    "PFI_time_days",
    "PFI.time.months",
    "PFI_time_months"
  )
)

if (is.na(pfi_time_col)) {
  pfi_time_col <- pick_first_grep(
    colnames(cdr),
    c("^PFI.*time")
  )
}

pfi_event_col <- pick_first(
  colnames(cdr),
  c(
    "PFI",
    "PFI_event",
    "PFI.status",
    "PFI_status"
  )
)

if (is.na(pfi_event_col)) {
  pfi_event_col <- pick_first_grep(
    colnames(cdr),
    c(
      "^PFI$",
      "PFI.*event",
      "PFI.*status"
    )
  )
}

if (is.na(cdr_patient_col)) {
  fail("No patient barcode column was found in TCGA-CDR.")
}

if (is.na(pfi_time_col)) {
  fail("No PFI time column was found in TCGA-CDR.")
}

if (is.na(pfi_event_col)) {
  fail("No PFI event column was found in TCGA-CDR.")
}

pfi_time_raw <- as_numeric_clean(
  cdr[[pfi_time_col]]
)

pfi_time_days <- if (
  grepl(
    "month",
    pfi_time_col,
    ignore.case = TRUE
  )
) {
  pfi_time_raw * 30.4375
} else {
  pfi_time_raw
}

pfi_rows <- tibble::tibble(
  patient_id = normalize_patient(
    cdr[[cdr_patient_col]]
  ),

  PFI_time_days = pfi_time_days,
  PFI_event = coerce_binary_event(
    cdr[[pfi_event_col]]
  )
) |>
  dplyr::filter(
    !is.na(patient_id)
  )

pfi_conflicts <- pfi_rows |>
  dplyr::group_by(patient_id) |>
  dplyr::summarise(
    distinct_times = dplyr::n_distinct(
      PFI_time_days[is.finite(PFI_time_days)]
    ),

    distinct_events = dplyr::n_distinct(
      PFI_event[!is.na(PFI_event)]
    ),

    .groups = "drop"
  ) |>
  dplyr::filter(
    distinct_times > 1L |
      distinct_events > 1L
  )

if (nrow(pfi_conflicts) > 0L) {
  fail(
    "TCGA-CDR contains conflicting duplicate PFI records for ",
    nrow(pfi_conflicts),
    " patient(s)."
  )
}

pfi_processed <- pfi_rows |>
  dplyr::group_by(patient_id) |>
  dplyr::summarise(
    PFI_time_days = first_finite(PFI_time_days),
    PFI_event = first_nonmissing_integer(PFI_event),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    PFI_time_months = PFI_time_days / 30.4375
  )

################################################################################
# 9. Canonical clinical master table
################################################################################

msg("Building canonical clinical master table")

clinical_master <- score_df |>
  dplyr::left_join(
    se_covariates,
    by = c(
      "sample_id",
      "patient_id"
    )
  ) |>
  dplyr::left_join(
    clin_endpoints,
    by = "patient_id"
  ) |>
  dplyr::left_join(
    pfi_processed,
    by = "patient_id"
  ) |>
  dplyr::mutate(
    platelet_score_median_group = factor(
      ifelse(
        Score_z <= canonical_score_median,
        "Low",
        "High"
      ),
      levels = c(
        "Low",
        "High"
      )
    )
  )

if (
  nrow(clinical_master) !=
  EXPECTED_CANONICAL_SAMPLES
) {
  fail(
    "Clinical master table contains ",
    nrow(clinical_master),
    " rows; expected ",
    EXPECTED_CANONICAL_SAMPLES,
    "."
  )
}

if (
  anyDuplicated(clinical_master$sample_id) ||
  anyDuplicated(clinical_master$patient_id)
) {
  fail(
    "Clinical master table contains duplicated canonical samples or patients."
  )
}

write_table(
  clinical_master,
  output_paths$master_csv
)

write_rds(
  clinical_master,
  output_paths$master_rds
)

availability_variables <- c(
  "age",
  "PSA",
  "gleason_primary",
  "gleason_secondary",
  "gleason_sum",
  "gleason_group",
  "T",
  "T_bin",
  "N",
  "N_bin",
  "LN_pos",
  "PFI_time_months",
  "PFI_event",
  "BCR_months",
  "bcr_flag",
  "RFS_months",
  "recurrence_flag",
  "OS_months",
  "death_flag"
)

clinical_availability <- tibble::tibble(
  variable = availability_variables,

  n_nonmissing = vapply(
    clinical_master[, availability_variables, drop = FALSE],
    function(x) sum(!is.na(x)),
    integer(1)
  ),

  n_missing = vapply(
    clinical_master[, availability_variables, drop = FALSE],
    function(x) sum(is.na(x)),
    integer(1)
  )
)

write_table(
  clinical_availability,
  output_paths$clinical_availability
)

################################################################################
# 10. Clinical association tests
################################################################################

msg("Running clinical association tests")

spearman_test <- function(df, variable, endpoint_label) {
  analysis_df <- df |>
    dplyr::transmute(
      Score_z = as_numeric_clean(Score_z),
      value = as_numeric_clean(.data[[variable]])
    ) |>
    dplyr::filter(
      is.finite(Score_z),
      is.finite(value)
    )

  if (nrow(analysis_df) < 10L) {
    return(
      tibble::tibble(
        endpoint = endpoint_label,
        test = "Spearman",
        n = nrow(analysis_df),
        estimate = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  test_result <- suppressWarnings(
    stats::cor.test(
      analysis_df$Score_z,
      analysis_df$value,
      method = "spearman",
      exact = FALSE
    )
  )

  tibble::tibble(
    endpoint = endpoint_label,
    test = "Spearman",
    n = nrow(analysis_df),
    estimate = unname(test_result$estimate),
    statistic = unname(test_result$statistic),
    p_value = test_result$p.value
  )
}

kruskal_test <- function(df, variable, endpoint_label) {
  analysis_df <- df |>
    dplyr::transmute(
      Score_z = as_numeric_clean(Score_z),
      group = as.character(.data[[variable]])
    ) |>
    dplyr::filter(
      is.finite(Score_z),
      !is.na(group),
      group != ""
    )

  if (
    nrow(analysis_df) < 10L ||
    dplyr::n_distinct(analysis_df$group) < 2L
  ) {
    return(
      tibble::tibble(
        endpoint = endpoint_label,
        test = "Kruskal-Wallis",
        n = nrow(analysis_df),
        estimate = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  test_result <- stats::kruskal.test(
    Score_z ~ group,
    data = analysis_df
  )

  tibble::tibble(
    endpoint = endpoint_label,
    test = "Kruskal-Wallis",
    n = nrow(analysis_df),
    estimate = NA_real_,
    statistic = unname(test_result$statistic),
    p_value = test_result$p.value
  )
}

wilcoxon_test <- function(df, variable, endpoint_label) {
  analysis_df <- df |>
    dplyr::transmute(
      Score_z = as_numeric_clean(Score_z),
      group = factor(.data[[variable]])
    ) |>
    dplyr::filter(
      is.finite(Score_z),
      !is.na(group)
    ) |>
    droplevels()

  if (
    nrow(analysis_df) < 10L ||
    nlevels(analysis_df$group) != 2L
  ) {
    return(
      tibble::tibble(
        endpoint = endpoint_label,
        test = "Wilcoxon rank-sum",
        n = nrow(analysis_df),
        estimate = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  test_result <- suppressWarnings(
    stats::wilcox.test(
      Score_z ~ group,
      data = analysis_df,
      exact = FALSE
    )
  )

  tibble::tibble(
    endpoint = endpoint_label,
    test = "Wilcoxon rank-sum",
    n = nrow(analysis_df),
    estimate = NA_real_,
    statistic = unname(test_result$statistic),
    p_value = test_result$p.value
  )
}

clinical_associations <- dplyr::bind_rows(
  spearman_test(
    clinical_master,
    "PSA",
    "PSA"
  ),

  spearman_test(
    clinical_master,
    "gleason_sum",
    "Gleason sum"
  ),

  spearman_test(
    clinical_master,
    "LN_pos",
    "Positive lymph nodes"
  ),

  kruskal_test(
    clinical_master,
    "T",
    "Pathologic T stage"
  ),

  wilcoxon_test(
    clinical_master,
    "T_bin",
    "T2 versus T3+"
  ),

  wilcoxon_test(
    clinical_master,
    "N_bin",
    "N0 versus N+"
  )
) |>
  dplyr::mutate(
    FDR = stats::p.adjust(
      p_value,
      method = "BH"
    ),

    p_label = vapply(
      p_value,
      format_p,
      character(1)
    ),

    FDR_label = vapply(
      FDR,
      format_p,
      character(1)
    )
  )

write_table(
  clinical_associations,
  output_paths$clinical_associations
)

################################################################################
# 11. TCGA Gleason figure
################################################################################

msg("Generating platelet score versus Gleason figure")

gleason_levels_figure <- c(
  "6",
  "3+4",
  "4+3",
  "8",
  ">=9"
)

gleason_df <- clinical_master |>
  dplyr::filter(
    is.finite(Score_z),
    as.character(gleason_group) %in%
      gleason_levels_figure
  ) |>
  dplyr::mutate(
    gleason_group = factor(
      as.character(gleason_group),
      levels = gleason_levels_figure,
      ordered = TRUE
    ),

    gleason_rank = as.numeric(
      gleason_group
    )
  )

if (
  nrow(gleason_df) < 10L ||
  dplyr::n_distinct(
    gleason_df$gleason_group
  ) < 2L
) {
  fail(
    "Insufficient structured Gleason data for the required TCGA Gleason figure."
  )
}

gleason_test <- suppressWarnings(
  stats::cor.test(
    gleason_df$Score_z,
    gleason_df$gleason_rank,
    method = "spearman",
    exact = FALSE
  )
)

gleason_rho <- unname(
  gleason_test$estimate
)

gleason_p <- gleason_test$p.value

gleason_summary <- gleason_df |>
  dplyr::group_by(
    gleason_group,
    gleason_rank
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_Score_z = stats::median(
      Score_z,
      na.rm = TRUE
    ),
    q1_Score_z = stats::quantile(
      Score_z,
      0.25,
      na.rm = TRUE
    ),
    q3_Score_z = stats::quantile(
      Score_z,
      0.75,
      na.rm = TRUE
    ),
    IQR_Score_z = q3_Score_z - q1_Score_z,
    .groups = "drop"
  )

gleason_counts <- gleason_summary |>
  dplyr::transmute(
    gleason_group = as.character(
      gleason_group
    ),
    n = n
  )

write_table(
  gleason_counts,
  output_paths$gleason_counts
)

write_table(
  gleason_summary,
  output_paths$gleason_summary
)

write_table(
  tibble::tibble(
    n = nrow(gleason_df),
    rho = gleason_rho,
    p_value = gleason_p,
    p_label = format_p(gleason_p)
  ),
  output_paths$gleason_spearman
)

gleason_x_labels <- setNames(
  paste0(
    as.character(gleason_summary$gleason_group),
    "\n(n=",
    gleason_summary$n,
    ")"
  ),
  gleason_summary$gleason_rank
)

gleason_y_quantiles <- stats::quantile(
  gleason_df$Score_z,
  c(0.01, 0.99),
  na.rm = TRUE
)

gleason_y_padding <- diff(
  gleason_y_quantiles
) * 0.35

gleason_y_limits <- c(
  gleason_y_quantiles[1L] - gleason_y_padding,
  gleason_y_quantiles[2L] + gleason_y_padding
)

gleason_x_min <- min(
  gleason_summary$gleason_rank
)

gleason_x_max <- max(
  gleason_summary$gleason_rank
)
gleason_p_label <- paste0(
  "italic(P) == '",
  format_p(gleason_p),
  "'"
)

gleason_plot <- ggplot2::ggplot(
  gleason_df,
  ggplot2::aes(
    x = gleason_rank,
    y = Score_z
  )
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.12,
      height = 0
    ),
    size = 0.45,
    alpha = 0.35,
    shape = 16,
    color = COL_GLEASON
  ) +
  ggplot2::geom_linerange(
    data = gleason_summary,
    ggplot2::aes(
      x = gleason_rank,
      ymin = q1_Score_z,
      ymax = q3_Score_z
    ),
    inherit.aes = FALSE,
    linewidth = 0.15,
    color = "black"
  ) +
  ggplot2::geom_point(
    data = gleason_summary,
    ggplot2::aes(
      x = gleason_rank,
      y = median_Score_z
    ),
    inherit.aes = FALSE,
    size = 0.60,
    shape = 16,
    color = "black"
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.20,
    color = COL_GLEASON_LINE
  ) +
  ggplot2::scale_x_continuous(
    breaks = gleason_summary$gleason_rank,
    labels = gleason_x_labels,
    limits = c(
      gleason_x_min - 0.35,
      gleason_x_max + 0.35
    ),
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(
    ylim = gleason_y_limits
  ) +
  ggplot2::annotate(
    "text",
    x = gleason_x_min - 0.22,
    y = gleason_y_limits[2L] -
      0.018 * diff(gleason_y_limits),
    label = paste0(
      "Spearman rho = ",
      sprintf("%.3f", gleason_rho)
    ),
    hjust = 0,
    vjust = 1,
    size = STAT_TEXT_PT / 2.845,
    family = FONT_FAMILY,
    color = "black"
  ) +
  ggplot2::annotate(
    "text",
    x = gleason_x_min - 0.22,
    y = gleason_y_limits[2L] -
      0.105 * diff(gleason_y_limits),
    label = gleason_p_label,
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = STAT_TEXT_PT / 2.845,
    family = FONT_FAMILY,
    color = "black"
  ) +
  ggplot2::labs(
    x = "Gleason group",
    y = "Platelet score (z)"
  ) +
  ggplot2::theme_bw(
    base_family = FONT_FAMILY,
    base_size = MAIN_TEXT_PT
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),

    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.22
    ),

    axis.text.x = ggplot2::element_text(
      family = FONT_FAMILY,
      color = "black",
      size = MAIN_TEXT_PT,
      lineheight = 0.75,
      margin = ggplot2::margin(t = 1)
    ),

    axis.text.y = ggplot2::element_text(
      family = FONT_FAMILY,
      color = "black",
      size = MAIN_TEXT_PT
    ),

    axis.title = ggplot2::element_text(
      family = FONT_FAMILY,
      color = "black",
      size = MAIN_TEXT_PT
    ),

    axis.ticks = ggplot2::element_line(
      color = "black",
      linewidth = 0.20
    ),

    plot.margin = ggplot2::margin(
      1,
      1,
      2,
      1,
      unit = "mm"
    )
  )

save_pdf(
  gleason_plot,
  output_paths$figure_gleason,
  width_cm = 4.0,
  height_cm = 3.5
)

################################################################################
# 12. TCGA T2 versus T3+ figure
################################################################################

msg("Generating platelet score T2 versus T3+ figure")

t_stage_df <- clinical_master |>
  dplyr::filter(
    is.finite(Score_z),
    T_bin %in% c(
      "T2",
      "T3+"
    )
  ) |>
  dplyr::mutate(
    T_bin = factor(
      T_bin,
      levels = c(
        "T2",
        "T3+"
      )
    )
  )

if (
  nrow(t_stage_df) < 10L ||
  nlevels(droplevels(t_stage_df$T_bin)) != 2L
) {
  fail(
    "Insufficient T2 versus T3+ data for the required TCGA T-stage figure."
  )
}

t_stage_test <- suppressWarnings(
  stats::wilcox.test(
    Score_z ~ T_bin,
    data = t_stage_df,
    exact = FALSE
  )
)

t_stage_p <- t_stage_test$p.value

write_table(
  tibble::tibble(
    n = nrow(t_stage_df),

    n_T2 = sum(
      t_stage_df$T_bin == "T2"
    ),

    n_T3plus = sum(
      t_stage_df$T_bin == "T3+"
    ),

    statistic = unname(
      t_stage_test$statistic
    ),

    p_value = t_stage_p,
    p_label = format_p(t_stage_p)
  ),
  output_paths$t_stage_wilcoxon
)

t_stage_y_quantiles <- stats::quantile(
  t_stage_df$Score_z,
  c(0.01, 0.99),
  na.rm = TRUE
)

t_stage_y_padding <- diff(
  t_stage_y_quantiles
) * 0.50

t_stage_y_limits <- c(
  t_stage_y_quantiles[1L] - t_stage_y_padding,
  t_stage_y_quantiles[2L] + t_stage_y_padding
)

t_stage_colors <- c(
  "T2" = COL_LOW,
  "T3+" = COL_HIGH
)

t_stage_plot <- ggplot2::ggplot(
  t_stage_df,
  ggplot2::aes(
    x = T_bin,
    y = Score_z
  )
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(
      color = T_bin,
      fill = T_bin
    ),
    width = 0.60,
    linewidth = 0.20,
    outlier.shape = NA,
    alpha = 0.35
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(
      color = T_bin
    ),
    width = 0.08,
    height = 0,
    size = 0.45,
    alpha = 0.30,
    shape = 16
  ) +
  ggplot2::scale_color_manual(
    values = t_stage_colors,
    guide = "none"
  ) +
  ggplot2::scale_fill_manual(
    values = t_stage_colors,
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    ylim = t_stage_y_limits
  ) +
  ggplot2::annotate(
    "text",
    x = 1.03,
    y = t_stage_y_limits[2L] -
      0.03 * diff(t_stage_y_limits),
    label = format_p_label(
      t_stage_p
    ),
    hjust = 0,
    vjust = 1,
    size = STAT_TEXT_PT / 2.845,
    family = FONT_FAMILY,
    color = "black"
  ) +
  ggplot2::labs(
    x = "Pathologic\nT stage",
    y = "Platelet score (z)"
  ) +
  ggplot2::theme_bw(
    base_family = FONT_FAMILY,
    base_size = MAIN_TEXT_PT
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),

    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.22
    ),

    axis.text = ggplot2::element_text(
      family = FONT_FAMILY,
      color = "black",
      size = MAIN_TEXT_PT
    ),

    axis.title = ggplot2::element_text(
      family = FONT_FAMILY,
      color = "black",
      size = MAIN_TEXT_PT
    ),

    axis.ticks = ggplot2::element_line(
      color = "black",
      linewidth = 0.20
    ),

    plot.margin = ggplot2::margin(
      2,
      2,
      2,
      2,
      unit = "mm"
    )
  )

save_pdf(
  t_stage_plot,
  output_paths$figure_t_stage,
  width_cm = 2.8,
  height_cm = 3.5
)

################################################################################
# 13. Survival helpers
################################################################################

make_endpoint_data <- function(
    df,
    time_column,
    event_column
) {
  endpoint_df <- df |>
    dplyr::transmute(
      patient_id = patient_id,
      Score_z = as_numeric_clean(Score_z),
      age = as_numeric_clean(age),
      gleason_sum = as_numeric_clean(gleason_sum),
      T = as.character(T),

      time = as_numeric_clean(
        .data[[time_column]]
      ),

      event_raw = as_numeric_clean(
        .data[[event_column]]
      )
    ) |>
    dplyr::mutate(
      event = dplyr::case_when(
        event_raw == 0 ~ 0L,
        event_raw > 0 ~ 1L,
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::filter(
      is.finite(Score_z),
      is.finite(time),
      time > 0,
      !is.na(event)
    )

  if (nrow(endpoint_df) > 0L) {
    endpoint_median <- stats::median(
      endpoint_df$Score_z
    )

    endpoint_df <- endpoint_df |>
      dplyr::mutate(
        median_cutpoint = endpoint_median,

        score_group = factor(
          ifelse(
            Score_z <= endpoint_median,
            "Low",
            "High"
          ),
          levels = c(
            "Low",
            "High"
          )
        ),

        age_per10years = age / 10,

        T_factor = factor(
          T,
          levels = c(
            "T2",
            "T3",
            "T4"
          )
        )
      )
  }

  endpoint_df
}

logrank_p_value <- function(endpoint_df) {
  survival_difference <- survival::survdiff(
    survival::Surv(time, event) ~ score_group,
    data = endpoint_df
  )

  stats::pchisq(
    survival_difference$chisq,
    df = length(survival_difference$n) - 1L,
    lower.tail = FALSE
  )
}

tidy_cox_model <- function(
    fit,
    endpoint,
    model_label,
    n,
    events
) {
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) |>
    dplyr::mutate(
      endpoint = endpoint,
      model = model_label,
      n = n,
      events = events,
      .before = 1L
    )
}

tidy_cox_ph <- function(
    fit,
    endpoint,
    model_label
) {
  ph_test <- survival::cox.zph(
    fit
  )

  ph_table <- as.data.frame(
    ph_test$table,
    stringsAsFactors = FALSE
  )

  ph_table$term <- rownames(
    ph_table
  )

  rownames(ph_table) <- NULL

  tibble::tibble(
    endpoint = endpoint,
    model = model_label,
    term = ph_table$term,
    chisq = as_numeric_clean(
      ph_table$chisq
    ),
    df = as_numeric_clean(
      ph_table$df
    ),
    p_value = as_numeric_clean(
      ph_table$p
    )
  )
}

make_km_plot <- function(
    fit,
    endpoint_df,
    endpoint,
    binary_cox_fit,
    output_path
) {
  palette <- c(
    "Low" = COL_LOW,
    "High" = COL_HIGH
  )

  logrank_p <- logrank_p_value(
    endpoint_df
  )

  binary_cox_table <- broom::tidy(
    binary_cox_fit,
    exponentiate = TRUE,
    conf.int = TRUE
  )

  if (nrow(binary_cox_table) != 1L) {
    fail(
      endpoint,
      ": binary median-split Cox model did not return one coefficient."
    )
  }

  hr_text <- paste0(
    "HR = ",
    sprintf(
      "%.2f",
      binary_cox_table$estimate[1L]
    ),
    " (95% CI ",
    sprintf(
      "%.2f",
      binary_cox_table$conf.low[1L]
    ),
    "-",
    sprintf(
      "%.2f",
      binary_cox_table$conf.high[1L]
    ),
    ")"
  )

  logrank_text <- paste0(
    "'Log-rank'~italic(P) == '",
    format_p(logrank_p),
    "'"
  )

  y_axis_labels <- c(
    BCR = "Biochemical recurrence-free survival (%)",
    PFI = "Progression-free interval (%)",
    RFS = "Recurrence-free survival (%)",
    OS = "Overall survival (%)"
  )

  y_axis_label <- y_axis_labels[[endpoint]]

  if (is.null(y_axis_label)) {
    y_axis_label <- "Event-free survival (%)"
  }

  km_object <- survminer::ggsurvplot(
    fit,
    data = endpoint_df,
    risk.table = FALSE,
    conf.int = FALSE,
    palette = unname(palette),
    legend = "none",
    pval = FALSE,
    size = KM_LINE_LW,
    censor.size = KM_CENSOR_SIZE,
    censor.shape = KM_CENSOR_SHAPE,
    break.time.by = KM_BREAK_BY,
    xlab = "Months",
    ylab = y_axis_label,
    ggtheme = ggplot2::theme_classic(
      base_family = FONT_FAMILY,
      base_size = MAIN_TEXT_PT
    )
  )


  risk_summary <- summary(
    fit,
    times = KM_TIMES,
    extend = TRUE
  )

  risk_table <- data.frame(
    time = risk_summary$time,
    n_risk = risk_summary$n.risk,
    strata = as.character(
      risk_summary$strata
    ),
    stringsAsFactors = FALSE
  )

  risk_table$strata <- sub(
    "^score_group=",
    "",
    risk_table$strata
  )

  risk_table$y <- ifelse(
    risk_table$strata == "Low",
    -0.225,
    -0.295
  )

  risk_low <- risk_table |>
    dplyr::filter(
      strata == "Low"
    )

  risk_high <- risk_table |>
    dplyr::filter(
      strata == "High"
    )

  n_low <- sum(
    endpoint_df$score_group == "Low"
  )

  n_high <- sum(
    endpoint_df$score_group == "High"
  )

  label_low <- paste0(
    "Platelet score <= median (n=",
    n_low,
    ")"
  )

  label_high <- paste0(
    "Platelet score > median (n=",
    n_high,
    ")"
  )

  x_annotation <- KM_XMAX * 0.020
  y_top <- 0.11
  y_step <- 0.070
  risk_title_y <- -0.155

  km_plot <- km_object$plot +
    ggplot2::coord_cartesian(
      xlim = c(0, KM_XMAX),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::scale_x_continuous(
      breaks = KM_TIMES,
      limits = c(0, KM_XMAX),
      expand = c(0, 0),
      oob = scales::squish
    )  +
    ggplot2::scale_y_continuous(
      breaks = c(
        0,
        0.25,
        0.50,
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
    ggplot2::theme(
      text = ggplot2::element_text(
        family = FONT_FAMILY,
        color = "black"
      ),

      panel.grid = ggplot2::element_blank(),

      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = PANEL_BORDER_LW
      ),

      axis.line = ggplot2::element_blank(),

      axis.text = ggplot2::element_text(
        size = AXIS_TEXT_PT,
        color = "black"
      ),

      axis.title = ggplot2::element_text(
        size = AXIS_TITLE_PT,
        color = "black"
      ),

      axis.ticks = ggplot2::element_line(
        linewidth = AXIS_TICK_LW,
        color = "black"
      ),

      axis.ticks.length = grid::unit(
        0.5,
        "mm"
      ),

      plot.margin = ggplot2::margin(
        2.2,
        1.8,
        8.0,
        1.8,
        unit = "mm"
      )
    ) +
    ggplot2::annotate(
      "text",
      x = x_annotation,
      y = y_top + 2 * y_step,
      label = label_low,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      family = FONT_FAMILY,
      color = COL_LOW
    ) +
    ggplot2::annotate(
      "text",
      x = x_annotation,
      y = y_top + y_step,
      label = label_high,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      family = FONT_FAMILY,
      color = COL_HIGH
    ) +
    ggplot2::annotate(
      "text",
      x = x_annotation,
      y = y_top,
      label = hr_text,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      family = FONT_FAMILY,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = x_annotation,
      y = y_top - 1.10 * y_step,
      label = logrank_text,
      parse = TRUE,
      hjust = 0,
      vjust = 0,
      size = 4.5 / 2.845,
      family = FONT_FAMILY,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = KM_XMAX * 0.5,
      y = risk_title_y,
      label = "Number at risk",
      hjust = 0.5,
      vjust = 0.5,
      size = 4.4 / 2.845,
      family = FONT_FAMILY,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = risk_low,
      ggplot2::aes(
        x = time,
        y = y,
        label = n_risk
      ),
      inherit.aes = FALSE,
      size = 4.5 / 2.845,
      family = FONT_FAMILY,
      color = COL_LOW,
      vjust = 0.5
    ) +
    ggplot2::geom_text(
      data = risk_high,
      ggplot2::aes(
        x = time,
        y = y,
        label = n_risk
      ),
      inherit.aes = FALSE,
      size = 4.5 / 2.845,
      family = FONT_FAMILY,
      color = COL_HIGH,
      vjust = 0.5
    )

  save_pdf(
    km_plot,
    output_path,
    width_cm = 4.2,
    height_cm = 4.6
  )

  invisible(
    list(
      plot = km_plot,
      logrank_p = logrank_p,
      hr_table = binary_cox_table,
      risk_table = risk_table
    )
  )
}

run_survival_endpoint <- function(
    endpoint,
    time_column,
    event_column,
    required = FALSE,
    km_output = NULL
) {
  endpoint_data <- make_endpoint_data(
    clinical_master,
    time_column,
    event_column
  )

  n_endpoint <- nrow(
    endpoint_data
  )

  endpoint_events <- sum(
    endpoint_data$event == 1L,
    na.rm = TRUE
  )

  cohort_summary <- tibble::tibble(
    endpoint = endpoint,
    n_total_canonical = nrow(clinical_master),
    n_analyzable = n_endpoint,
    n_excluded = nrow(clinical_master) - n_endpoint,
    events = endpoint_events,
    censored = n_endpoint - endpoint_events
  )

  if (
    n_endpoint < 10L ||
    endpoint_events < 1L ||
    dplyr::n_distinct(
      endpoint_data$score_group
    ) < 2L
  ) {
    error_text <- paste0(
      endpoint,
      ": insufficient analyzable survival data. n=",
      n_endpoint,
      "; events=",
      endpoint_events,
      "."
    )

    if (required) {
      fail(error_text)
    }

    add_warning(error_text)

    return(
      list(
        km = NULL,
        cox = NULL,
        ph = NULL,
        cohort = cohort_summary
      )
    )
  }

  if (
    endpoint == "OS" &&
    endpoint_events < 30L
  ) {
    add_warning(
      "OS has only ",
      endpoint_events,
      " events; OS estimates are exploratory and unstable."
    )
  }

  km_fit <- survival::survfit(
    survival::Surv(time, event) ~ score_group,
    data = endpoint_data
  )

  logrank_p <- logrank_p_value(
    endpoint_data
  )

  km_summary <- tibble::tibble(
    endpoint = endpoint,
    analysis = "KM_median_split",
    model = "logrank",
    n = n_endpoint,
    events = endpoint_events,
    term = "High versus Low",
    estimate = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    p_value = logrank_p,
    p_label = format_p(logrank_p),
    median_cutpoint = unique(
      endpoint_data$median_cutpoint
    )[1L],
    n_low = sum(
      endpoint_data$score_group == "Low"
    ),
    n_high = sum(
      endpoint_data$score_group == "High"
    )
  )

  cox_tables <- list()
  ph_tables <- list()

  if (endpoint_events >= 5L) {
    binary_cox <- survival::coxph(
      survival::Surv(time, event) ~ score_group,
      data = endpoint_data,
      x = TRUE
    )

    continuous_cox <- survival::coxph(
      survival::Surv(time, event) ~ Score_z,
      data = endpoint_data,
      x = TRUE
    )

    cox_tables[["binary"]] <- tidy_cox_model(
      binary_cox,
      endpoint,
      "binary_median",
      n_endpoint,
      endpoint_events
    )

    cox_tables[["continuous"]] <- tidy_cox_model(
      continuous_cox,
      endpoint,
      "continuous_canonical_Score_z",
      n_endpoint,
      endpoint_events
    )

    ph_tables[["binary"]] <- tidy_cox_ph(
      binary_cox,
      endpoint,
      "binary_median"
    )

    ph_tables[["continuous"]] <- tidy_cox_ph(
      continuous_cox,
      endpoint,
      "continuous_canonical_Score_z"
    )

    if (!is.null(km_output)) {
      make_km_plot(
        fit = km_fit,
        endpoint_df = endpoint_data,
        endpoint = endpoint,
        binary_cox_fit = binary_cox,
        output_path = km_output
      )
    }
  } else {
    add_warning(
      endpoint,
      ": Cox models were skipped because only ",
      endpoint_events,
      " events were available."
    )
  }

  multivariable_data <- endpoint_data |>
    dplyr::filter(
      is.finite(Score_z),
      is.finite(age_per10years),
      is.finite(gleason_sum),
      !is.na(T_factor)
    ) |>
    droplevels()

  multivariable_events <- sum(
    multivariable_data$event == 1L,
    na.rm = TRUE
  )

  t_levels <- dplyr::n_distinct(
    multivariable_data$T_factor
  )

  if (
    nrow(multivariable_data) >=
    MIN_MULTIVARIABLE_COMPLETE_CASES &&
    multivariable_events >=
    MIN_MULTIVARIABLE_EVENTS &&
    t_levels >= 2L
  ) {
    multivariable_cox <- survival::coxph(
      survival::Surv(time, event) ~
        Score_z +
        age_per10years +
        gleason_sum +
        T_factor,
      data = multivariable_data,
      x = TRUE
    )

    cox_tables[["multivariable"]] <- tidy_cox_model(
      multivariable_cox,
      endpoint,
      "multivariable",
      nrow(multivariable_data),
      multivariable_events
    )

    ph_tables[["multivariable"]] <- tidy_cox_ph(
      multivariable_cox,
      endpoint,
      "multivariable"
    )
  } else {
    add_warning(
      endpoint,
      ": multivariable Cox model skipped. Complete cases=",
      nrow(multivariable_data),
      "; events=",
      multivariable_events,
      "; observed T-stage levels=",
      t_levels,
      "."
    )
  }

  list(
    km = km_summary,
    cox = dplyr::bind_rows(
      cox_tables
    ),
    ph = dplyr::bind_rows(
      ph_tables
    ),
    cohort = cohort_summary
  )
}

################################################################################
# 14. Survival analyses
################################################################################

msg("Running PFI survival analysis")

survival_results <- list()

survival_results[["PFI"]] <- run_survival_endpoint(
  endpoint = "PFI",
  time_column = "PFI_time_months",
  event_column = "PFI_event",
  required = FALSE,
  km_output = NULL
)

msg("Running principal BCR survival analysis")

survival_results[["BCR"]] <- run_survival_endpoint(
  endpoint = "BCR",
  time_column = "BCR_months",
  event_column = "bcr_flag",
  required = TRUE,
  km_output = output_paths$figure_bcr_km
)

msg("Running complementary RFS survival analysis")

survival_results[["RFS"]] <- run_survival_endpoint(
  endpoint = "RFS",
  time_column = "RFS_months",
  event_column = "recurrence_flag",
  required = FALSE,
  km_output = NULL
)

msg("Running exploratory OS survival analysis")

survival_results[["OS"]] <- run_survival_endpoint(
  endpoint = "OS",
  time_column = "OS_months",
  event_column = "death_flag",
  required = FALSE,
  km_output = NULL
)

endpoint_counts <- dplyr::bind_rows(
  lapply(
    survival_results,
    function(x) x$cohort
  )
)

write_table(
  endpoint_counts,
  output_paths$endpoint_counts
)

km_summary_all <- dplyr::bind_rows(
  lapply(
    survival_results,
    function(x) x$km
  )
)

cox_summary_all <- dplyr::bind_rows(
  lapply(
    survival_results,
    function(x) x$cox
  )
)

cox_ph_all <- dplyr::bind_rows(
  lapply(
    survival_results,
    function(x) x$ph
  )
)

if (nrow(cox_ph_all) == 0L) {
  cox_ph_all <- tibble::tibble(
    endpoint = character(),
    model = character(),
    term = character(),
    chisq = numeric(),
    df = numeric(),
    p_value = numeric()
  )
}

cox_ph_all <- cox_ph_all |>
  dplyr::mutate(
    p_label = vapply(
      p_value,
      format_p,
      character(1)
    )
  )

write_table(
  cox_ph_all,
  output_paths$cox_ph_diagnostics
)

if (nrow(cox_summary_all) > 0L) {
  cox_summary_standardized <- cox_summary_all |>
    dplyr::transmute(
      endpoint = endpoint,

      analysis = dplyr::case_when(
        model == "binary_median" ~
          "Cox_median_split",

        model ==
          "continuous_canonical_Score_z" ~
          "Cox_continuous_score",

        model == "multivariable" ~
          "Cox_multivariable",

        TRUE ~ model
      ),

      model = model,
      n = n,
      events = events,
      term = term,
      estimate = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p_value = p.value,

      p_label = vapply(
        p.value,
        format_p,
        character(1)
      ),

      median_cutpoint = NA_real_,
      n_low = NA_integer_,
      n_high = NA_integer_
    )
} else {
  cox_summary_standardized <- tibble::tibble(
    endpoint = character(),
    analysis = character(),
    model = character(),
    n = integer(),
    events = integer(),
    term = character(),
    estimate = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    p_value = numeric(),
    p_label = character(),
    median_cutpoint = numeric(),
    n_low = integer(),
    n_high = integer()
  )
}

survival_summary_all <- dplyr::bind_rows(
  km_summary_all,
  cox_summary_standardized
)

write_table(
  survival_summary_all,
  output_paths$survival_summary
)

################################################################################
# 15. Final validation
################################################################################

required_outputs <- c(
  output_paths$master_csv,
  output_paths$master_rds,
  output_paths$clinical_associations,
  output_paths$clinical_availability,
  output_paths$gleason_counts,
  output_paths$gleason_summary,
  output_paths$gleason_spearman,
  output_paths$t_stage_wilcoxon,
  output_paths$endpoint_counts,
  output_paths$survival_summary,
  output_paths$cox_ph_diagnostics,
  output_paths$figure_gleason,
  output_paths$figure_t_stage,
  output_paths$figure_bcr_km
)

missing_outputs <- required_outputs[
  !file.exists(required_outputs)
]

if (length(missing_outputs) > 0L) {
  fail(
    "Required output(s) were not generated: ",
    paste(
      missing_outputs,
      collapse = " | "
    )
  )
}

empty_outputs <- required_outputs[
  file.info(required_outputs)$size <= 0
]

if (length(empty_outputs) > 0L) {
  fail(
    "Generated output(s) are empty: ",
    paste(
      empty_outputs,
      collapse = " | "
    )
  )
}

bcr_result <- survival_results[["BCR"]]

if (
  is.null(bcr_result$km) ||
  nrow(bcr_result$km) != 1L
) {
  fail(
    "Principal BCR Kaplan–Meier summary was not generated."
  )
}

bcr_continuous_result <- cox_summary_all |>
  dplyr::filter(
    endpoint == "BCR",
    model == "continuous_canonical_Score_z",
    term == "Score_z"
  )

if (nrow(bcr_continuous_result) != 1L) {
  fail(
    "Principal continuous BCR Cox result for canonical Score_z ",
    "was not generated."
  )
}

final_status <- if (
  length(warnings_log) == 0L
) {
  "PASS"
} else {
  "PASS_WITH_WARNINGS"
}

################################################################################
# 16. QC log
################################################################################

qc_lines <- c(
  "TCGA-PRAD clinical analysis using the canonical 41-gene score",
  "================================================================",
  "",
  paste("Date/time:", as.character(Sys.time())),
  paste("Final status:", final_status),
  paste("Project directory:", project_dir),
  paste("Script:", script_fp),
  "",
  "ANALYTICAL CONTRACT",
  "-------------------",
  "Clinical interpretation: association analysis; not biomarker validation",
  "Canonical cohort: 497 primary tumors; one sample per patient",
  "Median split: Kaplan-Meier visualization only",
  "Continuous Cox variable: canonical Score_z",
  "Endpoint-specific Score_z re-standardization: not performed",
  "Principal TCGA survival endpoint: BCR",
  "Complementary endpoints: PFI and RFS",
  "Exploratory endpoint: OS",
  "",
  "CANONICAL SCORE COHORT",
  "----------------------",
  paste("Score source:", score_fp),
  paste("Score rows:", nrow(score_df)),
  paste("Unique samples:", dplyr::n_distinct(score_df$sample_id)),
  paste("Unique patients:", dplyr::n_distinct(score_df$patient_id)),
  paste("Score_z mean:", signif(mean(score_df$Score_z), 8)),
  paste("Score_z SD:", signif(stats::sd(score_df$Score_z), 8)),
  paste("Canonical Score_z median:", signif(canonical_score_median, 8)),
  "Duplicate patient resolution inside clinical script: not performed",
  "Canonical sample IDs used directly: yes",
  "",
  "SE CLINICAL ALIGNMENT",
  "---------------------",
  paste("SE samples before alignment:", nrow(SummarizedExperiment::colData(se_tumor))),
  paste("Canonical samples aligned:", nrow(se_covariates)),
  "Alignment method: exact full sample_id",
  "Barcode truncation: not used",
  "",
  "DETECTED SE CLINICAL COLUMNS",
  "----------------------------",
  paste(
    capture.output(
      print(
        detected_columns,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "CLINICAL VARIABLE AVAILABILITY",
  "------------------------------",
  paste(
    capture.output(
      print(
        clinical_availability,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "CLINPLUS ENDPOINT HANDLING",
  "--------------------------",
  paste("Ambiguous BCR status patients:", ambiguous_bcr_patients),
  paste("BCR events without valid event time:", bcr_event_without_time),
  "Ambiguous BCR status converted to no-event: no",
  "BCR event without valid event date included in survival model: no",
  "OS unknown status converted to alive: no",
  "",
  "ENDPOINT ANALYSIS COUNTS",
  "------------------------",
  paste(
    capture.output(
      print(
        endpoint_counts,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "PRINCIPAL BCR RESULTS",
  "---------------------",
  paste(
    capture.output(
      print(
        bcr_result$km,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  paste(
    capture.output(
      print(
        bcr_continuous_result,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "COX PROPORTIONAL-HAZARDS DIAGNOSTICS",
  "------------------------------------",
  paste(
    capture.output(
      print(
        cox_ph_all,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "FROZEN KAPLAN-MEIER STYLE",
  "-------------------------",
  paste("Font:", FONT_FAMILY),
  paste("Low group color:", COL_LOW),
  paste("High group color:", COL_HIGH),
  paste("Curve linewidth:", KM_LINE_LW),
  paste("Censor size:", KM_CENSOR_SIZE),
  paste("Censor shape:", KM_CENSOR_SHAPE),
  paste("X-axis maximum:", KM_XMAX, "months"),
  paste("X-axis ticks:", paste(KM_TIMES, collapse = ", ")),
  paste("Axis text:", AXIS_TEXT_PT, "pt"),
  paste("Axis title:", AXIS_TITLE_PT, "pt"),
  paste("Panel border linewidth:", PANEL_BORDER_LW),
  "Confidence interval: not shown",
  "Manual number-at-risk table: yes",
  "Risk-table group labels: not shown",
  "",
  "REQUIRED FIGURES",
  "----------------",
  paste("Gleason:", output_paths$figure_gleason),
  paste("T2 versus T3+:", output_paths$figure_t_stage),
  paste("BCR Kaplan-Meier:", output_paths$figure_bcr_km),
  "",
  "GENERATED TABLES",
  "----------------",
  if (
    length(tables_generated) == 0L
  ) {
    "None"
  } else {
    paste(" -", tables_generated)
  },
  "",
  "GENERATED FIGURES",
  "-----------------",
  if (
    length(figures_generated) == 0L
  ) {
    "None"
  } else {
    paste(" -", figures_generated)
  },
  "",
  "GENERATED OBJECTS",
  "-----------------",
  if (
    length(objects_generated) == 0L
  ) {
    "None"
  } else {
    paste(" -", objects_generated)
  },
  "",
  "WARNINGS",
  "--------",
  if (
    length(warnings_log) == 0L
  ) {
    "None"
  } else {
    paste(" -", warnings_log)
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
  paste("FINAL STATUS:", final_status)
)

writeLines(
  qc_lines,
  con = output_paths$qc_log,
  useBytes = TRUE
)

msg("QC log written:", output_paths$qc_log)
msg("Final status:", final_status)
