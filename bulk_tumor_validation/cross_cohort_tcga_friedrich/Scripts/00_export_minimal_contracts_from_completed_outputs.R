#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 17)

EXPECTED_N <- c(
  TCGA_PRAD = 497L,
  Friedrich_GSE134051 = 164L
)

COHORT_LABELS <- c(
  TCGA_PRAD = "TCGA-PRAD",
  Friedrich_GSE134051 = "Friedrich/GSE134051"
)

REQUIRED_PATHWAYS <- c(
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

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

get_script_path <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(file_arg) != 1L) {
    fail("Could not determine adapter path from the Rscript invocation.")
  }

  normalizePath(
    sub("^--file=", "", file_arg),
    mustWork = TRUE
  )
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

  expected_script <- file.path(
    root,
    "bulk_tumor_validation",
    "cross_cohort_tcga_friedrich",
    "Scripts",
    basename(script_path)
  )

  if (
    !file.exists(expected_script) ||
    !identical(
      normalizePath(expected_script, mustWork = TRUE),
      script_path
    )
  ) {
    fail("Derived repository root is incompatible with the adapter location.")
  }

  root
}

parse_arguments <- function(args) {
  required <- c(
    "tcga-score",
    "friedrich-score",
    "tcga-gsea",
    "friedrich-gsea",
    "combined-emt"
  )

  force_count <- sum(args == "--force")
  if (force_count > 1L) {
    fail("--force may be supplied at most once.")
  }

  value_args <- args[args != "--force"]
  if (
    length(value_args) != length(required) ||
    any(!grepl("^--[a-z-]+=.+$", value_args))
  ) {
    fail(
      "Required arguments: ",
      paste0("--", required, "=<path>", collapse = " "),
      " [--force]"
    )
  }

  keys <- sub("^--([^=]+)=.*$", "\\1", value_args)
  values <- sub("^--[^=]+=", "", value_args)

  if (anyDuplicated(keys) || !setequal(keys, required)) {
    fail("Arguments are missing, duplicated or unsupported.")
  }

  names(values) <- keys

  missing_inputs <- values[!vapply(values, file.exists, logical(1))]
  if (length(missing_inputs) > 0L) {
    fail("One or more required input files do not exist.")
  }

  values <- vapply(
    values,
    normalizePath,
    character(1),
    mustWork = TRUE
  )

  list(paths = values, force = force_count == 1L)
}

read_csv <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

detect_column <- function(df, candidates, label) {
  found <- candidates[candidates %in% names(df)]

  if (length(found) == 0L) {
    fail("Missing supported ", label, " column.")
  }

  found[[1L]]
}

strict_numeric <- function(x, label) {
  if (length(x) == 0L || anyNA(x) || any(trimws(as.character(x)) == "")) {
    fail(label, " contains missing or empty values.")
  }

  out <- suppressWarnings(as.numeric(x))
  if (anyNA(out) || any(!is.finite(out))) {
    fail(label, " must contain only finite numeric values.")
  }

  out
}

validate_source_ids <- function(ids, cohort_id, expected_n) {
  ids <- trimws(as.character(ids))

  if (
    length(ids) != expected_n ||
    anyNA(ids) ||
    any(ids == "") ||
    anyDuplicated(ids)
  ) {
    fail(cohort_id, " sample IDs are missing, duplicated or have an unexpected count.")
  }

  valid_pattern <- if (identical(cohort_id, "TCGA_PRAD")) {
    grepl("^TCGA-", ids)
  } else {
    grepl("^GSM[0-9]+$", ids)
  }

  if (any(!valid_pattern)) {
    fail(cohort_id, " sample IDs do not match the expected technical-accession pattern.")
  }

  ids
}

validate_sample_contract <- function(df, cohort_id, expected_n) {
  expected_names <- c("cohort_id", "score_z", "emt_ssgsea")

  if (!identical(names(df), expected_names)) {
    fail(cohort_id, " sample_scores contract has an unexpected schema.")
  }

  if (
    nrow(df) != expected_n ||
    anyNA(df$cohort_id) ||
    !identical(unique(as.character(df$cohort_id)), cohort_id)
  ) {
    fail(cohort_id, " sample_scores contract has an unexpected cohort or row count.")
  }

  score_z <- strict_numeric(df$score_z, paste0(cohort_id, " score_z"))
  emt <- strict_numeric(df$emt_ssgsea, paste0(cohort_id, " emt_ssgsea"))

  if (abs(mean(score_z)) > 1e-8 || abs(stats::sd(score_z) - 1) > 1e-8) {
    fail(cohort_id, " score_z is not standardized within tolerance.")
  }

  if (!is.finite(stats::sd(emt)) || stats::sd(emt) <= 0) {
    fail(cohort_id, " emt_ssgsea must have positive finite sample standard deviation.")
  }

  forbidden_header <- grepl(
    "sample|barcode|gsm|rib",
    names(df),
    ignore.case = TRUE
  )

  text_values <- unlist(df, use.names = FALSE)
  forbidden_value <- grepl(
    "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}|GSM[0-9]+|RIB[[:alnum:]_-]*",
    as.character(text_values),
    ignore.case = TRUE
  )

  if (any(forbidden_header) || any(forbidden_value)) {
    fail(cohort_id, " sample_scores contract contains a forbidden identifier field or value.")
  }

  invisible(TRUE)
}

build_sample_contract <- function(
    score,
    combined,
    cohort_id,
    score_sample_candidates,
    expected_n
) {
  score_id_col <- detect_column(
    score,
    score_sample_candidates,
    paste0(cohort_id, " score sample-ID")
  )
  score_z_col <- detect_column(
    score,
    c("Score_z", "score_z", "platelet_score_z"),
    paste0(cohort_id, " score z")
  )

  combined_cohort_col <- detect_column(
    combined,
    c("cohort"),
    "combined EMT cohort"
  )
  combined_id_col <- detect_column(
    combined,
    c(
      "sample_id", "sample_name", "sample", "Sample", "SampleID",
      "sampleID", "barcode", "Barcode", "GSM", "gsm"
    ),
    "combined EMT sample-ID"
  )
  combined_score_col <- detect_column(
    combined,
    c("platelet_score_z", "Score_z", "score_z"),
    "combined EMT platelet score z"
  )
  combined_emt_col <- detect_column(
    combined,
    c(
      "EMT_score", "Hallmark_EMT_ssGSEA",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "emt_ssgsea"
    ),
    "combined EMT score"
  )

  if (nrow(score) != expected_n) {
    fail(cohort_id, " authoritative score table has an unexpected row count.")
  }

  score_ids <- validate_source_ids(score[[score_id_col]], cohort_id, expected_n)
  score_z <- strict_numeric(score[[score_z_col]], paste0(cohort_id, " authoritative Score_z"))

  expected_label <- unname(COHORT_LABELS[[cohort_id]])
  combined_cohort <- trimws(as.character(combined[[combined_cohort_col]]))
  sub <- combined[combined_cohort == expected_label, , drop = FALSE]

  if (nrow(sub) != expected_n) {
    fail(cohort_id, " combined EMT table has an unexpected cohort row count.")
  }

  combined_ids <- validate_source_ids(sub[[combined_id_col]], cohort_id, expected_n)

  if (!setequal(score_ids, combined_ids)) {
    fail(cohort_id, " score and combined EMT sample sets are not identical.")
  }

  idx <- match(score_ids, combined_ids)
  if (anyNA(idx)) {
    fail(cohort_id, " exact sample-ID matching failed.")
  }

  combined_score_z <- strict_numeric(
    sub[[combined_score_col]][idx],
    paste0(cohort_id, " combined-table platelet score z")
  )
  emt <- strict_numeric(
    sub[[combined_emt_col]][idx],
    paste0(cohort_id, " combined-table EMT score")
  )

  max_difference <- max(abs(score_z - combined_score_z))
  if (!is.finite(max_difference) || max_difference > 1e-12) {
    fail(cohort_id, " Score_z differs between authoritative and combined tables.")
  }

  out <- data.frame(
    cohort_id = rep(cohort_id, expected_n),
    score_z = score_z,
    emt_ssgsea = emt,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$score_z, out$emt_ssgsea), , drop = FALSE]
  rownames(out) <- NULL

  validate_sample_contract(out, cohort_id, expected_n)

  list(
    contract = out,
    maximum_score_difference = max_difference
  )
}

validate_hallmark_contract <- function(df, cohort_id) {
  expected_names <- c("cohort_id", "analysis_id", "pathway", "nes", "padj")

  if (!identical(names(df), expected_names)) {
    fail(cohort_id, " hallmark_gsea contract has an unexpected schema.")
  }

  if (
    nrow(df) != 50L ||
    !identical(unique(as.character(df$cohort_id)), cohort_id) ||
    !identical(unique(as.character(df$analysis_id)), "Q1Q4")
  ) {
    fail(cohort_id, " hallmark_gsea contract has an unexpected cohort, analysis or row count.")
  }

  pathways <- trimws(as.character(df$pathway))
  if (
    anyNA(pathways) ||
    any(pathways == "") ||
    any(!grepl("^HALLMARK_[A-Z0-9_]+$", pathways)) ||
    anyDuplicated(pathways) ||
    anyDuplicated(df[c("cohort_id", "analysis_id", "pathway")])
  ) {
    fail(cohort_id, " hallmark_gsea pathways or keys are malformed or duplicated.")
  }

  nes <- strict_numeric(df$nes, paste0(cohort_id, " NES"))
  padj <- strict_numeric(df$padj, paste0(cohort_id, " adjusted P value"))

  if (any(padj < 0 | padj > 1)) {
    fail(cohort_id, " hallmark_gsea padj values must be within [0, 1].")
  }

  if (!all(REQUIRED_PATHWAYS %in% pathways)) {
    fail(cohort_id, " hallmark_gsea is missing one or more required selected pathways.")
  }

  invisible(list(nes = nes, padj = padj, pathways = pathways))
}

build_hallmark_contract <- function(gsea, cohort_id) {
  required_source <- c("pathway", "NES", "padj")
  if (!all(required_source %in% names(gsea))) {
    fail(cohort_id, " GSEA source is missing pathway, NES or padj.")
  }

  out <- data.frame(
    cohort_id = rep(cohort_id, nrow(gsea)),
    analysis_id = rep("Q1Q4", nrow(gsea)),
    pathway = trimws(as.character(gsea$pathway)),
    nes = strict_numeric(gsea$NES, paste0(cohort_id, " source NES")),
    padj = strict_numeric(gsea$padj, paste0(cohort_id, " source padj")),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$pathway), , drop = FALSE]
  rownames(out) <- NULL

  validate_hallmark_contract(out, cohort_id)
  out
}

write_tsv <- function(df, path) {
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

script_path <- get_script_path()
repo_root <- derive_repo_root(script_path)
arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))

output_paths <- c(
  tcga_sample = file.path(
    repo_root, "bulk_tumor_validation", "tcga_prad", "Contracts", "sample_scores.tsv"
  ),
  tcga_gsea = file.path(
    repo_root, "bulk_tumor_validation", "tcga_prad", "Contracts", "hallmark_gsea.tsv"
  ),
  friedrich_sample = file.path(
    repo_root, "bulk_tumor_validation", "friedrich_gse134051", "Contracts", "sample_scores.tsv"
  ),
  friedrich_gsea = file.path(
    repo_root, "bulk_tumor_validation", "friedrich_gse134051", "Contracts", "hallmark_gsea.tsv"
  )
)

existing <- output_paths[file.exists(output_paths)]
if (length(existing) > 0L && !arguments$force) {
  fail("Contract output already exists; rerun with --force to replace authorized TSV files.")
}

score_tcga <- read_csv(arguments$paths[["tcga-score"]])
score_friedrich <- read_csv(arguments$paths[["friedrich-score"]])
gsea_tcga <- read_csv(arguments$paths[["tcga-gsea"]])
gsea_friedrich <- read_csv(arguments$paths[["friedrich-gsea"]])
combined_emt <- read_csv(arguments$paths[["combined-emt"]])

combined_cohort_col <- detect_column(combined_emt, c("cohort"), "combined EMT cohort")
observed_labels <- unique(trimws(as.character(combined_emt[[combined_cohort_col]])))
if (
  nrow(combined_emt) != sum(EXPECTED_N) ||
  !setequal(observed_labels, unname(COHORT_LABELS))
) {
  fail("Combined EMT table has unexpected rows or cohort labels.")
}

tcga_sample <- build_sample_contract(
  score_tcga,
  combined_emt,
  cohort_id = "TCGA_PRAD",
  score_sample_candidates = c("sample_id", "barcode", "Barcode", "sample_name"),
  expected_n = unname(EXPECTED_N[["TCGA_PRAD"]])
)

friedrich_sample <- build_sample_contract(
  score_friedrich,
  combined_emt,
  cohort_id = "Friedrich_GSE134051",
  score_sample_candidates = c(
    "sample_name", "sample_id", "sample", "Sample", "SampleID",
    "sampleID", "GSM", "gsm"
  ),
  expected_n = unname(EXPECTED_N[["Friedrich_GSE134051"]])
)

tcga_gsea <- build_hallmark_contract(gsea_tcga, "TCGA_PRAD")
friedrich_gsea <- build_hallmark_contract(gsea_friedrich, "Friedrich_GSE134051")

if (!setequal(tcga_gsea$pathway, friedrich_gsea$pathway)) {
  fail("TCGA and Friedrich Hallmark pathway sets are not identical.")
}

contracts <- list(
  tcga_sample = tcga_sample$contract,
  tcga_gsea = tcga_gsea,
  friedrich_sample = friedrich_sample$contract,
  friedrich_gsea = friedrich_gsea
)

for (path in unique(dirname(output_paths))) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    fail("Could not create an authorized Contracts directory.")
  }
}

temporary_paths <- vapply(
  names(contracts),
  function(name) tempfile(pattern = paste0(name, "_"), fileext = ".tsv"),
  character(1)
)
on.exit(unlink(temporary_paths[file.exists(temporary_paths)]), add = TRUE)

for (name in names(contracts)) {
  write_tsv(contracts[[name]], temporary_paths[[name]])
}

staged_tcga_sample <- utils::read.delim(
  temporary_paths[["tcga_sample"]], check.names = FALSE, stringsAsFactors = FALSE
)
staged_friedrich_sample <- utils::read.delim(
  temporary_paths[["friedrich_sample"]], check.names = FALSE, stringsAsFactors = FALSE
)
staged_tcga_gsea <- utils::read.delim(
  temporary_paths[["tcga_gsea"]], check.names = FALSE, stringsAsFactors = FALSE
)
staged_friedrich_gsea <- utils::read.delim(
  temporary_paths[["friedrich_gsea"]], check.names = FALSE, stringsAsFactors = FALSE
)

validate_sample_contract(staged_tcga_sample, "TCGA_PRAD", EXPECTED_N[["TCGA_PRAD"]])
validate_sample_contract(
  staged_friedrich_sample,
  "Friedrich_GSE134051",
  EXPECTED_N[["Friedrich_GSE134051"]]
)
validate_hallmark_contract(staged_tcga_gsea, "TCGA_PRAD")
validate_hallmark_contract(staged_friedrich_gsea, "Friedrich_GSE134051")

if (!setequal(staged_tcga_gsea$pathway, staged_friedrich_gsea$pathway)) {
  fail("Staged Hallmark contracts do not contain identical pathway sets.")
}

for (name in names(contracts)) {
  copied <- file.copy(
    temporary_paths[[name]],
    output_paths[[name]],
    overwrite = arguments$force
  )
  if (!isTRUE(copied)) {
    fail("Could not write an authorized contract output.")
  }
}

cat("Validated sample rows: TCGA_PRAD=497; Friedrich_GSE134051=164\n")
cat("Validated Hallmark rows: TCGA_PRAD=50; Friedrich_GSE134051=50\n")
cat(
  "Maximum Score_z difference: TCGA_PRAD=",
  format(tcga_sample$maximum_score_difference, scientific = FALSE),
  "; Friedrich_GSE134051=",
  format(friedrich_sample$maximum_score_difference, scientific = FALSE),
  "\n",
  sep = ""
)
cat("Contract outputs:\n")
cat(paste0("  ", output_paths, collapse = "\n"), "\n", sep = "")
cat("Final status: PASS\n")
