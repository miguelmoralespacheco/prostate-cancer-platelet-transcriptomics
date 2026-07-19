#!/usr/bin/env Rscript

# Deterministically export the canonical public platelet-signature resources.
options(stringsAsFactors = FALSE, scipen = 999)

canonical_signature_genes <- c(
  "TMSB4X",
  "PPBP",
  "FTH1",
  "PF4",
  "GNG11",
  "FTL",
  "SH3BGRL3",
  "NRGN",
  "CAVIN2",
  "CCL5",
  "TAGLN2",
  "NCOA4",
  "SERF2",
  "MAP3K7CL",
  "MYL6",
  "CLU",
  "TUBB1",
  "ITM2B",
  "PTMA",
  "MYL12A",
  "RGS18",
  "CALM3",
  "SPARC",
  "PGRMC1",
  "NAP1L1",
  "RGS10",
  "ARPC1B",
  "ACRBP",
  "PPDPF",
  "ACTG1",
  "TUBA4A",
  "TREML1",
  "MMD",
  "GP9",
  "FKBP1A",
  "LAMTOR1",
  "CTSA",
  "HMGB1",
  "TIMP1",
  "GPX4",
  "FERMT3"
)

broadly_expressed_exclusion_genes <- c(
  "B2M",
  "GAPDH",
  "ACTB",
  "PKM",
  "EIF1",
  "OAZ1",
  "SAT1",
  "OST4"
)

public_columns <- c(
  "signature_order",
  "gene_symbol",
  "bone_marrow_detection_pct",
  "blood_detection_pct"
)

get_script_path <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(file_arg) == 0) {
    return(NULL)
  }

  normalizePath(
    sub("^--file=", "", file_arg[1]),
    winslash = "/",
    mustWork = TRUE
  )
}

detect_repository_root <- function() {
  script_path <- get_script_path()

  if (!is.null(script_path)) {
    module_dir <- dirname(dirname(script_path))
  } else {
    current_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

    if (basename(current_dir) == "single_cell_platelet_signature") {
      module_dir <- current_dir
    } else if (dir.exists(file.path(current_dir, "single_cell_platelet_signature"))) {
      module_dir <- file.path(current_dir, "single_cell_platelet_signature")
    } else {
      stop(
        "Cannot detect the single_cell_platelet_signature module. Run this script with Rscript or set PUBLIC_RESOURCE_DIR.",
        call. = FALSE
      )
    }
  }

  module_dir <- normalizePath(module_dir, winslash = "/", mustWork = TRUE)

  if (basename(module_dir) != "single_cell_platelet_signature") {
    stop(
      "Detected module directory has an unexpected name: ",
      module_dir,
      call. = FALSE
    )
  }

  repository_root <- normalizePath(
    dirname(module_dir),
    winslash = "/",
    mustWork = TRUE
  )

  if (!dir.exists(file.path(repository_root, "single_cell_platelet_signature"))) {
    stop(
      "Cannot validate repository root: ",
      repository_root,
      call. = FALSE
    )
  }

  repository_root
}

resolve_score_creation_root <- function() {
  configured <- Sys.getenv("SCORE_CREATION_DIR", unset = "")
  candidate <- if (nzchar(configured)) configured else "."
  candidate <- path.expand(candidate)

  if (!dir.exists(candidate)) {
    stop(
      "SCORE_CREATION_DIR does not resolve to an existing directory: ",
      candidate,
      call. = FALSE
    )
  }

  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

resolve_public_resource_dir <- function() {
  configured <- Sys.getenv("PUBLIC_RESOURCE_DIR", unset = "")

  candidate <- if (nzchar(configured)) {
    path.expand(configured)
  } else {
    file.path(detect_repository_root(), "resources")
  }

  if (!dir.exists(candidate)) {
    created <- dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(created) && !dir.exists(candidate)) {
      stop(
        "Could not create PUBLIC_RESOURCE_DIR: ",
        candidate,
        call. = FALSE
      )
    }
  }

  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

sha256_file <- function(path) {
  commands <- list(
    list(command = "shasum", args = c("-a", "256")),
    list(command = "sha256sum", args = character()),
    list(command = "openssl", args = c("dgst", "-sha256"))
  )

  for (candidate in commands) {
    executable <- Sys.which(candidate$command)
    if (!nzchar(executable)) {
      next
    }

    output <- suppressWarnings(
      system2(
        executable,
        args = c(candidate$args, shQuote(path)),
        stdout = TRUE,
        stderr = TRUE
      )
    )
    status <- attr(output, "status")

    if (is.null(status) || identical(status, 0L)) {
      matches <- unlist(
        regmatches(output, gregexpr("[[:xdigit:]]{64}", output)),
        use.names = FALSE
      )
      if (length(matches) > 0) {
        return(tolower(matches[1]))
      }
    }
  }

  stop(
    "Unable to calculate SHA-256. Install or expose shasum, sha256sum, or openssl in PATH.",
    call. = FALSE
  )
}

prefix_module_path <- function(path) {
  if (!is.character(path) || length(path) != 1 || !nzchar(path)) {
    stop("Invalid internal derivation path in metadata.", call. = FALSE)
  }

  normalized <- gsub("\\\\", "/", path)

  if (grepl("^single_cell_platelet_signature/", normalized)) {
    return(normalized)
  }

  if (grepl("^/|^[A-Za-z]:/", normalized)) {
    stop(
      "Absolute internal path is not allowed in public metadata: ",
      normalized,
      call. = FALSE
    )
  }

  file.path("single_cell_platelet_signature", normalized)
}

validate_signature_source <- function(source_table) {
  required_columns <- c("signature_order", "gene", "pct_bm", "pct_blood")
  missing_columns <- setdiff(required_columns, names(source_table))

  if (length(missing_columns) > 0) {
    stop(
      "Canonical signature CSV is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(source_table) != 41L) {
    stop(
      "Canonical signature CSV must contain exactly 41 rows; observed ",
      nrow(source_table),
      ".",
      call. = FALSE
    )
  }

  genes <- trimws(source_table$gene)

  if (anyNA(genes) || any(!nzchar(genes))) {
    stop("Canonical signature CSV contains a missing or empty gene symbol.", call. = FALSE)
  }

  duplicated_genes <- unique(genes[duplicated(genes)])
  if (length(duplicated_genes) > 0) {
    stop(
      "Canonical signature CSV contains duplicated gene symbol(s): ",
      paste(duplicated_genes, collapse = ", "),
      call. = FALSE
    )
  }

  observed_order <- suppressWarnings(as.integer(source_table$signature_order))
  if (
    anyNA(observed_order) ||
      !identical(observed_order, seq_len(41L)) ||
      any(as.numeric(source_table$signature_order) != observed_order)
  ) {
    stop("signature_order must be exactly 1 through 41.", call. = FALSE)
  }

  if (!identical(genes, canonical_signature_genes)) {
    mismatch_index <- which(genes != canonical_signature_genes)[1]
    stop(
      "Canonical signature membership or order mismatch at position ",
      mismatch_index,
      ": expected '", canonical_signature_genes[mismatch_index],
      "', observed '", genes[mismatch_index], "'.",
      call. = FALSE
    )
  }

  excluded_present <- intersect(broadly_expressed_exclusion_genes, genes)
  if (length(excluded_present) > 0) {
    stop(
      "Excluded broadly expressed gene(s) found in canonical signature: ",
      paste(excluded_present, collapse = ", "),
      call. = FALSE
    )
  }

  for (column in c("pct_bm", "pct_blood")) {
    values <- source_table[[column]]
    numeric_values <- suppressWarnings(as.numeric(values))
    if (anyNA(numeric_values) || any(!is.finite(numeric_values))) {
      stop(
        "Canonical signature CSV contains a non-numeric or non-finite value in ",
        column,
        ".",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

build_public_table <- function(source_table) {
  public_table <- data.frame(
    signature_order = source_table$signature_order,
    gene_symbol = trimws(source_table$gene),
    bone_marrow_detection_pct = source_table$pct_bm,
    blood_detection_pct = source_table$pct_blood,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  names(public_table) <- public_columns
  public_table
}

validate_written_tsv <- function(path, expected_table) {
  observed <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    colClasses = "character",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!identical(names(observed), public_columns)) {
    stop("Written TSV header does not match the required public schema.", call. = FALSE)
  }

  if (!identical(observed, expected_table)) {
    stop("Written TSV content differs from the validated source values.", call. = FALSE)
  }

  invisible(TRUE)
}

build_public_metadata <- function(source_metadata, source_csv, public_tsv) {
  unscoped_method <- source_metadata$intended_score_method
  if (!is.null(unscoped_method) && !identical(unscoped_method, "Seurat::AddModuleScore")) {
    stop(
      "Unexpected intended_score_method in source metadata: ",
      paste(unlist(unscoped_method), collapse = " "),
      call. = FALSE
    )
  }

  source_metadata$intended_score_method <- NULL
  source_metadata$single_cell_reference_score_method <- list(
    method = "Seurat::AddModuleScore",
    scope = paste(
      "Used only for diagnostic scoring of the single-cell bone marrow",
      "and peripheral blood reference objects."
    )
  )
  source_metadata$bulk_tumor_score_method_scope <- paste(
    "This derivation module does not define one universal bulk-tumor scoring method.",
    "Downstream cohort pipelines must document expression scale, aggregation across",
    "available signature genes, missing-gene handling, and within-cohort standardization."
  )
  source_metadata$bulk_tumor_score_definitions <- NULL
  source_metadata$gene_identifier <- (
    "Human gene symbols as represented in the source expression matrices"
  )
  source_metadata$signature_order_definition <- paste(
    "Descending support_score, with gene symbol used to break ties;",
    "the resulting order must match the canonical 41-gene vector exactly."
  )
  source_metadata$support_score_formula <- (
    "(pct_bm + pct_blood) + 10 * (mean_bm + mean_blood)"
  )
  source_metadata$derivation_module <- "single_cell_platelet_signature"
  source_metadata$derivation_script <- (
    "single_cell_platelet_signature/Scripts/03_define_platelet_associated_signature.R"
  )

  if (!is.null(source_metadata$inputs)) {
    source_metadata$inputs <- lapply(
      source_metadata$inputs,
      prefix_module_path
    )
  }

  internal_output_paths <- source_metadata$output_files
  if (is.null(internal_output_paths)) {
    internal_output_paths <- source_metadata$generated_derivation_outputs$paths
  }
  if (is.null(internal_output_paths) || length(internal_output_paths) == 0) {
    stop("Source metadata does not identify generated derivation outputs.", call. = FALSE)
  }

  internal_output_paths <- lapply(internal_output_paths, prefix_module_path)
  source_metadata$output_files <- NULL
  source_metadata$generated_derivation_outputs <- list(
    tracked_in_repository = FALSE,
    note = paste(
      "These internal Results_* files are generated at runtime by the",
      "single-cell signature module and are intentionally not tracked in the repository."
    ),
    paths = internal_output_paths
  )

  source_metadata$expected_signature_genes_checked <- NULL
  source_metadata$canonical_signature_validation <- list(
    n_genes = 41L,
    exact_order_required = TRUE,
    ordered_gene_symbols = canonical_signature_genes,
    excluded_genes_absent = broadly_expressed_exclusion_genes
  )
  source_metadata$public_resource <- list(
    file = "resources/platelet_associated_transcriptional_signature.tsv",
    format = "TSV",
    n_signature_genes = 41L,
    columns = public_columns,
    detection_percentage_precision = (
      "Full precision retained from pct_bm and pct_blood"
    ),
    source_canonical_csv_sha256 = sha256_file(source_csv),
    public_tsv_sha256 = sha256_file(public_tsv)
  )

  source_metadata
}

commit_outputs <- function(temp_paths, final_paths, overwrite) {
  existing <- file.exists(final_paths)

  if (any(existing) && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing public resource(s): ",
      paste(final_paths[existing], collapse = ", "),
      ". Set PUBLIC_RESOURCE_OVERWRITE=true only after validating the export.",
      call. = FALSE
    )
  }

  backup_paths <- rep(NA_character_, length(final_paths))
  backup_moved <- rep(FALSE, length(final_paths))
  output_moved <- rep(FALSE, length(final_paths))

  rollback <- function() {
    for (index in which(output_moved)) {
      if (file.exists(final_paths[index])) {
        unlink(final_paths[index])
      }
    }
    for (index in rev(which(backup_moved))) {
      if (file.exists(backup_paths[index])) {
        file.rename(backup_paths[index], final_paths[index])
      }
    }
  }

  tryCatch(
    {
      for (index in which(existing)) {
        backup_paths[index] <- tempfile(
          pattern = paste0(".", basename(final_paths[index]), ".backup."),
          tmpdir = dirname(final_paths[index])
        )
        if (!file.rename(final_paths[index], backup_paths[index])) {
          stop("Could not stage existing output for replacement: ", final_paths[index])
        }
        backup_moved[index] <- TRUE
      }

      for (index in seq_along(final_paths)) {
        if (!file.rename(temp_paths[index], final_paths[index])) {
          stop("Could not atomically install validated output: ", final_paths[index])
        }
        output_moved[index] <- TRUE
      }
    },
    error = function(error) {
      rollback()
      stop(conditionMessage(error), call. = FALSE)
    }
  )

  for (index in which(backup_moved)) {
    if (file.exists(backup_paths[index])) {
      unlink(backup_paths[index])
    }
  }

  invisible(final_paths)
}

main <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Required package 'jsonlite' is unavailable. Install it outside this script and rerun.",
      call. = FALSE
    )
  }

  score_creation_root <- resolve_score_creation_root()
  public_resource_dir <- resolve_public_resource_dir()
  overwrite <- identical(
    tolower(trimws(Sys.getenv("PUBLIC_RESOURCE_OVERWRITE", unset = "false"))),
    "true"
  )

  source_csv <- file.path(
    score_creation_root,
    "Results_MergeSignature",
    "Signature",
    "platelet_associated_signature.csv"
  )
  source_metadata_json <- file.path(
    score_creation_root,
    "Results_MergeSignature",
    "Signature",
    "platelet_associated_signature_metadata.json"
  )

  missing_inputs <- c(source_csv, source_metadata_json)[
    !file.exists(c(source_csv, source_metadata_json))
  ]
  if (length(missing_inputs) > 0) {
    stop(
      "Missing required text-readable signature input(s): ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }

  output_tsv <- file.path(
    public_resource_dir,
    "platelet_associated_transcriptional_signature.tsv"
  )
  output_metadata_json <- file.path(
    public_resource_dir,
    "platelet_associated_transcriptional_signature_metadata.json"
  )

  existing_outputs <- c(output_tsv, output_metadata_json)[
    file.exists(c(output_tsv, output_metadata_json))
  ]
  if (length(existing_outputs) > 0 && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing public resource(s): ",
      paste(existing_outputs, collapse = ", "),
      ". Set PUBLIC_RESOURCE_OVERWRITE=true only after validating the export.",
      call. = FALSE
    )
  }

  source_table <- utils::read.csv(
    source_csv,
    header = TRUE,
    colClasses = "character",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  validate_signature_source(source_table)
  public_table <- build_public_table(source_table)

  temp_tsv <- tempfile(
    pattern = ".platelet_associated_transcriptional_signature.",
    tmpdir = public_resource_dir,
    fileext = ".tsv.tmp"
  )
  temp_metadata_json <- tempfile(
    pattern = ".platelet_associated_transcriptional_signature_metadata.",
    tmpdir = public_resource_dir,
    fileext = ".json.tmp"
  )
  temporary_paths <- c(temp_tsv, temp_metadata_json)
  on.exit(unlink(temporary_paths[file.exists(temporary_paths)]), add = TRUE)

  utils::write.table(
    public_table,
    file = temp_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    eol = "\n",
    fileEncoding = "UTF-8"
  )
  validate_written_tsv(temp_tsv, public_table)

  source_metadata <- jsonlite::fromJSON(
    source_metadata_json,
    simplifyVector = FALSE
  )
  public_metadata <- build_public_metadata(
    source_metadata,
    source_csv,
    temp_tsv
  )
  metadata_text <- jsonlite::toJSON(
    public_metadata,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
  writeLines(metadata_text, con = temp_metadata_json, useBytes = TRUE)

  metadata_json_text <- paste(
    readLines(temp_metadata_json, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  if (!jsonlite::validate(metadata_json_text)) {
    stop("Generated public metadata is not valid JSON.", call. = FALSE)
  }

  validated_metadata <- jsonlite::fromJSON(
    temp_metadata_json,
    simplifyVector = FALSE
  )
  if (!is.null(validated_metadata$intended_score_method)) {
    stop("Generated metadata retains unscoped intended_score_method.", call. = FALSE)
  }
  if (
    !identical(
      validated_metadata$public_resource$public_tsv_sha256,
      sha256_file(temp_tsv)
    )
  ) {
    stop("Generated metadata contains an incorrect public TSV SHA-256.", call. = FALSE)
  }

  commit_outputs(
    temp_paths = temporary_paths,
    final_paths = c(output_tsv, output_metadata_json),
    overwrite = overwrite
  )

  message("PASS: validated canonical ordered signature (41 genes).")
  message("PASS: public TSV written: ", output_tsv)
  message("PASS: public metadata written: ", output_metadata_json)
  message("PASS: deterministic public signature export completed.")

  invisible(list(tsv = output_tsv, metadata = output_metadata_json))
}

tryCatch(
  main(),
  error = function(error) {
    message("FAIL: ", conditionMessage(error))
    quit(save = "no", status = 1L, runLast = FALSE)
  }
)
