################################################################################
# Portable path configuration for the public Friedrich/GSE134051 module.
################################################################################

.resolve_repo_root <- function() {
  repo_override <- Sys.getenv("PLATELET_REPO_ROOT", unset = "")
  invocation_args <- commandArgs(trailingOnly = FALSE)
  file_args <- grep("^--file=", invocation_args, value = TRUE)

  if (nzchar(repo_override)) {
    candidate <- normalizePath(
      repo_override,
      winslash = "/",
      mustWork = TRUE
    )
  } else {
    if (length(file_args) != 1L) {
      stop(
        paste0(
          "Cannot resolve the repository outside Rscript execution. ",
          "Set PLATELET_REPO_ROOT to the public repository root before ",
          "sourcing this configuration."
        ),
        call. = FALSE
      )
    }

    config_file <- tryCatch(
      sys.frame(1)$ofile,
      error = function(e) NULL
    )

    if (is.null(config_file) || !nzchar(config_file)) {
      config_file <- sub("^--file=", "", file_args[[1]])
    }

    config_file <- normalizePath(
      config_file,
      winslash = "/",
      mustWork = TRUE
    )
    module_candidate <- normalizePath(
      file.path(dirname(config_file), ".."),
      winslash = "/",
      mustWork = TRUE
    )
    candidate <- normalizePath(
      file.path(module_candidate, "..", ".."),
      winslash = "/",
      mustWork = TRUE
    )
  }

  required_paths <- c(
    file.path(
      candidate,
      "resources",
      "platelet_associated_transcriptional_signature.tsv"
    ),
    file.path(candidate, "bulk_tumor_validation")
  )

  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0L) {
    stop(
      paste0(
        "Resolved repository root is invalid. Missing: ",
        paste(missing_paths, collapse = "; ")
      ),
      call. = FALSE
    )
  }

  candidate
}

REPO_ROOT <- .resolve_repo_root()
MODULE_DIR <- normalizePath(
  file.path(REPO_ROOT, "bulk_tumor_validation", "friedrich_gse134051"),
  winslash = "/",
  mustWork = TRUE
)
RESOURCE_DIR <- normalizePath(
  file.path(REPO_ROOT, "resources"),
  winslash = "/",
  mustWork = TRUE
)
INPUT_DIR <- normalizePath(
  file.path(MODULE_DIR, "Inputs"),
  winslash = "/",
  mustWork = TRUE
)
GENERATED_RESULTS_DIR <- file.path(MODULE_DIR, "Results", "generated")
GENERATED_FIGURES_DIR <- file.path(MODULE_DIR, "Figures", "generated")

read_canonical_platelet_signature <- function() {
  signature_path <- file.path(
    RESOURCE_DIR,
    "platelet_associated_transcriptional_signature.tsv"
  )

  if (!file.exists(signature_path)) {
    stop(
      paste0("Canonical platelet signature not found: ", signature_path),
      call. = FALSE
    )
  }

  signature_table <- utils::read.delim(
    signature_path,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(signature_table) != 41L) {
    stop("Canonical platelet signature must contain exactly 41 rows.", call. = FALSE)
  }
  if (!"gene_symbol" %in% names(signature_table)) {
    stop("Canonical platelet signature lacks the gene_symbol column.", call. = FALSE)
  }
  if (!"signature_order" %in% names(signature_table)) {
    stop("Canonical platelet signature lacks the signature_order column.", call. = FALSE)
  }

  gene_symbols <- trimws(as.character(signature_table$gene_symbol))
  signature_order <- suppressWarnings(as.integer(signature_table$signature_order))

  if (anyNA(gene_symbols) || any(!nzchar(gene_symbols))) {
    stop("Canonical platelet signature contains empty gene symbols.", call. = FALSE)
  }
  if (anyDuplicated(gene_symbols)) {
    stop("Canonical platelet signature contains duplicated gene symbols.", call. = FALSE)
  }
  if (!identical(signature_order, seq_len(41L))) {
    stop("Canonical platelet signature is not in exact signature_order 1:41.", call. = FALSE)
  }

  forbidden_symbols <- c(
    "B2M", "GAPDH", "ACTB", "PKM", "EIF1", "OAZ1", "SAT1", "OST4"
  )
  forbidden_present <- intersect(gene_symbols, forbidden_symbols)
  if (length(forbidden_present) > 0L) {
    stop(
      paste0(
        "Canonical platelet signature contains forbidden symbols: ",
        paste(forbidden_present, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  signature_table$gene_symbol <- gene_symbols
  signature_table$gene <- gene_symbols
  signature_table
}

rm(.resolve_repo_root)
