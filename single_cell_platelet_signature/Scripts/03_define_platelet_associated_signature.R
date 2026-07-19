#!/usr/bin/env Rscript

# 03_define_platelet_associated_signature.R
options(stringsAsFactors = FALSE, scipen = 999)

###############################################################################
# 1. Project directory
###############################################################################

get_project_dir <- function() {
  env_dir <- Sys.getenv("SCORE_CREATION_DIR", unset = "")
  
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = TRUE))
  }
  
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

project_dir <- get_project_dir()
set.seed(123)

###############################################################################
# 2. Run control
###############################################################################

overwrite <- TRUE

###############################################################################
# 3. Required packages
###############################################################################

required_packages <- c("jsonlite")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them outside this script and rerun.",
    call. = FALSE
  )
}

###############################################################################
# 4. Paths
###############################################################################

out_dir <- file.path(project_dir, "Results_MergeSignature")
tables_dir <- file.path(out_dir, "Tables")
signature_dir <- file.path(out_dir, "Signature")
reports_dir <- file.path(out_dir, "Reports")

for (d in c(out_dir, tables_dir, signature_dir, reports_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

bm_table <- file.path(
  project_dir,
  "Results_BoneMarrow",
  "Tables",
  "BM_platelet_associated_lowRNA_cluster_gene_universe_curated.csv"
)

blood_table <- file.path(
  project_dir,
  "Results_Blood",
  "Tables",
  "Blood_platelet_associated_lowRNA_cluster_gene_universe_curated.csv"
)

signature_txt <- file.path(
  signature_dir,
  "platelet_associated_signature.txt"
)

signature_csv <- file.path(
  signature_dir,
  "platelet_associated_signature.csv"
)

signature_rds <- file.path(
  signature_dir,
  "platelet_associated_signature.rds"
)

signature_metadata_json <- file.path(
  signature_dir,
  "platelet_associated_signature_metadata.json"
)

table_overlap <- file.path(
  tables_dir,
  "BM_vs_Blood_gene_overlap.csv"
)

table_shared_before_exclusion <- file.path(
  tables_dir,
  "shared_genes_before_exclusion.csv"
)

table_excluded_genes <- file.path(
  tables_dir,
  "excluded_broadly_expressed_genes.csv"
)

table_derivation <- file.path(
  tables_dir,
  "platelet_associated_signature_derivation_table.csv"
)

table_supplementary_signature <- file.path(
  tables_dir,
  "Supplementary_Table_X_platelet_associated_transcriptional_signature.csv"
)


###############################################################################
# 5. Parameters
###############################################################################

pct_min <- 10
mean_min <- 0.10

delta_pct_filter_used <- FALSE
sig_delta_pct_max <- NULL

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

expected_signature_n <- 41
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

###############################################################################
# 6. Safe write helpers
###############################################################################

safe_write_guard <- function(path, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing file: ",
      path,
      "\nSet overwrite <- TRUE only if this output should be regenerated.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

safe_write_csv <- function(x, path, row.names = FALSE) {
  safe_write_guard(path, overwrite = overwrite)
  utils::write.csv(x, path, row.names = row.names)
  invisible(path)
}

safe_write_lines <- function(text, path) {
  safe_write_guard(path, overwrite = overwrite)
  writeLines(text, con = path)
  invisible(path)
}

safe_write_rds <- function(x, path) {
  safe_write_guard(path, overwrite = overwrite)
  saveRDS(x, file = path)
  invisible(path)
}

###############################################################################
# 7. Logging
###############################################################################

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_file <- file.path(
  reports_dir,
  paste0("MergeSignature_", timestamp_tag, ".log")
)

session_info_file <- file.path(
  reports_dir,
  paste0("MergeSignature_sessionInfo_", timestamp_tag, ".txt")
)

sink(log_file, split = TRUE)
log_closed <- FALSE

close_log_and_write_session <- function() {
  if (isTRUE(log_closed)) {
    return(invisible(NULL))
  }
  
  cat("\n=== SCRIPT EXIT:", as.character(Sys.time()), "===\n")
  
  session_info_error <- tryCatch(
    {
      capture.output(sessionInfo(), file = session_info_file)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  
  if (!is.null(session_info_error)) {
    cat("WARNING: sessionInfo file was not written:", session_info_error, "\n")
  }
  
  if (sink.number() > 0) {
    sink()
  }
  
  log_closed <<- TRUE
  invisible(NULL)
}

on.exit(close_log_and_write_session(), add = TRUE)

cat("=== 03 DEFINE PLATELET-ASSOCIATED SIGNATURE ===\n")
cat("Project directory:", project_dir, "\n")
cat("BM input table:", bm_table, "\n")
cat("Blood input table:", blood_table, "\n")
cat("Output directory:", out_dir, "\n")
cat("Signature directory:", signature_dir, "\n")
cat("Tables directory:", tables_dir, "\n")
cat("Reports directory:", reports_dir, "\n")
cat("Signature TXT:", signature_txt, "\n")
cat("pct_min:", pct_min, "\n")
cat("mean_min:", mean_min, "\n")
cat("delta_pct_filter_used:", delta_pct_filter_used, "\n")
cat("sig_delta_pct_max:", ifelse(is.null(sig_delta_pct_max), "NULL", sig_delta_pct_max), "\n")
cat("Broadly expressed exclusion genes:", paste(broadly_expressed_exclusion_genes, collapse = ", "), "\n\n")

###############################################################################
# 8. Input validation
###############################################################################

if (!file.exists(bm_table)) {
  stop(
    "Missing BM curated universe table. Run script 01 first or check path: ",
    bm_table,
    call. = FALSE
  )
}

if (!file.exists(blood_table)) {
  stop(
    "Missing Blood curated universe table. Run script 02 first or check path: ",
    blood_table,
    call. = FALSE
  )
}

bm <- read.csv(bm_table, stringsAsFactors = FALSE)
blood <- read.csv(blood_table, stringsAsFactors = FALSE)

required_cols <- c("gene", "pct_cells_expr", "mean_counts")

missing_bm_cols <- setdiff(required_cols, colnames(bm))
missing_blood_cols <- setdiff(required_cols, colnames(blood))

if (length(missing_bm_cols) > 0) {
  stop(
    "BM table is missing required column(s): ",
    paste(missing_bm_cols, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_blood_cols) > 0) {
  stop(
    "Blood table is missing required column(s): ",
    paste(missing_blood_cols, collapse = ", "),
    call. = FALSE
  )
}

###############################################################################
# 9. Standardize BM and Blood tables
###############################################################################

bm <- bm[, required_cols]
blood <- blood[, required_cols]

colnames(bm) <- c("gene", "pct_bm", "mean_bm")
colnames(blood) <- c("gene", "pct_blood", "mean_blood")

bm$gene <- trimws(as.character(bm$gene))
blood$gene <- trimws(as.character(blood$gene))

bm <- bm[nzchar(bm$gene), , drop = FALSE]
blood <- blood[nzchar(blood$gene), , drop = FALSE]

if (anyDuplicated(bm$gene) > 0) {
  stop("Duplicated gene symbols found in BM curated universe.", call. = FALSE)
}

if (anyDuplicated(blood$gene) > 0) {
  stop("Duplicated gene symbols found in Blood curated universe.", call. = FALSE)
}

###############################################################################
# 10. Merge BM and Blood platelet-associated gene universes
###############################################################################

gene_overlap <- merge(bm, blood, by = "gene", all = TRUE)

gene_overlap$in_bm <- !is.na(gene_overlap$pct_bm)
gene_overlap$in_blood <- !is.na(gene_overlap$pct_blood)

gene_overlap$class <- ifelse(
  gene_overlap$in_bm & gene_overlap$in_blood,
  "shared",
  ifelse(gene_overlap$in_bm, "bm_only", "blood_only")
)

gene_overlap$pct_bm[is.na(gene_overlap$pct_bm)] <- 0
gene_overlap$mean_bm[is.na(gene_overlap$mean_bm)] <- 0
gene_overlap$pct_blood[is.na(gene_overlap$pct_blood)] <- 0
gene_overlap$mean_blood[is.na(gene_overlap$mean_blood)] <- 0

gene_overlap$delta_pct <- gene_overlap$pct_blood - gene_overlap$pct_bm
gene_overlap$delta_mean <- gene_overlap$mean_blood - gene_overlap$mean_bm

gene_overlap <- gene_overlap[
  order(gene_overlap$class, gene_overlap$gene),
  ,
  drop = FALSE
]

safe_write_csv(gene_overlap, table_overlap, row.names = FALSE)

cat("BM curated genes:", nrow(bm), "\n")
cat("Blood curated genes:", nrow(blood), "\n")
cat("Shared genes:", sum(gene_overlap$class == "shared"), "\n")
cat("BM-only genes:", sum(gene_overlap$class == "bm_only"), "\n")
cat("Blood-only genes:", sum(gene_overlap$class == "blood_only"), "\n\n")

###############################################################################
# 11. Define shared genes before exclusion
###############################################################################

shared_genes_before_exclusion <- subset(
  gene_overlap,
  class == "shared" &
    pct_bm >= pct_min &
    pct_blood >= pct_min &
    mean_bm >= mean_min &
    mean_blood >= mean_min
)

shared_genes_before_exclusion$support_score <- (
  shared_genes_before_exclusion$pct_bm +
    shared_genes_before_exclusion$pct_blood
) +
  10 * (
    shared_genes_before_exclusion$mean_bm +
      shared_genes_before_exclusion$mean_blood
  )

shared_genes_before_exclusion <- shared_genes_before_exclusion[
  order(-shared_genes_before_exclusion$support_score, shared_genes_before_exclusion$gene),
  ,
  drop = FALSE
]

safe_write_csv(
  shared_genes_before_exclusion,
  table_shared_before_exclusion,
  row.names = FALSE
)

cat("Shared genes before exclusion:", nrow(shared_genes_before_exclusion), "\n")

###############################################################################
# 12. Exclude broadly expressed genes
###############################################################################

signature_derivation_table <- shared_genes_before_exclusion

signature_derivation_table$excluded_broadly_expressed <- (
  signature_derivation_table$gene %in% broadly_expressed_exclusion_genes
)

signature_derivation_table$derivation_status <- ifelse(
  signature_derivation_table$excluded_broadly_expressed,
  "excluded_broadly_expressed",
  "included_in_signature"
)

excluded_broadly_expressed_genes <- signature_derivation_table[
  signature_derivation_table$excluded_broadly_expressed,
  c(
    "gene",
    "pct_bm",
    "pct_blood",
    "mean_bm",
    "mean_blood",
    "delta_pct",
    "delta_mean",
    "support_score",
    "derivation_status"
  ),
  drop = FALSE
]

safe_write_csv(
  excluded_broadly_expressed_genes,
  table_excluded_genes,
  row.names = FALSE
)

cat("Excluded broadly expressed genes:", nrow(excluded_broadly_expressed_genes), "\n")
cat(
  "Excluded gene symbols:",
  paste(excluded_broadly_expressed_genes$gene, collapse = ", "),
  "\n"
)

observed_excluded <- excluded_broadly_expressed_genes$gene

if (!setequal(observed_excluded, broadly_expressed_exclusion_genes)) {
  stop(
    "Broadly expressed exclusion validation failed.\n",
    "Expected: ",
    paste(broadly_expressed_exclusion_genes, collapse = ", "),
    "\nObserved: ",
    paste(observed_excluded, collapse = ", "),
    call. = FALSE
  )
}

###############################################################################
# 13. Final platelet-associated transcriptional signature
###############################################################################

platelet_associated_signature <- signature_derivation_table[
  !signature_derivation_table$excluded_broadly_expressed,
  ,
  drop = FALSE
]

platelet_associated_signature <- platelet_associated_signature[
  order(-platelet_associated_signature$support_score, platelet_associated_signature$gene),
  ,
  drop = FALSE
]

platelet_associated_signature$signature_order <- seq_len(nrow(platelet_associated_signature))

signature_derivation_table$signature_status <- ifelse(
  signature_derivation_table$excluded_broadly_expressed,
  "excluded",
  "included"
)

signature_derivation_table <- signature_derivation_table[
  order(
    signature_derivation_table$signature_status,
    -signature_derivation_table$support_score,
    signature_derivation_table$gene
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  signature_derivation_table,
  table_derivation,
  row.names = FALSE
)

###############################################################################
# 14. Signature validation
###############################################################################

observed_signature_n <- nrow(platelet_associated_signature)

if (observed_signature_n != expected_signature_n) {
  stop(
    "Platelet-associated signature must contain ",
    expected_signature_n,
    " genes; observed ",
    observed_signature_n,
    ".",
    call. = FALSE
  )
}

observed_signature_genes <- trimws(as.character(platelet_associated_signature$gene))

if (
  length(observed_signature_genes) != expected_signature_n ||
    anyNA(observed_signature_genes) ||
    any(!nzchar(observed_signature_genes))
) {
  stop(
    "Platelet-associated signature validation failed: all 41 gene symbols must be non-empty and non-missing.",
    call. = FALSE
  )
}

duplicated_signature_genes <- unique(
  observed_signature_genes[duplicated(observed_signature_genes)]
)

if (length(duplicated_signature_genes) > 0) {
  stop(
    "Platelet-associated signature validation failed: duplicated gene symbol(s): ",
    paste(duplicated_signature_genes, collapse = ", "),
    call. = FALSE
  )
}

observed_signature_order <- suppressWarnings(
  as.integer(platelet_associated_signature$signature_order)
)

if (
  length(observed_signature_order) != expected_signature_n ||
    anyNA(observed_signature_order) ||
    !identical(observed_signature_order, seq_len(expected_signature_n)) ||
    any(as.numeric(platelet_associated_signature$signature_order) != observed_signature_order)
) {
  stop(
    "Platelet-associated signature validation failed: signature_order must be exactly 1 through 41.",
    call. = FALSE
  )
}

if (!identical(observed_signature_genes, canonical_signature_genes)) {
  mismatch_index <- which(observed_signature_genes != canonical_signature_genes)[1]
  missing_canonical <- setdiff(canonical_signature_genes, observed_signature_genes)
  unexpected_observed <- setdiff(observed_signature_genes, canonical_signature_genes)

  stop(
    "Platelet-associated signature validation failed: membership or order differs from the canonical 41-gene vector.",
    " First mismatch at signature_order ", mismatch_index,
    ": expected '", canonical_signature_genes[mismatch_index],
    "', observed '", observed_signature_genes[mismatch_index], "'.",
    if (length(missing_canonical) > 0) {
      paste0(" Missing canonical gene(s): ", paste(missing_canonical, collapse = ", "), ".")
    } else {
      ""
    },
    if (length(unexpected_observed) > 0) {
      paste0(" Unexpected gene(s): ", paste(unexpected_observed, collapse = ", "), ".")
    } else {
      ""
    },
    call. = FALSE
  )
}

excluded_genes_in_signature <- intersect(
  broadly_expressed_exclusion_genes,
  observed_signature_genes
)

if (length(excluded_genes_in_signature) > 0) {
  stop(
    "Platelet-associated signature validation failed: excluded broadly expressed gene(s) remain in the final signature: ",
    paste(excluded_genes_in_signature, collapse = ", "),
    call. = FALSE
  )
}

if (!identical(delta_pct_filter_used, FALSE) || !is.null(sig_delta_pct_max)) {
  stop(
    "Invalid delta_pct configuration: delta_pct must be kept as an audit metric only, not as an active filter.",
    call. = FALSE
  )
}



###############################################################################
# 15. Write active signature
###############################################################################

signature_gene_vector <- platelet_associated_signature$gene

safe_write_lines(signature_gene_vector, signature_txt)
safe_write_csv(platelet_associated_signature, signature_csv, row.names = FALSE)
safe_write_rds(platelet_associated_signature, signature_rds)


###############################################################################
# 15A. Publication-ready supplementary signature table
###############################################################################

supplementary_signature_table <- platelet_associated_signature[
  ,
  c("gene", "pct_bm", "pct_blood"),
  drop = FALSE
]

supplementary_signature_table$pct_bm <- round(supplementary_signature_table$pct_bm, 1)
supplementary_signature_table$pct_blood <- round(supplementary_signature_table$pct_blood, 1)

colnames(supplementary_signature_table) <- c(
  "Gene symbol",
  "BM detection (%)",
  "Blood detection (%)"
)

safe_write_csv(
  supplementary_signature_table,
  table_supplementary_signature,
  row.names = FALSE
)

metadata <- list(
  signature_name = "platelet_associated_transcriptional_signature",
  score_name = "platelet_associated_transcriptional_score",
  signature_label = "Platelet-associated transcriptional signature",
  score_label = "Platelet-associated transcriptional score",
  signature_derivation_method = paste(
    "Shared bone marrow and blood platelet-associated gene universes",
    "followed by exact exclusion of broadly expressed genes"
  ),
  single_cell_reference_score_method = list(
    method = "Seurat::AddModuleScore",
    scope = paste(
      "Used only for diagnostic scoring of the single-cell bone marrow",
      "and peripheral blood reference objects."
    )
  ),
  bulk_tumor_score_method_scope = paste(
    "This derivation module does not define one universal bulk-tumor scoring method.",
    "Downstream cohort pipelines must document expression scale, aggregation across",
    "available signature genes, missing-gene handling, and within-cohort standardization."
  ),
  gene_identifier = "Human gene symbols as represented in the source expression matrices",
  signature_order_definition = paste(
    "Descending support_score, with gene symbol used to break ties;",
    "the resulting order must match the canonical 41-gene vector exactly."
  ),
  support_score_formula = "(pct_bm + pct_blood) + 10 * (mean_bm + mean_blood)",
  derivation_module = "single_cell_platelet_signature",
  derivation_script = file.path(
    "Scripts",
    "03_define_platelet_associated_signature.R"
  ),
  inputs = list(
    bm_table = file.path(
      "Results_BoneMarrow",
      "Tables",
      "BM_platelet_associated_lowRNA_cluster_gene_universe_curated.csv"
    ),
    blood_table = file.path(
      "Results_Blood",
      "Tables",
      "Blood_platelet_associated_lowRNA_cluster_gene_universe_curated.csv"
    )
  ),
  parameters = list(
    pct_min = pct_min,
    mean_min = mean_min,
    delta_pct_filter_used = delta_pct_filter_used,
    sig_delta_pct_max = sig_delta_pct_max
  ),
  exclusion = list(
    exclusion_type = "exact_match",
    exclusion_reason = "broadly_expressed_genes",
    excluded_genes = broadly_expressed_exclusion_genes
  ),
  output_files = list(
    signature_txt = file.path(
      "Results_MergeSignature",
      "Signature",
      "platelet_associated_signature.txt"
    ),
    signature_csv = file.path(
      "Results_MergeSignature",
      "Signature",
      "platelet_associated_signature.csv"
    ),
    signature_rds = file.path(
      "Results_MergeSignature",
      "Signature",
      "platelet_associated_signature.rds"
    ),
    metadata_json = file.path(
      "Results_MergeSignature",
      "Signature",
      "platelet_associated_signature_metadata.json"
    ),
    overlap_table = file.path(
      "Results_MergeSignature",
      "Tables",
      "BM_vs_Blood_gene_overlap.csv"
    ),
    shared_before_exclusion_table = file.path(
      "Results_MergeSignature",
      "Tables",
      "shared_genes_before_exclusion.csv"
    ),
    excluded_genes_table = file.path(
      "Results_MergeSignature",
      "Tables",
      "excluded_broadly_expressed_genes.csv"
    ),
    derivation_table = file.path(
      "Results_MergeSignature",
      "Tables",
      "platelet_associated_signature_derivation_table.csv"
    ),
    supplementary_signature_table = file.path(
      "Results_MergeSignature",
      "Tables",
      "Supplementary_Table_X_platelet_associated_transcriptional_signature.csv"
    )
  ),
  counts = list(
    n_bm_curated_genes = nrow(bm),
    n_blood_curated_genes = nrow(blood),
    n_shared_genes = sum(gene_overlap$class == "shared"),
    n_shared_genes_before_exclusion = nrow(shared_genes_before_exclusion),
    n_excluded_broadly_expressed_genes = nrow(excluded_broadly_expressed_genes),
    n_signature_genes = nrow(platelet_associated_signature)
  ),
  canonical_signature_validation = list(
    n_genes = expected_signature_n,
    exact_order_required = TRUE,
    ordered_gene_symbols = canonical_signature_genes,
    excluded_genes_absent = broadly_expressed_exclusion_genes
  ),
  created_at = as.character(Sys.time())
)

safe_write_lines(
  jsonlite::toJSON(
    metadata,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  ),
  signature_metadata_json
)


###############################################################################
# 16. Final report
###############################################################################

cat("\n=== OUTPUTS WRITTEN ===\n")
cat("Signature TXT:", signature_txt, "\n")
cat("Signature CSV:", signature_csv, "\n")
cat("Signature RDS:", signature_rds, "\n")
cat("Metadata JSON:", signature_metadata_json, "\n")
cat("Overlap table:", table_overlap, "\n")
cat("Shared genes before exclusion table:", table_shared_before_exclusion, "\n")
cat("Excluded genes table:", table_excluded_genes, "\n")
cat("Derivation table:", table_derivation, "\n")
cat("Supplementary signature table:", table_supplementary_signature, "\n")

cat("\n=== SIGNATURE COUNTS ===\n")
cat("BM curated genes:", nrow(bm), "\n")
cat("Blood curated genes:", nrow(blood), "\n")
cat("Shared genes:", sum(gene_overlap$class == "shared"), "\n")
cat("Shared genes before exclusion:", nrow(shared_genes_before_exclusion), "\n")
cat("Excluded broadly expressed genes:", nrow(excluded_broadly_expressed_genes), "\n")
cat("Platelet-associated signature genes:", nrow(platelet_associated_signature), "\n")

cat("\n=== ACTIVE SIGNATURE GENES ===\n")
cat(paste(signature_gene_vector, collapse = "\n"), "\n")

cat("\nActive platelet-associated transcriptional signature validation: PASS\n")

close_log_and_write_session()
