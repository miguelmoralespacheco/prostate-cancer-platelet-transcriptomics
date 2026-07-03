#!/usr/bin/env Rscript

###############################################################################
# 04_reactome_ORA_platelet_associated_signature.R

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

###############################################################################
# 2. Run control
###############################################################################

overwrite <- TRUE

###############################################################################
# 3. Required packages
###############################################################################

required_packages <- c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "ggplot2",
  "stringr",
  "jsonlite"
)

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

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
})

###############################################################################
# 4. Paths
###############################################################################

signature_file <- file.path(
  project_dir,
  "Results_MergeSignature",
  "Signature",
  "platelet_associated_signature.txt"
)

out_dir <- file.path(project_dir, "Results_Reactome")
tables_dir <- file.path(out_dir, "Tables")
figures_dir <- file.path(out_dir, "Figures")
reports_dir <- file.path(out_dir, "Reports")

for (d in c(out_dir, tables_dir, figures_dir, reports_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

mapping_csv <- file.path(
  tables_dir,
  "SYMBOL_to_ENTREZ_mapping.csv"
)

unmapped_txt <- file.path(
  tables_dir,
  "UNMAPPED_symbols.txt"
)

reactome_full_csv <- file.path(
  tables_dir,
  "ORA_Reactome_platelet_signature_full.csv"
)

reactome_top3_csv <- file.path(
  tables_dir,
  "ORA_Reactome_platelet_signature_top3.csv"
)

reactome_metadata_json <- file.path(
  tables_dir,
  "ORA_Reactome_platelet_signature_metadata.json"
)

reactome_top3_pdf <- file.path(
  figures_dir,
  "FIG_Reactome_ORA_TOP3.pdf"
)

###############################################################################
# 5. Parameters
###############################################################################

expected_signature_n <- 41
expected_signature_genes <- c("NRGN", "NCOA4")
topn <- 3

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

save_plot_pdf <- function(plot_obj, file_name, width = 6, height = 5) {
  safe_write_guard(file_name, overwrite = overwrite)
  
  grDevices::pdf(
    file = file_name,
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj)
  
  invisible(file_name)
}

write_unmapped_symbols <- function(symbols, path) {
  if (length(symbols) == 0) {
    symbols <- "NONE"
  }
  
  safe_write_lines(symbols, path)
  invisible(path)
}

###############################################################################
# 7. Logging
###############################################################################

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_file <- file.path(
  reports_dir,
  paste0("Reactome_ORA_", timestamp_tag, ".log")
)

session_info_file <- file.path(
  reports_dir,
  paste0("Reactome_ORA_sessionInfo_", timestamp_tag, ".txt")
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

cat("=== 04 REACTOME ORA ===\n")
cat("Project directory:", project_dir, "\n")
cat("Signature input:", signature_file, "\n")
cat("Results directory:", out_dir, "\n")
cat("Tables directory:", tables_dir, "\n")
cat("Figures directory:", figures_dir, "\n")
cat("Reports directory:", reports_dir, "\n\n")

###############################################################################
# 8. Input validation
###############################################################################

if (!file.exists(signature_file)) {
  stop(
    "Missing platelet-associated signature file. Run script 03 first or check path: ",
    signature_file,
    call. = FALSE
  )
}

genes <- unique(trimws(readLines(signature_file, warn = FALSE)))
genes <- genes[nzchar(genes)]

cat("Input signature genes:", length(genes), "\n")

if (length(genes) != expected_signature_n) {
  stop(
    "Expected ",
    expected_signature_n,
    " platelet-associated signature genes; observed ",
    length(genes),
    ". Check script 03 output.",
    call. = FALSE
  )
}

missing_expected_genes <- setdiff(expected_signature_genes, genes)

if (length(missing_expected_genes) > 0) {
  stop(
    "Input signature is missing expected gene(s): ",
    paste(missing_expected_genes, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(genes) > 0) {
  stop("Duplicated gene symbols found in input signature.", call. = FALSE)
}

###############################################################################
# 9. SYMBOL to ENTREZ mapping
###############################################################################

cat("\nMapping SYMBOL to ENTREZ...\n")

conversion <- clusterProfiler::bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)

conversion <- as.data.frame(conversion, stringsAsFactors = FALSE)

if (nrow(conversion) == 0) {
  stop("No genes could be mapped from SYMBOL to ENTREZID.", call. = FALSE)
}

conversion <- conversion[!duplicated(conversion[, c("SYMBOL", "ENTREZID")]), , drop = FALSE]

mapped_symbols <- unique(conversion$SYMBOL)
unmapped_symbols <- setdiff(genes, mapped_symbols)
entrez <- unique(conversion$ENTREZID)

safe_write_csv(conversion, mapping_csv, row.names = FALSE)
write_unmapped_symbols(unmapped_symbols, unmapped_txt)

cat("Mapped SYMBOL:", length(mapped_symbols), "\n")
cat("Unmapped SYMBOL:", length(unmapped_symbols), "\n")
cat("Mapped ENTREZ:", length(entrez), "\n")

if (length(entrez) < 10) {
  stop(
    "Too few mapped ENTREZ IDs for Reactome ORA: ",
    length(entrez),
    call. = FALSE
  )
}

###############################################################################
# 10. Reactome ORA
###############################################################################

cat("\nRunning Reactome ORA using Reactome annotated/default universe...\n")

reactome_result <- ReactomePA::enrichPathway(
  gene = entrez,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

reactome_df <- as.data.frame(reactome_result)

if (is.null(reactome_df)) {
  reactome_df <- data.frame()
}

expected_empty_cols <- c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "qvalue",
  "geneID",
  "Count"
)

if (nrow(reactome_df) == 0) {
  reactome_df <- as.data.frame(
    stats::setNames(
      replicate(length(expected_empty_cols), character(0), simplify = FALSE),
      expected_empty_cols
    )
  )
  
  safe_write_csv(reactome_df, reactome_full_csv, row.names = FALSE)
  safe_write_csv(reactome_df, reactome_top3_csv, row.names = FALSE)
  
  cat("Reactome ORA returned no enriched terms. Empty tables were written.\n")
} else {
  reactome_df <- reactome_df[
    order(reactome_df$p.adjust, reactome_df$pvalue),
    ,
    drop = FALSE
  ]
  
  safe_write_csv(reactome_df, reactome_full_csv, row.names = FALSE)
  
  top_df <- reactome_df[
    seq_len(min(topn, nrow(reactome_df))),
    ,
    drop = FALSE
  ]
  
  safe_write_csv(top_df, reactome_top3_csv, row.names = FALSE)
  
  cat("Reactome ORA rows:", nrow(reactome_df), "\n")
  cat("Top Reactome terms:\n")
  print(top_df[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "Count")])
}

###############################################################################
# 11. Plot top Reactome terms
###############################################################################

if (nrow(reactome_df) > 0) {
  plot_df <- top_df
  
  plot_df$term <- stringr::str_wrap(plot_df$Description, width = 28)
  
  plot_df <- plot_df[
    order(plot_df$Count, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  plot_df$term <- factor(plot_df$term, levels = rev(plot_df$term))
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Count, y = term)
  ) +
    ggplot2::geom_col(
      width = 0.65,
      fill = "#AAC4F5",
      color = "grey60",
      linewidth = 0.12
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 4, 8, 12),
      limits = c(-0.4, 12.5),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::theme_classic(
      base_family = "Helvetica",
      base_size = 5
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica", face = "plain"),
      axis.title.x = ggplot2::element_text(family = "Helvetica", size = 4),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(family = "Helvetica", size = 4),
      axis.text.y = ggplot2::element_text(family = "Helvetica", size = 4.5),
      axis.line = ggplot2::element_line(linewidth = 0.25),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      axis.ticks.length = grid::unit(0.08, "cm"),
      plot.title = ggplot2::element_text(
        family = "Helvetica",
        size = 5,
        face = "plain",
        hjust = 0.5
      ),
      plot.margin = grid::unit(c(0.4, 0.4, 0.3, 0.4), "mm")
    ) +
    ggplot2::labs(
      x = "Gene count",
      y = NULL,
      title = "Reactome pathway enrichment"
    )
  
  save_plot_pdf(
    p,
    reactome_top3_pdf,
    width = 5.0 / 2.54,
    height = 2.1 / 2.54
  )
}

###############################################################################
# 12. Metadata
###############################################################################

metadata <- list(
  analysis_name = "Reactome_ORA_platelet_associated_signature",
  input_signature = file.path(
    "Results_MergeSignature",
    "Signature",
    "platelet_associated_signature.txt"
  ),
  signature_name = "platelet_associated_transcriptional_signature",
  n_input_genes = length(genes),
  n_mapped_symbols = length(mapped_symbols),
  n_unmapped_symbols = length(unmapped_symbols),
  unmapped_symbols = if (length(unmapped_symbols) == 0) "NONE" else unmapped_symbols,
  n_mapped_entrez = length(entrez),
  ora_method = "ReactomePA::enrichPathway",
  organism = "human",
  p_adjust_method = "BH",
  pvalue_cutoff = 1,
  qvalue_cutoff = 1,
  background = "Reactome annotated/default universe",
  interpretation = paste(
    "Functional interpretation of the active 41-gene platelet-associated",
    "transcriptional signature; not tested against the BM/Blood candidate universe."
  ),
  output_files = list(
    mapping_csv = file.path(
      "Results_Reactome",
      "Tables",
      "SYMBOL_to_ENTREZ_mapping.csv"
    ),
    unmapped_txt = file.path(
      "Results_Reactome",
      "Tables",
      "UNMAPPED_symbols.txt"
    ),
    reactome_full_csv = file.path(
      "Results_Reactome",
      "Tables",
      "ORA_Reactome_platelet_signature_full.csv"
    ),
    reactome_top3_csv = file.path(
      "Results_Reactome",
      "Tables",
      "ORA_Reactome_platelet_signature_top3.csv"
    ),
    figure_top3_pdf = file.path(
      "Results_Reactome",
      "Figures",
      "FIG_Reactome_ORA_TOP3.pdf"
    ),
    metadata_json = file.path(
      "Results_Reactome",
      "Tables",
      "ORA_Reactome_platelet_signature_metadata.json"
    )
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
  reactome_metadata_json
)

###############################################################################
# 13. Final report
###############################################################################

cat("\n=== OUTPUTS WRITTEN ===\n")
cat("SYMBOL to ENTREZ mapping:", mapping_csv, "\n")
cat("Unmapped symbols:", unmapped_txt, "\n")
cat("Reactome full table:", reactome_full_csv, "\n")
cat("Reactome top3 table:", reactome_top3_csv, "\n")
cat("Reactome metadata:", reactome_metadata_json, "\n")

if (nrow(reactome_df) > 0) {
  cat("Reactome top3 figure:", reactome_top3_pdf, "\n")
}

cat("\n=== REACTOME COUNTS ===\n")
cat("Input signature genes:", length(genes), "\n")
cat("Mapped SYMBOL:", length(mapped_symbols), "\n")
cat("Unmapped SYMBOL:", length(unmapped_symbols), "\n")
cat("Mapped ENTREZ:", length(entrez), "\n")
cat("Reactome rows:", nrow(reactome_df), "\n")

cat("\nReactome ORA validation: PASS\n")

close_log_and_write_session()
