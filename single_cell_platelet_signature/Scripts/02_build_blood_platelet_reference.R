###############################################################################
# Build Blood platelet-associated single-cell reference
# Clean reproduction script for Score_creation. Does not install packages.
###############################################################################

get_project_dir <- function() {
  env_dir <- Sys.getenv("SCORE_CREATION_DIR", unset = "")
  if (nzchar(env_dir)) return(normalizePath(env_dir, mustWork = TRUE))
  normalizePath(".", mustWork = TRUE)
}

project_dir <- get_project_dir()
set.seed(123)

tissue_label <- "Blood"
prefix <- "BLOOD"
manifest_path <- file.path(project_dir, "Inputs", "Blood", "Manifests", "manifest_BL_HiSeq9.tsv")
output_dir <- file.path(project_dir, "Results_Blood")
tables_dir <- file.path(output_dir, "Tables")
figures_dir <- file.path(output_dir, "Figures")
reports_dir <- file.path(output_dir, "Reports")
objects_dir <- file.path(output_dir, "Objects")
for (d in c(tables_dir, figures_dir, reports_dir, objects_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

true_platelet_cluster <- "3"
min_features_all <- 30
min_cells_gene <- 3
plt_nFeature_min <- 30
plt_nFeature_max <- 200
plt_nCount_min <- 50
plt_nCount_max <- 2000
plt_percent_mt_max <- 25
n_pcs <- 20
resolution <- 0.5
min_pct_expr <- 20
min_mean_count <- 0.05
threshold_pct <- 10
overwrite <- TRUE

parameter_contract <- list(
  official_score_name = "platelet-associated transcriptional score",
  method = "Seurat::AddModuleScore",
  source = "single-cell bone marrow and peripheral blood platelet/MK-associated clusters",
  blood_platelet_cluster = true_platelet_cluster,
  nFeature_RNA = paste0(plt_nFeature_min, "-", plt_nFeature_max),
  nCount_RNA = paste0(plt_nCount_min, "-", plt_nCount_max),
  percent_mt = paste0("<=", plt_percent_mt_max),
  pct_cells_expr = paste0(">=", min_pct_expr),
  mean_counts = paste0(">", min_mean_count)
)

markers_platelet <- c("PPBP", "PF4", "NRGN", "SPARC", "SDPR", "ITGA2B", "GP9", "NBEAL2", "TUBB1", "CLU")
markers_myeloid <- c("LYZ", "S100A8", "S100A9", "MPO", "PRTN3", "FCGR3B")
markers_lymphoid <- c("CD3D", "CD3E", "MS4A1", "CD79A", "JCHAIN")
markers_erythro <- c("HBB", "HBA1", "HBA2", "ALAS2", "CA1")

required_packages <- c("Matrix", "hdf5r", "Seurat", "SeuratObject", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       ". Install them outside this script and rerun.")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(reports_dir, paste0(prefix, "_score_creation_", timestamp_tag, ".log"))
session_info_file <- file.path(reports_dir, paste0(prefix, "_sessionInfo_", timestamp_tag, ".txt"))
sink(log_file, split = TRUE)
log_closed <- FALSE
close_log_and_write_session <- function() {
  if (isTRUE(log_closed)) return(invisible(NULL))
  cat("\n=== SCRIPT EXIT:", as.character(Sys.time()), "===\n")
  session_info_error <- tryCatch({
    capture.output(sessionInfo(), file = session_info_file)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(session_info_error)) {
    cat("WARNING: sessionInfo file was not written:", session_info_error, "\n")
  }
  if (sink.number() > 0) sink()
  log_closed <<- TRUE
  invisible(NULL)
}
on.exit(close_log_and_write_session(), add = TRUE)

cat("Project dir:", project_dir, "\n")
cat("Manifest:", manifest_path, "\n")
cat("Output dir:", output_dir, "\n")

theme_pub <- function() {
  ggplot2::theme_classic(base_family = "Helvetica", base_size = 5) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica", size = 5, face = "plain"),
      axis.title = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain"),
      axis.text = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain"),
      plot.title = ggplot2::element_text(family = "Helvetica", size = 5, face = "plain"),
      strip.text = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain"),
      legend.title = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain"),
      legend.text = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain"),
      axis.line = ggplot2::element_line(linewidth = 0.25),
      axis.ticks = ggplot2::element_line(linewidth = 0.25)
    )
}

save_plot_pdf <- function(plot_obj, file_name, width = 6, height = 5) {
  safe_write_guard(file_name, overwrite = overwrite)
  grDevices::pdf(file_name, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  print(plot_obj)
  grDevices::dev.off()
  cat("Saved plot:", file_name, "\n")
  invisible(file_name)
}

save_plot_grid_pdf <- function(plot_list, file_name, width = 6, height = 5, ncol = 2) {
  safe_write_guard(file_name, overwrite = overwrite)
  nrow <- ceiling(length(plot_list) / ncol)
  grDevices::pdf(file_name, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
  for (i in seq_along(plot_list)) {
    row_i <- ceiling(i / ncol)
    col_i <- ((i - 1) %% ncol) + 1
    print(
      plot_list[[i]],
      vp = grid::viewport(layout.pos.row = row_i, layout.pos.col = col_i)
    )
  }
  cat("Saved plot:", file_name, "\n")
  invisible(file_name)
}

set_point_style <- function(plot_obj, size = 0.5, alpha = 0.35) {
  if (length(plot_obj$layers) > 0) {
    plot_obj$layers[[1]]$aes_params$shape <- 16
    plot_obj$layers[[1]]$aes_params$stroke <- 0
    plot_obj$layers[[1]]$aes_params$size <- size
    plot_obj$layers[[1]]$aes_params$alpha <- alpha
  }
  plot_obj
}

sort_cluster_levels <- function(x) {
  x <- as.character(x)
  if (all(grepl("^[0-9]+$", x))) return(as.character(sort(as.integer(x))))
  sort(x)
}

cluster_palette <- function(n) {
  base_cols <- c(
    "#004C94",   
    "#9151B8",   
    "#FF8F5C",  
    "#DB1A1A",  
    "#81912F",  
    "#00D4C6",
    "#9F431C",  
    "#3C93FA",  
    "#FF57B0",  
    "#ADDBFF"   
  )
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  grDevices::hcl.colors(n, palette = "Set 3")
}

platelet_focus_col <- "#DB1A1A"
context_cell_col <- "grey85"
dot_low_col <- "#FEE5D9" 
dot_high_col <- "#A80000"

get_last_score_col <- function(obj, prefix) {
  md <- colnames(obj@meta.data)
  hits <- md[grepl(paste0("^", prefix, "[0-9]+$"), md)]
  if (length(hits) == 0) return(NA_character_)
  ord <- order(as.integer(gsub(paste0("^", prefix), "", hits)))
  hits[ord][length(hits)]
}

squish_oob <- function(x, range = c(0, 1), only.finite = TRUE) {
  finite <- if (isTRUE(only.finite)) is.finite(x) else rep(TRUE, length(x))
  x[finite] <- pmin(pmax(x[finite], range[1]), range[2])
  x
}

resolve_input_path <- function(path) {
  if (grepl("^/", path)) return(path.expand(path))
  file.path(project_dir, path)
}

safe_write_guard <- function(path, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing file: ", path,
      "\nSet overwrite <- TRUE only if you intentionally want to regenerate this output.",
      call. = FALSE
    )
  }
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

safe_save_rds <- function(object, path) {
  safe_write_guard(path, overwrite = overwrite)
  saveRDS(object, path)
  invisible(path)
}

validate_parameter_contract <- function(contract) {
  expected <- list(
    blood_platelet_cluster = true_platelet_cluster,
    nFeature_RNA = paste0(plt_nFeature_min, "-", plt_nFeature_max),
    nCount_RNA = paste0(plt_nCount_min, "-", plt_nCount_max),
    percent_mt = paste0("<=", plt_percent_mt_max),
    pct_cells_expr = paste0(">=", min_pct_expr),
    mean_counts = paste0(">", min_mean_count)
  )
  for (nm in names(expected)) {
    if (!identical(contract[[nm]], expected[[nm]])) {
      stop("Internal parameter contract mismatch for ", nm, call. = FALSE)
    }
  }
  invisible(TRUE)
}

keep_present <- function(obj, genes) intersect(genes, rownames(obj))

safe_find_variable_features <- function(obj, n1 = 2000, n2 = 1000) {
  warning_text <- NULL
  obj2 <- withCallingHandlers(
    Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = n1),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warning_text) > 0) {
    cat("FindVariableFeatures warning; retrying with nfeatures=", n2, "\n", sep = "")
    obj2 <- Seurat::FindVariableFeatures(obj2, selection.method = "vst", nfeatures = n2)
  }
  obj2
}

join_layers_if_available <- function(obj) {
  if (exists("JoinLayers", mode = "function")) return(JoinLayers(obj))
  obj
}

get_counts_matrix <- function(obj) {
  out <- try(SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"), silent = TRUE)
  if (inherits(out, "try-error")) {
    out <- SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts")
  }
  out
}

curate_gene_universe <- function(df_stats) {
  df_stats$flag_MT <- grepl("^MT-", df_stats$gene)
  df_stats$flag_Ribo <- grepl("^RPL|^RPS", df_stats$gene)
  df_stats$flag_MALAT1 <- df_stats$gene == "MALAT1"
  df_stats$flag_HLA <- grepl("^HLA-", df_stats$gene)
  df_stats$flag_HIST <- grepl("^HIST", df_stats$gene) | grepl("^H[1-4]F", df_stats$gene)
  subset(df_stats, !flag_MT & !flag_Ribo & !flag_HLA & !flag_HIST)
}

if (!file.exists(manifest_path)) stop("Manifest not found: ", manifest_path)
validate_parameter_contract(parameter_contract)
manifest <- read.delim(manifest_path, stringsAsFactors = FALSE)
required_cols <- c("donor_id", "h5_path", "file_name")
if (!all(required_cols %in% colnames(manifest))) {
  stop("Manifest must contain columns: ", paste(required_cols, collapse = ", "))
}
manifest$h5_path_resolved <- vapply(manifest$h5_path, resolve_input_path, character(1))
missing_files <- manifest$h5_path_resolved[!file.exists(manifest$h5_path_resolved)]
if (length(missing_files) > 0) stop("Missing H5 inputs:\n", paste(missing_files, collapse = "\n"))

objs_all <- list()
for (i in seq_len(nrow(manifest))) {
  donor_id <- manifest$donor_id[i]
  h5_path <- manifest$h5_path_resolved[i]
  cat("\nLoading donor:", donor_id, "\n")
  counts_raw <- Seurat::Read10X_h5(h5_path)
  if (is.list(counts_raw)) {
    counts <- if ("Gene Expression" %in% names(counts_raw)) counts_raw[["Gene Expression"]] else counts_raw[[1]]
  } else {
    counts <- counts_raw
  }
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = donor_id,
    assay = "RNA",
    min.cells = min_cells_gene,
    min.features = min_features_all
  )
  obj$donor_id <- donor_id
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  objs_all[[donor_id]] <- Seurat::RenameCells(obj, add.cell.id = donor_id)
  cat("Cells:", ncol(obj), " Genes:", nrow(obj), "\n")
}

donor_ids <- names(objs_all)
obj_global_all <- objs_all[[1]]
if (length(objs_all) > 1) {
  for (j in 2:length(objs_all)) obj_global_all <- merge(obj_global_all, y = objs_all[[j]])
}
cat("Merged cells:", ncol(obj_global_all), " Genes:", nrow(obj_global_all), "\n")
safe_save_rds(obj_global_all, file.path(objects_dir, paste0(prefix, "_obj_global_all_raw.rds")))

obj_global_platelet <- subset(
  obj_global_all,
  subset = nFeature_RNA >= plt_nFeature_min &
    nFeature_RNA <= plt_nFeature_max &
    nCount_RNA >= plt_nCount_min &
    nCount_RNA <= plt_nCount_max &
    percent.mt <= plt_percent_mt_max
)
if (ncol(obj_global_platelet) < 200) stop("Too few cells after platelet-friendly filter.")

obj_global_platelet <- Seurat::NormalizeData(obj_global_platelet, normalization.method = "LogNormalize", scale.factor = 10000)
obj_global_platelet <- safe_find_variable_features(obj_global_platelet)
obj_global_platelet <- Seurat::ScaleData(obj_global_platelet, vars.to.regress = "percent.mt")
obj_global_platelet <- Seurat::RunPCA(obj_global_platelet, features = Seurat::VariableFeatures(obj_global_platelet), verbose = FALSE)

max_pcs_allowed <- min(n_pcs, ncol(obj_global_platelet) - 1)
if (max_pcs_allowed < 5) stop("Too few cells for PCA/UMAP.")
obj_global_platelet <- Seurat::FindNeighbors(obj_global_platelet, dims = 1:max_pcs_allowed)
obj_global_platelet <- Seurat::FindClusters(obj_global_platelet, resolution = resolution)
obj_global_platelet <- Seurat::RunUMAP(obj_global_platelet, dims = 1:max_pcs_allowed)

cluster_levels <- sort_cluster_levels(unique(obj_global_platelet$seurat_clusters))
obj_global_platelet$seurat_clusters <- factor(as.character(obj_global_platelet$seurat_clusters), levels = cluster_levels)
Seurat::Idents(obj_global_platelet) <- "seurat_clusters"

p_umap <- Seurat::DimPlot(
  obj_global_platelet,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 0.5,
  raster = FALSE
) +
  ggplot2::scale_color_manual(
    values = stats::setNames(cluster_palette(length(cluster_levels)), cluster_levels)
  ) +
  ggplot2::ggtitle("Blood UMAP clustering") +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") +
  ggplot2::scale_x_continuous(breaks = seq(-10, 25, by = 5)) +
  ggplot2::scale_y_continuous(breaks = seq(-10, 10, by = 5)) +
  ggplot2::labs(color = NULL) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = "Helvetica",
      size = 5,
      face = "plain",
      margin = ggplot2::margin(0, 0, 2, 0)
    ),
    axis.title = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    axis.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.position = "right",
    legend.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.key.height = grid::unit(0.13, "cm"),
    legend.key.width  = grid::unit(0.13, "cm"),
    legend.spacing.y  = grid::unit(-0.08, "cm"),
    legend.margin = ggplot2::margin(-4, 0, 0, 2),
    legend.box.margin = ggplot2::margin(-4, 0, 0, 0)
  )

p_umap <- set_point_style(p_umap, size = 0.35, alpha = 0.45)
save_plot_pdf(p_umap, file.path(figures_dir, "BLOOD_UMAP_GLOBAL_clusters.pdf"), width = 5.7 / 2.54, height = 4.5 / 2.54)

donor_levels <- sort(unique(as.character(obj_global_platelet$donor_id)))
obj_global_platelet$donor_id <- factor(
  as.character(obj_global_platelet$donor_id),
  levels = donor_levels
)

donor_cols <- stats::setNames(
  c(
    "#73CC80",
    "#3C93FA",
    "#028183",
    "#F24F4F",
    "#FF8F5C",
    "#EEB72B",
    "#FF57B0",
    "#9151B8"
  )[seq_along(donor_levels)],
  donor_levels
)



p_umap_donor <- Seurat::DimPlot(
  obj_global_platelet,
  reduction = "umap",
  group.by = "donor_id",
  label = FALSE,
  pt.size = 0.5,
  raster = FALSE
) +
  ggplot2::scale_color_manual(values = donor_cols) +
  ggplot2::ggtitle("Blood UMAP by donor") +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") +
  ggplot2::scale_x_continuous(breaks = seq(-10, 25, by = 5)) +
  ggplot2::scale_y_continuous(breaks = seq(-10, 10, by = 5)) +
  ggplot2::labs(color = NULL) +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = "Helvetica",
      size = 5,
      face = "plain",
      margin = ggplot2::margin(0, 0, 2, 0)
    ),
    axis.title = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    axis.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.position = "right",
    legend.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.key.height = grid::unit(0.13, "cm"),
    legend.key.width  = grid::unit(0.13, "cm"),
    legend.spacing.y  = grid::unit(-0.08, "cm"),
    legend.margin = ggplot2::margin(-4, 0, 0, -2),
    legend.box.margin = ggplot2::margin(-4, -4, 0, 0)
  )

p_umap_donor <- set_point_style(p_umap_donor, size = 0.35, alpha = 0.55)

save_plot_pdf(
  p_umap_donor,
  file.path(figures_dir, "BLOOD_UMAP_GLOBAL_byDonor.pdf"),
  width = 5.7 / 2.54,
  height = 4.5 / 2.54
)



p_qc <- Seurat::VlnPlot(obj_global_platelet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  ggplot2::ggtitle("Blood platelet-friendly QC") + theme_pub()
save_plot_pdf(p_qc, file.path(figures_dir, "BLOOD_QC_GLOBAL_platelet_filtered.pdf"), width = 10, height = 4)

platelet_panel_present <- keep_present(obj_global_platelet, markers_platelet)
if (length(platelet_panel_present) == 0) stop("No platelet marker genes found in object.")

lineage_marker_list <- list(
  Platelet = platelet_panel_present,
  Myeloid = keep_present(obj_global_platelet, markers_myeloid),
  Lymphoid = keep_present(obj_global_platelet, markers_lymphoid),
  Erythroid = keep_present(obj_global_platelet, markers_erythro)
)
lineage_marker_list <- lineage_marker_list[vapply(lineage_marker_list, length, integer(1)) > 0]
cat(
  "Lineage markers present (Platelet/Myeloid/Lymphoid/Erythroid): ",
  length(lineage_marker_list[["Platelet"]]), "/",
  length(lineage_marker_list[["Myeloid"]]), "/",
  length(lineage_marker_list[["Lymphoid"]]), "/",
  length(lineage_marker_list[["Erythroid"]]), "\n",
  sep = ""
)

p_dot_base <- Seurat::DotPlot(obj_global_platelet, features = lineage_marker_list, assay = "RNA")
dot_data <- p_dot_base$data
if (!("feature.groups" %in% colnames(dot_data))) dot_data$feature.groups <- "Markers"
dot_qlo <- unname(stats::quantile(dot_data$avg.exp.scaled, 0.05, na.rm = TRUE))
dot_qhi <- unname(stats::quantile(dot_data$avg.exp.scaled, 0.95, na.rm = TRUE))
dot_data$avg_clip <- pmin(pmax(dot_data$avg.exp.scaled, dot_qlo), dot_qhi)
dot_data$pct_plot <- ifelse(dot_data$pct.exp < 15, dot_data$pct.exp * 0.15, dot_data$pct.exp)


p_dot_lineage <- ggplot2::ggplot(
  dot_data,
  ggplot2::aes(
    x = features.plot,
    y = id,
    size = pct_plot,
    color = avg_clip
  )
) +
  ggplot2::geom_point(
    alpha = 0.95,
    shape = 16,
    stroke = 0
  ) +
  ggplot2::facet_grid(
    . ~ feature.groups,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_y_discrete(
    limits = rev(cluster_levels)
  ) +
  ggplot2::scale_size(
    range = c(0.01, 2.8),
    breaks = c(0, 25, 50, 75),
    name = "Percent\nexpressed"
  ) +
  ggplot2::scale_color_gradient(
    low = dot_low_col,
    high = dot_high_col,
    name = "Average\nexpression"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      ticks = FALSE,
      frame.colour = NA,
      barheight = grid::unit(1.20, "cm"),
      barwidth = grid::unit(0.22, "cm")
    ),
    size = ggplot2::guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(
        color = "black",
        alpha = 1
      )
    )
  ) +
  ggplot2::ggtitle("Lineage marker expression across blood clusters") +
  ggplot2::xlab(NULL) +
  ggplot2::ylab("Cluster") +
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = "Helvetica",
      size = 5,
      face = "plain",
      margin = ggplot2::margin(0, 0, 3, 0)
    ),
    axis.title.y = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    axis.text.y = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    axis.text.x = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    strip.background = ggplot2::element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.25
    ),
    strip.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    panel.spacing.x = grid::unit(0.18, "cm"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.key.height = grid::unit(0.22, "cm"),
    legend.key.width = grid::unit(0.22, "cm"),
    legend.ticks = ggplot2::element_blank(),
    legend.spacing.y = grid::unit(0.04, "cm"),
    legend.box.spacing = grid::unit(0.12, "cm"),
    legend.margin = ggplot2::margin(0, 0, 0, 2),
    legend.box.margin = ggplot2::margin(0, 0, 0, 2)
  )




save_plot_pdf(
  p_dot_lineage,
  file.path(figures_dir, "BLOOD_DotPlot_GLOBAL_markers_HYBRID.pdf"),
  width = 10 / 2.54,
  height = 5.5 / 2.54
)

dp <- Seurat::DotPlot(obj_global_platelet, features = platelet_panel_present)$data
genes_detected_by_cluster <- aggregate(pct.exp ~ id, data = dp, FUN = function(x) sum(x >= threshold_pct))
colnames(genes_detected_by_cluster) <- c("cluster", paste0("n_platelet_genes_pct>=", threshold_pct))
genes_detected_by_cluster <- genes_detected_by_cluster[order(-genes_detected_by_cluster[, 2]), ]
safe_write_csv(genes_detected_by_cluster, file.path(tables_dir, "BLOOD_platelet_genes_detected_by_cluster.csv"), row.names = FALSE)

obj_score_plot <- obj_global_platelet

score_panels <- list(
  PlateletMarkerModule = platelet_panel_present,
  MyeloidMarkerModule = keep_present(obj_global_platelet, c(
    "LYZ", "S100A8", "S100A9", "CTSS", "LGALS3", "FCER1G", "TYMP", "MNDA",
    "MS4A7", "CST3", "AIF1", "LST1", "FCGR3A", "TYROBP"
  )),
  ErythroidMarkerModule = keep_present(obj_global_platelet, c("HBB", "HBA1", "HBA2", "ALAS2", "CA1")),
  LymphoidMarkerModule = keep_present(obj_global_platelet, c(
    "CD3D", "CD3E", "TRAC", "IL7R", "LTB", "MS4A1", "CD79A", "CD74", "JCHAIN"
  ))
)
score_min_genes <- c(
  PlateletMarkerModule = 1,
  MyeloidMarkerModule = 5,
  ErythroidMarkerModule = 3,
  LymphoidMarkerModule = 5
)

for (score_name in names(score_panels)) {
  score_genes <- score_panels[[score_name]]
  if (length(score_genes) >= score_min_genes[[score_name]]) {
    obj_score_plot <- Seurat::AddModuleScore(
      object = obj_score_plot,
      features = list(score_genes),
      assay = "RNA",
      name = score_name
    )
  } else {
    cat("Skipping diagnostic ", score_name, ": too few genes present (", length(score_genes), ").\n", sep = "")
  }
}

score_cols <- vapply(names(score_panels), function(score_name) {
  get_last_score_col(obj_score_plot, score_name)
}, character(1), USE.NAMES = FALSE)
score_cols <- score_cols[!is.na(score_cols) & score_cols %in% colnames(obj_score_plot@meta.data)]
if (length(score_cols) >= 2) {
  score_values <- unlist(obj_score_plot@meta.data[, score_cols, drop = FALSE])
  score_lims <- as.numeric(stats::quantile(score_values, probs = c(0.01, 0.99), na.rm = TRUE))
  if (!all(is.finite(score_lims)) || score_lims[1] == score_lims[2]) {
    score_lims <- range(score_values, na.rm = TRUE)
  }
  names(score_lims) <- c("p1", "p99")
  score_breaks <- pretty(score_lims, n = 4)
  score_breaks <- score_breaks[score_breaks >= score_lims["p1"] & score_breaks <= score_lims["p99"]]
  if (0 >= score_lims["p1"] && 0 <= score_lims["p99"] && !any(abs(score_breaks) < 1e-9)) {
    score_breaks <- sort(c(score_breaks, 0))
  }
  blue_palette <- c("#EAEFEF", "#BBE0EF", "#5A82E8", "#123A9C", "#071B4A")
  p_feat4 <- Seurat::FeaturePlot(
    object = obj_score_plot,
    features = score_cols,
    ncol = 2,
    cols = blue_palette,
    min.cutoff = score_lims["p1"],
    max.cutoff = score_lims["p99"],
    order = TRUE,
    pt.size = 0.10,
    combine = FALSE
  )
  diagnostic_titles <- c(
    PlateletMarkerModule = "Platelet marker module",
    MyeloidMarkerModule = "Myeloid marker module",
    ErythroidMarkerModule = "Erythroid marker module",
    LymphoidMarkerModule = "Lymphoid marker module"
  )
  p_feat4 <- lapply(seq_along(p_feat4), function(i) {
    score_prefix <- sub("[0-9]+$", "", score_cols[[i]])
    p <- p_feat4[[i]] +
      ggplot2::ggtitle(diagnostic_titles[[score_prefix]]) +
      theme_pub() +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain"),
        legend.key.height = grid::unit(0.24, "cm"),
        legend.key.width = grid::unit(0.16, "cm")
      ) +
      ggplot2::guides(color = ggplot2::guide_colorbar(
        ticks = TRUE,
        frame.colour = NA,
        barheight = grid::unit(1.1, "cm"),
        barwidth = grid::unit(0.20, "cm")
      )) +
      ggplot2::scale_color_gradientn(
        colours = blue_palette,
        limits = c(score_lims["p1"], score_lims["p99"]),
        breaks = score_breaks,
        oob = squish_oob
      )
    set_point_style(p, size = 0.4, alpha = 0.6)
  })
  save_plot_grid_pdf(
    p_feat4,
    file.path(figures_dir, "Blood_marker_module_featureplots.pdf"),
    width = 10 / 2.54,
    height = 8 / 2.54,
    ncol = 2
  )
} else {
  cat("Skipping diagnostic FeaturePlot 4 marker modules: fewer than two diagnostic marker modules were generated.\n")
}



platelet_marker_module_col <- get_last_score_col(obj_score_plot, "PlateletMarkerModule")

if (!is.na(platelet_marker_module_col) && platelet_marker_module_col %in% colnames(obj_score_plot@meta.data)) {
  
  vln_cluster_cols <- stats::setNames(
    cluster_palette(length(cluster_levels)),
    cluster_levels
  )
  
  p_vln_platelet <- Seurat::VlnPlot(
    obj_score_plot,
    features = platelet_marker_module_col,
    group.by = "seurat_clusters",
    pt.size = 0,
    cols = vln_cluster_cols
  ) +
    ggplot2::ggtitle("Platelet marker module by cluster") +
    ggplot2::xlab("Cluster") +
    ggplot2::ylab("Marker module signal") +
    theme_pub() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "Helvetica",
        size = 5,
        face = "plain",
        margin = ggplot2::margin(0, 0, 4, 0)
      ),
      axis.title = ggplot2::element_text(
        family = "Helvetica",
        size = 4.5,
        face = "plain"
      ),
      axis.text.x = ggplot2::element_text(
        family = "Helvetica",
        size = 4.5,
        face = "plain",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica",
        size = 4.5,
        face = "plain"
      ),
      legend.position = "none"
    )
  
  if (length(p_vln_platelet$layers) >= 1) {
    p_vln_platelet$layers[[1]]$aes_params$linewidth <- 0.18
    p_vln_platelet$layers[[1]]$aes_params$colour <- "grey20"
    p_vln_platelet$layers[[1]]$aes_params$alpha <- 1
  }
  
  save_plot_pdf(
    p_vln_platelet,
    file.path(figures_dir, "BLOOD_SUPP_Violin_PlateletMarkerModule_byCluster.pdf"),
    width = 6 / 2.54,
    height = 4.5 / 2.54
  )
  
} else {
  cat("Skipping diagnostic platelet marker module violin: PlateletMarkerModule was not generated.\n")
}


clusters_exist <- levels(Seurat::Idents(obj_global_platelet))
if (!(true_platelet_cluster %in% clusters_exist)) {
  stop("Configured true platelet cluster ", true_platelet_cluster, " is absent. Available: ", paste(clusters_exist, collapse = ", "))
}


cells_cluster3 <- colnames(obj_global_platelet)[
  as.character(obj_global_platelet$seurat_clusters) == true_platelet_cluster
]

p_umap_true <- Seurat::DimPlot(
  obj_global_platelet,
  reduction = "umap",
  cells.highlight = cells_cluster3,
  cols = context_cell_col,
  cols.highlight = platelet_focus_col,
  sizes.highlight = 0.5,
  pt.size = 0.5,
  label = FALSE,
  raster = FALSE
) +
  ggplot2::ggtitle("Blood low-RNA UMAP") +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") +
  ggplot2::scale_x_continuous(breaks = seq(-10, 25, by = 5)) +
  ggplot2::scale_y_continuous(breaks = seq(-10, 10, by = 5)) +
  ggplot2::annotate(
    "text",
    x = 11.2,
    y = 4.5,
    label = "Platelet-associated\ncluster",
    hjust = 0.5,
    vjust = 1,
    family = "Helvetica",
    size = 1.582,
    fontface = "plain",
    color = "black",
    lineheight = 0.85
  ) + 
  theme_pub() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = "Helvetica",
      size = 5,
      face = "plain",
      margin = ggplot2::margin(0, 0, 2, 0)
    ),
    axis.title = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    axis.text = ggplot2::element_text(
      family = "Helvetica",
      size = 4.5,
      face = "plain"
    ),
    legend.position = "none"
  )

p_umap_true <- set_point_style(p_umap_true, size = 0.35, alpha = 0.45)

save_plot_pdf(
  p_umap_true,
  file.path(figures_dir, "Blood_platelet_associated_lowRNA_cluster_UMAP.pdf"),
  width = 5.1 / 2.54,
  height = 4.5 / 2.54
)



obj_true_plt <- subset(obj_global_platelet, idents = true_platelet_cluster)
obj_true_plt <- join_layers_if_available(obj_true_plt)
counts_true <- get_counts_matrix(obj_true_plt)
pct_cells_expr <- Matrix::rowMeans(counts_true > 0) * 100
mean_counts <- Matrix::rowMeans(counts_true)
gene_stats <- data.frame(
  gene = rownames(counts_true),
  pct_cells_expr = as.numeric(pct_cells_expr),
  mean_counts = as.numeric(mean_counts),
  stringsAsFactors = FALSE
)

gene_universe_raw <- subset(gene_stats, pct_cells_expr >= min_pct_expr & mean_counts > min_mean_count)
gene_universe_curated <- curate_gene_universe(gene_universe_raw)
gene_universe_raw <- gene_universe_raw[order(-gene_universe_raw$pct_cells_expr, -gene_universe_raw$mean_counts), ]
gene_universe_curated <- gene_universe_curated[order(-gene_universe_curated$pct_cells_expr, -gene_universe_curated$mean_counts), ]

safe_write_csv(gene_universe_raw, file.path(tables_dir, "Blood_platelet_associated_lowRNA_cluster_gene_universe_raw.csv"), row.names = FALSE)
safe_write_csv(gene_universe_curated, file.path(tables_dir, "Blood_platelet_associated_lowRNA_cluster_gene_universe_curated.csv"), row.names = FALSE)
safe_write_lines(gene_universe_raw$gene, file.path(tables_dir, "Blood_platelet_associated_lowRNA_cluster_genes_raw.txt"))
safe_write_lines(gene_universe_curated$gene, file.path(tables_dir, "Blood_platelet_associated_lowRNA_cluster_genes_curated.txt"))

obj_global_all <- join_layers_if_available(obj_global_all)
obj_global_all <- Seurat::NormalizeData(obj_global_all, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
genes_cur_present <- intersect(gene_universe_curated$gene, rownames(obj_global_all))
if (length(genes_cur_present) < 20) stop("Too few curated platelet-associated genes present for scoring.")
obj_global_all <- Seurat::AddModuleScore(
  object = obj_global_all,
  features = list(genes_cur_present),
  assay = "RNA",
  name = "platelet_associated_transcriptional_score_Blood_universe"
)
score_col <- "platelet_associated_transcriptional_score_Blood_universe1"
if (!(score_col %in% colnames(obj_global_all@meta.data))) stop("Score column not generated: ", score_col)
safe_write_csv(obj_global_all@meta.data, file.path(tables_dir, "Blood_global_all_metadata_with_platelet_associated_transcriptional_score.csv"), row.names = TRUE)
safe_save_rds(obj_global_platelet, file.path(objects_dir, "BLOOD_obj_global_platelet.rds"))
safe_save_rds(obj_global_all, file.path(objects_dir, "Blood_obj_global_all_with_platelet_associated_transcriptional_score.rds"))

p_vln_score_donor <- Seurat::VlnPlot(
  obj_global_all,
  features = score_col,
  group.by = "donor_id",
  pt.size = 0
) +
  ggplot2::ggtitle("Blood universe platelet score by donor") +
  ggplot2::xlab(NULL) +
  ggplot2::ylab("Module score") +
  theme_pub() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(family = "Helvetica", size = 4.5, face = "plain", angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )
save_plot_pdf(p_vln_score_donor, file.path(figures_dir, "Blood_platelet_associated_transcriptional_score_by_donor.pdf"), width = 5.2, height = 3.5)

score_summary <- data.frame(
  metric = c("min", "p25", "median", "mean", "p75", "max"),
  value = as.numeric(c(
    min(obj_global_all@meta.data[[score_col]], na.rm = TRUE),
    stats::quantile(obj_global_all@meta.data[[score_col]], 0.25, na.rm = TRUE),
    stats::median(obj_global_all@meta.data[[score_col]], na.rm = TRUE),
    mean(obj_global_all@meta.data[[score_col]], na.rm = TRUE),
    stats::quantile(obj_global_all@meta.data[[score_col]], 0.75, na.rm = TRUE),
    max(obj_global_all@meta.data[[score_col]], na.rm = TRUE)
  ))
)
safe_write_csv(score_summary, file.path(tables_dir, "Blood_platelet_associated_transcriptional_score_summary.csv"), row.names = FALSE)

cat("Completed Blood reference build.\n")
close_log_and_write_session()
