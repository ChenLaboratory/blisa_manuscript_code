# Run from the repository root. --main-only renders Figure 2; --supp-only S1-S3.
# Layout source: run_blisa_on_xenium.pdf pages 3, 8, 14, 22 and 27.
# Plotting reads processed inputs and complete CSV results; it never runs inference.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "helpers", "plotting_helpers.R"))
source(repo_path("codes", "helpers", "manuscript_plotting_helpers.R"))
require_packages(c("blisa", "scider", "SpatialExperiment", "S4Vectors", "Matrix",
  "sf", "ggplot2", "ggnewscale", "RColorBrewer", "ComplexHeatmap", "circlize",
  "viridisLite", "gridExtra", "scales"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(Matrix)
})
binned <- readRDS(repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds"))
read_result <- function(name) read.csv(repo_path("res", "xenium", name),
  check.names = FALSE, stringsAsFactors = FALSE)
summary <- read_result("lri_summary.csv")
local <- read_result("lri_by_bin.csv")
cci <- read_result("cci_scores.csv")
dominant <- read_result("xenium_dominant_celltype_pairs_by_bin.csv")
lr_ids <- summary$lr_id
selected <- "CXCL12_CXCR4"
stopifnot(selected %in% lr_ids, !anyDuplicated(lr_ids))
other_lr <- setdiff(lr_ids, selected)
if (is.null(binned$cells_metadata)) {
  stop("Run codes/analysis/xenium/01_prepare_xenium_50um_bins.R first: ",
       "cell-level metadata are required for Figure 2a.")
}
save_panel <- save_figure_panel
roi_colours <- c(Invasive_Tumor = "#377EB8", DCIS_1 = "#4DAF4A", DCIS_2 = "#006D2C")
aligned_rows <- function(table, id) {
  rows <- table[table$lr_id == id, , drop = FALSE]
  index <- match(binned$bins$bin_id, rows$bin_id)
  stopifnot(nrow(rows) == nrow(binned$bins), !anyDuplicated(rows$bin_id), !anyNA(index))
  rows[index, , drop = FALSE]
}

# Reassemble the native plotting object from complete CSVs. This is formatting,
# not a second significance calculation or an inference-result RDS dependency.
lr <- summary[, setdiff(names(summary), c("rank", "lr_id", "n_tested", "n_total_bins")),
  drop = FALSE]
rownames(lr) <- lr_ids
rows_by_lr <- lapply(lr_ids, function(id) aligned_rows(local, id))
lr$sig_index <- lapply(rows_by_lr, function(x) which(x$native_significant))
lr$sig_pval <- lapply(rows_by_lr, function(x) x$p_value[x$native_significant])
lr$all_pval <- lapply(rows_by_lr, function(x) x$p_value)
lr$all_lisa <- lapply(rows_by_lr, function(x) x$score)
lr$all_quadrant <- lapply(rows_by_lr, function(x) x$quadrant)
res <- structure(list(LR_results = lr, bins = binned$bins,
  spatial_weights = binned$spatial_weights, CCI_scores = cci, level = "lr"), class = "blisa")

# Same expression palette, titles and bottom colour bars as lrExpressionMap()
# in the source report. Expression panels do not have pathology outlines.
expression_map <- function(partner, role) {
  genes <- strsplit(gsub("\\s+", "", partner), "[,_]")[[1]]
  stopifnot(all(genes %in% rownames(binned$counts_matrix)))
  index <- match(paste0("bin_factor", binned$bins$bin_id), colnames(binned$counts_matrix))
  stopifnot(!anyNA(index))
  values <- binned$counts_matrix[genes, index, drop = FALSE]
  bins <- binned$bins
  bins$.bin_value <- log1p(if (length(genes) == 1L) as.numeric(values) else
    apply(as.matrix(values), 2, min))
  ggplot(bins) + geom_sf(aes(fill = .bin_value), colour = NA) +
    scale_fill_gradientn(colours = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
      "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F"),
      breaks = scales::breaks_pretty(n = 3), name = "log1p(count)") +
    labs(title = paste0(role, ": ", partner)) + theme_void() +
    spatialTheme() + guides(fill = compactColourbar())
}
hotspot_map <- function(id) {
  plotHotspotsWithPaperROIs(res, match(id, lr_ids), binned$paper_rois, roi_colours)
}
cci_heatmap <- function(id, limits = NULL) {
  values <- cciMatrix(cci, cci[[id]])
  # Supplementary heatmaps retain the agreed common score scale.
  colours <- if (is.null(limits)) viridisLite::viridis(10) else
    circlize::colorRamp2(seq(limits[1], limits[2], length.out = 10),
      viridisLite::viridis(10))
  plot <- squareCCI(values, cell_mm = min(2.5, 32.5 / max(dim(values))),
    colours = colours)
  heatmapGrob(plot, width = 3.6, height = 2.55)
}
dominant_map <- function(id) {
  values <- aligned_rows(dominant, id)
  ranked <- sort(table(values$cell_pair[values$native_significant]), decreasing = TRUE)
  shown <- names(head(ranked, 8))
  bins <- binned$bins
  # Preserve native labels, category palette and top-pair order. The full
  # assignments were already exported; this does not recalculate CCI scores.
  bins$cell_pair_plot <- ifelse(values$tested, "Non-Significant", "Empty")
  hotspots <- values$native_significant
  bins$cell_pair_plot[hotspots] <- ifelse(values$cell_pair[hotspots] %in% shown,
    values$cell_pair[hotspots], "rare pairs")
  bins$cell_pair_plot <- factor(bins$cell_pair_plot,
    levels = c(shown, "rare pairs", "Non-Significant", "Empty"))
  colours <- c(setNames(getFromNamespace("cols", "blisa")[seq_along(shown)], shown),
    "rare pairs" = "#404040", "Non-Significant" = "#D3D3D3", Empty = "#FFFFFF")
  p <- ggplot(bins) + geom_sf(aes(fill = cell_pair_plot), colour = NA) +
    scale_fill_manual(values = colours, drop = FALSE, labels = cellLabel) +
    guides(fill = guide_legend(title = if (length(ranked) <= 8) "All pairs" else "Top 8 pairs")) +
    labs(title = "Spatial CCI") + theme_void() + spatialTheme(7.5, "right")
  addPaperROIBorders(p, binned$paper_rois)
}

# Figure 2a, report page 3: native scider axes, original cell labels and ROI paths.
if (!"--supp-only" %in% commandArgs(trailingOnly = TRUE)) {
cells <- binned$cells_metadata
broad_levels <- levels(cells$cell_type_broad)
if (is.null(broad_levels)) broad_levels <- unique(cells$cell_type_broad)
cell_colours <- setNames(c(RColorBrewer::brewer.pal(7, "Set1"),
  RColorBrewer::brewer.pal(6, "Set2")), broad_levels)
cell_colours <- c(cell_colours[names(cell_colours) != "DCIS"],
  DCIS_1 = unname(cell_colours["DCIS"]), DCIS_2 = "#006D2C")
cell_colours <- cell_colours[levels(cells$cell_type)]
cell_colours["Stromal"] <- "#4D4D4D"
# Only coordinates and annotations are needed to call the original cell plotter.
cell_spe <- SpatialExperiment::SpatialExperiment(
  colData = S4Vectors::DataFrame(cells, row.names = cells$cell_id),
  spatialCoords = as.matrix(cells[, c("x", "y")]))
borders <- paperROIBorderPaths(binned$paper_rois)
p_cells <- scider::plotSpatial(cell_spe, group.by = "cell_type", pt.size = 0.2,
  pt.alpha = 0.9, cols = cell_colours) +
  geom_path(data = borders, aes(x = x, y = y, group = piece), inherit.aes = FALSE,
    colour = "white", linewidth = 1.8, show.legend = FALSE) +
  geom_path(data = borders, aes(x = x, y = y, group = piece, colour = region),
    inherit.aes = FALSE, linewidth = 1.2, show.legend = FALSE) +
  spatialTheme(9, "right") + labs(subtitle = "Outlines: paper Figure 4 ROIs")
save_panel(p_cells, "figure2a", 7.2, 4.2, pointsize = 10)

# Figure 2b, report page 8: native annotation colours and native axis range.
p_rank <- blisa::plotLRrank(res, top = 30, pt_size = 2.5) +
  labs(x = "Number of hotspots", y = NULL) + theme(
    axis.text = element_text(size = 9, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 8), panel.grid.minor = element_blank(),
    plot.margin = margin(3, 3, 3, 3))
save_panel(p_rank, "figure2b", 4.5, 3.8, pointsize = 10)

# Figure 2c, report page 22: fixed physical tile sizes and coloured ROI labels.
regions <- as.character(binned$bins_metadata$region[
  match(binned$bins$bin_id, binned$bins_metadata$bin_id)])
region_names <- names(roi_colours)
coverage <- t(vapply(rows_by_lr, function(rows) {
  vapply(region_names, function(region) {
    sum(rows$native_significant & regions == region) / sum(regions == region)
  }, numeric(1))
}, numeric(length(region_names))))
dimnames(coverage) <- list(lr_ids, region_names)
save_panel(heatmapGrob(regionalHotspotHeatmap(coverage, roi_colours),
  width = 7.2, height = 4.1), "figure2c", 7.2, 4.1, pointsize = 10)

# Figure 2d: page-14 expression maps and page-27 hotspot map. Keep exact source
# subpanels for assembly, and stack receptor above ligand to the left of hotspots.
receptor <- expression_map("CXCR4", "Receptor")
ligand <- expression_map("CXCL12", "Ligand")
hotspots <- hotspot_map(selected)
save_panel(receptor, "figure2d_receptor_expression", 2.4, 2)
save_panel(ligand, "figure2d_ligand_expression", 2.4, 2)
save_panel(hotspots, "figure2d_hotspots", 7.2, 4.3)
expression_column <- gridExtra::arrangeGrob(receptor, ligand, ncol = 1)
p_expression_hotspots <- gridExtra::arrangeGrob(expression_column, hotspots,
  nrow = 1, widths = c(2.4, 7.2))
save_panel(p_expression_hotspots, "figure2d", 9.6, 4.3)

# Figure 2e-f: the CXCL12-CXCR4 heatmap/spatial CCI block on report page 14.
save_panel(cci_heatmap(selected), "figure2e", 3.6, 2.55)
save_panel(dominant_map(selected), "figure2f", 3.6, 2.55)
}

# Supplementary panels retain native LR rank, excluding the main-text example.
if (!"--main-only" %in% commandArgs(trailingOnly = TRUE)) {
  source(repo_path("codes", "helpers", "xenium_supplementary_plotting.R"))
}
