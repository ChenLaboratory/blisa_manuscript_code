# S1-S3: assemble the eleven non-CXCL12-CXCR4 pairs in native LR rank.
# Called by 01_plot_xenium_results_Figure2_FigureS1_S3.R after reading full CSVs.
# Layout reference: manuscript_supp_figures/xenium_{hotspot_pvalues,
# cci_heatmaps,spatial_cci}.pdf. No inference is performed here.
supp_title <- function(id) sprintf("%s (%s significant hotspots)",
  gsub("_", "\u2013", id), format(summary$sig_numbers[match(id, lr_ids)], big.mark = ","))
supp_layout <- function(panels, legend, title) {
  gridExtra::arrangeGrob(grobs = c(panels, list(legend)), ncol = 3,
    top = grid::textGrob(title, x = grid::unit(2, "mm"), hjust = 0,
      gp = grid::gpar(fontsize = 14, fontface = "bold")))
}
roi_legend <- ComplexHeatmap::Legend(title = "ROI boundaries",
  labels = c("Invasive tumour", "DCIS 1", "DCIS 2"), type = "lines",
  legend_gp = grid::gpar(col = unname(roi_colours), lwd = 1.2),
  title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 8))
warm_colours <- c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")
hotspot_panels <- lapply(other_lr, function(id) {
  values <- aligned_rows(local, id)
  bins <- binned$bins
  bins$.background <- ifelse(values$tested, "#D3D3D3", "#FFFFFF")
  bins$.p <- pmin(4, -log10(pmax(values$p_value, 1e-4)))
  p <- ggplot(bins) + geom_sf(fill = bins$.background, colour = NA) +
    geom_sf(data = bins[values$native_significant, ], aes(fill = .p), colour = NA) +
    scale_fill_gradientn(colours = warm_colours, limits = c(-log10(0.05), 4),
      guide = "none") + theme_void() + labs(title = supp_title(id)) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.margin = margin(4, 5, 4, 5))
  addPaperROIBorders(p, binned$paper_rois, linewidth = 0.3)
})
hotspot_legend <- grid::grid.grabExpr({
  lgd <- ComplexHeatmap::packLegend(
    ComplexHeatmap::Legend(title = "Hotspot -log10(p)",
      col_fun = circlize::colorRamp2(seq(-log10(0.05), 4, length.out = length(warm_colours)), warm_colours),
      at = c(-log10(0.05), 2, 3, 4), labels = c("1.30", "2", "3", "4"),
      legend_height = grid::unit(28, "mm"),
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8)),
    ComplexHeatmap::Legend(title = "Other bins", labels = c("Non-significant", "Untested"),
      legend_gp = grid::gpar(fill = c("#D3D3D3", "#FFFFFF"), col = NA),
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"), labels_gp = grid::gpar(fontsize = 8)),
    roi_legend, gap = grid::unit(4, "mm"))
  ComplexHeatmap::draw(lgd)
}, width = 3.9, height = 3)
save_panel(supp_layout(hotspot_panels, hotspot_legend,
  "Xenium | LR hotspot p-values"), "figureS1", 11.7, 12.4)

# The original combined heatmap uses one explicit 0-18 scale for every LR.
stopifnot(max(as.matrix(cci[, other_lr, drop = FALSE]), na.rm = TRUE) <= 18)
cci_colours <- circlize::colorRamp2(seq(0, 18, length.out = 10), viridisLite::viridis(10))
cci_panels <- lapply(other_lr, function(id) {
  mat <- cciMatrix(cci, cci[[id]])
  ht <- ComplexHeatmap::Heatmap(mat, name = "CCI score", col = cci_colours,
    show_heatmap_legend = FALSE, width = grid::unit(45, "mm"), height = grid::unit(45, "mm"),
    row_title = "Receiver", column_title = "Sender",
    row_title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
    row_labels = cellLabel(rownames(mat)), column_labels = cellLabel(colnames(mat)),
    row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 6.5),
    column_names_rot = 60, row_dend_width = grid::unit(3, "mm"),
    column_dend_height = grid::unit(3, "mm"), rect_gp = grid::gpar(col = "#FFFFFF55", lwd = 0.25))
  gridExtra::arrangeGrob(heatmapGrob(ht, 3.9, 3),
    top = grid::textGrob(supp_title(id), gp = grid::gpar(fontsize = 9, fontface = "bold")))
})
cci_legend <- grid::grid.grabExpr(ComplexHeatmap::draw(
  ComplexHeatmap::Legend(title = "CCI score", col_fun = cci_colours, at = seq(0, 18, 3),
    legend_height = grid::unit(55, "mm"), title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 8))), width = 3.9, height = 3)
save_panel(supp_layout(cci_panels, cci_legend,
  "Xenium | LR-specific CCI heatmaps"), "figureS2", 11.7, 12.76)

dominant_panels <- lapply(other_lr, function(id) {
  values <- aligned_rows(dominant, id)
  ranked <- sort(table(values$cell_pair[values$native_significant]), decreasing = TRUE)
  shown <- names(head(ranked, 8))
  bins <- binned$bins
  bins$.pair <- ifelse(values$tested, "Non-significant", "Untested")
  bins$.pair[values$native_significant] <- ifelse(values$cell_pair[values$native_significant] %in% shown,
    values$cell_pair[values$native_significant], "Rare pairs")
  bins$.pair <- factor(bins$.pair, levels = c(shown, "Rare pairs", "Non-significant", "Untested"))
  pair_rows <- values[match(shown, values$cell_pair), ]
  pair_labels <- paste(cellLabel(pair_rows$sender), "\u2192", cellLabel(pair_rows$receiver))
  colours <- c(setNames(getFromNamespace("cols", "blisa")[seq_along(shown)], shown),
    "Rare pairs" = "#404040", "Non-significant" = "#D3D3D3", Untested = "#FFFFFF")
  p <- ggplot(bins) + geom_sf(aes(fill = .pair), colour = NA) +
    scale_fill_manual(values = colours, drop = FALSE,
      labels = c(pair_labels, "Rare pairs", "Non-significant", "Untested"),
      name = "Sender \u2192 receiver") + theme_void() + spatialTheme(7.5, "right") +
    labs(title = supp_title(id)) + theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 5.8), legend.title = element_text(size = 6.5, face = "bold"),
      legend.key.height = grid::unit(1.8, "mm"), legend.key.width = grid::unit(1.8, "mm"),
      plot.margin = margin(4, 4, 4, 4))
  addPaperROIBorders(p, binned$paper_rois, linewidth = 0.3)
})
dominant_legend <- grid::grid.grabExpr({
  ComplexHeatmap::draw(roi_legend, y = grid::unit(0.68, "npc"))
  grid::grid.text(paste("Top eight dominant cell-type pairs per LR",
    "Remaining pairs are shown as rare pairs", "Colours are defined separately in each panel", sep = "\n"),
    y = 0.18, gp = grid::gpar(fontsize = 8, lineheight = 1.6))
}, width = 4.3, height = 2.25)
save_panel(supp_layout(dominant_panels, dominant_legend,
  "Xenium | Spatial CCI by LR pair"), "figureS3", 12.9, 9.24)
