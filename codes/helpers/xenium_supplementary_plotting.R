# S1-S3: assemble the eleven non-CXCL12-CXCR4 pairs in native LR rank.
# Called by 01_plot_xenium_results_Figure2_FigureS1_S3.R after reading the BLISA RDS.
# Layout reference: manuscript_supp_figures/xenium_{hotspot_pvalues,
# cci_heatmaps,spatial_cci}.pdf. No inference is performed here.
supp_title <- function(id) sprintf("%s (%s significant hotspots)",
  gsub("_", "\u2013", id), format(lr$sig_numbers[match(id, lr_ids)], big.mark = ","))
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
  index <- match(id, lr_ids)
  bins <- res$bins
  # Display the saved native p-values and hotspot indices. The weight mask
  # identifies bins included in the corresponding nearby/distant LRI test.
  weight_key <- if (lr$ccc_mode[index] == "nearby") "keep_idx_queen" else "keep_idx_dist"
  tested <- seq_len(nrow(bins)) %in% res$spatial_weights[[weight_key]] &
    is.finite(lr$all_pval[[index]]) & !is.na(lr$all_quadrant[[index]])
  bins$.background <- ifelse(tested, "#D3D3D3", "#FFFFFF")
  bins$.p <- pmin(4, -log10(pmax(lr$all_pval[[index]], 1e-4)))
  p <- ggplot(bins) + geom_sf(fill = bins$.background, colour = NA) +
    geom_sf(data = bins[lr$sig_index[[index]], ], aes(fill = .p), colour = NA) +
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
  p <- dominant_map(id, roi_linewidth = 0.3) +
    guides(fill = guide_legend(title = "Sender \u2192 receiver")) +
    labs(title = supp_title(id)) + theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 5.8), legend.title = element_text(size = 6.5, face = "bold"),
      legend.key.height = grid::unit(1.8, "mm"), legend.key.width = grid::unit(1.8, "mm"),
      plot.margin = margin(4, 4, 4, 4))
  # Keep the supplementary legend wording while using native pair assignments.
  fill_scale <- p$scales$get_scales("fill")
  fill_scale$labels <- function(x) {
    labels <- gsub(" -> ", " \u2192 ", cellLabel(x), fixed = TRUE)
    labels[x == "rare pairs"] <- "Rare pairs"
    labels[x == "Non-Significant"] <- "Non-significant"
    labels[x == "Empty"] <- "Untested"
    labels
  }
  p
})
dominant_legend <- grid::grid.grabExpr({
  ComplexHeatmap::draw(roi_legend, y = grid::unit(0.68, "npc"))
  grid::grid.text(paste("Top eight dominant cell-type pairs per LR",
    "Remaining pairs are shown as rare pairs", "Colours are defined separately in each panel", sep = "\n"),
    y = 0.18, gp = grid::gpar(fontsize = 8, lineheight = 1.6))
}, width = 4.3, height = 2.25)
save_panel(supp_layout(dominant_panels, dominant_legend,
  "Xenium | Spatial CCI by LR pair"), "figureS3", 12.9, 9.24)
