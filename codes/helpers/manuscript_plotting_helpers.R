# Manuscript layout functions extracted from the original Xenium report.
# These functions only format plots; they do not prepare data or run inference.

# Source: blisa_on_xenium/plotManuscript.R
cellLabel <- function(x) gsub("_", " ", x, fixed = TRUE)

# Source: blisa_on_xenium/plotManuscript.R
spatialTheme <- function(font_size = 8, legend_position = "bottom") {
  ggplot2::theme(
    text = ggplot2::element_text(family = "sans", size = font_size),
    plot.title = ggplot2::element_text(size = font_size, face = "bold", hjust = 0),
    plot.subtitle = ggplot2::element_text(size = font_size - 1),
    legend.position = legend_position,
    legend.title = ggplot2::element_text(size = font_size - 1, face = "bold"),
    legend.text = ggplot2::element_text(size = font_size - 1),
    legend.key.height = grid::unit(2.3, "mm"),
    legend.key.width = grid::unit(2.3, "mm"),
    legend.spacing.y = grid::unit(0.3, "mm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(1, "mm"),
    plot.margin = ggplot2::margin(2, 3, 2, 3)
  )
}

# Source: blisa_on_xenium/plotManuscript.R
compactColourbar <- function() {
  ggplot2::guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(26, "mm"), barheight = grid::unit(1.8, "mm")
  )
}

# Source: blisa_on_xenium/plotManuscript.R
heatmapGrob <- function(heatmap, width, height) {
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap, heatmap_legend_side = "right", annotation_legend_side = "right",
      padding = grid::unit(rep(1, 4), "mm")
    ),
    width = width, height = height
  )
}

# Source: blisa_on_xenium/plotManuscript.R
cciMatrix <- function(scores, values) {
  mat <- matrix(NA_real_, nrow = length(unique(scores$Receiver)),
                ncol = length(unique(scores$Sender)),
                dimnames = list(unique(scores$Receiver), unique(scores$Sender)))
  mat[cbind(match(scores$Receiver, rownames(mat)),
            match(scores$Sender, colnames(mat)))] <- values
  mat
}

# Source: blisa_on_xenium/plotManuscript.R
squareCCI <- function(mat, cell_mm = 2.5, font_size = 7,
                      colours = viridisLite::viridis(10)) {
  ComplexHeatmap::Heatmap(
    mat, name = "CCI score", col = colours, na_col = "grey90",
    width = grid::unit(ncol(mat) * cell_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_mm, "mm"),
    row_title = "Receiver", column_title = "Sender",
    row_title_gp = grid::gpar(fontsize = font_size, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = font_size, fontface = "bold"),
    row_labels = cellLabel(rownames(mat)), column_labels = cellLabel(colnames(mat)),
    row_names_gp = grid::gpar(fontsize = font_size, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = font_size - 0.5, fontface = "bold"),
    column_names_rot = 60,
    row_dend_width = grid::unit(3, "mm"), column_dend_height = grid::unit(3, "mm"),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = font_size, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = font_size)
    )
  )
}

# Source: blisa_on_xenium/plotManuscript.R
addPaperROIBorders <- function(plot, rois, linewidth = 0.45) {
  if (is.null(rois)) return(plot)
  colours <- c(Invasive_Tumor = "#377EB8", DCIS_1 = "#4DAF4A", DCIS_2 = "#006D2C")
  stopifnot(all(as.character(rois$region) %in% names(colours)))
  plot <- plot + ggplot2::geom_sf(
    data = rois, inherit.aes = FALSE, fill = NA, colour = "white",
    linewidth = linewidth + 0.3, show.legend = FALSE
  )
  for (region in names(colours)) {
    plot <- plot + ggplot2::geom_sf(
      data = rois[as.character(rois$region) == region, ], inherit.aes = FALSE,
      fill = NA, colour = colours[[region]], linewidth = linewidth, show.legend = FALSE
    )
  }
  plot
}

# Source: blisa_on_xenium/regionSummary.R
paperROIBorderPaths <- function(rois) {
  do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
    xy <- sf::st_coordinates(sf::st_geometry(rois)[i])
    rings <- xy[, startsWith(colnames(xy), "L"), drop = FALSE]
    data.frame(
      x = xy[, "X"], y = xy[, "Y"], region = as.character(rois$region[i]),
      # Keep disconnected polygons and interior rings as separate paths.
      piece = paste(i, apply(rings, 1, paste, collapse = "_"), sep = "_")
    )
  }))
}

# Source: blisa_on_xenium/regionSummary.R
plotHotspotsWithPaperROIs <- function(res, index, rois, colours) {
  region_names <- as.character(rois$region)
  if (!all(region_names %in% names(colours)) || anyNA(colours[region_names])) {
    stop("Paper ROI colours are incomplete.")
  }
  lr_name <- rownames(res$LR_results)[index]
  n_hotspots <- length(res$LR_results$sig_index[[index]])

  # Retain BLISA's hotspot values, bin colours, and untested-bin handling.
  p <- blisa::plotHotspots(res, index = index) +
    ggplot2::labs(title = sprintf("%s (%d hotspots)", lr_name, n_hotspots),
                  subtitle = NULL, colour = "-log10(p)") +
    ggplot2::guides(colour = ggplot2::guide_colorbar(
      title.position = "top", order = 1,
      # ggnewscale renames this aesthetic when adding the ROI colour scale.
      available_aes = "any",
      barheight = grid::unit(20, "mm"), barwidth = grid::unit(3, "mm")
    ))
  pvalue_scale <- p$scales$get_scales("colour")
  pvalue_scale$breaks <- scales::breaks_pretty(n = 3)

  # Separate the p-value scale from the categorical ROI outline colours.
  p + ggnewscale::new_scale_colour() +
    ggplot2::geom_sf(
      data = rois, ggplot2::aes(colour = region), inherit.aes = FALSE,
      fill = NA, linewidth = 1.2, show.legend = "line"
    ) +
    ggplot2::scale_colour_manual(
      values = colours[region_names], breaks = region_names,
      labels = cellLabel, name = "Paper ROI",
      guide = ggplot2::guide_legend(order = 2)
    ) + spatialTheme(9, "right") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
                    legend.key.width = grid::unit(5, "mm"))
}

# Source: blisa_on_xenium/regionSummary.R
regionalHotspotHeatmap <- function(mat, colours, column_title = "Paper ROI") {
  if (!nrow(mat) || !ncol(mat)) stop("No LR-region proportions to plot.")
  label_colours <- unname(colours[colnames(mat)])
  if (anyNA(label_colours)) stop("Region label colours are incomplete.")
  ComplexHeatmap::Heatmap(
    mat, name = "Hotspot\nproportion",
    col = circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                              c("white", "#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")),
    na_col = "grey80", cluster_rows = FALSE, cluster_columns = FALSE,
    width = grid::unit(ncol(mat) * 10, "mm"),
    height = grid::unit(nrow(mat) * 4, "mm"),
    row_names_side = "left", row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    column_labels = cellLabel(colnames(mat)), column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8, fontface = "bold", col = label_colours),
    column_title = column_title, row_title = "LR pair",
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    # Outline every tile, including zero-proportion tiles with white fill.
    rect_gp = grid::gpar(col = "grey30", lwd = 0.6),
    border = TRUE, border_gp = grid::gpar(col = "grey30", lwd = 0.8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      proportion <- mat[i, j]
      # Match the benchmarking report's percentage labels and text contrast.
      label <- if (!is.finite(proportion)) {
        "NA"
      } else if (proportion > 0 && proportion < 0.001) {
        "<0.1%"
      } else {
        sprintf("%.1f%%", 100 * proportion)
      }
      text_colour <- if (is.finite(proportion) && proportion >= 0.72) "white" else "#222222"
      grid::grid.text(label, x, y, gp = grid::gpar(fontsize = 7, col = text_colour))
    },
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"),
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8)
    )
  )
}
