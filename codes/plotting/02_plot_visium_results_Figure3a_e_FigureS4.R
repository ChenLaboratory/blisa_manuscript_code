# Regenerate Figure 3a-e and Supplementary Figure S4 from the processed SPE and
# complete result CSVs. No inference, preparation, database download or result RDS.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
source(repo_path("codes", "helpers", "plotting_helpers.R"))
source(repo_path("codes", "helpers", "manuscript_plotting_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "S4Vectors",
  "ggplot2", "scales", "gridExtra", "data.table", "R.utils", "ComplexHeatmap", "circlize",
  "scider", "sf", "ggnewscale"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(SpatialExperiment)
  library(blisa)
})

# Shared image/spot functions retain the original report's calibration, palettes,
# display radius, ROI registration and distinction between untested and nonsignificant.


spatial_spot_radius <- function(spots, coverage = 0.88) {
  xy <- unique(as.matrix(spots[, c("plot_x", "plot_y")]))
  if (nrow(xy) < 2L || any(!is.finite(xy))) {
    stop("Spot sizing requires at least two distinct finite coordinates.")
  }
  if (length(coverage) != 1L || !is.finite(coverage) || coverage <= 0 || coverage >= 1) {
    stop("Spot coverage must lie strictly between zero and one.")
  }
  index <- unique(round(seq(1, nrow(xy), length.out = min(nrow(xy), 256L))))
  nearest <- vapply(index, function(i) {
    distance2 <- (xy[, 1L] - xy[i, 1L])^2 + (xy[, 2L] - xy[i, 2L])^2
    sqrt(min(distance2[distance2 > 0]))
  }, numeric(1))
  stats::median(nearest) * coverage / 2
}


geom_spatial_spot <- function(mapping = NULL, data = NULL, ..., radius,
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  if (length(radius) != 1L || !is.finite(radius) || radius <= 0) {
    stop("Supply a positive finite display radius for spatial spots.")
  }
  geom <- ggplot2::ggproto(NULL, ggplot2::Geom,
    required_aes = c("x", "y"),
    default_aes = ggplot2::aes(fill = "grey80", colour = "white", alpha = 1,
                               linewidth = 0.05, linetype = 1),
    draw_key = function(data, params, size) {
      data$shape <- 21; data$size <- 2.2; data$stroke <- data$linewidth
      ggplot2::draw_key_point(data, params, size)
    },
    draw_panel = function(data, panel_params, coord, radius, na.rm = FALSE) {
      if (!nrow(data)) return(grid::nullGrob())
      theta <- seq(0, 2 * pi, length.out = 25L)[-25L]
      polygon <- data[rep(seq_len(nrow(data)), each = length(theta)), , drop = FALSE]
      polygon$x <- polygon$x + radius * rep(cos(theta), times = nrow(data))
      polygon$y <- polygon$y + radius * rep(sin(theta), times = nrow(data))
      polygon$group <- rep(seq_len(nrow(data)), each = length(theta))
      ggplot2::GeomPolygon$draw_panel(polygon, panel_params, coord)
    })
  ggplot2::layer(geom = geom, stat = "identity", position = "identity",
    mapping = mapping, data = data, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(radius = radius, na.rm = na.rm, ...))
}


spatial_canvas <- function(spots, background = NULL, borders = NULL) {
  x <- spots
  if (!is.null(background)) {
    x$plot_x <- x$image_x * background$scale
    x$plot_y <- -x$image_y * background$scale
    plot <- ggplot2::ggplot() + ggplot2::annotation_raster(background$image,
      xmin = 0, xmax = background$width, ymin = -background$height, ymax = 0)
  } else {
    x$plot_x <- x$grid_x_um; x$plot_y <- -x$grid_y_um
    plot <- ggplot2::ggplot()
  }
  radius <- spatial_spot_radius(x)
  # Include complete outer disks while keeping the tissue panels compact.
  x_limits <- range(x$plot_x) + c(-radius, radius)
  y_limits <- range(x$plot_y) + c(-radius, radius)
  if (!is.null(borders)) {
    if (is.null(background)) stop("Registered image borders require the calibrated image background.")
    # All panels in a triptych include the same extent, even though only the
    # hotspot panel draws the outlines. Some ROIs extend beyond retained spots.
    x_limits <- range(x_limits, borders$x * background$scale)
    y_limits <- range(y_limits, -borders$y * background$scale)
  }
  plot <- plot + ggplot2::expand_limits(
    x = x_limits, y = y_limits)
  list(data = x, spot_radius = radius,
    plot = plot + ggplot2::coord_fixed() + ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 8, face = "bold"),
        legend.text = ggplot2::element_text(size = 7, colour = "black"),
        legend.key.height = grid::unit(2.5, "mm"),
        legend.key.width = grid::unit(8, "mm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(1, "mm"),
        plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
        plot.margin = ggplot2::margin(2, 2, 2, 2)))
}


plot_spot_values <- function(spots, values, title, background = NULL, discrete = FALSE,
                             point_size = NULL, point_stroke = 0,
                             legend_title = NULL, borders = NULL) {
  canvas <- spatial_canvas(spots, background, borders)
  canvas$data$value <- values
  point_colour <- if (point_stroke > 0) "white" else "transparent"
  layer <- if (is.null(point_size)) {
    geom_spatial_spot(data = canvas$data, ggplot2::aes(plot_x, plot_y, fill = value),
      radius = canvas$spot_radius, colour = point_colour, linewidth = point_stroke)
  } else {
    ggplot2::geom_point(data = canvas$data, ggplot2::aes(plot_x, plot_y, fill = value),
      shape = 21, colour = point_colour, size = point_size, stroke = point_stroke)
  }
  plot <- canvas$plot + layer + ggplot2::labs(title = title, fill = legend_title)
  if (discrete) plot + ggplot2::scale_fill_viridis_d(na.value = "grey80") else
    # Exact continuous defaults from scider::plotSpatial/plotGrid, also used
    # by Xenium's plotBinned: rev(scider:::col.spec), low purple to high red.
    # Copy the stops to avoid depending on an unexported package object.
    plot + ggplot2::scale_fill_gradientn(
      colours = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
                  "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F"),
      na.value = "grey80", breaks = scales::breaks_pretty(n = 3),
      guide = ggplot2::guide_colourbar(
        title.position = "top", barwidth = grid::unit(22, "mm"),
        barheight = grid::unit(2.5, "mm"), ticks.colour = "white"))
}


add_pathology_borders <- function(plot, borders, x_scale = 1, y_scale = 1,
                                  show_legend = FALSE, linewidth = 0.6) {
  if (is.null(borders)) return(plot)
  colours <- c(Invasive_Tumor = "#377EB8", DCIS_1 = "#4DAF4A", DCIS_2 = "#006D2C")
  if (!is.data.frame(borders) || !all(c("x", "y", "piece", "region") %in% names(borders)) ||
      !nrow(borders) || !is.numeric(borders$x) || !is.numeric(borders$y) ||
      any(!is.finite(borders$x)) || any(!is.finite(borders$y)) ||
      anyNA(borders$piece) || any(!nzchar(as.character(borders$piece))) ||
      any(!borders$region %in% names(colours))) {
    stop("Supply registered pathology borders with finite x/y, piece and region columns.")
  }
  if (length(x_scale) != 1L || length(y_scale) != 1L ||
      !is.finite(x_scale) || !is.finite(y_scale) || x_scale == 0 || y_scale == 0) {
    stop("Pathology borders require finite nonzero coordinate scales.")
  }
  borders$plot_x <- borders$x * x_scale
  borders$plot_y <- borders$y * y_scale
  borders$path_group <- interaction(borders$region, borders$piece, drop = TRUE)
  # Native pathway plots use colour for p-values; custom LR maps use fill.
  if (!is.null(plot$scales$get_scales("colour"))) {
    require_packages("ggnewscale")
    plot <- plot + ggnewscale::new_scale_colour()
  }
  plot +
    ggplot2::geom_path(data = borders,
      ggplot2::aes(plot_x, plot_y, group = path_group), inherit.aes = FALSE,
      colour = "white", linewidth = linewidth + 0.3, show.legend = FALSE) +
    ggplot2::geom_path(data = borders,
      ggplot2::aes(plot_x, plot_y, group = path_group, colour = region),
      inherit.aes = FALSE, linewidth = linewidth, show.legend = show_legend) +
    ggplot2::scale_colour_manual(values = colours, breaks = names(colours),
      limits = names(colours), drop = FALSE, name = "Pathology ROI",
      labels = c("Invasive tumour", "DCIS 1", "DCIS 2"),
      guide = ggplot2::guide_legend(order = 2, nrow = 1))
}


plot_hotspot_values <- function(spots, rows, title, background = NULL,
                                point_size = NULL, point_stroke = 0, borders = NULL) {
  if (anyDuplicated(rows$spot_id) || !setequal(spots$spot_id, rows$spot_id)) {
    stop("Hotspot maps require exactly one row per retained spot.")
  }
  index <- match(spots$spot_id, rows$spot_id)
  if (anyNA(index)) stop("Hotspot map requires an explicit row for every retained spot.")
  rows <- rows[index, , drop = FALSE]
  canvas <- spatial_canvas(spots, background, borders)
  d <- canvas$data
  d$tested <- rows$tested %in% TRUE; d$significant <- rows$native_significant %in% TRUE
  if (any(d$significant & (!d$tested | !is.finite(rows$p_value) | rows$p_value < 0 | rows$p_value > 1))) {
    stop("Hotspot maps require valid p-values for tested significant spots.")
  }
  d$value <- -log10(pmax(rows$p_value, .Machine$double.xmin))
  # Every retained spot has a background state; significant tested spots keep
  # their native p-values and use the warm colour scale.
  background_fill <- ifelse(d$tested, "grey80", "white")
  point_colour <- if (point_stroke > 0) "white" else "transparent"
  layer <- function(data, mapping, ...) {
    if (is.null(point_size)) {
      geom_spatial_spot(data = data, mapping = mapping, radius = canvas$spot_radius,
        colour = point_colour, linewidth = point_stroke, ...)
    } else {
      ggplot2::geom_point(data = data, mapping = mapping, shape = 21, colour = point_colour,
        size = point_size, stroke = point_stroke, ...)
    }
  }
  plot <- canvas$plot +
    layer(d, ggplot2::aes(plot_x, plot_y), fill = background_fill) +
    layer(d[d$significant, ], ggplot2::aes(plot_x, plot_y, fill = value)) +
    ggplot2::scale_fill_gradientn(colours = c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F"),
      breaks = scales::breaks_pretty(n = 3),
      guide = ggplot2::guide_colourbar(title.position = "top",
        barwidth = grid::unit(22, "mm"), barheight = grid::unit(2.5, "mm"),
        ticks.colour = "white")) +
    ggplot2::labs(title = title, fill = expression(-log[10](italic(p))))
  if (!is.null(borders)) {
    if (is.null(background)) stop("Registered image borders require the calibrated image background.")
    plot <- add_pathology_borders(plot, borders,
      x_scale = background$scale, y_scale = -background$scale)
  }
  plot
}


lr_triptych <- function(spe, manifest_row, lr_spot, background = NULL, point_size = NULL,
                        borders = NULL) {
  counts <- SummarizedExperiment::assay(spe, "counts")
  spots <- as.data.frame(SummarizedExperiment::colData(spe))
  id <- manifest_row$lr_id
  rows <- lr_spot[lr_spot$lr_id == id, , drop = FALSE]
  l <- minimum_partner_counts(counts, manifest_row$ligand.symbol)
  r <- minimum_partner_counts(counts, manifest_row$receptor.symbol)
  # Preserve report map styles; stack expression maps as in the manuscript.
  ligand <- plot_spot_values(spots, log1p(l), paste0("Ligand: ", manifest_row$ligand.symbol),
    background, point_size = point_size, legend_title = "log(1 + counts)", borders = borders)
  receptor <- plot_spot_values(spots, log1p(r), paste0("Receptor: ", manifest_row$receptor.symbol),
    background, point_size = point_size, legend_title = "log(1 + counts)", borders = borders)
  hotspots <- plot_hotspot_values(spots, rows, "PPIA-BSG hotspots", background,
    point_size = point_size, borders = borders)
  gridExtra::arrangeGrob(
    gridExtra::arrangeGrob(ligand, receptor, ncol = 1), hotspots,
    ncol = 2, widths = c(2.6, 4.6))
}


plot_visium_cell_types_with_rois <- function(spe, borders) {
  require_packages(c("SpatialExperiment", "SummarizedExperiment", "ggplot2"))
  cd <- SummarizedExperiment::colData(spe)
  if (!all(c("image_x", "image_y") %in% names(cd))) {
    stop("Prepare the export's image_x/image_y full-resolution pixel coordinates before plotting.")
  }
  type_column <- intersect(c("dominant_ct", "cell_type"), names(cd))[1]
  if (is.na(type_column)) stop("Add the dominant cell-type annotations before plotting.")
  cell_types <- as.character(cd[[type_column]])
  cell_types[is.na(cell_types) | !nzchar(cell_types)] <- "Unknown"
  spots <- data.frame(image_x = cd$image_x, image_y = cd$image_y,
                      cell_type = cell_types)
  if (!nrow(spots) || any(!is.finite(as.matrix(spots[, c("image_x", "image_y")])))) {
    stop("The cell-type map requires non-empty, finite spot coordinates.")
  }
  regions <- c("Invasive_Tumor", "DCIS_1", "DCIS_2")
  if (!is.data.frame(borders) || !all(c("x", "y", "piece", "region") %in% names(borders)) ||
      !nrow(borders) || !is.numeric(borders$x) || !is.numeric(borders$y) ||
      any(!is.finite(borders$x)) || any(!is.finite(borders$y)) ||
      anyNA(borders$piece) || any(!nzchar(as.character(borders$piece))) ||
      any(!borders$region %in% regions)) {
    stop("Supply registered pathology borders with finite x/y, piece and region columns.")
  }
  borders$path_group <- interaction(borders$region, borders$piece, drop = TRUE)

  # Use the same category order and colours as the Xenium report, including
  # categories with nonzero proportions but no dominant Visium spots.
  colours <- visium_cell_type_palette()
  unexpected <- setdiff(unique(cell_types), c(names(colours), "Unknown"))
  if (length(unexpected)) {
    stop("Harmonize Visium cell types with annotate_visium() before plotting: ",
         paste(unexpected, collapse = ", "))
  }
  if ("Unknown" %in% cell_types) colours <- c(colours, Unknown = "grey75")
  spots$cell_type <- factor(spots$cell_type, levels = names(colours))
  sample_ids <- unique(as.character(spe$sample_id))
  if (length(sample_ids) != 1L || is.na(sample_ids) || !nzchar(sample_ids)) {
    stop("The cell-type map requires one Visium sample.")
  }
  images <- SpatialExperiment::imgData(spe)
  candidates <- which(images$sample_id == sample_ids &
                        is.finite(images$scaleFactor) & images$scaleFactor > 0)
  if (!length(candidates)) stop("The cell-type map requires a calibrated embedded image.")
  image_index <- candidates[which.max(images$scaleFactor[candidates])]
  raster <- SpatialExperiment::imgRaster(spe, sample_id = sample_ids,
                                         image_id = images$image_id[image_index])
  background <- list(image = raster, width = ncol(raster), height = nrow(raster),
                     scale = images$scaleFactor[image_index])
  canvas <- spatial_canvas(spots, background)
  # Spots and borders receive the identical pixel-to-image transform. Negative
  # image y places row zero at the top without a second reversed scale.
  borders$plot_x <- borders$x * background$scale
  borders$plot_y <- -borders$y * background$scale

  canvas$plot +
    ggplot2::expand_limits(x = c(0, background$width), y = c(-background$height, 0)) +
    ggplot2::geom_point(data = canvas$data, ggplot2::aes(plot_x, plot_y, fill = cell_type),
                        shape = 21, colour = "transparent", stroke = 0, size = 1.35,
                        show.legend = c(fill = TRUE, colour = FALSE)) +
    ggplot2::geom_path(data = borders, ggplot2::aes(plot_x, plot_y, group = path_group),
                       inherit.aes = FALSE, colour = "white", linewidth = 1.8,
                       show.legend = FALSE) +
    ggplot2::geom_path(data = borders,
                       ggplot2::aes(plot_x, plot_y, group = path_group, colour = region),
                       inherit.aes = FALSE, linewidth = 1.2) +
    # Apply opacity only to the fill: layer-wide alpha would turn the
    # transparent point edge into a visible white stroke.
    ggplot2::scale_fill_manual(values = scales::alpha(colours, 0.9),
                               labels = function(x) gsub("_", " ", x),
                               limits = names(colours), drop = FALSE,
                               name = "Dominant inferred cell type") +
    ggplot2::scale_colour_manual(values = colours[regions], breaks = regions,
                                 labels = c("Invasive tumour", "DCIS 1", "DCIS 2"),
                                 name = "Pathology ROI") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1,
        override.aes = list(shape = 21, size = 3, alpha = 1, colour = NA, stroke = 0)),
      colour = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(title = "Visium after QC: dominant cell types",
                  subtitle = "Outlines: paper Figure 4 pathology ROIs registered to Visium") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "right", legend.title = ggplot2::element_text(face = "bold"),
                   legend.key.height = grid::unit(0.35, "cm"),
                   plot.title = ggplot2::element_text(face = "bold"))
}


regional_proportion_overview <- function(mat, colours, title, show_lr_names = FALSE,
                                         show_legend = TRUE) {
  ComplexHeatmap::Heatmap(
    mat, name = "Hotspot\nproportion",
    col = circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
      c("white", "#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")),
    na_col = "grey80", cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = show_lr_names, show_column_names = TRUE, use_raster = FALSE,
    show_heatmap_legend = show_legend,
    row_names_side = "left", row_names_gp = grid::gpar(fontsize = 5.3),
    width = grid::unit(ncol(mat) * 10, "mm"),
    height = grid::unit(196, "mm"),
    column_labels = cellLabel(colnames(mat)), column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8, fontface = "bold",
      col = unname(colours[colnames(mat)])),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    rect_gp = grid::gpar(col = if (nrow(mat) <= 100L) "grey70" else NA, lwd = 0.25),
    border = TRUE, border_gp = grid::gpar(col = "grey30", lwd = 0.6),
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9)))
}


minimum_partner_counts <- function(counts, symbols) {
  genes <- unique(trimws(strsplit(symbols, ",", fixed = TRUE)[[1]]))
  if (!length(genes) || any(!genes %in% rownames(counts))) stop("Incomplete partner: ", symbols)
  Reduce(pmin, lapply(genes, function(g) as.numeric(counts[g, ])))
}

# The original Visium R Markdown uses the standard PDF device (not Cairo).
# Preserve that device, default point size, and the original chunk dimensions.
save_panel <- function(plot, panel, width, height) {
  save_figure_panel(plot, panel, width, height, device = grDevices::pdf)
}


spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))
stopifnot(nrow(spe) == 18085L, ncol(spe) == 4412L)
summary <- utils::read.csv(repo_path("res", "visium", "lri_summary.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
spot_files <- sort(list.files(repo_path("res", "visium"),
  pattern = "^lri_by_spot_part[0-9]{3}\\.csv\\.gz$", full.names = TRUE))
stopifnot(length(spot_files) == 3L)
local <- as.data.frame(data.table::rbindlist(lapply(spot_files, function(path) {
  data.table::fread(path,
    select = c("lr_id", "spot_id", "p_value", "tested", "native_significant"),
    na.strings = "NA")
})))
summary <- summary[order(summary$rank), , drop = FALSE]
stopifnot(nrow(summary) == 1616L, !anyDuplicated(summary$lr_id),
          nrow(local) == nrow(summary) * ncol(spe),
          identical(local$lr_id, rep(summary$lr_id, each = ncol(spe))),
          identical(local$spot_id, rep(colnames(spe), times = nrow(summary))),
          !anyNA(local$tested), !anyNA(local$native_significant),
          all(!local$native_significant | local$tested))
spots <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
borders <- S4Vectors::metadata(spe)$pathology$borders_visium_pixels
stopifnot(is.data.frame(borders), nrow(borders) > 0L)
images <- SpatialExperiment::imgData(spe)
visium <- S4Vectors::metadata(spe)$visium
image_index <- which(images$sample_id == visium$sample_id & images$image_id == visium$image_id)
stopifnot(length(image_index) == 1L)
raster <- SpatialExperiment::imgRaster(spe, sample_id = visium$sample_id,
                                      image_id = visium$image_id)
background <- list(image = raster, width = ncol(raster), height = nrow(raster),
                    scale = images$scaleFactor[image_index])

# Rebuild the native plotting object from complete CSV results. This reconstructs
# significant lists only; it does not recompute weights, LRI or permutations.
significant <- matrix(local$native_significant, nrow = ncol(spe))
p_values <- matrix(local$p_value, nrow = ncol(spe))
stopifnot(identical(as.integer(colSums(significant)), as.integer(summary$sig_numbers)))
lr_native <- summary
rownames(lr_native) <- summary$lr_id
lr_native$sig_index <- I(lapply(seq_len(nrow(summary)), function(i) which(significant[, i])))
lr_native$sig_pval <- I(lapply(seq_len(nrow(summary)), function(i)
  p_values[significant[, i], i]))
# With an image background, native plotHotspots draws img_x/img_y and uses the
# sf geometry only for an invisible legend anchor. Retain the original SPE's
# micrometre coordinates rather than substituting the inference lattice.
image_xy <- SpatialExperiment::spatialCoords(spe)[,
  c("pxl_col_in_fullres", "pxl_row_in_fullres"), drop = FALSE]
stopifnot(identical(S4Vectors::metadata(spe)$coord_unit, "micron"))
bins <- sf::st_as_sf(data.frame(spot_id = colnames(spe),
  img_x = image_xy[, 1], img_y = image_xy[, 2],
  x = spe$grid_x_um, y = spe$grid_y_um), coords = c("x", "y"))
res <- structure(list(LR_results = lr_native, bins = bins,
  spatial_weights = NULL, CCI_scores = NULL, level = "lr"), class = c("blisa", "list"))
res_pathway <- blisa::blisaPathway(res)
stopifnot(inherits(res_pathway, "blisa"), identical(res_pathway$level, "pathway"),
          identical(rownames(res_pathway$LR_results)[1L], "LAMININ"))
# Keep the source chunk names below; these aliases refer to the same annotated
# native-feature SPE and already registered borders loaded above.
spe_qc <- spe
visium_borders <- borders

# Figure 3a: exact cell-types-pathology chunk (9 x 6.8 inches).
save_panel(plot_visium_cell_types_with_rois(spe, borders), "figure3a", 9, 6.8)

# Figure 3b: exact lr-rank chunk, preserving the native annotation colours,
# category ordering, margins, legend settings and PDF dimensions.
p_lr_rank <- plotLRrank(res, top = 30, flip = TRUE, pt_size = 2.2) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_pretty(n = 5),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "High–High hotspots", title = "Top 30 LR pairs", colour = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    text = element_text(family = "sans"),
    axis.text = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(size = 7.5, angle = 65, hjust = 1, vjust = 1,
                               margin = margin(t = 4)),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.height = grid::unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(3, 4, 2, 3)
  )
save_panel(p_lr_rank, "figure3b", 7.2, 5.4)

# Figure 3c: original PPIA_BSG maps, arranged as in the updated manuscript.
pair <- summary[summary$lr_id == "PPIA_BSG", , drop = FALSE]
stopifnot(nrow(pair) == 1L)
save_panel(lr_triptych(spe, pair, local, background, borders = borders),
           "figure3c", 7.2, 6.2)

# Figure 3d: exact pathway-rank chunk, including fixed 500-hotspot tick spacing.
p_pathway_rank <- plotLRrank(res_pathway, flip = TRUE) +
  scale_size_continuous(name = "LR pairs", range = c(1.5, 4.2),
                        breaks = scales::breaks_pretty(n = 3)) +
  guides(size = guide_legend(nrow = 1)) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 2000, by = 500),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "High–High hotspots", title = "Top 30 pathways") +
  theme(
    text = element_text(family = "sans"),
    axis.text = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(size = 8, angle = 65, hjust = 1, vjust = 1,
                               margin = margin(t = 4)),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8, face = "bold", margin = margin(r = 4)),
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7, face = "bold"),
    legend.key.height = grid::unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(3, 4, 2, 3)
  )
save_panel(p_pathway_rank, "figure3d", 8.2, 2.8)

# Figure 3e: exact pathway-hotspots chunk. Native pathway maps use their own
# colour range, title, grey background states and separate p-value/ROI legends.
# Match the image background to the native result's micrometre coordinates.
spe_image <- spe
imgData(spe_image)$scaleFactor <- imgData(spe_image)$scaleFactor / metadata(spe_image)$um_per_pixel
p_img <- scider::plotImage(spe_image, image_id = metadata(spe_qc)$visium$image_id)
p_pathway <- plotHotspots(res_pathway, index = 1, background = p_img)
# Retain BLISA's spot colours and p-value scale; size disks to the spot grid.
for (i in seq_along(p_pathway$layers)) {
  layer <- p_pathway$layers[[i]]
  if (is.data.frame(layer$data) && "fill_col" %in% names(layer$data)) {
    radius <- spatial_spot_radius(data.frame(plot_x = layer$data$x, plot_y = layer$data$y))
    p_pathway$layers[[i]] <- geom_spatial_spot(
      mapping = layer$mapping, data = layer$data, radius = radius,
      colour = "grey95", linewidth = 0.05
    )
  }
}
p_pathway <- p_pathway +
  labs(subtitle = NULL) + theme_void(base_size = 9) +
  spatialTheme(font_size = 9) +
  guides(colour = guide_colourbar(
    title = expression(-log[10](p)), title.position = "top", title.hjust = 0.5,
    available_aes = "any", order = 1,
    barwidth = grid::unit(26, "mm"), barheight = grid::unit(2.5, "mm")
  ))
# Native hotspot coordinates are micrometres, with y reversed by the plot scale.
# Include full ROI outlines beyond the native spot crop, preserving its orientation.
pathway_pixel_size <- metadata(spe)$um_per_pixel
p_pathway <- p_pathway + coord_fixed(
  xlim = range(p_pathway$coordinates$limits$x, visium_borders$x * pathway_pixel_size),
  ylim = rev(range(p_pathway$coordinates$limits$y, visium_borders$y * pathway_pixel_size))
)
p_pathway <- add_pathology_borders(p_pathway, visium_borders,
  x_scale = metadata(spe)$um_per_pixel, y_scale = metadata(spe)$um_per_pixel,
  show_legend = TRUE, linewidth = 0.8) + theme(legend.box = "vertical")
save_panel(p_pathway, "figure3e", 4.6, 4.8)

# S4 retains the native LR order and divides by ALL retained spots in each region,
# including untested spots. Non-tumour means the complement of the three ROIs.
regions <- as.character(spe$region)
regions[regions == "Outside ROI"] <- "Non-tumour"
region_order <- c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Non-tumour")
stopifnot(!anyNA(regions), all(regions %in% region_order))
proportions <- vapply(region_order, function(region) {
  keep <- regions == region
  if (!any(keep)) return(rep(NA_real_, nrow(summary)))
  colSums(significant[keep, , drop = FALSE]) / sum(keep)
}, numeric(nrow(summary)))
rownames(proportions) <- summary$lr_id
region_colours <- c(Invasive_Tumor = "#377EB8", DCIS_1 = "#4DAF4A",
  DCIS_2 = "#006D2C", "Non-tumour" = "#555555")
# Export the two panels separately at the original report dimensions.
top100_grob <- heatmapGrob(regional_proportion_overview(
  proportions[seq_len(100L), , drop = FALSE], region_colours,
  "Top 100 LR pairs", show_lr_names = TRUE, show_legend = FALSE),
  width = 4, height = 8.8)
all_grob <- heatmapGrob(regional_proportion_overview(proportions,
  region_colours, "All 1,616 LR pairs"), width = 3.2, height = 8.8)
save_figure_panel(top100_grob, "figureS4a", 4, 8.8)
save_figure_panel(all_grob, "figureS4b", 3.2, 8.8)
