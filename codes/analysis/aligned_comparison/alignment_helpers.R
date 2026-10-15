# Pure functions extracted from the manuscript alignment preparation and comparison.
# Sourcing this file does not read data, run inference or write outputs.
# The complete original Visium array defines assignments before overlap/QC selection.

alignment_visium_positions <- function(spe) {
  cd <- SummarizedExperiment::colData(spe)
  xy <- SpatialExperiment::spatialCoords(spe)
  fields <- c("in_tissue", "array_row", "array_col")
  coords <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

  ids <- colnames(spe)

  positions <- data.frame(barcode = ids, as.data.frame(cd[, fields]),
                         xy[, coords, drop = FALSE], row.names = NULL)
  # The example SPE's pixel-named columns have already been converted to
  # micrometres. Respect the recorded unit when recovering image pixels.
  unit <- S4Vectors::metadata(spe)$coord_unit
  if (identical(unit, "micron")) {
    um_per_pixel <- S4Vectors::metadata(spe)$um_per_pixel

    positions[, coords] <- positions[, coords] / um_per_pixel
  }

  positions
}

project_xy <- function(xy, transform) {
  z <- cbind(as.matrix(xy), 1) %*% t(transform)

  z[, 1:2, drop = FALSE] / z[, 3]
}

alignment_assign_cells <- function(xy, positions, diameter_fullres) {
  # Assign against all original spot centres, before overlap/QC subsetting.
  spot_xy <- as.matrix(positions[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  nearest <- RANN::nn2(spot_xy, xy, k = 1)
  centres <- sf::st_as_sf(data.frame(x = spot_xy[, 1], y = spot_xy[, 2]), coords = c("x", "y"))
  # The paper's estimated VISIUM capture area, independent of Xenium image bounds.
  capture <- sf::st_buffer(sf::st_convex_hull(sf::st_union(centres)), diameter_fullres)
  cells <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"))
  inside <- lengths(sf::st_intersects(cells, capture)) > 0
  list(nearest = nearest$nn.idx[, 1],
       distance = nearest$nn.dists[, 1], inside_capture = inside)
}

alignment_sum_counts <- function(counts, nearest, inside_capture, spot_ids) {
  mapping <- Matrix::sparseMatrix(i = which(inside_capture), j = nearest[inside_capture], x = 1,
    dims = c(ncol(counts), length(spot_ids)), dimnames = list(colnames(counts), spot_ids))
  aggregated <- counts %*% mapping

  list(counts = aggregated, mapping = mapping,
       n_cells = as.integer(Matrix::colSums(mapping)))
}

alignment_subset_spe <- function(spe, ids) {
  out <- spe[, ids, drop = FALSE]
  md <- S4Vectors::metadata(out)
  table_names <- c("cell_type_estimates", "cell_type_proportions",
    "cell_type_estimates_original", "cell_type_proportions_original", "cell_type_counts")
  for (name in intersect(table_names, names(md))) {
    value <- md[[name]]

    md[[name]] <- value[ids, , drop = FALSE]
  }
  if (!is.null(md$deconvolution$unmatched_spots))
    md$deconvolution$unmatched_spots <- intersect(md$deconvolution$unmatched_spots, ids)
  S4Vectors::metadata(out) <- md
  out
}

# Four explicit stages; source features (including repeated symbols) are preserved.
build_alignment_stages <- function(x, v, positions, scales, transform) {
  index <- match(colnames(v), positions$barcode)
  v_original <- v
  xxy <- SpatialExperiment::spatialCoords(x)
  vxy <- as.matrix(positions[index, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  aligned_xy <- project_xy(xxy, transform)

  colnames(aligned_xy) <- colnames(vxy)
  rownames(aligned_xy) <- colnames(x)
  rownames(vxy) <- colnames(v)
  x$xenium_x_um <- xxy[, 1]
  x$xenium_y_um <- xxy[, 2]
  x$image_x <- aligned_xy[, 1]
  x$image_y <- aligned_xy[, 2]
  v$source_spatial_x <- SpatialExperiment::spatialCoords(v)[, 1]
  v$source_spatial_y <- SpatialExperiment::spatialCoords(v)[, 2]
  v$image_x <- vxy[, 1]
  v$image_y <- vxy[, 2]
  v$in_tissue <- positions$in_tissue[index]
  v$array_row <- positions$array_row[index]
  v$array_col <- positions$array_col[index]
  v$grid_x_um <- v$array_col * 50
  v$grid_y_um <- v$array_row * 50 * sqrt(3)
  SpatialExperiment::spatialCoords(x) <- aligned_xy
  SpatialExperiment::spatialCoords(v) <- vxy
  if (nrow(SpatialExperiment::imgData(x))) SpatialExperiment::imgData(x) <- S4Vectors::DataFrame()
  S4Vectors::metadata(x)$coord_unit <- "pixel"
  S4Vectors::metadata(x)$um_per_pixel <- NULL
  S4Vectors::metadata(v)$source_coord_unit <- S4Vectors::metadata(v_original)$coord_unit
  S4Vectors::metadata(v)$coord_unit <- "pixel"

  aligned <- list(xenium_cells = x, visium = v)

  assignment <- alignment_assign_cells(aligned_xy, positions, scales$spot_diameter_fullres)
  sums <- alignment_sum_counts(SummarizedExperiment::assay(x, "counts"),
    assignment$nearest, assignment$inside_capture, positions$barcode)
  x$nearest_original_spot_id <- positions$barcode[assignment$nearest]
  x$inside_visium_capture <- assignment$inside_capture
  x$assigned_spot_id <- ifelse(assignment$inside_capture, x$nearest_original_spot_id, NA_character_)
  x$distance_to_spot_px <- assignment$distance
  v$n_xenium_cells <- sums$n_cells[index]
  v$observed_cell_support <- v$n_xenium_cells > 0L
  bcounts <- sums$counts[, colnames(v), drop = FALSE]
  bcd <- S4Vectors::DataFrame(spot_id = colnames(v), in_tissue = v$in_tissue,
    array_row = v$array_row, array_col = v$array_col, image_x = v$image_x, image_y = v$image_y,
    grid_x_um = v$grid_x_um, grid_y_um = v$grid_y_um, n_xenium_cells = v$n_xenium_cells,
    observed_cell_support = v$observed_cell_support,
    total_counts = as.numeric(Matrix::colSums(bcounts)),
    n_genes = as.integer(Matrix::colSums(bcounts > 0)), row.names = colnames(v))
  if ("sample_id" %in% names(SummarizedExperiment::colData(v))) bcd$sample_id <- v$sample_id
  bmeta <- list(coord_unit = "pixel", platform = "Xenium aggregated onto Visium spots")
  if ("cell_type" %in% names(SummarizedExperiment::colData(x))) {
    labels <- as.character(x$cell_type)
    if ("Cluster" %in% names(SummarizedExperiment::colData(x))) {
      dcis <- !is.na(labels) & labels == "DCIS" & x$Cluster %in% c("DCIS_1", "DCIS_2")
      labels[dcis] <- as.character(x$Cluster[dcis])
    }
    labels[is.na(labels) | !nzchar(labels)] <- "Unlabelled"
    labels <- factor(labels)
    indicator <- Matrix::sparseMatrix(i = as.integer(labels), j = seq_along(labels), x = 1,
      dims = c(nlevels(labels), length(labels)), dimnames = list(levels(labels), colnames(x)))
    ct <- as.matrix(Matrix::t(indicator %*% sums$mapping[, colnames(v), drop = FALSE]))
    prop <- ct / v$n_xenium_cells
    prop[v$n_xenium_cells == 0L, ] <- NA_real_
    dominant <- colnames(ct)[max.col(ct, ties.method = "first")]
    dominant[v$n_xenium_cells == 0L] <- NA_character_
    bcd$cell_type <- dominant
    bmeta$cell_type_counts <- ct
    bmeta$cell_type_proportions <- prop
  }
  b <- SpatialExperiment::SpatialExperiment(assays = list(counts = bcounts),
    rowData = SummarizedExperiment::rowData(x), colData = bcd, spatialCoords = vxy, metadata = bmeta)
  binned <- list(xenium_cells = x, visium = v, xenium_binned = b)
  subset_spots <- function(objects, spot_ids, rule) {
    cell_ids <- colnames(objects$xenium_cells)[objects$xenium_cells$assigned_spot_id %in% spot_ids]
    out <- list(xenium_cells = alignment_subset_spe(objects$xenium_cells, cell_ids),
      visium = alignment_subset_spe(objects$visium, spot_ids),
      xenium_binned = alignment_subset_spe(objects$xenium_binned, spot_ids))
    for (name in names(out)) S4Vectors::metadata(out[[name]])$selection <- rule

    out
  }
  overlap <- subset_spots(binned, colnames(v)[v$n_xenium_cells > 0L],
    "Stage 3: retain spots with at least one assigned source Xenium cell and their associated cells; no expression filter.")
  visium_total <- as.numeric(Matrix::colSums(SummarizedExperiment::assay(v, "counts")))
  keep_qc <- Matrix::colSums(SummarizedExperiment::assay(overlap$visium, "counts")) > 3000
  qc <- subset_spots(overlap, colnames(overlap$visium)[keep_qc],
    "Stage 4: within overlap, retain original-feature Visium raw count totals > 3000 and matched Xenium spots/cells.")
  stages <- list(aligned = aligned, binned = binned, overlap = overlap, qc = qc)
  summary <- do.call(rbind, lapply(names(stages), function(name) {
    z <- stages[[name]]
    data.frame(stage = name, xenium_cells = ncol(z$xenium_cells), visium_spots = ncol(z$visium),
      xenium_features = nrow(z$xenium_cells), visium_features = nrow(z$visium))
  }))
  qc_summary <- data.frame(metric = c("All input Visium spots", "All spots passing counts > 3000",
    "Overlap spots before count QC", "Overlap spots removed by count QC", "Overlap spots after count QC",
    "Xenium cells before spatial subsetting", "Xenium cells in overlap", "Xenium cells after Visium QC"),
    value = c(ncol(v), sum(visium_total > 3000), ncol(overlap$visium), sum(!keep_qc), ncol(qc$visium),
      ncol(x), ncol(overlap$xenium_cells), ncol(qc$xenium_cells)))
  cells <- data.frame(cell_id = colnames(x), xenium_x_um = x$xenium_x_um, xenium_y_um = x$xenium_y_um,
    image_x = x$image_x, image_y = x$image_y, assigned_spot_id = x$assigned_spot_id,
    inside_visium_capture = x$inside_visium_capture,
    in_overlap = colnames(x) %in% colnames(overlap$xenium_cells),
    after_visium_qc = colnames(x) %in% colnames(qc$xenium_cells))
  c(stages, list(cell_to_spot = cells,
    summary = summary, qc_summary = qc_summary))
}
