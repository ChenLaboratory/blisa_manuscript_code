# Pure functions extracted from the manuscript alignment preparation and comparison.
# Sourcing this file does not read data, run inference or write outputs.
# The complete original Visium array defines assignments before overlap/QC selection.

alignment_visium_positions <- function(spe) {
  if (!inherits(spe, "SpatialExperiment")) stop("Expected a SpatialExperiment.")
  cd <- SummarizedExperiment::colData(spe)
  xy <- SpatialExperiment::spatialCoords(spe)
  fields <- c("in_tissue", "array_row", "array_col")
  coords <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  if (!all(fields %in% names(cd)) || !all(coords %in% colnames(xy))) {
    stop("Expected the vignette SPE's tissue, array and full-resolution pixel coordinates.")
  }
  ids <- colnames(spe)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Unique, non-empty spot barcodes are required.")
  }
  positions <- data.frame(barcode = ids, as.data.frame(cd[, fields]),
                         xy[, coords, drop = FALSE], row.names = NULL)
  # The example SPE's pixel-named columns have already been converted to
  # micrometres. Respect the recorded unit when recovering image pixels.
  unit <- S4Vectors::metadata(spe)$coord_unit
  if (identical(unit, "micron")) {
    um_per_pixel <- S4Vectors::metadata(spe)$um_per_pixel
    if (length(um_per_pixel) != 1L || !is.numeric(um_per_pixel) ||
        !is.finite(um_per_pixel) || um_per_pixel <= 0) {
      stop("Micrometre coordinates require a positive metadata$um_per_pixel.")
    }
    positions[, coords] <- positions[, coords] / um_per_pixel
  } else if (!identical(unit, "pixel")) {
    stop("Expected SPE metadata$coord_unit to be 'micron' or 'pixel'.")
  }
  if (any(!is.finite(as.matrix(positions[, -1])))) stop("Invalid spot coordinates.")
  positions
}


project_xy <- function(xy, transform) {
  z <- cbind(as.matrix(xy), 1) %*% t(transform)
  if (any(!is.finite(z)) || any(abs(z[, 3]) < 1e-8))
    stop("Invalid homogeneous coordinates in published transform.")
  z[, 1:2, drop = FALSE] / z[, 3]
}

spot_geometry <- function(positions, diameter_fullres) {
  xy <- as.matrix(positions[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  # Fit orientation against the complete array, not a QC-thinned subset.
  lattice <- cbind(1, positions$array_col / 2, positions$array_row * sqrt(3) / 2)
  fit <- qr.solve(lattice, xy)
  residual <- sqrt(rowSums((lattice %*% fit - xy)^2))
  pitch <- sqrt(sum(fit[2, ]^2))
  if (max(residual) > max(3, 0.03 * pitch))
    stop("Visium image positions are not consistent with their full array lattice.")
  # CytAssist image registration can make image metadata inconsistent with the
  # nominal 100/55 µm lattice. Record this; never silently rescale one to fit.
  nominal_ratio_consistent <- abs(pitch / diameter_fullres - 100 / 55) <= 0.1
  # Exact image-space Euclidean Voronoi cells match the nearest-image-centre
  # assignment even under shear. Clip unbounded edge cells to the paper-style
  # capture boundary (convex hull padded by one image-space spot diameter).
  points <- sf::st_sfc(lapply(seq_len(nrow(xy)), function(i) sf::st_point(xy[i, ])))
  capture <- sf::st_buffer(sf::st_convex_hull(sf::st_union(points)), diameter_fullres)
  box <- sf::st_bbox(capture)
  box[c("xmin", "ymin")] <- box[c("xmin", "ymin")] - 2 * pitch
  box[c("xmax", "ymax")] <- box[c("xmax", "ymax")] + 2 * pitch
  polygons <- sf::st_collection_extract(sf::st_voronoi(sf::st_union(points),
    envelope = sf::st_as_sfc(box)), "POLYGON")
  hit <- sf::st_intersects(points, polygons)
  if (any(lengths(hit) != 1)) stop("Cannot assign Voronoi polygons to original Visium centres.")
  polygons <- polygons[vapply(hit, `[`, integer(1), 1)]
  polygons <- sf::st_sf(spot_id = positions$barcode, geometry = polygons)
  polygons <- suppressWarnings(sf::st_intersection(polygons, sf::st_sf(geometry = capture)))
  polygons <- polygons[match(positions$barcode, polygons$spot_id), ]
  if (anyNA(polygons$spot_id) || anyDuplicated(polygons$spot_id)) stop("Invalid original-array Voronoi footprints.")
  list(hex = polygons,
       xy = xy, pitch_fullres = pitch, fit_max_error_px = max(residual),
       nominal_ratio_consistent = nominal_ratio_consistent)
}


# Preparation only; no inference or downloads.

alignment_check_spe <- function(spe, label) {
  if (!inherits(spe, "SpatialExperiment")) stop(label, " must be a SpatialExperiment.")
  for (ids in dimnames(spe)) {
    if (!length(ids) || anyNA(ids) || any(!nzchar(ids)))
      stop(label, " requires nonempty feature and observation identifiers.")
  }
  if (anyDuplicated(colnames(spe))) stop(label, " requires unique observation identifiers.")
  if (!"counts" %in% SummarizedExperiment::assayNames(spe)) stop(label, " lacks raw counts.")
  counts <- SummarizedExperiment::assay(spe, "counts")
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values)) || any(values < 0) || any(values != round(values)))
    stop(label, " counts must be finite, nonnegative integers.")
  invisible(TRUE)
}

alignment_assign_cells <- function(xy, positions, diameter_fullres) {
  required <- c("barcode", "in_tissue", "array_row", "array_col",
                "pxl_col_in_fullres", "pxl_row_in_fullres")
  if (!all(required %in% names(positions)) || anyNA(positions$barcode) ||
      any(!nzchar(positions$barcode)) || anyDuplicated(positions$barcode) ||
      any(!is.finite(as.matrix(positions[, required[-1]]))) ||
      !all(positions$in_tissue %in% c(0, 1))) stop("Invalid original Visium position table.")
  if (length(diameter_fullres) != 1L || !is.finite(diameter_fullres) || diameter_fullres <= 0)
    stop("A positive Visium spot diameter in full-resolution pixels is required.")
  grid <- spot_geometry(positions, diameter_fullres)
  nearest <- RANN::nn2(grid$xy, xy, k = 1)
  centres <- sf::st_as_sf(data.frame(x = grid$xy[, 1], y = grid$xy[, 2]), coords = c("x", "y"))
  # The paper's estimated VISIUM capture area, independent of Xenium image bounds.
  capture <- sf::st_buffer(sf::st_convex_hull(sf::st_union(centres)), diameter_fullres)
  cells <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"))
  inside <- lengths(sf::st_intersects(cells, capture)) > 0
  list(grid = grid, nearest = nearest$nn.idx[, 1],
       distance = nearest$nn.dists[, 1], inside_capture = inside)
}

alignment_sum_counts <- function(counts, nearest, inside_capture, spot_ids) {
  if (length(nearest) != ncol(counts) || length(inside_capture) != ncol(counts) ||
      anyNA(nearest) || anyNA(inside_capture) ||
      any(!nearest %in% seq_along(spot_ids))) stop("Invalid cell-to-spot indices.")
  mapping <- Matrix::sparseMatrix(i = which(inside_capture), j = nearest[inside_capture], x = 1,
    dims = c(ncol(counts), length(spot_ids)), dimnames = list(colnames(counts), spot_ids))
  aggregated <- counts %*% mapping
  expected <- Matrix::rowSums(counts[, inside_capture, drop = FALSE])
  if (!isTRUE(all.equal(as.numeric(Matrix::rowSums(aggregated)), as.numeric(expected), tolerance = 0)))
    stop("Binning did not conserve per-gene raw counts.")
  list(counts = aggregated, mapping = mapping,
       n_cells = as.integer(Matrix::colSums(mapping)))
}

alignment_subset_spe <- function(spe, ids) {
  if (anyNA(ids) || anyDuplicated(ids) || !all(ids %in% colnames(spe)))
    stop("Invalid overlap observation identifiers.")
  out <- spe[, ids, drop = FALSE]
  md <- S4Vectors::metadata(out)
  table_names <- c("cell_type_estimates", "cell_type_proportions",
    "cell_type_estimates_original", "cell_type_proportions_original", "cell_type_counts")
  for (name in intersect(table_names, names(md))) {
    value <- md[[name]]
    if (length(dim(value)) != 2L || anyDuplicated(rownames(value)) ||
        !setequal(rownames(value), colnames(spe)))
      stop("Observation metadata lacks matching row identifiers: ", name)
    md[[name]] <- value[ids, , drop = FALSE]
  }
  if (!is.null(md$deconvolution$unmatched_spots))
    md$deconvolution$unmatched_spots <- intersect(md$deconvolution$unmatched_spots, ids)
  S4Vectors::metadata(out) <- md
  out
}

# Four explicit stages; source features (including repeated symbols) are preserved.
build_alignment_stages <- function(x, v, positions, scales, transform) {
  alignment_check_spe(x, "Xenium")
  alignment_check_spe(v, "Visium")
  stopifnot(setequal(colnames(v), positions$barcode))
  index <- match(colnames(v), positions$barcode)
  x_original <- x
  v_original <- v
  xxy <- SpatialExperiment::spatialCoords(x)
  vxy <- as.matrix(positions[index, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  aligned_xy <- project_xy(xxy, transform)
  stopifnot(max(abs(project_xy(aligned_xy, solve(transform)) - xxy)) < 1e-6)
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
  stopifnot(identical(SummarizedExperiment::assay(x, "counts"), SummarizedExperiment::assay(x_original, "counts")),
    identical(SummarizedExperiment::assay(v, "counts"), SummarizedExperiment::assay(v_original, "counts")),
    identical(dimnames(x), dimnames(x_original)), identical(dimnames(v), dimnames(v_original)))
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
    stopifnot(identical(colnames(out$visium), colnames(out$xenium_binned)),
      identical(SpatialExperiment::spatialCoords(out$visium), SpatialExperiment::spatialCoords(out$xenium_binned)),
      ncol(out$xenium_cells) == sum(out$xenium_binned$n_xenium_cells),
      isTRUE(all.equal(as.numeric(Matrix::rowSums(SummarizedExperiment::assay(out$xenium_cells, "counts"))),
        as.numeric(Matrix::rowSums(SummarizedExperiment::assay(out$xenium_binned, "counts"))), tolerance = 0)))
    out
  }
  overlap <- subset_spots(binned, colnames(v)[v$n_xenium_cells > 0L],
    "Stage 3: retain spots with at least one assigned source Xenium cell and their associated cells; no expression filter.")
  visium_total <- as.numeric(Matrix::colSums(SummarizedExperiment::assay(v, "counts")))
  keep_qc <- Matrix::colSums(SummarizedExperiment::assay(overlap$visium, "counts")) > 3000
  qc <- subset_spots(overlap, colnames(overlap$visium)[keep_qc],
    "Stage 4: within overlap, retain original-feature Visium raw count totals > 3000 and matched Xenium spots/cells.")
  stages <- list(aligned = aligned, binned = binned, overlap = overlap, qc = qc)
  for (stage in stages) for (spe in stage) methods::validObject(spe)
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
  c(stages, list(spot_geometry = assignment$grid$hex, cell_to_spot = cells,
    summary = summary, qc_summary = qc_summary))
}

blisa_spot_table <- function(result, analysis) {
  lr <- result$LR_results
  ids <- as.character(lr$interaction_name)
  spots <- as.character(result$bins$spot_id)
  n <- nrow(result$bins)
  stopifnot(length(ids) == nrow(lr), !anyNA(ids), all(nzchar(ids)),
            !anyDuplicated(ids), length(spots) == n, !anyNA(spots),
            !anyDuplicated(spots))
  rows <- lapply(seq_len(nrow(lr)), function(i) {
    mode <- as.character(lr$ccc_mode[i])
    stopifnot(mode %in% c("nearby", "diffuse"))
    key <- if (mode == "nearby") "keep_idx_queen" else "keep_idx_dist"
    keep <- result$spatial_weights[[key]]
    p <- lr$all_pval[[i]]
    has_lisa <- "all_lisa" %in% names(lr)
    score <- if (has_lisa) lr$all_lisa[[i]] else rep(NA_real_, n)
    quadrant <- as.character(lr$all_quadrant[[i]])
    sig_index <- lr$sig_index[[i]]
    stopifnot(!is.null(keep), !anyNA(keep), all(keep %in% seq_len(n)),
              length(p) == n, length(score) == n, length(quadrant) == n, !anyNA(sig_index),
              !anyDuplicated(sig_index), all(sig_index %in% seq_len(n)))
    # Isolated bins can have finite placeholder P = 1: exclude those explicitly.
    tested <- seq_len(n) %in% keep & is.finite(p)
    if (has_lisa) tested <- tested & is.finite(score)
    significant <- seq_len(n) %in% sig_index
    stopifnot(!any(significant & !tested))
    p[!tested] <- NA_real_
    score[!tested] <- NA_real_
    quadrant[!tested] <- NA_character_
    significant[!tested] <- NA
    data.frame(analysis = analysis, database_id = ids[i],
      lr_label = paste(lr$ligand.symbol[i], lr$receptor.symbol[i], sep = " - "),
      pathway = as.character(lr$pathway_name[i]), mode = mode, spot_id = spots,
      p_value = p, lisa = score, quadrant = quadrant,
      tested = tested, significant = significant,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

summarize_blisa <- function(result, analysis) {
  spots <- blisa_spot_table(result, analysis)
  groups <- split(seq_len(nrow(spots)),
                  factor(spots$database_id, levels = unique(spots$database_id)))
  rows <- lapply(groups, function(index) {
    x <- spots[index, ]
    n_tested <- sum(x$tested)
    n_sig <- sum(x$significant, na.rm = TRUE)
    data.frame(x[1, c("analysis", "database_id", "lr_label", "pathway", "mode")],
      n_sig = n_sig, n_tested = n_tested, n_bins = nrow(x),
      n_untested = nrow(x) - n_tested,
      prop_sig = if (n_tested) n_sig / n_tested else NA_real_,
      prop_sig_all = n_sig / nrow(x), row.names = NULL)
  })
  do.call(rbind, unname(rows))
}
