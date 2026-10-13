# Pathology annotation functions extracted from the manuscript analysis.
# Janesick et al. (2023), DOI: 10.1038/s41467-023-43458-x.
# Read the unmodified source workbook; do not save intermediate ROI CSV files.

paperRegionLevels <- function() c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI")

paperBarcode <- function(x) {
  value <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.finite(value)) || any(abs(value - round(value)) > 1e-6)) {
    stop("Paper ROI barcodes must be finite integer identifiers.")
  }
  # Excel writes one barcode in scientific notation; compare numeric IDs.
  as.character(round(value))
}

paperMembershipRegion <- function(x) {
  # Fig. 4b-d uses the opposite DCIS numbering from Fig. 4a and spe$Cluster.
  crosswalk <- c(
    InvasiveA = "Invasive_Tumor",
    DCIS1A = "DCIS_2", DCIS1B = "DCIS_2", DCIS1C = "DCIS_2",
    DCIS1D = "DCIS_2", DCIS1E = "DCIS_2",
    DCIS2A = "DCIS_1", DCIS2B = "DCIS_1", DCIS2C = "DCIS_1"
  )
  if (anyNA(x) || any(!x %in% names(crosswalk))) stop("Unknown paper ROI membership.")
  unname(crosswalk[x])
}

readPaperROIWorkbook <- function(workbook) {
  coordinates <- as.data.frame(readxl::read_excel(workbook,
    sheet = "Fig. 4a ROIs", skip = 2, col_names = FALSE, col_types = "numeric"))
  column_pairs <- paste0(LETTERS[seq(1, 19, 2)], ":", LETTERS[seq(2, 20, 2)])
  vertices <- do.call(rbind, lapply(seq_along(column_pairs), function(i) {
    x <- coordinates[[2L * i - 1L]]
    y <- coordinates[[2L * i]]
    keep <- !is.na(x) | !is.na(y)
    data.frame(column_pair = column_pairs[i], vertex = which(keep),
      x = x[keep], y = y[keep])
  }))
  cells <- as.data.frame(readxl::read_excel(workbook, sheet = "Fig. 4b-d ",
    range = readxl::cell_cols("A:E"),
    col_types = c("text", "text", "text", "numeric", "numeric")))
  list(vertices = vertices, cells = cells)
}

loadPaperROIs <- function(workbook = repo_path("Data", "downloaded",
    "janesick_2023_springer", "janesick_2023_source_data.xlsx")) {
  source <- readPaperROIWorkbook(workbook)
  vertices <- source$vertices
  # O:P duplicates C:D and Q:R duplicates E:F in the source workbook.
  source_columns <- c("S:T", "A:B", "C:D", "E:F", "G:H", "I:J", "K:L", "M:N")
  ids <- c("Invasive", "DCIS_1_A", "DCIS_1_B", "DCIS_1_C",
           "DCIS_2_A", "DCIS_2_B", "DCIS_2_C", "DCIS_2_D")
  regions <- c("Invasive_Tumor", rep("DCIS_1", 3), rep("DCIS_2", 4))
  geometry <- lapply(source_columns, function(column) {
    v <- vertices[vertices$column_pair == column, ]
    v <- v[order(v$vertex), ]
    if (nrow(v) < 3L || any(!is.finite(as.matrix(v[, c("x", "y")])))) {
      stop("Incomplete paper ROI polygon: ", column)
    }
    # Registered against the published cell memberships, not inferred cell types.
    xy <- cbind(v$x, v$y - 3692.59)
    if (!identical(xy[1, ], xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    sf::st_polygon(list(xy))
  })
  components <- sf::st_make_valid(sf::st_sf(
    roi_id = ids, region = regions, source_columns = source_columns,
    geometry = sf::st_sfc(geometry)
  ))
  region_names <- paperRegionLevels()[1:3]
  pooled <- sf::st_sfc(lapply(region_names, function(region) {
    sf::st_union(sf::st_geometry(components)[components$region == region])[[1]]
  }))
  # The source cell-membership table resolves the invasive/DCIS overlap this way.
  pooled[2] <- sf::st_difference(pooled[2], pooled[1])
  pooled[3] <- sf::st_difference(pooled[3], pooled[1])
  if (sum(as.numeric(sf::st_area(sf::st_intersection(pooled[2], pooled[3])))) > 1e-6) {
    stop("DCIS 1 and DCIS 2 ROIs overlap after registration.")
  }
  rois <- sf::st_sf(region = factor(region_names, levels = paperRegionLevels()),
                    geometry = pooled)
  list(rois = rois, components = components, y_offset = -3692.59,
       source_cells = source$cells)
}

classifyPaperRegions <- function(xy, rois) {
  xy <- as.matrix(xy)
  if (ncol(xy) != 2L || any(!is.finite(xy))) stop("Finite X/Y coordinates are required.")
  points <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"))
  hit <- sf::st_intersects(points, rois)
  region <- vapply(hit, function(index) {
    # Invasive is first, including any shared boundary points.
    if (length(index)) as.character(rois$region[index[1]]) else "Outside ROI"
  }, character(1))
  factor(region, levels = paperRegionLevels())
}

validatePaperROIs <- function(paper, spe) {
  source <- paper$source_cells
  source_id <- paperBarcode(source$Barcode)
  current_id <- paperBarcode(spe$Barcode)
  if (anyDuplicated(source_id) || anyDuplicated(current_id)) stop("Duplicate cell barcodes.")
  expected_source <- paperMembershipRegion(source$ROI_Membership)
  predicted_source <- classifyPaperRegions(source[, c("x_centroid", "y_centroid")], paper$rois)
  if (any(as.character(predicted_source) != expected_source)) {
    stop("Registered polygons do not reproduce the published ROI memberships.")
  }

  index <- match(current_id, source_id)
  matched <- !is.na(index)
  if (!any(matched)) stop("No paper ROI cells match this Xenium object.")
  xy <- SpatialExperiment::spatialCoords(spe)
  source_xy <- as.matrix(source[index[matched], c("x_centroid", "y_centroid")])
  coordinate_error <- max(abs(xy[matched, ] - source_xy))
  if (coordinate_error > 1e-5) stop("Paper ROI barcodes and current spatial coordinates disagree.")

  expected <- rep("Outside ROI", length(current_id))
  expected[matched] <- expected_source[index[matched]]
  predicted <- classifyPaperRegions(xy, paper$rois)
  if (any(as.character(predicted) != expected)) {
    stop("ROI polygons and source memberships disagree on current cells.")
  }
  fine_source <- source$celltype[index[matched]]
  fine_current <- as.character(spe$Cluster[matched])
  swap <- c(DCIS_1 = "DCIS_2", DCIS_2 = "DCIS_1")
  is_dcis <- fine_source %in% names(swap)
  if (any(fine_current[is_dcis] != unname(swap[fine_source[is_dcis]]))) {
    stop("The documented DCIS naming crosswalk does not match this dataset.")
  }
  list(
    cells = data.frame(cell_id = colnames(spe), Barcode = current_id,
                       source_ROI = source$ROI_Membership[index], region = predicted),
    audit = data.frame(
      source_roi_cells = nrow(source), current_cells = length(current_id),
      matched_roi_cells = sum(matched),
      source_membership_agreement = mean(as.character(predicted_source) == expected_source),
      current_membership_agreement = mean(as.character(predicted) == expected),
      max_coordinate_error = coordinate_error, polygon_y_offset = paper$y_offset
    )
  )
}

annotatePaperBins <- function(bins, rois) {
  # Region membership is geometric and includes every cell type within an ROI.
  xy <- sf::st_coordinates(sf::st_centroid(bins))
  classifyPaperRegions(xy, rois)
}

splitDCISLabels <- function(cell_type, cluster) {
  broad <- as.character(cell_type)
  fine <- as.character(cluster)
  if (length(broad) != length(fine) || anyNA(broad)) {
    stop("Cell-type and cluster labels must align and cell types cannot be missing.")
  }
  is_dcis <- broad == "DCIS"
  if (any(!fine[is_dcis] %in% c("DCIS_1", "DCIS_2"))) {
    stop("Every DCIS cell must have a DCIS_1 or DCIS_2 Cluster annotation.")
  }
  old_levels <- if (is.factor(cell_type)) levels(cell_type) else unique(broad)
  new_levels <- unlist(lapply(old_levels, function(x) {
    if (x == "DCIS") c("DCIS_1", "DCIS_2") else x
  }), use.names = FALSE)
  broad[is_dcis] <- fine[is_dcis]
  factor(broad, levels = unique(new_levels))
}


# Add cell-level plot metadata without re-binning or running inference.
attach_xenium_cell_metadata <- function(binned, spe, paper = loadPaperROIs()) {
  spe$cell_type_broad <- spe$cell_type
  spe$cell_type <- splitDCISLabels(spe$cell_type_broad, spe$Cluster)
  validation <- validatePaperROIs(paper, spe)
  xy <- SpatialExperiment::spatialCoords(spe)
  binned$cells_metadata <- data.frame(cell_id = colnames(spe),
    x = xy[, 1], y = xy[, 2], cell_type = spe$cell_type,
    cell_type_broad = spe$cell_type_broad, cluster = spe$Cluster,
    region = validation$cells$region)
  binned$paper_rois <- paper$rois
  binned$roi_alignment_audit <- validation$audit
  binned
}

# Export saved or newly computed results. This function never calls blisa().
# Return updated input geometry so callers can save it under Data/generated/.
write_xenium_results <- function(result, binned, output = repo_path("res", "xenium")) {
  stopifnot(identical(result$bins$bin_id, binned$bins$bin_id))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  lr <- result$LR_results
  n_bins <- nrow(result$bins)
  ids <- rownames(lr)
  summary <- lr[, !vapply(lr, is.list, logical(1)), drop = FALSE]
  summary <- data.frame(rank = seq_len(nrow(lr)), lr_id = ids, summary,
    row.names = NULL, check.names = FALSE)
  local_results <- lapply(seq_len(nrow(lr)), function(i) {
    key <- if (lr$ccc_mode[i] == "nearby") "keep_idx_queen" else "keep_idx_dist"
    tested <- seq_len(n_bins) %in% result$spatial_weights[[key]] &
      is.finite(lr$all_pval[[i]]) & !is.na(lr$all_quadrant[[i]])
    native_hh <- seq_len(n_bins) %in% lr$sig_index[[i]]
    stopifnot(all(tested[native_hh]))
    data.frame(lr_id = ids[i], bin_index = seq_len(n_bins),
      bin_id = result$bins$bin_id, p_value = lr$all_pval[[i]],
      score = lr$all_lisa[[i]], quadrant = lr$all_quadrant[[i]],
      tested = tested, native_significant = native_hh)
  })
  summary$n_tested <- vapply(local_results, function(x) sum(x$tested), integer(1))
  summary$n_total_bins <- n_bins
  write_csv(summary, file.path(output, "lri_summary.csv"))
  write_csv(do.call(rbind, local_results), file.path(output, "lri_by_bin.csv"))
  write_csv(result$CCI_scores, file.path(output, "cci_scores.csv"))

  # The following calculation is extracted from blisa 1.0.2 plotCCIspatial().
  # For each native High-High hotspot, sum ligand expression over the hotspot
  # and its neighbours; use receptor expression at that hotspot. Preserve the
  # native first-maximum tie rule and native complex-subunit expression helper.
  get_min_expr <- getFromNamespace("get_min_expr", "blisa")
  ct_names <- names(binned$counts_by_group)
  dominant <- lapply(seq_len(nrow(lr)), function(i) {
    rows <- local_results[[i]][, c("lr_id", "bin_index", "bin_id", "tested",
                                   "native_significant")]
    rows$sender <- rows$receiver <- rows$cell_pair <- NA_character_
    rows$cci_score <- NA_real_
    neighbours <- if (lr$ccc_mode[i] == "nearby")
      result$spatial_weights$queen_nb_full else result$spatial_weights$dist_nb_full
    receptor <- lapply(binned$counts_by_group,
      function(m) get_min_expr(lr$receptor.symbol[i], m))
    ligand <- lapply(binned$counts_by_group,
      function(m) get_min_expr(lr$ligand.symbol[i], m))
    for (h in lr$sig_index[[i]]) {
      neighbourhood <- unique(c(h, neighbours[[h]]))
      r_scores <- vapply(ct_names, function(ct) sum(receptor[[ct]][h]), numeric(1))
      l_scores <- vapply(ct_names, function(ct) sum(ligand[[ct]][neighbourhood]), numeric(1))
      scores <- 0.5 * log2(outer(r_scores, l_scores) + 1)
      best <- which(scores == max(scores), arr.ind = TRUE)[1L, ]
      rows$sender[h] <- ct_names[best[2L]]
      rows$receiver[h] <- ct_names[best[1L]]
      rows$cell_pair[h] <- paste(rows$sender[h], rows$receiver[h], sep = " -> ")
      rows$cci_score[h] <- scores[best[1L], best[2L]]
    }
    rows
  })
  write_csv(do.call(rbind, dominant),
    file.path(output, "xenium_dominant_celltype_pairs_by_bin.csv"))
  # Spatial weights are reusable model-input geometry. Keep them with processed
  # inputs; res/ contains only complete, portable analysis-result CSV tables.
  binned$spatial_weights <- result$spatial_weights
  invisible(binned)
}
