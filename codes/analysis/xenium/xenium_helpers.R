# Pathology annotation functions extracted from the manuscript analysis.
# Janesick et al. (2023), DOI: 10.1038/s41467-023-43458-x.
# Read the unmodified source workbook; do not save intermediate ROI CSV files.

paperRegionLevels <- function() c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI")

paperBarcode <- function(x) {
  value <- suppressWarnings(as.numeric(as.character(x)))

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

  rois <- sf::st_sf(region = factor(region_names, levels = paperRegionLevels()),
                    geometry = pooled)
  list(rois = rois, components = components, y_offset = -3692.59,
       source_cells = source$cells)
}

classifyPaperRegions <- function(xy, rois) {
  xy <- as.matrix(xy)

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

  expected_source <- paperMembershipRegion(source$ROI_Membership)
  predicted_source <- classifyPaperRegions(source[, c("x_centroid", "y_centroid")], paper$rois)

  index <- match(current_id, source_id)
  matched <- !is.na(index)

  xy <- SpatialExperiment::spatialCoords(spe)
  source_xy <- as.matrix(source[index[matched], c("x_centroid", "y_centroid")])
  coordinate_error <- max(abs(xy[matched, ] - source_xy))

  expected <- rep("Outside ROI", length(current_id))
  expected[matched] <- expected_source[index[matched]]
  predicted <- classifyPaperRegions(xy, paper$rois)

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

  is_dcis <- broad == "DCIS"

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
