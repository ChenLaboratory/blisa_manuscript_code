#!/usr/bin/env Rscript
# Run from the repository root after making the two published SPEs available.
# No download or BLISA inference occurs here. Only final overlap data are saved.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "aligned_comparison", "alignment_helpers.R"))
suppressPackageStartupMessages(library(SpatialExperiment))

xenium <- read_source_data("xenium")
visium <- read_source_data("visium")
stopifnot(identical(dim(xenium), c(313L, 163797L)),
          identical(dim(visium), c(18085L, 4992L)))

# Preserve the original 4,992-spot array for nearest-centre assignment. Do not
# substitute the standalone post-QC Visium object or refit after subsetting.
positions <- alignment_visium_positions(visium)
parameter_dir <- repo_path("Data", "parameters", "from_10x_genomics")
matrix_file <- file.path(parameter_dir, "xenium_to_visium_transform.csv")
calibration_file <- file.path(parameter_dir, "xenium_coordinate_calibration.json")
homography <- as.matrix(read.csv(matrix_file, header = FALSE, check.names = FALSE))
storage.mode(homography) <- "double"
calibration <- jsonlite::fromJSON(calibration_file)
stopifnot(identical(dim(homography), c(3L, 3L)), all(is.finite(homography)),
          abs(det(homography)) > 1e-12, calibration$pixel_size_um == 0.2125)
transform <- homography %*% diag(c(1 / calibration$pixel_size_um,
                                  1 / calibration$pixel_size_um, 1))
scales <- list(spot_diameter_fullres = 55 / metadata(visium)$um_per_pixel)

# Transform cell centroids, sum whole-cell counts into nearest original spots,
# select cell-supported overlap, then retain original-feature Visium counts >3,000.
# Earlier stages are held in memory to check counts; they are not written out.
prepared <- build_alignment_stages(xenium, visium, positions, scales, transform)
v <- prepared$qc$visium
x <- prepared$qc$xenium_binned
cells <- prepared$cell_to_spot[prepared$cell_to_spot$after_visium_qc,
  c("cell_id", "xenium_x_um", "xenium_y_um", "image_x", "image_y", "assigned_spot_id")]
stopifnot(ncol(v) == 3762L, nrow(cells) == 136932L,
          identical(colnames(v), colnames(x)),
          identical(spatialCoords(v), spatialCoords(x)))

# This is the exact union of retained polygons computed against ALL original
# spots. Store it once for the Figure 3f outline; never rebuild it from QC centres.
polygons <- prepared$spot_geometry[
  match(colnames(v), prepared$spot_geometry$spot_id), ]
stopifnot(identical(as.character(polygons$spot_id), colnames(v)))
metadata(v)$overlap_support_geometry <- sf::st_sf(
  geometry = sf::st_union(sf::st_geometry(polygons)))
metadata(v)$overlap_support_units <- "Visium full-resolution image pixels; x right, y down"
provenance <- list(
  source = "ChenLaboratory/example_data release v1.0.0",
  matrix_sha256 = digest::digest(file = matrix_file, algo = "sha256"),
  calibration_sha256 = digest::digest(file = calibration_file, algo = "sha256"),
  xenium_um_to_visium_pixels = transform,
  assignment = "Nearest original Visium centre inside convex hull buffered by one spot diameter",
  overlap = "At least one assigned Xenium cell, followed by original-feature Visium counts >3000",
  counts = "Whole-cell Xenium counts summed by centroid; all source feature rows retained",
  stage_counts = prepared$summary)
metadata(v)$alignment <- metadata(x)$alignment <- provenance

out <- repo_path("Data", "generated", "aligned_overlap")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base::saveRDS(v, file.path(out, "visium_aligned_overlap_qc.rds"))
base::saveRDS(x, file.path(out, "xenium_aggregated_to_visium_overlap.rds"))
write_csv(cells, file.path(out, "xenium_overlap_cells.csv"))
