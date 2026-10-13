#!/usr/bin/env Rscript
# Import saved matched-platform inputs/results and the original overlap support.
# Usage from repository root: Rscript codes/analysis/aligned_comparison/04_import_saved_results.R --source-root /path/to/original/workspace
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "aligned_comparison", "alignment_helpers.R"))
suppressPackageStartupMessages(library(SpatialExperiment))
require_packages(c("sf", "Matrix", "digest"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1L] != "--source-root") {
  stop("Supply --source-root followed by the original manuscript workspace.")
}
source_root <- normalizePath(args[2L], mustWork = TRUE)
source_dir <- file.path(source_root, "blisa_on_aligned_xenium_visium")
prepared_dir <- file.path(source_dir, "prepare_alignment", "output", "04_visium_qc")
result_dir <- file.path(source_dir, "blisa_analysis", "output")
v <- readRDS(file.path(prepared_dir, "visium_aligned.rds"))
x <- readRDS(file.path(prepared_dir, "xenium_binned_visium.rds"))
xc <- readRDS(file.path(prepared_dir, "xenium_cells_aligned.rds"))
stopifnot(ncol(v) == 3762L, ncol(xc) == 136932L,
  identical(colnames(v), colnames(x)), identical(spatialCoords(v), spatialCoords(x)))
spot_geometry <- readRDS(file.path(source_dir, "prepare_alignment", "output",
  "audit", "visium_spot_geometry.rds"))
polygons <- spot_geometry[match(colnames(v), spot_geometry$spot_id), ]
stopifnot(identical(as.character(polygons$spot_id), colnames(v)))
metadata(v)$overlap_support_geometry <- sf::st_sf(
  geometry = sf::st_union(sf::st_geometry(polygons)))
metadata(v)$overlap_support_units <- "Visium full-resolution image pixels; x right, y down"
cells <- data.frame(cell_id = colnames(xc), xenium_x_um = xc$xenium_x_um,
  xenium_y_um = xc$xenium_y_um, image_x = xc$image_x, image_y = xc$image_y,
  assigned_spot_id = xc$assigned_spot_id)
stopifnot(nrow(cells) == 136932L, all(cells$assigned_spot_id %in% colnames(v)),
  all(lengths(sf::st_intersects(sf::st_as_sf(cells, coords = c("image_x", "image_y")),
    polygons)) > 0L))
data_dir <- repo_path("Data", "generated", "aligned_overlap")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
metadata(v)$import_source <- "blisa_on_aligned_xenium_visium/prepare_alignment/output/04_visium_qc/visium_aligned.rds"
metadata(v)$import_source_sha256 <- digest::digest(file = file.path(prepared_dir, "visium_aligned.rds"), algo = "sha256")
base::saveRDS(v, file.path(data_dir, "visium_aligned_overlap_qc.rds"))
base::saveRDS(x, file.path(data_dir, "xenium_aggregated_to_visium_overlap.rds"))
write_csv(cells, file.path(data_dir, "xenium_overlap_cells.csv"))
shared_lr <- readRDS(file.path(result_dir, "shared", "lr_input.rds"))
settings <- readRDS(file.path(result_dir, "shared", "settings.rds"))
stopifnot(nrow(shared_lr) == 15L, length(settings$shared_genes) == 307L,
  settings$settings$nsim == 9999L, settings$settings$cpu_threads == 4L)
lr_file <- repo_path("Data", "generated", "ligand_receptor_tables", "aligned_shared_lr_pairs.csv")
write_csv(data.frame(database_id = rownames(shared_lr), shared_lr,
  row.names = NULL, check.names = FALSE), lr_file)

out <- repo_path("res", "aligned_overlap")
for (name in c("visium", "xenium")) {
  arm <- if (name == "visium") "visium" else "xenium_visium"
  result_file <- file.path(result_dir, arm, "result.rds")
  result <- readRDS(result_file)
  inference_input <- readRDS(file.path(result_dir, arm, "input.rds"))
  source_counts <- if (name == "visium") counts(v[settings$shared_genes, ]) else counts(x[settings$shared_genes, ])
  stopifnot(identical(counts(inference_input), source_counts),
    identical(colnames(inference_input), colnames(v)),
    nrow(result$bins) == 3762L, nrow(result$LR_results) == 15L,
    identical(as.integer(result$bins$bin_id), seq_len(ncol(v))))
  if (!is.null(result$bins$spot_id)) stopifnot(identical(result$bins$spot_id, colnames(v)))
  result$bins$spot_id <- colnames(v)
  summary <- summarize_blisa(result, name)
  stopifnot(identical(as.numeric(summary$n_sig), as.numeric(result$LR_results$sig_numbers)))
  summary$rank <- seq_len(nrow(summary))
  for (field in names(settings$settings)) summary[[field]] <- settings$settings[[field]]
  summary$fastlisa_seed <- settings$fastlisa_default_seed
  summary$import_source <- paste0("blisa_on_aligned_xenium_visium/blisa_analysis/output/", arm, "/result.rds")
  summary$import_source_sha256 <- digest::digest(file = result_file, algo = "sha256")
  summary$lr_input_sha256 <- digest::digest(file = lr_file, algo = "sha256")
  write_csv(summary, file.path(out, paste0(name, "_lri_summary.csv")))
  write_csv(blisa_spot_table(result, name), file.path(out, paste0(name, "_lri_by_spot.csv")))
}
message("Imported two saved 15-LR x 3,762-spot arms and 136,932 aligned cells; no inference performed.")
