# Run from the repository root after obtaining the release RDS and references.
# This stage preserves native feature rows, the embedded image and coordinate calibration.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
source(repo_path("codes", "analysis", "xenium", "xenium_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "S4Vectors",
                   "Matrix", "readxl", "digest", "jsonlite", "sf"))

# 1. Retain original feature rows and filter spots with raw total counts >3000.
spe <- read_source_data("visium")

original_positions <- visium_position_table(spe)
database <- load_visium_cellchat()
before <- visium_qc_dimensions(spe, database)
total_counts <- as.numeric(Matrix::colSums(SummarizedExperiment::assay(spe, "counts")))
qc <- data.frame(spot_id = colnames(spe), total_counts = total_counts,
                 keep = total_counts > 3000)
spe <- spe[, qc$keep]

after <- visium_qc_dimensions(spe, database)

# 2. Annotate the retained spots using the published deconvolution workbook.
# Add annotations to this native-count object; never aggregate duplicate symbols.
positions <- visium_position_table(spe)
um_per_pixel <- S4Vectors::metadata(spe)$um_per_pixel

fields <- list(spot_id = colnames(spe),
  image_x = positions$pxl_col_in_fullres, image_y = positions$pxl_row_in_fullres,
  x_um = positions$pxl_col_in_fullres * um_per_pixel,
  y_um = positions$pxl_row_in_fullres * um_per_pixel,
  grid_x_um = positions$array_col * 50,
  grid_y_um = positions$array_row * 50 * sqrt(3),
  total_counts = total_counts[qc$keep],
  n_genes = as.numeric(Matrix::colSums(SummarizedExperiment::assay(spe, "counts") > 0)))
for (field in names(fields)) {
  SummarizedExperiment::colData(spe)[[field]] <- fields[[field]]
}
workbook <- repo_path("Data", "downloaded", "janesick_2023_springer",
                      "janesick_2023_source_data.xlsx")
deconvolution <- as.data.frame(readxl::read_xlsx(workbook,
  sheet = "Sup. Fig. 10 Visium Deconv."))
spe <- annotate_visium(spe, deconvolution, workbook)

# 3. Register published pathology ROIs into Visium image coordinates.
# Published Xenium ROIs are registered to Visium image pixels. Region assignment
# uses spot centres and retains the original invasive-tumour overlap precedence.
paper <- loadPaperROIs(workbook = workbook)
transform <- as.matrix(utils::read.csv(repo_path("Data", "parameters",
  "from_10x_genomics", "xenium_to_visium_transform.csv"), header = FALSE))
storage.mode(transform) <- "double"
calibration <- jsonlite::fromJSON(repo_path("Data", "parameters",
  "from_10x_genomics", "xenium_coordinate_calibration.json"))

xenium_um_to_visium <- transform %*%
  diag(c(1 / calibration$pixel_size_um, 1 / calibration$pixel_size_um, 1))
centres <- project_xy(cbind(spe$image_x, spe$image_y), solve(xenium_um_to_visium))
spe$region <- classifyPaperRegions(centres, paper$rois)
borders <- paperROIBorderPaths(paper$rois)
borders[, c("x", "y")] <- project_xy(borders[, c("x", "y")], xenium_um_to_visium)

# 4. Preserve calibrated images and source provenance in the processed object.
images <- SpatialExperiment::imgData(spe)
sample_id <- unique(as.character(spe$sample_id))

image_candidates <- which(images$sample_id == sample_id &
                            is.finite(images$scaleFactor) & images$scaleFactor > 0)

image_index <- image_candidates[which.max(images$scaleFactor[image_candidates])]
md <- S4Vectors::metadata(spe)
md$qc <- list(rule = "raw total counts > 3000", before = nrow(qc), after = ncol(spe),
              gene_feature_rows = nrow(spe), feature_aggregation = "none",
              spot_filter = qc, summary = data.frame(metric = names(before),
                before_qc = unname(before), after_qc = unname(after)))
md$source_url <- "https://github.com/ChenLaboratory/example_data/releases/download/v1.0.0/spe_visium_bc_s1.rds"
md$visium <- list(source = "BLISA release SpatialExperiment, native feature rows",
  source_coordinate_unit = md$coord_unit, source_positions_before_qc = original_positions,
  image_id = as.character(images$image_id[image_index]), sample_id = sample_id,
  spot_pitch_um = 100, spot_diameter_um = 55, um_per_fullres_pixel = um_per_pixel,
  scalefactors = list(spot_diameter_fullres = 55 / um_per_pixel,
    tissue_hires_scalef = images$scaleFactor[image_index]),
  coordinate_note = paste("Inference uses the ideal 100-micrometre array lattice.",
    "Original spatialCoords, imgData and calibration are retained for image overlays."))
md$pathology <- list(source = "Janesick et al. Figure 4a polygons",
  assignment = "Retained spot image centres; pooled components; invasive tumour takes precedence",
  regions = paperRegionLevels(), xenium_um_to_visium_pixels = xenium_um_to_visium,
  borders_visium_pixels = borders)
md$deconvolution$source_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-43458-x/MediaObjects/41467_2023_43458_MOESM4_ESM.xlsx"
md$deconvolution$sheet <- "Sup. Fig. 10 Visium Deconv."
md$deconvolution$region_definition <- "Published pathology ROIs; dominant cell types remain in dominant_ct and cell_type"
S4Vectors::metadata(spe) <- md

out <- repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
base::saveRDS(spe, out)
