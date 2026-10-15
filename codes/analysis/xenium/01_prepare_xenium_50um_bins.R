# Run from the repository root. Original release data stay in the external cache.
# This script performs preparation only; 02_run_xenium_lri_cci.R runs inference.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "xenium", "xenium_helpers.R"))
suppressPackageStartupMessages({
  library(blisa)
  library(SpatialExperiment)
  library(Matrix)
})

# 1. Read the release cells and harmonize DCIS labels with the published ROIs.
spe <- read_source_data("xenium")
spe$cell_type_broad <- spe$cell_type
spe$cell_type <- splitDCISLabels(spe$cell_type_broad, spe$Cluster)
paper <- loadPaperROIs()
roi_validation <- validatePaperROIs(paper, spe)

coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
counts <- SummarizedExperiment::assay(spe, "counts")
# 2. Aggregate expression into 50-micrometre hexagons.
# Keep the native hexBinCells filtering defaults used in the manuscript.
binned <- blisa::hexBinCells(coords_df = coords, counts_matrix = counts,
  bin_size = 50, group = spe$cell_type)

# 3. Preserve cell counts, proportions and region labels with the processed input.
# Count annotated cells separately from cell-type-specific gene expression.
cell_type <- droplevels(factor(spe$cell_type))
indicators <- Matrix::sparseMatrix(
  i = as.integer(cell_type), j = seq_along(cell_type), x = 1,
  dims = c(nlevels(cell_type), length(cell_type)),
  dimnames = list(levels(cell_type), colnames(counts)))
cell_bins <- blisa::hexBinCells(coords_df = coords, counts_matrix = indicators,
  bin_size = 50, min_cells = 1, min_total_counts = 0)
sample_id <- colnames(binned$counts_matrix)
index <- match(sample_id, colnames(cell_bins$counts_matrix))

numbers <- cell_bins$counts_matrix[, index, drop = FALSE]

proportions <- numbers %*% Matrix::Diagonal(x = 1 / Matrix::colSums(numbers))
dimnames(proportions) <- dimnames(numbers)
number_matrix <- as.matrix(Matrix::t(numbers))
proportion_matrix <- as.matrix(Matrix::t(proportions))
centres <- sf::st_coordinates(sf::st_centroid(binned$bins))
binned$bins_metadata <- data.frame(
  sample_id = sample_id, bin_id = binned$bins$bin_id,
  x = centres[, "X"], y = centres[, "Y"],
  n_cells = binned$bins$n_cells, total_counts = binned$bins$total_counts,
  dominant_ct = colnames(number_matrix)[max.col(number_matrix, ties.method = "first")],
  row.names = sample_id, check.names = FALSE)
binned$bins_metadata[paste0("cell_count::", colnames(number_matrix))] <-
  as.data.frame(number_matrix, check.names = FALSE)
binned$bins_metadata[paste0("cell_proportion::", colnames(proportion_matrix))] <-
  as.data.frame(proportion_matrix, check.names = FALSE)
binned$bins_metadata$region <- annotatePaperBins(binned$bins, paper$rois)
binned$paper_rois <- paper$rois
binned$roi_alignment_audit <- roi_validation$audit
# Figure 2a uses all retained cells, rather than hexagon composition.
binned$cells_metadata <- data.frame(cell_id = colnames(spe),
  x = coords[[1]], y = coords[[2]], cell_type = spe$cell_type,
  cell_type_broad = spe$cell_type_broad, cluster = spe$Cluster,
  region = roi_validation$cells$region)
binned$preparation <- list(bin_size_um = 50,
  source = "spe_xenium_bc_s1rep1.rds", roi_y_offset = paper$y_offset,
  blisa_version = as.character(utils::packageVersion("blisa")))

# 4. Save the processed input used by LRI/CCI analysis and expression plots.
destination <- repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds")
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
saveRDS(binned, destination)
