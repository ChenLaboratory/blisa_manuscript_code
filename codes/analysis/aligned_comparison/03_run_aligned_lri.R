#!/usr/bin/env Rscript
# Two matched-platform analyses for Figure 3f-h and Figure S5.
# No 50-micrometre third arm, plotting, downloads, or native-result RDS exports.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "aligned_comparison", "alignment_helpers.R"))
suppressPackageStartupMessages({library(SpatialExperiment); library(blisa)})

# 1. Restrict both overlap inputs to shared genes and identical spot coordinates.
data_dir <- repo_path("Data", "generated", "aligned_overlap")
v <- readRDS(file.path(data_dir, "visium_aligned_overlap_qc.rds"))
xb <- readRDS(file.path(data_dir, "xenium_aggregated_to_visium_overlap.rds"))
shared_genes <- intersect(rownames(xb), rownames(v))

v <- v[shared_genes, ]
spatial_fields <- c("in_tissue", "array_row", "array_col", "sample_id",
  "image_x", "image_y", "grid_x_um", "grid_y_um", "n_xenium_cells")
colData(v) <- colData(v)[, spatial_fields]
assays(v) <- S4Vectors::SimpleList(counts = counts(v))
imgData(v) <- S4Vectors::DataFrame()
metadata(v) <- list(coord_unit = "pixel")
v$total_counts <- as.numeric(Matrix::colSums(counts(v)))
v$n_genes <- as.integer(Matrix::colSums(counts(v) > 0))
xv <- v
assays(xv) <- S4Vectors::SimpleList(counts = counts(xb[shared_genes, ]))
rowData(xv) <- rowData(xb[shared_genes, ])
xv$total_counts <- as.numeric(Matrix::colSums(counts(xv)))
xv$n_genes <- as.integer(Matrix::colSums(counts(xv) > 0))

# 2. Read the shared CellChat LR candidate table.
lr_file <- repo_path("Data", "generated", "ligand_receptor_tables", "aligned_shared_lr_pairs.csv")
shared_lr <- read.csv(lr_file, check.names = FALSE, stringsAsFactors = FALSE)

rownames(shared_lr) <- shared_lr$database_id
shared_lr$database_id <- NULL

# 3. Run the two original BLISA calls: 9,999 permutations, four threads,
# P <=0.05, 250-micrometre diffuse neighbourhood, no further abundance filters.
# The source uses fastLISA's native default seed 123456789, not R's set.seed().
# Pin BLISA/fastLISA to the recorded manuscript revisions when reproducing.
run_arm <- function(spe) {
  result <- blisa(spe, platform = "visium", bin_size = 100,
    min_total_counts = 0, group = NULL, LR_df = shared_lr,
    nsim = 9999, cpu_threads = 4L, p_cutoff = 0.05, dmax = 250,
    min_ligand = 0, min_receptor = 0)

  result$bins$spot_id <- colnames(spe)
  result
}
visium_result <- run_arm(v)
xenium_result <- run_arm(xv)

# 4. Save complete native results for each platform.
# Long tables preserve every LR x retained-spot observation. Native sig_index
# defines High-High hotspots; isolated/unavailable observations remain untested.
out <- repo_path("res", "aligned_overlap")
for (name in c("visium", "xenium")) {
  result <- if (name == "visium") visium_result else xenium_result
  summary <- summarize_blisa(result, name)

  summary$rank <- seq_len(nrow(summary))
  write_csv(summary, file.path(out, paste0(name, "_lri_summary.csv")))
  local <- blisa_spot_table(result, name)
  local <- local[, c("database_id", "spot_id", "p_value", "lisa", "quadrant",
    "tested", "significant")]
  write_csv(local, file.path(out, paste0(name, "_lri_by_spot.csv")))
}
