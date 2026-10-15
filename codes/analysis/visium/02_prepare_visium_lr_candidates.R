# Prepare the complete native CellChat candidate table from the annotated QC input.
# External thresholds apply to every subunit; BLISA applies its maximum-count filter later.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "Matrix"))
# 1. Read the native-feature Visium QC object and downloaded CellChat snapshot.
spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))

database <- load_visium_cellchat()
# 2. Keep genes with at least 200 counts across at least 200 nonempty spots.
raw_counts <- SummarizedExperiment::assay(spe, "counts")
eligible <- Matrix::rowSums(raw_counts) >= 200 & Matrix::rowSums(raw_counts > 0) >= 200
eligible_genes <- unique(rownames(spe)[eligible])
# Keep complete ligand/receptor definitions whose every subunit passes the filter.
filtered_lr <- getFromNamespace("getLRpairs", "blisa")(
  gene_panel = eligible_genes, LR_df = database)

# 3. Save the candidates consumed by the native LRI analysis.
write_csv(data.frame(native_lr_id = rownames(filtered_lr), filtered_lr,
  row.names = NULL, check.names = FALSE), repo_path("Data", "generated",
  "ligand_receptor_tables", "visium_expression_filtered_lr_candidates.csv"))
