# Prepare the complete native CellChat candidate table from the annotated QC input.
# External thresholds apply to every subunit; BLISA applies its maximum-count filter later.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "Matrix",
                   "digest"))
spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))
stopifnot(nrow(spe) == 18085L, ncol(spe) == 4412L)
database <- load_visium_cellchat()
coverage <- lr_panel_coverage(database, rownames(spe))
preparation <- system.time({
  raw_counts <- SummarizedExperiment::assay(spe, "counts")
  gene_audit <- data.frame(feature_id = as.character(SummarizedExperiment::rowData(spe)$ID),
    gene = rownames(spe), total_count = as.numeric(Matrix::rowSums(raw_counts)),
    nonempty_spots = as.numeric(Matrix::rowSums(raw_counts > 0)))
  gene_audit$passes_external <- with(gene_audit, total_count >= 200 & nonempty_spots >= 200)
  eligible_genes <- unique(gene_audit$gene[gene_audit$passes_external])
  filtered_lr <- getFromNamespace("getLRpairs", "blisa")(
    gene_panel = eligible_genes, LR_df = database)
  predicted_lr <- getFromNamespace("filterLRpairs", "blisa")(
    raw_counts, min_ligand = 3, min_receptor = 3, LR_df = filtered_lr)
  lr_audit <- data.frame(native_lr_id = rownames(database), database,
                         row.names = NULL, check.names = FALSE)
  lr_audit$panel_covered <- coverage$panel_covered
  lr_audit$passes_external <- lr_audit$native_lr_id %in% rownames(filtered_lr)
  lr_audit$custom_eligible <- lr_audit$native_lr_id %in% rownames(predicted_lr)
})
stopifnot(nrow(filtered_lr) == 1917L, nrow(predicted_lr) == 1616L)
write_csv(data.frame(native_lr_id = rownames(filtered_lr), filtered_lr,
  row.names = NULL, check.names = FALSE), repo_path("Data", "generated",
  "ligand_receptor_tables", "visium_expression_filtered_lr_candidates.csv"))
