# Reuse the QC and annotated Visium input from the main Figure 3 analysis.
# The historical benchmark sums duplicate gene symbols (18,082 genes), whereas
# the main Figure 3 analysis retains all 18,085 feature rows.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(Matrix)
})

# 1. Read the retained spots, coordinates, deconvolution and registered pathology.
# QC and annotations are produced once by visium/01_prepare_visium_qc_annotations.R.
spe <- readRDS(repo_path("Data", "generated", "visium",
                         "visium_qc_counts_gt3000_annotated.rds"))
counts <- SummarizedExperiment::assay(spe, "counts")

# 2. Sum duplicate symbols, preserving first-occurrence gene order and all counts.
# This is the aggregation used in the historical R/Python exchange input.
features <- SummarizedExperiment::rowData(spe)
symbols <- as.character(features$Symbol)
genes <- unique(symbols)
mapping <- Matrix::sparseMatrix(i = match(symbols, genes), j = seq_along(symbols),
  x = 1, dims = c(length(genes), length(symbols)))
aggregated <- methods::as(mapping %*% counts, "dgCMatrix")
dimnames(aggregated) <- list(genes, colnames(spe))
spe <- spe[match(genes, symbols), ]
SummarizedExperiment::assay(spe, "counts") <- aggregated
feature_ids <- split(as.character(features$ID), factor(symbols, levels = genes))
SummarizedExperiment::rowData(spe)$ID <- vapply(feature_ids,
  function(ids) paste(unique(ids), collapse = ";"), character(1))
SummarizedExperiment::rowData(spe)$gene_symbol <- genes
SummarizedExperiment::rowData(spe)$gene_id <- SummarizedExperiment::rowData(spe)$ID
spe$source_n_genes <- as.numeric(Matrix::colSums(counts > 0))
spe$n_genes <- as.numeric(Matrix::colSums(aggregated > 0))

# 3. Record the benchmark gene aggregation and reuse image and ROI calibration.
# Original QC, deconvolution and pathology metadata remain attached to the input.
md <- S4Vectors::metadata(spe)
md$feature_aggregation <- list(rule = "Sum raw counts by gene symbol",
  source_features = length(symbols), exported_genes = length(genes),
  duplicated_symbols = unique(symbols[duplicated(symbols)]))
scales <- md$visium$scalefactors
borders <- md$pathology$borders_visium_pixels
S4Vectors::metadata(spe) <- md

# 4. Select input-present, single-subunit CellChat pairs before method filtering.
# Use primary ligand/receptor names for all methods, retaining source symbols.
database_file <- repo_path("Data", "downloaded", "cellchat_github", "CellChatDB.human.rda")
database <- new.env()
load(database_file, envir = database)
db <- database$CellChatDB.human
interactions <- db$interaction
keep <- !interactions$ligand %in% rownames(db$complex) &
  !interactions$receptor %in% rownames(db$complex) &
  interactions$ligand %in% genes & interactions$receptor %in% genes
lr <- interactions[keep, , drop = FALSE]
lr <- lr[!duplicated(lr[, c("ligand", "receptor")]), , drop = FALSE]
lr$cellchat_id <- rownames(lr)
lr$lr_id <- paste(lr$ligand, lr$receptor, sep = "__")
lr$source_ligand.symbol <- lr$ligand.symbol
lr$source_receptor.symbol <- lr$receptor.symbol
lr$ligand.symbol <- as.character(lr$ligand)
lr$receptor.symbol <- as.character(lr$receptor)
lr$cellchat_source_url <- "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.human.rda"
lr$cellchat_source_sha256 <- digest::digest(file = database_file, algo = "sha256")
lr$native_lr_id <- paste(lr$ligand, lr$receptor, sep = "_")
lr$single_gene <- TRUE

# 5. Save benchmark gene counts, candidates and reused annotations when executed.
# The benchmark border CSV uses negative image y, as in the report panels.
out <- repo_path("Data", "generated", "benchmarks", "visium")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base::saveRDS(spe, file.path(out, "visium_qc_benchmark_18082_genes.rds"))
write_csv(lr, file.path(out, "visium_benchmark_single_subunit_lr_candidates.csv"))
spots <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
spots <- spots[, vapply(spots, function(x) is.atomic(x) && is.null(dim(x)), logical(1))]
spots$x <- spots$grid_x_um
spots$y <- spots$grid_y_um
write_csv(spots, file.path(out, "visium_benchmark_spot_annotations.csv"))
borders$y <- -borders$y
write_csv(borders, file.path(out, "visium_benchmark_registered_roi_borders.csv"))
jsonlite::write_json(scales, file.path(out, "visium_benchmark_image_scalefactors.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
