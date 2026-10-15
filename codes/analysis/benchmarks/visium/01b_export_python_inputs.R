# Export the prepared benchmark counts for stLearn and SpatialDM.
# Run 01a_prepare_visium_benchmark.R first, then 01c_prepare_python_inputs.py.
# 1. Read the duplicate-symbol-collapsed input and its complete spot metadata.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
suppressPackageStartupMessages(library(SpatialExperiment))
data_dir <- repo_path("Data", "generated", "benchmarks", "visium")
spe <- readRDS(file.path(data_dir, "visium_qc_benchmark_18082_genes.rds"))
counts <- SummarizedExperiment::assay(spe, "counts")
spots <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
spots <- spots[, vapply(spots, function(x) is.atomic(x) && is.null(dim(x)), logical(1)), drop = FALSE]

# 2. Write sparse raw counts, genes and coordinates without normalization.
out <- file.path(data_dir, "python_inputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Matrix::writeMM(counts, file.path(out, "counts.mtx"))
write_csv(data.frame(gene_symbol = rownames(spe)), file.path(out, "genes.csv"))
write_csv(spots, file.path(out, "spots.csv"))
message("Raw counts and spot metadata saved; next run 01c_prepare_python_inputs.py.")
