# Export the main Xenium analysis's 50-um bins for stLearn and SpatialDM.
# Run after analysis/xenium/01_prepare_xenium_50um_bins.R, then run 01b.
# No new binning, filtering, normalization, or method inference is performed.
source("codes/helpers/source_data_helpers.R")
# 1. Reuse the main analysis's counts and metadata, including separate DCIS labels.
binned <- readRDS(repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds"))
input <- repo_path("Data", "generated", "benchmarks", "xenium")
dir.create(input, recursive = TRUE, showWarnings = FALSE)
counts <- binned$counts_matrix
metadata <- binned$bins_metadata
# 2. Write one raw-count row per bin, preserving gene and bin order.
connection <- gzfile(file.path(input, "xenium_benchmark_raw_counts.csv.gz"), "wt")
write.csv(data.frame(sample_id = colnames(counts), as.matrix(Matrix::t(counts)),
  check.names = FALSE), connection, row.names = FALSE)
close(connection)
write_csv(metadata, file.path(input, "xenium_benchmark_bin_metadata.csv"))
cat("Exported", ncol(counts), "bins x", nrow(counts), "genes without changing counts.\n")
