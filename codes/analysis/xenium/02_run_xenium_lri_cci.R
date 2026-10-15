# Run after 01_prepare_xenium_50um_bins.R, from the repository root.
# Faithful manuscript settings: 9,999 permutations and four CPU threads.
# No random seed was recorded for the original Xenium run; no new seed is imposed.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "xenium", "xenium_helpers.R"))
suppressPackageStartupMessages(library(blisa))

# 1. Load the prepared 50-micrometre bin counts and cell-type expression.
input <- repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds")
binned <- readRDS(input)
# 2. Run the native LRI permutation test and cell-type communication scores.
# Retain the original native CellChat loading and LR-filtering defaults.
# The historical run used blisa 1.0.2; see software_versions.csv.
result <- blisa::blisa(binned$counts_matrix, bins = binned$bins,
  n_cells_col = "n_cells", counts_by_group = binned$counts_by_group,
  nsim = 9999, cpu_threads = 4L)

# 3. Export full CSV results and retain reusable spatial geometry with the input.
binned <- write_xenium_results(result, binned)
saveRDS(binned, input)
