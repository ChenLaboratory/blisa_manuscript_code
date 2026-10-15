# Summarize the saved Xenium LRI/CCI results for Figure 2 and Figures S1-S3.
source("codes/helpers/source_data_helpers.R")

# 1. Read complete native hotspot calls, cell-pair assignments and region labels.
binned <- readRDS(repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds"))
out <- repo_path("res", "xenium")
summary <- read.csv(file.path(out, "lri_summary.csv"), check.names = FALSE)
local <- read.csv(file.path(out, "lri_by_bin.csv"), check.names = FALSE)
dominant <- read.csv(file.path(out, "xenium_dominant_celltype_pairs_by_bin.csv"),
  check.names = FALSE)
regions <- as.character(binned$bins_metadata$region[
  match(binned$bins$bin_id, binned$bins_metadata$bin_id)])
region_order <- c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI")

# 2. Divide native High-High hotspot counts by every retained bin in each region.
# Untested bins remain in this denominator, as in the original regional heatmap.
coverage <- do.call(rbind, lapply(summary$lr_id, function(id) {
  rows <- local[local$lr_id == id, ]
  rows <- rows[match(binned$bins$bin_id, rows$bin_id), ]
  do.call(rbind, lapply(region_order, function(region) {
    keep <- regions == region
    n_hotspots <- sum(rows$native_significant[keep])
    data.frame(lr_id = id, region = region, n_hotspots = n_hotspots,
      n_bins = sum(keep), hotspot_proportion = n_hotspots / sum(keep))
  }))
}))
write_csv(coverage, file.path(out, "xenium_lri_regional_coverage.csv"))

# 3. Rank all dominant sender-receiver pairs by their number of native hotspots.
# Preserve table()'s alphabetical order for ties, matching the original maps.
frequencies <- do.call(rbind, lapply(summary$lr_id, function(id) {
  rows <- dominant[dominant$lr_id == id & dominant$native_significant, ]
  counts <- sort(table(rows$cell_pair), decreasing = TRUE)
  pairs <- rows[match(names(counts), rows$cell_pair), c("sender", "receiver")]
  data.frame(lr_id = rep(id, length(counts)), rank = seq_along(counts), pairs,
    cell_pair = names(counts), n_hotspots = as.integer(counts),
    hotspot_proportion = as.integer(counts) / nrow(rows), row.names = NULL)
}))
write_csv(frequencies, file.path(out, "xenium_dominant_celltype_pair_frequencies.csv"))
