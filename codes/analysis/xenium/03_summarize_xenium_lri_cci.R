# Calculate regional hotspot proportions from the saved native Xenium result.
source("codes/helpers/source_data_helpers.R")

# 1. Match processed region labels to the native result's bin order.
binned <- readRDS(repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds"))
out <- repo_path("res", "xenium")
result <- readRDS(file.path(out, "xenium_blisa_result.rds"))
lr <- result$LR_results
regions <- as.character(binned$bins_metadata$region[
  match(result$bins$bin_id, binned$bins_metadata$bin_id)])
region_order <- c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI")

# 2. Divide native High-High hotspot counts by every retained bin in each region.
# Untested bins remain in this denominator, as in the original regional heatmap.
coverage <- do.call(rbind, lapply(seq_len(nrow(lr)), function(i) {
  hotspots <- seq_len(nrow(result$bins)) %in% lr$sig_index[[i]]
  do.call(rbind, lapply(region_order, function(region) {
    keep <- regions == region
    n_hotspots <- sum(hotspots[keep])
    data.frame(lr_id = rownames(lr)[i], region = region, n_hotspots = n_hotspots,
      n_bins = sum(keep), hotspot_proportion = n_hotspots / sum(keep))
  }))
}))
write_csv(coverage, file.path(out, "xenium_lri_regional_coverage.csv"))
