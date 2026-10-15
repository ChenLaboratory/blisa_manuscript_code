# Compare native 10,000-permutation Xenium results without rerunning inference.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))

# 1. Read each method's saved calls and the processed bin-region annotations.
result_dir <- repo_path("res", "benchmarks", "xenium")
summary <- read.csv(file.path(result_dir, "xenium_benchmark_lr_summary.csv"))
local <- read.csv(file.path(result_dir, "xenium_benchmark_lri_by_bin.csv.gz"))
# Older cached CSVs also contain BLISA rows; use the native object for BLISA.
local <- local[local$method != "blisa", ]
result <- readRDS(file.path(result_dir, "xenium_benchmark_blisa_result.rds"))
bins <- readRDS(repo_path("Data", "generated", "xenium",
                         "xenium_50um_binned.rds"))$bins_metadata
methods <- c("blisa", "stlearn", "spatialdm")
ids <- unique(summary$lr_id[order(summary$rank)])

# 2. Shared eligibility requires an available local test in all three methods.
# It is independent of whether any bin is significant.
common_ids <- Reduce(intersect, lapply(methods, function(method) {
  summary$lr_id[summary$method == method & summary$returned & summary$n_available > 0]
}))
sets <- split(local$sample_id[local$significant],
              paste(local$method[local$significant], local$lr_id[local$significant], sep = "::"))
sample_ids <- bins$sample_id[match(result$bins$bin_id, bins$bin_id)]
lr_ids <- paste(result$LR_results$ligand.symbol, result$LR_results$receptor.symbol, sep = "__")
blisa_sets <- setNames(lapply(result$LR_results$sig_index, function(index) sample_ids[index]),
                      paste("blisa", lr_ids, sep = "::"))
sets <- c(blisa_sets, sets)

# 3. Compare significant-bin sets for every LR and ordered method pair.
# An empty union has an undefined Jaccard score; unavailable tests remain missing.
overlap <- do.call(rbind, lapply(ids, function(id) {
  do.call(rbind, lapply(methods, function(a) {
    do.call(rbind, lapply(methods, function(b) {
      x <- summary[summary$method == a & summary$lr_id == id, ]
      y <- summary[summary$method == b & summary$lr_id == id, ]
      available <- x$returned && y$returned && x$n_available > 0 && y$n_available > 0
      sa <- sets[[paste(a, id, sep = "::")]]
      sb <- sets[[paste(b, id, sep = "::")]]
      numerator <- if (available) length(intersect(sa, sb)) else NA_integer_
      denominator <- if (available) length(union(sa, sb)) else NA_integer_
      data.frame(lr_id = id, method_a = a, method_b = b,
        n_a = x$significant_bins, n_b = y$significant_bins,
        intersection = numerator, union = denominator,
        jaccard = if (available && denominator > 0) numerator / denominator else NA_real_,
        status = if (!x$returned || !y$returned) "NR" else if (!available) "NA" else
          if (denominator == 0) "no_significant_bins" else "available",
        common_three_eligible = id %in% common_ids)
    }))
  }))
}))

# 4. Pool intersections and unions across common LRs, preserving LR-bin identity.
pooled <- aggregate(cbind(intersection, union) ~ method_a + method_b,
                    overlap[overlap$common_three_eligible, ], sum)
pooled$n_lr_eligible <- length(common_ids)
pooled$jaccard <- ifelse(pooled$union > 0, pooled$intersection / pooled$union, NA_real_)
pooled$status <- ifelse(pooled$union > 0, "available", "no_significant_bins")

# 5. Count hotspots in every annotated region, including the area outside ROIs.
# The denominator is every bin in that region, not just bins with an available test.
regions <- c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI")
regional <- do.call(rbind, lapply(seq_len(nrow(summary)), function(i) {
  row <- summary[i, ]
  available <- row$returned && row$n_available > 0
  selected <- sets[[paste(row$method, row$lr_id, sep = "::")]]
  do.call(rbind, lapply(regions, function(region) {
    region_ids <- bins$sample_id[bins$region == region]
    hits <- if (available) sum(region_ids %in% selected) else NA_integer_
    data.frame(method = row$method, lr_id = row$lr_id, region = region,
      n_total = length(region_ids), n_sig = hits, proportion = hits / length(region_ids),
      status = if (!row$returned) "NR" else if (available) "available" else "NA")
  }))
}))
write_csv(overlap, file.path(result_dir, "xenium_benchmark_jaccard_by_lr.csv"))
write_csv(pooled, file.path(result_dir, "xenium_benchmark_jaccard_pooled.csv"))
write_csv(regional, file.path(result_dir, "xenium_benchmark_regional_coverage.csv"))
