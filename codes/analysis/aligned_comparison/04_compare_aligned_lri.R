# Compare the two saved native BLISA results on the same aligned Visium spots.
source("codes/helpers/source_data_helpers.R")
out <- repo_path("res", "aligned_overlap")

# 1. Match LR pairs and spots explicitly, retaining the native Visium LR order.
visium <- readRDS(file.path(out, "visium_blisa_result.rds"))
xenium <- readRDS(file.path(out, "xenium_blisa_result.rds"))
lr_v <- visium$LR_results
lr_x <- xenium$LR_results
lr_index <- match(lr_v$interaction_name, lr_x$interaction_name)
spot_index <- match(visium$bins$spot_id, xenium$bins$spot_id)
paired <- data.frame(database_id = lr_v$interaction_name,
  lr_label = paste(lr_v$ligand.symbol, lr_v$receptor.symbol, sep = " - "),
  mode = lr_v$ccc_mode, visium = lr_v$sig_numbers,
  xenium = lr_x$sig_numbers[lr_index])
write_csv(paired, file.path(out, "aligned_lri_hotspot_counts.csv"))

# 2. Compare LR-level hotspot counts with the manuscript's Spearman correlation.
correlation <- data.frame(method = "spearman", n_lr_pairs = nrow(paired),
  rho = stats::cor(paired$visium, paired$xenium, method = "spearman"))
write_csv(correlation, file.path(out, "aligned_lri_hotspot_count_correlation.csv"))

# 3. Classify spot-level agreement only where both platform tests are available.
# Native weight indices exclude isolated spots even if their placeholder P is 1.
local_flags <- function(result, i) {
  lr <- result$LR_results
  key <- if (lr$ccc_mode[i] == "nearby") "keep_idx_queen" else "keep_idx_dist"
  spots <- seq_len(nrow(result$bins))
  list(tested = spots %in% result$spatial_weights[[key]] &
    is.finite(lr$all_pval[[i]]) & is.finite(lr$all_lisa[[i]]),
    significant = spots %in% lr$sig_index[[i]])
}
agreement <- do.call(rbind, lapply(seq_len(nrow(lr_v)), function(i) {
  v <- local_flags(visium, i)
  x <- local_flags(xenium, lr_index[i])
  common <- v$tested & x$tested[spot_index]
  state <- rep("untested", length(common))
  state[common] <- ifelse(v$significant[common],
    ifelse(x$significant[spot_index][common], "both", "visium_only"),
    ifelse(x$significant[spot_index][common], "xenium_only", "neither"))
  data.frame(database_id = lr_v$interaction_name[i],
    spot_id = visium$bins$spot_id, state = state)
}))
write_csv(agreement, file.path(out, "aligned_lri_agreement_by_spot.csv"))
