# Compare the two saved BLISA analyses on the same aligned Visium spots.
source("codes/helpers/source_data_helpers.R")
out <- repo_path("res", "aligned_overlap")

# 1. Match LR pairs and spots explicitly, retaining the native Visium LR order.
summary_v <- read.csv(file.path(out, "visium_lri_summary.csv"), check.names = FALSE)
summary_x <- read.csv(file.path(out, "xenium_lri_summary.csv"), check.names = FALSE)
spots_v <- read.csv(file.path(out, "visium_lri_by_spot.csv"), check.names = FALSE)
spots_x <- read.csv(file.path(out, "xenium_lri_by_spot.csv"), check.names = FALSE)
key_v <- paste(spots_v$database_id, spots_v$spot_id, sep = "::")
key_x <- paste(spots_x$database_id, spots_x$spot_id, sep = "::")
spots_x <- spots_x[match(key_v, key_x), ]
summary_x <- summary_x[match(summary_v$database_id, summary_x$database_id), ]
paired <- data.frame(database_id = summary_v$database_id,
  lr_label = summary_v$lr_label, mode = summary_v$mode,
  visium = summary_v$n_sig, xenium = summary_x$n_sig)
write_csv(paired, file.path(out, "aligned_lri_hotspot_counts.csv"))

# 2. Compare LR-level hotspot counts with the manuscript's Spearman correlation.
correlation <- data.frame(method = "spearman", n_lr_pairs = nrow(paired),
  rho = stats::cor(paired$visium, paired$xenium, method = "spearman"))
write_csv(correlation, file.path(out, "aligned_lri_hotspot_count_correlation.csv"))

# 3. Classify spot-level agreement only where both platform tests are available.
common <- spots_v$tested & spots_x$tested
state <- rep("untested", nrow(spots_v))
state[common] <- ifelse(spots_v$significant[common],
  ifelse(spots_x$significant[common], "both", "visium_only"),
  ifelse(spots_x$significant[common], "xenium_only", "neither"))
agreement <- data.frame(database_id = spots_v$database_id,
  spot_id = spots_v$spot_id, state = state)
write_csv(agreement, file.path(out, "aligned_lri_agreement_by_spot.csv"))
