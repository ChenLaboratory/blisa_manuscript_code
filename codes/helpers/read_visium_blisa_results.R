# Read a saved native BLISA object for the Python comparison and plot scripts.
# Arguments: result RDS, spot annotations CSV, temporary output CSV, optional LR IDs.
# This only reshapes saved p-values and hotspot calls; it runs no inference.
args <- commandArgs(trailingOnly = TRUE)
result <- readRDS(args[1])
spots <- read.csv(args[2], stringsAsFactors = FALSE)
spot_ids <- spots$spot_id[as.integer(result$bins$bin_id)]
lr <- result$LR_results
ids <- paste(lr$ligand.symbol, lr$receptor.symbol, sep = "__")
selected <- seq_len(nrow(lr))
if (length(args) > 3L) selected <- which(ids %in% args[-seq_len(3)])

# Preserve the native keep mask and significant High-High indices.
connection <- file(args[3], "wt")
for (j in seq_along(selected)) {
  i <- selected[j]
  p <- lr$all_pval[[i]]
  keep <- result$spatial_weights[[if (lr$ccc_mode[i] == "nearby")
    "keep_idx_queen" else "keep_idx_dist"]]
  tested <- seq_along(p) %in% keep & is.finite(p)
  significant <- seq_along(p) %in% lr$sig_index[[i]]
  local <- data.frame(lr_id = ids[i], spot_id = spot_ids,
    p_value = ifelse(tested, sprintf("%.17g", p), NA_character_),
    tested = as.integer(tested),
    significant = ifelse(tested, as.integer(significant), NA_integer_))
  write.table(local, connection, sep = ",", row.names = FALSE,
    col.names = j == 1L, na = "NA")
}
close(connection)
