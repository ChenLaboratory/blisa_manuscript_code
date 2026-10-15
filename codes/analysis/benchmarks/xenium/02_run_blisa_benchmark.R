# Run BLISA on the main Xenium analysis's prepared 50-um bins.
# Replaces BLISA rows while retaining existing stLearn/SpatialDM results.
source("codes/helpers/source_data_helpers.R")
library(blisa)
# 1. Read raw binned counts and any already completed method results.
output <- repo_path("res", "benchmarks", "xenium")
binned <- readRDS(repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds"))
read_result <- function(name) {
  path <- file.path(output, name)
  if (!file.exists(path)) return(data.frame(method = character()))
  table <- read.csv(path, check.names = FALSE)
  for (column in intersect(c("returned", "local_available", "significant"), names(table))) {
    table[[column]] <- as.logical(table[[column]])
  }
  table
}
bind_results <- function(a, b) {
  columns <- union(names(a), names(b))
  for (name in setdiff(columns, names(a))) a[[name]] <- rep(NA, nrow(a))
  for (name in setdiff(columns, names(b))) b[[name]] <- rep(NA, nrow(b))
  rbind(a[, columns, drop = FALSE], b[, columns, drop = FALSE])
}
local <- read_result("xenium_benchmark_lri_by_bin.csv.gz")
summary <- read_result("xenium_benchmark_lr_summary.csv")
performance <- read_result("xenium_benchmark_computation_resources.csv")
local_columns <- c("method", "lr_id", "sample_id", "p_value", "local_available",
  "significant", "lr_score", "quadrant", "bin_id")
summary_columns <- c("method", "lr_id", "rank", "returned", "n_available",
  "significant_bins", "status")
resource_columns <- c("method", "condition", "nperm", "ncpu", "elapsed_seconds",
  "peak_rss_mib", "rss_peak_reset", "rss_scope", "timing_scope")
ncpu <- 10L
nperm <- 10000L

# 2. Start Linux process-level resource measurement. This includes allocations and
# threads, but excludes separate worker processes. Imports, input loading and
# export are outside the measured stage; native database loading is inside it.
process_memory_mib <- function(field) {
  tryCatch({
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]][2]) / 1024
  }, error = function(e) NA_real_)
}
gc(full = TRUE, reset = TRUE)
rss_peak_reset <- tryCatch({
  writeLines("5", "/proc/self/clear_refs")
  TRUE
}, warning = function(w) FALSE, error = function(e) FALSE)

# 3. Run the original native call: 10,000 permutations, 10 threads, no R seed.
# Default LR database handling remains inside the timed native call.
# LRI uses pooled bin counts; the unchanged DCIS cell-type labels are not inputs.
timing <- system.time({
  result <- blisa(binned$counts_matrix, bins = binned$bins,
    n_cells_col = "n_cells", counts_by_group = NULL,
    nsim = nperm, cpu_threads = ncpu)
}, gcFirst = FALSE)
peak_rss_mib <- if (rss_peak_reset) process_memory_mib("VmHWM") else NA_real_
gc(full = TRUE)
lr <- result$LR_results

# 4. Export every LR/bin test and summarize the native High-High hotspots.
metadata <- binned$bins_metadata
# BLISA canonicalizes native row names (CLEC2C_KLRB1 becomes CD69_KLRB1),
# so construct comparison IDs from the returned ligand/receptor symbols.
new_local <- do.call(rbind, lapply(seq_len(nrow(lr)), function(i) {
  # Native vectors include finite placeholders for bins excluded by spatial
  # weights. Use the native keep indices to distinguish these from tested bins.
  keep <- result$spatial_weights[[if (lr$ccc_mode[i] == "nearby")
    "keep_idx_queen" else "keep_idx_dist"]]
  available <- seq_len(nrow(metadata)) %in% keep & is.finite(lr$all_pval[[i]]) &
    is.finite(lr$all_lisa[[i]])
  data.frame(method = "blisa", lr_id = paste(lr$ligand.symbol[i],
    lr$receptor.symbol[i], sep = "__"), sample_id = metadata$sample_id,
    p_value = ifelse(available, lr$all_pval[[i]], NA_real_), local_available = available,
    significant = seq_len(nrow(metadata)) %in% lr$sig_index[[i]],
    lr_score = NA_real_,
    quadrant = lr$all_quadrant[[i]], bin_id = metadata$bin_id)
}))
local <- bind_results(new_local, local[local$method != "blisa", , drop = FALSE])
new_summary <- do.call(rbind, lapply(seq_len(nrow(lr)), function(i) {
  id <- paste(lr$ligand.symbol[i], lr$receptor.symbol[i], sep = "__")
  rows <- new_local[new_local$lr_id == id, ]
  data.frame(method = "blisa", lr_id = id, rank = NA_integer_, returned = TRUE,
    n_available = sum(rows$local_available),
    significant_bins = if (any(rows$local_available)) sum(rows$significant) else NA_integer_,
    status = if (any(rows$local_available)) "available" else "no_local_results")
}))
# 5. Merge BLISA into the three-method comparison and rank by hotspot counts.
old_order <- if (nrow(summary)) unique(summary$lr_id[order(summary$rank)]) else character()
union_ids <- unique(c(old_order, new_summary$lr_id))
measured_methods <- union(summary$method[summary$status != "not_run"], "blisa")
summary <- bind_results(new_summary, summary[summary$method != "blisa", , drop = FALSE])
for (id in union_ids) for (method in c("blisa", "stlearn", "spatialdm")) {
  if (!any(summary$method == method & summary$lr_id == id)) {
    missing <- new_summary[NA_integer_, , drop = FALSE]
    missing$method <- method
    missing$lr_id <- id
    missing$returned <- FALSE
    missing$status <- if (method %in% measured_methods) "pair_not_returned" else "not_run"
    summary <- bind_results(summary, missing)
  }
}
# Re-rank the full LR union by the new BLISA result; absent pairs remain last.
rank_rows <- summary[summary$method == "blisa", ]
rank_ids <- rank_rows$lr_id[order(-rank_rows$significant_bins,
  match(rank_rows$lr_id, union_ids), na.last = TRUE)]
summary$rank <- match(summary$lr_id, rank_ids)

# 6. Save resource measurements and complete unsplit CSV result tables.
measured <- list(method = "blisa", condition = "blisa",
  ncpu = ncpu, nperm = nperm, elapsed_seconds = unname(timing["elapsed"]),
  peak_rss_mib = peak_rss_mib, rss_peak_reset = rss_peak_reset,
  rss_scope = "Stage process peak; includes threads, excludes separate workers",
  timing_scope = "Native spatial preprocessing, database loading and LR inference; excludes input loading and exports")
performance <- bind_results(performance[performance$method != "blisa", , drop = FALSE],
  as.data.frame(measured, stringsAsFactors = FALSE))
write_csv(summary[, summary_columns], file.path(output, "xenium_benchmark_lr_summary.csv"))
write_csv(performance[, resource_columns], file.path(output, "xenium_benchmark_computation_resources.csv"))
local <- local[, local_columns]
local$p_value <- ifelse(is.na(local$p_value), NA_character_, sprintf("%.17g", local$p_value))
local$lr_score <- ifelse(is.na(local$lr_score), NA_character_, sprintf("%.17g", local$lr_score))
connection <- gzfile(file.path(output, "xenium_benchmark_lri_by_bin.csv.gz"), "wt")
write.csv(local, connection, row.names = FALSE, na = "NA")
close(connection)
