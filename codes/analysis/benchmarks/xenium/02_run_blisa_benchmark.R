# Optional computational rerun: BLISA only. Not needed for the supplied figures.
# Run manually after 01_import_saved_xenium_benchmark_results.R has staged the input.
# Replaces BLISA rows in the generic result tables; retains stLearn/SpatialDM rows.
source("codes/helpers/source_data_helpers.R")
require_packages("blisa")
library(blisa)
input_dir <- repo_path("Data", "generated", "benchmarks", "xenium")
output <- repo_path("res", "benchmarks", "xenium")
binned <- readRDS(file.path(input_dir, "xenium_benchmark_50um_binned.rds"))
reference_lr <- read.csv(file.path(input_dir, "xenium_benchmark_blisa_lr_input.csv"),
  check.names = FALSE)
local <- read.csv(file.path(output, "xenium_benchmark_lri_by_bin.csv.gz"), check.names = FALSE)
summary <- read.csv(file.path(output, "xenium_benchmark_lr_summary.csv"), check.names = FALSE)
performance <- read.csv(file.path(output, "xenium_benchmark_computation_resources.csv"),
  check.names = FALSE)
ncpu <- 10L
nperm <- 10000L

# The original Linux process-level measurements include native allocations and
# threads, but exclude separate worker processes. Imports, input loading and
# export are outside the measured stage; native database loading is inside it.
process_memory_mib <- function(field) {
  tryCatch({
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    if (length(line) != 1L) return(NA_real_)
    as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]][2]) / 1024
  }, error = function(e) NA_real_)
}
gc_before <- gc(full = TRUE, reset = TRUE)
baseline_rss_mib <- process_memory_mib("VmRSS")
rss_peak_reset <- tryCatch({
  writeLines("5", "/proc/self/clear_refs")
  TRUE
}, warning = function(w) FALSE, error = function(e) FALSE)

# Faithful to bechmark_methods_on_xenium/run_blisa.R: no seed was specified.
# Default LR database handling remains inside the timed native call. The frozen
# LR input reference is checked afterwards, so database drift is not silent.
timing <- system.time({
  result <- blisa(binned$counts_matrix, bins = binned$bins,
    n_cells_col = "n_cells", counts_by_group = NULL,
    nsim = nperm, cpu_threads = ncpu)
}, gcFirst = FALSE)
peak_rss_mib <- if (rss_peak_reset) process_memory_mib("VmHWM") else NA_real_
gc_after <- gc(full = TRUE)
lr <- result$LR_results
stopifnot(setequal(rownames(lr), reference_lr$interaction_name),
  identical(result$bins$bin_id, binned$bins$bin_id))

# Export complete native local tests, with unavailable values kept distinct.
metadata <- binned$bins_metadata
new_local <- do.call(rbind, lapply(seq_len(nrow(lr)), function(i) {
  row <- summary[which(summary$method == "blisa" & summary$native_lr_id == rownames(lr)[i]), ]
  stopifnot(nrow(row) == 1L)
  available <- is.finite(lr$all_pval[[i]])
  data.frame(method = "blisa", lr_id = row$lr_id, sample_id = metadata$sample_id,
    p_value = lr$all_pval[[i]], local_available = available,
    significant = seq_len(nrow(metadata)) %in% lr$sig_index[[i]],
    status = ifelse(available, "tested", "excluded_spatial_weights"),
    native_lr_id = rownames(lr)[i], score = lr$all_lisa[[i]],
    quadrant = lr$all_quadrant[[i]], bin_id = metadata$bin_id)
}))
local <- rbind(new_local[, names(local)], local[local$method != "blisa", ])
for (i in which(summary$method == "blisa" & summary$returned)) {
  rows <- new_local[new_local$lr_id == summary$lr_id[i], ]
  summary$n_available[i] <- sum(rows$local_available)
  summary$significant_bins[i] <- if (any(rows$local_available)) sum(rows$significant) else NA_integer_
  summary$status[i] <- if (any(rows$local_available)) "available" else "no_local_results"
}
# Re-rank the full LR union by the new BLISA result; absent pairs remain last.
rank_rows <- summary[summary$method == "blisa", ]
rank_ids <- rank_rows$lr_id[order(-rank_rows$significant_bins, rank_rows$rank, na.last = TRUE)]
summary$rank <- match(summary$lr_id, rank_ids)

measured <- list(package_version = as.character(packageVersion("blisa")),
  stage = "lr_analysis", ncpu = ncpu, nperm = nperm,
  n_genes = nrow(binned$counts_matrix), n_bins = ncol(binned$counts_matrix),
  n_lr_pairs = nrow(lr), elapsed_seconds = unname(timing["elapsed"]),
  user_cpu_seconds = unname(timing["user.self"]), system_cpu_seconds = unname(timing["sys.self"]),
  baseline_rss_mib = baseline_rss_mib, peak_rss_mib = peak_rss_mib,
  extra_rss_mib = max(0, peak_rss_mib - baseline_rss_mib), rss_peak_reset = rss_peak_reset,
  baseline_r_heap_mib = sum(gc_before[, 2]),
  peak_r_heap_upper_bound_mib = sum(gc_after[, ncol(gc_after)]),
  extra_r_heap_upper_bound_mib = sum(gc_after[, ncol(gc_after)]) - sum(gc_before[, 2]),
  result_size_mib = as.numeric(object.size(result)) / 1024^2,
  source_dir = "codes/analysis/benchmarks/xenium/02_run_blisa_benchmark.R",
  reused = FALSE, input_dir = "Data/generated/benchmarks/xenium",
  run_dir = "res/benchmarks/xenium",
  elapsed_minutes = unname(timing["elapsed"]) / 60, peak_rss_gib = peak_rss_mib / 1024)
for (name in names(measured)) performance[performance$method == "blisa", name] <- measured[[name]]
write_csv(summary, file.path(output, "xenium_benchmark_lr_summary.csv"))
write_csv(performance, file.path(output, "xenium_benchmark_computation_resources.csv"))
local$p_value <- ifelse(is.na(local$p_value), NA_character_, sprintf("%.17g", local$p_value))
local$score <- ifelse(is.na(local$score), NA_character_, sprintf("%.17g", local$score))
connection <- gzfile(file.path(output, "xenium_benchmark_lri_by_bin.csv.gz"), "wt")
write.csv(local, connection, row.names = FALSE, na = "NA")
close(connection)
