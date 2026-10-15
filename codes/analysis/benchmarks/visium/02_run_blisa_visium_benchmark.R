# Run the Visium BLISA benchmark after 01a_prepare_visium_benchmark.R.
# Use 10,000 permutations, 40 threads and seed 20260915, as in the recorded run.
# 1. Read the prepared counts and shared single-subunit CellChat candidates.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
suppressPackageStartupMessages(library(SpatialExperiment))
data_dir <- repo_path("Data", "generated", "benchmarks", "visium")
out <- repo_path("res", "benchmarks", "visium")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ncpu <- as.integer(Sys.getenv("NCPU", "40"))
nperm <- 10000L
seed <- 20260915L
radius <- 250
spe <- readRDS(file.path(data_dir, "visium_qc_benchmark_18082_genes.rds"))
lr_df <- read.csv(file.path(data_dir, "visium_benchmark_single_subunit_lr_candidates.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
rownames(lr_df) <- lr_df$cellchat_id

# 2. Retain genes with total counts >=200 AND detection in at least 200 spots.
# External total-count/detection filtering precedes the measured stage.
# The native maximum-count >=3 rule stays inside BLISA's measured call.
counts <- SummarizedExperiment::assay(spe, "counts")
eligible <- rownames(counts)[Matrix::rowSums(counts) >= 200 & Matrix::rowSums(counts > 0) >= 200]
lr_for_blisa <- getFromNamespace("getLRpairs", "blisa")(
  gene_panel = eligible, LR_df = lr_df, species = "human")

# 3. Run native BLISA; measure runtime and Linux process peak RAM.
# Process RSS includes native allocations and threads, excluding separate
# workers. Imports, loading, shared LR preparation and exports are not timed.
process_memory_mib <- function(field) {
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]][2]) / 1024
}
set.seed(seed)
gc(full = TRUE)
rss_peak_reset <- tryCatch({ writeLines("5", "/proc/self/clear_refs"); TRUE },
                           warning = function(w) FALSE, error = function(e) FALSE)
timing <- system.time({
  result <- blisa::blisa(spe, platform = "visium", bin_size = 100, min_total_counts = 0,
    LR_df = lr_for_blisa, group = NULL, nsim = nperm, cpu_threads = ncpu,
    p_cutoff = 0.05, dmax = radius, min_ligand = 3, min_receptor = 3)
}, gcFirst = FALSE)
peak_rss_mib <- if (rss_peak_reset) process_memory_mib("VmHWM") else NA_real_

# 4. Preserve the native object, including local p-values and High-High calls.
saveRDS(result, file.path(out, "visium_benchmark_blisa_result_10000_samples.rds"))
lr <- result$LR_results
lr_ids <- paste(lr$ligand.symbol, lr$receptor.symbol, sep = "__")
summary <- data.frame(lr_id = lr_df$lr_id, returned = lr_df$lr_id %in% lr_ids,
                       n_available = NA_integer_, n_sig = NA_integer_)
summary$native_status <- ifelse(summary$returned, "tested",
  ifelse(lr_df$ligand.symbol %in% eligible & lr_df$receptor.symbol %in% eligible,
          "native_filtered", "custom_filtered"))
# Summarize all candidates for the method comparison, including filtered pairs.
for (i in seq_len(nrow(lr))) {
  p <- lr$all_pval[[i]]
  keep <- result$spatial_weights[[if (lr$ccc_mode[i] == "nearby") "keep_idx_queen" else "keep_idx_dist"]]
  tested <- seq_along(p) %in% keep & is.finite(p)
  summary$n_available[summary$lr_id == lr_ids[i]] <- sum(tested)
  summary$n_sig[summary$lr_id == lr_ids[i]] <- if (any(tested)) length(lr$sig_index[[i]]) else NA_integer_
}
summary$status <- ifelse(!summary$returned, "NR",
                         ifelse(summary$n_available > 0, "available", "NA"))
summary$budget <- nperm
summary$condition <- "blisa"
counts_path <- file.path(out, "visium_lri_hotspot_counts_and_test_status_10000_samples.csv")
summary <- summary[, c("budget", "condition", "lr_id", "returned", "n_available", "n_sig",
                       "status", "native_status")]
all_counts <- if (file.exists(counts_path)) {
  read.csv(counts_path, stringsAsFactors = FALSE)
} else summary[FALSE, ]
all_counts$status[is.na(all_counts$status)] <- "NA"
write_csv(rbind(all_counts[all_counts$condition != "blisa", names(summary)], summary), counts_path)
# 5. Save this method's resource measurements alongside the other methods.
performance <- data.frame(method = "blisa", condition = "blisa", budget = nperm,
  nperm = nperm, ncpu = ncpu, elapsed_seconds = unname(timing["elapsed"]),
  peak_rss_mib = peak_rss_mib, rss_peak_reset = rss_peak_reset,
  rss_scope = "Stage process peak; includes threads, excludes separate workers",
  timing_scope = "Native spatial preprocessing, maximum-count filtering and LR inference; excludes shared database preparation, input loading, external total-count/detection filtering and exports")
performance_path <- file.path(out, "visium_benchmark_computation_resources_10000_samples.csv")
all_performance <- if (file.exists(performance_path)) {
  read.csv(performance_path, stringsAsFactors = FALSE, check.names = FALSE)
} else performance[FALSE, ]
write_csv(rbind(all_performance[!(all_performance$condition == "blisa" & all_performance$budget == nperm),
                                names(performance)], performance), performance_path)
message("BLISA results updated; stLearn/SpatialDM saved files were unchanged.")
