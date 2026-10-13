# Optional BLISA inference rerun, extracted from the historical Visium runner.
# This was NOT executed while reproducing the manuscript figures.
# Default: 10,000 permutations, 40 threads, seed 20260915 (historical run).
# Explicitly running this script replaces BLISA's saved results/measurements and
# plotting derives comparisons from those results. Imported stLearn/SpatialDM files stay unchanged.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "Matrix"))
suppressPackageStartupMessages(library(SpatialExperiment))
data_dir <- repo_path("Data", "generated", "benchmarks", "visium")
out <- repo_path("res", "benchmarks", "visium")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ncpu <- as.integer(Sys.getenv("NCPU", "40"))
nperm <- 10000L
seed <- 20260915L
radius <- 250
stopifnot(ncpu > 0L)
spe <- readRDS(file.path(data_dir, "visium_qc_benchmark_18082_genes.rds"))
spots <- as.data.frame(SummarizedExperiment::colData(spe), optional = TRUE)
stopifnot(identical(as.character(spots$spot_id), colnames(spe)), nrow(spe) == 18082L,
          ncol(spe) == 4412L, !anyDuplicated(spots$spot_id))
lr_df <- read.csv(file.path(data_dir, "visium_benchmark_single_subunit_lr_candidates.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(lr_df) == 1154L, !anyDuplicated(lr_df$lr_id))
rownames(lr_df) <- lr_df$cellchat_id

# External total-count/detection filtering precedes the measured stage.
# The native maximum-count >=3 rule stays inside BLISA's measured call.
filter_start <- proc.time()
counts <- SummarizedExperiment::assay(spe, "counts")
eligible <- rownames(counts)[Matrix::rowSums(counts) >= 200 & Matrix::rowSums(counts > 0) >= 200]
lr_for_blisa <- getFromNamespace("getLRpairs", "blisa")(
  gene_panel = eligible, LR_df = lr_df, species = "human")
filter_seconds <- unname((proc.time() - filter_start)["elapsed"])

# Linux process RSS includes native allocations and threads, excluding separate
# workers. Imports, loading, shared LR preparation and exports are not timed.
process_memory_mib <- function(field) {
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]][2]) / 1024
}
set.seed(seed)
gc(full = TRUE)
baseline_rss_mib <- process_memory_mib("VmRSS")
rss_peak_reset <- tryCatch({ writeLines("5", "/proc/self/clear_refs"); TRUE },
                           warning = function(w) FALSE, error = function(e) FALSE)
timing <- system.time({
  result <- blisa::blisa(spe, platform = "visium", bin_size = 100, min_total_counts = 0,
    LR_df = lr_for_blisa, group = NULL, nsim = nperm, cpu_threads = ncpu,
    p_cutoff = 0.05, dmax = radius, min_ligand = 3, min_receptor = 3)
}, gcFirst = FALSE)
peak_rss_mib <- if (rss_peak_reset) process_memory_mib("VmHWM") else NA_real_

bin_index <- as.integer(result$bins$bin_id)
stopifnot(identical(sort(bin_index), seq_len(nrow(spots))))
spot_id <- spots$spot_id[bin_index]
lr <- result$LR_results
lr_ids <- paste(lr$ligand.symbol, lr$receptor.symbol, sep = "__")
stopifnot(!anyDuplicated(lr_ids), all(lr_ids %in% lr_df$lr_id))
summary <- data.frame(lr_id = lr_df$lr_id, returned = lr_df$lr_id %in% lr_ids,
                       n_available = NA_integer_, n_sig = NA_integer_)
summary$native_status <- ifelse(summary$returned, "tested",
  ifelse(lr_df$ligand.symbol %in% eligible & lr_df$receptor.symbol %in% eligible,
          "native_filtered", "custom_filtered"))
summary$native_reason <- ifelse(summary$returned, "Returned by native BLISA",
  ifelse(summary$native_status == "native_filtered", "maximum count <3",
          "total count <200 or non-empty spots <200"))

# Full local p-values/scores/quadrants are retained in one compressed CSV, in the
# format used for imported results. No native result RDS is written to res/.
order_ids <- sort(lr_df$lr_id)
current <- order_ids[order_ids %in% lr_ids]
path <- file.path(out, "blisa_permutation_local_lri_results_10000_samples.csv.gz")
con <- gzfile(path, open = "wt")
if (!length(current)) writeLines("lr_id,spot_id,p_value,score,tested,significant,quadrant", con)
for (j in seq_along(current)) {
  i <- match(current[j], lr_ids)
  p <- lr$all_pval[[i]]
  score <- lr$all_lisa[[i]]
  quadrant <- lr$all_quadrant[[i]]
  keep <- result$spatial_weights[[if (lr$ccc_mode[i] == "nearby") "keep_idx_queen" else "keep_idx_dist"]]
  tested <- seq_len(nrow(spots)) %in% keep & is.finite(p)
  significant <- seq_len(nrow(spots)) %in% lr$sig_index[[i]]
  stopifnot(all(tested[significant]), all(quadrant[significant] == "High-High"),
            all(p[significant] <= .05))
  p[!tested] <- score[!tested] <- NA_real_
  local <- data.frame(lr_id = lr_ids[i], spot_id = spot_id, p_value = p, score = score,
                       tested = as.integer(tested), significant = ifelse(tested, as.integer(significant), NA_integer_),
                       quadrant = as.character(quadrant))
  write.table(local, con, sep = ",", row.names = FALSE, col.names = j == 1L, na = "NA")
  summary$n_available[summary$lr_id == lr_ids[i]] <- sum(tested)
  summary$n_sig[summary$lr_id == lr_ids[i]] <- if (any(tested)) sum(significant) else NA_integer_
}
close(con)
summary$status <- ifelse(!summary$returned, "NR",
                         ifelse(summary$n_available > 0, "available", "NA"))
summary$budget <- nperm
summary$condition <- "blisa"
counts_path <- file.path(out, "visium_lri_hotspot_counts_and_test_status_10000_samples.csv")
all_counts <- read.csv(counts_path, stringsAsFactors = FALSE)
all_counts$status[is.na(all_counts$status)] <- "NA"
write_csv(rbind(all_counts[all_counts$condition != "blisa", ], summary[, names(all_counts)]), counts_path)
performance <- data.frame(method = "blisa", condition = "blisa", stage = "lr_analysis",
  package_version = as.character(utils::packageVersion("blisa")), ncpu = ncpu,
  nperm = nperm, nperm_requested = nperm, seed = seed, radius_um = radius,
  n_genes = nrow(spe), n_bins = ncol(spe), n_lr_pairs = nrow(lr), n_input_lr_pairs = nrow(lr_df),
  n_prefilter_lr_pairs = nrow(lr_for_blisa), min_gene_total = 200, min_gene_max = 3, min_gene_spots = 200,
  filter_prep_elapsed_seconds = filter_seconds, elapsed_seconds = unname(timing["elapsed"]),
  user_cpu_seconds = unname(timing["user.self"]), system_cpu_seconds = unname(timing["sys.self"]),
  baseline_rss_mib = baseline_rss_mib, peak_rss_mib = peak_rss_mib,
  extra_rss_mib = max(0, peak_rss_mib-baseline_rss_mib), rss_peak_reset = rss_peak_reset,
  rss_scope = "Stage process peak; excludes separate workers",
  timing_scope = "Native spatial preprocessing, maximum-count filtering and LR inference; excludes shared database preparation, input loading, external total-count/detection filtering and exports")
performance$budget <- nperm
performance$folder <- "blisa"
performance_path <- file.path(out, "visium_benchmark_computation_resources_10000_samples.csv")
all_performance <- read.csv(performance_path, stringsAsFactors = FALSE, check.names = FALSE)
for (name in setdiff(names(all_performance), names(performance))) performance[[name]] <- NA
write_csv(rbind(all_performance[!(all_performance$condition == "blisa" & all_performance$budget == nperm), ],
                performance[, names(all_performance)]), performance_path)
message("BLISA results updated; stLearn/SpatialDM saved files were unchanged.")
