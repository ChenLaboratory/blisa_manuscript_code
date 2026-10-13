# Import completed benchmark runs; never run BLISA, stLearn or SpatialDM.
# Rscript codes/analysis/benchmarks/xenium/01_import_saved_xenium_benchmark_results.R --source-root /path/to/original/workspace
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "xenium", "xenium_helpers.R"))
require_packages(c("blisa", "sf", "readxl", "SpatialExperiment"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != "--source-root") {
  stop("Supply --source-root followed by the original workspace directory.")
}
source_root <- normalizePath(args[2], mustWork = TRUE)
source_path <- function(...) file.path(source_root, "bechmark_methods_on_xenium", ...)
tables <- source_path("manuscript_update_20261012", "report_tables")
read_table <- function(name) read.csv(file.path(tables, name), check.names = FALSE)
read_native <- function(run, method, name) read.csv(source_path("results", run, method, name),
  check.names = FALSE)
flag <- function(x) as.character(x) %in% c("TRUE", "True", "1")
methods <- c("blisa", "stlearn", "spatialdm")

# The cache spans every returned LR identity, not just the figure selections.
local <- read_table("selected_lr_spot_results.csv.gz")
local$local_available <- flag(local$local_available)
local$significant <- flag(local$significant)
summary <- read_table("linked_dot_counts.csv")
summary$returned <- flag(summary$returned)
mapping <- read.csv(source_path("results", "29983663", "report_tables", "lr_summary.csv"),
  check.names = FALSE)
mapping <- unique(mapping[mapping$method %in% methods,
  c("method", "lr_id", "native_lr_id", "ligand", "receptor")])
key <- function(x) paste(x$method, x$lr_id, sep = "::")
stopifnot(!anyDuplicated(key(mapping)), !anyDuplicated(key(summary)))
mi <- match(key(local), key(mapping))
stopifnot(!anyNA(mi))
local$native_lr_id <- mapping$native_lr_id[mi]
si <- match(key(summary), key(mapping))
summary$native_lr_id <- mapping$native_lr_id[si]
summary$ligand <- mapping$ligand[si]
summary$receptor <- mapping$receptor[si]

binned <- readRDS(source_path("input", "binned", "binned.rds"))
metadata <- read.csv(source_path("input", "binned", "metadata.csv"), check.names = FALSE)
regions <- read.csv(source_path("input", "binned", "bin_regions.csv"), check.names = FALSE)
bin_index <- match(metadata$bin_id, binned$bins$bin_id)
stopifnot(nrow(metadata) == 17215L, !anyNA(bin_index),
  identical(metadata$sample_id, colnames(binned$counts_matrix)),
  identical(metadata$bin_id, regions$bin_id))
binned$bins <- binned$bins[bin_index, ]
metadata$region <- regions$region
binned$bins_metadata <- metadata
binned$paper_rois <- loadPaperROIs()$rois
stopifnot(identical(as.character(classifyPaperRegions(metadata[, c("x", "y")],
  binned$paper_rois)), as.character(metadata$region)))
input_path <- repo_path("Data", "generated", "benchmarks", "xenium", "xenium_benchmark_50um_binned.rds")
dir.create(dirname(input_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(binned, input_path)

# Restore p-values directly from completed native outputs. In particular,
# values immediately below 0.05 must not be rounded across the stLearn BH cutoff.
native_blisa <- readRDS(source_path("results", "29957811", "blisa", "result.rds"))
stopifnot(identical(native_blisa$bins$bin_id, binned$bins$bin_id))
# Frozen native database annotations of the 12 returned, gene-matched inputs.
# The optional benchmark run retains default database loading inside its timing.
lr_columns <- setdiff(names(native_blisa$LR_results), c("ccc_mode", "sig_numbers",
  "sig_index", "sig_pval", "all_pval", "all_lisa", "all_quadrant"))
write_csv(native_blisa$LR_results[, lr_columns], file.path(dirname(input_path),
  "xenium_benchmark_blisa_lr_input.csv"))
stlearn_p <- read_native("29957811", "stlearn", "p_adjs.csv")
stlearn_score <- read_native("29957811", "stlearn", "lr_scores.csv")
spatialdm_p <- read_native("29983663", "spatialdm_permutation", "local_perm_p.csv")
spatialdm_score <- read_native("29983663", "spatialdm_permutation", "local_I.csv")
spatialdm_global <- read_native("29983663", "spatialdm_permutation", "global_res.csv")
local$score <- NA_real_
local$quadrant <- NA_character_
local$bin_id <- metadata$bin_id[match(local$sample_id, metadata$sample_id)]
for (group in split(seq_len(nrow(local)), key(local))) {
  method <- local$method[group[1]]
  id <- local$native_lr_id[group[1]]
  ids <- local$sample_id[group]
  available <- local$local_available[group]
  if (method == "blisa") {
    i <- match(id, rownames(native_blisa$LR_results))
    j <- match(ids, metadata$sample_id)
    p <- native_blisa$LR_results$all_pval[[i]][j]
    score <- native_blisa$LR_results$all_lisa[[i]][j]
    local$quadrant[group] <- native_blisa$LR_results$all_quadrant[[i]][j]
    expected <- j %in% native_blisa$LR_results$sig_index[[i]]
  } else if (method == "stlearn") {
    p <- stlearn_p[[id]][match(ids, stlearn_p$sample_id)]
    score <- stlearn_score[[id]][match(ids, stlearn_score$sample_id)]
    expected <- available & score > 0 & p < 0.05
  } else {
    p <- spatialdm_p[[id]][match(ids, spatialdm_p$sample_id)]
    score <- spatialdm_score[[id]][match(ids, spatialdm_score$sample_id)]
    gate <- spatialdm_global$fdr[match(id, spatialdm_global$cellchat_id)] < 0.10
    expected <- available & gate & p < 0.05
  }
  stopifnot(length(p) == length(group),
    all(abs(p[available] - local$p_value[group][available]) < 1e-13),
    identical(as.logical(expected), local$significant[group]))
  local$p_value[group] <- ifelse(available, p, NA_real_)
  local$score[group] <- score
}
summary$global_permutation_p <- summary$global_fdr <- NA_real_
sd <- which(summary$method == "spatialdm")
index <- match(summary$native_lr_id[sd], spatialdm_global$cellchat_id)
summary$global_permutation_p[sd] <- spatialdm_global$perm_pval[index]
summary$global_fdr[sd] <- spatialdm_global$fdr[index]
rules <- read_table("comparison_rules.csv")
summary$significance_rule <- rules$rule[match(summary$method, rules$method)]
for (i in seq_len(nrow(summary))) {
  rows <- local[key(local) == key(summary[i, ]), ]
  count <- if (any(rows$local_available)) sum(rows$significant) else NA_integer_
  stopifnot(isTRUE(all.equal(as.numeric(count), summary$significant_bins[i])))
}

performance <- read_table("runtime_memory_data.csv")
stopifnot(setequal(performance$method, methods), nrow(performance) == 3L)
# Retain run provenance without embedding workspace-specific absolute paths.
for (name in intersect(c("source_dir", "input_dir", "run_dir"), names(performance))) {
  performance[[name]] <- sub(paste0(source_root, "/"), "", performance[[name]], fixed = TRUE)
}
output <- repo_path("res", "benchmarks", "xenium")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_csv(summary, file.path(output, "xenium_benchmark_lr_summary.csv"))
write_csv(performance, file.path(output, "xenium_benchmark_computation_resources.csv"))
# Base write.csv formats numeric doubles to 15 digits. Format these two columns
# explicitly so reading the CSV round-trips native scores and threshold p-values.
local$p_value <- ifelse(is.na(local$p_value), NA_character_, sprintf("%.17g", local$p_value))
local$score <- ifelse(is.na(local$score), NA_character_, sprintf("%.17g", local$score))
connection <- gzfile(file.path(output, "xenium_benchmark_lri_by_bin.csv.gz"), "wt")
write.csv(local, connection, row.names = FALSE, na = "NA")
close(connection)
cat("Imported completed Xenium benchmark results; no method inference was run.\n")
