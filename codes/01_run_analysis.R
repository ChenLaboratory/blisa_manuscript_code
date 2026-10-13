#!/usr/bin/env Rscript
# Run from CCC_manuscript_code/. Each step starts a separate R session.
# Usage: Rscript codes/01_run_analysis.R [all|xenium|visium|aligned_comparison]
# Benchmark import and optional BLISA reruns are separate explicit commands;
# this wrapper never runs stLearn or SpatialDM.

source("codes/helpers/source_data_helpers.R")
steps <- list(
  xenium = c("01_prepare_xenium_50um_bins.R", "02_run_xenium_lri_cci.R"),
  visium = c("01_prepare_visium_qc_annotations.R", "02_prepare_visium_lr_candidates.R",
             "03_run_visium_lri.R"),
  aligned_comparison = c("01_prepare_aligned_overlap.R", "02_prepare_shared_lr_pairs.R",
                         "03_run_aligned_lri.R")
)
args <- commandArgs(trailingOnly = TRUE)
workflow <- if (length(args)) args[1L] else "all"
if (length(args) > 1L || !workflow %in% c("all", names(steps))) {
  stop("Choose all, xenium, visium, or aligned_comparison.")
}
selected <- if (workflow == "all") names(steps) else workflow
Sys.setenv(CCC_REPO_ROOT = repo_root())
for (name in selected) {
  for (filename in steps[[name]]) {
    script <- repo_path("codes", "analysis", name, filename)
    message("Running ", name, "/", filename)
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(script)))
    if (status != 0L) stop("Stage failed: ", script)
  }
}
