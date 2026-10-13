#!/usr/bin/env Rscript
# Optional recovery from the original workspace, never permutation inference.
# Usage: Rscript codes/00_import_saved_results.R /path/to/original/workspace [all|visium|aligned_comparison|xenium_benchmark|visium_benchmark]
# Already staged CSVs are sufficient for plotting; this command is not required.
source("codes/helpers/source_data_helpers.R")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) stop("Supply the original workspace path and optionally one workflow.")
source_root <- normalizePath(args[1], mustWork = TRUE)
workflow <- if (length(args) == 2L) args[2] else "all"
scripts <- c(
  visium = "visium/04_import_saved_results.R",
  aligned_comparison = "aligned_comparison/04_import_saved_results.R",
  xenium_benchmark = "benchmarks/xenium/01_import_saved_xenium_benchmark_results.R",
  visium_benchmark = "benchmarks/visium/01_import_saved_visium_benchmark_results.R")
if (!workflow %in% c("all", names(scripts))) stop("Unknown import workflow: ", workflow)
selected <- if (workflow == "all") names(scripts) else workflow
Sys.setenv(CCC_REPO_ROOT = repo_root())
for (name in selected) {
  script <- repo_path("codes", "analysis", scripts[[name]])
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), "--source-root", shQuote(source_root)))
  if (status != 0L) stop("Import failed: ", name)
}
