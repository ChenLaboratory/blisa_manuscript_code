#!/usr/bin/env Rscript
# Read processed inputs and CSV results; never run permutation inference.
# Usage: Rscript codes/02_make_figures.R [all|xenium|visium|aligned_comparison|xenium_benchmark|visium_benchmark]

source("codes/helpers/source_data_helpers.R")
scripts <- c(
  xenium = "01_plot_xenium_results_Figure2_FigureS1_S3.R",
  visium = "02_plot_visium_results_Figure3a_e_FigureS4.R",
  aligned_comparison = "03_plot_aligned_results_Figure3f_h_FigureS5.R",
  xenium_benchmark = "04_plot_xenium_benchmark_Figure4a_c_g_FigureS6_S8a.R",
  visium_benchmark = "05_plot_visium_benchmark_Figure4d_f_h_i_FigureS7_S8b_c.R"
)
args <- commandArgs(trailingOnly = TRUE)
workflow <- if (length(args)) args[1L] else "all"
if (length(args) > 1L || !workflow %in% c("all", names(scripts))) {
  stop("Choose all, xenium, visium, aligned_comparison, xenium_benchmark, or visium_benchmark.")
}
selected <- if (workflow == "all") names(scripts) else workflow
Sys.setenv(CCC_REPO_ROOT = repo_root())
for (name in selected) {
  script <- repo_path("codes", "plotting", scripts[[name]])
  message("Plotting ", name)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(script)))
  if (status != 0L) stop("Plotting failed: ", script)
}
