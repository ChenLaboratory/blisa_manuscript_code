# Figure 4d-f,h-i and Supplementary S7,S8b-c, using saved 10,000-sample results.
# The Python companion preserves the original matplotlib spatial/heatmap layout.
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "helpers", "plotting_helpers.R"))
source(repo_path("codes", "helpers", "benchmark_jaccard_helpers.R"))
source(repo_path("codes", "helpers", "visium_benchmark_plotting_helpers.R"))
require_packages(c("ggplot2", "scales"))
suppressPackageStartupMessages(library(ggplot2))
read_result <- function(name) read.csv(repo_path("res", "benchmarks", "visium", name),
                                      stringsAsFactors = FALSE, check.names = FALSE)
methods <- c("blisa", "stlearn_bh", "spatialdm_permutation")
labels <- setNames(c("blisa", "stLearn", "SpatialDM"), methods)
colours <- setNames(c("#222222", "#0072B2", "#D55E00"), methods)
shapes <- setNames(c(21, 24, 22), methods)
performance <- read_result("visium_benchmark_computation_resources_10000_samples.csv")
performance$peak_rss_gib <- performance$peak_rss_mib / 1024
performance$elapsed_minutes <- performance$elapsed_seconds / 60
performance$condition <- factor(performance$condition, levels = methods)
counts <- read_result("visium_lri_hotspot_counts_and_test_status_10000_samples.csv")
overlap <- read_result("visium_benchmark_jaccard_by_lr.csv")
pooled <- read_result("visium_benchmark_jaccard_pooled.csv")
# Python CSV exports use True/False rather than R's uppercase logical literals.
counts$returned <- tolower(as.character(counts$returned)) == "true"
overlap$common_three_eligible <- tolower(as.character(overlap$common_three_eligible)) == "true"
# 1. Select the saved common LR set in the native BLISA hotspot-count order.
blisa <- counts[counts$condition == "blisa", ]
blisa <- blisa[order(!blisa$returned, is.na(blisa$n_sig), -blisa$n_sig, blisa$lr_id), ]
common_ids <- blisa$lr_id[blisa$lr_id %in% overlap$lr_id[overlap$common_three_eligible]]
save_panel <- function(plot, name, width, height) {
  save_figure_panel(plot + labs(title = NULL, caption = NULL), name, width, height,
                    device = grDevices::pdf, useDingbats = FALSE)
}
# 2. Draw runtime/RAM, native hotspot counts and saved Jaccard scores.
save_panel(runtime_plot(), "figure4d", 3.35, 3.35)
save_panel(linked_plot(common_ids[1:30], NULL), "figure4e", 6.9, 4.4)
save_panel(linked_plot(common_ids, NULL, FALSE), "figureS7a", 6.9, 4.4)
p <- jaccard_summary(overlap, pooled, common_ids[1:10], methods, labels, colours,
                     NULL, sprintf("Pooled (%d common LRs)", length(common_ids)), NULL)
save_panel(p, "figure4f", 6.9, 5.04)
# 3. Draw spatial maps and heatmaps from saved local and regional results.
python <- Sys.getenv("CCC_PYTHON", unset = "python3")
status <- system2(python, c(shQuote(repo_path("codes", "plotting",
  "05_plot_visium_benchmark_Figure4d_f_h_i_FigureS7_S8b_c.py")),
  "--repo-root", shQuote(repo_root()),
  vapply(commandArgs(trailingOnly = TRUE), shQuote, character(1))))
if (status != 0L) stop("Visium benchmark spatial/heatmap plotting failed.")
