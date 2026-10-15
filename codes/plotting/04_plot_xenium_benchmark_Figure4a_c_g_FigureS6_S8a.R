# Draw the Xenium benchmark panels from saved native and comparison results.
# Set CCC_PYTHON to a Python executable with matplotlib, pandas, numpy and scikit-image.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "helpers", "benchmark_jaccard_helpers.R"))
source(repo_path("codes", "helpers", "xenium_benchmark_plot_helpers.R"))
source(repo_path("codes", "helpers", "manuscript_plotting_helpers.R"))
require_packages(c("ggplot2", "scales", "sf"))
library(ggplot2)
methods <- c("blisa", "stlearn", "spatialdm")
labels <- c(blisa = "blisa", stlearn = "stLearn", spatialdm = "SpatialDM")
colours <- c(blisa = "#222222", stlearn = "#0072B2", spatialdm = "#D55E00")
shapes <- c(blisa = 21, stlearn = 24, spatialdm = 22)
result_dir <- repo_path("res", "benchmarks", "xenium")
read_result <- function(name) read.csv(file.path(result_dir, paste0("xenium_benchmark_", name)),
  check.names = FALSE)
performance <- read_result("computation_resources.csv")
performance$elapsed_minutes <- performance$elapsed_seconds / 60
performance$peak_rss_gib <- performance$peak_rss_mib / 1024
dot_data <- read_result("lr_summary.csv")
overlap <- read_result("jaccard_by_lr.csv")
pooled <- read_result("jaccard_pooled.csv")
binned <- readRDS(repo_path("Data", "generated", "xenium", "xenium_50um_binned.rds"))
figures <- repo_path("Figures")
dir.create(figures, showWarnings = FALSE)
save_plot <- function(plot, panel, width, height) {
  ggsave(file.path(figures, paste0(panel, ".pdf")), plot, device = grDevices::pdf,
    width = width, height = height, units = "in", useDingbats = FALSE, bg = "white")
}

# 1. Use the common LR set already determined by the comparison analysis.
ids <- unique(dot_data$lr_id[order(dot_data$rank)])
common_ids <- ids[ids %in% overlap$lr_id[overlap$common_three_eligible]]
common <- dot_data[dot_data$lr_id %in% common_ids, ]

# 2. Draw runtime/RAM and native significant-bin counts with the manuscript layout.
scatter <- performance[match(methods, performance$method), ]
scatter$method <- factor(scatter$method, levels = methods)
memory_max <- max(scatter$peak_rss_gib)
runtime_max <- max(scatter$elapsed_minutes, 1)
memory_breaks <- pretty(c(0, memory_max * 1.25), n = 4)
runtime_breaks <- pretty(c(0, runtime_max * 1.5), n = 4)
runtime_plot <- ggplot(scatter,
  aes(x = peak_rss_gib, y = elapsed_minutes, colour = method,
    shape = method, fill = method)) +
  geom_point(size = 3.8, stroke = 0.8) +
  geom_text(aes(label = sprintf("%.2f min\n%.2f GiB", elapsed_minutes, peak_rss_gib),
    x = peak_rss_gib + ifelse(method == "spatialdm", -1, 1) * max(memory_breaks) * 0.055,
    hjust = ifelse(method == "spatialdm", 1, 0)),
    size = 2.6, lineheight = 0.95, vjust = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = colours, breaks = methods, labels = labels, drop = FALSE) +
  scale_shape_manual(values = shapes, breaks = methods, labels = labels, drop = FALSE) +
  scale_fill_manual(values = colours, breaks = methods, labels = labels, drop = FALSE) +
  scale_x_continuous(limits = range(memory_breaks),
    breaks = memory_breaks, labels = scales::label_comma(),
    expand = expansion(mult = c(0.035, 0.02))) +
  scale_y_continuous(limits = range(runtime_breaks),
    breaks = runtime_breaks, labels = scales::label_comma(),
    expand = expansion(mult = c(0.035, 0.02))) +
  labs(x = "Peak RAM (GiB)", y = "Runtime (min)", colour = NULL, shape = NULL, fill = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 11, base_family = "sans") +
  theme(aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(face = "bold", colour = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "bottom", legend.text = element_text(size = 10),
    legend.key.size = grid::unit(0.18, "in"), legend.spacing.x = grid::unit(0.05, "in"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0.04, "in"),
    plot.margin = margin(4, 4, 4, 4))
save_plot(runtime_plot, "figure4a", 3.35, 3.35)
save_plot(linked_dot(common, NULL), "figure4b", 6.9, 4.4)

# 3. Display the saved per-LR and pooled Jaccard scores.
overlap_plot <- jaccard_summary(overlap, pooled, common_ids, methods, labels, colours,
  NULL, sprintf("Pooled (%d common LRs)", length(common_ids)), NULL) + labs(caption = NULL)
save_plot(overlap_plot, "figure4c", 6.9, 5.52)

# 4. Pass saved polygons, annotations and native BLISA tests to Python.
# Temporary plotting tables are derived inputs, not additional repository results.
geometry_dir <- tempfile("xenium_benchmark_geometry_")
dir.create(geometry_dir)
metadata <- binned$bins_metadata
write_csv(metadata, file.path(geometry_dir, "bins.csv"))
vertices <- do.call(rbind, lapply(seq_len(nrow(binned$bins)), function(i) {
  xy <- sf::st_coordinates(sf::st_geometry(binned$bins)[i])
  data.frame(sample_id = metadata$sample_id[i], bin_id = metadata$bin_id[i],
    vertex = seq_len(nrow(xy)), x = xy[, "X"], y = xy[, "Y"])
}))
write_csv(vertices, file.path(geometry_dir, "hexagons.csv"))
write_csv(paperROIBorderPaths(binned$paper_rois), file.path(geometry_dir, "roi_paths.csv"))
result <- readRDS(file.path(result_dir, "xenium_benchmark_blisa_result.rds"))
lr <- result$LR_results
sample_ids <- metadata$sample_id[match(result$bins$bin_id, metadata$bin_id)]
local <- do.call(rbind, lapply(seq_len(nrow(lr)), function(i) {
  keep <- result$spatial_weights[[if (lr$ccc_mode[i] == "nearby")
    "keep_idx_queen" else "keep_idx_dist"]]
  available <- seq_len(nrow(result$bins)) %in% keep & is.finite(lr$all_pval[[i]]) &
    is.finite(lr$all_lisa[[i]])
  data.frame(method = "blisa", lr_id = paste(lr$ligand.symbol[i],
    lr$receptor.symbol[i], sep = "__"), sample_id = sample_ids,
    p_value = ifelse(available, sprintf("%.17g", lr$all_pval[[i]]), NA_character_),
    local_available = available,
    significant = seq_len(nrow(result$bins)) %in% lr$sig_index[[i]])
}))
write_csv(local, file.path(geometry_dir, "blisa_local.csv"))
python <- Sys.getenv("CCC_PYTHON", unset = Sys.which("python3"))
status <- system2(python, c(shQuote(repo_path("codes", "plotting",
  "04_plot_xenium_benchmark_spatial_Figure4g_FigureS6_S8a.py")),
  "--repo", shQuote(repo_root()), "--geometry", shQuote(geometry_dir)))
if (status != 0L) stop("The Xenium benchmark spatial plotting script failed.")
