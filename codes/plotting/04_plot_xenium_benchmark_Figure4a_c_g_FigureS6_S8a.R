# Saved-result plots for the three-method Xenium benchmark; no inference.
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
dot_data <- read_result("lr_summary.csv")
local <- read_result("lri_by_bin.csv.gz")
binned <- readRDS(repo_path("Data", "generated", "benchmarks", "xenium", "xenium_benchmark_50um_binned.rds"))
figures <- repo_path("Figures")
dir.create(figures, showWarnings = FALSE)
save_plot <- function(plot, panel, width, height) {
  ggsave(file.path(figures, paste0(panel, ".pdf")), plot, device = grDevices::pdf,
    width = width, height = height, units = "in", useDingbats = FALSE, bg = "white")
}

# Eligibility means an available local result in every method, not significance.
ids <- unique(dot_data$lr_id[order(dot_data$rank)])
common_ids <- ids[vapply(ids, function(id) {
  rows <- dot_data[dot_data$lr_id == id, ]
  nrow(rows) == 3L && all(rows$returned & rows$n_available > 0)
}, logical(1))]
stopifnot(length(common_ids) == 12L)
common <- dot_data[dot_data$lr_id %in% common_ids, ]

# Match the manuscript report's measured runtime, pseudo-log axes and symbols.
scatter <- performance[match(methods, performance$method), ]
scatter$method <- factor(scatter$method, levels = methods)
memory_max <- max(scatter$peak_rss_gib)
runtime_max <- max(scatter$elapsed_minutes, 1)
runtime_plot <- ggplot(scatter,
  aes(x = peak_rss_gib, y = elapsed_minutes, colour = method,
    shape = method, fill = method)) +
  geom_point(size = 3.8, stroke = 0.8) +
  scale_colour_manual(values = colours, breaks = methods, labels = labels, drop = FALSE) +
  scale_shape_manual(values = shapes, breaks = methods, labels = labels, drop = FALSE) +
  scale_fill_manual(values = colours, breaks = methods, labels = labels, drop = FALSE) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5, base = 10),
    limits = c(0, memory_max * 1.25), breaks = c(0, 1, 10, 50),
    expand = expansion(mult = c(0.035, 0.02))) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.05, base = 10),
    limits = c(0, runtime_max * 1.5), breaks = c(0, 1, 5, 20),
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

# Keep LR-bin assignments distinct: pooled Jaccard sums intersections and unions.
sets <- split(local$sample_id[local$significant],
  paste(local$method[local$significant], local$lr_id[local$significant], sep = "::"))
overlap <- do.call(rbind, lapply(common_ids, function(id) {
  do.call(rbind, lapply(methods, function(a) {
    do.call(rbind, lapply(methods, function(b) {
      sa <- sets[[paste(a, id, sep = "::")]]
      sb <- sets[[paste(b, id, sep = "::")]]
      numerator <- length(intersect(sa, sb))
      denominator <- length(union(sa, sb))
      data.frame(lr_id = id, method_a = a, method_b = b,
        intersection = numerator, union = denominator,
        jaccard = if (denominator) numerator / denominator else NA_real_,
        status = if (denominator) "available" else "no_significant_bins")
    }))
  }))
}))
pooled <- aggregate(cbind(intersection, union) ~ method_a + method_b, overlap, sum)
pooled$jaccard <- ifelse(pooled$union > 0, pooled$intersection / pooled$union, NA_real_)
overlap_plot <- jaccard_summary(overlap, pooled, common_ids, methods, labels, colours,
  NULL, "Pooled (12 common LRs)", NULL) + labs(caption = NULL)
save_plot(overlap_plot, "figure4c", 6.9, 5.52)
print(pooled[pooled$method_a != pooled$method_b, ])

# Pass exact saved polygons and annotations to the original Python plotting style.
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
python <- Sys.getenv("CCC_PYTHON", unset = Sys.which("python3"))
if (!nzchar(python)) stop("Set CCC_PYTHON to a Python executable with the plotting dependencies.")
status <- system2(python, c(shQuote(repo_path("codes", "plotting",
  "04_plot_xenium_benchmark_spatial_Figure4g_FigureS6_S8a.py")),
  "--repo", shQuote(repo_root()), "--geometry", shQuote(geometry_dir)))
if (status != 0L) stop("The Xenium benchmark spatial plotting script failed.")
