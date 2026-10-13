#!/usr/bin/env Rscript
# Plot saved CSV results only. Does not run BLISA or change processed inputs.
# Execute from the repository root after importing saved results or running the analysis.
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "helpers", "plotting_helpers.R"))
suppressPackageStartupMessages({library(SpatialExperiment); library(ggplot2); library(sf)})
require_packages(c("ggrastr", "ragg", "ggrepel", "gridExtra", "scales"))

data_dir <- repo_path("Data", "generated", "aligned_overlap")
result_dir <- repo_path("res", "aligned_overlap")
v <- readRDS(file.path(data_dir, "visium_aligned_overlap_qc.rds"))
cell_file <- file.path(data_dir, "xenium_overlap_cells.csv")
if (!file.exists(cell_file)) {
  stop("Run codes/analysis/aligned_comparison/01_prepare_aligned_overlap.R first: ",
       "the aligned cell coordinates are required for Figure 3f.")
}
cells <- read.csv(cell_file, stringsAsFactors = FALSE)
read_result <- function(name) read.csv(file.path(result_dir, name),
  check.names = FALSE, stringsAsFactors = FALSE)
summary_v <- read_result("visium_lri_summary.csv")
summary_x <- read_result("xenium_lri_summary.csv")
spots_v <- read_result("visium_lri_by_spot.csv")
spots_x <- read_result("xenium_lri_by_spot.csv")
stopifnot(nrow(summary_v) == 15L, nrow(summary_x) == 15L,
          setequal(summary_v$database_id, summary_x$database_id),
          all(cells$assigned_spot_id %in% colnames(v)))

# Keep native Visium LR order; match every counterpart explicitly by LR and spot.
key_v <- paste(spots_v$database_id, spots_v$spot_id, sep = "::")
key_x <- paste(spots_x$database_id, spots_x$spot_id, sep = "::")
stopifnot(!anyDuplicated(key_v), !anyDuplicated(key_x), setequal(key_v, key_x),
          nrow(spots_v) == nrow(summary_v) * ncol(v))
spots_x <- spots_x[match(key_v, key_x), ]
stopifnot(identical(spots_v$spot_id, spots_x$spot_id),
          identical(spots_v$database_id, spots_x$database_id),
          identical(spots_v$tested, spots_x$tested))
summary_x_matched <- summary_x[match(summary_v$database_id, summary_x$database_id), ]
paired <- data.frame(database_id = summary_v$database_id,
  lr_label = gsub(" - ", "–", summary_v$lr_label, fixed = TRUE), mode = summary_v$mode,
  visium = summary_v$n_sig, xenium = summary_x_matched$n_sig)
for (i in seq_len(nrow(paired))) {
  keep <- spots_v$database_id == paired$database_id[i]
  stopifnot(sum(spots_v$significant[keep], na.rm = TRUE) == paired$visium[i],
            sum(spots_x$significant[keep], na.rm = TRUE) == paired$xenium[i])
}
common <- spots_v$tested & spots_x$tested
state <- rep("untested", nrow(spots_v))
state[common] <- ifelse(spots_v$significant[common],
  ifelse(spots_x$significant[common], "both", "visium_only"),
  ifelse(spots_x$significant[common], "xenium_only", "neither"))
spot_index <- match(spots_v$spot_id, colnames(v))
stopifnot(!anyNA(spot_index))
maps <- data.frame(database_id = spots_v$database_id,
  image_x = v$image_x[spot_index], image_y = -v$image_y[spot_index],
  state = factor(state, levels = c("neither", "visium_only", "xenium_only", "both", "untested")))

figure_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) + ggplot2::theme(
    text = ggplot2::element_text(colour = "black"),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(face = "bold", colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.7),
    panel.border = ggplot2::element_rect(linewidth = 0.9),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "grey92"),
    strip.text = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "white", linewidth = 0.8),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold"),
    legend.key.size = grid::unit(4, "mm"))
}

map_theme <- function(base_size = 11) {
  figure_theme(base_size) + ggplot2::theme(
    axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}

plot_lr_ranking <- function(summary, label, colour, unit) {
  # Retain native order for tied counts and rank each analysis independently.
  d <- summary[order(-summary$n_sig, seq_len(nrow(summary))), ]
  d$lr_label <- factor(d$lr_label, levels = rev(d$lr_label))
  ggplot2::ggplot(d, ggplot2::aes(n_sig, lr_label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = n_sig, yend = lr_label),
      colour = "grey80", linewidth = 0.5) +
    ggplot2::geom_point(colour = colour, size = 2.7) +
    ggplot2::geom_text(ggplot2::aes(label = scales::comma(n_sig, accuracy = 1)),
      hjust = -0.6, size = 3.1) +
    ggplot2::scale_x_continuous(limits = c(0, NA),
      labels = scales::label_comma(), breaks = scales::breaks_pretty(n = 4),
      expand = ggplot2::expansion(mult = c(0, 0.17))) +
    figure_theme(base_size = 11) +
    ggplot2::labs(title = label, subtitle = "blisa LR ranking", y = NULL,
      x = paste("Significant High–High", unit)) +
    ggplot2::theme(axis.text = ggplot2::element_text(face = "plain"),
      axis.text.y = ggplot2::element_text(size = 9.5),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.5),
      plot.title = ggplot2::element_text(size = 12))
}

save_aligned_panel <- save_figure_panel
# Figure 3f shows the saved aligned cells and spots without a border.
# The manuscript outline is added manually in BioRender.
overlap_cells <- data.frame(image_x = cells$image_x, image_y = -cells$image_y)
overlap_spots <- data.frame(image_x = v$image_x, image_y = -v$image_y)
p <- ggplot() +
  ggrastr::rasterise(geom_point(data = overlap_cells, aes(image_x, image_y),
    colour = "#E69F00", shape = 16, size = 0.12, alpha = 0.12),
    dpi = 450, dev = "ragg") +
  geom_point(data = overlap_spots, aes(image_x, image_y),
    colour = "#0072B2", shape = 1, size = 0.65, stroke = 0.2) +
  coord_sf(datum = NA) + map_theme(10) +
  theme(panel.border = element_blank()) +
  labs(title = "Aligned Xenium and Visium overlap",
    subtitle = paste0("Gold: ", format(nrow(overlap_cells), big.mark = ","),
      " Xenium cells; blue: ", format(nrow(overlap_spots), big.mark = ","),
      " Visium spots"))
save_aligned_panel(p, "figure3f", 5.4, 5.2)

scatter_theme <- figure_theme(10) + theme(
  axis.text = element_text(face = "plain", colour = "black"),
  axis.title = element_text(face = "plain"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_rect(linewidth = 0.6),
  axis.ticks = element_line(linewidth = 0.5), legend.position = "top",
  legend.title = element_text(face = "plain"),
  plot.margin = margin(5, 6, 5, 5))
lr_limit <- ceiling(max(paired$visium, paired$xenium) / 100) * 100
selected_lr <- c("CDH1_CDH1", "CXCL12_CXCR4", "CEACAM8_CEACAM6")
stopifnot(all(selected_lr %in% paired$database_id))
selected_labels <- paired
selected_labels$lr_label <- ifelse(selected_labels$database_id %in% selected_lr,
  gsub(" - ", "–", selected_labels$lr_label, fixed = TRUE), "")
p_selected <- ggplot(paired, aes(visium, xenium, colour = mode)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
    colour = "grey45", linewidth = 0.5) +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(data = selected_labels, aes(label = lr_label),
    size = 3.1, max.overlaps = Inf, seed = 1, box.padding = 0.4,
    point.padding = 0.3, min.segment.length = 0, segment.size = 0.3,
    max.time = 3, max.iter = 50000, show.legend = FALSE) +
  scale_colour_manual(values = c(nearby = "#0072B2", diffuse = "#D55E00"),
    labels = c(nearby = "Nearby", diffuse = "Diffuse")) +
  scale_x_continuous(limits = c(0, lr_limit), breaks = seq(0, lr_limit, 200),
    expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, lr_limit), breaks = seq(0, lr_limit, 200),
    expand = expansion(mult = 0)) +
  coord_equal() + scatter_theme +
  labs(x = "Visium significant spots", y = "Xenium significant spots",
    colour = NULL, subtitle = paste("Spearman rho =",
      round(stats::cor(paired$visium, paired$xenium, method = "spearman"), 3)))
save_aligned_panel(p_selected, "figure3g", 4, 3.8)

agreement_map <- function(ids, ncol = 3) {
  d <- maps[maps$database_id %in% ids, ]
  info <- paired[match(ids, paired$database_id), ]
  titles <- paste0(info$lr_label, "\nVisium: ", info$visium, "; Xenium: ", info$xenium)
  d$panel <- factor(d$database_id, levels = ids, labels = titles)
  ggplot(d, aes(image_x, image_y, colour = state)) +
    geom_point(size = 0.5) + coord_equal() + facet_wrap(~panel, ncol = ncol) +
    scale_colour_manual(values = c(neither = "grey85", visium_only = "#0072B2",
      xenium_only = "#E69F00", both = "#7B3294", untested = "grey35"), drop = FALSE,
      labels = c("Neither", "Visium only", "Xenium only", "Both", "Untested")) +
    map_theme() + theme(legend.position = "bottom",
      strip.text = element_text(size = 9, face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size = 3))) + labs(colour = NULL)
}
save_aligned_panel(agreement_map(selected_lr), "figure3h", 8, 3.2)

# S5a ranks each arm independently; S5b includes every one of the 15 shared pairs.
rankings <- gridExtra::arrangeGrob(
  plot_lr_ranking(summary_v, "Visium", "#0072B2", "spots"),
  plot_lr_ranking(summary_x, "Xenium at Visium spots", "#E69F00", "spots"), ncol = 2)
all_maps <- agreement_map(paired$database_id, ncol = 5)
save_aligned_panel(rankings, "figureS5a", 11.6, 4.8)
save_aligned_panel(all_maps, "figureS5b", 11.6, 9.6)
