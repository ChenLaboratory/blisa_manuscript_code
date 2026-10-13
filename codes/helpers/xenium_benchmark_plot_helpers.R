# Native linked-dot layout from the saved three-method Xenium report.
linked_dot <- function(data, title) {
  ids <- unique(data$lr_id[order(data$rank)])
  shown <- data
  shown$method <- factor(shown$method, levels = methods)
  shown$x <- match(shown$lr_id, ids)
  dot_max <- max(c(1, dot_data$significant_bins), na.rm = TRUE)
  dot_na_y <- -0.11 * dot_max
  dot_breaks <- pretty(c(0, dot_max), n = 5)
  dot_breaks <- dot_breaks[dot_breaks >= 0 & dot_breaks <= dot_max * 1.07]
  offsets <- c(blisa = -0.18, stlearn = 0, spatialdm = 0.18)
  available <- shown[is.finite(shown$significant_bins), , drop = FALSE]
  missing <- shown[!is.finite(shown$significant_bins), , drop = FALSE]
  has_missing <- nrow(missing) > 0L
  segments <- do.call(rbind, lapply(methods, function(method) {
    x <- shown[shown$method == method, , drop = FALSE]
    x <- x[order(x$x), , drop = FALSE]
    x$x_next <- c(x$x[-1], NA_real_)
    x$y_next <- c(x$significant_bins[-1], NA_real_)
    x[is.finite(x$significant_bins) & is.finite(x$y_next), , drop = FALSE]
  }))
  ggplot(available, aes(x = x, y = significant_bins, colour = method)) +
    (if (has_missing) annotate("rect", xmin = 0.45, xmax = length(ids) + 0.55,
      ymin = -0.17 * dot_max, ymax = -0.035 * dot_max, fill = "#F2F2F2", colour = NA)) +
    geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.35) +
    geom_segment(data = segments, aes(xend = x_next, yend = y_next, linetype = method), linewidth = 0.65) +
    geom_point(aes(shape = method, fill = method), size = 2.6, stroke = 0.75) +
    geom_point(data = missing,
      aes(x = x + offsets[as.character(method)], y = dot_na_y,
        colour = method, shape = method, fill = method), inherit.aes = FALSE,
      size = 2.1, stroke = 0.65, show.legend = FALSE) +
    scale_colour_manual(values = colours, breaks = methods, labels = labels, drop = FALSE) +
    scale_shape_manual(values = shapes, breaks = methods, labels = labels, drop = FALSE) +
    scale_fill_manual(values = colours, breaks = methods, labels = labels, drop = FALSE) +
    scale_linetype_manual(values = setNames(rep("solid", 3), methods),
      breaks = methods, labels = labels, drop = FALSE) +
    scale_x_continuous(breaks = seq_along(ids), labels = gsub("__", "-", ids, fixed = TRUE),
      limits = c(0.45, length(ids) + 0.55), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = if (has_missing) c(dot_na_y, dot_breaks) else dot_breaks,
      labels = if (has_missing) c("NA", scales::label_comma()(dot_breaks)) else scales::label_comma()(dot_breaks),
      limits = c((if (has_missing) -0.18 else -0.035) * dot_max, 1.07 * dot_max),
      expand = expansion(mult = 0)) +
    labs(title = title, x = NULL, y = "Significant bins", colour = NULL,
      shape = NULL, fill = NULL, linetype = NULL) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1),
      fill = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
    theme_bw(base_size = 11) +
    theme(text = element_text(face = "bold", colour = "#222222"),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = 0.3),
      panel.border = element_rect(colour = "#555555", linewidth = 0.5),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1,
        size = 11, face = "bold", colour = "#222222"),
      axis.text.y = element_text(size = 11, face = "bold", colour = "#222222"),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 5)),
      legend.position = "bottom", legend.justification = "left",
      legend.text = element_text(size = 11, face = "bold"),
      legend.key.width = grid::unit(0.8, "cm"), legend.key.height = grid::unit(0.4, "cm"),
      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"),
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
      plot.margin = margin(6, 8, 4, 5))
}
