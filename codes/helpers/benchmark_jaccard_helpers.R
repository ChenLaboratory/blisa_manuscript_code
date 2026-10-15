jaccard_summary <- function(overlap, pooled, lr_ids, methods, labels, colours,
                            title, pooled_label, caption) {
  # Three unique method comparisons; pooled values retain the saved denominator.
  pair_a <- methods[c(2, 3, 3)]
  pair_b <- methods[c(1, 1, 2)]
  ids <- c("Pooled", lr_ids)
  shown <- do.call(rbind, lapply(seq_along(ids), function(i) {
    source <- if (i == 1L) pooled else overlap[overlap$lr_id == ids[i], ]
    do.call(rbind, lapply(seq_along(pair_a), function(j) {
      value <- source[source$method_a == pair_a[j] & source$method_b == pair_b[j], ]
      status <- as.character(value$status)
      label <- if (is.finite(value$jaccard)) sprintf("%.2f", value$jaccard) else "NA"
      if (!is.na(status) && status %in% c("pair_not_returned", "NR")) label <- "NR"
      if (!is.na(status) && status %in% c("local_result_unavailable", "NA")) label <- "NA"
      if ((!is.na(status) && status == "no_significant_bins") ||
          (is.finite(value$union) && value$union == 0)) label <- "Empty"
      data.frame(lr_id = ids[i], method_a = pair_a[j], method_b = pair_b[j],
        intersection = value$intersection, union = value$union,
        jaccard = value$jaccard, status = status, label = label,
        x = j, y = length(ids) + 1L - i, stringsAsFactors = FALSE)
    }))
  }))
  n <- length(ids)
  row_labels <- c(pooled_label, gsub("__", "-", lr_ids, fixed = TRUE))
  key <- character()
  if ("NR" %in% shown$label) key <- c(key, "NR: pair not returned")
  if ("NA" %in% shown$label) key <- c(key, "NA: local result unavailable")
  if ("Empty" %in% shown$label) key <- c(key, "Empty: no significant bins in either method")
  if (length(key)) caption <- paste(c(caption, paste(key, collapse = "; ")), collapse = "\n")
  ggplot2::ggplot(shown, ggplot2::aes(x = x, y = y, fill = jaccard)) +
    ggplot2::geom_tile(colour = "#D9D9D9", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = label), fontface = "bold", size = 11 / 2.845,
      colour = ifelse(is.finite(shown$jaccard) & shown$jaccard > 0.65, "white", "#222222")) +
    ggplot2::annotate("rect", xmin = 0.5, xmax = 3.5, ymin = n - 0.5,
      ymax = n + 0.5, fill = NA, colour = "#333333", linewidth = 0.7) +
    ggplot2::annotate("text", x = 1:3, y = n + 2.1,
      label = unname(labels[pair_b]), colour = unname(colours[pair_b]),
      fontface = "bold", size = 12 / 2.845) +
    ggplot2::annotate("text", x = 1:3, y = n + 1.55,
      label = "vs", colour = "#555555", size = 8.5 / 2.845) +
    ggplot2::annotate("text", x = 1:3, y = n + 1,
      label = unname(labels[pair_a]), colour = unname(colours[pair_a]),
      fontface = "bold", size = 12 / 2.845) +
    ggplot2::scale_fill_gradientn(
      colours = grDevices::colorRampPalette(c("white", "#6B1740"), space = "rgb")(256),
      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1"), na.value = "#E6E6E6",
      name = "Jaccard similarity") +
    ggplot2::scale_x_continuous(breaks = NULL, limits = c(0.5, 3.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(n), labels = rev(row_labels),
      limits = c(0.5, n + 2.65), expand = c(0, 0)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, caption = caption) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(direction = "horizontal", raster = FALSE,
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(2.8, "in"), barheight = grid::unit(0.14, "in"))) +
    ggplot2::theme_minimal(base_size = 11, base_family = "sans") +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", size = 11, colour = "#222222",
        margin = ggplot2::margin(r = 8)),
      axis.text.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "#222222",
        margin = ggplot2::margin(b = 2)),
      plot.title.position = "plot", plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(size = 9, colour = "#444444", hjust = 0,
        lineheight = 1.15, margin = ggplot2::margin(t = 7)),
      legend.position = "bottom", legend.justification = "center",
      legend.title = ggplot2::element_text(face = "bold", size = 11),
      legend.text = ggplot2::element_text(face = "bold", size = 10),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(5, "pt"),
      plot.margin = ggplot2::margin(7, 8, 6, 5))
}
