# Shared figure export. Domain-specific plotting functions stay with their figures.
# Sourcing this helper does not create a graphics device or write a file.
if (!interactive()) {
  # Heatmap layout may request a temporary device before the explicit PDF export.
  # Keep that device in memory instead of leaving Rplots.pdf in the repository.
  options(device = function(...) grDevices::pdf(file = NULL, ...))
}
save_figure_panel <- function(plot, panel, width, height,
                              device = grDevices::cairo_pdf, ...) {
  directory <- repo_path("Figures")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(directory, paste0(panel, ".pdf")), plot = plot,
    width = width, height = height, device = device, bg = "white", ...)
}
