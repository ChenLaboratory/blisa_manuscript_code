# Shared paths and file I/O. Source from the repository root.
# This file does not download data or run analysis.

repo_root <- function() {
  normalizePath(Sys.getenv("CCC_REPO_ROOT", unset = "."), mustWork = TRUE)
}

repo_path <- function(...) file.path(repo_root(), ...)

require_packages <- function(packages) {
  invisible(lapply(packages, loadNamespace))
}

# Release objects remain outside the repository. Supply individual paths or
# a directory containing the two release filenames listed in the manifest.
source_data_path <- function(id) {
  manifest <- utils::read.csv(repo_path("Data", "data_sources.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  entry <- manifest[manifest$dataset_id == paste0(id, "_release"), , drop = FALSE]
  filename <- basename(entry$source_url)
  explicit <- Sys.getenv(paste0("CCC_", toupper(id), "_RDS"), unset = "")
  if (nzchar(explicit)) return(explicit)
  data_dir <- Sys.getenv("CCC_SOURCE_DATA_DIR", unset = "")
  if (nzchar(data_dir)) return(file.path(data_dir, filename))
  stop("Set CCC_", toupper(id), "_RDS or CCC_SOURCE_DATA_DIR for ", entry$source_url)
}

read_source_data <- function(id) readRDS(source_data_path(id))

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA", quote = TRUE)
}
