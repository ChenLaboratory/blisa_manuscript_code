# Shared paths and file I/O. Source from the repository root.
# This file does not download data or run analysis.

repo_root <- function() {
  root <- normalizePath(Sys.getenv("CCC_REPO_ROOT", unset = "."), mustWork = TRUE)
  if (!file.exists(file.path(root, "Data", "data_sources.csv"))) {
    stop("Run from CCC_manuscript_code/, or set CCC_REPO_ROOT to that directory.")
  }
  root
}

repo_path <- function(...) file.path(repo_root(), ...)

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the required packages first: ", paste(missing, collapse = ", "))
}

# Release objects remain outside the repository. Explicit file paths take
# precedence over a shared data directory and the original BLISA cache.
source_data_path <- function(id) {
  if (!id %in% c("xenium", "visium")) stop("Unknown release dataset: ", id)
  manifest <- utils::read.csv(repo_path("Data", "data_sources.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  entry <- manifest[manifest$dataset_id == paste0(id, "_release"), , drop = FALSE]
  stopifnot(nrow(entry) == 1L)
  filename <- basename(entry$source_url)
  explicit <- Sys.getenv(paste0("CCC_", toupper(id), "_RDS"), unset = "")
  if (nzchar(explicit)) {
    if (!file.exists(explicit)) stop("The configured ", id, " input does not exist: ", explicit)
    return(normalizePath(explicit))
  }
  data_dir <- Sys.getenv("CCC_SOURCE_DATA_DIR", unset = "")
  candidates <- c(if (nzchar(data_dir)) file.path(data_dir, filename),
                  file.path(tools::R_user_dir("blisa", "cache"), filename))
  available <- candidates[file.exists(candidates)]
  if (length(available)) return(normalizePath(available[1L]))
  stop("Missing external release data: ", filename,
       "\nDownload it from ", entry$source_url,
       "\nSet CCC_SOURCE_DATA_DIR or CCC_", toupper(id), "_RDS to its location.",
       "\nNo automatic download was attempted.")
}

read_source_data <- function(id) readRDS(source_data_path(id))

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA", quote = TRUE)
}
