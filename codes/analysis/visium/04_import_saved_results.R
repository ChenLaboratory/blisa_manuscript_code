#!/usr/bin/env Rscript
# Import the original saved Visium inference without running permutations.
# Usage from repository root: Rscript codes/analysis/visium/04_import_saved_results.R --source-root /path/to/original/workspace
source("codes/helpers/source_data_helpers.R")
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "Matrix", "data.table", "digest"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1L] != "--source-root") {
  stop("Supply --source-root followed by the original manuscript workspace.")
}
source_root <- normalizePath(args[2L], mustWork = TRUE)
source_dir <- file.path(source_root, "blisa_on_visium")
input <- repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds")
if (!file.exists(input)) {
  Sys.setenv(CCC_VISIUM_RDS = file.path(source_dir, "input", "spe_visium_bc_s1.rds"))
  source(repo_path("codes", "analysis", "visium", "01_prepare_visium_qc_annotations.R"))
}
spe <- readRDS(input)
raw <- readRDS(file.path(source_dir, "input", "spe_visium_bc_s1.rds"))
raw <- raw[, Matrix::colSums(SummarizedExperiment::assay(raw, "counts")) > 3000]
stopifnot(identical(dim(spe), c(18085L, 4412L)),
  identical(dimnames(spe), dimnames(raw)),
  identical(SummarizedExperiment::assay(spe, "counts"), SummarizedExperiment::assay(raw, "counts")),
  identical(SpatialExperiment::spatialCoords(spe), SpatialExperiment::spatialCoords(raw)))
result_file <- file.path(source_dir, "lri_results", "result.rds")
result <- readRDS(result_file)
performance <- utils::read.csv(file.path(source_dir, "lri_results", "blisa_performance.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(result$bins) == 4412L, nrow(result$LR_results) == 1616L,
  performance$nsim == 9999L, performance$cpu_threads == 4L,
  performance$database_sha256 == "582e99db2bb1b1e05c0364553893d804952bd007fef870720a9cc659a7bd3ad8")
summary <- lri_result_summary(result)
for (field in names(performance)) summary[[paste0("run_", field)]] <- performance[[field]]
summary$import_source <- "blisa_on_visium/lri_results/result.rds"
summary$import_source_sha256 <- digest::digest(file = result_file, algo = "sha256")
tables <- native_report_tables(result, spe)
stopifnot(nrow(tables$lr_spot) == 1616L * 4412L,
  sum(tables$lr_spot$native_significant) == sum(summary$sig_numbers))
out <- repo_path("res", "visium")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv(summary, file.path(out, "lri_summary.csv"))
write_visium_spot_csv(tables$lr_spot, out)
message("Imported 1,616 LR pairs x 4,412 Visium spots; no inference performed.")
