# Import completed 10,000-sample benchmark outputs; no inference is rerun.
# Rscript codes/analysis/benchmarks/visium/01_import_saved_visium_benchmark_results.R \
#   --source-root /path/to/original/manuscript/workspace
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
python <- Sys.getenv("CCC_PYTHON", unset = "python3")
script <- repo_path("codes", "analysis", "benchmarks", "visium",
                    "01_import_saved_visium_benchmark_results.py")
status <- system2(python, c(shQuote(script), "--repo-root", shQuote(repo_root()),
                           vapply(commandArgs(trailingOnly = TRUE), shQuote, character(1))))
if (status != 0L) stop("Saved Visium benchmark import failed.")
