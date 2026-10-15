# Run the recorded native Visium analysis; output every returned LR pair and spot as CSV.
# The original main analysis used nsim=9999 (distinct from benchmark nperm=10000).
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "Matrix",
                   "callr", "ps", "data.table"))
# 1. Load the annotated count input and expression-filtered CellChat candidates.
spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))

filtered_lr <- utils::read.csv(repo_path("Data", "generated", "ligand_receptor_tables",
  "visium_expression_filtered_lr_candidates.csv"), check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")

rownames(filtered_lr) <- filtered_lr$native_lr_id
filtered_lr$native_lr_id <- NULL
invisible(load_visium_cellchat())

# 2. Run the recorded native BLISA call.
# Preserve the source call, including the package's remaining default arguments.
# The source notebook did not set an explicit seed here; no new seed is imposed.
measurement <- measure_blisa_call(blisa::blisa(spe, LR_df = filtered_lr,
  min_ligand = 3, min_receptor = 3, nsim = 9999, cpu_threads = 4L))
res <- measurement$value

# 3. Save computation resources once, separately from the LR results.
performance <- data.frame(method = "BLISA", nsim = 9999L, cpu_threads = 4L,
  elapsed_seconds = measurement$performance$elapsed_seconds,
  peak_rss_mib = measurement$performance$sampled_peak_rss_mib,
  timing_scope = measurement$performance$timing_scope,
  memory_scope = measurement$performance$memory_scope)
write_csv(performance, repo_path("res", "visium", "visium_blisa_computation_resources.csv"))

# 4. Export native LR summaries and every spot-level test, including untested spots.
local <- visium_spot_table(res, spe)
summary <- lri_result_summary(res)
write_csv(summary, repo_path("res", "visium", "lri_summary.csv"))
write_visium_spot_csv(local, repo_path("res", "visium"))
