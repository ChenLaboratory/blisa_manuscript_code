# Run the recorded native Visium analysis; output every returned LR pair and spot as CSV.
# The original main analysis used nsim=9999 (distinct from benchmark nperm=10000).
source(file.path(Sys.getenv("CCC_REPO_ROOT", unset = "."),
                 "codes", "helpers", "source_data_helpers.R"))
source(repo_path("codes", "analysis", "visium", "visium_helpers.R"))
require_packages(c("blisa", "SpatialExperiment", "SummarizedExperiment", "Matrix",
                   "digest", "callr", "ps", "data.table"))
spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))
stopifnot(nrow(spe) == 18085L, ncol(spe) == 4412L)
filtered_lr <- utils::read.csv(repo_path("Data", "generated", "ligand_receptor_tables",
  "visium_expression_filtered_lr_candidates.csv"), check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
stopifnot(nrow(filtered_lr) == 1917L, !anyDuplicated(filtered_lr$native_lr_id))
rownames(filtered_lr) <- filtered_lr$native_lr_id
filtered_lr$native_lr_id <- NULL
invisible(load_visium_cellchat())

# Preserve the source call, including the package's remaining default arguments.
# The source notebook did not set an explicit seed here; no new seed is imposed.
measurement <- measure_blisa_call(blisa::blisa(spe, LR_df = filtered_lr,
  min_ligand = 3, min_receptor = 3, nsim = 9999, cpu_threads = 4L))
res <- measurement$value
stopifnot(nrow(res$bins) == 4412L, nrow(res$LR_results) == 1616L)
performance <- measurement$performance
performance$blisa_version <- as.character(utils::packageVersion("blisa"))
performance$nsim <- 9999L
performance$cpu_threads <- 4L
performance$min_ligand <- performance$min_receptor <- 3
performance$min_total_count <- performance$min_nonempty_spots <- 200
performance$n_input_lr_pairs <- nrow(filtered_lr)
performance$database_sha256 <- "582e99db2bb1b1e05c0364553893d804952bd007fef870720a9cc659a7bd3ad8"
performance$database_source <- "https://github.com/jinworks/CellChat/raw/main/data/CellChatDB.human.rda"
performance$n_predicted_lr_pairs <- 1616L
performance$n_gene_features <- nrow(spe)
performance$n_spots <- nrow(res$bins)
performance$n_lr_pairs <- nrow(res$LR_results)
performance$database_state <- "Fixed native CellChat snapshot; external total/detection filter; native max-count-3 filter"
performance$timing_scope <- "blisa() only, including native max-count-3 filter; excludes preparation, exports and plotting"

tables <- native_report_tables(res, spe)
summary <- lri_result_summary(res)
# Repeat scalar run metadata so the public result CSVs are self-contained.
for (field in names(performance)) summary[[paste0("run_", field)]] <- performance[[field]]
write_csv(summary, repo_path("res", "visium", "lri_summary.csv"))
write_visium_spot_csv(tables$lr_spot, repo_path("res", "visium"))

# These summaries support the manuscript's pathway and partner-sharing text.
tested <- res$LR_results[, c("pathway_name", "ligand", "receptor"), drop = FALSE]
stopifnot(!anyNA(tested), all(vapply(tested,
  function(x) all(nzchar(trimws(x))), logical(1))))
pathways <- as.data.frame(table(tested$pathway_name), stringsAsFactors = FALSE)
names(pathways) <- c("pathway_name", "n_lr_pairs")
ligands <- table(tested$ligand)
receptors <- table(tested$receptor)
sharing <- data.frame(lr_id = rownames(tested), tested, row.names = NULL,
  n_pairs_with_ligand = as.integer(ligands[tested$ligand]),
  n_pairs_with_receptor = as.integer(receptors[tested$receptor]))
sharing$shares_ligand <- sharing$n_pairs_with_ligand > 1L
sharing$shares_receptor <- sharing$n_pairs_with_receptor > 1L
sharing$shares_either <- sharing$shares_ligand | sharing$shares_receptor
network_summary <- data.frame(n_tested_lr_pairs = nrow(tested), n_pathways = nrow(pathways),
  median_lr_pairs_per_pathway = median(pathways$n_lr_pairs),
  n_pairs_sharing_ligand = sum(sharing$shares_ligand),
  n_pairs_sharing_receptor = sum(sharing$shares_receptor),
  n_pairs_sharing_either = sum(sharing$shares_either),
  percent_pairs_sharing_either = 100 * mean(sharing$shares_either),
  sharing_definition = "Same CellChatDB entity in the same role in another tested pair; complexes are single entities")
print(network_summary)
