# Summarize saved Visium LR tests; no BLISA permutations are rerun.
source("codes/helpers/source_data_helpers.R")
suppressPackageStartupMessages(library(SpatialExperiment))
out <- repo_path("res", "visium")

# 1. Read the native result and processed spot annotations.
spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))
result <- readRDS(file.path(out, "visium_blisa_result.rds"))
lr <- result$LR_results
# Native bins index tissue spots before the package's count filter.
spot_index <- which(spe$in_tissue == 1)[as.integer(result$bins$bin_id)]
significant <- vapply(lr$sig_index, function(index) {
  seq_len(ncol(spe)) %in% spot_index[index]
}, logical(ncol(spe)))

# 2. Apply the report's native minP pathway aggregation: each pathway uses the
# minimum P value among its LR High-High hotspots at that spot. No extra P-value
# adjustment is applied, matching blisaPathway()'s original defaults.
pathway <- blisa::blisaPathway(result)
saveRDS(pathway, file.path(out, "visium_blisa_pathway_result.rds"))

# 3. Divide LR hotspot counts by all retained spots in each pathology region.
# The non-tumour region is the complement of the three registered tumour ROIs.
regions <- as.character(spe$region)
regions[regions == "Outside ROI"] <- "Non-tumour"
region_order <- c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Non-tumour")
coverage <- do.call(rbind, lapply(region_order, function(region) {
  keep <- regions == region
  n_hotspots <- colSums(significant[keep, , drop = FALSE])
  data.frame(lr_id = rownames(lr), region = region, n_hotspots = n_hotspots,
    n_spots = sum(keep), hotspot_proportion = n_hotspots / sum(keep))
}))
write_csv(coverage, file.path(out, "visium_lri_regional_coverage.csv"))

# 4. Count pathway membership and ligand/receptor sharing for the manuscript text.
# Complexes are single CellChatDB entities, and sharing is assessed in the same role.
tested <- data.frame(lr_id = rownames(lr),
  lr[, c("pathway_name", "ligand", "receptor"), drop = FALSE], row.names = NULL)
pathways <- as.data.frame(table(tested$pathway_name), stringsAsFactors = FALSE)
names(pathways) <- c("pathway_name", "n_lr_pairs")
ligands <- table(tested$ligand)
receptors <- table(tested$receptor)
sharing <- data.frame(tested,
  n_pairs_with_ligand = as.integer(ligands[tested$ligand]),
  n_pairs_with_receptor = as.integer(receptors[tested$receptor]))
sharing$shares_ligand <- sharing$n_pairs_with_ligand > 1L
sharing$shares_receptor <- sharing$n_pairs_with_receptor > 1L
sharing$shares_either <- sharing$shares_ligand | sharing$shares_receptor
network <- data.frame(n_tested_lr_pairs = nrow(tested), n_pathways = nrow(pathways),
  median_lr_pairs_per_pathway = median(pathways$n_lr_pairs),
  n_pairs_sharing_ligand = sum(sharing$shares_ligand),
  n_pairs_sharing_receptor = sum(sharing$shares_receptor),
  n_pairs_sharing_either = sum(sharing$shares_either),
  percent_pairs_sharing_either = 100 * mean(sharing$shares_either))
write_csv(sharing, file.path(out, "visium_lri_partner_sharing.csv"))
write_csv(network, file.path(out, "visium_lri_network_summary.csv"))
