# Summarize saved Visium LR tests; no BLISA permutations are rerun.
source("codes/helpers/source_data_helpers.R")
suppressPackageStartupMessages(library(SpatialExperiment))
out <- repo_path("res", "visium")

# 1. Read the native LR ranking and complete spot results in their saved order.
spe <- readRDS(repo_path("Data", "generated", "visium", "visium_qc_counts_gt3000_annotated.rds"))
summary <- read.csv(file.path(out, "lri_summary.csv"), check.names = FALSE)
summary <- summary[order(summary$rank), ]
local <- data.table::fread(file.path(out, "lri_by_spot.csv.gz"),
  select = c("lr_id", "spot_id", "p_value", "tested", "native_significant"),
  na.strings = "NA")
significant <- matrix(local$native_significant, nrow = ncol(spe))
p_values <- matrix(local$p_value, nrow = ncol(spe))

# 2. Apply the report's native minP pathway aggregation: each pathway uses the
# minimum P value among its LR High-High hotspots at that spot. No extra P-value
# adjustment is applied, matching blisaPathway()'s original defaults.
lr <- summary
rownames(lr) <- lr$lr_id
lr$sig_index <- I(lapply(seq_len(nrow(lr)), function(i) which(significant[, i])))
lr$sig_pval <- I(lapply(seq_len(nrow(lr)), function(i) p_values[significant[, i], i]))
pathway <- blisa::blisaPathway(lr)$LR_results
pathway_summary <- data.frame(rank = seq_len(nrow(pathway)),
  pathway = rownames(pathway), sig_numbers = pathway$sig_numbers,
  n_LR_pairs = pathway$n_LR_pairs)
write_csv(pathway_summary, file.path(out, "pathway_summary.csv"))
membership <- do.call(rbind, lapply(seq_len(nrow(pathway)), function(i) {
  data.frame(pathway = rownames(pathway)[i], lr_id = pathway$LR_pairs[[i]])
}))
write_csv(membership, file.path(out, "pathway_lr_pairs.csv"))
hotspots <- do.call(rbind, lapply(seq_len(nrow(pathway)), function(i) {
  data.frame(pathway = rep(rownames(pathway)[i], length(pathway$sig_index[[i]])),
    spot_id = colnames(spe)[pathway$sig_index[[i]]], p_value = pathway$sig_pval[[i]],
    top_lr_pairs = vapply(pathway$top_LR_pair[[i]], paste, collapse = ";", character(1)))
}))
data.table::fwrite(hotspots, file.path(out, "pathway_hotspots_by_spot.csv.gz"),
  compress = "gzip", na = "NA", quote = "auto")

# 3. Divide LR hotspot counts by all retained spots in each pathology region.
# The non-tumour region is the complement of the three registered tumour ROIs.
regions <- as.character(spe$region)
regions[regions == "Outside ROI"] <- "Non-tumour"
region_order <- c("Invasive_Tumor", "DCIS_1", "DCIS_2", "Non-tumour")
coverage <- do.call(rbind, lapply(region_order, function(region) {
  keep <- regions == region
  n_hotspots <- colSums(significant[keep, , drop = FALSE])
  data.frame(lr_id = summary$lr_id, region = region, n_hotspots = n_hotspots,
    n_spots = sum(keep), hotspot_proportion = n_hotspots / sum(keep))
}))
write_csv(coverage, file.path(out, "visium_lri_regional_coverage.csv"))

# 4. Count pathway membership and ligand/receptor sharing for the manuscript text.
# Complexes are single CellChatDB entities, and sharing is assessed in the same role.
tested <- summary[, c("lr_id", "pathway_name", "ligand", "receptor")]
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
