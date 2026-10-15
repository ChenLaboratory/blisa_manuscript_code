#!/usr/bin/env Rscript
# Prepare the exact shared-panel LR input; no inference or downloads.
source("codes/helpers/source_data_helpers.R")
suppressPackageStartupMessages({library(SpatialExperiment); library(blisa)})

# 1. Read both overlap inputs and identify their shared gene panel.
data_dir <- repo_path("Data", "generated", "aligned_overlap")
visium <- readRDS(file.path(data_dir, "visium_aligned_overlap_qc.rds"))
xenium <- readRDS(file.path(data_dir, "xenium_aggregated_to_visium_overlap.rds"))
# Preserve Xenium feature order, matching the original preparation's gene index.
shared_genes <- intersect(rownames(xenium), rownames(visium))

# 2. Load the downloaded CellChat database and retain complete shared-panel pairs.
database_file <- repo_path("Data", "downloaded", "cellchat_github", "CellChatDB.human.rda")
database_env <- new.env()
load(database_file, envir = database_env)
cache <- getFromNamespace(".blisa_cache", "blisa")
cache[["CellChatDB.human"]] <- database_env$CellChatDB.human$interaction
database <- getFromNamespace(".load_cellchat_db", "blisa")("human")
# Retain all complete ligand/receptor definitions present in the shared panel.
# Unlike the cross-method benchmark, this is not a singleton-only database.
shared_lr <- getFromNamespace("getLRpairs", "blisa")(
  gene_panel = shared_genes, LR_df = database)

# 3. Save the LR candidates used identically by both platform analyses.
write_csv(data.frame(database_id = rownames(shared_lr), shared_lr,
                     row.names = NULL, check.names = FALSE),
  repo_path("Data", "generated", "ligand_receptor_tables", "aligned_shared_lr_pairs.csv"))
