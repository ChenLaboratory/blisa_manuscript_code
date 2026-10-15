#!/usr/bin/env python3
"""Run the Visium SpatialDM permutation benchmark (SpatialDM 0.2.0).

Run preparation steps 01a, 01b and 01c first, then:
    NCPU=40 python codes/analysis/benchmarks/visium/04_run_spatialdm_permutation_visium_benchmark.py

Extracted from benchmark_methods_on_visium/run_spatialdm_permutation.py.
Uses all staged genes, the frozen CellChat candidates and 10,000 permutations.
Global BH<0.1 gates local tests. Native sig_spots uses p<0.1; the manuscript
comparison exports the historical raw local permutation p<0.05 selection.
Only SpatialDM's results and summary/resource rows are replaced when executed.
"""
import os
import random

# Match the historical job wrapper before importing numerical libraries.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = os.environ.get("NCPU", "40")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

from visium_native_workflow_helpers import (
    NPERM, RADIUS_UM, SEED, export_local_results, export_performance,
    finish_measurement, load_inputs, start_measurement, spatialdm_local_frames,
)


import scanpy as sc
import spatialdm as sdm

# 1. Use the shared single-subunit pairs and native min_cell=3 prevalence rule.
# NCPU=40 reproduces the recorded allocation and divides 10,000 permutations.
adata, candidates, output, ncpu, _ = load_inputs()
candidates["native_lr_id"] = candidates.cellchat_id
ordered = candidates.copy()
ordered["source_annotation"] = ordered.annotation
ordered["annotation"] = ordered.annotation.replace({"Non-protein Signaling": "Secreted Signaling"})
ordered = ordered.sort_values("annotation", kind="stable")
# Native extract_lr min_cell=3 applied to raw counts, without a database
# download or BLISA's separate total/max/detected-spot gene filters.
expressed_spots = pd.Series(np.asarray((adata.X > 0).sum(axis=0)).ravel(), index=adata.var_names)
ordered = ordered.loc[(ordered.ligand.map(expressed_spots) >= 3)
                      & (ordered.receptor.map(expressed_spots) >= 3)]
adata.uns["ligand"] = ordered[["ligand"]].rename(columns={"ligand": "Ligand0"})
adata.uns["receptor"] = ordered[["receptor"]].rename(columns={"receptor": "Receptor0"})
adata.uns["geneInter"] = ordered.drop(columns=["ligand", "receptor"])
adata.uns["num_pairs"], adata.uns["mean"] = len(ordered), "algebra"
np.random.seed(SEED)
random.seed(SEED)

# 2. Log-transform counts and construct the native contact/diffuse weights.
# Time dense conversion, log transformation, weights and both inference stages.
measurement = start_measurement()
adata.X = adata.X.toarray().astype(float)
adata.raw = adata.copy()
sc.pp.log1p(adata)
sdm.weight_matrix(adata, l=100, cutoff=np.exp(-RADIUS_UM**2 / (2 * 100**2)),
                  n_nearest_neighbors=7, single_cell=False)
# 3. Select globally supported pairs (BH<0.1), then run local permutations.
sdm.spatialdm_global(adata, n_perm=NPERM, method="permutation", nproc=ncpu)
sdm.sig_pairs(adata, method="permutation", fdr=True, threshold=0.1)
print(f"Global LR pairs tested: {len(adata.uns['global_res'])}; "
      f"selected for local testing: {int(adata.uns['global_res'].selected.sum())}", flush=True)
sdm.spatialdm_local(adata, n_perm=NPERM, method="permutation", specified_ind=None, nproc=ncpu)
sdm.sig_spots(adata, method="permutation", fdr=False, threshold=0.1)
measured = finish_measurement(measurement)

# 4. Export local p<0.05 results and retain the global gate in pair status.
global_results = adata.uns["global_res"]
local_ids = global_results.index[global_results.selected]
status = pd.DataFrame({"returned": False, "status": "expression_filtered"}, index=candidates.index)
status.loc[global_results.index, ["returned", "status"]] = [True, "global_filtered"]
for pair in local_ids:
    tested = np.isfinite(adata.uns["local_perm_p"].loc[pair, adata.obs_names].to_numpy())
    status.loc[pair, "status"] = "tested" if tested.any() else "local_untestable"
export_local_results(output, "spatialdm_permutation", "spatialdm_permutation", "p_value",
                     candidates, adata.obs_names, spatialdm_local_frames(adata, candidates), status)
# 5. Save this method's runtime, process peak RAM and analysis settings.
performance = dict(method="spatialdm", nperm=NPERM, ncpu=ncpu,
    timing_scope="Dense count conversion, raw copy, log1p, spatial weights and global/local inference; excludes shared LR preparation/loading, min_cell=3 filtering, native table formatting and exports",
    **measured)
export_performance(output, "spatialdm_permutation", performance)
print("Updated SpatialDM permutation results; other saved methods unchanged.", flush=True)
