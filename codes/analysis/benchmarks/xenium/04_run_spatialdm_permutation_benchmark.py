#!/usr/bin/env python
"""Run SpatialDM permutation tests on prepared Xenium bins.

Run the main Xenium preparation and benchmark 01a/01b first, then execute
in the SpatialDM 0.2.0 environment.
Source: bechmark_methods_on_xenium/spatialdm_modes/run_spatialdm_permutation.py,
saved run 29983663. The original run set no random seed. Defaults: 10,000
permutations, 10 CPUs. Native extract_lr downloads its human Figshare tables
inside the timed stage, as in the original run; frozen GitHub candidate decoding
and result exports are outside it. All original 313 expression genes are kept.
"""
import argparse
import os

# 1. Set the recorded permutation budget and thread limits before imports.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ncpu", type=int, default=int(os.environ.get("NCPU", "10")))
parser.add_argument("--nperm", type=int, default=int(os.environ.get("NPERM", "10000")))
args = parser.parse_args()
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[name] = str(args.ncpu)
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdm as sdm
from xenium_benchmark_helpers import (INPUT,
    cellchat_candidates, start_measurement, finish_measurement, result_tables, update_outputs)

# 2. Load raw bins and the frozen panel-present CellChat singleton candidates.
adata = anndata.read_h5ad(INPUT / "xenium_benchmark_raw_counts.h5ad")
candidates = cellchat_candidates(adata.var_names)

# 3. Time log transformation, spatial weights and native LR preparation.
started = start_measurement()
adata.raw = adata.copy()
sc.pp.log1p(adata)
# n_nearest_neighbors remains SpatialDM 0.2.0's native default (6).
sdm.weight_matrix(adata, l=100, cutoff=np.exp(-(250**2) / (2 * 100**2)),
    n_neighbors=None, single_cell=False)
sdm.extract_lr(adata, "human", min_cell=3)

# 4. Align native pairs to the shared candidates and retain min-cell eligibility.
native_pairs = {(tuple(adata.uns["ligand"].loc[pair_id].dropna()),
    tuple(adata.uns["receptor"].loc[pair_id].dropna()))
    for pair_id in adata.uns["geneInter"].index}
retained = []
for pair_id, pair in candidates.iterrows():
    if (((pair.ligand,), (pair.receptor,)) in native_pairs or
        (int((adata[:, [pair.ligand]].X > 0).sum()) >= 3 and
         int((adata[:, [pair.receptor]].X > 0).sum()) >= 3)):
        retained.append(pair_id)
# Contact/ECM rows must precede secreted rows for SpatialDM's positional
# choice of spatial weights. All three tables share this exact order.
pairs = candidates.loc[retained].sort_values("annotation")
adata.uns["ligand"] = pairs[["ligand"]].rename(columns={"ligand": "Ligand0"})
adata.uns["receptor"] = pairs[["receptor"]].rename(columns={"receptor": "Receptor0"})
adata.uns["geneInter"] = pairs.drop(columns=["ligand", "receptor"])
adata.uns["num_pairs"] = len(pairs)

# 5. Select pairs with global permutation BH < 0.1, then run local tests.
sdm.spatialdm_global(adata, n_perm=args.nperm, specified_ind=None,
    method="permutation", nproc=args.ncpu)
sdm.sig_pairs(adata, method="permutation", fdr=True, threshold=0.1)
global_results = adata.uns["global_res"]
local_ids = global_results.index[global_results.selected]
if len(local_ids):
    sdm.spatialdm_local(adata, n_perm=args.nperm, method="permutation",
        specified_ind=None, nproc=args.ncpu)
    # Retain this timed native call from the historical run. Manuscript
    # export below applies local p < 0.05, not this native 0.1 mask.
    sdm.sig_spots(adata, method="permutation", fdr=False, threshold=0.1)
measured = finish_measurement(started)

# 6. Export local permutation p < 0.05 calls after the global selection gate.
p_values = adata.uns["local_perm_p"].T if len(local_ids) else pd.DataFrame(index=adata.obs_names)
local, summary = result_tables("spatialdm", adata.obs, pairs, p_values,
    global_results=global_results)

# 7. Save complete results and the measured process-level computation resources.
performance = dict(method="spatialdm", condition="spatialdm_permutation",
    ncpu=args.ncpu, nperm=args.nperm,
    timing_scope="Raw copy, log1p, spatial weights, native database download/preparation and global/local inference; excludes input loading, frozen candidate decoding and exports",
    **measured)
update_outputs("spatialdm", local, summary, performance)
