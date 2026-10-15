#!/usr/bin/env python
"""Run the stLearn BH benchmark on prepared Xenium bins.

Run the main Xenium preparation and benchmark 01a/01b first, then execute
in the stLearn 0.4.12 environment.
Source: bechmark_methods_on_xenium/run_stlearn.py, saved run 29957811.
Defaults reproduce the executed settings: all 15 singleton candidates in one
call, 10,000 background pairs, 10 CPUs. The original run set no random seed.
Normalization and inference are timed; input loading, CellChat decoding and
CSV export are excluded. Existing other-method results are retained.
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
import pandas as pd
import stlearn as st
from xenium_benchmark_helpers import (INPUT,
    cellchat_candidates, start_measurement, finish_measurement, result_tables, update_outputs)

# 2. Read raw counts and select CellChat singleton pairs present in the panel.
adata = anndata.read_h5ad(INPUT / "xenium_benchmark_raw_counts.h5ad")
adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]
adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]
candidates = cellchat_candidates(adata.var_names)
native_ids = candidates.ligand + "_" + candidates.receptor
lrs = native_ids.to_numpy(dtype=str)
# Keep the same baseline as the original runner: candidate metadata is
# released before measuring, then decoded again only for result export.
del candidates, native_ids

# 3. Time normalization and the native 250-um stLearn LR test.
started = start_measurement()
st.pp.normalize_total(adata)
st.tl.cci.run(adata, lrs, min_spots=20, distance=250,
    n_pairs=args.nperm, n_cpus=args.ncpu)
measured = finish_measurement(started)

# 4. Use BH-adjusted p-values and positive LR scores for local significance.
summary = adata.uns["lr_summary"]
candidates = cellchat_candidates(adata.var_names)
candidates.index = candidates.ligand + "_" + candidates.receptor
pairs = candidates.loc[summary.index]
p_values = pd.DataFrame(adata.obsm["p_adjs"], index=adata.obs_names,
    columns=summary.index)
scores = pd.DataFrame(adata.obsm["lr_scores"], index=adata.obs_names,
    columns=summary.index)
local, exported_summary = result_tables("stlearn", adata.obs, pairs, p_values, scores)

# 5. Save complete local results, hotspot counts and computation resources.
performance = dict(method="stlearn", condition="stlearn", ncpu=args.ncpu,
    nperm=args.nperm,
    timing_scope="Normalization and native LR inference; excludes input loading, CellChat decoding and exports",
    **measured)
update_outputs("stlearn", local, exported_summary, performance)
