#!/usr/bin/env python3
"""Run the Visium stLearn BH benchmark (stLearn 0.4.12).

Run preparation steps 01a, 01b and 01c first, then:
    NCPU=40 python codes/analysis/benchmarks/visium/03_run_stlearn_visium_benchmark.py

Extracted from benchmark_methods_on_visium/run_stlearn.py. Uses the frozen
CellChat candidates, all 18,082 input genes, 10,000 backgrounds and seed 20260915.
Replaces only stLearn's result file and its summary/resource rows when executed.
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
    CONTACT_RADIUS_UM, NPERM, RADIUS_UM, SEED, export_local_results,
    export_performance, finish_measurement, load_inputs, start_measurement, stlearn_local_frames,
)


import numba
import stlearn as st
from statsmodels.stats.multitest import multipletests

# 1. Divide CellChat pairs into contact (120 um) and diffuse (250 um) groups.
adata, candidates, output, ncpu, _ = load_inputs()
candidates["native_lr_id"] = candidates.ligand + "_" + candidates.receptor
candidates["signaling_group"] = np.where(
    candidates.annotation.isin(["Secreted Signaling", "Non-protein Signaling"]), "diffuse", "contact")
candidates["distance_um"] = np.where(candidates.signaling_group == "contact", CONTACT_RADIUS_UM, RADIUS_UM)
xy = np.asarray(adata.obsm["spatial"], float)
adata.obs["imagecol"], adata.obs["imagerow"] = xy[:, 0], xy[:, 1]

@numba.njit
def seed_numba(seed):
    np.random.seed(seed)

# 2. Normalize counts and run each signaling group with 10,000 backgrounds.
# Time normalization, both native calls and pooled BH; exclude input loading.
measurement = start_measurement()
st.pp.normalize_total(adata)
lr_ids, score_blocks, pvalue_blocks = [], [], []
for group, distance_um in (("contact", CONTACT_RADIUS_UM), ("diffuse", RADIUS_UM)):
    lrs = candidates.loc[candidates.signaling_group == group, "native_lr_id"].to_numpy(dtype=str)
    print(f"stLearn {group}: {len(lrs)} pairs, {distance_um:g} um, {NPERM} backgrounds", flush=True)
    group_adata = adata.copy()
    # stLearn 0.4.12 has no random_state parameter; reset all three samplers
    # for each group. Native backgrounds exclude that group's retained LR genes.
    np.random.seed(SEED)
    random.seed(SEED)
    seed_numba(SEED)
    st.tl.cci.run(group_adata, lrs, min_spots=20, distance=distance_um, n_pairs=NPERM,
                 n_cpus=ncpu, adj_method="fdr_bh", pval_adj_cutoff=0.05,
                 min_expr=0.0, save_bg=False, neg_binom=False)
    group_summary = group_adata.uns["lr_summary"]
    lr_ids.extend(group_summary.index)
    score_blocks.append(np.asarray(group_adata.obsm["lr_scores"]))
    pvalue_blocks.append(np.asarray(group_adata.obsm["p_vals"]))
    del group_adata

# 3. Apply BH across both groups within each spot, using positive-score tests.
scores = np.concatenate(score_blocks, axis=1)
pvalues = np.concatenate(pvalue_blocks, axis=1)
del score_blocks, pvalue_blocks
tested = (scores > 0) & np.isfinite(pvalues)
adjusted = np.ones(pvalues.shape)
for spot in range(adata.n_obs):
    eligible = tested[spot]
    if eligible.any():
        adjusted[spot, eligible] = multipletests(pvalues[spot, eligible], method="fdr_bh")[1]
significant = tested & (adjusted < 0.05)
summary = pd.DataFrame({"n_spots": (scores > 0).sum(axis=0),
    "n_spots_sig": significant.sum(axis=0),
    "n_spots_sig_pval": (tested & (pvalues < 0.05)).sum(axis=0)}, index=pd.Index(lr_ids))
pair_metadata = candidates.set_index("native_lr_id")
summary["signaling_group"] = pair_metadata.signaling_group.reindex(summary.index)
summary["distance_um"] = pair_metadata.distance_um.reindex(summary.index)
order = np.argsort(-summary.n_spots_sig.to_numpy(), kind="stable")
summary = summary.iloc[order]
adata.uns["lr_summary"] = summary
adata.obsm["lr_scores"] = scores[:, order]
adata.obsm["p_vals"] = pvalues[:, order]
adata.obsm["p_adjs"] = adjusted[:, order]
adata.obsm["lr_sig_scores"] = np.where(significant, scores, 0)[:, order]
measured = finish_measurement(measurement)

# 4. Export adjusted p-values and hotspot calls for every candidate and spot.
returned = candidates.native_lr_id.isin(summary.index)
status = pd.DataFrame({"returned": returned,
    "status": np.where(returned, "tested", "native_filtered")}, index=candidates.index)
export_local_results(output, "stlearn_bh", "stlearn_bh_adjusted", "bh_adjusted_p_value",
                     candidates, adata.obs_names, stlearn_local_frames(adata, candidates), status)
# 5. Save this method's runtime, process peak RAM and analysis settings.
performance = dict(method="stlearn", nperm=NPERM, ncpu=ncpu,
    timing_scope="Normalization, group copies, both native LR calls, combination and pooled BH; excludes input, shared LR preparation/loading, coordinate setup and exports",
    **measured)
export_performance(output, "stlearn_bh", performance)
print("Updated stLearn BH results; other saved methods unchanged.", flush=True)
