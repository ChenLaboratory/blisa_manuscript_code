"""Read shared inputs, measure method resources and export complete result CSVs."""
import gc
import gzip
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd

NPERM = 10000
SEED = 20260915
RADIUS_UM = 250.0
CONTACT_RADIUS_UM = 120.0
RESOURCE_COLUMNS = ["method", "condition", "budget", "nperm", "ncpu",
    "elapsed_seconds", "peak_rss_mib", "rss_peak_reset", "rss_scope", "timing_scope"]


def load_inputs():
    import anndata

    repo = Path(os.environ.get("CCC_REPO_ROOT", Path(__file__).resolve().parents[4])).resolve()
    data = repo / "Data/generated/benchmarks/visium"
    output = repo / "res/benchmarks/visium"
    ncpu = int(os.environ.get("NCPU", "40"))
    adata = anndata.read_h5ad(data / "python_inputs/visium_qc.h5ad")
    adata.X = adata.layers["counts"].copy()
    lr_table = data / "visium_benchmark_single_subunit_lr_candidates.csv"
    candidates = pd.read_csv(lr_table, dtype=str, keep_default_na=False).set_index("cellchat_id", drop=False)
    output.mkdir(parents=True, exist_ok=True)
    return adata, candidates, output, ncpu, lr_table


def process_memory_mib(field):
    fields = dict(line.split(":", 1) for line in Path("/proc/self/status").read_text().splitlines())
    return int(fields[field].split()[0]) / 1024


def start_measurement():
    gc.collect()
    # Linux may deny peak reset; record that rather than reporting an old peak.
    try:
        Path("/proc/self/clear_refs").write_text("5\n")
        reset = True
    except OSError:
        reset = False
    return reset, time.perf_counter()


def finish_measurement(start):
    reset, clock_start = start
    elapsed = time.perf_counter() - clock_start
    peak = process_memory_mib("VmHWM") if reset else float("nan")
    return dict(elapsed_seconds=elapsed, peak_rss_mib=peak, rss_peak_reset=reset,
                rss_scope="Stage process peak; includes threads, excludes separate workers")


def update_method_rows(path, replacement, condition):
    """Replace only this method/budget, preserving the other saved method rows."""
    if path.exists():
        old = pd.read_csv(path, converters={"status": str}, float_precision="round_trip")
        keep = ~((old.condition == condition) & (old.budget == NPERM))
        replacement = pd.concat([old.loc[keep, replacement.columns], replacement], ignore_index=True)
    replacement.to_csv(path, index=False, na_rep="NA", float_format="%.17g")


def export_local_results(output, condition, stem, p_column, candidates, spots,
                         local_frames, native_status):
    """Write all candidates × spots, including missing/native-filtered tests.

    local_frames yields (cellchat_id, p_values, lr_scores, tested, significant).
    Only stLearn supplies lr_scores; SpatialDM supplies its p-values and flags.
    Significance is stLearn pooled BH<0.05 or SpatialDM permutation p<0.05.
    The SpatialDM native 0.1 selection is intentionally not the report flag.
    """
    target = output / f"{stem}_local_lri_results_10000_samples.csv.gz"
    summaries = []
    with gzip.open(target, "wt", compresslevel=6) as handle:
        first = True
        for pair, pvalues, scores, tested, significant in local_frames:
            candidate = candidates.loc[pair]
            state = native_status.loc[pair]
            ok = np.asarray(tested, bool) & np.isfinite(pvalues)
            hits = ok & np.asarray(significant, bool)
            frame = pd.DataFrame(dict(lr_id=candidate.lr_id, spot_id=spots))
            frame[p_column] = np.where(ok, pvalues, np.nan)
            # Retain stLearn's native LR score; SpatialDM exports its local test.
            if condition == "stlearn_bh":
                frame["lr_score"] = scores
            frame["tested"] = np.asarray(tested, int)
            frame["significant"] = np.where(ok, hits.astype(int), np.nan)
            frame.to_csv(handle, index=False, header=first, na_rep="NA", float_format="%.17g")
            first = False
            returned = bool(state.returned)
            n_available = int(ok.sum())
            summaries.append(dict(budget=NPERM, condition=condition, lr_id=candidate.lr_id,
                returned=returned, n_available=n_available if returned else np.nan,
                n_sig=int(hits.sum()) if returned and n_available else np.nan,
                status="NR" if not returned else ("available" if n_available else "NA"),
                native_status=state.status))
    summary = pd.DataFrame(summaries).sort_values("lr_id")
    update_method_rows(output / "visium_lri_hotspot_counts_and_test_status_10000_samples.csv",
                       summary, condition)


def export_performance(output, condition, performance):
    performance.update(condition=condition, budget=NPERM)
    update_method_rows(output / "visium_benchmark_computation_resources_10000_samples.csv",
                       pd.DataFrame([performance], columns=RESOURCE_COLUMNS), condition)


def stlearn_local_frames(adata, candidates):
    """Convert the native pooled matrices without dropping untested candidates."""
    positions = {pair: i for i, pair in enumerate(adata.uns["lr_summary"].index)}
    for pair, candidate in candidates.iterrows():
        if candidate.native_lr_id in positions:
            column = positions[candidate.native_lr_id]
            score = np.asarray(adata.obsm["lr_scores"][:, column]).reshape(-1)
            p = np.asarray(adata.obsm["p_vals"][:, column]).reshape(-1)
            adjusted = np.asarray(adata.obsm["p_adjs"][:, column]).reshape(-1)
            tested = (score > 0) & np.isfinite(p)
            yield pair, adjusted, score, tested, tested & (adjusted < 0.05)
        else:
            missing = np.full(adata.n_obs, np.nan)
            no_test = np.zeros(adata.n_obs, bool)
            yield pair, missing, missing, no_test, no_test


def spatialdm_local_frames(adata, candidates):
    global_results = adata.uns["global_res"]
    local_ids = global_results.index[global_results.selected]
    for pair in candidates.index:
        if pair in local_ids:
            p = adata.uns["local_perm_p"].loc[pair, adata.obs_names].to_numpy()
            tested = np.isfinite(p)
            # The manuscript comparison uses local p<0.05 after global selection.
            yield pair, p, None, tested, tested & (p < 0.05)
        else:
            missing = np.full(adata.n_obs, np.nan)
            no_test = np.zeros(adata.n_obs, bool)
            yield pair, missing, missing, no_test, no_test
