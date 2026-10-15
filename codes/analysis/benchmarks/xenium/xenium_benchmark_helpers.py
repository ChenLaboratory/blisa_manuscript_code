"""Shared input, measurement and CSV export for the two explicit Python reruns.

No method inference runs on import. Existing rows belonging to other methods
are retained; run the numbered methods sequentially when sharing an output path.
"""
import gc
import io
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("CCC_REPO_ROOT", Path(__file__).resolve().parents[4]))
INPUT = ROOT / "Data/generated/benchmarks/xenium"
OUTPUT = ROOT / "res/benchmarks/xenium"
CELLCHAT = ROOT / "Data/downloaded/cellchat_github/CellChatDB.human.rda"
METHODS = ("blisa", "stlearn", "spatialdm")
LOCAL_COLUMNS = ["method", "lr_id", "sample_id", "p_value", "local_available",
    "significant", "lr_score", "quadrant", "bin_id"]
SUMMARY_COLUMNS = ["method", "lr_id", "rank", "returned", "n_available",
    "significant_bins", "status"]
RESOURCE_COLUMNS = ["method", "condition", "nperm", "ncpu", "elapsed_seconds",
    "peak_rss_mib", "rss_peak_reset", "rss_scope", "timing_scope"]


def cellchat_candidates(genes):
    """Decode the frozen database with R, excluding complexes on either side."""
    # Read the downloaded CellChat snapshot, excluding multi-subunit complexes.
    r_code = '''
        env <- new.env()
        load(commandArgs(trailingOnly = TRUE)[1], envir = env)
        db <- env$CellChatDB.human
        x <- db$interaction
        keep <- !x$ligand %in% rownames(db$complex) & !x$receptor %in% rownames(db$complex)
        x$cellchat_id <- rownames(x)
        write.csv(x[keep, , drop = FALSE], stdout(), row.names = FALSE, na = "")
    '''
    result = subprocess.run([os.environ.get("RSCRIPT_BIN", "Rscript"), "--vanilla",
        "-e", r_code, str(CELLCHAT)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    table = pd.read_csv(io.StringIO(result.stdout.decode()), index_col="cellchat_id",
        keep_default_na=False)
    table = table[table.ligand.isin(genes) & table.receptor.isin(genes)]
    table = table.drop_duplicates(["ligand", "receptor"])
    return table


def process_memory_mib(field):
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith(field + ":"):
                return int(line.split()[1]) / 1024
    except (OSError, ValueError):
        pass
    return float("nan")


def start_measurement():
    # Same process-level Linux RSS accounting as the original runs: includes
    # threads, excludes separate workers; imports/input/database decoding excluded.
    gc.collect()
    try:
        Path("/proc/self/clear_refs").write_text("5\n")
        reset = True
    except OSError:
        reset = False
    return reset, time.perf_counter()


def finish_measurement(started):
    reset, clock_start = started
    elapsed = time.perf_counter() - clock_start
    peak = process_memory_mib("VmHWM") if reset else float("nan")
    return dict(elapsed_seconds=elapsed, peak_rss_mib=peak, rss_peak_reset=reset,
        rss_scope="Stage process peak; includes threads, excludes separate workers")


def result_tables(method, metadata, pairs, p_values, scores=None, global_results=None):
    """Convert complete native results using the manuscript's native test rules."""
    # 1. Map each native LR pair to its complete vector of local tests.
    local, summary = [], []
    for native_id, pair in pairs.iterrows():
        lr_id = pair.ligand + "__" + pair.receptor
        p = (p_values[native_id].reindex(metadata.index).to_numpy(float)
            if native_id in p_values else np.full(len(metadata), np.nan))
        available = np.isfinite(p)
        score = np.full(len(metadata), np.nan)
        # 2. Preserve method-specific eligibility before counting hotspots.
        if method == "stlearn":
            score = scores[native_id].reindex(metadata.index).to_numpy(float)
            available &= np.isfinite(score) & (score > 0)
        else:
            global_fdr = global_results.loc[native_id, "fdr"]
            available &= np.isfinite(global_fdr) & (global_fdr < 0.10)
        significant = available & (p < 0.05)
        # 3. Keep unavailable bins distinct from tested bins with no hotspot.
        local.append(pd.DataFrame(dict(method=method, lr_id=lr_id,
            sample_id=metadata.index, p_value=np.where(available, p, np.nan),
            local_available=available, significant=significant,
            lr_score=score, quadrant=np.nan,
            bin_id=metadata.bin_id.to_numpy())))
        summary.append(dict(method=method, lr_id=lr_id, returned=True,
            n_available=int(available.sum()),
            significant_bins=int(significant.sum()) if available.any() else np.nan,
            status="available" if available.any() else "no_local_results"))
    return (pd.concat(local, ignore_index=True) if local else pd.DataFrame(columns=LOCAL_COLUMNS),
        pd.DataFrame(summary, columns=[name for name in SUMMARY_COLUMNS if name != "rank"]))


def update_outputs(method, local, summary, performance, output=OUTPUT):
    """Replace only this method; support both cached and initially empty results."""
    # 1. Combine this run with any methods that have already completed.
    output.mkdir(parents=True, exist_ok=True)
    local_path = output / "xenium_benchmark_lri_by_bin.csv.gz"
    summary_path = output / "xenium_benchmark_lr_summary.csv"
    performance_path = output / "xenium_benchmark_computation_resources.csv"
    read = lambda path: pd.read_csv(path, float_precision="round_trip", low_memory=False)
    old_summary = read(summary_path) if summary_path.exists() else pd.DataFrame(columns=SUMMARY_COLUMNS)
    order = old_summary.sort_values("rank").lr_id.drop_duplicates().tolist()
    order += [lr for lr in summary.lr_id if lr not in order]
    old_other = old_summary[old_summary.method != method]
    combined = pd.concat([old_other, summary], ignore_index=True)
    # 2. Explicitly mark methods not run on a new input; distinguish this from a
    # completed method that did not return an LR pair.
    measured_methods = set(old_summary.loc[old_summary.status != "not_run", "method"]) | {method}
    present = set(zip(combined.method, combined.lr_id))
    missing = [dict(method=m, lr_id=lr, returned=False,
        status="pair_not_returned" if m in measured_methods else "not_run")
        for lr in order for m in METHODS if (m, lr) not in present]
    combined = pd.concat([combined, pd.DataFrame(missing)], ignore_index=True)
    combined["rank"] = combined.lr_id.map({lr: i + 1 for i, lr in enumerate(order)})
    combined = combined.reindex(columns=SUMMARY_COLUMNS).sort_values(["rank", "method"])
    old_local = read(local_path) if local_path.exists() else pd.DataFrame(columns=LOCAL_COLUMNS)
    merged_local = pd.concat([old_local[old_local.method != method], local], ignore_index=True)
    old_performance = read(performance_path) if performance_path.exists() else pd.DataFrame(columns=["method"])
    merged_performance = pd.concat([old_performance[old_performance.method != method],
        pd.DataFrame([performance])], ignore_index=True)
    # 3. Save unsplit tables; 17 digits retain p-values close to the cutoff.
    merged_local.reindex(columns=LOCAL_COLUMNS).to_csv(local_path, index=False,
        na_rep="NA", float_format="%.17g", compression="gzip")
    combined.to_csv(summary_path, index=False, na_rep="NA", float_format="%.17g")
    merged_performance.reindex(columns=RESOURCE_COLUMNS).to_csv(performance_path,
        index=False, na_rep="NA", float_format="%.17g")
    print(f"Updated {method} in {output}; other available methods retained")
