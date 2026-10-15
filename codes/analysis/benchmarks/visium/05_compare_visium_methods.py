"""Compare saved 10,000-sample Visium calls without rerunning any method."""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

# 1. Read the native method summaries and processed spot-region annotations.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo-root", type=Path, default=Path(os.environ.get("CCC_REPO_ROOT", ".")))
args = parser.parse_args()
results = args.repo_root / "res/benchmarks/visium"
data = args.repo_root / "Data/generated/benchmarks/visium"
methods = ["blisa", "stlearn_bh", "spatialdm_permutation"]
stems = dict(zip(methods, ["blisa_permutation", "stlearn_bh_adjusted", "spatialdm_permutation"]))
regions = ["Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI"]
spots = pd.read_csv(data / "visium_benchmark_spot_annotations.csv", usecols=["spot_id", "region"])
counts = pd.read_csv(results / "visium_lri_hotspot_counts_and_test_status_10000_samples.csv",
                     converters={"status": str})
ids = sorted(counts.lr_id.unique())
lr_index = pd.Series(np.arange(len(ids)), index=ids)
spot_index = pd.Series(np.arange(len(spots)), index=spots.spot_id)
lookup = counts.set_index(["condition", "lr_id"])

# 2. Read the saved native significance calls in bounded chunks.
# Availability preserves failed/unreturned tests separately from nonsignificance.
# BLISA's High-High classification is already encoded in its saved significant flag.
selections, regional = {}, []
for method in methods:
    selected = np.zeros((len(ids), len(spots)), dtype=bool)
    pcol = "bh_adjusted_p_value" if method == "stlearn_bh" else "p_value"
    path = results / f"{stems[method]}_local_lri_results_10000_samples.csv.gz"
    for chunk in pd.read_csv(path, usecols=["lr_id", "spot_id", pcol, "tested", "significant"],
                             chunksize=200000):
        i = chunk.lr_id.map(lr_index).to_numpy()
        j = chunk.spot_id.map(spot_index).to_numpy()
        available = chunk.tested.eq(1) & np.isfinite(chunk[pcol]) & chunk.significant.notna()
        selected[i, j] = available & chunk.significant.eq(1)
    selections[method] = selected

    # Coverage divides the selected spots by all spots in the annotated region.
    status = lookup.loc[method].reindex(ids)
    usable = status.returned.to_numpy(bool) & (status.n_available.to_numpy() > 0)
    for region in regions:
        roi = spots.region.eq(region).to_numpy()
        hits = np.where(usable, selected[:, roi].sum(axis=1), np.nan)
        regional.append(pd.DataFrame(dict(method=method, lr_id=ids, region=region,
            n_total=int(roi.sum()), n_sig=hits, proportion=hits / roi.sum(),
            status=status.status.to_numpy())))

# 3. Define shared eligibility independently of hotspot counts.
common = set.intersection(*[set(counts.loc[
    (counts.condition == method) & (counts.status == "available") & (counts.n_available > 0), "lr_id"])
    for method in methods])

# 4. Compare all LR candidates for each ordered method pair.
# Unreturned/unavailable tests and empty significant sets retain their distinct meaning.
overlap = []
for a in methods:
    for b in methods:
        intersections = (selections[a] & selections[b]).sum(axis=1)
        unions = (selections[a] | selections[b]).sum(axis=1)
        for i, lr in enumerate(ids):
            x, y = lookup.loc[(a, lr)], lookup.loc[(b, lr)]
            state = "NR" if not (x.returned and y.returned) else (
                "available" if x.status == y.status == "available" else "NA")
            valid = state == "available"
            overlap.append(dict(lr_id=lr, method_a=a, method_b=b,
                n_a=x.n_sig if valid else np.nan, n_b=y.n_sig if valid else np.nan,
                intersection=intersections[i] if valid else np.nan,
                union=unions[i] if valid else np.nan,
                jaccard=intersections[i] / unions[i] if valid and unions[i] else np.nan,
                status=state, common_three_eligible=lr in common))
overlap = pd.DataFrame(overlap)

# 5. Pool intersections and unions across common LRs, preserving LR-spot identity.
pooled = overlap[overlap.common_three_eligible].groupby(
    ["method_a", "method_b"], as_index=False)[["intersection", "union"]].sum()
pooled["n_lr_eligible"] = len(common)
pooled["jaccard"] = pooled.intersection / pooled.union.replace(0, np.nan)
pooled["status"] = np.where(pooled.union > 0, "available", "no_significant_spots")
overlap.to_csv(results / "visium_benchmark_jaccard_by_lr.csv", index=False, na_rep="NA")
pooled.to_csv(results / "visium_benchmark_jaccard_pooled.csv", index=False, na_rep="NA")
pd.concat(regional, ignore_index=True).to_csv(
    results / "visium_benchmark_regional_coverage.csv", index=False, na_rep="NA")
