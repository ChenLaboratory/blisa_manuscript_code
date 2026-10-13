"""Derive plotting summaries in memory from all saved Visium local CSVs."""

import numpy as np
import pandas as pd

METHODS = ["blisa", "stlearn_bh", "spatialdm_permutation"]
STEMS = dict(zip(METHODS, ["blisa_permutation", "stlearn_bh_adjusted", "spatialdm_permutation"]))
REGIONS = ["Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI"]


def summarize_saved_results(repo, map_ids=()):
    results = repo / "res/benchmarks/visium"
    data = repo / "Data/generated/benchmarks/visium"
    spots = pd.read_csv(data / "visium_benchmark_spot_annotations.csv", usecols=["spot_id", "region"])
    counts_path = results / "visium_lri_hotspot_counts_and_test_status_10000_samples.csv"
    native = pd.read_csv(counts_path, converters={"status": str})
    native["returned"] = native.returned.astype(str).str.lower().isin(["true", "1"])
    ids = sorted(native.lr_id.unique())
    lr_index = pd.Series(np.arange(len(ids)), index=ids)
    spot_index = pd.Series(np.arange(len(spots)), index=spots.spot_id)
    selections, counts, regional, maps = {}, [], [], {}
    for method in METHODS:
        map_frames = []
        selected = np.zeros((len(ids), len(spots)), dtype=bool)
        available = np.zeros_like(selected)
        pcol = "bh_adjusted_p_value" if method == "stlearn_bh" else "p_value"
        path = results / f"{STEMS[method]}_local_lri_results_10000_samples.csv.gz"
        for chunk in pd.read_csv(path, usecols=["lr_id", "spot_id", pcol, "tested", "significant"], chunksize=200000):
            i = chunk.lr_id.map(lr_index).to_numpy()
            j = chunk.spot_id.map(spot_index).to_numpy()
            ok = chunk.tested.eq(1) & np.isfinite(chunk[pcol]) & chunk.significant.notna()
            available[i, j] = ok
            selected[i, j] = ok & chunk.significant.eq(1)
            shown = chunk[chunk.lr_id.isin(map_ids)].copy()
            if not shown.empty:
                shown["available"] = shown.tested.eq(1) & np.isfinite(shown[pcol]) & shown.significant.notna()
                shown["selected"] = shown.available & shown.significant.eq(1)
                map_frames.append(shown.rename(columns={pcol: "display_p_value"}))
        if map_frames:
            for lr, frame in pd.concat(map_frames).groupby("lr_id"):
                assert frame.spot_id.is_unique
                maps[method, lr] = frame.set_index("spot_id").reindex(spots.spot_id)
        selections[method] = selected
        status = native[native.condition == method].set_index("lr_id").reindex(ids)
        returned = status.returned.to_numpy(bool)
        n_available = available.sum(axis=1)
        usable = returned & (n_available > 0)
        status["n_available"] = np.where(returned, n_available, np.nan)
        status["n_sig"] = np.where(usable, selected.sum(axis=1), np.nan)
        status["status"] = np.where(~returned, "NR", np.where(usable, "available", "NA"))
        counts.append(status.reset_index())
        for region in REGIONS:
            roi = spots.region.eq(region).to_numpy()
            hits = np.where(usable, selected[:, roi].sum(axis=1), np.nan)
            regional.append(pd.DataFrame(dict(condition=method, lr_id=ids, region=region,
                n_total=int(roi.sum()), n_sig=hits, proportion=hits/roi.sum(), status=status.status.to_numpy())))
    counts = pd.concat(counts, ignore_index=True)
    common = set.intersection(*[set(counts.loc[(counts.condition == m) & (counts.status == "available"), "lr_id"]) for m in METHODS])
    # Eligibility follows deterministic input/filtering, independent of hotspot counts.
    assert len(common) == 679
    lookup = counts.set_index(["condition", "lr_id"])
    overlap = []
    for a in METHODS:
        for b in METHODS:
            intersections = (selections[a] & selections[b]).sum(axis=1)
            unions = (selections[a] | selections[b]).sum(axis=1)
            for i, lr in enumerate(ids):
                x, y = lookup.loc[(a, lr)], lookup.loc[(b, lr)]
                state = "NR" if not (x.returned and y.returned) else ("available" if x.status == y.status == "available" else "NA")
                valid = state == "available"
                overlap.append(dict(lr_id=lr, condition_a=a, condition_b=b,
                    n_a=x.n_sig if valid else np.nan, n_b=y.n_sig if valid else np.nan,
                    intersection=intersections[i] if valid else np.nan,
                    union=unions[i] if valid else np.nan,
                    jaccard=intersections[i]/unions[i] if valid and unions[i] else np.nan,
                    status=state, common_three_eligible=lr in common))
    overlap = pd.DataFrame(overlap)
    pooled = overlap[overlap.common_three_eligible].groupby(["condition_a", "condition_b"], as_index=False)[["intersection", "union"]].sum()
    pooled["n_lr_eligible"], pooled["jaccard"] = len(common), pooled.intersection/pooled.union.replace(0, np.nan)
    # Full summaries are verified against native LR-level exports, but are not
    # written as redundant plot-specific results.
    old = native.set_index(["condition", "lr_id"])
    new = counts.set_index(["condition", "lr_id"]).reindex(old.index)
    assert np.allclose(new[["n_available", "n_sig"]], old[["n_available", "n_sig"]], equal_nan=True)
    assert np.array_equal(new.status, old.status)
    return counts, pd.concat(regional, ignore_index=True), overlap, pooled, maps
