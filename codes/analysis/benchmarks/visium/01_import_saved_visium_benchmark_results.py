"""Import completed Visium benchmarks and summarize every LR/spot result.

Only BLISA permutation, stLearn BH-adjusted and SpatialDM permutation results
are retained. This script does not import inference packages or run a test.
Native local results are split by LR blocks to keep individual CSVs manageable;
all LR pairs and spot rows are retained, including unavailable/untested results.
"""
import argparse
import hashlib
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

METHODS = ["blisa", "stlearn_bh", "spatialdm_permutation"]
FOLDERS = dict(zip(METHODS, ["blisa", "stlearn", "spatialdm_permutation"]))
STEMS = dict(zip(METHODS, ["blisa_permutation", "stlearn_bh_adjusted", "spatialdm_permutation"]))
REGIONS = ["Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI"]


def read(path, **kwargs):
    return pd.read_csv(path, float_precision="round_trip", **kwargs)


def write(frame, path):
    frame.to_csv(path, index=False, na_rep="NA", float_format="%.17g")


def sha256(path):
    with path.open("rb") as handle:
        digest = hashlib.sha256()
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    source, repo = args.source_root.resolve(), args.repo_root.resolve()
    run = source / "benchmark_methods_on_visium/results/final_total200_max3_spots200"
    reference = source / "benchmark_methods_on_visium/manuscript_update_20261012/report_tables"
    data = repo / "Data/generated/benchmarks/visium"
    results = repo / "res/benchmarks/visium"
    data.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    provenance = []
    source_hashes = {}

    def record(original, target, detail):
        provenance.append(dict(repository_path=str(target.relative_to(repo)),
            source_file=str(original.relative_to(source)), source_sha256=source_hashes.setdefault(original, sha256(original)) if original not in source_hashes else source_hashes[original],
            sha256=sha256(target), bytes=target.stat().st_size, description=detail))

    def stage(original, name, detail):
        target = data / name
        shutil.copyfile(original, target)
        record(original, target, detail)

    stage(source / "blisa_on_visium/input/visium_qc.rds",
          "visium_qc_benchmark_18082_genes.rds",
          "Historical benchmark input: 18,082 gene rows and 4,412 QC spots; distinct from native Figure 3 preprocessing.")
    stage(run / "perm_10000/blisa/lr_candidates.csv", "visium_benchmark_single_subunit_lr_candidates.csv",
          "Complete 1,154 candidate single-subunit CellChat LR pairs, including native annotations.")
    stage(run / "perm_10000/blisa/spots.csv", "visium_benchmark_spot_annotations.csv",
          "Complete spot annotations, coordinates, pathology and cell-type proportions.")
    stage(run / "report_tables/spatial_roi_borders.csv", "visium_benchmark_registered_roi_borders.csv",
          "Registered published ROI border paths used in the completed report.")
    stage(source / "blisa_on_visium/input/spatial/scalefactors_json.json", "visium_benchmark_image_scalefactors.json",
          "Native full-resolution Visium spot footprint and image calibration.")
    spots = read(data / "visium_benchmark_spot_annotations.csv")
    spot_index = pd.Series(np.arange(len(spots)), index=spots.spot_id)
    candidates = read(data / "visium_benchmark_single_subunit_lr_candidates.csv")
    ids = sorted(candidates.lr_id)
    lr_index = pd.Series(np.arange(len(ids)), index=ids)
    assert len(ids) == 1154 and len(spots) == 4412
    selections, summaries, regional = {}, [], []
    for method in METHODS:
        folder = run / "perm_10000" / FOLDERS[method]
        status = read(folder / "pair_status.csv", keep_default_na=False).set_index("lr_id").reindex(ids)
        returned = status.returned.to_numpy(dtype=bool)
        available = np.zeros((len(ids), len(spots)), dtype=bool)
        selected = np.zeros_like(available)
        pcol = "adjusted_p_value" if method == "stlearn_bh" else "p_value"
        scol = "native_significant" if method == "stlearn_bh" else "raw_significant"
        fields = ["lr_id", "spot_id", pcol, "score", "tested", scol]
        if method == "blisa":
            fields.append("quadrant")
        paths = {}
        native = folder / "lr_spot.csv"
        for chunk in read(native, usecols=fields, chunksize=200000):
            i = chunk.lr_id.map(lr_index).to_numpy()
            j = chunk.spot_id.map(spot_index).to_numpy()
            assert not np.isnan(i).any() and not np.isnan(j).any()
            ok = chunk.tested.eq(True).to_numpy() & np.isfinite(chunk[pcol]) & chunk[scol].notna()
            hit = ok & chunk[scol].eq(1).to_numpy()
            available[i, j], selected[i, j] = ok, hit
            chunk = chunk.rename(columns={scol: "significant", pcol: "bh_adjusted_p_value" if method == "stlearn_bh" else "p_value"})
            # Significance is defined only for available local tests; unavailable
            # flags remain NA rather than becoming observed zero results.
            chunk["significant"] = np.where(ok, hit.astype(int), np.nan)
            chunk["tested"] = chunk.tested.astype(int)
            for block, part in chunk.groupby(i // 100, sort=True):
                path = results / f"{STEMS[method]}_local_lri_results_10000_samples_part{int(block)+1:02d}.csv.gz"
                first = path not in paths
                part.to_csv(path, mode="w" if first else "a", header=first, index=False,
                            na_rep="NA", float_format="%.17g", compression={"method": "gzip", "compresslevel": 6})
                paths[path] = True
        selections[method] = selected
        n_available, n_sig = available.sum(axis=1), selected.sum(axis=1)
        local_available = returned & (n_available > 0)
        display_status = np.where(~returned, "NR", np.where(local_available, "available", "NA"))
        summary = pd.DataFrame(dict(budget=10000, condition=method, lr_id=ids,
            returned=returned, n_available=np.where(returned, n_available, np.nan),
            n_sig=np.where(local_available, n_sig, np.nan), status=display_status,
            native_status=status.status.to_numpy(), native_reason=status.reason.to_numpy()))
        summaries.append(summary)
        for region in REGIONS:
            roi = spots.region.eq(region).to_numpy()
            hits = np.where(local_available, selected[:, roi].sum(axis=1), np.nan)
            regional.append(pd.DataFrame(dict(condition=method, lr_id=ids, region=region,
                n_total=int(roi.sum()), n_sig=hits, proportion=hits / roi.sum(), status=display_status)))
        for target in paths:
            record(native, target, f"All native {method} local tests; 100 LR identities per part, 10,000 samples.")
        print(f"Imported {method}: {n_available.sum():,} available local tests, {n_sig.sum():,} hotspots", flush=True)

    counts = pd.concat(summaries, ignore_index=True)
    regions = pd.concat(regional, ignore_index=True)
    returned_ids = set(counts.loc[counts.returned, "lr_id"])
    common = set.intersection(*[set(counts.loc[(counts.condition == m) & (counts.status == "available"), "lr_id"]) for m in METHODS])
    assert len(returned_ids) == 1152 and len(common) == 679
    write(counts, results / "visium_lri_hotspot_counts_and_test_status_10000_samples.csv")

    # Compare selections on the full shared QC spot grid. Pooled Jaccard is the
    # sum of intersections divided by the sum of unions across the common LRs,
    # including locally available LR pairs with zero hotspots.
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
    pooled["n_lr_eligible"], pooled["jaccard"] = len(common), pooled.intersection / pooled.union.replace(0, np.nan)

    # Only the manuscript's 10,000-sample resource measurements are retained.
    performance = read(reference / "performance.csv")
    performance = performance[performance.budget == 10000].copy()
    assert set(performance.condition) == set(METHODS)
    write(performance, results / "visium_benchmark_computation_resources_10000_samples.csv")
    global_source = run / "perm_10000/spatialdm_permutation/global_res.csv"
    target = results / "spatialdm_permutation_global_lri_results_10000_samples.csv"
    shutil.copyfile(global_source, target)
    record(global_source, target, "Saved SpatialDM global permutation p-values, BH FDR and gate selections.")

    # Assert that every derived count and overlap equals the published report.
    for frame, filename, keys, numeric in [
        (counts, "significant_spot_counts.csv", ["budget", "condition", "lr_id"], ["n_available", "n_sig"]),
        (regions, "regional_hotspot_proportions.csv", ["condition", "lr_id", "region"], ["n_total", "n_sig", "proportion"]),
        (overlap, "spatial_overlap_per_lr.csv", ["lr_id", "condition_a", "condition_b"], ["intersection", "union", "jaccard"]),
        (pooled, "spatial_overlap_pooled.csv", ["condition_a", "condition_b"], ["intersection", "union", "jaccard"]),
    ]:
        original = read(reference / filename, converters={"status": str})
        if "budget" in original:
            original = original[original.budget == 10000]
        current = frame.set_index(keys).reindex(pd.MultiIndex.from_frame(original[keys]))
        assert np.allclose(current[numeric], original[numeric], rtol=1e-12, atol=1e-12, equal_nan=True), filename
        if "status" in original:
            assert np.array_equal(current.status, original.status), filename
    write(pd.DataFrame(provenance), data / "visium_benchmark_import_provenance.csv")
    print("Verified all 10k counts, regional denominators and per-LR/pooled Jaccards against the saved report. No inference was run.", flush=True)


if __name__ == "__main__":
    main()
