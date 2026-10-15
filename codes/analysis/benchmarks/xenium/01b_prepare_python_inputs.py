#!/usr/bin/env python
"""Build raw AnnData after 01a_export_python_inputs.R; no inference is run.

The exported counts and metadata come from Data/generated/xenium/xenium_50um_binned.rds.
Gene/bin order, separate DCIS labels, and micrometre centroids are retained
without normalization.
Run in either the stLearn or SpatialDM environment, from any working directory.
"""
import hashlib
import json
import anndata
import pandas as pd

from xenium_benchmark_helpers import INPUT, CELLCHAT

# 1. Read the matched raw-count and metadata tables written by R.
counts_path = INPUT / "xenium_benchmark_raw_counts.csv.gz"
metadata_path = INPUT / "xenium_benchmark_bin_metadata.csv"
counts = pd.read_csv(counts_path, index_col="sample_id")
metadata = pd.read_csv(metadata_path, index_col="sample_id")
values = counts.to_numpy()
coordinates = metadata[["x", "y"]].to_numpy()

# 2. Store raw counts and micrometre centroids without normalization.
adata = anndata.AnnData(X=values, obs=metadata.copy(),
    var=pd.DataFrame(index=counts.columns), obsm={"spatial": coordinates},
    dtype=values.dtype)
output = INPUT / "xenium_benchmark_raw_counts.h5ad"
adata.write_h5ad(output, compression="gzip")

# 3. Record input versions and checksums alongside the generated AnnData.
files = [counts_path, metadata_path, output, CELLCHAT]
checksums = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
provenance = dict(source="Data/generated/xenium/xenium_50um_binned.rds", n_bins=adata.n_obs,
    n_genes=adata.n_vars, coordinate_units="micrometres", anndata_version=anndata.__version__,
    preparation="Raw counts, coordinates, identifiers and metadata retained", sha256=checksums)
(INPUT / "xenium_benchmark_python_input_provenance.json").write_text(
    json.dumps(provenance, indent=2) + "\n")
print(f"Prepared unchanged raw AnnData: {adata.n_obs} bins x {adata.n_vars} genes")
