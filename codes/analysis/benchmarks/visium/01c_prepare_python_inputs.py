#!/usr/bin/env python3
"""Convert 01b_export_python_inputs.R outputs to the benchmark AnnData input.

This retains every staged gene/spot and performs no normalization or inference.
Only count/coordinate inputs used by the two native methods are required.
"""
import argparse
import os
from pathlib import Path

import anndata
import pandas as pd
import scipy.io
import scipy.sparse


parser = argparse.ArgumentParser(description=__doc__)
repo = Path(os.environ.get("CCC_REPO_ROOT", Path(__file__).resolve().parents[4])).resolve()
parser.add_argument("--input-dir", type=Path,
                    default=repo / "Data/generated/benchmarks/visium/python_inputs")
input_dir = parser.parse_args().input_dir

# 1. Read genes, spots and sparse raw counts in their exported order.
genes = pd.read_csv(input_dir / "genes.csv", dtype=str, keep_default_na=False)
spots = pd.read_csv(input_dir / "spots.csv", dtype={"spot_id": str})
counts = scipy.sparse.csr_matrix(scipy.io.mmread(input_dir / "counts.mtx").T)
# 2. Use the ideal 100-um Visium lattice as the inference coordinates.
xy = spots[["grid_x_um", "grid_y_um"]].to_numpy(float)
spots.index = pd.Index(spots.spot_id, name="barcode")
genes.index = pd.Index(genes.gene_symbol, name="gene")
adata = anndata.AnnData(X=counts, obs=spots, var=genes)
adata.layers["counts"] = counts.copy()
adata.obsm["spatial"] = xy
adata.uns["coordinate_contract"] = "Ideal 100-micrometre Visium lattice; not image pixels"
# 3. Store unnormalized counts; each method applies its own normalization.
output = input_dir / "visium_qc.h5ad"
adata.write_h5ad(output, compression="gzip")
print(f"Saved {output}: {adata.n_obs} spots, {adata.n_vars} genes")
