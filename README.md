# BLISA manuscript code

Analysis and plotting code for the Xenium, Visium, aligned-platform and
three-method benchmark analyses. Plotting reads processed inputs and saved
results and does not run permutation inference.

This GitHub repository contains code, documentation, software versions and
selected data/results. Published data includes `Data/data_sources.csv`,
`Data/parameters/`, `Data/downloaded/`, and CSV/JSON files in `Data/generated/`.
Published `res/` files are individually smaller than 10,000,000 bytes, using
compressed file sizes for `.csv.gz` files. Generated RDS inputs and larger
result files remain local, together with figures and the private methods
reference. To reproduce all figures, the omitted inputs and results must first
be supplied in the working directory shown below.

## Working directory layout

```text
Data/
  downloaded/              # CellChat GitHub snapshot; Janesick publisher workbook
  parameters/              # Published registration and coordinate calibration
  generated/               # Processed expression, geometry, annotations and LR inputs
codes/
  00_import_saved_results.R
  01_run_analysis.R
  02_make_figures.R
  analysis/
    xenium/
    visium/
    aligned_comparison/
    benchmarks/
      xenium/
      visium/
  plotting/                # Main and supplementary panel scripts
  helpers/                 # Shared data and plotting functions
res/
  xenium/
  visium/
  aligned_overlap/
  benchmarks/
    xenium/
    visium/
Figures/                   # Regenerated PDFs; excluded from Git
software_versions.csv
```

`Data/data_sources.csv` records sources, sizes and SHA-256 checksums, including
files retained only locally. Original release objects are linked there and are
not duplicated. `Data/generated/`
contains the actual processed inputs; spatial annotations are inputs, not LRI
results. `res/` contains CSV analysis results and computation measurements.
Large complete local-test tables are gzip-compressed CSVs, with one file per
analysis or method. Results are never split to meet the upload cutoff: files
of 10,000,000 bytes or more remain local. No rows are dropped to reduce their
size, and no result RDS is required for plotting. Regional coverage and
agreement summaries are calculated from the complete results when plotting.

The standalone Visium spot results are stored in
`res/visium/lri_by_spot.csv.gz`. Each row is one
LR pair at one spot, including its local test statistics and native hotspot call.
The three Visium benchmark methods each have a single
`res/benchmarks/visium/*_local_lri_results_10000_samples.csv.gz` file. These four
complete tables exceed the upload cutoff and must be supplied locally for
plotting. `res/visium/lri_summary.csv` contains the separate per-LR summary and
cannot replace the spot-level results.

## Reproduce the figures

Run commands from this directory. `software_versions.csv` lists the main software
and versions by analysis. Each R script checks its required packages at startup.
R scripts call Python companions for the
original benchmark spatial and heatmap styles. Set `CCC_PYTHON` to a Python
executable with numpy, pandas, matplotlib, scipy and scikit-image installed.

```sh
export CCC_PYTHON=/path/to/python
Rscript codes/02_make_figures.R all
```

On the original cluster, the R environment also requires:

```sh
module load gdal geos proj/9.4.0 mariadb-connector-c/3.3.10 glpk/5.0 ImageMagick/7.1.1 R/4.4.2
```

| Plotting script | Figure coverage |
| --- | --- |
| `01_plot_xenium_results_Figure2_FigureS1_S3.R` | Figure 2; S1–S3 |
| `02_plot_visium_results_Figure3a_e_FigureS4.R` | Figure 3a–e; S4a–b |
| `03_plot_aligned_results_Figure3f_h_FigureS5.R` | Figure 3f–h; S5a–b |
| `04_plot_xenium_benchmark_Figure4a_c_g_FigureS6_S8a.R` | Figure 4a–c,g; S6; S8a |
| `05_plot_visium_benchmark_Figure4d_f_h_i_FigureS7_S8b_c.R` | Figure 4d–f,h–i; S7a–b; S8b–c |

Plots are written separately as `Figures/figure3a.pdf`, `figure3b.pdf`, etc.
Figures 3, 4, S4, S5, S7 and S8 are not assembled into whole figures, and no
supplementary document is generated. S1–S3 and S6 retain their multi-map layouts.
Figure 3f has no overlap outline or rectangular frame; the manuscript outline
is added manually in BioRender. S8c follows page 28 of the updated Visium
benchmark report, using the 679 shared LR pairs.

To render one workflow, replace `all` with `xenium`, `visium`,
`aligned_comparison`, `xenium_benchmark` or `visium_benchmark`.

The Xenium plotting file also accepts `--main-only` or `--supp-only`.
The included plotting scripts generate each plot from `Data/` and `res/`.
They do not read existing figure PDFs or extract images from the manuscript.

**Layout provenance:** S1–S3 have been reconstructed from their existing combined
PDF references. The original assembly code is still being located for exact
layout and dendrogram-order matching. Their underlying data and native hotspot
calls are preserved.

## Results and analysis settings

Figures use native BLISA calls, stLearn **BH-adjusted** local results and
SpatialDM **permutation** results. Raw-stLearn and SpatialDM z-score modes are
excluded. Benchmark inference used 10,000 samples; common-LR comparisons use 12
Xenium pairs and 679 Visium pairs. Availability/untested states remain distinct
from a tested pair with no hotspots. SpatialDM's local p-value comparison also
retains the original global selection gate.

Standalone Xenium, standalone Visium and aligned BLISA runs preserve their
executed `nsim = 9999` and four-thread settings. The benchmark BLISA settings
are 10,000 permutations with 10 CPUs for Xenium and 40 for Visium. These are
separate analyses, including distinct native versus duplicate-collapsed Visium
inputs. Historical computation measurements retain their original timing and
memory scope; a fresh run on other hardware will not reproduce those timings.

## Import saved results or rerun BLISA

The locally saved CSVs are sufficient for figures. To recover them from the original
workspace, without rerunning any method:

```sh
Rscript codes/00_import_saved_results.R /path/to/original/workspace
```

This imports saved Visium/aligned BLISA and both benchmark workflows. In
particular, stLearn and SpatialDM are always imported from their completed runs.
It never starts either inference method. Original source paths and hashes are
recorded in the imported metadata.

For fresh standalone BLISA analysis, provide the original release files through
`CCC_XENIUM_RDS` and `CCC_VISIUM_RDS` (or `CCC_SOURCE_DATA_DIR`) and run the
numbered workflow. The helper also checks the existing BLISA cache and never
downloads data automatically.

```sh
export CCC_XENIUM_RDS=/path/to/spe_xenium_bc_s1rep1.rds
export CCC_VISIUM_RDS=/path/to/spe_visium_bc_s1.rds
Rscript codes/01_run_analysis.R xenium
Rscript codes/01_run_analysis.R visium
Rscript codes/01_run_analysis.R aligned_comparison
```

Optional benchmark BLISA reruns are explicit commands:

```sh
Rscript codes/analysis/benchmarks/xenium/02_run_blisa_benchmark.R
Rscript codes/analysis/benchmarks/visium/02_run_blisa_visium_benchmark.R
```

They replace only BLISA rows and measurements, preserving imported stLearn and
SpatialDM results. The Xenium native database loader did not record its original
snapshot or an explicit R seed; the rerun checks returned LR identities against
the saved run. A fresh inference run has not been validated in this extraction.

## Validation and provenance

The current figures were rendered from saved results, without permutation
reruns. Imports checked native hotspot calls, LR/spot identities, availability,
regional coverage and benchmark overlap against the original completed analyses.
The complete panel command was also run with an initially empty `Figures/`
directory and access only to the repository code, local processed data and CSV results.
It generated all 39 plot PDFs; the data and result files remained unchanged.
The Visium export contains all 7,129,792 LR–spot rows; the aligned comparison
preserves both 15-LR × 3,762-spot arms. Benchmark selections and common-pair
Jaccard values match the updated report tables.

`software_versions.csv` has three columns: analysis, software and version.
It lists the main analysis tools and plotting runtimes, not their full dependency
trees or environment settings.
`METHODS_writing_reference.md` is a private, Gitignored writing reference.
