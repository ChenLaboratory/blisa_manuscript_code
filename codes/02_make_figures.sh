#!/usr/bin/env bash
# Generate panels from existing processed inputs and CSV results.
# Usage: bash codes/02_make_figures.sh [all|xenium|visium|aligned_comparison|xenium_benchmark|visium_benchmark]
set -euo pipefail

# 1. Set the repository and interpreters; CCC_PYTHON runs plotting companions.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
export CCC_REPO_ROOT="$PWD"
export RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
export CCC_PYTHON="${CCC_PYTHON:-python3}"
workflow="${1:-all}"
case "$workflow" in
    all) workflows=(xenium visium aligned_comparison xenium_benchmark visium_benchmark) ;;
    xenium|visium|aligned_comparison|xenium_benchmark|visium_benchmark) workflows=("$workflow") ;;
    *) echo "Choose all, xenium, visium, aligned_comparison, xenium_benchmark, or visium_benchmark." >&2; exit 1 ;;
esac

# 2. Render the matching main and supplementary panels without inference.
for analysis in "${workflows[@]}"; do
    echo "Plotting: $analysis"
    case "$analysis" in
        xenium) script=01_plot_xenium_results_Figure2_FigureS1_S3.R ;;
        visium) script=02_plot_visium_results_Figure3a_e_FigureS4.R ;;
        aligned_comparison) script=03_plot_aligned_results_Figure3f_h_FigureS5.R ;;
        xenium_benchmark) script=04_plot_xenium_benchmark_Figure4a_c_g_FigureS6_S8a.R ;;
        visium_benchmark) script=05_plot_visium_benchmark_Figure4d_f_h_i_FigureS7_S8b_c.R ;;
    esac
    "$RSCRIPT_BIN" --vanilla "codes/plotting/$script"
done
