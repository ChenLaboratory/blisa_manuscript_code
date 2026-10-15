#!/usr/bin/env bash
# Prepare inputs and run the selected analysis with its R/Python interpreter.
# Usage: bash codes/01_run_analysis.sh [workflow] [--reuse-stlearn-spatialdm]
set -euo pipefail

# 1. Use the repository directory for all relative input and output paths.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
export CCC_REPO_ROOT="$PWD"
export RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
workflow="${1:-all}"
reuse_python_results="${2:-}"
case "$reuse_python_results" in
    ""|--reuse-stlearn-spatialdm) ;;
    *) echo "Optional second argument: --reuse-stlearn-spatialdm" >&2; exit 1 ;;
esac
case "$workflow" in
    all) workflows=(xenium visium aligned_comparison xenium_benchmark visium_benchmark) ;;
    xenium|visium|aligned_comparison|xenium_benchmark|visium_benchmark) workflows=("$workflow") ;;
    *) echo "Choose all, xenium, visium, aligned_comparison, xenium_benchmark, or visium_benchmark." >&2; exit 1 ;;
esac

# These environments also perform AnnData conversion and saved-result comparison.
case "$workflow" in
    all|xenium_benchmark|visium_benchmark)
        : "${CCC_STLEARN_PYTHON:?Set CCC_STLEARN_PYTHON to the stLearn Python executable}"
        : "${CCC_SPATIALDM_PYTHON:?Set CCC_SPATIALDM_PYTHON to the SpatialDM Python executable}"
        ;;
esac

# 2. Run each preparation and analysis step in its own process, in order.
for analysis in "${workflows[@]}"; do
    echo "Running analysis: $analysis"
    case "$analysis" in
        xenium)
            # Bin cells, add annotations, then calculate native LRI and CCI.
            "$RSCRIPT_BIN" --vanilla codes/analysis/xenium/01_prepare_xenium_50um_bins.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/xenium/02_run_xenium_lri_cci.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/xenium/03_summarize_xenium_lri_cci.R
            ;;
        visium)
            # Apply QC, prepare LR candidates, then calculate native LRI.
            "$RSCRIPT_BIN" --vanilla codes/analysis/visium/01_prepare_visium_qc_annotations.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/visium/02_prepare_visium_lr_candidates.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/visium/03_run_visium_lri.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/visium/04_summarize_visium_lri_pathways.R
            ;;
        aligned_comparison)
            # Align the two platforms, select shared pairs, then test both arms.
            "$RSCRIPT_BIN" --vanilla codes/analysis/aligned_comparison/01_prepare_aligned_overlap.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/aligned_comparison/02_prepare_shared_lr_pairs.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/aligned_comparison/03_run_aligned_lri.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/aligned_comparison/04_compare_aligned_lri.R
            ;;
        xenium_benchmark)
            # Reuse standalone Xenium bins and export their raw counts to AnnData.
            "$RSCRIPT_BIN" --vanilla codes/analysis/benchmarks/xenium/01a_export_python_inputs.R
            "$CCC_STLEARN_PYTHON" codes/analysis/benchmarks/xenium/01b_prepare_python_inputs.py
            "$RSCRIPT_BIN" --vanilla codes/analysis/benchmarks/xenium/02_run_blisa_benchmark.R
            if [[ "$reuse_python_results" != --reuse-stlearn-spatialdm ]]; then
                "$CCC_STLEARN_PYTHON" codes/analysis/benchmarks/xenium/03_run_stlearn_benchmark.py
                "$CCC_SPATIALDM_PYTHON" codes/analysis/benchmarks/xenium/04_run_spatialdm_permutation_benchmark.py
            fi
            # Compare fresh BLISA calls with saved or newly computed method calls.
            "$RSCRIPT_BIN" --vanilla codes/analysis/benchmarks/xenium/05_compare_xenium_methods.R
            ;;
        visium_benchmark)
            # Reuse standalone Visium QC/annotations; collapse duplicate symbols.
            "$RSCRIPT_BIN" --vanilla codes/analysis/benchmarks/visium/01a_prepare_visium_benchmark.R
            "$RSCRIPT_BIN" --vanilla codes/analysis/benchmarks/visium/01b_export_python_inputs.R
            "$CCC_STLEARN_PYTHON" codes/analysis/benchmarks/visium/01c_prepare_python_inputs.py
            "$RSCRIPT_BIN" --vanilla codes/analysis/benchmarks/visium/02_run_blisa_visium_benchmark.R
            if [[ "$reuse_python_results" != --reuse-stlearn-spatialdm ]]; then
                "$CCC_STLEARN_PYTHON" codes/analysis/benchmarks/visium/03_run_stlearn_visium_benchmark.py
                "$CCC_SPATIALDM_PYTHON" codes/analysis/benchmarks/visium/04_run_spatialdm_permutation_visium_benchmark.py
            fi
            "$CCC_SPATIALDM_PYTHON" codes/analysis/benchmarks/visium/05_compare_visium_methods.py
            ;;
    esac
done
