#!/usr/bin/env bash
# Run preparation, analysis and figure generation for the selected workflow.
# Usage: bash codes/run_workflow.sh [workflow] [--reuse-stlearn-spatialdm]
# The optional flag reuses saved stLearn/SpatialDM inference and resource rows.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
workflow="${1:-all}"

# 1. Generate processed inputs and analysis results from the starting data.
bash codes/01_run_analysis.sh "$workflow" "${2:-}"

# 2. Use those inputs and results to generate the figure panels.
bash codes/02_make_figures.sh "$workflow"
