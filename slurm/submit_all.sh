#!/bin/bash
# slurm/submit_all.sh — Submit the DeSurv pipeline to SLURM
#
# Usage:
#   cd /work/users/a/y/ayoung31/DeSurv-paper-rl
#   bash slurm/submit_all.sh
#
# Submits the main pipeline (steps 1-5, 8-9).
# Steps 6-7 (sensitivity analysis, simulations) use precomputed results.

set -euo pipefail

PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJ_DIR"

mkdir -p "$PROJ_DIR/slurm/logs"

echo "=== Submitting DeSurv pipeline ==="

JOB_MAIN=$(sbatch --parsable slurm/run_main_pipeline.sh)
echo "Main pipeline: job $JOB_MAIN (steps 1-5, 8-9)"

echo ""
echo "=== Job submitted ==="
echo "Monitor with: squeue -u $(whoami)"
echo "Cancel with:  scancel $JOB_MAIN"
