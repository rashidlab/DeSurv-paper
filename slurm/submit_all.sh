#!/bin/bash
# slurm/submit_all.sh — Submit the full DeSurv pipeline to SLURM
#
# Workflow (jobs 1-3 start immediately in parallel):
#   1. Main pipeline    (steps 01-05, 08, 09a)        ─┐
#   2. Simulation array (step 07, 600 tasks)           ─┤→ 4. Merge + render (steps 07b, 09c, 10)
#   3. CV grid search   (step 06)                     ─┘
#
# Job 4 waits for all three to succeed (afterok).
#
# Usage:
#   cd /path/to/DeSurv-paper-rl
#   bash slurm/submit_all.sh

set -euo pipefail

PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJ_DIR"
mkdir -p "$PROJ_DIR/slurm/logs"

echo "=== Submitting DeSurv pipeline from: $PROJ_DIR ==="

MAIN_JOB=$(sbatch --parsable slurm/run_main_pipeline.sh)
echo "Main pipeline:    job $MAIN_JOB  (steps 01-05, 08, 09a)"

SIM_JOB=$(sbatch --parsable slurm/run_simulations.sh)
echo "Simulation array: job $SIM_JOB  (step 07, 600 tasks)"

CV_JOB=$(sbatch --parsable slurm/run_cv_grid.sh)
echo "CV grid search:   job $CV_JOB  (step 06, independent)"

MERGE_JOB=$(sbatch --parsable \
  --dependency=afterok:${MAIN_JOB}:${SIM_JOB}:${CV_JOB} \
  slurm/run_merge_simulations.sh)
echo "Merge + render:   job $MERGE_JOB (steps 07b, 09c, 10) — waits for $MAIN_JOB, $SIM_JOB, $CV_JOB"

echo ""
echo "=== All jobs submitted ==="
echo "Monitor: squeue -u \$(whoami)"
echo "Cancel:  scancel $MAIN_JOB $SIM_JOB $CV_JOB $MERGE_JOB"
