#!/bin/bash
#SBATCH --job-name=desurv-sim
#SBATCH --partition=general
#SBATCH --array=0-599
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=slurm/logs/sim_%A_%a.log

# 600 tasks = 3 scenarios x 2 analyses x 100 replicates
# Layout: task_id = scenario_idx * 200 + analysis_idx * 100 + (replicate - 1)
SCENARIOS=(R0_easy R00_null R_mixed)
ANALYSES=(bo_tune_ntop bo_tune_ntop_alpha0)

scenario_idx=$(( SLURM_ARRAY_TASK_ID / 200 ))
analysis_idx=$(( (SLURM_ARRAY_TASK_ID % 200) / 100 ))
replicate=$(( (SLURM_ARRAY_TASK_ID % 100) + 1 ))

export DESURV_SIM_SCENARIO=${SCENARIOS[$scenario_idx]}
export DESURV_SIM_ANALYSIS=${ANALYSES[$analysis_idx]}
export DESURV_SIM_REPLICATE=$replicate

# rev 3 preamble: R 4.5.0, rev 3 library, DeSurv >= 1.2.0, isolated output dirs,
# DESURV_RECOMPUTE=TRUE, DESURV_NCORES from the allocation. Simulations are rerun in
# rev 3 only if their hyperparameters come from the BO (owner stop, plan section 5).
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

echo "=== Simulation task $SLURM_ARRAY_TASK_ID ==="
echo "Scenario: $DESURV_SIM_SCENARIO  Analysis: $DESURV_SIM_ANALYSIS  Replicate: $DESURV_SIM_REPLICATE"
echo "Cores: $SLURM_CPUS_PER_TASK  Node: $SLURM_NODELIST  Start: $(date)"
echo "========================================="

Rscript code/07_simulations.R

echo "=== Task $SLURM_ARRAY_TASK_ID finished: $(date) ==="
