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
export DESURV_NCORES=$SLURM_CPUS_PER_TASK

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR

echo "=== Simulation task $SLURM_ARRAY_TASK_ID ==="
echo "Scenario: $DESURV_SIM_SCENARIO  Analysis: $DESURV_SIM_ANALYSIS  Replicate: $DESURV_SIM_REPLICATE"
echo "Cores: $SLURM_CPUS_PER_TASK  Node: $SLURM_NODELIST  Start: $(date)"
echo "========================================="

Rscript code/07_simulations.R

echo "=== Task $SLURM_ARRAY_TASK_ID finished: $(date) ==="
