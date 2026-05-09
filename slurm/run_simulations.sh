#!/bin/bash
#SBATCH --job-name=desurv-sim
#SBATCH --partition=general
#SBATCH --array=0-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=slurm/logs/sim_%A_%a.log

# Map array task ID to (scenario, analysis)
SCENARIOS=(R0_easy R0_easy R00_null R00_null R_mixed R_mixed)
ANALYSES=(bo_tune_ntop bo_tune_ntop_alpha0 bo_tune_ntop bo_tune_ntop_alpha0 bo_tune_ntop bo_tune_ntop_alpha0)

export DESURV_SIM_SCENARIO=${SCENARIOS[$SLURM_ARRAY_TASK_ID]}
export DESURV_SIM_ANALYSIS=${ANALYSES[$SLURM_ARRAY_TASK_ID]}
export DESURV_RECOMPUTE=TRUE
export DESURV_NCORES=$SLURM_CPUS_PER_TASK

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR

echo "=== Simulation array task $SLURM_ARRAY_TASK_ID ==="
echo "Scenario:  $DESURV_SIM_SCENARIO"
echo "Analysis:  $DESURV_SIM_ANALYSIS"
echo "Cores:     $SLURM_CPUS_PER_TASK"
echo "Node:      $SLURM_NODELIST"
echo "Start:     $(date)"
echo "========================================="

Rscript code/07_simulations.R

echo "=== Task $SLURM_ARRAY_TASK_ID finished: $(date) ==="
