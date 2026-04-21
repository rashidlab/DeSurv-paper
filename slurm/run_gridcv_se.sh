#!/bin/bash
#SBATCH --job-name=desurv-gridcv-se
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=slurm/logs/gridcv_se_%j.log

# Regenerates results/precomputed/gridcv/desurv_grid_history_*.rds
# (and dependent desurv_bo_results_*, tar_k_selection_*, tar_params_best_*)
# with the new se_cindex column.
#
# Before submitting, delete:
#   results/precomputed/gridcv/desurv_grid_history_{,_alpha0,_elbowk}_tcgacptac.rds
#   results/precomputed/gridcv/desurv_bo_results_{,_alpha0,_elbowk}_tcgacptac.rds
#   results/precomputed/gridcv/tar_k_selection_tcgacptac.rds
#   results/precomputed/gridcv/tar_params_best_{,_alpha0,_elbowk}_tcgacptac.rds
#
# fit_std_tcgacptac.rds, std_nmf_selected_k_tcgacptac.rds, and
# tar_data_filtered_*.rds are reused from cache.

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv Grid CV rerun (adds se_cindex) ==="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "Cores:    $NCORES"
echo "Start:    $(date)"
echo "=============================================="

export DESURV_NCORES=$NCORES
# DESURV_RECOMPUTE intentionally NOT set — relies on deleted files to
# trigger recompute and keeps fit_std / tar_data_filtered cached.

Rscript code/03_gridcv.R || { echo "FAILED at step 3"; exit 1; }

echo "=== done: $(date) ==="
