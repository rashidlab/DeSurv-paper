#!/bin/bash
#SBATCH --job-name=desurv-cv-grid
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/cv_grid_%j.log

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

echo "=== DeSurv CV Grid Search ==="
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURM_NODELIST"
echo "Cores:     $SLURM_CPUS_PER_TASK"
echo "Start:     $(date)"
echo "Directory: $(pwd)"
echo "R version: $(Rscript --version 2>&1)"
echo "============================="

make cv-grid NCORES=$SLURM_CPUS_PER_TASK

echo "=== CV Grid finished: $(date) ==="
