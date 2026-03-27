#!/bin/bash
#SBATCH --job-name=desurv-pipeline
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/logs/desurv_%j.log

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR

mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

echo "=== DeSurv full pipeline ==="
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURM_NODELIST"
echo "Cores:     $SLURM_CPUS_PER_TASK"
echo "Start:     $(date)"
echo "Directory: $(pwd)"
echo "R version: $(Rscript --version 2>&1)"
echo "==========================="

make all NCORES=$SLURM_CPUS_PER_TASK

echo "=== Pipeline finished: $(date) ==="
