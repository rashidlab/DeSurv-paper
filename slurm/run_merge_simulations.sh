#!/bin/bash
#SBATCH --job-name=desurv-merge
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=slurm/logs/merge_%j.log

# Merges partial simulation results and compiles the manuscript.
# Should run after both the main pipeline and simulation array complete.

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR

echo "=== Merging simulation results ==="
echo "Start: $(date)"

Rscript code/07b_merge_simulations.R

echo "=== Compiling manuscript ==="
Rscript code/09_render_paper.R

echo "=== Merge + paper complete: $(date) ==="
