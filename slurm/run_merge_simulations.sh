#!/bin/bash
#SBATCH --job-name=desurv-merge
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=slurm/logs/merge_%j.log

# Runs after main pipeline, simulation array, and CV grid all complete.
# Merges simulation results, generates all SI figures (including CV grid figure),
# generates Fig 2, and renders the paper.
# Use slurm/submit_all.sh to set up job dependencies automatically.

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR

echo "=== Merging simulation results ==="
echo "Start: $(date)"

Rscript code/07b_merge_simulations.R

echo "=== Generating SI figures ==="
Rscript code/09b_si_figures.R

echo "=== Generating simulation figures ==="
Rscript code/09c_sim_figures.R

echo "=== Compiling manuscript ==="
make paper

echo "=== Merge + paper complete: $(date) ==="
