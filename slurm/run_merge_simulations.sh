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

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

echo "=== Merging simulation results ==="
echo "Start: $(date)"

Rscript code/07b_merge_simulations.R

echo "=== Generating SI figures ==="
Rscript code/09b_si_figures.R

echo "=== Generating simulation figures ==="
Rscript code/09c_sim_figures.R

# Rendering writes the tracked paper/paper.pdf and paper/si_appendix.pdf, needs LaTeX and
# every static figure inside DESURV_FIGURE_DIR, so rev 3 does not render here by default.
if [ "${REV3_RENDER:-0}" = 1 ]; then
  echo "=== Compiling manuscript ==="
  DESURV_RECOMPUTE=FALSE Rscript code/10_render_paper.R
else
  echo "=== Render skipped (set REV3_RENDER=1 to render from DESURV_RESULTS_DIR) ==="
fi

echo "=== Merge + paper complete: $(date) ==="
