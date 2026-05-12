#!/bin/bash
#SBATCH --job-name=desurv-main
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=18:00:00
#SBATCH --output=slurm/logs/main_%j.log

# Runs steps 01-05, 08, 09a, 09b (BO, fitting, validation, main/SI figures).
# Simulations run in parallel via slurm/run_simulations.sh + run_merge_simulations.sh.
# Use slurm/submit_all.sh to coordinate the full workflow with job dependencies.

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

echo "=== DeSurv main pipeline ==="
echo "Job ID: $SLURM_JOB_ID  Node: $SLURM_NODELIST  Cores: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "============================"

make main NCORES=$SLURM_CPUS_PER_TASK

echo "=== Main pipeline finished: $(date) ==="
