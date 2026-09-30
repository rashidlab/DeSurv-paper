#!/bin/bash
#SBATCH --job-name=desurv-gate3a
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work/users/n/u/nur2/desurv-rerun/logs/gate3a_%j.out
#SBATCH --error=/work/users/n/u/nur2/desurv-rerun/logs/gate3a_%j.err

# Track 3a convergence probe, fixed-R rule (section 3 amendment d120e56): four R = 50 ensembles,
# and the R = 100 stage in the same job only if R = 50 is borderline. Resubmitting after a kill resumes from the checkpointed run chunks
# in $GATE_OUT_DIR/runs. Submit from the directory holding gate_3a.R (see README.md).
set -euo pipefail

module load r/4.5.0

export R_LIBS=/work/users/n/u/nur2/desurv-rerun/lib/r4.5
export GATE_DATA_DIR=/work/users/n/u/nur2/desurv-rerun/data
export GATE_OUT_DIR=/work/users/n/u/nur2/desurv-rerun/gate3a
export GATE_DRY=FALSE
unset GATE_FORCE_BORDERLINE   # dry-run-only override; the script refuses it outside GATE_DRY=TRUE
# one BLAS thread per forked worker
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p /work/users/n/u/nur2/desurv-rerun/logs "$GATE_OUT_DIR"

echo "Job $SLURM_JOB_ID on $(hostname), $SLURM_CPUS_PER_TASK cpus, started $(date)"
Rscript "${SLURM_SUBMIT_DIR}/gate_3a.R"
echo "Finished $(date)"
