#!/bin/bash
#SBATCH --job-name=qual-grid
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --array=1-180
#SBATCH --output=slurm/logs/qual_%x_%A_%a.log
# Stage 1 (prespec section 9.1): one grid/reference fit per array task (180 = S1, S2, S3 x replicates 1-5 x g 1-12,
# manifest order), the 30 starts in parallel on 30 cores. The job refuses to run unless the verifier controls passed
# and, for S3, the stage 1 realized-ratio check passed.
# Resources (timing run 5e55295): a 30-start fit at p = 3,000 on 30 cores takes about the slowest start, at most ~700 s
# (n = 218); 4 h leaves ample margin for oversubscription and recording.
# this script was written. Interim timing_run.log (15 cores, oversubscribed machine): 551 s (n = 175) and 848 s
# (n = 218) per 30-start fit; the values above are placeholders.
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/slurm/qual/qual_env.sh"
i=$(( ${SLURM_ARRAY_TASK_ID:?} - 1 ))
SCEN=(S1 S2 S3); s=${SCEN[$(( i / 60 ))]}; r=$(( (i % 60) / 12 + 1 )); g=$(( i % 12 + 1 ))
Rscript analysis/estimator-gate/qualification/qual_grid_job.R "$s" "$r" "$g"
