#!/bin/bash
#SBATCH -p short-serial
#SBATCH --job-name=createSummary
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=00:10:00
#SBATCH --array=RANGES

module add jaspy

./createSum.sh ${SLURM_ARRAY_TASK_ID} FILEDIR USER

