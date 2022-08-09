#!/bin/bash
#SBATCH -p short-serial
#SBATCH --job-name=RUNCODE
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=TIME
#SBATCH --array=RANGES
#SBATCH --mem=16000

./RUNCODE.sh ${SLURM_ARRAY_TASK_ID} FILEDIR

