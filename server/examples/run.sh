#!/bin/bash
#SBATCH --job-name=CpH
# normal-32g or normal-64g
#SBATCH --partition=normal-64g
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28gb

module swap gnu12 intel/2023.2.1

srun ./ANNT.exe $SLURM_ARRAY_TASK_ID
