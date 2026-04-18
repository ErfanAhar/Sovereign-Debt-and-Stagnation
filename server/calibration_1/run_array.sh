#!/bin/bash
#SBATCH --job-name=calib
#SBATCH --partition=normal-64g
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --array=1-4000%150
#SBATCH --output=slurm/calib_%A_%a.out
#SBATCH --error=slurm/calib_%A_%a.err

set -euo pipefail

module load julia

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${ROOT_DIR}"
mkdir -p slurm result/calibration_grid

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun julia --project="${ROOT_DIR}" \
    "${ROOT_DIR}/run_one_case.jl" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${ROOT_DIR}/result/calibration_grid"
