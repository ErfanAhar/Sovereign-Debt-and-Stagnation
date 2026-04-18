#!/bin/bash
#SBATCH --job-name=cal_pre
#SBATCH --partition=normal-64g
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --output=slurm/precompile_%j.out
#SBATCH --error=slurm/precompile_%j.err

set -euo pipefail

module load julia

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${ROOT_DIR}"
mkdir -p slurm result/calibration_grid

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

julia --project="${ROOT_DIR}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
