#!/bin/bash
set -euo pipefail

array_job_id=$(sbatch --parsable run_array.sh)
echo "array_job_id=${array_job_id}"

manifest_job_id=$(sbatch --parsable --dependency=afterany:${array_job_id} build_manifest.sh)
echo "manifest_job_id=${manifest_job_id}"
