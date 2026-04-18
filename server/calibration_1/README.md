# calibration_1

Self-contained Slurm bundle for running the 4000-point sunspot calibration on the server.

## Files

- `src/`: copied model source files
- `Project.toml`: Julia environment
- `precompile.sh`: one-time package instantiate/precompile job
- `run_array.sh`: 4000-point Slurm array job
- `build_manifest.sh`: rebuilds `result/calibration_grid/manifest.csv`
- `submit_all.sh`: submits the array job, then the manifest job with a dependency
- `run_one_case.jl`: one Slurm task = one calibration point

## Server Usage

From inside `calibration_1`:

```bash
sbatch precompile.sh
sbatch run_array.sh
```

Or submit the array and manifest together:

```bash
bash submit_all.sh
```

## Notes

- The array is `1-4000%100`: all 4000 points are scheduled, with at most 100 concurrent tasks.
- Each task writes one `.jld2` file into `result/calibration_grid/`.
- `manifest.csv` is built after the array finishes, not during the array run.
- `cpus-per-task=8` is a reasonable starting point because the solver already uses Julia threads internally.
- If the cluster is busy, lower `%100`. If the nodes are underused, try `%150` or `%200`.
