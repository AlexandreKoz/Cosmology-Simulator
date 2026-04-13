# Initial release benchmark summary (v0.1.0-initial)

Benchmarks provide performance expectations and regression visibility. They are not correctness evidence.

## Environment disclosure template (required for release reports)

- Build preset: `cpu-only-release`
- Build type: `Release`
- Host CPU model: `<fill from run host>`
- Logical threads used: `<fill>`
- MPI ranks used: `<fill>`
- Enabled features (`mpi/hdf5/fftw/cuda/python`): `<fill from build metadata>`

## Required release benchmark hooks

- `bench_reference_workflow`
- `bench_mini_run_pipeline`
- `bench_config_parser`
- `bench_docs_reference_scan`

## Workstation-class expectation envelope (guidance)

These are intentionally conservative and used for trend checks:

- `bench_reference_workflow`: should complete in seconds on modern desktop CPUs for default mini workload.
- `bench_mini_run_pipeline`: should remain within a single-digit second envelope for baseline sizes.
- `bench_config_parser`: should remain sub-millisecond for canonical example configs.
- `bench_docs_reference_scan`: should remain sub-100ms on SSD-backed developer environments.

## Reporting guidance

- Separate setup and steady-state timing when the benchmark supports it.
- Include at least one throughput or rate metric where available.
- Keep baseline comparisons preset-stable (`cpu-only-release` vs previous `cpu-only-release`).
