# Benchmarks

Benchmark targets follow `bench_<topic>` naming and remain independently runnable.

## Benchmark suite layout

- `bench/gravity/bench_gravity_kernels.cpp`: Tree gravity and PM kernel timings with stable particle/grid sizes.
- `bench/hydro/bench_hydro_kernels.cpp`: finite-volume face sweep kernel with scratch/cache reuse.
- `bench/amr/bench_amr_regrid_kernel.cpp`: AMR regrid + reflux synchronization kernel costs.
- `bench/io/bench_io_restart_kernel.cpp`: restart checkpoint write/read latency and bandwidth proxies.
- `bench/mini/bench_mini_run_pipeline.cpp`: representative small end-to-end mini-run phase timing.
- `bench/reporting/bench_report.hpp`: common reporter and execution metadata helper.
- `bench/baselines/benchmark_sizes_v1.txt`: reproducibility contract for benchmark input sizes.

## Reporting conventions

All suite benchmarks emit single-line key-value records with:

- `build_type`, `threads`, `mpi_ranks`, and `device`.
- Warmup and measurement iteration counts.
- Wall-time breakdown for major phases.
- Throughput counters (e.g., updates/s, faces/s).
- `bytes_moved` or explicit `bytes_moved_proxy` with effective bandwidth proxy.

Use environment overrides when comparing scaling variants:

- `COSMOSIM_BENCH_THREADS`
- `COSMOSIM_BENCH_MPI_RANKS`
- `COSMOSIM_BENCH_DEVICE`
