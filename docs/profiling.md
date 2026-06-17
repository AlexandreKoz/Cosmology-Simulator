# Profiling and benchmark workflow

CosmoSim benchmarks are lightweight hooks for performance trend visibility, not proofs of correctness.

## Benchmark directory structure

- Core/system hooks: `bench/bench_*.cpp`
- Module-focused hooks: `bench/gravity`, `bench/hydro`, `bench/amr`, `bench/io`, `bench/mini`
- Reporting utility: `bench/reporting/bench_report.hpp`
- Baseline size presets: `bench/baselines/benchmark_sizes_v1.txt`

## Build and run

```bash
cmake --preset cpu-only-release
cmake --build --preset build-cpu-release --target bench_config_parser bench_profiling_overhead bench_pm_solver bench_tree_pm_coupling bench_parallel_decomposition_exchange
./build/cpu-only-release/bench_config_parser
./build/cpu-only-release/bench_profiling_overhead
./build/cpu-only-release/bench_pm_solver
./build/cpu-only-release/bench_tree_pm_coupling
./build/cpu-only-release/bench_parallel_decomposition_exchange
```

## Reporting expectations

Each benchmark report should include:

- build type/preset,
- hardware summary and thread count,
- enabled feature flags,
- setup vs steady-state timing where practical,
- throughput/rate metrics or effective bandwidth proxies for memory-sensitive paths.
- cache/build reuse counters for persistent infrastructure where applicable (for example PM FFT plan cache build count, warmup vs measured solve counts).

Do not claim scientific correctness from benchmark throughput.

## Distributed TreePM phase-2 performance hardening notes

- `bench_pm_solver` now reports `plan_cache_size` and `plan_build_count`; in a stable slab layout `plan_build_count` should stay flat across measured solves.
- `bench_tree_pm_coupling` runs repeated measured iterations and exercises PM cadence reuse (`refresh_long_range_field=false` on selected iterations) so PM reuse impacts are visible in timing deltas.
- `bench_parallel_decomposition_exchange` remains the communication planning baseline for ownership/decomposition costs; compare alongside TreePM coupling output when evaluating communication overhead.

## Profiling discipline

- Keep hot-path instrumentation low-overhead and explicitly gated.
- Preserve deterministic behavior when profiling is disabled.
- Compare against prior baselines using the same preset and workload size.

## Hydro runtime events

The reference workflow emits `hydro.conservation` events from the Godunov hydro stage when a profiler session is
present. Payload values are volume-integrated local diagnostics over the cells updated by that stage. The event includes
before/after mass, flux/source/floor mass deltas, residuals for mass, momentum x/y/z, total energy, and derived internal
energy, total-energy source/floor deltas, `internal_energy_floor_count`, and the tolerances used by the CPU closed-box
regression. Source-term deltas are reported separately from flux deltas so gravity/expansion work is not classified as
a flux-conservation error.

The solver also records these same totals in `HydroProfileEvent::conservation`, alongside the existing face count,
fallback counters, bytes moved, and stage timings.

For production AMR hydro, the reference workflow emits `hydro.amr_production_stage` after the AMR hydro synchronization
point. The payload includes patch and active cell/face counts, `flux_register_entry_count`, reflux corrected cell count,
corrected mass, corrected momentum x/y/z, corrected total energy, corrected internal energy, complete register count,
and skipped register counts for incomplete, area-mismatched, or missing-target registers. These diagnostics are
observational only; flux-register ownership remains in the AMR hydro synchronization path and is not persisted as
restart truth.

## Workflow hooks for documentation/scaffolding changes

Documentation changes still need auditable developer workflow checks:

- integration doc-scaffold smoke test (`integration_docs_scaffold`)
- docs reference scan benchmark (`bench_docs_reference_scan`)

These hooks ensure core docs remain present and internally referenced as the codebase evolves.
