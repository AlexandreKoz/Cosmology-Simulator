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
cmake --build --preset build-cpu-release --target bench_config_parser bench_profiling_overhead
./build/cpu-only-release/bench_config_parser
./build/cpu-only-release/bench_profiling_overhead
```

## Reporting expectations

Each benchmark report should include:

- build type/preset,
- hardware summary and thread count,
- enabled feature flags,
- setup vs steady-state timing where practical,
- throughput/rate metrics or effective bandwidth proxies for memory-sensitive paths.

Do not claim scientific correctness from benchmark throughput.

## Profiling discipline

- Keep hot-path instrumentation low-overhead and explicitly gated.
- Preserve deterministic behavior when profiling is disabled.
- Compare against prior baselines using the same preset and workload size.

## Workflow hooks for documentation/scaffolding changes

Documentation changes still need auditable developer workflow checks:

- integration doc-scaffold smoke test (`integration_docs_scaffold`)
- docs reference scan benchmark (`bench_docs_reference_scan`)

These hooks ensure core docs remain present and internally referenced as the codebase evolves.
