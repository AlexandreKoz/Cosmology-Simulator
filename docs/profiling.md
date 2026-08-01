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

## Initial-condition ingestion events

After initial state construction, every rank emits one
`io.ic_ingestion.summary` event in subsystem `io.initial_conditions`. Its
payload records the single `provenance_authority` selected for the import and
includes:

- `files_assigned`
- `chunks_assigned`
- `metadata_bytes_read`
- `hash_bytes_read`
- `payload_bytes_read`
- `converted_payload_bytes`
- `bytes_serialized`
- `bytes_sent`
- `bytes_received`
- `manifest_metadata_bytes_communicated`
- `records_read`, `records_converted`, and `records_routed`
- `source_file_open_count`, `source_dataset_open_count`
- `full_file_hash_pass_count`, `source_identity_validation_count`
- `routing_batch_count`, `reader_batches_assigned`, `reader_records_assigned`, `reader_record_imbalance`
- `main_exchange_count`, `exact_audit_exchange_count`, `distributed_id_audit_round_count`
- `logical_consensus_phase_count`, `routing_logical_consensus_phase_count`
- compatibility aliases `collective_phase_count`, `routing_collective_phase_count`
  (logical phases only)
- `mpi_collective_call_count`, `routing_mpi_collective_call_count`,
  `nonrouting_mpi_collective_call_count`
- `mpi_allreduce_call_count`, `mpi_bcast_call_count`,
  `mpi_gather_call_count`, `mpi_gatherv_call_count`,
  `mpi_alltoall_call_count`, `mpi_alltoallv_call_count`
- `collectives_per_million_records`
- `wall_time_nanoseconds`
- `peak_staging_bytes`
- final local particle, gas, star, black-hole, and tracer counts
- `already_partitioned`

`bytes_read` is retained as the checked sum of metadata, SHA-256, and particle
payload reads. Metadata is the decoded logical HDF5 header payload; hashing
counts every source byte read by the assigned hashing rank; payload counts only
datasets actually read, never merely inspected or explicitly dropped fields.

Peak staging is a capacity-based high-water mark for the actual bounded
import/routing workspace, not the final authoritative owner-local state. It
includes all simultaneously live coordinate, velocity, mass, ID, gas, star,
black-hole, tracer, `ParticleRecord`, nested per-peer, flattened exchange,
count/displacement, decode, coverage, and ID-reconciliation buffers. The exact global duplicate-ID audit
uses rank-local sorted temporary runs and a bounded-memory external merge; disk
run bytes are not RAM staging and are removed before import returns. For a distributed fixture, the required evidence is that each nonempty source
file has one stable payload reader/session, complete-file SHA-256 work is bounded
by three passes per source file (inspection, session start, session completion)
independent of batch count, each source chunk is assigned once, reader-record imbalance is reported as the maximum minus minimum assigned-record count across ranks, main exchanges
scale with routing batches rather than source chunks, each source ID balances
against one final ID, and no rank allocates authoritative arrays sized to the
global particle count merely because MPI is enabled. `main_exchange_count` and the compatibility
`routing_collective_phase_count` are global protocol counters recorded on rank
zero; byte counters remain rank-local and may be reduced by the caller. The
compatibility collective fields count logical rank-consistent protocol phases,
not raw MPI calls. Actual production communicator calls are counted centrally by
the `mpi_*_call_count` fields. Their per-kind sum equals
`mpi_collective_call_count`; `routing_mpi_collective_call_count` is exactly 30
calls per successful routing batch in protocol version 1 (22 consensus votes,
three coverage reductions, two `Alltoall`/`Alltoallv` pairs, and one exact
reconciliation reduction). Fixed discovery, manifest, final audit, and
finalization calls are reported by `nonrouting_mpi_collective_call_count`.
For protocol version 1 the successful-path identity is:

```text
routing_mpi_collective_call_count
  = 30 * routing_batch_count

nonrouting_mpi_collective_call_count
  = 40
  + (validate_runtime_cosmology ? 1 : 0)
  + source_file_count
  + 10 * distributed_id_audit_round_count
  + mpi_bcast_call_count
```

The `mpi_bcast_call_count` term is explicit because length-prefixed metadata
broadcasts may require more than one 64 MiB payload chunk. The distributed MPI
acceptance test checks both identities rather than applying an arbitrary loose
ceiling. `collectives_per_million_records` uses the globally routed record
count and is zero when no records are routed.

These counters are scalability evidence, not a substitute for scientific
validation. Exact distributed duplicate-ID, count, mass, ownership, provenance,
and sidecar checks run separately before the workflow accepts the state.

## Workflow hooks for documentation/scaffolding changes

Documentation changes still need auditable developer workflow checks:

- integration doc-scaffold smoke test (`integration_docs_scaffold`)
- docs reference scan benchmark (`bench_docs_reference_scan`)

These hooks ensure core docs remain present and internally referenced as the codebase evolves.

## Star-formation profiling

`bench_star_formation_spawn` profiles three explicit adaptive-model phases: `eligibility_no_births`, `eligibility_sparse_births`, and `dense_birth_plan_and_batch_append`. Each line reports `cells_per_second`, `births_per_second`, `bytes_allocated_per_step`, `allocation_count_per_step`, and `peak_temporary_memory_upper_bound_bytes`. The peak is an upper bound because the tracker includes storage retained by the appended authoritative particle and stellar sidecars as well as temporary planning storage. The implementation uses one particle resize, one star-sidecar resize, one species-index rebuild, no collective inside the cell loop, no full existing-particle-ID scan, and no allocation node per candidate or newborn ID. Compare allocation counts across cell-count sweeps: a fixed number of vector/sidecar growth allocations is acceptable, while allocation growth proportional to candidate count is a regression.
