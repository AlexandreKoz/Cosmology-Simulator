# Profiling and benchmark workflow

## OpenMP concurrency contract

`ProfilerSession` and `CounterRegistry` use per-thread shards. Worker threads
record scopes, events, byte counts, and counters without a global mutex on each
operation. The first use by a thread registers its shard under a short mutex;
report generation and reset are quiescent-boundary operations that merge shards
deterministically after worker activity has joined. `AllocatorStats` uses atomic
accounting. Callers must not request a report/reset while worker threads are
still mutating the same session.


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

## TreePM residual counters, timers, and geometry telemetry

The gravity step event derived from `TreePmDiagnostics` / `TreePmProfileEvent`
reports:

- **Pair counters.** `local_pair_evaluations` and
  `incoming_remote_pair_evaluations` are non-overlapping residual bundles.
  `remote_pair_evaluations` is a compatibility alias of the incoming-remote
  value; `total_pair_evaluations` is the residual total and equals the exact
  sum of the two local/incoming components for the same solve.
- **Wall timers.** `tree_wall_ms_recent` is short-range tree wall time alone;
  `pm_wall_ms_recent` is `pm_profile.total_ms` alone (PM phase entry through
  tree short-range start, including long-range refresh when taken).
  `incoming_request_decode_validation_ms`,
  `incoming_remote_target_compute_ms`, and
  `incoming_response_encode_pack_ms` split incoming request
  decode/validation/hash, pure incoming tree compute (force evaluation of
  validated targets against this rank's tree only), and response
  encode/pack respectively;
  `protocol_validation_ms`, `protocol_consensus_ms`, and
  `response_exchange_ms` split response layout/payload sizing, failure
  consensus, and the response `MPI_Neighbor_alltoallv` call respectively.
  `let_remote_traversal_ms` is a compatibility alias of
  `incoming_remote_target_compute_ms`. OpenMP provenance
   fields `openmp_compiled`, `openmp_configured_workers`, and
   `openmp_observed_workers` come from `TreePmDiagnostics`
   (configured is the planned configured/runtime maximum; observed is the
   maximum actual team size recorded across the local, distributed overlap,
   and incoming-target regions during the last short-range solve; serial
   execution reports one).

  `residual_local_target_count` / `residual_incoming_target_count` count
  targets evaluated in each bundle for that solve;
   `residual_worker_scratch_high_water_bytes` is the retained contiguous
   worker-stack capacity (`worker_count * (1 + 7 * kMaximumTreeDepth) *
   sizeof(TreeLocalIndex)`); worker counter bundles are `O(worker_count)` and
   deterministic block floating diagnostics are `O(ceil(active_targets/64))`.
   The matching MemoryGovernor estimates are
   `gravity.estimate.treepm_residual_worker_scratch` and
   `gravity.estimate.treepm_residual_block_diagnostics`.

- **LET high-water.** `let_wire_buffer_high_water_bytes` is the four reusable
  short-range payload buffers; `let_known_workspace_high_water_bytes` is the
   known workspace peak (wire + structured + counts/masks/accumulators +
   metadata + transient codec). CHUÍ-owned retained storage is read from

  actual `vector.capacity()`; the transient codec term remains a modeled
  conservative upper envelope. These are capacity high-waters, not bytes
  communicated and not preflight estimates.
- **Domain geometry.** `domain_geometry_source_generation`,
  `current_gravity_source_generation`, `domain_geometry_fresh`,
  `domain_geometry_refreshed`, `domain_geometry_fallback_used`,
  `domain_geometry_fallback_reason` (string name of
  `TreePmDomainGeometryFallbackReason`),
  `domain_geometry_refreshed_leaf_count`,
  `domain_geometry_out_of_seed_range_source_count`, and
  `domain_geometry_uncovered_source_count` record the derived top-domain
  freshness/refit/fallback path. The published refit is seeded from the
  stable decomposition-local seed leaf set (not the previous published
  result), and an ownership commit invalidates freshness until reinstall.
  They are observational; restart/snapshot schema is unchanged.

## PM routing runtime events

The `gravity.treepm` runtime event keeps total communication traffic separate from
resident routing memory. Its PM payload includes logical density/force/potential
route counts and peer counts, total MPI bytes sent/received, MPI wait time, and
four distinct capacity high-waters: `pm_routed_send_buffer_high_water_bytes`,
`pm_routed_receive_buffer_high_water_bytes`,
`pm_routed_combined_buffer_high_water_bytes`, and
`pm_routed_workspace_high_water_bytes`. The combined value is the simultaneous
capacity of the two reusable wire buffers; the workspace value additionally
includes retained rank-scale routing count/displacement metadata and the reused
per-peer packing/response cursor. The capacity model intentionally leaves 64 KiB
below the 128 MiB/rank M1A engineering ceiling. These are capacity/high-water
metrics, not aliases for bytes communicated over the phase.

The DMO process preflight also reports hierarchical scheduler memory as three
separate values: logical live bytes, retained owned capacity, and historical
retained-capacity high-water. Candidate labels remain transient compatibility
inputs and do not contribute population-scale scheduler storage.

## M1B memory-governor telemetry

The existing profiler `memory_report` JSON object may contain a `governor`
object. It reports the authoritative process hard limit, reconciled baseline
owned bytes, configured external/planned reserves, governed committed bytes,
pending reserved bytes, raw and safety-adjusted accounted demand, raw headroom,
pressure (`green`, `amber`, `red`, or `trip`), historical committed/reserved/
accounted high-waters, and the reservation rejection count. This is an additive
view in the existing profiler format, not a second profiling stream.

Do not add `memory_report.totals` to governor `committed_bytes`: governed
physical blocks remain visible as ordinary `MemoryEntry` capacity for ownership
reporting, but are excluded from `baseline_owned_bytes` using the
`governed_commitment` reconciliation marker. The governor snapshot is the policy
view; the ordinary memory entries remain the ownership/capacity view. M1B does
not claim RSS/PSS reconciliation.

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
- `logical_metadata_bytes_read`
- `hash_bytes_read`
- `logical_payload_bytes_read`
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
file has one stable payload reader/session. In the default verified-identity mode, complete-file
SHA-256 work is one inspection pass per source file plus stable file-identity validation around
ingestion; explicit strict-full-rehash mode performs one additional completion SHA-256 pass.
The hash-pass count is independent of batch count, each source chunk is assigned once, reader-record imbalance is reported as the maximum minus minimum assigned-record count across ranks, main exchanges
scale with routing batches rather than source chunks, each source ID balances
against one final ID, and no rank allocates authoritative arrays sized to the
global particle count merely because MPI is enabled. `main_exchange_count` and the compatibility
`routing_collective_phase_count` are global protocol counters recorded on rank
zero; byte counters remain rank-local and may be reduced by the caller. The
compatibility collective fields count logical rank-consistent protocol phases,
not raw MPI calls. Actual production communicator calls are counted centrally by
the `mpi_*_call_count` fields. Their per-kind sum equals
`mpi_collective_call_count`; `routing_mpi_collective_call_count` is exactly 20
calls per successful routing batch in routing protocol version 3 (12 consensus votes,
three coverage reductions, two `Alltoall`/`Alltoallv` pairs, and one exact
reconciliation reduction). Fixed discovery, manifest, final audit, and
finalization calls are reported by `nonrouting_mpi_collective_call_count`.
For routing protocol version 3 the successful-path identity is:

```text
routing_mpi_collective_call_count
  = 20 * routing_batch_count

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


## Star-formation and effective-ISM benchmarks

`bench_star_formation_spawn` reports no-birth, sparse-birth, and dense plan/append throughput plus total allocations, allocated bytes, and temporary-memory upper bounds. Exact ID precommit uses sorted contiguous batches, avoiding one heap node per birth. `bench_effective_multiphase_ism` separately reports table initialization, lookup, direct-equilibrium evaluation, and hydro-closure throughput. Table lookup reconstructs neither cooling tables nor EOS state and allocates nothing per cell.

## M1C-1 process-memory reconciliation fields

When the runtime memory report is attached, profiling emits a `process_memory`
object with `current_declared_residency_bytes`, the compatibility alias
`known_accounted_bytes`, optional current RSS, optional peak RSS, optional PSS,
non-negative `unexplained_resident_bytes`, and an optional observed/known ratio.
Current declared residency excludes inactive future output/restart policy
headroom and uncommitted reservations; governor `accounted_bytes` remains the
separate conservative admission/policy quantity. Unavailable OS measurements
are JSON `null`, not zero.

A `distributed_process_memory` object records rank count and per-metric local,
global-sum, rank-max, rank-mean, and max/mean imbalance for governor-accounted
demand, current RSS, peak RSS, and communication high-water. OS-derived metrics
are valid only when all ranks provide a value. These fields are observational;
the deterministic memory governor remains the allocation authority.

## M2D runtime topology telemetry

The reference workflow emits a `runtime.topology` profiler event containing MPI
world size/rank, MPI node-local size/rank where available, OpenMP compile
status, requested thread count, and configured thread count. These fields are
intended for rank/thread topology comparisons and memory/work imbalance
interpretation. They are diagnostic evidence only and do not alter numerical
configuration or restart state.

Optional science-diagnostic pressure handling emits
`analysis.memory_pressure_deferral` with current pressure plus light/heavy due
and pending state. When a previously deferred diagnostic later executes on a
non-Red analysis hook outside its original cadence, the profiler emits
`analysis.memory_pressure_catchup` with light/heavy catch-up flags. These events
make bounded deferral and starvation recovery auditable without promoting the
pending bits into restart or scientific state.

## M2D-1 optional cadence evidence

`analysis.memory_pressure_deferral` includes pending state and missed counts. `analysis.memory_pressure_catchup` records class, first/latest due step, missed/coalesced count, actual execution step, and `historical_state_replayed=false`. `analysis.optional_cadence_dropped` records pending work discarded at segment termination; `analysis.optional_cadence_restart_policy` records the previous checkpoint's cadence summary and the nonpersistent restart policy. These are operational provenance events, not scientific outputs or an exact-cadence guarantee.

## M2D diagnostic admission evidence

The existing governor snapshot/rejection counters and `RuntimeEvent` stream
are the only memory-reporting authority. `analysis.memory_pressure_deferral`
records class, requested bytes, headroom and reason; catch-up/terminal events
retain actual versus due epochs. The model distinguishes FFT mesh coexistence
from reduction workspace, and the existing process RSS/PSS report remains the
measurement authority. Compare identical numerical configurations and thread
counts; do not equate a declared reservation with measured RSS. A required
health preflight is owner-managed and does not constitute a second physical
allocation.

### PM spectral-operator cache telemetry

`PmProfileEvent::spectral_operator_rebuilds` counts reconstructions of the
scale-free periodic PM Poisson/gradient operator arrays. Cosmological
`scale_factor` remains part of PM collective-entry consensus and force-state
coordination, but it is not an invariant spectral-operator cache dependency.
Changing only `scale_factor` therefore must leave this counter at zero after an
operator has already been built for the same mesh, box, split scale, gravity
constant, assignment scheme, and deconvolution policy.

## MPI global diagnostics ownership

Production `AnalysisRuntime` distinguishes rank-local operational state from shared science/run-health diagnostics. A due shared diagnostic is first prepared locally on every rank; local preparation failures are coordinated before any diagnostic reduction. Global health counters and SFR bins use communicator sums, invariant booleans use communicator-wide AND semantics with failing-rank counts, and angular-momentum **vectors** are summed component-wise before any norm is reported. Slice quicklooks reduce density sums plus sample counts before forming the global average; projection grids are element-wise sums.

Only rank 0 publishes ordinary shared diagnostic JSON/CSV files after those reductions complete, using `.part` then rename transactional publication. Memory reporting preserves publisher-rank local ownership while also exposing distributed rank-sum, rank-max, and imbalance fields; RSS/local ownership is never relabeled as one global process value. Large particle state is not gathered to the publisher.

A correct MPI power spectrum requires a global density field and global/distributed FFT. Until that contract exists, distributed diagnostics record `unsupported_under_mpi_requires_global_density_fft` and publish no averaged rank-local spectrum. Serial power-spectrum behavior is unchanged. Floating reductions retain the repository tolerance-equivalence contract rather than promising bitwise rank-count invariance; integer counts and boolean outcomes are exact.
