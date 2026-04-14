# Repair open issues (P01–P19 freeze ledger)

_Date captured: 2026-04-07 (UTC)_

## Current blocker ledger after dependency-enabled validation

No open P19 blockers remain after the PM/FFTW stabilization repair and full preset revalidation.

## Recently closed items

| ID | Status | Area | Prior symptom | Closure evidence |
|---|---|---|---|---|
| P19-STATE-OWNERSHIP-005 | Closed | Core state ownership/invariants decomposition | `src/core/simulation_state.cpp` concentrated ownership rules, species indexing, transfer packing, metadata parsing, and active-view hot/cold logic in one implementation unit. | Split into focused implementation units; explicit gravity/hydro hot-field contract notes added; targeted invariant/contract checks added in `tests/unit/test_simulation_state.cpp`; verified by `cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state` and `./build/cpu-only-debug/test_unit_simulation_state`. |
| P19-IO-CONTRACT-006 | Closed | Snapshot vs restart contract separation | Restart continuation-critical metadata/scheduler requirements and failure behavior were under-specified relative to snapshot interoperability boundaries. | Added shared continuation metadata contract helper, explicit restart completeness checklist, and negative tests for schema mismatch/missing required scheduler lane/finalize failure behavior; verified via `test_unit_snapshot_hdf5_schema`, `test_unit_restart_checkpoint_schema`, and `test_integration_restart_checkpoint_roundtrip`. |
| P19-GATE-FFTW-TEST-001 | Closed | PM FFTW analytic mode response | `unit_pm_solver` failed at `tests/unit/test_pm_solver.cpp:82` (`cosine_similarity > 0.98`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"` now passes. |
| P19-GATE-FFTW-TEST-002 | Closed | TreePM periodic PM long-range coupling (FFTW path) | `integration_tree_pm_coupling_periodic` failed with `rel_l2=18129.9` (required `<= 0.75`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` now passes (`36/36`). |
| P19-ARCH-CORE-BOUNDARY-003 | Closed | Core dependency direction | `include/cosmosim/core/reference_workflow.hpp` assembled analysis/I/O/physics workflow concerns inside `core/`. | Workflow assembly moved to `workflows/` and `integration_core_dependency_direction` guard added to fail on forbidden upward includes. |
| P19-CONFIG-CONTRACT-004 | Closed | Frozen configuration contract typing | Policy keys (`gravity_solver`, `hydro_solver`, `coordinate_frame`, boundaries, feedback mode/variant) remained string contracts after parse. | `unit_config_parser`, `unit_simulation_mode`, and `integration_simulation_mode_toy_runs` now validate enum-backed freeze behavior with deterministic normalized text/hash and centralized key+alias registry handling. |
| P19-OBS-EVENTS-007 | Closed | Operational observability | Runtime diagnostics lacked a structured run-health/event surface linked to provenance; failures were harder to audit from CI artifacts. | Added compact runtime event model in `core::ProfilerSession`, new `reference_operational_events.json` report, and test coverage in `test_unit_profiling` + `test_integration_reference_workflow`. |
| P19-DIAGNOSTICS-TIER-008 | Closed | Diagnostics maturity/scalability labeling | Diagnostics bundles mixed run-health outputs with heavy reference algorithms and did not encode maturity tier in code/output policy. | Added typed analysis execution policy, per-diagnostic tier/maturity/scalability metadata in bundle JSON, and tests proving provisional heavy diagnostics are non-default unless explicitly opted in (`all_including_provisional`). |
| P19-PARALLEL-CONTRACT-009 | Closed | Distributed-memory ownership/determinism contracts | Owned-vs-ghost residency and reduction/config consensus expectations were partially implicit in ghost-owner vectors and ad-hoc checks. | Added typed ghost residency descriptors + validation, deterministic reduction agreement helpers, rank config-freeze consensus checks, and targeted unit/integration pseudo-multi-rank tests in `test_parallel_distributed_memory` and `test_parallel_two_rank_restart`. |

## Verified non-blocking evidence (this run)

- PM HDF5+FFTW path passes:
  - `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"`
  - `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` (`36/36` passed)
- HDF5 path passes:
  - `ctest --preset test-hdf5-debug --output-on-failure` (`36/36` passed)
- CPU-only path passes:
  - `ctest --preset test-cpu-debug --output-on-failure` (`36/36` passed)

## Outcome tags

- **CPU-only preserved**
- **HDF5 path preserved**
- **PM HDF5+FFTW path restored**
- **P19 blocker ledger cleared**
