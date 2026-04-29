# Stage 0 runtime-truth freeze gate (P0-10)

_Date: 2026-04-26 (UTC)_

## Stage 0 status

**Status: CLOSURE CANDIDATE AFTER P0-09/P0-10 REPAIR; full local certification still required.**

This repair pass removes the documented P0-09/P0-10 closure blockers found in the prior gate. Full Stage 0 closure still requires successful local completion of the listed CPU, HDF5, and PM+HDF5+FFTW preset suites on a normal development machine.

## Scope and ownership freeze summary

Stage 0 consolidation reviewed and retained the runtime ownership authorities documented in:

- `docs/architecture/runtime_truth_map.md`
- `docs/architecture/adr_runtime_truth_ownership.md`

No new ownership authority lanes were introduced in this gate pass. Runtime-truth owners remain discoverable in those two architecture artifacts.

## Ownership ADRs and gate references

- ADR-INFRA-AUDIT-012 — Stage 0 runtime-truth freeze initiated.
- ADR-INFRA-OWNERSHIP-013 — single-source-of-truth runtime ownership policy.
- ADR-INFRA-STAGE0-GATE-014 — this Stage 0 consolidation gate result (see `docs/architecture/decision_log.md`).

## Stage 0 required domains checklist

The following domains have documented source-of-truth owner, mutating authority, mirror/cache/view policy, restart policy, and reorder/resize/migration policy (where applicable) in the ownership ADR + runtime truth map:

- timestep/bin authority + scheduler active-set construction/invalidation
- particle ordering/resize/reorder sidecar movement
- species migration/grouping
- gas-cell identity mapping
- softening override/default priority and persistence
- config raw/normalized/derived/runtime/provenance ownership
- restart/reload round-trip continuation contracts

## Test evidence (commands run in this gate)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "runtime|truth|timestep|scheduler|active|particle|sidecar|species|gas|hydro|softening|config|restart|checkpoint|provenance"
ctest --preset test-cpu-debug --output-on-failure

cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure

cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

### Result summary

- Stage 0 targeted runtime-truth suite on `cpu-only-debug`: **pass**.
- Full CPU debug suite: **fail/not fully complete in this environment** due:
  - `unit_snapshot_hdf5_schema` assertion mismatch (`gadget_arepo_v2` expected),
  - extremely long `integration_tree_pm_coupling_periodic` runtime in this environment preventing timely complete-cycle closure.
- HDF5 debug suite: **fail/not fully complete in this environment** due:
  - `integration_softening_ownership_invariants` restart metadata failure (`requires provenance.config_hash_hex`),
  - `unit_snapshot_hdf5_schema` schema-name mismatch,
  - same long-running `integration_tree_pm_coupling_periodic` timing blocker.
- PM+HDF5+FFTW preset suite: **fails fast with explicit blockers** (8 failing tests), including:
  - `integration_softening_ownership_invariants`
  - `unit_pm_solver`
  - `unit_snapshot_hdf5_schema`
  - `integration_tree_pm_coupling_periodic`
  - `integration_docs_scaffold`
  - `integration_release_readiness_artifacts`
  - `integration_runtime_app_smoke`
  - `validation_phase2_mpi_gravity_single_rank`

## Invariants enforced by this gate pass

- Stage 0 runtime-truth regex suite successfully re-ran across timestep/bin, active-set, particle/sidecar, species, gas, softening, config/provenance, and restart/checkpoint domains on CPU debug.
- No invariant tests were weakened or removed.
- Added a dedicated grouped test preset:
  - `test-stage0-runtime-truth-cpu-debug`
  to keep Stage 0 runtime-truth tests runnable as one named suite.

## Remaining repair tickets / blockers

See `docs/repair_open_issues.md` for explicit ticketed blockers recorded from this gate pass (including Stage 0 closure blockers).

## Physics-scope confirmation

This gate introduced no new physics models, no new solver/numerics features, and no intended solver-behavior changes. Changes are limited to test grouping and repair-state/decision documentation.

## Reproducibility impact

Determinism/reproducibility policy is unchanged. This gate pass only records command-backed evidence and closure status.

## Known limitations

- Full CPU/HDF5 suite completion in this container is constrained by the runtime cost of `integration_tree_pm_coupling_periodic`.
- PM+HDF5+FFTW path remains red due explicit failing tests listed above.

## Blocked future work and permitted next stage

- **Blocked:** new physics/performance feature work while Stage 0 remains not closed.
- **Permitted next-stage work:** targeted infrastructure repairs for Stage 0 blockers, then re-run full required preset commands and re-evaluate closure.

## Commit context (at gate run)

- Evaluated commit: `dc989e07111786621e4dc8e1e095f509fdd13c1f`

## Follow-up repair pass — P0-00 through P0-04 hardening

This follow-up pass addresses the first audit batch only: runtime-truth documentation, timestep/bin mirror authority, particle ordering/sidecar invariants, species migration invariants, and active-view lifetime safety. It does not close later Stage 0 restart/HDF5 gates.

Changes made:

- Added `syncTimeBinMirrorsFromScheduler(...)` with an explicit `TimeBinMirrorDomain` so timestep/bin mirror refreshes use a named authoritative synchronization path instead of ad-hoc local copy loops.
- Replaced pointer-keyed global active-view generation registries with generation fields carried directly by mutable kernel views (`GravityParticleKernelView::source_particle_index_generation` and `HydroCellKernelView::source_cell_index_generation`). Scatter now validates the view against the owning `SimulationState` generation without hidden process-global pointer maps.
- Added an explicit `ParticleSidecar::has_gravity_softening_override` mask plus helper methods so materialized species/default softening values are not confused with true per-particle override authority.
- Updated particle reorder and migration paths to preserve both softening values and override-mask semantics across reorder, resize, and species migration.
- Added/updated invariant coverage so species migration preserves materialized softening values without promoting them to override authority, true override flags survive migration/reorder, and migration records distinguish numeric softening-value presence from override authority.

Validation note:

- Static edits were completed in the repository tree. The execution environment used for this repair pass repeatedly timed out during CMake/compiler invocations before producing actionable diagnostics, so this archive should be rebuilt with the standard presets in a normal development environment.
- Required local verification commands remain:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "time|scheduler|active|particle|sidecar|species|migration|reorder"
ctest --preset test-cpu-debug --output-on-failure
```

## Follow-up repair pass — P0-09 through P0-10 restart/gate closure hardening

This pass targets the remaining Stage 0 restart/reload and consolidation-gate blockers only. It does not add new physics, solvers, or performance features.

Changes made:

- Restored the missing `tests/unit` and `tests/validation` source files referenced by `CMakeLists.txt`, eliminating a configure-time test-registry/source-tree drift blocker in the uploaded archive.
- Repaired `integration_softening_ownership_invariants` restart payload construction so the restart writer receives a valid normalized config hash, matching `ProvenanceRecord::config_hash_hex`, and a valid single-rank `DistributedRestartState` instead of tripping the continuation-metadata guard.
- Hardened the legacy softening compatibility round-trip by clearing both the optional softening value lane and the authoritative override-mask lane when constructing the legacy payload.
- Aligned docs, release manifest, and gate checks with the current runtime schemas: snapshot `gadget_arepo_v4`, restart `cosmosim_restart_v5`, and provenance `provenance_v4`.
- Updated the runtime app smoke gate to check the canonical axis-aware normalized key `treepm_pm_grid_nx = 24` rather than the deprecated scalar alias `treepm_pm_grid = 24`.

Validation evidence from this sandbox:

- `cmake --preset cpu-only-debug` completed successfully after restoring the missing test sources.
- `cmake --preset hdf5-debug` reached the normal configure/generate phase with HDF5 detected (`1.14.5`), but the sandbox command timed out after configure output and before a complete build/test cycle could be certified.
- Targeted build attempts for schema tests began compiling normally but exceeded the sandbox wall-time limit without producing compiler diagnostics.

Required local certification commands remain:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "runtime|truth|timestep|scheduler|active|particle|sidecar|species|gas|hydro|softening|config|restart|checkpoint|provenance"
ctest --preset test-cpu-debug --output-on-failure

cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure

cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```
