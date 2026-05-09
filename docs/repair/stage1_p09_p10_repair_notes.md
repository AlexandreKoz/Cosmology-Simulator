# Stage 1 Prompt 09-10 Repair Notes

Date: 2026-05-08
Scope: AMR/moving-mesh gas-identity seam and runtime-truth CI gate consolidation.

## Prompt 09 - GasCellIdentityMap seam hardening

Implemented repairs and upgrades:

- Upgraded `core::GasCellIdentityMap` from a thin vector wrapper into a validated identity mapping seam.
- Added atomic lookup-table rebuild during `assign(...)` for stable `gas_cell_id` and transient `local_cell_row` lookups.
- Reserved `gas_cell_id == 0` as an uninitialized sentinel and reject it during assignment.
- Added a monotonic `generation()` counter so future production hydro/restart views can reject stale row mappings.
- Added dense local row coverage validation through `coversDenseLocalRows(...)` and `requireCoversDenseLocalRows(...)`.
- Added explicit lineage/ownership query helpers: `rowsForParentParticleId(...)` and `rowsForPatch(...)`.
- Added `buildGasCellNewToOldRowMap(...)`, which remaps by stable `gas_cell_id` and marks genuinely new cells with `kInvalidGasCellRow` so callers cannot silently inherit old hydro state.
- Fixed a duplicated `gas_cell_id.resize(count)` call in `GasCellSidecar::resize(...)`.
- Kept the seam out of production hydro/restart state; the current `requireParticleBoundGasCellContract(...)` remains the production gate until a deliberate schema migration is implemented.

Tests added/strengthened:

- `unit_gas_cell_identity_invariants` now validates generation behavior, nonzero ID enforcement, dense local row coverage, parent/patch queries, lookup helpers, and stable-ID row remapping.

## Prompt 10 - CI/runtime-truth gate hardening

Implemented repairs and upgrades:

- Added `stage1_runtime_truth_targets`, an aggregate CMake build target for the exact compiled executables used by the Stage 1 runtime-truth P0 gate.
- Upgraded `scripts/ci/run_preset_pipeline.sh` to honor `COSMOSIM_CI_BUILD_TARGETS`, allowing CI jobs to build a scoped target set instead of compiling every benchmark, validation executable, and unrelated test binary.
- Updated `scripts/ci/run_stage1_runtime_truth_gate.sh` to default to `stage1_runtime_truth_targets` while still allowing callers to override the target scope.
- Updated `scripts/ci/enforce_infra_gates.sh` so the manifest-level Stage 1 gate also uses the scoped runtime-truth build target.
- Updated CI documentation with the exact scoped build command and the target-scope semantics.
- Reduced `unit_gas_cell_identity_invariants` linkage from `cosmosim_all` to `cosmosim_core`, matching its actual dependency boundary and avoiding unnecessary gravity/hydro/PM compilation for this core identity test.

## Validation performed

Commands executed in the repair environment:

```bash
bash -n scripts/ci/run_preset_pipeline.sh scripts/ci/run_stage1_runtime_truth_gate.sh scripts/ci/enforce_infra_gates.sh
cmake --fresh --preset cpu-only-debug
cmake --build --preset build-cpu-debug --parallel 2 --target test_unit_gas_cell_identity_invariants
ctest --preset test-cpu-debug --output-on-failure -R '^unit_gas_cell_identity_invariants$'
ctest --preset test-cpu-debug --output-on-failure -R '^integration_runtime_truth_ctest_labels$'
COSMOSIM_CI_BUILD_TARGETS=test_unit_gas_cell_identity_invariants COSMOSIM_CI_BUILD_PARALLEL=2 \
  bash ./scripts/ci/run_preset_pipeline.sh cpu-only-debug build-cpu-debug test-cpu-debug \
  unit_gas_cell_identity_invariants /tmp/cosmosim_p9_pipeline 0
```

Observed results:

- CPU-only configure completed successfully.
- `unit_gas_cell_identity_invariants` built and passed.
- `integration_runtime_truth_ctest_labels` passed.
- The targeted `run_preset_pipeline.sh` invocation produced a passing JSON report and confirmed `COSMOSIM_CI_BUILD_TARGETS` is recorded in the artifact metadata.

Limitations:

- A full `stage1_runtime_truth_targets` build was started but exceeded the interactive container budget while compiling heavier gravity/PM-linked targets. The scoped CI mechanism is now present specifically to keep the Stage 1 gate exact and controllable, but full CI should still be run on a normal runner before release closure.
