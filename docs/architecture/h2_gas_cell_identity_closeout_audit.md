# H2 gas-cell identity closeout audit

Date: 2026-06-14
Prompt ID: H2.7
Scope: audit/documentation only; no solver numerics, AMR hydro implementation, config keys, output names, or restart datasets are changed by this closeout.

## Verdict

The H2 gas-cell identity design criterion is satisfied by code inspection: production hydro/runtime no longer requires local gas-cell count to equal local gas-particle count.

The remaining `requireParticleBoundGasCellContract(...)` executable callers are legacy/import compatibility and tests only. The production runtime now validates dense `SimulationState::gas_cell_identity` coverage, gathers/scatters hydro state by stable `gas_cell_id`, and treats `parent_particle_id` as optional lineage/mirror metadata rather than identity authority.

Command-backed H2 closure is blocked in this worktree by a CPU debug build failure in `src/core/simulation_state_species.cpp`: four calls to `requireGasCellIdentityMapCoversDenseRows(*this, "...")` resolve to the one-argument `SimulationState` member on MSVC instead of the intended free helper. H3 should not begin from this worktree until that compile blocker is fixed and the targeted H2 validation suite passes.

After that compile blocker is fixed, H3 may begin only on the APIs and blockers listed below. H3 must not claim mature AMR gas ownership until scheduler cell authority, patch exchange, ghost synchronization, and parentless/multi-cell production paths are promoted beyond the current compatibility bridges.

## Evidence basis

Inspected files required by H2.7:

- `src/workflows/reference_workflow.cpp`
- `src/core/simulation_state_active_views.cpp`
- `include/cosmosim/core/simulation_state.hpp`
- `src/core/simulation_state.cpp`
- `src/io/restart_checkpoint.cpp`
- `tests/unit/test_gas_cell_identity_invariants.cpp`
- `tests/integration/test_hydro_decoupled_gas_cells.cpp`

Additional H2 surfaces inspected:

- `docs/architecture/h2_gas_cell_identity_migration_plan.md`
- `src/core/simulation_state_ownership.cpp`
- `src/core/simulation_state_species.cpp`
- `src/core/time_integration.cpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `tests/integration/test_gas_cell_identity_migration.cpp`
- `tests/integration/test_gas_cell_split_merge_remap.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `CMakeLists.txt`
- `docs/output_schema.md`
- `docs/restart_checkpointing.md`

Discovery commands:

```bash
rg -n "requireParticleBoundGasCellContract" src include tests docs
rg -n "parentParticleRowForGasCellRow|parentParticleIdForGasCellRow|gasParticleIndexForCellRow|gasCellRowForParticleId" src include tests
rg -n "count\\(ParticleSpecies::kGas\\)|gas_particle_count|globalIndices\\(ParticleSpecies::kGas\\)|cells\\.size\\(\\)\\s*[!=]=|gas_cells\\.size\\(\\)\\s*[!=]=" src include tests
rg -n -C 3 "test_unit_gas_cell_identity_invariants|test_integration_gas_cell_identity_migration|test_integration_gas_cell_split_merge_remap|test_integration_hydro_decoupled_gas_cells|test_integration_restart_equivalence_hydro_toy" CMakeLists.txt
```

## `requireParticleBoundGasCellContract(...)` classification

| Call site | Classification | Closeout status |
|---|---|---|
| `src/core/simulation_state.cpp::legacyRequireParticleBoundGasCellContract` | Legacy/import compatibility implementation | Allowed. This is the named compatibility gate for old particle-bound layouts. |
| `src/core/simulation_state.cpp::requireParticleBoundGasCellContract` | Backward-compatible wrapper over the legacy gate | Allowed. Kept for tests/import callers and older API compatibility. |
| `include/cosmosim/core/simulation_state.hpp` declarations | Public compatibility API | Allowed. Production H2 paths should prefer identity-map coverage and ID lookups. |
| `tests/unit/test_gas_cell_identity_invariants.cpp` | Test-only, including negative legacy rejection | Allowed. Proves legacy particle-bound states still fail loudly when broken. |
| `tests/integration/test_species_migration_invariants.cpp` | Test-only legacy fixture coverage | Allowed. Exercises old particle-bound migration invariants. |
| `docs/**` mentions | Historical/audit/planning documentation | Allowed. Some older docs still describe the pre-H2 state. |

No production source caller outside the helper implementation remains in `src` after the H2.7 grep.

## Production hydro/runtime independence

Confirmed production hydro no longer depends on `cells.size() == local gas particle count`:

- `buildHydroCellKernelView(...)` in `src/core/simulation_state_active_views.cpp` requires dense `GasCellIdentityMap` coverage, gathers `gas_cell_id`, and captures `source_gas_cell_identity_generation`.
- `scatterHydroCellKernelView(...)` rejects stale identity generations and resolves destinations with `rowForGasCellId(...)`, not with parent particle row or active cell row.
- `HydroStageCallback::onStage` in `src/workflows/reference_workflow.cpp` advances hydro cells through cell-local conserved/primitive state. Parent-particle lookup is optional and used only to refresh compatibility particle mirrors when a local parent exists.
- `test_integration_hydro_decoupled_gas_cells` constructs 3 gas cells with only 1 gas particle parent, includes one parentless cell and two cells sharing one parent, confirms the legacy contract rejects the state, and then runs hydro plus gas-row reorder using stable `gas_cell_id`.

The remaining runtime parent lookups are H3 surfaces, not H2 stop-condition blockers:

- adaptive timestep metadata can emit a sentinel when no parent particle exists;
- drift/gravity/hydro callbacks update parent particle mirrors only when a parent is present;
- flux-correction and ghost paths still need H3 promotion to `gas_cell_id`-native exchange before AMR production maturity;
- scheduler cell mirror validation still has parent-backed compatibility behavior where a local parent exists.

## Parentless and multi-cell coverage

Confirmed tests exist and are registered:

- `tests/unit/test_gas_cell_identity_invariants.cpp`
  - parentless active hydro scatter by `gas_cell_id`;
  - duplicate `gas_cell_id` rejection;
  - duplicate local-row rejection;
  - zero `gas_cell_id` rejection;
  - multi-cell parent lookup through `rowsForParentParticleId(...)`;
  - `buildGasCellNewToOldRowMap(...)` row-remap by stable ID.
- `tests/integration/test_hydro_decoupled_gas_cells.cpp`
  - production-style hydro over `cells.size() != gas particle count`;
  - one parentless cell and two cells sharing one parent;
  - gas-row reorder while hydro state follows stable `gas_cell_id`.
- `tests/integration/test_gas_cell_identity_migration.cpp`
  - parentless gas-cell migration without a particle;
  - multi-cell same-parent migration keyed by `gas_cell_id`, not `parent_particle_id`;
  - stale gas ghost rejection by identity generation and hydro epoch.
- `tests/integration/test_gas_cell_split_merge_remap.cpp`
  - split case with two cells sharing one parent and conserved totals checked;
  - merge result with parentless identity;
  - row reorder scatter by `gas_cell_id`.

CTest registration was found in `CMakeLists.txt` for `unit_gas_cell_identity_invariants`, `integration_gas_cell_identity_migration`, `integration_gas_cell_split_merge_remap`, and `integration_hydro_decoupled_gas_cells`.

## Restart schema and migration coverage

Restart coverage is sufficient for H2 closeout:

- `RestartSchema` is v15 in `include/cosmosim/io/restart_checkpoint.hpp`.
- `src/io/restart_checkpoint.cpp` writes `/state/gas_cell_identity` with `gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`, `local_cell_row`, `local_row_reconstruction_policy=explicit_dense_local_cell_row`, and `identity_generation_at_write`.
- The restart reader accepts the documented v14 particle-bound import path only by materializing identity records from legacy mirrors.
- New v15 reads reject malformed identity records through `GasCellIdentityMap::assign(...)`, explicit parent-mask validation, dense-row coverage validation, sidecar mirror validation, and patch ownership validation.
- Restart integrity hashing includes sorted identity records for current schema payloads.
- `tests/unit/test_restart_checkpoint_schema.cpp` covers parentless identity hashing and rejects an invalid parent flag shape.
- `tests/integration/test_restart_checkpoint_roundtrip.cpp` checks identity record roundtrip and sidecar-lane agreement.

## H3 API map

H3 should use these existing APIs and data lanes rather than introducing a second identity system:

- `SimulationState::gas_cell_identity`: authoritative in-memory gas-cell identity map.
- `GasCellIdentityRecord::{gas_cell_id,parent_particle_id,owning_patch_id,local_cell_row}`: stable identity, optional lineage, stable patch owner, transient dense local row.
- `GasCellIdentityMap::rowForGasCellId(...)`: stable ID to current local row.
- `GasCellIdentityMap::gasCellIdForLocalRow(...)`: local row to stable ID.
- `GasCellIdentityMap::rowsForParentParticleId(...)`: optional lineage query; never uniqueness authority.
- `GasCellIdentityMap::rowsForPatch(...)`: patch-local row discovery for AMR ownership.
- `SimulationState::requireGasCellIdentityMapCoversDenseRows(...)`: production dense-row precondition.
- `SimulationState::requireGasCellIdentityMapFresh(...)`: stale-generation guard for cached views/remap plans.
- `SimulationState::gasCellIdentityMapMatchesSidecarLanes(...)`: mirror-lane validation before persistence or structural transforms.
- `buildGasCellNewToOldRowMap(...)`: row remap by stable `gas_cell_id`; absent old rows return `kInvalidGasCellRow` and require explicit initialization/conservation handling.
- `HydroCellKernelView::{gas_cell_id,local_cell_row,source_gas_cell_identity_generation}` plus `buildHydroCellKernelView(...)` / `scatterHydroCellKernelView(...)`: hydro hot-view identity/scatter contract.
- `SimulationState::packGasCellMigrationRecords(...)` and `SimulationState::commitGasCellMigration(...)`: gas-cell migration payload/commit surface keyed by `gas_cell_id`.
- Restart `/state/gas_cell_identity/*`: persistence contract for exact continuation.

## H3 blockers and guardrails

These are not H2 closeout blockers, but they are H3 prerequisites before AMR hydro integration can claim production readiness:

1. Scheduler cell authority must be keyed by `gas_cell_id` for parentless and split cells. Current parent-backed cell `time_bin` mirror validation is compatibility behavior when a local parent exists.
2. Runtime ghost and conservative flux exchange must address cells by `gas_cell_id` end to end. Parent IDs may remain diagnostics only.
3. Patch migration/rebalance must move gas-cell ownership with `GasCellMigrationRecord`/`GasCellMigrationCommit`, not by inferring patch cells from parent particle rows.
4. Parent-particle mirror writes in drift/gravity/hydro callbacks must remain optional and must not become conserved gas-cell truth.
5. Split/merge implementation must define conservation rules for `kInvalidGasCellRow` new cells and consumed old cells before mutating hydro state.
6. Multi-rank H3 work must exchange full gas-cell identity, scheduler identity, patch ownership, ghost epoch, and hydro fields. `parent_particle_id` and `CellSoa::time_bin` are not sufficient continuation authority.

## Validation status

Commands attempted:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "gas_cell|hydro_decoupled|migration|restart_equivalence_hydro_toy"
bash ./scripts/ci/check_repo_hygiene.sh
```

Outcomes:

- `cmake --preset cpu-only-debug` passed when invoked through Visual Studio's bundled CMake because plain `cmake` is not on this PowerShell `PATH`.
- The first build attempt without `VsDevCmd.bat` failed before repo code because `cl.exe` could not find standard library headers such as `<algorithm>` and `<array>`.
- The build attempt inside `VsDevCmd.bat` reached repo code and failed before linking/tests:
  - `src/core/simulation_state_species.cpp:365`
  - `src/core/simulation_state_species.cpp:496`
  - `src/core/simulation_state_species.cpp:1066`
  - `src/core/simulation_state_species.cpp:1088`
  - MSVC error: `SimulationState::requireGasCellIdentityMapCoversDenseRows` does not take two arguments.
- The requested CTest regex was not run because `build-cpu-debug` did not complete.
- `bash ./scripts/ci/check_repo_hygiene.sh` did not run because this Windows environment routes `bash` to WSL and WSL reports no installed Linux distributions.

## Reproducibility impact

This closeout audit is documentation-only. It changes no runtime behavior, solver numerics, restart schema version, HDF5 dataset names, config normalization, provenance payloads, output naming, rank coordination, or scheduler semantics.
