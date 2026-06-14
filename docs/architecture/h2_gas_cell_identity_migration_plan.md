# H2 gas-cell identity migration plan

Date: 2026-06-13
Prompt ID: H2.0
Scope: production migration plan only; no runtime behavior, restart schema, hydro solver, AMR, or moving-mesh changes.

## Verdict

H2 may proceed. `docs/architecture/h1_hydro_closeout_audit.md` marks H1 green for the targeted hydro, CFL, validation, HDF5 restart, and HDF5 reference-workflow gates, with the explicit caveat that production hydro remains particle-bound. This plan treats that caveat as the H2 migration target.

Current production authority is still the temporary particle-bound gas-cell contract:

- `GasCellSidecar::{gas_cell_id,parent_particle_id}` persist identity lanes, but `gas_cell_id == parent_particle_id == gas particle_id` is required whenever local gas cells exist.
- `core::requireParticleBoundGasCellContract(...)` is the production gate.
- `GasCellIdentityMap` exists as an isolated RFC seam and is not stored in `SimulationState`, restart payloads, hydro active views, scheduler mirror refresh, or workflow migration authority.

## H2 authority decision

H2 should promote one core-owned `GasCellIdentityMap` to authoritative runtime state. The owner should be `core::SimulationState`, unless H2.1 finds a compelling reason to introduce a narrow `core` ownership component that is still stored by `SimulationState` and serialized through restart as part of state.

The promoted map owns:

- stable `gas_cell_id`, unique and nonzero;
- optional `parent_particle_id`, where absence is represented explicitly and zero is not overloaded as absent;
- stable `owning_patch_id`, naming the patch owner rather than a local patch row;
- transient `local_cell_row`, dense over local `CellSoa`/`GasCellSidecar` rows for a given generation;
- generation counter advanced on any row, ownership, split, merge, or patch remap.

Hydro, restart, scheduler mirror refresh, migration, and diagnostics must consume gas-cell rows through this one map. No hydro-local identity map is allowed.

## Current particle-bound call-site classification

| File / symbol | Current call | Classification | H2 disposition |
|---|---|---|---|
| `src/core/simulation_state.cpp::requireParticleBoundGasCellContract` | Defines temporary contract | Legacy import compatibility | Keep through H2.6 as a compatibility validator for materializing legacy particle-bound maps. Do not remove in H2. |
| `src/core/simulation_state.cpp::parentParticleIdForGasCellRow` | Requires particle-bound map | Legacy compatibility helper | Replace production callers with identity-map lookups; keep helper for legacy tests and fallback imports. |
| `src/core/simulation_state.cpp::gasCellRowForParticleId` | Linear parent-particle lookup | Legacy compatibility helper | Stop using for production hydro/flux correction. Add `rowForGasCellId` and optional `rowsForParentParticleId` call sites. |
| `src/core/simulation_state.cpp::gasParticleIndexForCellRow` | Row-to-gas-particle lookup | Legacy compatibility helper | Remove from production hydro/scheduler paths after gas velocity/mass storage is decoupled from particle rows. Keep only legacy fixtures/imports. |
| `src/core/simulation_state.cpp::reorderParticles` | Rejects relative gas-particle reorder | Migration guard | Replace with gas-cell row remap by stable `gas_cell_id`; keep rejection only while state is in legacy particle-bound mode. |
| `src/core/simulation_state_active_views.cpp::buildHydroCellKernelView` | Requires particle-bound state before active hydro view | Production hydro | Migrate to `state.gas_cell_identity.requireCoversDenseLocalRows(...)` plus generation capture. |
| `src/core/simulation_state_active_views.cpp::scatterHydroCellKernelView` | Requires particle-bound state before scatter | Production hydro | Migrate to identity-map generation validation plus dense-row validation. |
| `src/core/simulation_state_ownership.cpp::validateOwnershipInvariants` | Requires particle-bound when cells exist | Migration | Replace with identity-map consistency, dense row coverage, gas sidecar count, patch ownership, and legacy fallback validation. |
| `src/core/simulation_state_species.cpp::packParticleMigrationRecords` | Packs gas cell fields by parent particle ID | Migration | Split particle migration from gas-cell migration. Gas-cell migration records must be keyed by `gas_cell_id`, with optional parent metadata. |
| `src/core/simulation_state_species.cpp::commitParticleMigration` | Rebuilds gas cells from kept/inbound gas particles | Migration | Replace with gas-cell remap/commit keyed by `gas_cell_id`; particle species migration must not imply one gas cell per gas particle. |
| `src/core/time_integration.cpp::syncGasCellTimeBinMirrorsFromParticleScheduler` | Maps cell bins through parent gas particle | Scheduler mirror | H2 must define scheduler element identity for gas cells. Until then this stays legacy-only; production cell mirror refresh must use scheduler records keyed by `gas_cell_id`. |
| `src/core/time_integration.cpp::timeBinMirrorsMatchScheduler` | Validates cell mirrors through parent particle | Scheduler mirror | Migrate to gas-cell scheduler identity records keyed by `gas_cell_id`, with legacy fallback only for old payloads. |
| `src/core/time_integration.cpp::debugAssertTimeBinMirrorAuthorityInvariant` | Asserts cell mirrors through parent particle | Scheduler mirror | Same migration as mirror matching; no mirror may become authority. |
| `src/io/restart_checkpoint.cpp::validateRestartTimeBinMirrorsAgainstScheduler` | Requires particle-bound cells to validate cell mirrors | Restart | Keep for schema v14 legacy. New schema must validate serialized gas identity map and scheduler gas-cell identity records. |
| `src/io/restart_checkpoint.cpp::rebuildRestartTimeBinMirrorsFromScheduler` | Rebuilds cell mirrors through parent particle | Restart | H2 schema upgrade must rebuild cell mirrors by `gas_cell_id` scheduler authority. |
| `src/workflows/reference_workflow.cpp::buildAdaptiveTimeStepCriteriaView` | Stores gas particle index per cell | Production hydro / scheduler mirror | Replace with gas-cell IDs and gas-cell velocity/mass lanes or explicit legacy velocity source adapter. |
| `src/workflows/reference_workflow.cpp::collectLocalGasCellRecords` | Collects records by gas particle index/ID | Migration | Replace with records keyed by `gas_cell_id`; parent particle is optional metadata only. |
| `src/workflows/reference_workflow.cpp::gasParticleIdByOldCellIndex` | Builds old row parent-particle list | Migration | Replace with `gasCellIdByOldCellIndex` from identity map. |
| `src/workflows/reference_workflow.cpp::rebuildLocalGasStateFromParticleIds` | Rebuilds rows by parent particle ID | Migration | Replace with `rebuildLocalGasStateFromGasCellIds` using `buildGasCellNewToOldRowMap`. |
| `src/workflows/reference_workflow.cpp::estimateDecompositionItemMemoryBytesForDecomposition` and load-balance helpers | Use parent gas particle to infer patch cost | Migration / AMR readiness | Use gas-cell identity and patch ownership directly; particle entity cost should not imply cell ownership. |
| `src/workflows/reference_workflow.cpp::buildLocalAmrPatchCellPayloadRecords` | Derives gas identity from parent particle | Migration / AMR readiness | Emit stable `gas_cell_id`, optional parent, and patch ID from authoritative map. |
| `src/workflows/reference_workflow.cpp::runtime rebalance AMR patch loop` | Moves particles for patch cell ownership | Migration / AMR readiness | Replace with gas-cell ownership transfer keyed by `gas_cell_id`; parent particles migrate only when required by a documented policy. |
| `src/workflows/reference_workflow.cpp::snapshotHydroGhostConservedCells` | Ghost records use parent particle IDs | Production hydro / ghost exchange | Change flux correction and ghost identity to `gas_cell_id`; keep parent only as compatibility metadata. |
| `src/workflows/reference_workflow.cpp::initializeReferenceSchedulerBins` | Sets scheduler bins for cell row and parent particle row | Scheduler mirror | Introduce scheduler identity records for gas cells; no positional cell row equals scheduler row assumption. |
| `src/workflows/reference_workflow.cpp::DriftStageCallback::onStage` | Mirrors cell centers from parent particles | Production hydro | Remove as production rule after gas cells own centers. Legacy particle-bound adapter may still refresh generated toy states. |
| `src/workflows/reference_workflow.cpp::GravityStageCallback::onStage` | Mirrors gas-cell acceleration from parent active particle | Production hydro / gravity coupling | Define gravity source coupling by gas-cell center/mass or explicit parent adapter; do not require one active particle per cell. |
| `src/workflows/reference_workflow.cpp::HydroStageCallback::onStage` | Reads/writes gas velocity and mass through parent particles | Production hydro | Migrate hydro primitive/conserved load/store to gas-cell-owned velocity/momentum/mass contract or explicit adapter. |
| `src/workflows/reference_workflow.cpp::HydroStageCallback::flux correction` | Maps correction parent_particle_id to cell row | Production hydro / ghost exchange | Change correction records to `gas_cell_id`; parent ID remains optional diagnostic/legacy field. |
| `src/workflows/reference_workflow.cpp::hydroCflDiagnosticsForCell` | Gets velocity through parent particle | Production hydro | Use gas-cell velocity/momentum contract after H2.3. |
| `tests/unit/test_gas_cell_identity_invariants.cpp` | Tests both temporary contract and isolated map | Test-only | Keep legacy tests; add production map tests as H2 stages wire authority. |
| `tests/integration/test_species_migration_invariants.cpp` | Tests gas rebuild by particle ID | Test-only | Add gas-cell-ID split/merge/remap migration tests before replacing behavior. |
| `tests/integration/test_reorder_compaction_sidecars.cpp` | Tests stale hydro view after gas identity rebuild | Test-only | Extend to identity-map generation invalidation and gas-cell row remap. |
| `tests/integration/test_restart_checkpoint_roundtrip.cpp` | Tests restart gas density by particle ID and particle-bound mirror validation | Test-only / restart | Add schema-versioned gas identity map roundtrip and legacy materialization tests. |

## Restart and schema plan

Do not change restart schema in H2.0. H2.3 or H2.4 must introduce a versioned restart schema bump after the in-memory authority is green.

Required new restart state:

- `/state/gas_cell_identity/gas_cell_id`: nonzero stable IDs, unique.
- `/state/gas_cell_identity/has_parent_particle_id`: explicit optional mask.
- `/state/gas_cell_identity/parent_particle_id`: value lane interpreted only where the mask is true.
- `/state/gas_cell_identity/owning_patch_id`: stable owner ID, not local row.
- `/state/gas_cell_identity/local_cell_row` or a documented deterministic row reconstruction rule.
- `/state/gas_cell_identity/generation`: optional diagnostic mirror only; reader may assign a fresh generation after import.

Compatibility behavior:

- Legacy schema v14 materializes `GasCellIdentityMap` with `gas_cell_id == parent_particle_id == particle_sidecar.particle_id[gas species row]`, `owning_patch_id` from `patches.patch_id[cells.patch_index[row]]` when patches exist, and dense `local_cell_row`.
- Legacy materialization must reject `parent_particle_id == 0` as absent. Absence exists only through `has_parent_particle_id == 0` in the new schema.
- New schema must not require unique `parent_particle_id`; split cells may share a parent and parentless cells must roundtrip.
- Payload integrity hashing must include all identity-map lanes after the schema bump.

## Hot-path lookup strategy

H2 production code must avoid repeated linear scans:

- Keep `GasCellIdentityMap::rowForGasCellId(gas_cell_id)` and `gasCellIdForLocalRow(local_cell_row)` as O(1) lookup APIs backed by the existing unordered maps.
- Add dense arrays beside or inside the map when profiling shows unordered lookups in hot loops:
  - `gas_cell_id_by_local_row[cell_row]`;
  - `local_row_by_dense_identity_index` if a sorted/packed gas-cell iteration order is needed;
  - `scheduler_element_by_gas_cell_row` or `scheduler_identity_record_by_gas_cell_id` for cell timestep mirror refresh.
- Hydro loops must gather once per stage into compact arrays and then index dense local rows, not call map lookups per face.
- Split/merge/remap operations must build one `new_row -> old_row` vector with `buildGasCellNewToOldRowMap(...)` and then move all hydro/cell/patch lanes in one commit boundary.

## H2 staged patch plan

### H2.1: Promote identity map into core state without changing production semantics

Files/symbols:

- `include/cosmosim/core/simulation_state.hpp`: add `GasCellIdentityMap gas_cell_identity` to `SimulationState`; add `refreshGasCellIdentityMapFromParticleBoundState`, `gasCellIdentityMapMatchesParticleBoundState`, and dense-row invariant helpers.
- `src/core/simulation_state.cpp`: implement legacy materialization from current particle-bound lanes.
- `src/core/simulation_state_ownership.cpp`: validate both old lanes and the materialized map while in legacy mode.
- `tests/unit/test_gas_cell_identity_invariants.cpp`: prove map materialization from legacy state, dense coverage, generation update, parentless rejection in legacy mode, and no duplicate `gas_cell_id`.
- `docs/state_model_memory_layout.md` and `docs/architecture/adr_runtime_truth_ownership.md`: record `SimulationState::gas_cell_identity` as the one future authority while particle-bound mode remains enabled.

Stop conditions:

- H1 audit no longer green.
- Materialized map can drift from `gas_cells.{gas_cell_id,parent_particle_id}` without an invariant failure.
- Any production caller starts accepting parentless/multi-parent cells before restart and hydro storage contracts exist.

Validation floor:

- `cmake --preset cpu-only-debug`
- `cmake --build --preset build-cpu-debug --target test_unit_gas_cell_identity_invariants`
- `ctest --preset test-cpu-debug --output-on-failure -R "gas_cell_identity"`

### H2.2: Replace migration and reorder row remaps with gas_cell_id keys

Files/symbols:

- `src/core/simulation_state.cpp::reorderParticles`: allow gas-particle relative reorder only when gas-cell rows are remapped by `gas_cell_id` in the same commit.
- `src/core/simulation_state_species.cpp::{packParticleMigrationRecords,commitParticleMigration}`: stop deriving cell continuity from parent particle uniqueness; introduce gas-cell migration payloads keyed by `gas_cell_id`.
- `src/workflows/reference_workflow.cpp::{collectLocalGasCellRecords,gasParticleIdByOldCellIndex,rebuildLocalGasStateFromParticleIds}`: replace with gas-cell-ID keyed helpers.
- `tests/integration/test_species_migration_invariants.cpp`: add split and parentless gas-cell cases; assert conserved hydro fields follow `gas_cell_id`, not particle row or parent ID.
- `tests/integration/test_reorder_compaction_sidecars.cpp`: add gas-cell row reorder with stable IDs and stale active-view generation rejection.

Stop conditions:

- A kept or migrated gas cell can be lost because its parent particle is absent or non-unique.
- Any path uses `parent_particle_id` as a map key for conserved hydro state.
- Reorder/remap mutates `cells` and `gas_cells` in separate commits.

Validation floor:

- `cmake --build --preset build-cpu-debug --target test_integration_species_migration_invariants test_integration_reorder_compaction_sidecars`
- `ctest --preset test-cpu-debug --output-on-failure -R "migration|reorder|gas|sidecar"`

### H2.3: Define gas-cell-owned hydro runtime lanes and active-view generation

Files/symbols:

- `include/cosmosim/core/simulation_state.hpp`: add or formalize gas-cell velocity/momentum ownership needed to stop reading velocities from parent particles. If adding lanes, use explicit `_peculiar` / `_comov` / `_code` suffixes and document restart impact.
- `src/core/simulation_state_active_views.cpp::{buildHydroCellKernelView,scatterHydroCellKernelView}`: validate dense identity-map coverage and capture identity-map generation.
- `src/workflows/reference_workflow.cpp::{buildAdaptiveTimeStepCriteriaView,DriftStageCallback::onStage,GravityStageCallback::onStage,HydroStageCallback::onStage,hydroCflDiagnosticsForCell}`: route hydro load/store and CFL through gas-cell-owned lanes or a clearly named legacy adapter.
- `tests/unit/test_gas_cell_identity_invariants.cpp`: active hydro view rejects stale identity-map generations.
- Hydro integration tests: add a parentless cell and a two-cells-one-parent fixture that exercises active hydro without particle-row velocity authority.

Stop conditions:

- Hydro writes particle mass/velocity as the only owner of gas-cell conserved state.
- Active hydro view can scatter after identity-map generation changes.
- CFL diagnostics still require `gasParticleIndexForCellRow`.

Validation floor:

- `cmake --build --preset build-cpu-debug --target test_unit_gas_cell_identity_invariants test_integration_hydro_sod_like test_integration_hydro_conservation_periodic`
- `ctest --preset test-cpu-debug --output-on-failure -R "gas_cell_identity|integration_hydro_|hydro_conservation"`

### H2.4: Version restart schema and legacy materialization

Files/symbols:

- `include/cosmosim/io/restart_checkpoint.hpp`: bump `RestartSchema` only after H2.1-H2.3 in-memory tests are green.
- `src/io/restart_checkpoint.cpp`: write/read `/state/gas_cell_identity/*`; include lanes in integrity hash; keep schema v14 legacy materialization path if compatibility is required.
- `docs/output_schema.md` and `docs/restart_checkpointing.md`: document new fields, optional-parent mask, legacy fallback, and deterministic row reconstruction.
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`: roundtrip split, parentless, and patch-owned cells; verify legacy particle-bound files materialize an equivalent map.
- `tests/unit/test_restart_checkpoint_schema.cpp`: add missing/duplicate/zero gas ID, absent-parent mask, and stale local-row negative cases.

Stop conditions:

- New schema uses `parent_particle_id == 0` to mean absent.
- Payload hash omits any identity lane.
- Reader silently accepts duplicate `gas_cell_id` or non-dense local rows.

Validation floor:

- `cmake --preset hdf5-debug`
- `cmake --build --preset build-hdf5-debug --target test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip`
- `ctest --preset test-hdf5-debug --output-on-failure -R "restart_checkpoint|restart.*gas|gas_cell_identity"`

### H2.5: Scheduler mirror migration to gas-cell identity

Files/symbols:

- `include/cosmosim/core/time_integration.hpp` and `src/core/time_integration.cpp`: add scheduler identity records for gas cells keyed by `gas_cell_id`; update cell mirror refresh and mirror validation to stop routing through parent particle rows.
- `src/workflows/reference_workflow.cpp::initializeReferenceSchedulerBins`: register gas-cell scheduler elements by stable identity rather than setting both `cell_index` and parent particle index opportunistically.
- `tests/unit/test_time_integration.cpp`: add gas-cell scheduler mirror tests for split cells and parentless cells.
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`: verify exact scheduler continuation by `gas_cell_id` across restart.

Stop conditions:

- `CellSoa::time_bin` is treated as scheduler authority.
- Scheduler identity exchange still requires one gas particle per gas cell.
- New cell scheduler records cannot be rebuilt after split/merge.

Validation floor:

- `cmake --build --preset build-cpu-debug --target test_unit_time_integration test_integration_restart_checkpoint_roundtrip`
- `ctest --preset test-cpu-debug --output-on-failure -R "time_integration|scheduler|restart.*gas"`

### H2.6: Retire production particle-bound call sites

Files/symbols:

- `src/workflows/reference_workflow.cpp`: remove production uses of `requireParticleBoundGasCellContract`, `gasParticleIndexForCellRow`, `gasCellRowForParticleId`, and `parentParticleIdForGasCellRow` except in explicit legacy import/materialization code.
- `src/core/simulation_state_active_views.cpp`: require identity-map coverage/generation only.
- `src/io/restart_checkpoint.cpp`: keep old contract only in legacy schema reader.
- `docs/architecture/gas_cell_identity_map_rfc.md`: mark RFC as promoted and point to production ownership docs.
- `docs/repair_open_issues.md`: update only if a real blocker remains.

Stop conditions:

- `rg -n "requireParticleBoundGasCellContract|gasParticleIndexForCellRow|gasCellRowForParticleId|parentParticleIdForGasCellRow" src include` shows production callers outside legacy compatibility, tests, or deprecation wrappers.
- H2 tests pass only because parent-particle uniqueness is restored in fixtures.
- AMR patch payloads still derive cell ownership from particle rows.

Validation floor:

- `cmake --preset cpu-only-debug`
- `cmake --build --preset build-cpu-debug`
- `ctest --preset test-cpu-debug --output-on-failure`
- HDF5 restart gate from H2.4.

## Test plan before implementation

Required tests before removing the particle-bound production gate:

- Split: one parent particle produces two gas cells with distinct `gas_cell_id`; both preserve density, pressure, mass, velocity/momentum, patch ownership, and scheduler mirror state by `gas_cell_id`.
- Merge: two gas cells merge into one new `gas_cell_id`; old rows are consumed through an explicit conservation rule rather than accidental row overwrite.
- Parentless: a gas cell with no parent particle roundtrips through runtime ownership, hydro active view, restart, and scheduler mirror refresh.
- Remap: row order changes while `gas_cell_id` stays stable; all hydro fields and host-cell sidecars follow the stable ID.
- Restart legacy fallback: schema v14 particle-bound state materializes a map with `gas_cell_id == parent_particle_id`, while new schema preserves absent parents through a mask.
- Scheduler mirror: cell `time_bin` mirrors refresh from scheduler gas-cell records; stale mirrors are rejected without consulting parent particle rows.
- Ghost/flux correction: conservative correction records address cells by `gas_cell_id`, not parent particle ID.

## Reproducibility impact

H2.0 is documentation-only. It changes no solver numerics, config keys, restart schema, HDF5 dataset names, output naming, scheduler behavior, or rank coordination. Future H2.4 schema work must explicitly update restart/output docs and validation evidence before any compatibility claim.
