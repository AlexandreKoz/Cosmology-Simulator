# Stage 1 runtime transform audit: mutation-map proof

_Date captured: 2026-05-06 (UTC)_

Scope: baseline audit freeze for particle, gas-cell, species-sidecar, restart, snapshot, migration, and reorder mutation paths. This note is documentation-only evidence for later repair prompts; it does not authorize physics, solver, data-layout, or schema changes.

## Audit commands and files inspected

- Baseline archive did not include `.git/` metadata in the reviewed zip; mutation-map proof therefore used source/diff inspection rather than repository status claims.
- Source/doc sweep used `rg` for `SimulationState`, sidecar, ownership, runtime-truth, reorder, migration, restart, snapshot, gas, and softening terms across `tests`, `include`, `src`, and architecture/repair docs.
- Targeted source reads covered:
  - `include/cosmosim/core/simulation_state.hpp`
  - `src/core/simulation_state.cpp`
  - `src/core/simulation_state_active_views.cpp`
  - `src/core/simulation_state_ownership.cpp`
  - `src/core/simulation_state_species.cpp`
  - `src/core/simulation_state_structures.cpp`
  - `src/workflows/reference_workflow.cpp`
  - `src/io/restart_checkpoint.cpp`
  - `src/io/snapshot_hdf5.cpp`
  - `docs/architecture/adr_runtime_truth_ownership.md`
  - `docs/architecture/runtime_truth_map.md`
- Sidecar/ownership tests inspected by name/discovery: `test_unit_simulation_state`, `test_unit_hot_cold_sidecar_layout`, `test_unit_gas_cell_identity_invariants`, `test_integration_reorder_compaction_sidecars`, `test_integration_species_migration_invariants`, `test_integration_softening_ownership_invariants`, and `test_integration_restart_checkpoint_roundtrip`.

## Authoritative owners and derived mirrors

| Domain | Authoritative owner/storage | Derived mirror/cache/view | Authorized mutators observed | Current audit classification |
|---|---|---|---|---|
| Particle hot SoA rows | `core::SimulationState::particles` row order, paired with `particle_sidecar` IDs/species/rank metadata | Active views, gravity kernel compact views, transfer/migration records, snapshot/restart payloads | `SimulationState::resizeParticles`, `reorderParticles`, `SimulationState::commitParticleMigration`, restart/snapshot readers for new state construction, solver scatter through `scatterGravityParticleKernelView` for hot physics lanes only | Implemented, with tests for consistency, active view derivation, reorder, and migration invariants |
| Particle sidecar metadata | `core::SimulationState::particle_sidecar` (`particle_id`, `sfc_key`, `species_tag`, flags, owning rank, softening lanes) | `ParticleSpeciesIndex`, `SpeciesContainer`, transfer packets, snapshot groups, restart groups | `resizeParticles`, `reorderParticles`, `commitParticleMigration`, `rebuildSpeciesIndex`; workflow decomposition may set `owning_rank` before compaction | Implemented, but workflow direct `owning_rank` assignment is an acknowledged ad-hoc path requiring owner-aware framing in later repair |
| Species ledger/index | `particle_sidecar.species_tag` as row truth; `species.count_by_species` and `particle_species_index` as maintained mirrors | Species-local index spans and I/O type grouping | `ParticleSpeciesIndex::rebuild`, `SimulationState::rebuildSpeciesIndex`, `commitParticleMigration`, snapshot/restart read reconstruction | Implemented; direct test fixture writes exist but production mutation should continue to go through rebuild/commit APIs |
| Species sidecars | `star_particles`, `black_holes`, `tracers` sidecar rows keyed by `particle_index` | Migration records and snapshot tracer datasets | `resize` for sidecar containers, `reorderParticles` via remapped parent indirection or physical row move, `commitParticleMigration`; gas rebuild remaps black-hole/tracer `host_cell_index` after cell reconstruction | Implemented for star/BH/tracer migration and reorder; snapshot support is partial/observer-like because only tracer sidecars round-trip through snapshot |
| Gas cell row identity | Stable gas-particle IDs (`particle_sidecar.particle_id` for gas species) reflected into `gas_cells.gas_cell_id` and `parent_particle_id`; cell rows are local/transient | `cells`/`gas_cells` row order keyed to canonical gas species order; active hydro views | `refreshGasCellIdentityFromParticleOrder`, `rebuildLocalGasStateFromParticleIds`, guarded `reorderParticles`, restart read/write | Implemented for temporary 1:1 gas-particle/gas-cell contract; still unsafe where tests or initial-condition builders directly write gas arrays before identity refresh |
| Restart continuation truth | HDF5 restart schema payload plus typed `RestartReadResult` | Integrity hash, normalized config/provenance datasets, reconstructed scheduler | `writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`, private `writeStateGroup`/`readStateGroup` | Implemented for full continuation; schema changes remain forbidden without version/docs/tests |
| Snapshot interchange truth | Gadget/AREPO-style HDF5 snapshot payload; by-type particle groups | Snapshot read state re-derived from groups and config defaults | `writeGadgetArepoSnapshotHdf5`, `readGadgetArepoSnapshotHdf5` | Partially implemented by design: snapshots are not restart continuation truth and currently reconstruct particles/species/tracers/softening but not full gas-cell sidecar or scheduler state |
| Reorder/compaction/migration transforms | `SimulationState` transform APIs and workflow gas rebuild helpers | Migration packets, gas-cell migration records, sidecar row-order maps | `buildParticleReorderMap`, `reorderParticles`, `packParticleMigrationRecords`, `commitParticleMigration`, `compactStateToCurrentOwner`, `rebuildLocalGasStateFromParticleIds` | Implemented with targeted tests; gas ownership compaction is split between core particle commit and workflow gas rebuild |

## Authoritative mutators by subsystem

### Core `SimulationState` structure and view mutators

- `SimulationState::resizeParticles` and `SimulationState::resizeCells` resize the paired SoA/sidecar domains and bump the corresponding generation counters.
- `ParticleSidecar::{setGravitySofteningOverride,clearGravitySofteningOverride}` are the only explicit per-particle override-mask mutators; the value vector alone is a materialized/default mirror unless the mask marks an override.
- `SimulationState::rebuildSpeciesIndex` rebuilds the derived species index and opportunistically initializes gas identity only when gas-cell rows are 1:1 and identity lanes are uninitialized.
- `buildParticleActiveView` and `buildCellActiveView` create read-only compact views; they are derived views only.
- `buildGravityParticleKernelView` and `buildHydroCellKernelView` create mutable compact kernel workspaces; only their scatter functions are authorized to write back hot lanes, and generation counters reject stale scatter after resize/reorder.
- `validateOwnershipInvariants`, `validateUniqueParticleIds`, `debugAssertNoStaleParticleIndices`, `gasCellIdentityMatchesParticleOrder`, and `debugAssertGasCellIdentityContract` are guards/observers, not mutators.

### Particle reorder paths

- `buildParticleReorderMap` builds explicit old-to-new/new-to-old maps for `kByTimeBin`, `kBySfcKey`, and `kBySpecies`.
- `reorderParticles` applies the permutation to particle SoA lanes and particle sidecar lanes, optionally moves species sidecar rows with parents or remaps sidecar `particle_index` indirection, rebuilds species index, and bumps the particle generation.
- `reorderParticles` rejects gas-cell reorders that would violate the current 1:1 gas identity contract unless the reorder leaves gas order compatible with gas-cell identity.

### Particle migration paths

- `packSpeciesTransferPacket` is a read-only packer for species-local transfer data.
- `packParticleMigrationRecords` serializes all particle hot lanes, common sidecar lanes, softening value/mask state, and required star/BH/tracer sidecar fields for selected local indices.
- `commitParticleMigration` is the core authority for removing outbound/stale particles, appending inbound particles, rebuilding star/BH/tracer sidecars, preserving softening override semantics, recomputing species counts/index, and bumping the particle generation.
- `ParticleMigrationCommit` explicitly separates `outbound_local_indices`, `inbound_records`, and `stale_local_ghost_indices`; duplicate IDs and mismatched inbound sidecar contracts are rejected.

### Workflow gas-cell and rank-ownership paths

- `collectLocalGasCellRecords` captures gas-cell/cell thermodynamic and geometry lanes keyed by gas particle ID after enforcing the gas identity contract.
- `gasParticleIdByOldCellIndex` captures old cell index to gas-particle ID mapping for host-cell remaps.
- `rebuildLocalGasStateFromParticleIds` rebuilds `cells` and `gas_cells` from gas-particle IDs after particle compaction/migration, writes gas identity lanes from particle IDs, bumps the cell generation, and remaps black-hole/tracer `host_cell_index` values.
- `compactStateToCurrentOwner` collects outbound particle indices from `particle_sidecar.owning_rank`, invokes `commitParticleMigration`, then rebuilds gas state from captured gas records.
- `applyInitialGravityAwareDecomposition` directly writes `particle_sidecar.owning_rank` from a decomposition plan before compaction. This is currently an ad-hoc but known rank-ownership mutation path.
- `initializeSchedulerBins`/`syncTimeBinsFromScheduler` own scheduler bin setup/sync; `particles.time_bin` and `cells.time_bin` remain scheduler mirrors during stepping and must not become scheduling authority.

### Restart mutation paths

- `writeRestartCheckpointHdf5` validates continuation metadata and `state.validateOwnershipInvariants()` before writing.
- `writeStateGroup` serializes particle SoA, particle sidecars, cell/gas sidecars including gas identity lanes, patches, species counts, star/BH/tracer sidecars, metadata, and module sidecars.
- `readStateGroup` reconstructs those same state lanes, requires gas identity datasets, reads optional legacy softening lanes with documented behavior, rebuilds species index, and validates the reconstructed state.
- `readRestartCheckpointHdf5` reads schema/integrity/provenance/integrator/scheduler/distributed-gravity payloads and recomputes the restart payload hash before returning typed state.

### Snapshot mutation paths

- `writeGadgetArepoSnapshotHdf5` observes `SimulationState`, groups particles by mapped species/PartType, writes coordinates/velocities/masses/IDs, optional softening values/masks, tracer sidecar datasets, config, and provenance.
- `readGadgetArepoSnapshotHdf5` constructs a new `SimulationState` from snapshot PartType groups, sets particle hot lanes, IDs, species tags, owner rank default, time-bin default, optional softening lanes/masks, species counts, and tracer sidecar lanes; it calls `rebuildSpeciesIndex` afterward.
- Snapshot read/write is not a restart mutator: scheduler state, integrator continuation truth, full gas-cell thermodynamics/identity, star sidecars, and black-hole sidecars are not complete continuation payloads in the snapshot path.

## Forbidden mutators and anti-drift rules for later prompts

- Do not add direct writes to `particle_sidecar.species_tag` in production paths without same-call species ledger/index rebuild and tests.
- Do not add direct writes to `particles.time_bin` or `cells.time_bin` as scheduling authority; these are mirrors refreshed from scheduler state.
- Do not mutate gas-cell identity lanes except through `refreshGasCellIdentityFromParticleOrder`, restart read reconstruction, or `rebuildLocalGasStateFromParticleIds`.
- Do not make snapshot read/write a silent restart schema surrogate; any continuation-payload expansion belongs in restart schema/version/docs/tests.
- Do not add module-local active-set builders, sidecar caches, or rank ownership stores that duplicate scheduler/state authority.
- Do not silently infer per-particle softening override authority from a populated value vector without `has_gravity_softening_override`.
- Do not change solver numerics, SoA layout, HDF5 schema names, or restart compatibility behavior in response to this audit.

## Still-unsafe or ad-hoc paths found by `rg`

- Test fixtures directly fill state lanes (`particle_sidecar`, `species.count_by_species`, gas/cell arrays, sidecars) before invoking invariant/rebuild APIs. This is acceptable test setup but should not be copied into runtime code.
- `applyInitialGravityAwareDecomposition` directly assigns `particle_sidecar.owning_rank` before `compactStateToCurrentOwner`; later repair may want an explicit rank-ownership mutator or naming contract, but this prompt only records the path.
- Snapshot read constructs new state from PartType groups and defaults missing continuation lanes (`time_bin = 0`, owner rank `0`, missing tracer optional fields to zero). This is intentional interchange behavior, not restart continuation.
- Gas-cell state mutation is split between core particle migration and workflow gas rebuild; later repairs must cite both authorities rather than treating either one as complete by itself.
- `SimulationState::rebuildSpeciesIndex` may initialize gas identity if the state is 1:1 and identity lanes are empty/zero. This is a controlled compatibility convenience but should not become a general gas-cell mutation surface.

## P0 target classification for this Stage 1 audit gate

Because no separate in-repository Stage 1 target list was found by `rg "1\\.1|1\\.2|1\\.3|Stage 1|stage1|runtime transform"`, the audit classifies the three requested P0 targets by their prompt scope:

| P0 target | Scope used for this audit | Classification | Evidence basis |
|---|---|---|---|
| 1.1 | Particle/canonical-order mutation authority and sidecar synchronization map | Implemented | Core resize/reorder/migration APIs and ownership tests are present; direct test setup writes are not runtime mutators. |
| 1.2 | Gas-cell identity/rebuild and restart/snapshot transform boundaries | Partially implemented | Gas/restart authority is implemented; snapshot remains a partial interchange path by design and not a full gas/restart continuation transform. |
| 1.3 | Species-sidecar, migration, and reorder transform invariants | Implemented | Star/BH/tracer sidecar migration/reorder paths and targeted tests are present; gas rebuild remains workflow-owned and separately recorded under 1.2. |

## Validation report for this prompt

Required command-backed validation for this documentation-only prompt:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure -R "simulation_state|ownership|runtime|truth"
```


Observed result:

- `cmake --preset cpu-only-debug` passed.
- `cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state` passed.
- The first `ctest --preset test-cpu-debug --output-on-failure -R "simulation_state|ownership|runtime|truth"` selected two integration tests whose executables had not been built by the narrower unit target and failed as `Not Run`; this was a validation-order issue, not a source/test-registration change.
- After `cmake --build --preset build-cpu-debug -j4`, the same `ctest --preset test-cpu-debug --output-on-failure -R "simulation_state|ownership|runtime|truth"` passed 3/3.

Reproducibility impact: documentation-only audit. No normalized config dumps, provenance payloads, restart/snapshot schemas, output naming, deterministic scheduling behavior, solver numerics, or SoA/hot-cold layouts are changed.

## CI gate consolidation update: runtime-truth repair suite

The Stage 1 P0 runtime-truth tests are now registered with explicit `runtime_truth`/`p0` CTest labels and subsystem labels for softening, sidecar, gas, migration, restart, provenance, scheduler, and active-view failures. The dependency-free local/CI closeout command is:

```bash
./scripts/ci/run_stage1_runtime_truth_gate.sh ci_artifacts/stage1_runtime_truth
```

The script expands to `cmake --fresh --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, and `ctest --preset test-stage1-runtime-truth-cpu-debug` over the exact P0 test list; it also validates the build metadata contract. The `integration_runtime_truth_ctest_labels` audit test verifies that P0 tests are registered and labeled and that feature-gated MPI/CUDA/Python/HDF5 app-smoke tests are absent from CPU-only builds rather than silently passing in a disabled feature configuration.

Reproducibility impact: CI/test registration only. No normalized config dumps, provenance payloads, restart/snapshot schemas, output naming, deterministic scheduling behavior, solver numerics, or SoA/hot-cold layouts are changed.

## Follow-up repair note for prompt 0-4 verification pass

The post-Codex verification pass found three concrete Stage-1 gaps and this repair applies them:

1. **TreePM target softening closure.** `TreeSofteningView` now has active-target species lanes and a target resolver that can fall back through active-target species, source-indexed active target overrides/species, and finally the scalar global fallback. `ReferenceWorkflow` now builds compact owned source softening/species sidecars and active-target sidecars instead of passing full-state particle sidecars into compact TreePM source views.
2. **Gas host-cell remap fail-fast behavior.** `rebuildLocalGasStateFromParticleIds` no longer maps invalid or removed black-hole/tracer host cells to cell zero. Invalid host remaps now throw, preventing silent attachment to the wrong gas cell after ownership compaction.
3. **Migration softening completeness.** `commitParticleMigration` now rejects inbound records that omit a softening value when the destination state carries a softening sidecar or override mask. Ownership migration records must carry the runtime-truth softening lane into such states.

The temporary particle-bound gas-cell contract remains intentionally named and asserted rather than replaced with AMR/moving-mesh semantics in this Stage-1 repair.

## Prompt 5-8 repair addendum — 2026-05-08

This follow-up pass tightened the post-P0 hardening layer from prompts 5 through 8.

- Restart exactness: reference workflow restart verification now uses a full runtime-state equality check across particle hot lanes, particle sidecar lanes, gas-cell identity and hydro fields, patch state, star/BH/tracer sidecars, species ledgers, module sidecar payloads, scheduler state, and distributed gravity metadata. A restart is no longer considered roundtrip-ok merely because counts, particle IDs, and rank ownership match.
- Payload integrity: module sidecar payload length is now included explicitly in the restart integrity hash before hashing payload bytes, preventing boundary ambiguity between adjacent module payloads.
- Restart tests: the restart roundtrip fixture now contains a tracer particle in addition to dark matter, gas, star, and black-hole state, and asserts the softening override mask and tracer sidecar lanes after restart and migration/restart.
- Transform fuzzing and active views: the hot/cold layout test now verifies workspace capacity preservation across `TransientStepWorkspace::clear()` and stale scatter rejection after particle and cell index-space mutation.
- Distributed ownership floor: local/global ownership identity summaries now include a particle-ID square-sum invariant in addition to count, sum, XOR, and local uniqueness. Reference workflow compares reduced rank-local identity against the pre-partition generated IC identity, not against its own reduced values.
- Ghost-vs-migration payload separation: `GhostExchangeBuffer` is explicitly restricted to ghost-refresh payloads when used with transfer descriptors. Ownership migration must use `ParticleMigrationRecord`; attempts to pack/unpack ghost payloads under ownership-migration intent now throw.

Validation performed in this pass:

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j1 --target test_unit_parallel_distributed_memory test_unit_hot_cold_sidecar_layout test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip test_integration_transform_fuzz_invariants
ctest --preset test-cpu-debug --output-on-failure -R "parallel_distributed_memory|hot_cold|transform_fuzz|restart_checkpoint"
```

The targeted CPU tests passed. The full Stage-1 CI helper was started with a fresh build but exceeded the interactive container timeout before completing the full repository build; no full green all-target CI claim is made here.
