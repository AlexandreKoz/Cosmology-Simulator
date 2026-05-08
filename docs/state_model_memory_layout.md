# Core Data Model and Memory Layout

This document defines the persistent and transient memory contract for simulation state.

## Persistent ownership (`SimulationState`)

`SimulationState` is the single ownership root for persistent run state:

- `ParticleSoa`: shared gravity-hot particle skeleton (`pos_*_comoving`, `vel_*_peculiar`, `mass_code`, `time_bin`).
- `ParticleSidecar`: shared metadata (`particle_id`, species tags, flags, rank ownership).
- `CellSoa`: gas-cell gravity skeleton (`center_*_comoving`, `mass_code`, `time_bin`, `patch_index`).
- `GasCellSidecar`: thermodynamics and hydro reconstruction sidecar (`density_code`, `pressure_code`, gradients, etc.) plus the temporary particle-bound gas identity lanes (`gas_cell_id`, `parent_particle_id`). Gas particle ID is the stable identity anchor; gas cell row indices are local/transient and must not be serialized or consumed as universal identity.
- `StarParticleSidecar`, `BlackHoleParticleSidecar`, `TracerParticleSidecar`: species-cold metadata blocks keyed by global particle index.
- `PatchSoa`: AMR patch descriptors and contiguous cell ranges.
- `GasCellIdentityMap`: proposed, isolated identity seam for future AMR/moving-mesh decoupling. It is not yet stored as production state; see `docs/architecture/gas_cell_identity_map_rfc.md`.
- `SpeciesContainer`: explicit species accounting ledger.
- `ParticleSpeciesIndex`: explicit local/global species indexing map for branch-light species iteration.
- `StateMetadata`: schema/provenance-sensitive run metadata.
- `ModuleSidecarRegistry`: module-specific persistent payload blocks.

The ownership invariant is validated by `SimulationState::validateOwnershipInvariants()`.

## Species indexing and transfer packing

`ParticleSpeciesIndex::rebuild()` creates explicit species-local to global index mappings.
This supports:

- species-aware loops that avoid branching on irrelevant species fields,
- sidecar consistency checks,
- explicit `packSpeciesTransferPacket(...)` staging for MPI/device transfer paths.

For distributed ownership migration boundaries:

- `packParticleMigrationRecords(...)` exports all required authoritative hot lanes plus metadata sidecars (`particle_id`, `sfc_key`, species, flags, ownership) and species-specific sidecar payloads.
- `commitParticleMigration(...)` is the only state mutation point for ownership transfer:
  - removes outbound-owned rows and stale ghost/import rows,
  - appends inbound authoritative rows owned by `world_rank`,
  - rebuilds species sidecar indices and species counts.
- This keeps sidecar ownership auditable and prevents ambiguous partial ownership updates.

## Transient ownership (`TransientStepWorkspace`)

`TransientStepWorkspace` owns temporary compact arrays for active solver subsets and a
`MonotonicScratchAllocator` for scratch bytes. This workspace is explicitly resettable and does
not retain persistent simulation data.

## Active-set views

`buildParticleActiveView` and `buildCellActiveView` materialize compact, contiguous active buffers
from index lists. Returned views are `std::span`-based and stable for the lifetime of the workspace
storage.

`buildGravityParticleKernelView` / `scatterGravityParticleKernelView` and
`buildHydroCellKernelView` / `scatterHydroCellKernelView` provide explicit read/write compact views
for gravity-hot and hydro-active loops. Hydro active-view build/scatter paths validate the temporary
particle-bound gas-cell contract before solver kernels run; they may carry transient `cell_index` rows
for scatter, but they must not infer stable particle identity from those rows without the named contract
helpers. This keeps solver kernels independent from full-state layout and centralizes scatter/update semantics.

## Temporary gas-cell ownership contract

Stage-0 gas cells are particle-bound finite-volume carriers. Until AMR or moving-mesh decoupling is
implemented, any path that consumes local gas-cell rows must pass `requireParticleBoundGasCellContract(...)`:

- local `CellSoa` rows, `GasCellSidecar` rows, and local gas-particle rows have equal counts;
- `GasCellSidecar::gas_cell_id[cell]` and `parent_particle_id[cell]` equal the parent gas particle ID;
- `parentParticleIdForGasCellRow(...)` maps transient cell rows to stable gas particle IDs;
- `gasCellRowForParticleId(...)` maps stable gas particle IDs back to transient local rows after rebuild;
- migration/compaction must rebuild hydro fields by gas particle ID before hydro kernels resume.

This is an explicit quarantine of the current 1:1 assumption, not an AMR design. Future many-cells-per-particle
or cell-without-particle layouts should replace the helper implementations and contract, not add positional
assumptions to workflow or solver code.

## Future decoupled gas-cell identity seam

`GasCellIdentityMap` is the documented next-step seam for AMR/moving-mesh readiness. It can represent multiple
gas cells with the same optional parent particle, gas cells without a parent particle, stable patch ownership, and
a transient `local_cell_row` mapping. Its current use is intentionally limited to isolated validation tests so it
does not duplicate production authority or alter hydro/restart behavior. Promotion to production requires the
restart-schema migration plan in `docs/architecture/gas_cell_identity_map_rfc.md`, plus tests proving hydro state
remaps by stable `gas_cell_id` rather than particle index or row position.

## Layout policy and reorder contract

The solver-facing state is split into three ownership classes:

1. **Persistent hot arrays** (`ParticleSoa`, `CellSoa`) used by dominant gravity/hydro loops.
2. **Persistent cold sidecars** (`ParticleSidecar`, species sidecars, gas thermodynamics, module sidecars).
3. **Transient compact arenas** (`TransientStepWorkspace`) for active-set kernels and temporary staging.

`ParticleReorderMap` is the auditable reorder contract. It stores both `old_to_new_index` and
`new_to_old_index` maps and can be generated by:

- `ParticleReorderMode::kByTimeBin`
- `ParticleReorderMode::kBySfcKey`
- `ParticleReorderMode::kBySpecies`

`reorderParticles(...)` is the only allowed particle reorder API. It applies the same
permutation to all hot arrays and parent-aligned metadata lanes, then synchronizes species sidecars
through the explicit `SidecarSyncPolicy`:

- `kUseParentIndirection` (default): preserve sidecar row order and remap only `particle_index`
  through `old_to_new_index`, so the row remains attached to the same physical parent ID without
  moving cold payload lanes.
- `kMoveWithParent`: move the complete species sidecar row with the parent order, including all
  scalar payload lanes and fixed channel arrays, then remap `particle_index` exactly once to the
  post-reorder global index.

Ad-hoc partial sidecar reorders are forbidden. New built-in species sidecar fields must be added to
the typed row-movement helpers in `src/core/simulation_state.cpp` and covered by payload-identity
regression tests. Future module sidecars must expose a small, auditable reorder contract equivalent
to "move full rows" or "remap parent indices" rather than mutating `particle_index` metadata in
isolation.

`debugAssertNoStaleParticleIndices(...)` is a debug-oriented guard that throws if any species sidecar
index is stale after reorder or compaction passes. `debugAssertSpeciesSidecarOwnershipInvariants(...)`
adds the stricter row-ownership check: exactly one sidecar row for each eligible star/black-hole/tracer
particle and no species sidecar rows attached to ineligible particles.

## Reproducibility and schema implications

`StateMetadata` provides a deterministic key/value serialization surface for schema version, run
identity, snapshot/restart naming, step index, and normalized config hash fields.

Conservative assumptions:

1. Metadata serialization format is line-based key/value text and intentionally minimal.
2. Species tags are encoded as a bounded integer enum (0..4).
3. Patch-to-cell mapping uses contiguous ranges (`first_cell`, `cell_count`) for locality and future
   MPI packing.
4. Species sidecars reference global particle indices rather than duplicating IDs.
5. `ParticleSidecar::sfc_key` is sidecar metadata owned with parent particle rows and is valid for
   reorder/grouping even when no space-filling curve module is linked.

## Reusable SoA substrate (`soa_storage.hpp`)

`core/soa_storage.hpp` adds a reusable SoA substrate built around aligned contiguous field arrays,
field-keyed typed span access, and gather/scatter helpers for active kernels:

- `SoaFieldArray<T>`: aligned contiguous field storage with explicit `size`, `capacity`,
  `reserve`, `resize`, `swapErase`, and stable compaction.
- `ParticleSoaStorage`: canonical particle-oriented field pack (`pos_*`, `vel_*`, `mass`, `id`,
  `rho`, `u_int`) with typed access via `ParticleSoaField`.
- `gatherSpan` / `scatterSpan`: index-driven data movement for active kernels.

Schema/provenance implications:

1. In-memory schema version advances to `2` to reflect explicit species sidecars.
2. Canonical external naming remains unchanged in configuration and restart metadata.
3. The substrate keeps host-side semantics compatible with future device mirrors by using
   per-field contiguous arrays and explicit logical sizes.
