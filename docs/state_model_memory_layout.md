# Core Data Model and Memory Layout

This document defines the persistent and transient memory contract for simulation state.

## Persistent ownership (`SimulationState`)

`SimulationState` is the single ownership root for persistent run state:

- `ParticleSoa`: shared gravity-hot particle skeleton (`pos_*_comoving`, `vel_*_peculiar`, `mass_code`, `time_bin`).
  `ParticleSoa` is the only authoritative persistent owner of particle hot lanes in runtime state.
  Persistent acceleration arrays are intentionally absent; acceleration/force outputs are transient
  kernel/workspace products and are not restart truth.
- `ParticleSidecar`: shared metadata (`particle_id`, species tags, flags, rank ownership).
- `CellSoa`: gas-cell gravity skeleton (`center_*_comoving`, `mass_code`, `time_bin`, `patch_index`).
- `GasCellSidecar`: persistent gas-cell sidecar (`velocity_[xyz]_peculiar`, `density_code`, `pressure_code`, `internal_energy_code`, `temperature_code`, `sound_speed_code`) plus compatibility identity mirrors (`gas_cell_id`, `parent_particle_id`). `parent_particle_id == 0` mirrors an authoritative parentless identity record; it is not a gas-cell ID sentinel. Hydro transient reconstruction gradients are transient scratch, not restart truth.
- `GasCellIdentityMap`: authoritative in-memory gas-cell identity owner. It stores nonzero stable `gas_cell_id`, optional lineage-only `parent_particle_id`, stable `owning_patch_id`, and transient dense `local_cell_row` mappings with a generation counter. `SimulationState::refreshGasCellIdentityMapFromParticleBoundState()` materializes this map from legacy one-cell-per-gas-particle lanes, while general production identity validation accepts parentless cells and multiple rows with the same parent. Mutable hydro cell views capture the map generation so stale row mappings fail before scatter.
- `StarParticleSidecar`, `BlackHoleParticleSidecar`, `TracerParticleSidecar`: species-cold metadata blocks keyed by global particle index.
- `PatchSoa`: AMR patch descriptors and contiguous cell ranges.
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
for gravity-hot and hydro-active loops. Hydro active-view build/scatter paths validate dense
`GasCellIdentityMap` coverage, carry stable `gas_cell_id` plus transient `local_cell_row`, and scatter
by resolving the stable gas-cell ID through the current identity map. They do not require
`cell_index == gas particle row`; parent-particle mirrors are optional workflow compatibility mirrors.
This keeps solver kernels independent from full-state layout and
centralizes scatter/update semantics.

## H2 gas-cell ownership transition

H2 promotes `SimulationState::gas_cell_identity` as the authoritative in-memory identity map. Production identity
validation now requires only dense local-cell coverage, unique nonzero `gas_cell_id`, unique `local_cell_row`, and
sidecar mirror agreement. Parent lineage is optional and non-unique: AMR split rows may share one old parent, and
merged or newly created cells may be parentless.

Legacy particle-bound adapters remain available behind `legacyRequireParticleBoundGasCellContract(...)`
and the compatibility wrapper `requireParticleBoundGasCellContract(...)`:

- local `CellSoa` rows, `GasCellSidecar` rows, and local gas-particle rows have equal counts;
- `GasCellIdentityMap` covers dense local rows and matches the compatibility mirror lanes;
- `GasCellSidecar::gas_cell_id[cell]` and `parent_particle_id[cell]` mirror the parent gas particle ID;
- `parentParticleIdForGasCellRow(...)` returns `std::optional<std::uint64_t>` from the identity map;
- `gasCellRowForParticleId(...)` remains a legacy unique inverse and must only be used after the particle-bound contract;
- migration/compaction paths that still operate by gas particle ID must rebuild hydro fields by particle ID before those legacy adapters resume.

Hydro kernels consume cell-local velocity from `GasCellSidecar::velocity_[xyz]_peculiar`. When a gas cell has a
local parent particle, hydro writeback may update the particle mass and velocity compatibility mirrors through
explicit optional-parent lookup. If multiple cells share one parent, the parent mirror is written at most once during
the hydro store pass; cell-local gas state remains authoritative for every row. Parentless cells update only
cell-local gas state, so restart can preserve their hydro state without particle velocity access.

## Future decoupled gas-cell identity seam

`GasCellIdentityMap` can represent multiple gas cells with the same optional parent particle, gas cells without a
parent particle, stable patch ownership, and a transient `local_cell_row` mapping. H2.6 makes those layouts valid
production state for local identity validation and hydro active-view scatter without requiring local gas-particle
count to equal local gas-cell count. Remaining future work is distributed
gas-cell scheduler identity exchange and schema-versioned compatibility import for older restart files.

## Hot/cold ownership contract (Stage 6.2)

| Ownership class | Authoritative structures | Allowed contents | Forbidden in hot kernel views |
| --- | --- | --- | --- |
| Particle hot lanes | `ParticleSoa`, `GravityParticleKernelView` | `position_[xyz]_comoving`, `velocity_[xyz]_peculiar`, `mass_code`, transient `particle_index` scatter key | `particle_id`, species tag, owning rank, flags, SFC key, module state |
| Particle cold metadata | `ParticleSidecar` + species sidecars | IDs, species tags, ownership, flags, SFC, provenance-adjacent bookkeeping, optional per-species/module payloads | mutation from gravity/hydro hot-view scatter paths |
| Gas/cell hot lanes | `CellSoa`, `HydroCellKernelView` | `center_[xyz]_comoving`, `mass_code`, `density_code`, `pressure_code`, stable `gas_cell_id`, transient `local_cell_row` | parent-particle metadata, patch metadata, transient reconstruction gradients |
| Transient scratch | `TransientStepWorkspace`, reconstruction/work arrays | active-list mirrors, kernel gather/scatter buffers, temporary reconstruction/face work | restart truth or persistent ownership |

The unit contract tests now assert that gravity/hydro hot-view scatter updates only allowed hot lanes and leaves cold sidecars untouched.

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
- `ParticleSoaStorage`: reusable particle-oriented field pack (`pos_*`, `vel_*`, `mass`, `id`,
  `rho`, `u_int`) with typed access via `ParticleSoaField`.
- `gatherSpan` / `scatterSpan`: index-driven data movement for active kernels.

Schema/provenance implications:

1. In-memory schema version advances to `2` to reflect explicit species sidecars.
2. Canonical external naming remains unchanged in configuration and restart metadata.
3. The substrate keeps host-side semantics compatible with future device mirrors by using
   per-field contiguous arrays and explicit logical sizes.
4. `ParticleSoaStorage` is utility/test substrate only; it is not authoritative runtime
   particle truth and must not replace `SimulationState::particles` + sidecars.
