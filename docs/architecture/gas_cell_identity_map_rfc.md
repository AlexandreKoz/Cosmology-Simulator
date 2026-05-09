# RFC: GasCellIdentityMap seam for decoupled gas ownership

Status: Proposed infrastructure seam (P2 design-enabling upgrade)  
Scope: identity/mapping API only; no hydro solver, AMR, moving-mesh, or restart schema behavior changes.

## Context

The current production gas contract is intentionally transitional: each local gas cell is a finite-volume
carrier bound 1:1 to a local gas particle, and `GasCellSidecar::gas_cell_id` plus
`GasCellSidecar::parent_particle_id` are refreshed from canonical gas-particle order. That contract is safe
for current hydro and restart paths because row counts and row order are validated before hydro active views
or migration remaps consume gas-cell rows.

AMR and moving-mesh evolution need a different identity model. Gas cells must be able to split, merge, move
between patches, or exist without an active parent particle while preserving stable cell identity and a local
row mapping for SoA storage. This RFC defines the next API seam without replacing the current P0 production
contract.

## Proposed seam

`core::GasCellIdentityMap` is an isolated, validated map of `core::GasCellIdentityRecord` rows. Each record
contains:

- `gas_cell_id`: stable gas-cell identity. This is the primary key and must not be interpreted as a particle
  index or local row.
- `parent_particle_id`: optional parent particle identity. It can be absent for mesh-generated cells and can
  repeat across multiple gas cells during refinement, so it is not a uniqueness authority.
- `owning_patch_id`: stable patch identity for AMR/moving-mesh ownership decisions. It names the owner; it is
  not a patch row index.
- `local_cell_row`: transient local row into `CellSoa`/`GasCellSidecar` storage for the current rank/layout.

The current implementation validates nonzero unique `gas_cell_id` values and unique `local_cell_row` values,
rebuilds lookup tables atomically during `assign(...)`, exposes read-only records, supports lookup by stable gas-cell
ID or local row, and carries a monotonic `generation()` counter. The generation counter is already present so a
future production promotion can make hydro views and restart writers reject stale row mappings exactly as particle and
cell active views do today. It is deliberately not stored inside `SimulationState` yet, so it cannot shadow the
authoritative production gas contract or affect solver paths.

The seam also provides `coversDenseLocalRows(cell_count)` / `requireCoversDenseLocalRows(...)` and
`buildGasCellNewToOldRowMap(old_map, new_map)`. These utilities are intentionally small but not toy behavior: they
let tests represent cell split/reorder/new-cell cases and prove that remapping is keyed by stable `gas_cell_id`, with
new cells marked by `kInvalidGasCellRow` so callers must deliberately initialize conserved hydro state.

## Ownership and invalidation rules

The future production owner should be a single identity map owned by `SimulationState` or a dedicated cell
ownership component under `core/`. No hydro module should maintain an independent duplicate map.

When the seam is promoted to production, these invalidation rules must be enforced before any solver use:

1. Any operation that changes cell count, local row ordering, AMR patch ownership, moving-mesh connectivity, or
   rank ownership invalidates `local_cell_row` mappings and must advance the identity-map generation.
2. Such operations must rebuild the identity map before hydro kernels or restart writers consume gas rows.
3. Hydro remapping must key conserved state by `gas_cell_id` first and only use `local_cell_row` for local SoA
   scatter/gather after `requireCoversDenseLocalRows(...)` succeeds.
4. Parent-particle relationships are lineage metadata only. Multiple gas cells may share one parent particle,
   and a gas cell may have no parent particle.
5. Patch ownership changes must update `owning_patch_id`, row maps, and rank transfer payloads in the same
   auditable commit as row remapping.
6. `gas_cell_id == 0` is reserved as an uninitialized sentinel and must not be accepted into a production map.

## Migration plan

1. **Current state (this RFC):** keep `requireParticleBoundGasCellContract(...)` as the production gate. Add a
   compile/unit-tested `GasCellIdentityMap` shape that can represent non-1:1 mappings in isolation.
2. **Schema planning:** before serializing the map, define a new restart schema version and compatibility policy
   in `docs/output_schema.md` and migration notes. No restart payload changes are made by this RFC.
3. **Dual-read transition:** when implemented, restart readers should accept the legacy particle-bound lanes and
   materialize a `GasCellIdentityMap` where `gas_cell_id == parent_particle_id` for legacy files. New files should
   write explicit `gas_cell_id`, optional `parent_particle_id`, `owning_patch_id`, and `local_cell_row` or a
   restart-safe row-order equivalent.
4. **Hydro remap transition:** introduce remap tests that split one parent particle into multiple gas cells,
   move a cell between patches, and verify conserved hydro fields are associated by `gas_cell_id` rather than
   particle index or row position.
5. **Contract replacement:** only after restart compatibility, migration tests, and hydro remap tests are green
   should production paths replace the temporary 1:1 gas contract.

## Restart and reproducibility implications

This RFC does not change snapshot, restart, or provenance schema. Therefore existing restart determinism and HDF5
schema stability are preserved.

A future schema change must be explicit because decoupled gas identity affects restart reproducibility:

- stable gas-cell IDs must be serialized deterministically;
- optional parent-particle IDs must preserve absence vs zero-valued IDs;
- patch ownership must be serialized by stable patch ID, not local patch row;
- local row order must either be reproduced deterministically or reconstructed from a documented stable ordering;
- migration readers must document how legacy particle-bound restarts materialize the new map.

## Non-goals

- No hydro solver rewrite.
- No AMR or moving-mesh implementation.
- No snapshot/restart schema mutation.
- No change to current production gas-particle/gas-cell equality checks.
