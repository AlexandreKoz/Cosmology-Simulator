# M2A Hydro Memory Architecture — 2026-09-03

## Status

Campaign M2A converts the production hydro/AMR execution path from full-population
primitive/reconstruction residency to deterministic active-batch staging governed by
the existing process `MemoryGovernor`.  It does not change the Euler equations,
adiabatic index, MUSCL-Hancock reconstruction, HLLC solver, source ordering, AMR
reflux semantics, floor policy, CFL policy, or restart schema.

## Authoritative state classification

The current `SimulationState` schema remains the restart and migration authority.
M2A does not create a competing conserved-state authority.

| State | Class | Lifetime / rationale |
| --- | --- | --- |
| `CellSoa.center_{x,y,z}_comoving`, `mass_code`, `patch_index` | canonical persistent | Cell geometry/mass and patch ownership/restart truth. |
| `CellSoa.time_bin` | diagnostic/compatibility mirror | Derived from scheduler authority; not timestep truth. |
| `GasCellSidecar.gas_cell_id`, `parent_particle_id` | canonical persistent identity/lineage | Restart/migration compatibility mirrors around the authoritative identity map. |
| `GasCellSidecar.velocity_{x,y,z}_peculiar` | canonical persistent | Cell-local restart truth, including parentless/split gas cells. |
| `density_code`, `pressure_code`, `internal_energy_code`, `metal_mass_code`, `temperature_code`, `sound_speed_code` | current restart-authoritative hydro state | Preserved because current restart/source/closure contracts consume these lanes. Removing derivable thermodynamic lanes requires a separately versioned restart-schema migration and is not hidden inside M2A. |
| `HydroConservedStateSoa` | phase-resident | Patch-local solver representation; never a second persistent authority. Six `double` lanes = 48 logical bytes/staged cell. |
| Primitive stencil state | active-batch scratch | Derived from step-start conserved state only for cells needed by the current face batch. |
| MUSCL predictor / left-right face states / Riemann fluxes | active-batch scratch | Bounded by selected active batch and recycled between face batches. |
| Sparse conservative deltas | phase-resident scratch | One delta per actually touched cell; committed only after the full face set so all reconstruction observes the same step-start state. |
| Source-term vectors | patch-local phase state | Production AMR derives metallicity/temperature/density inputs for the current patch instead of four full-population vectors. |
| Prepared AMR ghost shell | communication/synchronization | Surface-scaled step-start snapshot retained only through the hydro phase to preserve old-state neighbor semantics while patches scatter sequentially. |
| AMR topology / patch geometry | topology | Existing patch topology only. M2A introduces no permanent full-population hydro neighbor table; sparse/topology redesign is M2B scope. |
| Profile counters / high-water fields | diagnostic | Subsystem-scale metadata only. |

### Persistent logical width note

The current dense gas arrays contain 125 logical bytes/gas row before allocator
capacity and identity-map/hash overhead: 37 bytes in `CellSoa` (including the
1-byte non-authoritative time-bin mirror) plus 88 bytes in `GasCellSidecar`.
Excluding the time-bin mirror, the current canonical/compatibility lanes are
124 logical bytes/gas row.  This is a schema inventory, not a claim that every
thermodynamic lane is irreducible.  M2A deliberately avoids a restart format
change while eliminating the mandatory full-population solver-side duplicate.

## Active-batch policy

`numerics.hydro_active_batch_max_cells` is the single typed user-facing control.

- `0` means deterministic automatic maximum: 16,384 active cells.
- A positive value is a hard configured maximum.
- The existing process `MemoryGovernor` may reduce the selected batch using its
  current headroom and a conservative 4,096-byte/active-cell scratch reservation.
- Selection aligns to 8 cells when at least 8 cells fit; a one-cell batch remains a
  valid low-headroom fallback.
- No AIMD/predictive feedback loop is introduced in M2A.
- Reservation owner: `hydro.active_batch`, memory class `ScratchArena`.
- Distributed preparation failure is coordinated before the collective AMR path.

The 4,096-byte coefficient intentionally upper-bounds the dominant bounded
workspace terms: up to six face records per active cell, four primitive stencil
states per face, left/right face states, Riemann fluxes, and sparse batch indices.
It is admission policy, not a claim that every reserved byte is committed.

## Numerical-order contract

MUSCL-Hancock face reconstruction must see one common step-start state.  M2A
therefore batches only reconstruction/Riemann workspace.  Flux contributions are
accumulated into sparse touched-cell deltas and **no conserved cell is committed
until all active faces have been processed**.  This avoids the scientifically
invalid alternative of allowing later batches to reconstruct from already-updated
neighbors.

For production AMR, every active local patch first captures only its filled
step-start ghost shell.  The solver then loads, advances, and scatters one target
patch at a time.  This keeps the previous synchronized-old-state ghost contract
without retaining conserved volumes for all local patches simultaneously.

## Memory/scaling contract after M2A

Production hydro peak is intended to follow

```text
persistent SimulationState gas truth
+ max(one target patch conserved staging,
      target + geometrically touching source patches during ghost preparation)
+ surface-scaled prepared ghost snapshots
+ bounded active-batch primitive/face/Riemann scratch
+ sparse touched-cell conservative deltas
+ bounded distributed AMR payloads
```

rather than

```text
persistent SimulationState gas truth
+ all-local-patch conserved copies
+ full-population primitive cache
+ full-active-face left/right/flux arrays
+ four full-population source-property vectors
```

`HydroScratchBuffers` reports capacity/high-water for active-cell batches, face
batches, primitive stencil states, touched-cell metadata, and aggregate owned
capacity. `HydroAmrRuntime::memoryReport()` contributes hydro scratch, patch
workspace, prepared ghost snapshots, active metadata, and communication
high-water to the existing central memory report, so MPI max/mean reduction uses
the same report authority when MPI is available.

## Scope handoff

M2B owns broader sparse AMR topology/adjacency compaction. A later versioned
restart-schema campaign may investigate removal/rematerialization of persistent
`pressure_code`, `temperature_code`, and `sound_speed_code` only with explicit
restart compatibility, source/closure derivation, performance measurements, and
scientific equivalence evidence. M2A does not claim that schema compaction.

## M2A.1 acceptance-closure addendum — 2026-09-04

The post-M2A adversarial audit found that the local hydro batching architecture was sound but the directed distributed-AMR hydro handoff still materialized one `AmrPatchCellPayloadRecord` for every local cell, sent that full population-scale vector to each candidate peer, used a one-shot classic-MPI byte-count contract, and reported transferred bytes as communication memory high-water. It also found four full-grid derived source-property vectors on the supported single-rank fixed-grid hydro path.

M2A.1 closes those defects without changing hydro equations, reconstruction, Riemann semantics, CFL/floor policy, AMR reflux ordering, or restart schema:

- Directed AMR exchange now sends compact patch descriptors first, derives actual face-sharing remote interfaces, and asks a lazy provider to materialize only requested boundary-face cells. Edge/corner-only patch contacts do not request hydro cells.
- The old full-owned-cell directed-hydro payload builder is removed. Ownership migration remains a distinct commit protocol and is not replaced by sparse ghost state.
- Directed peer record transport uses deterministic bounded `MPI_Sendrecv` rounds. The default round ceiling is 16 MiB and every individual classic-MPI byte count is checked before narrowing to `int`. Logical payloads may therefore exceed `INT_MAX` without requiring a single illegal MPI count.
- Receive allocation and lazy payload preparation have pairwise readiness handshakes before payload traffic so a local preparation failure does not leave the peer waiting in the payload exchange.
- Remote AMR patch storage carries an availability mask. Unadvertised interior cells are not treated as valid hydro state; ghost fill throws if geometry selection ever requests a cell omitted by the sparse interface protocol.
- Communication telemetry separates traffic (`patch_*_payload_bytes`) from resident capacity high-water (`patch_cell_*_capacity_high_water_bytes` and `communication_workspace_high_water_bytes`).
- The fixed-grid cooling path no longer builds four `N_gas` derived vectors. Cooling can obtain per-cell density, hydrogen density, metallicity, and temperature through a read-only property provider derived on demand from authoritative state. The dense-span API remains as a compatibility surface and is regression-tested for numerical equivalence.
- `pressure_code`, `temperature_code`, and `sound_speed_code` remain in the restart-compatible schema for this campaign, but memory accounting now classifies them explicitly as `PersistentCache` rather than independent canonical conserved truth. Their eventual removal requires a versioned schema/performance campaign.

### Scaling consequence

For a rectangular patch with dimensions `nx × ny × nz`, a single face request now materializes `ny*nz`, `nx*nz`, or `nx*ny` cell records rather than `nx*ny*nz`. Multiple requested faces are unioned without duplicate edge/corner records. Thus the directed hydro-cell wire representation scales with advertised interface surface for ordinary patch adjacency, not the entire local gas population. Patch descriptors remain patch-scale control metadata; broader AMR-topology compaction remains M2B scope.

### New acceptance regressions

- `unit_parallel_distributed_memory` checks face-only interface planning and plans a logical AMR cell payload larger than `INT_MAX` into bounded 16 MiB rounds without allocating the multi-GiB payload.
- `integration_amr_distributed_remote_patch_boundary` runs with an unavailable remote interior cell and only the advertised boundary cell resident, proving the ghost/reflux path consumes sparse remote state.
- `unit_cooling_heating` checks on-demand source-property provider results against the dense compatibility path.
- `unit_memory_accounting` asserts the three derivable restart-compatible thermodynamic lanes report as `PersistentCache`.
