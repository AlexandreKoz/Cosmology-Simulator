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

## M2A.2 bounded remote hydro ghost acceptance closure — 2026-09-04

The post-M2A.1 adversarial audit found that M2A.1 had made the directed AMR
*wire* representation interface-scaled and classic-MPI-count safe, but had not
yet made the complete resident live set sparse. The receiver still allocated a
full logical receive vector, retained all-peer wire records, expanded sparse
boundary state into full remote patch-volume arrays and a full remote hydro
geometry/SoA, and deep-copied geometry to form a combined tracing collection.
Those ownership/lifetime defects could restore O(remote patch volume) memory
even though the transmitted records were only O(interface area).

M2A.2 changes the contract from merely bounded MPI calls to bounded physical
residency:

- directed AMR patch-cell transport streams record-aligned chunks through a
  producer/consumer interface; the complete logical receive payload is never
  allocated as one vector;
- the sender and receiver reuse bounded chunk buffers whose simultaneous
  capacity is limited by the negotiated transport window and admitted through
  the existing `MemoryGovernor`;
- received wire records are consumed directly into compact
  `AmrHydroSparseRemoteCell` stores keyed by patch-local cell offset, with no
  all-peer owning wire-record aggregate;
- sparse remote ghost state contains only advertised boundary cells. Missing
  remote interiors have no placeholder `gas_cell_id`, conserved-state, or
  availability lanes;
- ghost selection computes remote patch-local offsets from the immutable patch
  descriptor and reads the compact sparse store directly. Same-level,
  coarse-to-fine, and fine-to-coarse source selection do not construct remote
  solver geometry, remote internal faces, or a full remote
  `HydroConservedStateSoa`;
- the local AMR solver continues to build full geometry only for locally owned
  patches that it actually advances. The previous combined `trace_geometries`
  owning-copy stage is removed from distributed ghost-source use;
- sparse-store growth is admitted through `MemoryGovernor` as communication
  memory for the lifetime of the distributed hydro ghost phase. Transport
  workspace and sparse ghost residency therefore participate in the same
  process ceiling as hydro active-batch work;
- diagnostics keep traffic volume separate from resident capacity and report
  bounded transport workspace plus sparse ghost cell/metadata high-waters.

### M2A.2 scaling contract

For a remote patch with volume `nx*ny*nz` and advertised interface cell count
`N_interface`, production ghost-source residency is now approximately:

```
2 * bounded_transport_chunk
+ N_interface * sizeof(AmrHydroSparseRemoteCell)
+ O(remote_patch_count) descriptors/reservations
```

It no longer allocates `O(nx*ny*nz)` remote conserved state or remote solver
geometry solely to fill local ghosts. A regression compares patches with the
same 64x64 advertised face but depths 4 and 128 and requires identical sparse
cell-store capacity. Missing advertised source cells fail deterministically
rather than becoming default/zero state.

The restart-compatible `pressure_code`, `temperature_code`, and
`sound_speed_code` lanes remain `PersistentCache` in M2A.2. Their 24 B/gas-cell
persistent cost and any future rematerialization/schema migration remain a
separate performance/restart campaign and are not silently changed here.

## M2A.3 final distributed-hydro correctness addendum — 2026-09-04

Post-M2A.2 adversarial review found that the sparse remote boundary protocol
requested exactly one source layer per face. That is sufficient for same-level
and the current coarse-to-fine interpolation path, but it is insufficient when
a fine remote source supplies a coarse local ghost: one coarse ghost-cell width
covers multiple fine layers. The sparse selector also accepted whichever fine
records happened to be present, allowing an incomplete restriction stencil to
produce a plausible but incorrect average.

M2A.3 closes that correctness gap while preserving the M2A.2 bounded-memory
contract:

- `AmrPatchBoundaryCellRequest` now carries a compact per-face source depth.
  The planner derives the depth from the actual target/source cell-width ratio;
  a normal 2:1 fine-to-coarse interface requests two fine source layers, while
  same-level and coarse-to-fine requests remain one layer.
- Boundary payload counting, streamed production, and receive validation honor
  the per-face depth and union overlapping face slabs without constructing a
  full remote patch representation.
- Fine-to-coarse ghost selection independently derives the complete aligned
  source-cell index range for the target coarse ghost and requires every sparse
  offset before averaging. Missing second-layer state now fails deterministically
  instead of performing partial restriction.
- The transient distributed-hydro cell record now carries canonical
  `metal_mass_code`. The sparse receiver reconstructs the passive conserved
  metal lane using the same `metal_mass_code / mass_code` mass-fraction contract
  as local AMR hydro. This changes only the MPI transient wire record; HDF5
  snapshot/restart schema is unchanged.

The resulting residency law remains bounded communication workspace plus
`O(interface_area * required_normal_source_depth)` sparse source state and
`O(remote interfaces)` metadata. M2A.3 does not reintroduce dense remote patch
hydro state, remote solver topology, or full logical receive buffers.

Live multi-rank certification still requires the repository's authoritative MPI
preset on a host with MPI C++ development support. The source contract and
serial dense-vs-sparse regressions are intended to fail closed until that
external runtime evidence is available.

## M2A.4 constant-space directed-AMR transport planning addendum — 2026-09-04

Post-M2A.3 adversarial review found one remaining memory-scaling violation in
its otherwise bounded patch-cell stream. The production peer exchange called
`planDirectedAmrPatchCellTransferRounds()` for both directions and retained one
`std::size_t` entry per physical transport round. Under normal 16 MiB rounds
that metadata is small, but memory-governor pressure is explicitly allowed to
shrink the negotiated window. At a one-record window, planner residency grew
linearly with the complete logical interface population even though the actual
send/receive chunks were bounded.

M2A.4 removes that hidden scaling term without changing AMR geometry, ghost
selection, reconstruction, Riemann fluxes, reflux, passive-metal transport,
MPI wire records, restart schema, or numerical tolerances:

- production directed-AMR patch-cell streaming now uses the constant-space
  `DirectedAmrPatchCellTransferPlan`, containing only the logical record count,
  records per round, and arithmetic round count for each direction;
- per-round send and receive record counts are derived from logical offsets as
  the stream advances, so planner residency is O(1) with respect to interface
  population and transport-round count;
- the legacy/materialized round-size helper is retained for focused diagnostics
  and compatibility tests, but refuses to allocate more than 65,536 round
  entries and directs large callers to the constant-space planner;
- sender and receiver logical offsets are both checked at stream completion so
  the arithmetic plan remains fail-closed against dropped/overrun coverage.

### M2A.4 scaling evidence

On a 64-bit platform, the removed production metadata had the asymptotic form

```
2 * N_rounds * sizeof(std::size_t)
```

for simultaneous send/receive round vectors, before allocator over-capacity.
Thus a one-record window would require about 160 MB of planner entries at
10 million logical records and about 1.6 GB at 100 million records. The new
production plan keeps two fixed-size descriptors regardless of logical record
count. `unit_parallel_distributed_memory` includes a synthetic
100,000,000-record / one-record-per-round plan and completes without
materializing those 100 million entries; it also verifies that the diagnostic
materializer fails closed above its metadata cap.

The physical communication high-water remains the M2A.2/M2A.3 contract:
bounded reusable send/receive chunks plus sparse interface source state and
O(remote-interface) metadata. M2A.4 specifically removes planner metadata from
being another population-scaled live object under the smallest transport
windows.

### M2A.4 post-artifact scientific closure

The first broader CPU-debug inventory after the artifact-first checkpoint
exposed a pre-existing M2A.3 regression in the same acceptance boundary.
M2A.3 correctly required a complete fine-to-coarse restriction stencil, but
its selector assumed the entire coarse ghost volume was covered by one fine
source patch. Production octree refinement normally partitions that volume
across several child patches; for a 2:1 Cartesian refinement a coarse ghost may
require eight fine cells distributed across four tangential child patches.
The overly narrow source assumption therefore rejected valid shock-tube and
refine/derefine synchronization runs.

The closure now:

- intersects the coarse ghost volume with every finer local and sparse-remote
  candidate patch rather than requiring one patch to contain the whole volume;
- verifies that all contributors use one fine AMR level and consistent cell
  widths, maps every selected cell onto the expected refinement lattice, and
  rejects overlaps, gaps, off-lattice cells, or incomplete sparse payloads;
- orders contributors by geometric stencil ordinal before arithmetic averaging,
  so patch-vector/rank discovery order does not silently change the
  fine-to-coarse accumulation order;
- preserves the existing synchronization-time requirement for fine-to-coarse
  reads and does not change reflux, reconstruction, or Riemann semantics;
- updates the temporal-history negative fixture to stale the saved history
  fingerprint directly instead of corrupting live patch geometry, which would
  now (correctly) fail the stricter fine-to-coarse geometry gate first.

Focused regression coverage includes a coarse ghost tiled by four fine child
patches (eight fine source cells), both dense-local and sparse-remote source
forms, plus a missing-cell rejection. The previously red AMR temporal-history,
shock-tube, and synchronization-stress tests pass with the corrected selector.
The temporary coverage map is only one coarse-ghost restriction stencil in
size (`O(product(refinement_ratio_axis))`), not patch- or population-scale
persistent state.
