# H3 AMR/Hydro Coupling Plan

Date: 2026-06-14
Prompt ID: H3.0
Stage: H3 - AMR hydro integration
Scope: design map only. No solver rewrite, moving mesh backend, subgrid coupling, config key, output name, or restart schema change is authorized by this document.

## Verdict

H3 can proceed from the H2 gas-cell identity closeout, but production AMR hydro is not implemented yet. The current code has two disconnected representations:

- `core::SimulationState::{cells,gas_cells,gas_cell_identity,patches}` is the restartable production owner for gas-cell identity, geometry attachment, thermodynamic lanes, patch ownership, and dense local cell rows.
- `amr::PatchHierarchy` / `amr::AmrPatch` currently owns a standalone `std::vector<amr::ConservedState>` and synthetic `CellMetrics`. That storage is useful for scaffold tests and benchmarks, but it must not become a second production owner for live hydro state.

H3 production work must make `SimulationState` the persistent owner for AMR hydro cells. AMR patch objects may provide topology, patch lifecycle, geometry views, and transient transfer/reflux workspaces, but persistent conserved hydro truth must map back to stable `gas_cell_id` rows in `SimulationState`.

H3.2 scaffold status: `PatchHierarchy::refinePatch` initializes child
`AmrPatch::conservedView()` values with a conservative piecewise-constant
prolongation from the parent patch. This preserves volume-integrated mass,
momentum components, and total energy for isolated AMR algorithm tests, while
leaving slope-limited conservative prolongation and production
`SimulationState` gas-cell remapping to later H3 stages.

H3.3 scaffold status: `PatchHierarchy::derefinePatch` restricts child
`AmrPatch::conservedView()` values back into parent cells before erasing child
patches. Parent scaffold gas-cell IDs remain stable, retired child gas-cell IDs
are archived in hierarchy diagnostics, parent metrics are refreshed from the
restricted state, and volume-integrated mass, momentum components, and total
energy are preserved across derefine. This remains an isolated AMR algorithm
contract; production derefine against `SimulationState` rows by stable
`gas_cell_id` remains future H3 work.

## Evidence Basis

Inspected files required by H3.0:

- `docs/architecture/h2_gas_cell_identity_closeout_audit.md`
- `include/cosmosim/amr/amr_framework.hpp`
- `src/amr/amr_framework.cpp`
- `src/amr/internal/amr_boundary.hpp`
- `include/cosmosim/hydro/hydro_core_solver.hpp`
- `src/hydro/hydro_core_solver.cpp`
- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/core/simulation_state.hpp`
- `src/io/restart_checkpoint.cpp`
- `tests/unit/test_amr_refinement.cpp`
- `tests/integration/test_amr_static_refinement_sync.cpp`
- `bench/amr/bench_amr_regrid_kernel.cpp`
- `bench/bench_amr_regrid_sync.cpp`

Additional targeted inspection:

- `src/core/simulation_state_active_views.cpp`
- `include/cosmosim/hydro/hydro_cartesian_patch.hpp`

Discovery commands used:

```bash
rg -n "class HydroStageCallback|HydroStageCallback|advancePatch|FluxRegisterEntry|RefluxSynchronizer|GasCellIdentity|refreshGasCellIdentityMapFromSidecarLanes|owning_patch_id|rowsForPatch|HydroCartesianPatchSpec|HydroPatchGeometry" src/workflows/reference_workflow.cpp src/core src/io include/cosmosim tests
rg -n "writeStateGroup|readStateGroup|gas_cell_identity|patches|gas_cells" src/io/restart_checkpoint.cpp
```

## Stop Conditions for H3 Follow-up Prompts

Stop immediately before production code changes if any of these are true:

1. `docs/architecture/h2_gas_cell_identity_closeout_audit.md` no longer says stable `gas_cell_id` is production-authoritative.
2. A proposed patch makes `AmrPatch::m_conserved` persistent restart truth or lets it diverge from `SimulationState` without an explicit view/import/export boundary.
3. A proposed live hydro sweep emits test-only `FluxRegisterEntry` records rather than records generated from production `HydroCoreSolver` face fluxes.
4. A proposed AMR regrid changes gas-cell rows without updating `GasCellIdentityMap`, `GasCellSidecar` mirror lanes, `PatchSoa`, scheduler identity/remap state, and restart validation in one commit boundary.
5. A proposed MPI patch migration moves cells using `parent_particle_id`, `CellSoa::time_bin`, or local row alone instead of full `gas_cell_id`, patch ownership, scheduler identity, and hydro fields.

## File and Symbol Map

| Area | Current symbol/file | Current role | H3 production role |
|---|---|---|---|
| AMR lifecycle | `amr::PatchHierarchy`, `amr::PatchDescriptor`, `amr::AmrPatch` in `include/cosmosim/amr/amr_framework.hpp` | Standalone patch tree with owned conserved/metric vectors | Patch topology and refinement lifecycle only; production cell state must be projected to/from `SimulationState` rows by `gas_cell_id` |
| AMR transfer | `amr::ConservativeTransfer` | Conservative prolong/restrict over `amr::ConservedState` | Reused as an algorithm after adding adapters between volume-integrated AMR conserved quantities and hydro density SoA rows |
| AMR reflux | `amr::FluxRegisterEntry`, `amr::RefluxSynchronizer` | Manually supplied scaffold corrections into `AmrPatch::m_conserved` | Flux-register data model should be extended/generated from live hydro face fluxes; final correction must apply to `SimulationState` gas rows by `gas_cell_id` |
| Hydro state | `hydro::HydroConservedStateSoa` | Density-style conserved lanes consumed by `HydroCoreSolver` | Transient patch-local SoA view/cache built from `SimulationState::gas_cells` and scattered back by `gas_cell_id` |
| Hydro geometry | `hydro::HydroPatchGeometry`, `hydro::HydroFace`, `hydro::HydroGhostCell` | Cartesian patch geometry, faces, ghosts, mutation rights | Per-AMR-patch geometry surface; same-level and coarse-fine ghost fill must populate these ghost slots before solver execution |
| Workflow coupling | `HydroStageCallback::onStage` in `src/workflows/reference_workflow.cpp` | Builds one row-order geometry, advances active cells, writes back gas sidecar state | Split active cells by `owning_patch_id`, advance each leaf patch, generate reflux records at coarse-fine faces, then reflux before parent/coarse state is made visible |
| Runtime owner | `core::SimulationState` in `include/cosmosim/core/simulation_state.hpp` | Persistent particles, gas cells, identity map, patch table | Single persistent owner for gas-cell hydro truth and patch row ownership |
| Identity | `core::GasCellIdentityMap` and `GasCellIdentityRecord` | Stable `gas_cell_id`, optional parent lineage, `owning_patch_id`, dense `local_cell_row` | Required key for AMR row remap, patch membership, split/merge, MPI migration, active views, restart |
| Restart | `src/io/restart_checkpoint.cpp` | Persists `/state/gas_cells`, `/state/gas_cell_identity`, `/state/patches` and validates patch coverage | Must remain the exact continuation owner; any H3 schema expansion requires versioning and docs updates |
| Existing tests | `test_amr_refinement`, `test_amr_static_refinement_sync` | Scaffold-level AMR lifecycle/manual reflux tests | Keep as AMR unit floor; add integration tests only when production hydro emits real flux registers |
| Existing benches | `bench_amr_regrid_kernel`, `bench_amr_regrid_sync` | Synthetic regrid/reflux work | Later update to measure production patch-row projection, ghost fill, flux register generation, and reflux apply costs |

## Persistent Ownership Decision

Persistent hydro truth belongs to `core::SimulationState`, not `amr::AmrPatch`.

Production H3 should treat:

- `SimulationState::cells` as the gravity-facing gas-cell skeleton: centers, mass mirror, timestep mirror, and patch row index.
- `SimulationState::gas_cells` as persistent hydro/cold gas lanes: stable `gas_cell_id` mirror, optional `parent_particle_id` mirror, velocity, density, pressure, internal energy, temperature, and sound speed.
- `SimulationState::gas_cell_identity` as the authoritative stable identity map: `gas_cell_id`, optional parent lineage, `owning_patch_id`, and dense local row.
- `SimulationState::patches` as the restartable local patch table: `patch_id`, level, contiguous `[first_cell, first_cell + cell_count)`, and `owning_rank`.

`AmrPatch::m_conserved` may remain for isolated AMR algorithm tests, but production hydro must not evolve it as persistent truth. If H3 keeps `AmrPatch::conservedView()` for algorithms, it must be a transient projection from `SimulationState` rows or explicitly documented as scaffold-only until replaced.

## Conserved State Mapping

`amr::ConservedState` is volume-integrated:

- `mass_code`
- `momentum_{x,y,z}_code`
- `total_energy_code`

`hydro::HydroConservedStateSoa` is density-style:

- `mass_density_comoving`
- `momentum_density_{x,y,z}_comoving`
- `total_energy_density_comoving`

For a cell with volume `cell_volume_comoving`:

```text
amr.mass_code = hydro.mass_density_comoving * cell_volume_comoving
amr.momentum_i_code = hydro.momentum_density_i_comoving * cell_volume_comoving
amr.total_energy_code = hydro.total_energy_density_comoving * cell_volume_comoving
```

Production H3 transfer/reflux code must state which representation it accepts. Recommended boundary:

- hydro sweeps consume density-style `HydroConservedStateSoa`;
- AMR transfer/reflux diagnostics store volume-integrated deltas because conservation checks are volume-integrated;
- adapters convert at the patch boundary using `HydroPatchGeometry::cell_volume_comoving`.

The source of row identity during all conversions is `GasCellIdentityRecord::gas_cell_id`, never raw dense row alone.

## Patch Geometry Contract

Each leaf AMR patch must expose a hydro geometry view equivalent to `hydro::HydroPatchGeometry`:

1. `PatchDescriptor::{origin_comov,extent_comov,cell_dims}` or `SimulationState::patches` plus cell centers define `HydroCartesianPatchSpec`.
2. Real cells are the dense local rows whose `GasCellIdentityRecord::owning_patch_id == patch_id`.
3. `HydroPatchGeometry::cellCount()` is the number of real cells on the patch.
4. `HydroPatchGeometry::totalCellStorageCount()` may include real plus ghost storage.
5. `HydroFace` lists same-level real-real faces, physical-boundary faces, same-level ghost faces, coarse-fine ghost faces, and imported-MPI ghost faces.
6. `HydroGhostCell` records source/owner cell references and `HydroGhostMutationRights`; imported MPI and coarse-fine ghosts are read-only scratch, not authority.

The current `HydroStageCallback` builds one geometry from row-order cell centers. H3 must replace this with a patch-local geometry builder that groups rows by `GasCellIdentityMap::rowsForPatch(patch_id)` / `PatchSoa` ranges and rejects stale identity generations before scatter.

## Refinement and Derefinement Order

### Refine / Prolong

Required order:

1. Freeze a regrid boundary at a restart-safe and scheduler-safe synchronization point.
2. Snapshot old identity map generation, scheduler identity records, and patch table.
3. Choose leaf patches/cells to refine using `CellMetrics` derived from `SimulationState` gas rows.
4. Allocate child gas-cell IDs as new stable nonzero `gas_cell_id` values.
5. Build a new `PatchSoa` and dense cell row layout.
6. Use `buildGasCellNewToOldRowMap(old_map,new_map)` to identify carried rows and `kInvalidGasCellRow` new rows.
7. Prolong parent volume-integrated conserved state into child rows using a conservative law; initialize cold gas sidecar mirrors from derived primitive state.
8. Assign the new `GasCellIdentityMap`, update `gas_cells.{gas_cell_id,parent_particle_id}`, `cells.patch_index`, patch ownership, and bump cell/identity generations.
9. Rebuild/remap scheduler cell entries by stable `gas_cell_id`; do not use old `CellSoa::time_bin` as authority.
10. Validate dense identity coverage and patch coverage before hydro resumes.

### Derefine / Restrict

Required order:

1. Derefine only at a synchronization point where child fluxes have been refluxed and no stale ghost correction is pending.
2. Sum child volume-integrated conserved values into the coarse target using `ConservativeTransfer::restrictToCoarse` semantics.
3. Choose the coarse `gas_cell_id` according to a documented rule: retain an existing parent/coarse ID when present, otherwise allocate a new ID and record consumed child IDs in diagnostics if such diagnostics are added.
4. Build the new dense row map, with removed child rows absent and coarse rows initialized by restriction.
5. Refresh `GasCellIdentityMap`, sidecar mirrors, `PatchSoa`, scheduler identity, and active views in the same commit boundary.

Derefinement must not silently average density lanes without volume weighting.
For the current scaffold hierarchy, `PatchHierarchy::derefinePatch` implements
this order locally: it restricts child `ConservedState` values first, keeps the
parent gas-cell IDs as the coarse identity, archives retired child gas-cell IDs
as non-persistent diagnostics, refreshes parent metrics, then erases the child
patches and marks the parent as a leaf.

## Ghost Fill and Boundary Ownership

### Same-level Ghost Fill

Same-level neighbors must fill `HydroGhostCell` slots from the authoritative owner row:

- same rank: copy from `SimulationState` row resolved by neighboring `gas_cell_id`;
- cross rank: import read-only ghost payload keyed by `gas_cell_id`, source patch ID, owner rank, and ghost epoch.

Ghost writes are scratch-only. They must never mutate `SimulationState` on the receiving side.

### Coarse-fine Ghost Fill

Coarse-fine ghost fill must be explicit:

- fine ghost from coarse: prolong coarse primitive/conserved state into fine ghost slots;
- coarse ghost from fine: restrict fine boundary state into coarse ghost slots where needed for reconstruction;
- record interpolation order, conservation expectation, and source epoch in the patch-boundary structure before live use.

Coarse-fine ghosts are read-only inputs to reconstruction and Riemann solves. The conservative correction across a coarse-fine interface belongs to the flux register/reflux path, not ghost storage mutation.

### Physical Boundaries

Existing `HydroBoundaryKind::{kPeriodic,kOpen,kReflective,kImportedMpi}` remains the starting boundary vocabulary. H3 may add AMR-specific boundary metadata, but it must not replace hydro boundary kinds without a public interface migration note.

## Flux Register Generation and Reflux Timing

Production `FluxRegisterEntry` records must be generated inside live hydro face accumulation, specifically in or immediately adjacent to the `HydroCoreSolver::advancePatchActiveSetWithScratch` face loop where `scratch.fluxes[active_face_slot]` is known and `geometry.faces[face_index]` identifies the face.

Required generated fields:

- coarse patch ID;
- coarse cell `gas_cell_id` or a stable coarse cell lookup that resolves to one before apply;
- coarse face flux from the coarse sweep;
- summed fine face fluxes from child sweeps covering the same coarse face;
- face area in comoving units;
- `dt_code`;
- level pair and face orientation diagnostics;
- scheduler/integrator step index for auditability.

The existing `amr::FluxRegisterEntry` only stores `coarse_patch_id` and `coarse_cell_index`. H3 production work should either extend it with stable `gas_cell_id`/face metadata or create a production companion type while keeping the scaffold type for current AMR tests.

Reflux timing:

1. Fill ghosts for the level about to advance.
2. Advance hydro leaf patches at the correct scheduler time boundary.
3. During face accumulation, capture real coarse-fine face fluxes into level-pair registers.
4. After all fine substeps covering the coarse step are complete, apply reflux once to authoritative coarse `SimulationState` rows by `gas_cell_id`.
5. Refresh primitive/cold mirrors and invalidate any active/hydro views that captured old cell or identity generations.
6. Only then allow derefine, restart write, output, or MPI migration that observes corrected coarse state.

## Restart and Schema Impact

Current restart already persists and validates:

- `/state/gas_cells/*`
- `/state/gas_cell_identity/{gas_cell_id,has_parent_particle,parent_particle_id,owning_patch_id,local_cell_row}`
- `/state/patches/{patch_id,level,first_cell,cell_count,owning_rank}`

H3.0 makes no schema change. Future H3 code changes require schema/version/docs updates if any of the following become persistent:

- patch geometry beyond current `PatchSoa`;
- AMR hierarchy parent/child links not derivable from `patch_id`, level, and cell ranges;
- pending flux registers across restart boundaries;
- ghost epochs needed for exact restart continuation;
- split/merge provenance needed to explain new `gas_cell_id` allocation.

If flux registers are not restartable, checkpoint writes must be blocked while any register is pending. That block belongs with the restart boundary contract, not in solver code.

## MPI Migration and Ownership Order

Patch/gas migration must move a complete ownership packet:

1. `GasCellMigrationRecord` fields, including `gas_cell_id`, optional parent flag/ID, `owning_patch_id`, destination row, identity generation, ghost hydro epoch, centers, mass, time-bin mirror, patch index, velocities, density, pressure, internal energy, temperature, and sound speed.
2. Patch table updates for moved or rebalanced patches, including `patch_id`, level, cell range, and `owning_rank`.
3. Scheduler identity records keyed by `gas_cell_id` for gas-cell activation and by particle ID for particles.
4. Ghost invalidation epoch and stale ghost records keyed by `gas_cell_id`.
5. Post-commit dense identity map refresh and sidecar mirror validation.

MPI migration must not infer gas ownership from parent particles. Parent IDs are optional lineage only; parentless and split cells are production-supported after H2.

## Prompt-by-prompt Implementation Sequence

### H3.1: State-backed AMR Patch View

Goal: introduce a production patch view that maps `PatchSoa` / `GasCellIdentityMap` rows into patch-local hydro state without changing solver numerics.

Stop unless: identity coverage, patch range coverage, and sidecar mirrors are valid.

Tests:

- patch view rejects stale identity generation;
- patch view resolves parentless and multi-cell-parent rows by `gas_cell_id`;
- scaffold `AmrPatch::m_conserved` remains non-production or adapter-only.

### H3.2: Patch-local Geometry and Ghost Fill

Goal: build `HydroPatchGeometry` per leaf patch, including same-level and physical ghost slots.

Stop unless: ghost mutation rights and owner/source epochs are represented.

Tests:

- same-level ghost fill copies from owner rows without mutating receiver truth;
- physical boundary ghosts match existing periodic/open/reflective behavior;
- row reorder invalidates cached geometry/view generations.

### H3.3: Coarse-fine Transfer

Goal: implement refine/prolong and derefine/restrict against `SimulationState` rows by `gas_cell_id`.

Stop unless: total mass, momentum, and total energy are conserved to tolerance across split/merge.

Tests:

- refine creates stable nonzero child IDs and updates patch ownership;
- derefine volume-restricts child conserved state;
- restart roundtrip preserves the new identity/patch layout.

### H3.4: Live Flux Register Production

Goal: generate flux register records from real `HydroCoreSolver` face fluxes.

Stop unless: records are emitted only for real coarse-fine interfaces and contain stable coarse/fine cell identity.

Tests:

- no fake/manual register path is used in production integration;
- summed fine flux and coarse flux produce expected reflux correction;
- conservation report includes reflux delta separately from local patch flux/source/floor deltas.

### H3.5: Reflux Apply to SimulationState

Goal: apply reflux corrections to authoritative coarse rows in `SimulationState` and invalidate stale views.

Stop unless: correction resolves target rows by `gas_cell_id`, not stale local row.

Tests:

- stale identity generation rejects reflux apply;
- corrected coarse state survives restart;
- active hydro scatter after reflux sees the new generation.

### H3.6: MPI Patch/Gas Ownership

Goal: move patches and gas cells across ranks with full scheduler, identity, ghost, and hydro payloads.

Stop unless: multi-rank exchange carries full gas-cell and scheduler identity records.

Tests:

- parentless gas cells migrate without particles;
- split cells sharing one parent migrate independently;
- stale ghost records cannot mutate owner truth;
- restart after migration preserves patch ownership and scheduler activation.

## Validation Floor

For this H3.0 plan-only patch:

```bash
./scripts/ci/check_repo_hygiene.sh
```

For future CPU code patches:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

If a future patch touches restart/HDF5 behavior, also run the HDF5 preset from `docs/build_instructions.md`. If it touches MPI patch/gas ownership, run the matching MPI preset; if `fftw3_mpi` or MPI C++ tooling is unavailable, report the exact configure command and missing dependency.

## Reproducibility Impact

This H3.0 document changes no runtime behavior, solver numerics, config semantics, output names, or restart schema. It records the required ownership and sequencing contract for future AMR hydro integration so that production work does not introduce duplicate conserved-state ownership or fake flux-register evidence.
