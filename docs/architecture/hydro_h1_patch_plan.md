# H1 hydro baseline audit and patch plan

Date: 2026-06-10

Mode: Audit, with a documentation-only repository artifact requested by prompt H1.0.

## Verdict

The live hydro path is a production-wired Godunov finite-volume chain, but the reference workflow currently presents it with a one-dimensional, x-axis geometry and a particle-bound gas-cell row contract. The hydro core API already carries face normals, owner/neighbor cells, full vector momenta, active face sets, scratch buffers, and source hooks; the limiting H1 risk is not the Riemann/update kernel shape, but the workflow geometry/reconstruction/CFL bridge and the gas-cell identity assumptions around `SimulationState`.

No prerequisite invariant from the prior runtime-truth work is missing for H1.0. The temporary production contract is still explicit: hydro consumers must pass `requireParticleBoundGasCellContract(...)` before mapping a gas cell row to a gas particle row, and restart validates gas-cell timestep mirrors through the parent gas particle. H1 follow-up patches must not bypass that contract until a deliberate schema/identity migration is designed and tested.

## Evidence basis

Inspected files requested by H1.0:

- `AGENTS.md`
- `docs/hydro_core_solver.md`
- `docs/reference_workflow.md`
- `docs/time_integration.md`
- `docs/restart_checkpointing.md`
- `docs/output_schema.md`
- `include/cosmosim/hydro/hydro_core_solver.hpp`
- `src/hydro/hydro_core_solver.cpp`
- `include/cosmosim/hydro/hydro_reconstruction.hpp`
- `src/hydro/hydro_reconstruction.cpp`
- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/core/time_integration.hpp`
- `src/core/time_integration.cpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/io/restart_checkpoint.cpp`
- `tests/unit/test_hydro_core_solver.cpp`
- `tests/unit/test_hydro_reconstruction.cpp`
- `tests/integration/test_hydro_sod_like.cpp`
- `CMakeLists.txt`
- `CMakePresets.json`

Additional orientation files inspected under the repo contract:

- `README.md`
- `CONTRIBUTING.md`
- `docs/architecture/overview.md`
- `docs/architecture/developer_workflow_contract.md`
- `docs/architecture/runtime_truth_map.md`
- `docs/architecture/adr_runtime_truth_ownership.md`
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`

Additional hydro evidence discovered by search:

- `tests/validation/test_validation_integration.cpp`
- `tests/validation/test_validation_convergence.cpp`
- `tests/integration/test_restart_equivalence_hydro_toy.cpp`

## Current H1 data flow

`HydroStageCallback::onStage(...)` in `src/workflows/reference_workflow.cpp` is the only live reference-workflow hydro update callback. It declares only `IntegrationStage::kHydroUpdate` through `k_hydro_stage` and its `StageContract` requires `GasCells | Particles | GhostCells`, mutates `GasCells | Particles`, uses `StageActiveSetFamily::kGasCells`, and is marked restart/output unsafe for the in-step stage.

Current data flow from `SimulationState` to `HydroConservedStateSoa` and back:

1. `HydroStageCallback::onStage(...)` returns early when `context.state.cells.size() == 0`.
2. It refreshes particle ghosts for `"hydro.godunov"` through `refreshParticleGhostsForSolver(...)`.
3. It calls `core::requireParticleBoundGasCellContract(context.state, "hydro callback")`.
4. It calls `rebuildGeometryIfNeeded(context.state.cells.size(), context.integrator_state.dt_time_code)`.
5. It builds active hydro faces from `context.active_set.cell_indices` with `buildActiveFaceView(...)`.
6. It resizes `m_conserved` to `context.state.cells.size()`.
7. For each cell row, it maps to a gas particle row with `core::gasParticleIndexForCellRow(context.state, cell_index)`.
8. It reads `gas_cells.density_code`, `gas_cells.pressure_code`, and `gas_cells.internal_energy_code`; pressure falls back to `(gamma - 1) * rho * internal_energy` when the stored pressure is nonpositive.
9. It reads particle velocities from `particles.velocity_{x,y,z}_peculiar[gas_particle_index]`.
10. It converts the primitive state with `HydroCoreSolver::conservedFromPrimitive(...)` and stores into `m_conserved`.
11. It snapshots imported ghost-cell conserved state with `snapshotHydroGhostConservedCells(...)`.
12. It builds `HydroUpdateContext` from `IntegratorState::dt_time_code`, `IntegratorState::current_scale_factor`, and an optional cosmological `hubble_rate_code`.
13. It builds `HydroSourceContext` with gravity cell acceleration spans from `GravityStageCallback::cellAccel{X,Y,Z}()`. Metallicity, temperature, and hydrogen number density are currently zero-filled vectors sized to `state.cells.size()`.
14. It invokes `HydroCoreSolver::advancePatchActiveSetWithScratch(...)` with `m_conserved`, `m_geometry`, active cells/faces, `MusclHancockReconstruction`, `HllcRiemannSolver`, `ComovingGravityExpansionSource`, `m_scratch`, and `m_primitive_cache`.
15. It restores imported hydro ghosts and exchanges conservative owner-side flux corrections with `executeBlockingHydroConservativeFluxCorrectionExchange(...)`, keyed by `parent_particle_id`.
16. For owned gas rows only, it converts updated conserved state back to primitive with `primitiveFromConserved(...)`.
17. It writes back `gas_cells.density_code`, `gas_cells.pressure_code`, `gas_cells.internal_energy_code`, `cells.mass_code`, `particles.mass_code`, and particle velocities for the parent gas particle.

Important ownership implication: `m_conserved`, `m_scratch`, `m_primitive_cache`, and `m_geometry` are callback-local transient state. Restart serialization is explicitly rooted at `RestartPersistentStateView::simulation_state`, so hydro scratch/geometry are not restart truth.

## Current 1D geometry chain

The reference workflow geometry is built in `HydroStageCallback::rebuildGeometryIfNeeded(...)`:

- Cache invalidation only keys on `cell_count`. It does not rebuild when `dt_time_code`, box geometry, or hydro boundary changes while cell count is unchanged.
- `HydroPatchGeometry::cell_volume_comoving` is set to `max(1.0, box_volume / cell_count)`.
- Interior faces are `cell_index -> cell_index + 1` only.
- Periodic closure adds `cell_count - 1 -> 0` only.
- Nonperiodic/open boundaries add one left ghost-like boundary on cell 0 with `normal_x = -1` and one right boundary on the last cell with `normal_x = 1`.
- All face areas are `1.0`.
- All non-boundary normals are `(1, 0, 0)`.
- H1.3 update: reconstruction CFL policy now carries per-axis `{dt/dx, dt/dy, dt/dz}` values while retaining
  `HydroReconstructionPolicy::dt_over_dx_code` as a scalar compatibility fallback.

The hydro core does validate general face normals and applies fluxes through `normal_{x,y,z}`, `area_comoving`, and `cell_volume_comoving` in `HydroCoreSolver::advancePatchActiveSetWithScratch(...)`, but the reference workflow currently feeds only an x-line face graph.

## H1.3 reconstruction status

`MusclHancockReconstruction` no longer infers a one-dimensional contiguous row. Reconstruction consumes explicit
`HydroFace::axis`, `owner_minus_cell`, and `neighbor_plus_cell` metadata populated by Cartesian geometry or hand-built
test geometry.

The current H1.3 contract is:

- MUSCL slopes are computed for `rho_comoving`, all three `vel_*_peculiar` components, and `pressure_comoving`.
- The predictor uses the face-normal velocity component and the axis-specific CFL entry from
  `dt_over_cell_width_code`.
- Missing explicit stencil rows are treated as boundary/degenerate zero-slope sides; production code does not recover
  stencils by row-difference inference.

This closes the x-only reconstruction caveat from the H1.0 audit without changing snapshot/restart schema or config keys.

## Current CFL hooks

Hydro CFL candidate submission is in `updateAdaptiveTimeBinsFromView(...)` in `src/workflows/reference_workflow.cpp`.

Current chain:

- `buildAdaptiveTimeStepCriteriaView(...)` calls `requireParticleBoundGasCellContract(...)` when cells exist.
- It stores `gas_particle_index_by_cell` by calling `gasParticleIndexForCellRow(...)` for each cell.
- `updateAdaptiveTimeBinsFromView(...)` estimates `cell_width = cbrt(box_volume / cell_count)` when cells exist.
- For each scheduler element with `element_index < cell_count`, it maps to `gas_index = gas_particle_index_by_cell[element_index]`.
- It computes `flow_speed = sqrt(vx^2 + vy^2 + vz^2)` from the gas particle velocity.
- It reads `sound_speed_code` from `view.gas_cells.sound_speed_code[element_index]`.
- It registers a CFL hook through `computeCflTimeStep({cell_width_code, flow_speed_code, sound_speed_code}, 0.4)`.
- It submits the CFL candidate both for the cell element (`"cell_hydro_cfl"`) and, when covered by scheduler size, for the parent gas particle (`"gas_particle_hydro_cfl"`).

Risk for H1.1-H1.3: solver geometry currently uses `box_size_x / cell_count` for reconstruction predictor while adaptive CFL uses `cbrt(box_volume / cell_count)`. A production patch must choose and document a single geometric cell-length authority for each H1 geometry mode, or deliberately separate face-normal predictor spacing from volume-equivalent CFL spacing with tests.

## Cell row == gas particle row assumptions and contract call sites

The current production contract does not allow naked `cell row == gas particle row` use. The intended seam is:

- `requireParticleBoundGasCellContract(...)`
- `parentParticleIdForGasCellRow(...)`
- `gasCellRowForParticleId(...)`
- `gasParticleIndexForCellRow(...)`

Relevant call sites found in the H1 path:

- `src/workflows/reference_workflow.cpp:305`: adaptive time-bin view construction checks the particle-bound gas contract before building `gas_particle_index_by_cell`.
- `src/workflows/reference_workflow.cpp:310`: adaptive time-bin view maps each cell row to its gas particle row.
- `src/workflows/reference_workflow.cpp:596`: `collectLocalGasCellRecords(...)` checks the contract before migration/compaction gas record collection.
- `src/workflows/reference_workflow.cpp:630`: `gasParticleIdByOldCellIndex(...)` checks the contract before old row to parent-particle-ID mapping.
- `src/workflows/reference_workflow.cpp:2578`: scheduler initialization checks the contract before assigning cell and parent gas-particle bins.
- `src/workflows/reference_workflow.cpp:2582`: scheduler initialization explicitly sets the parent gas particle bin through `gasParticleIndexForCellRow(...)`.
- `src/workflows/reference_workflow.cpp:2675`: drift callback checks before copying gas-particle positions to cell centers.
- `src/workflows/reference_workflow.cpp:3240`: gravity callback gas-cell acceleration sync checks before mapping cell rows to active gas-particle slots.
- `src/workflows/reference_workflow.cpp:3731`: hydro callback checks before converting `SimulationState` gas rows to conserved hydro rows.
- `src/workflows/reference_workflow.cpp:3738`, `3812`, and `3848`: hydro callback maps cell rows to parent gas particle rows for readback, flux-correction ownership, and writeback.
- `src/workflows/reference_workflow.cpp:3811`: hydro conservative flux correction records map `parent_particle_id` back to local gas-cell row with `gasCellRowForParticleId(...)`.
- `src/core/simulation_state_active_views.cpp:252` and `302`: hydro active-view build/scatter paths require the contract.
- `src/core/time_integration.cpp:1899`, `1962`, and `2029`: gas-cell time-bin mirror sync/match/assert paths require the contract before mapping cell mirrors through parent particles.
- `src/io/restart_checkpoint.cpp:136`: restart reader requires the contract before rebuilding `CellSoa::time_bin` mirrors from scheduler authority.

Conclusion: H1 should preserve this seam. If a future patch decouples cell rows from gas particle rows, it must introduce a new identity/remap contract and update restart/schema docs in the same patch.

## Geometry, reconstruction, CFL, restart, and test registration

Geometry:

- Public hydro geometry types live in `include/cosmosim/hydro/hydro_core_solver.hpp`: `HydroFace`, `HydroPatchGeometry`, and `HydroActiveSetView`.
- Core update validation and conservative flux accumulation live in `src/hydro/hydro_core_solver.cpp`.
- Reference-workflow geometry registration is local to `HydroStageCallback::rebuildGeometryIfNeeded(...)` in `src/workflows/reference_workflow.cpp`.

Reconstruction:

- Public reconstruction policy lives in `include/cosmosim/hydro/hydro_reconstruction.hpp`.
- Current production second-order reconstruction lives in `src/hydro/hydro_reconstruction.cpp`.
- The reference workflow owns the runtime policy instance as `HydroStageCallback::m_reconstruction`, initially constructed with
  zero CFL values and then replaced in `rebuildGeometryIfNeeded(...)` with per-axis `dt_over_cell_width_code`.

CFL:

- Public CFL API is `core::computeCflTimeStep(...)` in `include/cosmosim/core/time_integration.hpp`.
- Candidate source label is `TimeStepCandidateSource::kHydroCfl`.
- Reference-workflow submission is in `updateAdaptiveTimeBinsFromView(...)`.
- `docs/time_integration.md` documents CFL as `dt_CFL <= C_cfl * (Delta x / (|u| + c_s))`.

Restart:

- Restart schema is `cosmosim_restart_v14`.
- Hydro persistent lanes are not a separate hydro solver schema. They are part of `SimulationState`: cells, gas-cell identity, and gas thermodynamic fields.
- `src/io/restart_checkpoint.cpp` requires and writes `/state/gas_cells/{gas_cell_id,parent_particle_id,density_code,pressure_code,internal_energy_code,temperature_code,sound_speed_code}`.
- Restart integrity hashing includes the gas-cell hydro lanes and cell/particle time-bin mirrors.
- Hydro scratch, primitive cache, face geometry, profile counters, and ghost snapshots are transient and intentionally outside `RestartPersistentStateView`.

Tests and build targets:

- `cosmosim_hydro` is registered in `CMakeLists.txt` from `src/hydro/hydro_module.cpp`, `hydro_core_solver.cpp`, `hydro_reconstruction.cpp`, and `hydro_riemann.cpp`.
- `unit_hydro_reconstruction` is registered from `tests/unit/test_hydro_reconstruction.cpp`.
- `unit_hydro_riemann` is registered from `tests/unit/test_hydro_riemann.cpp`.
- `unit_hydro_core_solver` is registered from `tests/unit/test_hydro_core_solver.cpp`.
- `integration_hydro_sod_like` is registered from `tests/integration/test_hydro_sod_like.cpp`.
- `integration_hydro_axis_symmetry` is registered from `tests/integration/test_hydro_axis_symmetry.cpp`.
- `integration_restart_equivalence_hydro_toy` is registered and labeled `restart scheduler hydro gas equivalence integration`.
- Hydro benchmarks include `bench_hydro_kernels` and `bench_hydro_face_sweep`.
- Discovered validation paths include 1D hydro cases in `tests/validation/test_validation_integration.cpp` and `tests/validation/test_validation_convergence.cpp`.

## Validation gaps

Confirmed gaps for H1 follow-ups:

- No live reference-workflow test proves multidimensional hydro geometry.
- H1.3 adds unit coverage for non-x face reconstruction with transverse velocity slopes and z-axis predictor CFL.
- The geometry cache key includes `dt_time_code`, so reconstruction CFL is rebuilt when the timestep changes.
- No hydro-specific test proves face area/volume consistency for non-unit area faces.
- H1.3 reconstruction predictor policy uses per-axis Cartesian cell widths in the workflow.
- Restart equivalence covers hydro toy state persistence, but H1 changes that alter persistent gas fields or identity semantics will need explicit restart roundtrip/equivalence tests.
- Current validation decks are 1D hydro-oriented; they are useful regression anchors but not multidimensional production evidence.

Risk hypotheses for H1 follow-ups:

- Active-set face inclusion currently updates faces adjacent to active cells. With local timestepping, inactive neighbors can receive conservative deltas through active faces. This may be intended for conservative FV coupling, but H1.2 should state the active-neighbor update law and prove it with a focused test.
- Future CFL work should still reconcile scheduler CFL estimates with anisotropic Cartesian widths; reconstruction itself now uses
  per-axis widths and the geometry cache key includes `dt_time_code`.
- `cell_volume_comoving = max(1.0, box_volume / cell_count)` may be a defensive floor, but it is not a documented scientific geometry law for physical small volumes.

## H1.1-H1.7 patch order

### H1.1: Freeze the current hydro workflow contract

Expected files:

- `docs/hydro_core_solver.md`
- `docs/reference_workflow.md`
- `tests/unit/test_hydro_core_solver.cpp`
- `tests/integration/test_hydro_sod_like.cpp`

Patch intent:

- Add focused assertions around the existing 1D geometry/update assumptions without changing solver behavior.
- Document that the reference workflow currently feeds an x-axis 1D chain into a more general face-normal hydro core.
- Add a small regression for geometry/reconstruction cache staleness if production code changes are made.

Test target:

- `unit_hydro_core_solver`
- `integration_hydro_sod_like`

Stop conditions:

- Stop if tests cannot reproduce the current 1D Sod-like conservation baseline.
- Stop if the patch requires changing restart schema or gas-cell identity to explain the current behavior.

### H1.2: Repair geometry authority and cache invalidation

Expected files:

- `src/workflows/reference_workflow.cpp`
- possibly `include/cosmosim/hydro/hydro_core_solver.hpp`
- `tests/unit/test_hydro_core_solver.cpp`
- `tests/integration/test_hydro_sod_like.cpp`
- docs only if behavior or geometry semantics change.

Patch intent:

- Make geometry rebuild depend on all inputs that affect `HydroPatchGeometry` or reconstruction predictor policy, including `dt_time_code` and relevant box/boundary values.
- Replace undocumented `max(1.0, volume/cell_count)` behavior only if justified as a bug fix and covered by tests.
- Keep the topology 1D unless H1.2 explicitly scopes a documented multidimensional geometry.

Test target:

- A focused geometry-cache test.
- Existing hydro unit/integration tests.

Stop conditions:

- Stop if fixing geometry would change output schema or restart payloads.
- Stop if a multidimensional cell layout is needed but no ownership contract exists for row-to-neighbor mapping.

### H1.3: Make reconstruction geometry contract explicit

Expected files:

- `include/cosmosim/hydro/hydro_reconstruction.hpp`
- `src/hydro/hydro_reconstruction.cpp`
- `tests/unit/test_hydro_reconstruction.cpp`
- `tests/unit/test_hydro_core_solver.cpp`
- `docs/hydro_core_solver.md`

Patch intent:

- Either harden the existing MUSCL 1D fallback contract or add a deliberate axis/face-aware reconstruction policy.
- If adding multidimensional reconstruction, define how slopes are selected per face normal and how transverse velocities are predicted.
- Preserve conservative update and floors.

Test target:

- `unit_hydro_reconstruction`
- `unit_hydro_core_solver`

Stop conditions:

- Stop if the patch starts designing a moving mesh or Voronoi backend.
- Stop if non-contiguous cells need ghost-fill ownership that is not represented by current `HydroPatchGeometry`.

### H1.4: Align hydro CFL with geometry

Expected files:

- `include/cosmosim/core/time_integration.hpp`
- `src/core/time_integration.cpp`
- `src/workflows/reference_workflow.cpp`
- `tests/unit/test_time_integration.cpp`
- hydro integration tests as needed
- `docs/time_integration.md`

Patch intent:

- Define the hydro CFL length scale used by the live workflow.
- Ensure the reconstruction predictor and scheduler CFL candidate derive from compatible geometry data or document why they differ.
- Keep scheduler ownership intact; hydro submits candidates, scheduler remains authority.

Test target:

- `unit_time_integration`
- hydro unit/integration tests
- scheduler runtime-truth preset if scheduler behavior changes.

Stop conditions:

- Stop if a solver tries to mutate `state.*.time_bin` directly.
- Stop if CFL changes require scheduler schema changes without restart docs/tests.

### H1.5: Harden SimulationState hydro pack/unpack ownership

Expected files:

- `src/workflows/reference_workflow.cpp`
- `src/core/simulation_state_active_views.cpp`
- `tests/unit/test_gas_cell_identity_invariants.cpp`
- `tests/integration/test_restart_equivalence_hydro_toy.cpp`
- docs only if ownership behavior changes.

Patch intent:

- Keep `requireParticleBoundGasCellContract(...)` as the production gate.
- Make pack/unpack from `SimulationState` to conserved SoA auditable and testable.
- Ensure owned-vs-ghost writeback and conservative correction records remain keyed by parent particle ID.

Test target:

- gas identity unit tests
- hydro toy restart equivalence
- hydro workflow/integration test if available.

Stop conditions:

- Stop if a patch assumes `cell_index == particle_index` without the named contract/helper.
- Stop if a decoupled gas-cell identity migration is needed; that belongs to a separate schema/compatibility prompt.

### H1.6: Restart and output schema impact review

Expected files:

- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/io/restart_checkpoint.cpp`
- `docs/restart_checkpointing.md`
- `docs/output_schema.md`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_restart_equivalence_hydro_toy.cpp`

Patch intent:

- Only change restart/schema if H1.2-H1.5 introduced persistent hydro fields or changed gas identity semantics.
- Keep hydro scratch, geometry, profile counters, and primitive caches out of restart truth unless a versioned schema case is explicitly made.
- Preserve GADGET/HDF5 canonical snapshot dataset names.

Test target:

- restart schema unit tests
- restart roundtrip/equivalence tests
- HDF5 preset if HDF5 paths are touched.

Stop conditions:

- Stop if schema versioning, compatibility behavior, and docs cannot be updated in the same patch.
- Stop if a restart reader would silently accept ambiguous legacy hydro identity.

### H1.7: Validation deck and benchmark closeout

Expected files:

- `tests/validation/test_validation_integration.cpp`
- `tests/validation/test_validation_convergence.cpp`
- `validation/` decks if applicable
- `bench/hydro/bench_hydro_kernels.cpp`
- `bench/bench_hydro_face_sweep.cpp`
- `docs/validation_plan.md`
- `docs/profiling.md`

Patch intent:

- Add or update validation targets only after the production behavior is implemented.
- Preserve 1D Sod/contact/convergence baselines while adding the smallest multidimensional or geometry-specific validation needed by H1.
- Add profiling hooks only when hot-path geometry/reconstruction behavior changes.

Test target:

- relevant validation tests
- hydro benchmarks where performance-sensitive loops change.

Stop conditions:

- Stop if validation expectations are generic prose without numerical tolerances.
- Stop if the benchmark measures a path different from the production workflow path changed in H1.

## H1 closeout checklist

Before claiming H1 production fixed/patch hydro foundation:

- The live workflow geometry contract is documented and tested.
- Reconstruction behavior is either deliberately 1D with fail-fast/fallback evidence, or deliberately multidimensional with face/axis tests.
- CFL candidates use a documented geometry length scale and submit only through scheduler-owned candidate APIs.
- `SimulationState` hydro pack/unpack uses named gas identity helpers and never relies on naked row equality.
- Active-cell/active-face update semantics are documented and tested under local active sets.
- Restart impact is explicitly classified as no schema change or implemented as a versioned, documented schema change with compatibility behavior.
- Existing hydro unit/integration tests, restart hydro toy equivalence, and repo hygiene pass or have exact named blockers.

## Minimal validation for H1.0

H1.0 changes documentation only. Required validation floor:

```bash
./scripts/ci/check_repo_hygiene.sh
```
