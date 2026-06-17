# H3 AMR Hydro Closeout Audit

Date: 2026-06-16
Prompt ID: H3.9
Mode: Audit

## Verdict

H3 is **locally conservative at the implemented single-rank production AMR-hydro level**, but it is **not a full galaxy-formation credibility closeout** and should not be advertised as distributed/restart-complete AMR hydro.

The implemented production path is real: it uses `core::SimulationState` as persistent gas truth, builds patch-local hydro views, fills AMR ghosts, advances `HydroCoreSolver` over production patch geometry, generates flux-register records from actual Riemann face fluxes, scatters by stable `gas_cell_id`, and applies reflux automatically before the hydro stage returns.

The remaining closeout blocker is evidence breadth, not the local conservative mechanism itself: there is no dedicated production AMR hydro run/restart/run equivalence test. Restart schema v16 persists and validates patch geometry and gas identity lanes, and AMR patch migration has atomic payload tests, but H3 should keep restart-equivalence closure open until the production AMR hydro path itself is exercised across restart.

## Evidence Basis

Inspected before writing this audit:

- Repository contracts and ownership docs: `README.md`, `CONTRIBUTING.md`, `docs/architecture/overview.md`, `docs/architecture/developer_workflow_contract.md`, `docs/architecture/runtime_truth_map.md`, `docs/architecture/adr_runtime_truth_ownership.md`, `docs/repair_state_recap.md`, `docs/repair_open_issues.md`.
- H3 architecture/docs: `docs/architecture/h3_amr_hydro_coupling_plan.md`, `docs/architecture/h3_amr_hydro_production_closeout.md`, `docs/hydro_core_solver.md`, `docs/validation_ladder.md`, `docs/parallel_distributed_memory_contracts.md`, `docs/restart_checkpointing.md`, `docs/output_schema.md`.
- AMR/hydro implementation: `include/cosmosim/amr/amr_framework.hpp`, `src/amr/amr_framework.cpp`, `include/cosmosim/amr/amr_hydro_geometry.hpp`, `src/amr/amr_hydro_geometry.cpp`, `include/cosmosim/amr/amr_ghost_fill.hpp`, `src/amr/amr_ghost_fill.cpp`, `include/cosmosim/amr/amr_hydro_orchestrator.hpp`, `src/amr/amr_hydro_orchestrator.cpp`, `include/cosmosim/hydro/hydro_core_solver.hpp`, `src/hydro/hydro_core_solver.cpp`, `src/workflows/reference_workflow.cpp`.
- Restart/MPI impact surfaces: `include/cosmosim/io/restart_checkpoint.hpp`, `src/io/restart_checkpoint.cpp`, `include/cosmosim/core/simulation_state.hpp`, `src/core/simulation_state_species.cpp`.
- Tests and registration: `CMakeLists.txt`, `tests/unit/test_amr_*.cpp`, `tests/integration/test_amr_*.cpp`, `tests/integration/amr_hydro_validation_helpers.hpp`, `tests/integration/test_restart_equivalence_hydro_toy.cpp`, `tests/validation/`.

Discovery commands included:

```bash
rg -n "H3|AMR|amr|FluxRegister|Reflux|prolong|restrict|ghost|Sedov|shock|sync|subcycling|moving mesh|moving-mesh" src include tests docs CMakeLists.txt
rg --files tests | rg "test_amr_|validation|hydro|restart_equivalence_hydro_toy"
rg -n "advanceProductionAmrHydro|refineProductionPatchInSimulationState|derefineProductionPatchInSimulationState|fillAmrHydroGhostCells|populateAmrHydroFluxRegisterFaces|applyFluxRegistersToSimulationState" src include tests docs
rg -n "integration_amr_|unit_amr_|validation_hydro_classics|restart_equivalence_hydro_toy" CMakeLists.txt
rg -n "restart|checkpoint|AmrPatch|gas_cell_id|patch" src/io include/cosmosim/io tests/unit/test_restart_checkpoint_schema.cpp tests/integration/test_restart_checkpoint_roundtrip.cpp
```

## Conservative AMR Transfer

Status: **implemented and tested at scaffold and production-state levels**.

- `ConservativeTransfer::prolongateFromCoarse` splits volume-integrated `ConservedState` equally into fine cells, preserving mass, momentum x/y/z, and total energy for piecewise-constant refinement.
- `ConservativeTransfer::restrictToCoarse` sums child volume-integrated states back to a coarse cell.
- `PatchHierarchy::refinePatch` and `PatchHierarchy::derefinePatch` use those transfers for standalone AMR scaffold patches.
- `refineProductionPatchInSimulationState` and `derefineProductionPatchInSimulationState` repeat the conservation contract against persistent `SimulationState` gas rows, rebuild `GasCellIdentityMap`, update patch descriptors, and bump the cell index generation.

Coverage:

- `integration_amr_conservative_refine`
- `integration_amr_conservative_derefine`
- `integration_amr_production_hydro_integration`
- `integration_amr_synchronization_stress`

These tests check conservation of volume-integrated mass, momentum x/y/z, and total energy across refine/derefine. The implemented prolongation is piecewise-constant; there is no slope-limited conservative prolongation claim.

## Ghost Fill

Status: **production path, not test-only helper**.

`advanceProductionAmrHydro` builds `AmrHydroPatchGeometry` for every explicit patch descriptor, loads patch-local conserved storage from `SimulationState`, and calls `fillAmrHydroGhostCells` before solving. The ghost fill path handles:

- H1 physical boundary rules for periodic/open/reflective ghosts;
- same-level AMR ghost copy;
- coarse-to-fine ghost injection;
- fine-to-coarse coarse ghost averaging;
- stale remote ghost epoch rejection for imported/read-only ghost surfaces.

Coverage:

- `unit_amr_ghost_fill`
- `integration_amr_patch_boundary_consistency`
- `integration_amr_production_hydro_integration`
- `integration_amr_hydro_shock_tube`
- `integration_amr_hydro_sedov`
- `integration_amr_synchronization_stress`

Ghost rows are scratch-only and remain outside restart truth. Source same-level, coarse-fine, and remote imported ghosts are read-only from the perspective of the receiving patch.

## Flux Registers and Reflux

Status: **production path, automatically applied**.

`HydroCoreSolver::advancePatchActiveSetWithScratch` accepts a `HydroFluxRegisterSink`. During the Riemann face loop it calls the sink only when `HydroPatchGeometry::flux_register_faces` marks a face as coarse or fine. `populateAmrHydroFluxRegisterFaces` marks coarse-fine ghost faces before the sweep; `FluxRegisterAccumulator` area-weights coarse/fine contributions by register key; `advanceProductionAmrHydro` scatters patch-local states to `SimulationState` and then applies `applyFluxRegistersToSimulationState` automatically.

The reflux target is resolved through stable coarse patch/gas-cell metadata. Stale target mappings fail before correction. Incomplete or area-mismatched registers are skipped and counted; missing coarse or fine sides are not treated as implicit zero flux.

Coverage:

- `unit_amr_flux_register_generation` proves records come from actual hydro sweeps and configured Riemann fluxes.
- `integration_amr_reflux_conservation` proves automatic production reflux changes the coarse owner and preserves total state within tolerance.
- `integration_amr_production_hydro_integration` proves geometry, ghost fill, hydro sweep, flux-register generation, reflux, scatter, and identity consistency through one production path.
- `integration_amr_hydro_shock_tube`, `integration_amr_hydro_sedov`, and `integration_amr_synchronization_stress` exercise the production AMR path with CI-scale physical guards.

## Validation Registration

Status: **registered in CMake; command execution blocked in this environment**.

Registered H3-relevant CPU tests include:

- `unit_amr_refinement`
- `unit_amr_hydro_geometry`
- `unit_amr_ghost_fill`
- `unit_amr_flux_register_generation`
- `integration_amr_static_refinement_sync`
- `integration_amr_patch_boundary_consistency`
- `integration_amr_conservative_refine`
- `integration_amr_conservative_derefine`
- `integration_amr_reflux_conservation`
- `integration_amr_production_hydro_integration`
- `integration_amr_hydro_shock_tube`
- `integration_amr_hydro_sedov`
- `integration_amr_synchronization_stress`
- `integration_amr_patch_migration`

The validation ladder describes AMR hydro shock/Sedov/synchronization stress as CI-scale production guards. They are not publication-grade convergence decks.

## Restart and MPI Impact

Status: **documented and partially tested; H3 restart equivalence remains open**.

Restart:

- Restart schema is `cosmosim_restart_v16` in code.
- Patch geometry lanes, cell `patch_index`, gas-cell identity records, gas thermodynamic lanes, and patch ownership are serialized and integrity-hashed.
- Restart validation rejects invalid patch geometry coverage and malformed gas identity records.
- `test_restart_checkpoint_roundtrip` covers patch geometry and parentless gas-cell restart cases.

Blocker:

- `test_restart_equivalence_hydro_toy` is an H1 Cartesian hydro restart-equivalence test, not a production AMR hydro restart-equivalence test.
- No dedicated test currently proves `advanceProductionAmrHydro` plus AMR patch geometry/ghost/reflux state survives a run/restart/run split with equivalence to uninterrupted evolution.

MPI/patch migration:

- `AmrPatchMigrationRecord` carries patch descriptors plus complete gas-cell migration records.
- `commitAmrPatchMigration` updates patch descriptors and gas-cell rows atomically, rebuilds `GasCellIdentityMap`, and rejects stale gas ghosts by identity generation/ghost epoch.
- `integration_amr_patch_migration` tests atomic descriptor/gas sidecar payload behavior and stale ghost epoch rejection.

Limit:

- True distributed AMR ghost exchange and multi-rank production AMR hydro are not complete. Current remote ghost epoch machinery is a correctness seam and local/test-facing contract, not a full MPI AMR hydro closure.

## Limitations

- No AMR time subcycling is implemented or evidenced here. The production H3 path advances the selected active patch cells for the current hydro update and skips incomplete reflux registers rather than claiming subcycled fine/coarse synchronization.
- No cooling, star formation, stellar feedback, AGN feedback, metal enrichment, or other baryonic/subgrid coupling is closed by H3.
- No MHD, radiation transport, cosmic rays, conduction, viscosity, or chemistry is implemented by this stage.
- No moving mesh is implemented, validated, or recommended by this closeout.
- CI-scale shock/Sedov/stress tests are credible regression guards, not paper-grade validation or convergence campaigns.
- MPI AMR patch migration is represented by local/pseudo-rank payload contracts; full distributed AMR ghost exchange and multi-rank restart equivalence remain future work.
- Restart schema docs have some surrounding historical wording drift in older sections, but the implementation-facing H3 evidence depends on current code identity `cosmosim_restart_v16` and the AMR v16 patch geometry note.

## Future Moving-Mesh Seams

Moving mesh remains a future backend, not an H3 recommendation. If considered later, the required seams are:

- Keep `gas_cell_id` as stable cell identity independent of dense row order.
- Keep `parent_particle_id` as optional lineage metadata, not gas ownership truth.
- Preserve geometry adapters that feed `HydroPatchGeometry`-like solver inputs without promoting mesh-specific scratch to restart authority.
- Define face ownership, neighbor topology, ghost import epochs, and conservative correction ownership before solver coupling.
- Version restart/schema changes for any new topology, face, or mesh-generating-point fields.
- Add restart equivalence and migration tests before claiming a new backend is production-ready.

## Commands and Outcomes

Attempted:

```bash
cmake --preset cpu-only-debug
```

Outcome: **blocked**. PowerShell reported:

```text
cmake : O termo 'cmake' nao e reconhecido como nome de cmdlet, funcao, arquivo de script ou programa operavel.
```

Additional environment check:

```bash
Get-Command cmake -ErrorAction SilentlyContinue
```

Outcome: **blocked/not found**. No bundled workspace dependency runtime was configured in Codex Desktop.

Not run because configure was blocked by missing `cmake`:

```bash
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "hydro|amr|gas_cell|restart_equivalence_hydro_toy"
ctest --preset test-hdf5-debug --output-on-failure -R "restart_checkpoint|restart_equivalence_hydro_toy"
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R "amr|migration|parallel"
```

## Risk Register

| Risk | Status | Required next evidence |
| --- | --- | --- |
| Production AMR hydro local conservation | Closure candidate | Run CPU AMR/hydro targeted CTest bundle in an environment with CMake available. |
| AMR hydro restart equivalence | Open blocker | Add/run a dedicated production AMR hydro run/restart/run equivalence test that exercises patch geometry, ghost fill, flux registers, reflux, and stable gas-cell identity. |
| Distributed AMR hydro ghost exchange | Open future work | Implement multi-rank AMR ghost exchange with full identity/epoch contracts and MPI tests before any distributed AMR closure claim. |
| AMR subcycling | Not implemented | Define scheduler/reflux timing law and tests before claiming subcycled AMR. |
| Galaxy-formation physics coupling | Out of H3 scope | Separate validated stages for cooling, star formation, feedback, metals, and related source terms. |

## Allowed Claims

- Production AMR hydro has a real single-rank path through `SimulationState`, patch-local geometry, scratch ghost fill, live hydro sweeps, flux-register generation, scatter, and automatic reflux.
- Conservative refine/derefine is implemented for piecewise-constant transfers at scaffold and production `SimulationState` levels.
- CI-scale AMR hydro shock, Sedov, and synchronization stress guards are registered.

## Forbidden Claims

- H3 is fully closed for distributed AMR hydro.
- AMR hydro restart equivalence is proven for the production AMR path.
- AMR subcycling is implemented.
- The H3 evidence is publication-grade hydro validation.
- Cooling/star formation/feedback/MHD/radiation or moving-mesh production readiness is implied by H3.

## Reproducibility Impact

This audit changes documentation only. No solver numerics, config keys, restart schema fields, output names, HDF5 dataset names, scheduler behavior, MPI rank coordination, or persistent runtime ownership semantics were changed.
