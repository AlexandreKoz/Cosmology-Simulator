# H1 hydro closeout audit

Date: 2026-06-12
Scope: Stage H1 production fixed/patch hydro repair and hardening pass after Codex implementation.

## Verdict

H1 is repaired to the point where the targeted hydro, CFL, validation, HDF5 restart, and HDF5 reference-workflow gates are green in this repository snapshot. H2 may begin only with the known caveat that production hydro is still particle-bound at the workflow/state boundary; that dependency is explicitly H2 scope and is not hidden by this audit.

This is a credible fixed/patch Cartesian hydro foundation, not a claim of paper-grade hydrodynamic validation. The classic tests remain CI-level guards.

## Repairs applied

- Fixed the directional CFL unit test so the expected limiting axis matches the H1 contract: the timestep is the minimum of `cfl * dx_axis / (abs(v_axis) + c_s)` over x/y/z.
- Fixed production periodic-boundary conservation for ghost-backed Cartesian patches. Periodic ghost cells remain reconstruction scratch, but flux deltas now update the wrapped authoritative real cell instead of being lost into non-authoritative ghost storage.
- Hardened active-set flux application: the solver now applies flux deltas to every authoritative real cell touched by an active face, while source terms remain restricted to the explicitly active cells. This prevents active-face updates from dropping conservative flux into inactive real neighbors or periodic wrapped cells.
- Added `testActivePeriodicGhostFluxUpdatesWrappedRealCell`, a regression that advances a single periodic ghost face with only the boundary owner active and verifies that the wrapped real cell changes and that mass, momentum x/y/z, and total energy are conserved over the touched real cells.
- Strengthened the Noh CI validation setup so the converging inflow evolves long enough to measure actual compression instead of a nearly unevolved initial state.
- Fixed restart reader diagnostics so malformed hydro patch geometry is rejected with a specific `/state/cells/patch_index` or `/state/patches` error before the generic ownership-invariant error masks the root cause.
- Materialized a single root hydro patch for generated/legacy states that have gas cells but no patch descriptors, so restart payload validation and HDF5 reference workflow output see explicit `/state/patches` coverage instead of an implicit empty-patch state.
- Tightened `HydroStageCallback` Cartesian geometry reconstruction: if gas-cell centers form a complete 3D Cartesian product but the dense rows are not row-major, the workflow now fails loudly instead of silently falling back to a generated near-cubic patch geometry. The near-cubic fallback remains available for generated/toy states whose centers do not describe a full rectilinear patch.

## Evidence commands

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug --target \
  test_unit_hydro_core_solver \
  test_unit_hydro_reconstruction \
  test_unit_hydro_riemann \
  test_unit_hydro_boundary_conditions \
  test_unit_time_integration \
  test_integration_hierarchical_timestep_regression \
  test_integration_hydro_sod_like \
  test_integration_hydro_conservation_periodic \
  test_integration_hydro_axis_symmetry \
  test_validation_hydro_classics \
  test_integration_reference_workflow -j2
ctest --preset test-cpu-debug --output-on-failure \
  -R "unit_hydro_|integration_hydro_|validation_hydro|unit_time_integration|hierarchical_timestep"
```

Result: 10/10 passed.

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug --target \
  test_unit_restart_checkpoint_schema \
  test_integration_restart_checkpoint_roundtrip \
  test_integration_restart_equivalence_hydro_toy \
  test_integration_reference_workflow -j2
ctest --preset test-hdf5-debug --output-on-failure \
  -R "restart_checkpoint|restart_equivalence_hydro_toy|^integration_reference_workflow$"
```

Result: 4/4 passed.

```bash
bash ./scripts/ci/check_repo_hygiene.sh
```

Result: passed.

## Full-suite note

`integration_reference_workflow` still aborts under the CPU-only preset when it is run directly because that preset has `COSMOSIM_ENABLE_HDF5=OFF` while the config-driven workflow requests snapshot/restart output. The same test passes under the HDF5 preset, and the HDF5 restart-specific gate is green.

## Remaining H2 blockers

These are not H1 failures, but they remain hard prerequisites for H2:

- `HydroStageCallback::onStage` still calls `core::requireParticleBoundGasCellContract`.
- Hydro workflow load/store still uses `core::gasParticleIndexForCellRow` and mirrors gas-cell velocity/mass through parent gas particles.
- Flux-correction lookup remains parent-particle-id based.
- Parentless and multi-cell-parent gas-cell cases are not production runtime paths yet.

## H3 blockers

These are intentionally outside H1:

- AMR patches are not yet first-class hydro geometry providers.
- Conservative prolongation/restriction for live AMR hydro remains H3 work.
- Coarse-fine ghost fill, production flux-register generation, and automatic refluxing are not complete H1 responsibilities.

## H1 status checklist

- 3D Cartesian hydro geometry: present and tested.
- Ghost zones and boundary conditions: present and tested, including periodic ghost-backed conservation repair.
- Axis-aware reconstruction: present and tested for x/y/z symmetry.
- CFL enforcement: targeted unit/integration gates repaired and passing.
- Conservation accounting: present and exercised by periodic and active periodic ghost hydro tests.
- Classical validation guards: Sedov, Noh, Gresho, Kelvin-Helmholtz, and Evrard toy guard registered and passing.
- Hydro restart equivalence and patch-geometry rejection: HDF5 gate passing.
- HDF5 reference workflow output: passing with explicit root patch materialization.
