# H1 hydro closeout audit

Date: 2026-06-11
Scope: Stage H1 production fixed/patch hydro repair pass after Codex implementation.

## Verdict

H1 is repaired to the point where the targeted hydro, CFL, validation, and HDF5 restart gates are green in this repository snapshot. H2 may begin only with the known caveat that production hydro is still particle-bound at the workflow/state boundary; that dependency is explicitly H2 scope and is not hidden by this audit.

## Repairs applied

- Fixed the directional CFL unit test so the expected limiting axis matches the H1 contract: the timestep is the minimum of `cfl * dx_axis / (abs(v_axis) + c_s)` over x/y/z. The previous test data made x, not y, the limiting axis while expecting y.
- Fixed production periodic-boundary conservation for ghost-backed Cartesian patches. Periodic ghost cells remain reconstruction scratch, but flux deltas now update the wrapped authoritative real cell instead of being lost into non-authoritative ghost storage.
- Added a boundary-condition regression that advances a periodic ghost-backed patch and checks real-cell mass, momentum x/y/z, and total energy conservation.
- Strengthened the Noh CI validation setup so the converging inflow evolves long enough to measure actual compression instead of a nearly unevolved initial state.
- Fixed restart reader diagnostics so malformed hydro patch geometry is rejected with a specific `/state/cells/patch_index` or `/state/patches` error before the generic ownership-invariant error masks the root cause.

## Evidence commands

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug --target \
  test_unit_hydro_core_solver \
  test_unit_hydro_reconstruction \
  test_unit_hydro_riemann \
  test_unit_hydro_boundary_conditions \
  test_unit_time_integration \
  test_integration_hydro_sod_like \
  test_integration_hydro_conservation_periodic \
  test_integration_hydro_axis_symmetry \
  test_validation_hydro_classics -- -j2
ctest --preset test-cpu-debug --output-on-failure \
  -R "unit_hydro_|integration_hydro_|validation_hydro|unit_time_integration"
```

Result: 9/9 passed.

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug --target \
  test_unit_restart_checkpoint_schema \
  test_integration_restart_checkpoint_roundtrip \
  test_integration_restart_equivalence_hydro_toy -- -j2
ctest --preset test-hdf5-debug --output-on-failure \
  -R "restart_checkpoint|restart_equivalence_hydro_toy"
```

Result: 3/3 passed.

```bash
ctest --preset test-cpu-debug --output-on-failure \
  -R "validation_(integration|convergence|regression|phase2)|validation_hydro"
```

Result: 5/5 passed.

## Full-suite note

A full `ctest --preset test-cpu-debug --output-on-failure` was also attempted after a full CPU build. All executed tests through the H1/Hydro and runtime-truth sections passed except `integration_reference_workflow`, which aborts because the CPU-only preset has `COSMOSIM_ENABLE_HDF5=OFF` while that config-driven workflow requests snapshot/restart output. The HDF5 restart-specific gate above was run under `hdf5-debug` and passed.

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
- CFL enforcement: targeted unit gate repaired and passing.
- Conservation accounting: present and exercised by periodic hydro tests.
- Classical validation guards: Sedov, Noh, Gresho, Kelvin-Helmholtz, and Evrard toy guard registered and passing.
- Hydro restart equivalence and patch-geometry rejection: HDF5 gate passing.
