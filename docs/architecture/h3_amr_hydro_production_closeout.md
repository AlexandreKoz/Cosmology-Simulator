# H3 AMR Hydro Production Closeout

Date: 2026-06-17

## Scope of the accepted local production path

The accepted H3 production claim is intentionally narrow. CHUÍ/CosmoSim now has a local, synchronized AMR hydro path where `core::SimulationState` owns restart-authoritative gas and hydro state, patch descriptors provide explicit geometry, gas cells are addressed through stable identity records, ghost fill carries auditable status, and reflux applies only complete area-consistent records after stable-ID and patch-local validation.

The production path is selected only when production AMR coverage exists. In that case the workflow must not fall back to the old single global Cartesian toy hydro patch.

## AMR hydro restart equivalence

A dedicated HDF5-gated integration test now exists:

- source: `tests/integration/test_restart_equivalence_amr_hydro.cpp`
- CTest name: `integration_restart_equivalence_amr_hydro`
- CMake alias target: `integration_restart_equivalence_amr_hydro`

The test compares a direct AMR hydro continuation against a checkpoint/reload/continue path using the existing restart-equivalence harness. It exercises explicit patch geometry, gas identity coverage, ghost-fill/reflux safety, scheduler/integrator persistence, and final gas state comparison.

This proves local HDF5 AMR hydro continuation equivalence for the exercised deterministic synchronized-sweep scenario. It does not prove MPI restart, restart after AMR migration, subcycling restart, or deferred flux-register replay.

## H3.8 distributed boundary

H3.8 remains **local-only contract coverage** unless and until a real MPI test is added. The local migration test checks patch descriptor geometry lanes, gas sidecars, identity rebuild, row-order changes, stale ghost epoch rejection, and production AMR coverage after local commit. It does not send a patch/gas payload between two MPI ranks and does not run distributed AMR hydro after migration.

Do not describe the current H3.8 state as production distributed AMR hydro.

## AMR subcycling status

AMR hydro currently runs synchronized local sweeps. There is no recursive fine-level time subcycling, no Berger-Colella level hierarchy, and no level-by-level dt authority integrated with the production scheduler.

A real subcycling implementation will need a separate design and test stage covering:

1. scheduler-owned level dt hierarchy;
2. coarse/fine temporal interpolation and ghost timing rules;
3. deferred flux-register accumulation with area and dt coverage;
4. restart-safe pending register state;
5. MPI-safe synchronization and rollback behavior.

No placeholder flag should be treated as subcycling acceptance.

## Flux-register persistence status

Current reflux behavior is deliberately conservative: incomplete, missing-target, stale, or area-mismatched records are skipped/rejected instead of being applied. Skipped records are **not** persisted and are **not** replayed later.

This is safe for synchronized complete local sweeps where all required coarse/fine contributions are present by the end of the sweep. It is insufficient for subcycling or multi-step deferred synchronization. Future work must introduce a restart-safe pending register owner with patch IDs, gas-cell IDs, area coverage, dt coverage, generation/epoch metadata, and replay tests.

## Validation posture

The H3.7 shock tube, Sedov, and synchronization-stress tests are CI-scale guards. They check finite/positive state, AMR path usage, identity coverage, conservation bounds, and qualitative behavior. They are not scientific validation decks.

Future validation ladder:

1. fast CI smoke tests for AMR shock tube, AMR Sedov, and synchronization stress;
2. fixed-seed resolution/convergence studies with published tolerances;
3. cross-code comparison against established hydro/AMR solvers or analytic benchmarks;
4. conservation/restart/MPI equivalence on multi-rank runs;
5. larger science-readiness benchmarks.

## Latest tested evidence

Commands run in this hardening pass:

```bash
cmake --preset cpu-debug
cmake --build --preset build-cpu-debug --target cosmosim_amr -j2
ctest --preset test-cpu-debug \
  -R "amr_hydro|amr_synchronization|amr_patch_migration|restart_equivalence_amr_hydro" \
  --output-on-failure

cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug --target \
  test_integration_restart_equivalence_harness \
  test_integration_restart_equivalence_dm_only \
  test_integration_restart_equivalence_treepm \
  test_integration_restart_equivalence_hydro_toy \
  test_integration_restart_equivalence_multirate_bins \
  test_integration_restart_equivalence_output_enabled \
  test_integration_restart_equivalence_stochastic_sources \
  test_integration_restart_equivalence_amr_hydro -j2
ctest --preset test-hdf5-debug -R "restart_equivalence" --output-on-failure
```

Result: the targeted CPU AMR tests passed, and the HDF5 restart-equivalence suite passed 8/8 including `integration_restart_equivalence_amr_hydro`.

---

## Post-H3 update: local AMR subcycling and restart-safe pending registers

The production AMR hydro path now has two explicitly selectable sweep modes:

- synchronized local sweep, preserving the H3 behavior;
- local level-subcycled sweep through the AMR hydro orchestrator.

The subcycling mode is intentionally scoped. It is real in the sense that finer AMR levels are advanced with smaller timesteps and expected fine-substep coverage is tracked before reflux can be applied. It is not a claim that the global production scheduler owns an arbitrary Berger-Colella AMR timeline.

Pending flux-register state is now restart-authoritative. Incomplete records can be merged into `core::SimulationState::pending_flux_registers`, serialized through HDF5 restart schema v17, restored, validated by stable gas-cell identity, and completed after restart. The new `integration_restart_equivalence_amr_flux_registers` test exercises this direct-vs-restart path.

The validation posture is unchanged in the scientific sense: shock/Sedov/synchronization tests remain fast CI/regression guards. They are stronger than smoke tests, but they are not convergence studies or cross-code validation.
