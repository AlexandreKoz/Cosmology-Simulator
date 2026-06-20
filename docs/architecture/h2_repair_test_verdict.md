# H2 gas-cell identity repair — test verdict

Date: 2026-06-19  
Suggested PR title: `fix(h2): make gas-cell IDs authoritative across workflow scheduling and restart`

## Verdict

**Accepted for the exercised local CPU and single-rank HDF5 H2 contract.**

The repaired workflow uses `gas_cell_id` as the stable gas identity, keeps the
`GasCellIdentityMap` authoritative over compatibility mirrors, schedules gas cells through
an independent gas-cell scheduler, and persists that scheduler in restart schema v19.

**Not accepted as an MPI migration claim.** The MPI+HDF5+FFTW configure preset could not be
run because the environment lacks an MPI C++ installation (`mpi-cxx` / `MPI_CXX`). Therefore
real multi-rank migration, post-migration geometry rebind, and restart-after-MPI-migration
remain unexecuted gates rather than assumed capabilities.

## What changed

- Added controlled `SimulationState` identity replacement/restoration APIs that derive legacy
  gas-cell and patch-index mirrors from canonical identity records.
- Removed normal root-patch materialization dependence on reverse-building identity truth from
  mutable sidecar lanes.
- Split particle and gas-cell adaptive scheduler submission and lifecycle ownership; gas-cell
  time-bin lanes are derived mirrors keyed by stable `gas_cell_id`.
- Added restart schema v19 group `/gas_cell_scheduler`, with records keyed by
  `gas_cell_id`, integrity coverage, diagnostics, strict modern validation, and an explicit
  v14--v18 compatibility reconstruction path.
- Allowed legal patchless non-AMR/legacy cell states with zero patch sentinels while retaining
  strict coverage validation for AMR-backed states.
- Made the non-AMR Cartesian workflow path support 1D, 2D, and 3D uniform Cartesian layouts
  in arbitrary dense-row order without sorting canonical state.
- Added a full reference-workflow regression that compares canonical and shuffled storage order
  by `gas_cell_id`, including one parentless gas cell and two gas cells sharing one parent
  particle.
- Updated release manifest and current-schema documentation from stale v18 declarations to v19.

## Executed evidence

### CPU-only Debug

Built focused H2 targets and ran:

```text
ctest --test-dir build/cpu-only-debug -R \
'^(unit_gas_cell_identity_invariants|unit_time_integration|integration_gas_cell_identity_migration|integration_gas_cell_split_merge_remap|integration_hydro_decoupled_gas_cells|integration_reference_workflow)$' \
--output-on-failure
```

Result: **6 / 6 passed**.

### HDF5 Debug

Configured with `cmake --preset hdf5-debug`, built the exact targets, then ran:

```text
ctest --test-dir build/hdf5-debug -R \
'^(unit_gas_cell_identity_invariants|unit_time_integration|unit_restart_checkpoint_schema|integration_gas_cell_identity_migration|integration_gas_cell_split_merge_remap|integration_hydro_decoupled_gas_cells|integration_amr_production_hydro_integration|integration_amr_hydro_subcycling|integration_reference_workflow|integration_restart_checkpoint_roundtrip|integration_reorder_compaction_sidecars|integration_species_migration_invariants|integration_transform_fuzz_invariants|integration_restart_equivalence_harness|integration_restart_equivalence_multirate_bins|integration_restart_equivalence_amr_hydro|integration_restart_equivalence_amr_flux_registers|integration_restart_equivalence_amr_temporal_ghosts|integration_docs_scaffold|integration_release_readiness_artifacts)$' \
--output-on-failure
```

Result: **20 / 20 passed**.

Coverage includes identity migration, split/merge remapping, decoupled gas cells, scheduler
invariants, reordered/compacted sidecars, restart schema and round trip, direct-vs-restart
multirate continuation, local AMR hydro/subcycling, pending flux-register and temporal-ghost
restart equivalence, and release/documentation consistency.

## Deliberately not claimed

- MPI rank migration or live multi-rank restart equivalence.
- Regridding/migration during an open temporal AMR subcycle.
- Scientific convergence or cross-code validation beyond the existing CI-scale regressions.
- A full unrestricted CTest run: the project registers tests whose executables are intentionally
  not always built by focused target builds, so this verdict uses explicitly built, named
  acceptance tests rather than treating `Not Run` inventory entries as solver failures.
