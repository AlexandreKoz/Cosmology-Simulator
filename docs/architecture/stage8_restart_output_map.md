# Stage 8 restart/output maturity map

This note is the current Stage 8 restart/output implementation map after the schema-v14 repair pass. It supersedes earlier reconnaissance notes that referenced older restart schema versions.

## Current restart and snapshot entry points

- Restart writer: `cosmosim::io::writeRestartCheckpointHdf5()` in `src/io/restart_checkpoint.cpp`.
- Restart reader: `cosmosim::io::readRestartCheckpointHdf5()` in `src/io/restart_checkpoint.cpp`.
- Snapshot writer: `cosmosim::io::writeGadgetArepoSnapshotHdf5()` in `src/io/snapshot_hdf5.cpp`.
- Snapshot reader: `cosmosim::io::readGadgetArepoSnapshotHdf5()` in `src/io/snapshot_hdf5.cpp`.
- Workflow output/restart orchestration: `maybeWriteOutputs()` in `src/workflows/reference_workflow.cpp`.

Restart files use `cosmosim_restart_v14` and root `cosmosim_file_kind=restart_checkpoint`. Ordinary science snapshots use `gadget_arepo_v4` and root `cosmosim_file_kind=science_snapshot`. Readers reject the wrong file kind, so a particle-only snapshot cannot be treated as a complete restart.

## Persistent restart truth

Authoritative restart truth includes:

- `SimulationState` persistent lanes and sidecars: particles, particle sidecar, gas cells, cells, patches, species ledger, stars, black holes, tracers, module sidecars, and metadata.
- `IntegratorState`: timeline values, scale factor/redshift, KDK factors, boundary kind state, open-half-step flag, restart-safe flag, time-bin context, PM force-validity flag, and PM synchronization persistent state.
- `HierarchicalTimeBinScheduler`: current tick, max bin, bin assignment, next activation ticks, active flags, and pending bin transitions.
- `parallel::DistributedRestartState`: rank/decomposition metadata, PM grid/slab metadata, long-range field cadence metadata, and deterministic rebuild policy.
- `OutputCadencePersistentState`: output enabled flag, restart policy, due flags, last completed step, next snapshot step, and file stems.
- `StochasticPersistentState`: deterministic module contracts for current stochastic source modules, seeds, rank-local offsets, and last committed step index.
- Normalized config, provenance, restart diagnostics, and payload integrity hash.

`SimulationState::particles.time_bin` and `SimulationState::cells.time_bin` remain compatibility mirrors. Exact scheduler continuation authority is the serialized scheduler state, not those mirrors.

## Boundary and validation contract

Restart writes call `core::assertCanWriteCheckpointAtBoundary()` and reject unsupported open KDK half-steps, unsafe local substeps, and illegal PM synchronization boundaries. Unsupported half-step restart is deliberately not encoded in schema v14.

The loader validates file kind, schema name/version, required groups/datasets/attributes, scheduler array shape, output cadence state, stochastic state, distributed gravity state, diagnostics consistency, time-bin mirrors, ownership invariants, and the payload integrity hash. Missing scheduler or restart-critical fields are fatal.

## Test coverage map

Focused Stage 8 coverage includes:

- `unit_restart_checkpoint_schema`
- `unit_snapshot_hdf5_schema`
- `integration_restart_checkpoint_roundtrip`
- `integration_snapshot_hdf5_roundtrip`
- `integration_restart_equivalence_harness`
- `integration_restart_equivalence_dm_only`
- `integration_restart_equivalence_treepm`
- `integration_restart_equivalence_hydro_toy`
- `integration_restart_equivalence_multirate_bins`
- `integration_restart_equivalence_output_enabled`
- `integration_restart_equivalence_stochastic_sources`
- `integration_runtime_truth_ctest_labels`
- `integration_reference_workflow`

The TreePM, hydro, output, and stochastic equivalence tests now exercise production modules instead of only comparing metadata: TreePM uses `TreePmCoordinator::solveActiveSet()`, hydro uses `HydroCoreSolver::advancePatch()`, stochastic sources use `StarFormationModel`, and output-enabled equivalence compares due event-step sequences.

## Remaining risks and follow-up items

- FFTW-enabled TreePM restart equivalence should be run in an environment with FFTW3 development libraries.
- Stateful RNG engine support is intentionally absent; future modules using stateful engines must add explicit engine serialization.
- Full production resume-from-file workflow entry points remain future work; the current equivalence harness exercises production restart I/O and selected production physics modules.
