# Reference workflow runtime path

This document describes the live config-driven runtime path assembled by `workflows::ReferenceWorkflowRunner` and exercised by the shipped executable target `cosmosim_harness`.

## What the runtime now does

A valid `param.txt` run goes through the following real path:

1. `cosmosim_harness <config.param.txt>` loads the typed frozen config with `loadFrozenConfigFromFile(...)`.
2. The runner validates the mode/build contract before stepping.
3. Initial conditions are dispatched from the authoritative config contract:
   - `mode.ic_file=generated` uses the in-repo generated IC path.
   - any other `mode.ic_file` is resolved relative to the config file and read through the HDF5 IC reader.
4. The live run loop uses `HierarchicalTimeBinScheduler` as the execution driver.
5. The orchestrator executes the canonical KDK stage sequence with real gravity and hydro callbacks.
6. Outputs, restart checkpoints, diagnostics, normalized config, and operational reports are written into the config-driven run directory.

## Runtime contract

The run directory is:

- `<output.output_directory>/<output.run_name>` in normal CLI usage.
- `<override_root>/<output.run_name>` only for tests/benchmarks that intentionally override the output root.

The runner honors these config fields directly:

- `mode.mode`
- `mode.ic_file`
- `output.output_directory`
- `output.run_name`
- `output.output_stem`
- `output.restart_stem`
- `output.snapshot_interval_steps`
- `output.write_restarts`
- `numerics.treepm_pm_grid`
- `numerics.treepm_asmth_cells`
- `numerics.treepm_rcut_cells`
- `numerics.treepm_assignment_scheme`
- `numerics.treepm_enable_window_deconvolution`
- `numerics.treepm_update_cadence_steps` (authoritative PM long-range refresh cadence in gravity kick opportunities; integer `>= 1`)

TreePM Phase-1 runtime mapping is explicit and auditable:

- `PmGridShape{N,N,N}` uses `N = numerics.treepm_pm_grid`
- `Δmesh = box_size / N`
- `r_s = treepm_asmth_cells * Δmesh`
- `r_cut = treepm_rcut_cells * Δmesh`
- Assignment/deconvolution are wired from typed config, not hidden workflow constants
- `treepm_assignment_scheme` maps directly to runtime PM assignment+gather (`cic` or `tsc`)
- `treepm_update_cadence_steps` is consumed by runtime PM cadence logic:
  - refresh long-range PM field every `N` gravity kick opportunities (`N = cadence_steps`)
  - reuse last long-range PM field between refreshes
  - default remains conservative (`N = 1`, refresh every kick)
- Reused PM long-range fields carry explicit provenance in runtime metadata (field version, build step, build scale factor), so each kick can be audited against the field it used.
- `r_cut` is resolved from typed config and drives explicit residual-traversal pruning in the TreePM residual path

## Canonical stage ordering

The workflow uses `core::StepOrchestrator` and preserves the canonical order:

1. `gravity_kick_pre`
2. `drift`
3. `hydro_update`
4. `source_terms`
5. `gravity_kick_post`
6. `analysis_hooks`
7. `output_check`

## Gravity, hydro, and scheduler ownership

- Gravity is executed through the live `gravity::TreePmCoordinator` callback.
- Hydro is executed through the live Godunov finite-volume callback using the hydro core solver, MUSCL-Hancock reconstruction, and HLLC fluxes.
- Active sets come from `HierarchicalTimeBinScheduler::beginSubstep()` and completed substeps are committed through `endSubstep()` before restart output is written.
- Restart payloads therefore serialize the scheduler state that actually drove the run.

## Artifacts emitted by the runtime

Every successful run writes the following core artifacts into the run directory:

- `normalized_config.param.txt`
- `profile.json`
- `profile.csv`
- `operational_events.json`

When output cadence conditions are met and HDF5 is enabled, the runtime also writes:

- `<output.output_stem>_NNN.hdf5`
- `<output.restart_stem>_NNN.hdf5` when `output.write_restarts=true`

## Failure transparency

The runner flushes the normalized config snapshot as soon as config loading succeeds. A top-level runtime wrapper then attempts to flush `operational_events.json` and profiler reports on both successful completion and runtime failure, so early first-run failures still leave a diagnostic trail whenever the filesystem remains writable.

## API migration note

`ReferenceWorkflowReport` includes PM cadence observability fields:

- `treepm_pm_grid`
- `treepm_update_cadence_steps`
- `treepm_long_range_refresh_count`
- `treepm_long_range_reuse_count`
- `treepm_cadence_records` (per-kick record including step, stage, field version, and refresh/reuse decision)
- `final_state_digest` (deterministic run-result checksum for reproducibility checks under fixed config/ICs)
