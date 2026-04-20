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
- `numerics.treepm_pm_grid_nx`
- `numerics.treepm_pm_grid_ny`
- `numerics.treepm_pm_grid_nz`
- `numerics.treepm_asmth_cells`
- `numerics.treepm_rcut_cells`
- `numerics.treepm_assignment_scheme`
- `numerics.treepm_enable_window_deconvolution`
- `numerics.treepm_update_cadence_steps` (authoritative PM long-range refresh cadence in gravity kick opportunities; integer `>= 1`)

TreePM runtime mapping is explicit and auditable:

- `PmGridShape{Nx,Ny,Nz}` uses `numerics.treepm_pm_grid_nx/ny/nz`
- `Δx = box_size_x / Nx`, `Δy = box_size_y / Ny`, `Δz = box_size_z / Nz`
- split-spacing scalar used for current TreePM split policy is `Δmesh = cbrt(Δx*Δy*Δz)`
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
- zoom-focused long-range correction (when enabled):
  - authoritative high-resolution membership is computed from the configured zoom center/radius every gravity kick.
  - long-range decomposition is explicit:
    - global coarse PM from all sources,
    - focused PM correction on high-resolution targets only:
      `a_PM_focused(high-res sources) - a_PM_coarse(high-res sources)`,
    - tree short-range residual from all sources.
  - contamination diagnostics count low-resolution source particles/mass inside the configured contamination radius and are emitted as runtime operational events.
- distributed cadence coherence:
  - all ranks advance `gravity_kick_opportunity` together for every gravity kick stage;
  - long-range PM refresh/reuse is checked as a rank-consensus decision, so rank-local counters cannot silently drift;
  - per-kick cadence records (`gravity_kick_opportunity`, `field_version`, refresh flag) are expected to match across ranks.

TreePM source/target ownership in the workflow path is now explicit:

- PM refresh uses all local particle masses on every rank (global mass sources are distributed by slab ownership).
- PM force interpolation and velocity kick updates target only the active particle subset for that kick.
- Distributed tree export/import is only for short-range residual target requests (active targets), not for PM source deposition.

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
- In MPI mode the workflow report records both local owned counts/checksums and reduced global counts/checksums (`particle_count`, `cell_count`, particle-ID sum/xor). This is part of the operational honesty contract for distributed runs: tests must be able to distinguish a true partitioned state from an accidentally replicated one.
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

- `treepm_pm_grid` (legacy report field: currently mirrors `treepm_pm_grid_nx` for compatibility)
- `treepm_update_cadence_steps`
- `treepm_long_range_refresh_count`
- `treepm_long_range_reuse_count`
- `treepm_cadence_records` (per-kick record including step, stage, field version, and refresh/reuse decision)
- `final_state_digest` (deterministic run-result checksum for reproducibility checks under fixed config/ICs)
- operational events include `gravity.zoom_force_diagnostics` with PM/tree/total force norms plus low-resolution contamination counters.
- operational events include `gravity.health_summary` plus targeted `gravity.health_check` warning/fatal events:
  - cheap checks are always-on (PM/force finiteness, sync-state invariants, decomposition/zoom sanity),
  - heavy reference checks are opt-in only when `analysis.diagnostics_execution_policy = all_including_provisional`.
  - fatal gravity-state failures are escalated as runtime errors and are never demoted to warnings.
