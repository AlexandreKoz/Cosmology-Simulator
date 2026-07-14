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
- `numerics.treepm_update_cadence_steps` (authoritative PM long-range refresh
  cadence; production value is exactly `1`)
- `numerics.hierarchical_max_rung` (production value is exactly `0`)

TreePM runtime mapping is explicit and auditable:

- `PmGridShape{Nx,Ny,Nz}` uses `numerics.treepm_pm_grid_nx/ny/nz`
- `Δx = box_size_x / Nx`, `Δy = box_size_y / Ny`, `Δz = box_size_z / Nz`
- split-spacing scalar used for current TreePM split policy is `Δmesh = cbrt(Δx*Δy*Δz)`
- `r_s = treepm_asmth_cells * Δmesh`
- `r_cut = treepm_rcut_cells * Δmesh`
- periodic typed validation and solve entry require
  `r_cut < 0.5*min(box_size_x,box_size_y,box_size_z)` because the residual tree
  evaluates one minimum image per source
- Assignment/deconvolution are wired from typed config, not hidden workflow constants
- `treepm_assignment_scheme` maps directly to runtime PM assignment+gather (`cic` or `tsc`)
- For an isolated/open gravity boundary, typed config validation requires
  `treepm_enable_window_deconvolution=false` before runtime construction. The
  periodic default is not silently overridden, so isolated decks must state
  the disabled policy explicitly.
- Production validation requires `treepm_update_cadence_steps = 1`, so every
  authorized rank-coordinated force-refresh rebuilds the PM field. Lower-level
  explicit reuse remains available for tests/future integration and fails if
  any rank lacks a compatible cache; it never silently changes into refresh.
- Field provenance includes force epoch and build step/scale-factor metadata.
  Particle decomposition epoch is deliberately excluded from PM cache
  compatibility because the field is owned by fixed FFT slabs.
- PM and tree both consume `core::newtonGravitationalConstantCode(...)` from
  frozen units and return a scale-free comoving kernel `A`. Collisionless KDK
  and the gas conservative source apply Hubble drag/expansion plus the `A/a²`
  response to their respective state; both use integrator-owned `a` and
  `H_code`.
- Hydro runs post-drift before step commit, so fixed-grid and AMR source
  contexts take `a,H_code` from
  `timeline_step.scale_factor_end/hubble_end_code`, not the still-step-begin
  `IntegratorState` values.
- Adaptive gravity timestep proposals run after commit. With comoving
  softening they call `computeComovingGravityTimeStep` with `eps_com`,
  scale-free `|A|`, and committed `a`. The public helper validates `a`, converts
  to comoving-coordinate acceleration `|A|/a^3`, and uses
  `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`.
- `r_cut` is resolved from typed config and drives explicit residual-traversal pruning in the TreePM residual path
- zoom-focused long-range correction (when enabled):
  - authoritative high-resolution membership is loaded from `mode.zoom_region_file` when provided (plain-text IDs or HDF5 `ParticleIDs`/`Region/ParticleIDs`); otherwise the configured zoom center/radius fallback is used.
  - long-range decomposition is explicit:
    - global coarse PM from all sources,
    - focused PM correction on high-resolution targets only:
      `a_PM_focused(high-res sources) - a_PM_coarse(high-res sources)`,
    - the focused correction uses isolated/open convolution and forces PM
      window deconvolution off; isolated deconvolution is not certified,
    - tree short-range residual from all sources.
  - contamination diagnostics count low-resolution source particles/mass inside the configured contamination radius and are emitted as runtime operational events.
- distributed cadence coherence:
  - all ranks advance `gravity_kick_opportunity` together for every
    integrator-issued PM synchronization event;
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
3. `force_refresh`
4. `hydro_update`
5. `source_terms`
6. `gravity_kick_post`
7. `analysis_hooks`
8. `output_check`

## Gravity, hydro, and scheduler ownership

- Gravity is executed through the live `gravity::TreePmCoordinator` callback.
- Hydro is executed through the live Godunov finite-volume callback using the hydro core solver, MUSCL-Hancock reconstruction, and HLLC fluxes.
- Active sets come from `HierarchicalTimeBinScheduler::beginSubstep()` and completed substeps are committed through `endSubstep()` before restart output is written. The production workflow currently requires rung zero; nonzero rungs are rejected because per-element kick/drift epochs are not yet authoritative.
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
- `treepm_pm_grid_nx`, `treepm_pm_grid_ny`, `treepm_pm_grid_nz`, `treepm_pm_grid_shape`
- `treepm_update_cadence_steps`
- `treepm_long_range_refresh_count`
- `treepm_long_range_reuse_count`
- `treepm_cadence_records` (per-kick record including step, stage, field version, and refresh/reuse decision)
- `final_state_digest` (deterministic run-result checksum for reproducibility checks under fixed config/ICs)
- operational events include `gravity.zoom_force_diagnostics` with PM/tree/total force norms plus low-resolution contamination counters.
- `gravity.pm_long_range_field` records the current pending refresh directive's
  opportunity, field version, build step, and build scale factor before cadence
  commit. It must not report the previous committed `pm_sync_state` values.
- Gravity diagnostic doubles in operational-event payloads use scientific
  `max_digits10` formatting. Small `G_code`, relative-MAC floors, derived
  split/cutoff/softening scales, force norms, and scale factors therefore
  remain distinguishable from zero.
- operational events include `gravity.health_summary` plus targeted `gravity.health_check` warning/fatal events:
  - cheap checks are always-on (PM/force finiteness, sync-state invariants, decomposition/zoom sanity),
  - heavy reference checks are opt-in only when `analysis.diagnostics_execution_policy = all_including_provisional`.
  - fatal gravity-state failures are escalated as runtime errors and are never demoted to warnings.
