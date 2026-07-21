# Reference workflow runtime path

This document describes the live config-driven runtime path assembled by `workflows::ReferenceWorkflowRunner` and exercised by the shipped executable target `cosmosim_harness`.

## What the runtime now does

A valid `param.txt` run goes through the following real path:

1. `cosmosim_harness <config.param.txt>` loads the typed frozen config with `loadFrozenConfigFromFile(...)`.
2. The runner validates the mode/build contract before stepping.
3. Initial conditions are dispatched from the authoritative typed convention:
   - `mode.ic_convention=generated` requires `mode.ic_file=generated`.
   - `chui_canonical_v1` reads the canonical CHUÍ v1 contract.
   - `gadget_arepo_bridge_v1` reads the explicitly versioned external bridge contract.
   - `manifest_v1` loads the strict audit manifest named by `mode.ic_manifest_file`.
   External paths without an explicit convention fail before runtime. Relative IC
   paths are resolved from the config; relative manifest source members are resolved
   from the manifest directory.
4. Serial input is chunked and multifile-aware. With MPI, source chunks are assigned
   once globally, converted once, and routed directly to deterministic initial owners
   with bounded staging. The result is marked already partitioned, so the legacy
   replicated-state ownership compactor is not applied.
5. The live run loop uses `HierarchicalTimeBinScheduler` as the execution driver.
5. The orchestrator executes the canonical KDK stage sequence through the
   gravity owner and the current hydro callback.
6. Outputs, restart checkpoints, diagnostics, normalized config, and operational reports are written into the config-driven run directory.

## Runtime contract

The run directory is:

- `<output.output_directory>/<output.run_name>` in normal CLI usage.
- `<override_root>/<output.run_name>` only for tests/benchmarks that intentionally override the output root.

The runner honors these config fields directly:

- `mode.mode`
- `mode.ic_file`
- `mode.ic_convention`
- `mode.ic_manifest_file`
- `mode.ic_chunk_particle_count`
- `mode.ic_staging_particle_count`
- `mode.ic_part_type2_policy` / `mode.ic_part_type3_policy`
- `output.output_directory`
- `output.run_name`
- `output.output_stem`
- `output.restart_stem`
- `output.snapshot_interval_steps`
- `output.snapshot_interval_time_code`
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

### Initial-condition ownership boundary

`InitialConditionRuntime` owns convention dispatch and borrows MPI/profiler services
from the process `RuntimeServices` bundle. The `io` layer owns file-set inspection,
conversion, wire records, and collective validation; it does not initialize MPI or
create a second runtime context. Distributed import finalizes rank-local state only
and returns `already_partitioned=true`. Generated or caller-supplied replicated state
continues through the established ownership initialization path.

Before stepping, distributed import verifies a canonical manifest digest, exact global
ID uniqueness, source/chunk coverage, ownership completeness and exclusivity, species
counts and mass totals, finite/domain-valid fields, and sidecar invariants. Rank-local
failures are coordinated before later collective phases so malformed input fails rather
than stranding peers in MPI.

The validated audit manifest is written by rank zero as `ic_manifest.json` in the run
directory. Every rank emits `io.ic_ingestion.summary` counters through the shared
profiler.

TreePM runtime mapping is explicit and auditable:

Output dispatch supports both step-modulo cadence and code-time cadence. A positive
`snapshot_interval_time_code` creates ordered events anchored at
`numerics.time_begin_code`; a larger proposed KDK interval is clipped at the event,
and the next event plus the pre-clip next-step proposal are persisted in a v21
checkpoint. Restart restores that authority instead of recomputing cadence from the
resumed step number. Schema v20 reads explicitly materialize the time-event lane as
disabled.

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

`workflows::TimeCoordinator` uses the dependency-safe dispatcher entry point
of `core::StepOrchestrator` and preserves the canonical order:

1. `gravity_kick_pre`
2. `drift`
3. `force_refresh`
4. `hydro_update`
5. `source_terms`
6. `gravity_kick_post`
7. `analysis_hooks`
8. `output_check`

## Gravity, hydro, and scheduler ownership

- `workflows::GravityRuntime` owns the live `gravity::TreePmCoordinator`, PM
  cadence decisions, force-cache invalidation/import/export, and gravity stage
  contribution. Hydro consumes only its read-only acceleration provider; the
  output/restart owner consumes only its persistence provider.
- `workflows::HydroAmrRuntime` owns the live Godunov finite-volume and supported
  AMR paths using the hydro core solver, MUSCL-Hancock reconstruction, HLLC
  fluxes, blocking ghost freshness, patch exchange, and existing reflux rules.
- `workflows::SourceRuntime` preserves star-formation then black-hole model
  ordering, while `workflows::AnalysisRuntime` owns stage-audit and diagnostics
  contributions.
- `workflows::RungZeroTimeState` is the sole workflow owner of particle and
  gas-cell schedulers, integrator truth, and output cadence. Active sets come
  from `HierarchicalTimeBinScheduler::beginSubstep()` and completed substeps
  are committed through `endSubstep()` before restart output is written. The
  production workflow currently requires rung zero; nonzero rungs are rejected
  because per-element kick/drift epochs are not yet authoritative.
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

External-IC runs additionally write `ic_manifest.json`, the exact validated
unit/frame/species/header contract used by ingestion. Generated-IC runs omit
this external-import artifact.

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

The additive `final_time_code` and `final_scale_factor` report fields expose
the committed integration endpoint without changing restart or snapshot
schemas. Before KDK preparation, the workflow clips an interval that would
cross `numerics.time_end_code`; the run therefore commits the configured final
time exactly and records a `time.endpoint_clip` operational event containing
the original and accepted intervals. Existing callers that do not read these
new fields require no migration.

The production substep loop consumes scheduler-owned compact active arrays as
read-only spans. It does not build duplicate particle/cell index vectors. One
`TransientStepWorkspace` is owned for the run segment and cleared between
executed KDK steps; clearing resets sizes and the monotonic scratch offset while
preserving vector and arena capacity. Profiling exposes
`workflow_workspace_reuses` and `scheduler_active_index_copy_bytes` (zero on
this path) so allocation/copy regressions remain visible.

## Runtime composition and resource views

`ReferenceWorkflowRunner` constructs frozen config, one `RuntimeServices`
bundle, IC/migration ownership, `RungZeroTimeState`, and the descriptor-backed
runtime composition. It then makes one high-level call to
`TimeCoordinator::runRungZeroSegment(...)`; numerical stages and the detailed
substep loop are not implemented in the runner.

Built-in `RuntimeModuleDescriptor` factories construct the real stage owners.
`RuntimeExecutionPlan` dispatches only typed `DriftParticleStageView`,
`GravityStageView`, `HydroAmrStageView`, `SourceMutationStageView`,
`AnalysisStageView`, and `OutputRestartStageView` capabilities. Each carries a
lease over the relevant particle/cell/gas-identity generation, scheduler tick,
and/or step epoch. Leases are transient and never restart truth. Reorder,
migration, compaction, scheduler-tick advance, or step advance invalidates a
captured view when its guarded epoch changes; callers must rebuild it.

`ReferenceWorkflowOptions::register_runtime_modules` is an additive
test/embedding seam. It runs after built-ins are declared and before the
registry freezes, so a caller can add an independent typed descriptor without
editing `ReferenceWorkflowRunner` or a central callback list. It is not a
configuration key, is not serialized, and does not affect restart schema or
provenance authority when unused.
