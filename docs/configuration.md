# Configuration reference

CosmoSim uses a GADGET/AREPO-style `param.txt` workflow and normalizes it into a typed, validated `SimulationConfig` before execution.

Authoritative structures: `include/cosmosim/core/config.hpp` and `src/core/config.cpp`.

## CLI runtime entry point

The shipped executable is now a real config-driven runtime application:

```bash
cosmosim_harness <config.param.txt>
```

That path:

- loads the config via `loadFrozenConfigFromFile(...)`
- validates the typed runtime contract
- constructs the reference workflow
- runs the live scheduler-driven execution path
- writes runtime artifacts into the config-driven run directory

The concrete run directory is:

```text
<output.output_directory>/<output.run_name>
```

## Parsing model

- Input supports `key=value` and `key value` forms.
- Comments: `#`, `;`, `//`.
- Sections (for example `[physics]`) prefix unqualified keys (`cooling_model` becomes `physics.cooling_model`).
- Unknown keys fail by default; compatibility opt-out is explicit (`compatibility.allow_unknown_keys=true`).

## Reproducibility invariants

- Config is normalized to canonical text.
- Canonical text is hashed (stable FNV-1a) for provenance.
- The normalized config snapshot is written to `normalized_config.param.txt` inside the run directory once config loading succeeds.
- Stable naming constraints apply to `output.output_stem`, `output.restart_stem`, and diagnostics/halo stems.
- Mode policy is validated before runtime (`mode.mode`, boundary selection, zoom requirements).

## External strings vs internal typed contract

- User-facing param files remain string-based (`gravity_solver=treepm`, `fb_mode=momentum`, etc.).
- During freeze, policy-like values are converted to enums in `core::SimulationConfig`:
  - `units.coordinate_frame` -> `core::CoordinateFrame`
  - `numerics.gravity_solver` / `numerics.hydro_solver` -> `core::GravitySolver` / `core::HydroSolver`
  - `mode.hydro_boundary` / `mode.gravity_boundary` -> `core::ModeHydroBoundary` / `core::ModeGravityBoundary`
  - `physics.fb_mode` / `physics.fb_variant` -> `core::FeedbackMode` / `core::FeedbackVariant`
  - `physics.uv_background_model` / `physics.self_shielding_model` / `physics.cooling_model` -> `core::UvBackgroundModel` / `core::SelfShieldingModel` / `core::CoolingModel`
- Free-form strings remain strings only where the value is path-like, label-like, or intentionally open-ended (`ic_file`, output stems, table paths, `reionization_model`).

## Central key registry

- `src/core/config.cpp` contains a single authoritative key registry that records canonical keys and defaults.
- The same source file also defines deprecated alias mappings.
- Unknown-key checks and alias translation both use this registry path during freeze.

## Key groups

## `schema_version`

- `schema_version` (int, required to be `1` in this build)

## `[mode]`

- `mode` (`cosmo_cube`, `zoom_in`, `isolated_galaxy`, `isolated_cluster`)
- `ic_file`
- `zoom_high_res_region` (bool)
- `zoom_region_file` (required when `zoom_high_res_region=true`)
- `zoom_long_range_strategy` (`disabled`, `global_coarse_plus_focused_highres_correction`)
- `zoom_region_center_x`, `zoom_region_center_y`, `zoom_region_center_z` (comoving Mpc)
- `zoom_region_radius` (comoving Mpc; required `>0` when `zoom_high_res_region=true`)
- `zoom_focused_pm_grid_nx`, `zoom_focused_pm_grid_ny`, `zoom_focused_pm_grid_nz` (required when focused correction is enabled)
- `zoom_contamination_radius` (comoving Mpc; defaults to `zoom_region_radius` when `<=0`)
- `hydro_boundary` (`auto`, `periodic`, `open`, `reflective`)
- `gravity_boundary` (`auto`, `periodic`, `isolated_monopole`)
  - `isolated_monopole` activates non-periodic/open PM long-range gravity.
    Single-rank execution is the normal path. Multi-rank execution is accepted
    only by the bounded root gather/solve/scatter compatibility route governed
    by `parallel.isolated_pm_root_workspace_limit_bytes`; it is not a scalable
    distributed-open-PM claim.

## `[cosmology]`

- `omega_matter`, `omega_lambda`, `omega_baryon`
- `hubble_param`, `sigma8`, `scalar_index_ns`
- canonical axis-aware box lengths:
  - `box_size_x`, `box_size_y`, `box_size_z` (recommended; normalized output always emits these keys)
- backward-compatible scalar input:
  - `box_size` / `box_size_mpc_comoving` (legacy alias; when provided alone it maps to `box_size_x=y=z`)

## `[numerics]`

- `t_code_begin`, `t_code_end`
- `a_begin`, `a_end`, `z_begin`, `z_end`, `t_phys_begin`, `t_phys_end`, `integrator_time_variable`
- `max_global_steps`, `hierarchical_max_rung`, `amr_max_level`
  - `hierarchical_max_rung` defaults to and currently requires `0` for the
    production `ReferenceWorkflow`. Nonzero values fail validation because the
    production KDK state does not yet carry per-element kick/drift epochs needed
    to advance mixed rungs without skipping elapsed time. Standalone scheduler
    APIs/tests remain available; they are not production multirate evidence.
- `gravity_softening` / `gravity_softening_kpc_comoving`
  - the production gravity timestep combines this comoving length with
    comoving-coordinate acceleration `|A|/a^3`, not directly with scale-free
    TreePM `|A|` or the peculiar-velocity response `|A|/a^2`; after step commit
    it evaluates `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`.
- `gravity_solver`, `hydro_solver`

Time/scale semantics (anti-ambiguity contract):

- `t_code_begin` / `t_code_end` are code-time boundaries for integration bookkeeping.
- Legacy `time_begin_code`/`time_end_code` and `initial_scale_factor`/`initial_redshift` are accepted only as user-input aliases.
- They are not redshift keys and not SI physical-time keys.
- Committed/restart cosmological scale-factor authority is
  `IntegratorState.current_scale_factor`, while redshift remains derived
  (`z=1/a-1`). Within an in-flight step, `StepContext.timeline_step` owns stage
  epochs; post-drift hydro uses its `scale_factor_end/hubble_end_code` before
  `IntegratorState` commits.
- TreePM runtime controls (typed + normalized, no hidden workflow defaults):
  - canonical axis-aware PM grid:
    - `treepm_pm_grid_nx`, `treepm_pm_grid_ny`, `treepm_pm_grid_nz` (normalized output always emits these keys)
  - backward-compatible scalar input:
    - `treepm_pm_grid` (legacy alias; when provided alone it maps to `treepm_pm_grid_nx=ny=nz`)
  - `treepm_asmth_cells` (float, default `1.25`; split scale in mesh-cell units)
  - `treepm_rcut_cells` (finite positive float, default `6.25`; user-facing short-range cutoff control in mesh-cell units; this drives explicit residual-traversal pruning in the TreePM coupling path)
  - `treepm_tree_opening_criterion` (`geometric`, `com_distance`, or `relative_force_error`; default `com_distance`)
  - `treepm_tree_opening_theta` (finite positive float, default `0.7`; Barnes-Hut threshold and deterministic first-evaluation fallback threshold for the relative-force criterion)
  - `treepm_tree_relative_force_tolerance` (finite positive float, default `0.005`; relative-force MAC tolerance)
  - `treepm_tree_relative_force_acceleration_floor` (finite positive float in code-acceleration units, default `1e-30`; lower bound for the previous-acceleration scale used by the relative-force MAC)
  - `treepm_assignment_scheme` (`cic` or `tsc`; default `tsc`; the independent Ewald gate certifies TSC, while CIC remains available as a diagnostic/compatibility variant)
  - `treepm_enable_window_deconvolution` (bool, default `true`; applies
    scheme-aware PM transfer deconvolution for both `cic` and `tsc` on the
    periodic solver path). Isolated/open gravity requires this value to be
    explicitly `false`; typed config validation rejects `true` because that
    solver path has no periodic assignment window to deconvolve.
  - `treepm_update_cadence_steps` (int, default and currently required value
    `1`; every integrator-issued production force-refresh surface rebuilds the
    PM long-range field)
  - `treepm_pm_decomposition_mode` (`slab` or `pencil`; default `slab`; `pencil` selects FFTW-MPI transposed spectral ownership while real-space particle deposition/interpolation retains x-slab ownership)
  - `treepm_tree_exchange_batch_bytes` (uint64 bytes, default `4194304`; cap for tree export/import payload chunking in distributed gravity paths)

TreePM split/cutoff semantics in this phase:

- `Δx = cosmology.box_size_x / numerics.treepm_pm_grid_nx`
- `Δy = cosmology.box_size_y / numerics.treepm_pm_grid_ny`
- `Δz = cosmology.box_size_z / numerics.treepm_pm_grid_nz`
- representative split spacing for TreePM scalar split controls remains
  `Δmesh = cbrt(Δx * Δy * Δz)` in this phase (pencil/isotropic split redesign deferred)
- `r_s = numerics.treepm_asmth_cells * Δmesh`
- `r_cut = numerics.treepm_rcut_cells * Δmesh`

For periodic gravity, typed validation additionally requires
`r_cut < 0.5 * min(Lx,Ly,Lz)`. The short-range residual evaluates one periodic
minimum image per source, so equality or a larger cutoff would make that
single-image traversal ambiguous. Invalid decks fail with the derived cutoff,
half-shortest-axis bound, and PM grid shape. Increase PM grid resolution or
reduce `treepm_rcut_cells`; the restriction does not apply to the isolated/open
operator.

Normalization emits the dimensionless controls exactly as provided and these are the same values consumed by runtime mapping.

Tree opening-criterion semantics are:

- `geometric`: accept an internal node when `l / r < theta`, where `l` is node width and `r` is target-to-node-COM distance.
- `com_distance` (default): tighten the Barnes-Hut test to account for the node center-to-COM offset, `(l + delta_com) / r < theta`.
- `relative_force_error`: accept when `G M l^2 <= alpha max(|a_previous|, a_floor) r^4`, using `treepm_tree_relative_force_tolerance` for `alpha`. When a finite previous-acceleration magnitude is unavailable, traversal deterministically falls back to the `com_distance` rule with `treepm_tree_opening_theta`.

All three criteria, the opening threshold, relative tolerance, and acceleration floor are typed fields in the frozen config and are emitted in normalized config text. This makes criterion selection and its numerical thresholds part of the reproducibility hash.

`treepm_update_cadence_steps` is production-restricted to `1`:

- Values other than `1` fail config validation with an explicit production-safety error.
- Together with required `hierarchical_max_rung=0`, every
  integrator-issued rank-coordinated production force-refresh surface rebuilds
  PM.
- The coordinator's explicit lower-level reuse API requires a transient cache
  signature match for force epoch, force-evaluation scale factor, `G_code`,
  split scale, all three box axes, assignment scheme, boundary condition, PM
  decomposition mode, and window-deconvolution policy. A missing or
  incompatible cache makes reuse fail coherently; it is not silently converted
  into a solve. Divergent explicit refresh votes fail before PM collectives.
- Particle decomposition epoch is deliberately not part of PM compatibility:
  particles may change owner while the PM field remains owned by the same fixed
  FFT slabs. Tree packets still carry the actual workflow decomposition epoch,
  and dense-row force history is invalidated on an ownership change.
- That signature prevents accidental reuse across known solver/runtime changes,
  but it is not a long-range time predictor/interpolator for arbitrary
  multi-kick evolution.

`treepm_rcut_cells` is normalized, carried into provenance, and consumed by runtime residual traversal pruning (node-level AABB pruning + acceptance guard + pair culling beyond `r_cut`).

The certified default combines `treepm_asmth_cells=1.25` with
`treepm_rcut_cells=6.25`, so `r_cut/r_s=5`. The former `4.5`-cell cutoff and
other positive explicit values remain available for compatibility and
parameter sweeps when the periodic half-shortest-axis geometry guard is also
satisfied, but are not covered by the default Ewald certification: the
4.5-cell profile showed about nine-percent force error at the cutoff transition. Raising
the cutoff increases short-range traversal radius, pair work, and distributed
target/response volume; this is a deliberate accuracy/cost tradeoff and is part
of normalized config/provenance.
Zoom long-range strategy semantics:

- `disabled`: standard TreePM split force only.
- `global_coarse_plus_focused_highres_correction`:  
  `a_total = a_PM_global_coarse(all sources) + a_PM_zoom_correction(high-res targets) + a_tree_short_residual(all sources)`  
  where `a_PM_zoom_correction = a_PM_focused(high-res sources) - a_PM_coarse(high-res sources)` to avoid double counting while retaining global tidal content from the coarse PM solve.
  The focused correction uses the isolated/open convolution and forces window
  deconvolution off; isolated deconvolution is unsupported even when the global
  periodic PM profile enables TSC deconvolution.
`treepm_assignment_scheme` now maps directly to matched PM assignment+gather behavior:

- `cic`: 2-point/axis stencil, lower cost, stronger smoothing.
- `tsc`: 3-point/axis stencil, higher cost, smoother transfer and typically lower anisotropy.

`treepm_enable_window_deconvolution=true` deconvolves the matched deposit/gather
transfer window for the selected scheme with a safety floor in k-space and is
enabled by default for the certified periodic TSC profile. It is invalid with
an isolated/open gravity boundary. There is no silent mode-dependent override:
tracked isolated decks set the key to `false`, and user isolated decks must do
the same.
`treepm_pm_decomposition_mode` and `treepm_tree_exchange_batch_bytes` are live
distributed runtime controls, not ignored Phase 2 placeholders. Both slab and
pencil FFTW-MPI paths are implemented; their current np1--np4 correctness
evidence remains a small-cluster validation boundary, not a large-scale
performance claim.

## Migration notes (axis-aware TreePM geometry)

- Existing scalar configs remain valid:
  - `cosmology.box_size = L` maps to `box_size_x=box_size_y=box_size_z=L`.
  - `numerics.treepm_pm_grid = N` maps to `treepm_pm_grid_nx=ny=nz=N`.
- Mixed partial axis input is rejected:
  - either provide all three axis keys or provide the scalar compatibility key.
- Normalized config dumps now always emit canonical axis-aware keys (`box_size_{x,y,z}`, `treepm_pm_grid_n{xyz}`).

## Migration note (isolated/open window deconvolution)

The global default for `numerics.treepm_enable_window_deconvolution` remains
`true` for the periodic TSC profile. Consequently, an isolated/open deck that
omits this key inherits `true` and now fails during typed config loading with a
key-specific error. Set
`numerics.treepm_enable_window_deconvolution = false` explicitly. Runtime code
does not silently rewrite the frozen config, so normalized dumps and
reproducibility hashes continue to describe the solver policy actually used.

## Migration note (periodic cutoff geometry)

Periodic decks must satisfy
`treepm_rcut_cells * cbrt((Lx/Nx)(Ly/Ny)(Lz/Nz)) < min(Lx,Ly,Lz)/2`.
Older coarse-mesh decks that were accepted at or above this boundary now fail
during typed loading, and direct API callers are rejected again at TreePM solve
entry. Increase `treepm_pm_grid_n{x,y,z}` as appropriate or explicitly lower
`treepm_rcut_cells`; either choice changes the normalized numerical profile and
therefore requires fresh force-accuracy evidence.

## `[physics]`

Core toggles:

- `enable_cooling`, `enable_star_formation`, `enable_feedback`, `enable_stellar_evolution`

Cooling/heating:

- `reionization_model`, `uv_background_model`, `self_shielding_model`
- `cooling_model`, `metal_line_table_path`, `temperature_floor_k`

Star formation:

- `sf_density_threshold_code`, `sf_temperature_threshold_k`
- `sf_min_converging_flow_rate_code`, `sf_epsilon_ff`
- `sf_min_star_particle_mass_code`, `sf_stochastic_spawning`, `sf_random_seed`

Stellar feedback:

- `fb_mode` (`thermal`, `kinetic`, `momentum`, `thermal_kinetic_momentum`)
- `fb_variant` (`none`, `delayed_cooling`, `stochastic`)
- `fb_use_returned_mass_budget`
- `fb_epsilon_thermal`, `fb_epsilon_kinetic`, `fb_epsilon_momentum`
- `fb_sn_energy_erg_per_mass_code`, `fb_momentum_code_per_mass_code`
- `fb_neighbor_count`, `fb_delayed_cooling_time_code`
- `fb_stochastic_event_probability`, `fb_random_seed`

Stellar evolution + AGN:

- `stellar_evolution_table_path`, `stellar_evolution_hubble_time_years`
- `enable_black_hole_agn`
- `bh_seed_halo_mass_threshold_code`, `bh_seed_mass_code`, `bh_seed_max_per_cell`
- `bh_alpha_bondi`, `bh_use_eddington_cap`, `bh_epsilon_r`, `bh_epsilon_f`
- `bh_feedback_coupling_efficiency`, `bh_duty_cycle_active_edd_ratio_threshold`
- `bh_proton_mass_si`, `bh_thomson_cross_section_si`, `bh_newton_g_si`, `bh_speed_of_light_si`

Tracers:

- `enable_tracers`, `tracer_track_mass`, `tracer_min_host_mass_code`

## `[output]`

- `run_name`
- `output_directory` (output root; the runtime appends `run_name` to form the concrete run directory)
- `output_stem`, `restart_stem` (stable character set only)
- `snapshot_interval_steps` (zero disables step-modulo events)
- `snapshot_interval_time_code` (zero disables code-time events; positive values are
  anchored at `numerics.time_begin_code` and steps are clipped rather than crossing them)
- `write_restarts`

At least one snapshot interval must be positive. Code-time cadence is part of the
normalized config and its next ordered event is restart-authoritative in restart
schema v21.

## `[analysis]`

- `enable_diagnostics`, `enable_halo_workflow`, `halo_on_the_fly`
- `diagnostics_execution_policy` (`run_health_only`, `run_health_and_light_science`, `all_including_provisional`)
- `run_health_interval_steps`, `science_light_interval_steps`, `science_heavy_interval_steps`
- `retention_bundle_count`
- `power_spectrum_mesh_n`, `power_spectrum_bin_count`, `sf_history_bin_count`, `quicklook_grid_n`
- `diagnostics_stem`, `halo_catalog_stem`, `merger_tree_stem`
- `halo_fof_linking_length_factor`, `halo_fof_min_group_size`
- `halo_include_gas`, `halo_include_stars`, `halo_include_black_holes`

Diagnostics maturity policy:

- `run_health_only`: only infrastructure health counters run.
- `run_health_and_light_science` (default): run-health + validated lightweight science diagnostics.
- `all_including_provisional`: also enables provisional/reference heavy diagnostics (currently the direct-summation power spectrum).

## `[parallel]`

- `mpi_ranks_expected`, `omp_threads`, `gpu_devices`

GPU execution contract:

- `gpu_devices = 0` keeps PM execution on the host.
- `gpu_devices > 0` explicitly requests the CUDA PM path for CIC assignment/interpolation. The runtime validates that CUDA devices are visible and assigns ranks to devices round-robin across the requested device pool.
- `gpu_devices` must be `>= 0`; requesting more devices than are visible is a configuration/runtime error rather than a silent fallback.
- `deterministic_reduction`
- `decomposition_runtime_rebalance_enabled`
- `decomposition_debug_exact_ownership_audit` (default `false`): keep the routine runtime rebalance path on distributed compact SFC-cut metadata, but run the exact global ownership partition audit after a migration commit in debug/small MPI runs. This exact audit is read-only evidence and is not a production ownership authority.
- `isolated_pm_root_workspace_limit_bytes` (default `268435456`): maximum rank-0 transient workspace allowed for the current isolated/open PM root-gather path. Multi-rank isolated PM is accepted only inside this explicit small-grid budget and diagnostics must not describe it as scalable.
- `zoom_high_res_allgather_limit_bytes` (default `268435456`): maximum transient four-field high-resolution source payload allowed for the current focused zoom PM correction all-gather. Exceeding it fails before allocation/communication.

Distributed-memory determinism contract (infrastructure scope):

- `mpi_ranks_expected` must match across all ranks for a valid distributed run contract.
- `deterministic_reduction=true` means rank-local contributions are compared against a rank-ordered deterministic reference sum in tests and smoke paths.
- `deterministic_reduction=false` is allowed only when reduction error is evaluated explicitly against a documented policy mode in code/tests.
- The normalized frozen config hash (`provenance.config_hash_hex`) is computed over the normalized config text itself and must match both the stored normalized-config hash field and the provenance hash field in continuation artifacts.

## `[units]`

- `length_unit`, `mass_unit`, `velocity_unit`, `coordinate_frame`

The three unit strings define one `core::UnitSystem` and its derived code-time
unit `T_code=L_code/V_code`. Production gravity derives

```text
G_code = G_SI * M_SI_per_code * T_SI_per_code^2 / L_SI_per_code^3
```

through `core::newtonGravitationalConstantCode(...)`; it does not assume
`G_code=1`. PM and the TreePM residual receive the same value. The normalized
unit strings and config hash remain reproducibility authority; the derived
value is also reported in the `gravity.treepm_setup` operational event. This is
an additive public helper, not a new config key or snapshot/restart schema
field.

## `[compatibility]`

- `allow_unknown_keys`

## Example configs

Canonical examples are in `configs/`:

- `minimal_cosmosim.param.txt`
- `cosmo_cube.param.txt`
- `zoom_in.param.txt`
- `isolated_galaxy.param.txt`
- `isolated_cluster.param.txt`
- `cooling_relaxation.param.txt`
- `configs/release/release_smoke_*.param.txt`
