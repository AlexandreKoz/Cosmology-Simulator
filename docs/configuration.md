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
  - `isolated_monopole` now activates non-periodic/open PM long-range gravity; current implementation requires `parallel.mpi_ranks_expected = 1` and fails fast otherwise.

## `[cosmology]`

- `omega_matter`, `omega_lambda`, `omega_baryon`
- `hubble_param`, `sigma8`, `scalar_index_ns`
- canonical axis-aware box lengths:
  - `box_size_x`, `box_size_y`, `box_size_z` (recommended; normalized output always emits these keys)
- backward-compatible scalar input:
  - `box_size` / `box_size_mpc_comoving` (legacy alias; when provided alone it maps to `box_size_x=y=z`)

## `[numerics]`

- `time_begin_code`, `time_end_code`
- `max_global_steps`, `hierarchical_max_rung`, `amr_max_level`
- `gravity_softening` / `gravity_softening_kpc_comoving`
- `gravity_solver`, `hydro_solver`

Time/scale semantics (anti-ambiguity contract):

- `time_begin_code` / `time_end_code` are code-time boundaries for integration bookkeeping.
- They are not redshift keys and not SI physical-time keys.
- Runtime cosmological scale-factor authority is `IntegratorState.current_scale_factor` (restart-continuation lane), while redshift remains a derived diagnostic (`z = 1/a - 1` when `a>0`).
- TreePM runtime controls (typed + normalized, no hidden workflow defaults):
  - canonical axis-aware PM grid:
    - `treepm_pm_grid_nx`, `treepm_pm_grid_ny`, `treepm_pm_grid_nz` (normalized output always emits these keys)
  - backward-compatible scalar input:
    - `treepm_pm_grid` (legacy alias; when provided alone it maps to `treepm_pm_grid_nx=ny=nz`)
  - `treepm_asmth_cells` (float, default `1.25`; split scale in mesh-cell units)
  - `treepm_rcut_cells` (float, default `4.5`; user-facing short-range cutoff control in mesh-cell units; this now drives explicit residual-traversal pruning in the TreePM coupling path)
  - `treepm_assignment_scheme` (`cic` or `tsc`; default `cic`)
  - `treepm_enable_window_deconvolution` (bool, default `false`; applies scheme-aware PM transfer deconvolution for both `cic` and `tsc`)
  - `treepm_update_cadence_steps` (int, default `1`; authoritative PM long-range refresh cadence in gravity-kick opportunities)
  - `treepm_pm_decomposition_mode` (`slab` only in Phase 2 freeze; default `slab`)
  - `treepm_tree_exchange_batch_bytes` (uint64 bytes, default `4194304`; cap for tree export/import payload chunking in distributed gravity paths)

TreePM split/cutoff semantics in this phase:

- `Δx = cosmology.box_size_x / numerics.treepm_pm_grid_nx`
- `Δy = cosmology.box_size_y / numerics.treepm_pm_grid_ny`
- `Δz = cosmology.box_size_z / numerics.treepm_pm_grid_nz`
- representative split spacing for TreePM scalar split controls remains
  `Δmesh = cbrt(Δx * Δy * Δz)` in this phase (pencil/isotropic split redesign deferred)
- `r_s = numerics.treepm_asmth_cells * Δmesh`
- `r_cut = numerics.treepm_rcut_cells * Δmesh`

Normalization emits the dimensionless controls exactly as provided and these are the same values consumed by runtime mapping.
`treepm_update_cadence_steps` is consumed directly by runtime:

- refresh PM long-range field every `N` gravity-kick opportunities (`N = treepm_update_cadence_steps`)
- reuse the most recent PM long-range field between refreshes
- record per-kick refresh/reuse metadata in the reference-workflow report and operational diagnostics events
`treepm_rcut_cells` is normalized, carried into provenance, and consumed by runtime residual traversal pruning (node-level AABB pruning + acceptance guard + pair culling beyond `r_cut`).
Zoom long-range strategy semantics:

- `disabled`: standard TreePM split force only.
- `global_coarse_plus_focused_highres_correction`:  
  `a_total = a_PM_global_coarse(all sources) + a_PM_zoom_correction(high-res targets) + a_tree_short_residual(all sources)`  
  where `a_PM_zoom_correction = a_PM_focused(high-res sources) - a_PM_coarse(high-res sources)` to avoid double counting while retaining global tidal content from the coarse PM solve.
`treepm_assignment_scheme` now maps directly to matched PM assignment+gather behavior:

- `cic`: 2-point/axis stencil, lower cost, stronger smoothing.
- `tsc`: 3-point/axis stencil, higher cost, smoother transfer and typically lower anisotropy.

`treepm_enable_window_deconvolution=true` deconvolves the matched deposit/gather transfer window for the selected scheme with a safety floor in k-space (disabled by default).
`treepm_pm_decomposition_mode` and `treepm_tree_exchange_batch_bytes` freeze the distributed TreePM contract surface for Phase 2 infrastructure work; they do not by themselves imply that distributed PM/tree algorithms are already implemented.

## Migration notes (axis-aware TreePM geometry)

- Existing scalar configs remain valid:
  - `cosmology.box_size = L` maps to `box_size_x=box_size_y=box_size_z=L`.
  - `numerics.treepm_pm_grid = N` maps to `treepm_pm_grid_nx=ny=nz=N`.
- Mixed partial axis input is rejected:
  - either provide all three axis keys or provide the scalar compatibility key.
- Normalized config dumps now always emit canonical axis-aware keys (`box_size_{x,y,z}`, `treepm_pm_grid_n{xyz}`).

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
- `snapshot_interval_steps`, `write_restarts`

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

Distributed-memory determinism contract (infrastructure scope):

- `mpi_ranks_expected` must match across all ranks for a valid distributed run contract.
- `deterministic_reduction=true` means rank-local contributions are compared against a rank-ordered deterministic reference sum in tests and smoke paths.
- `deterministic_reduction=false` is allowed only when reduction error is evaluated explicitly against a documented policy mode in code/tests.
- The normalized frozen config hash (`provenance.config_hash_hex`) is computed over the normalized config text itself and must match both the stored normalized-config hash field and the provenance hash field in continuation artifacts.

## `[units]`

- `length_unit`, `mass_unit`, `velocity_unit`, `coordinate_frame`

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
