# Configuration reference

CosmoSim uses a GADGET/AREPO-style `param.txt` workflow and normalizes it into a typed, validated `SimulationConfig` before execution.

Authoritative structures: `include/cosmosim/core/config.hpp` and `src/core/config.cpp`.

## Parsing model

- Input supports `key=value` and `key value` forms.
- Comments: `#`, `;`, `//`.
- Sections (for example `[physics]`) prefix unqualified keys (`cooling_model` becomes `physics.cooling_model`).
- Unknown keys fail by default; compatibility opt-out is explicit (`compatibility.allow_unknown_keys=true`).

## Reproducibility invariants

- Config is normalized to canonical text.
- Canonical text is hashed (stable FNV-1a) for provenance.
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
- `hydro_boundary` (`auto`, `periodic`, `open`, `reflective`)
- `gravity_boundary` (`auto`, `periodic`, `isolated_monopole`)

## `[cosmology]`

- `omega_matter`, `omega_lambda`, `omega_baryon`
- `hubble_param`, `sigma8`, `scalar_index_ns`
- `box_size` / `box_size_mpc_comoving`

## `[numerics]`

- `time_begin_code`, `time_end_code`
- `max_global_steps`, `hierarchical_max_rung`, `amr_max_level`
- `gravity_softening` / `gravity_softening_kpc_comoving`
- `gravity_solver`, `hydro_solver`

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

- `run_name`, `output_directory`
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
- `deterministic_reduction`

Distributed-memory determinism contract (infrastructure scope):

- `mpi_ranks_expected` must match across all ranks for a valid distributed run contract.
- `deterministic_reduction=true` means rank-local contributions are compared against a rank-ordered deterministic reference sum in tests and smoke paths.
- `deterministic_reduction=false` is allowed only when reduction error is evaluated explicitly against a documented policy mode in code/tests (for example via `parallel::ReductionAgreementMode::{kAbsoluteOnly,kRelativeOnly,kAbsoluteAndRelative,kAbsoluteOrRelative}` used by `parallel::satisfiesReductionAgreement`; no implicit “close enough” claims).
- Current pseudo-multi-rank infrastructure smoke/integration checks use `kAbsoluteOrRelative` by default and document that choice at call sites.
- The normalized frozen config hash (`provenance.config_hash_hex`) is treated as a cross-rank freeze identity and must not diverge across ranks in the same run.
- Cross-rank config consensus diagnostics should emit property-level mismatch records (property name + baseline/offending rank/value pairs), not only pass/fail booleans.

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
