# Output schema guide

This guide documents snapshot, restart, and provenance payload contracts currently implemented in CosmoSim.

Authoritative interfaces:

- `include/cosmosim/io/snapshot_hdf5.hpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `include/cosmosim/core/provenance.hpp`
- `include/cosmosim/core/profiling.hpp` (operational run-event report)

## 1) Snapshot schema (GADGET/AREPO interoperability)

Current schema identity:

- `schema_name = gadget_arepo_v4`
- `schema_version = 4`

Logical groups:

- `/Header`
- `/Config`
- `/Provenance`
- `/PartType0` ... `/PartType5`

Canonical fields and accepted read aliases:

- `Coordinates` (`Coordinates`, `Position`, `POS`)
- `Velocities` (`Velocities`, `Velocity`, `VEL`)
- `Masses` (`Masses`, `Mass`)
- `ParticleIDs` (`ParticleIDs`, `ParticleID`, `ID`)

`SnapshotIoPolicy` controls compression, chunk size, and whether alias groups are emitted.

## 2) Restart schema

Current restart identity:

- `name = cosmosim_restart_v20`

## H1 workflow force-cache restart note

Restart schema v20 adds `/gravity_force_cache`, containing the particle and gas-cell acceleration triplets consumed by the next KDK pre-kick, plus the persisted `IntegratorState::pm_refresh_enabled` policy bit. The ReferenceWorkflow writes a valid cache only at a restart-safe completed boundary and verifies it on read. Direct low-level checkpoint callers may write an explicitly invalid/empty cache, but those artifacts do not constitute exact `ReferenceWorkflow` continuation proof.
- `version = 19`

Restart payload includes:

- `SimulationState`
- hydro geometry state under `/state/cells`, `/state/gas_cells`, and `/state/patches`
  (`/state/gas_cells/velocity_[xyz]_peculiar` is persistent cell-local hydro velocity truth)
  (`patch_id`, `level`, `first_cell`, `cell_count`, `owning_rank`)
- authoritative gas-cell identity records under `/state/gas_cell_identity`:
  `gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`,
  `local_cell_row`, `@local_row_reconstruction_policy`, and
  `@identity_generation_at_write`
- `IntegratorState`
- `HierarchicalTimeBinScheduler` persistent state (`TimeBinPersistentState`)
- distributed TreePM continuation metadata (`DistributedRestartState`)
- normalized config text and normalized config hash
- provenance record
- payload integrity hash
- `/restart_diagnostics` audit metadata (schema/boundary/scheduler/PM/output/stochastic summaries; non-authoritative)

Compatibility rule is explicit through `isRestartSchemaCompatible(version)`.

## 2.1) Field ownership table (snapshot vs restart)

| Ownership | Fields |
|---|---|
| Shared metadata contract | normalized config text/hash, provenance payload, schema identity |
| Snapshot-only (interoperable science output) | `/Header` cosmology attrs, `/PartTypeN` particle datasets, read aliases (`Position`, `VEL`, `ID`, etc.) |
| Restart-only (exact continuation state) | compact `/restart_diagnostics` audit metadata plus full `SimulationState` hot/cold lanes, `StateMetadata`, hydro patch geometry lanes (`CellSoa::patch_index`, gas-cell mirror lanes, authoritative `/state/gas_cell_identity` records, gas-cell velocity/thermodynamic lanes, `PatchSoa` descriptors and ownership), module sidecars + schema versions, `IntegratorState`, scheduler persistent state (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`), distributed TreePM restart state (`decomposition_epoch`, owning-rank table, PM slab layout, cadence/long-range metadata, restart policy), payload integrity hashes. `ParticleSoa::time_bin` and `CellSoa::time_bin` are retained only as derived mirrors/diagnostics; exact continuation imports scheduler state, rejects stale mirror conflicts, and rebuilds mirrors from scheduler authority. For gas cells with a local parent, the cell mirror is validated/rebuilt against the parent gas particle's scheduler entry; parentless cells retain their cell-local mirror until a future cell scheduler authority is introduced. H1 Cartesian CFL widths are derived from persisted cell centers/config and are not serialized as scratch caches. |

Additive softening sidecar persistence:
- Snapshot `/PartTypeN/GravitySofteningComoving` (`float64`, optional; per-particle, comoving units) and `/PartTypeN/HasGravitySofteningOverride` (`uint8`, optional) are diagnostics/interchange mirrors.
- Restart `/state/particle_sidecar/gravity_softening_comoving` (`float64`, optional; per-particle, comoving units) and `/state/particle_sidecar/has_gravity_softening_override` (`uint8`, optional) are exact-continuation lanes. The mask is the authority for per-particle override truth; a value lane without a mask is a materialized default/diagnostic mirror and does not create overrides.

Snapshot and restart intentionally remain separate contracts: snapshot is analysis/interchange oriented;
restart is execution-resume oriented.


## Stage 2 timestep-authority schema note (2026-05-11)

Historical Stage 2 scheduler-authority documentation did not change snapshot/restart/provenance schemas. H2.4 historical material referenced `cosmosim_restart_v17`; v19 introduced authoritative gas-cell scheduling and the active schema is now `cosmosim_restart_v20`, retaining those identity records while adding checkpoint-authoritative gravity force caches. The compatibility behavior is explicit: restart payloads retain `ParticleSoa::time_bin` and `CellSoa::time_bin` as mirrors for corruption detection, reject stale mirror conflicts against scheduler truth, and rebuild valid parent-backed mirrors from scheduler state on import. Gas-cell parent lineage is optional metadata; parentless cells keep cell-local hydro velocity and timestep mirror lanes without particle velocity access.

## H1 Hydro Restart Geometry Note

Schema v15 persists the restart-authoritative hydro geometry inputs needed to rebuild
Cartesian gas-cell geometry deterministically: cell centers, gas-cell mirror lanes,
authoritative `/state/gas_cell_identity` records, cell `patch_index`, and `PatchSoa`
descriptors including patch ownership. Restart validation rejects patch ranges that
overlap, omit cells, exceed the cell arrays, disagree with per-cell `patch_index`, or
disagree with identity-record `owning_patch_id`. Transient hydro reconstruction scratch,
ghost-fill buffers, face flux arrays, and derived CFL width caches remain outside restart
payloads.

## 3) Provenance payload

`ProvenanceRecord` persists:

- schema tag (`provenance_v5`)
- config schema identity (`config_schema_name`, `config_schema_version`)
- audit payloads (`raw_input_config`, `normalized_config`, `derived_runtime_state`)
- deterministic normalized config hash (`normalized_config_hash_hex`)
- build identity (`git_sha`, compiler id/version, build preset, feature flags)
- deterministic config hash
- UTC timestamp and hardware summary
- rank attribution (`author_rank`)
- gravity/TreePM reproducibility contract:
  - controls: `gravity_treepm_pm_grid`, `gravity_treepm_assignment_scheme`,
    `gravity_treepm_window_deconvolution`, `gravity_treepm_asmth_cells`,
    `gravity_treepm_rcut_cells`, `gravity_treepm_update_cadence_steps`,
    `gravity_treepm_pm_decomposition_mode`, `gravity_treepm_tree_exchange_batch_bytes`
  - derived scales: `gravity_treepm_mesh_spacing_mpc_comoving` (`Δmesh`),
    `gravity_treepm_split_scale_mpc_comoving` (`r_s`),
    `gravity_treepm_cutoff_radius_mpc_comoving` (`r_cut`)
  - softening/backend: `gravity_softening_policy`, `gravity_softening_kernel`,
    `gravity_softening_epsilon_kpc_comoving`, `gravity_pm_fft_backend`
  - restart/debug continuation metadata:
    - `gravity_treepm_decomposition_epoch`
    - `gravity_treepm_restart_world_size`
    - `gravity_treepm_restart_pm_grid`
    - `gravity_treepm_restart_slab_signature`
    - `gravity_treepm_restart_kick_opportunity`
    - `gravity_treepm_restart_field_version`
    - `gravity_treepm_long_range_restart_policy`
  - zoom-gravity metadata:
    - `zoom_long_range_strategy`
    - `zoom_region_center_{x,y,z}_mpc_comoving`
    - `zoom_region_radius_mpc_comoving`
    - `zoom_focused_pm_grid`
    - `zoom_contamination_radius_mpc_comoving`

## 4) Naming and stability conventions

- Output stems are restricted to `[A-Za-z0-9_-]` for deterministic naming.
- Normalized config snapshots are written with run outputs.
- Schema names/versions are part of compatibility behavior and must be updated intentionally.
- Reference workflow operational report (`reference_operational_events.json`) is schema-tagged (`schema_version = 1`) and includes:
  - run label,
  - provenance config hash linkage (`provenance_config_hash_hex`),
  - severity-count summary and status,
  - structured events (`event_kind`, `severity`, `subsystem`, optional step/time/scale context, message, payload map).
  - gravity runtime events include:
    - `gravity.treepm_setup` (one setup event with PMGRID, assignment, deconvolution,
      ASMTH/RCUT controls, derived `Δmesh`/`r_s`/`r_cut`, cadence, softening policy/kernel,
      FFT backend),
    - `gravity.pm_long_range_field` (refresh/reuse event per gravity kick opportunity carrying
      the same PM contract plus cadence state/version payload).
    - `gravity.zoom_force_diagnostics` (per-kick zoom decomposition norms and low-resolution contamination counters).
    - `gravity.health_check` (targeted warning/fatal event for explicit gravity-state violations;
      fatal events are not downgraded into generic diagnostics).
    - `gravity.health_summary` (per-kick gravity health counters summarizing cheap always-on checks
      and heavy reference-only checks when policy enabled).

## 4.1) Diagnostics bundle maturity metadata

Diagnostics bundles (`<diagnostics_stem>_<diagnostic_class>_step_*.json`) now include explicit maturity metadata:

- `diagnostic_class`: cadence bucket (`run_health`, `science_light`, `science_heavy`)
- `diagnostics_execution_policy`: active policy (`run_health_only`, `run_health_and_light_science`, `all_including_provisional`)
- `diagnostic_records[]`: per-output metadata
  - `name`
  - `tier` (`infrastructure_health`, `validated_science`, `reference_science`)
  - `maturity` (`production`, `validated`, `provisional`)
  - `scalability` (`cheap`, `moderate`, `heavy_reference`)
  - `executed`
  - `policy_note`

Current intended classification:

- Production infrastructure health: `run_health_counters`
- Production infrastructure health: `gravity_health_summary`
- Validated lightweight science: `star_formation_history`, `angular_momentum_budget`, `gas_xy_slice_density`, `gas_xy_projection_density`
- Provisional heavy reference-only: `power_spectrum` (disabled unless `diagnostics_execution_policy = all_including_provisional`)

Run-health payload now includes additive gravity integrity counters:

- `gravity_softening_sidecar_size_ok`
- `non_finite_gravity_softening`
- `non_positive_particle_mass`

## 5) Change procedure for schema-affecting work

When changing snapshot/restart/provenance fields:

1. Update the corresponding interface headers and implementation.
2. Update `docs/configuration.md` if config keys or normalized text semantics changed.
3. Update validation expectations in `docs/validation_plan.md`.
4. Add/update tests in `tests/unit` + `tests/integration` + `tests/validation` as applicable.
5. Record rationale in `docs/architecture/decision_log.md`.

## Compatibility notes (2026-04-20)

- Snapshot schema was intentionally bumped to `gadget_arepo_v4` (`schema_version = 4`)
  to add optional per-particle softening sidecar dataset (`GravitySofteningComoving`) per particle group.
- No external `/PartType*` dataset names were changed.
- Restart schema version/name are now `cosmosim_restart_v20`, version `20`. It retains the v19 cell-local gas velocity, stable gas-cell identity, and authoritative gas-cell scheduler state; it additionally persists `/gravity_force_cache` and `IntegratorState::pm_refresh_enabled` so a workflow restart can apply the same next KDK pre-kick force state.
- Restart contract enforcement was tightened: missing continuation-critical metadata, a missing or wrong root file kind, or missing output-cadence state now fails fast with explicit path-aware errors instead of producing weak checkpoints.
- Restart schema is `cosmosim_restart_v20`; distributed TreePM state and the restart-authoritative gravity force cache are persisted under restart-only data and covered by restart integrity hashing.
- The reader accepts the documented legacy `cosmosim_restart_v14` particle-bound import path
  by materializing `/state/gas_cell_identity` from
  `/state/gas_cells/{gas_cell_id,parent_particle_id}` with
  `has_parent_particle=true`; it requires `gas_cell_id == parent_particle_id != 0` and
  does not reinterpret `parent_particle_id=0` as absent in old files.
- New v15 files require `/state/gas_cell_identity/has_parent_particle`; malformed identity
  records with duplicate or zero `gas_cell_id`, sparse `local_cell_row`, invalid parent
  flags, or invalid patch ownership are rejected.
- Restart v6 compatibility behavior is tightened without changing payload fields: non-empty
  `CellSoa::time_bin` mirrors must match the scheduler `bin_index` of each gas cell's
  parent particle, and restart import rebuilds cells from that parent-particle scheduler
  mapping. Mixed states where `cells.size() != particles.size()` are no longer an
  exact-size validation bypass.
- Diagnostics maturity metadata is additive to diagnostics JSON bundles and does not alter snapshot/restart/provenance schema compatibility.
- Provenance payload now includes additive zoom-gravity metadata and contamination-radius contract keys;
  this is a backward-compatible extension to `provenance_v5`.


### AMR patch geometry restart lanes (v16)

PatchSoa now persists restart-authoritative AMR patch geometry lanes: `parent_patch_id`, `morton_key`, `origin_x/y/z_comoving`, `extent_x/y/z_comoving`, and `cell_dim_x/y/z`. Production AMR hydro requires these lanes to be explicit; legacy restart inputs without them are accepted only as non-AMR/legacy patch states and do not silently enter the production AMR hydro path.


## AMR temporal restart state (v18)

`cosmosim_restart_v19` introduced `/state/amr_temporal_boundary_history` for active local AMR coarse temporal intervals. The current `cosmosim_restart_v20` retains it and adds `/gravity_force_cache` for exact workflow KDK continuation. Neither is part of analysis snapshots. v17 pending flux-register restart state remains supported as a legacy read path.
