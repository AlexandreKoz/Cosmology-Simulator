# Output schema guide

This guide documents snapshot, restart, and provenance payload contracts currently implemented in CosmoSim.

Authoritative interfaces:

- `include/cosmosim/io/snapshot_hdf5.hpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `include/cosmosim/core/provenance.hpp`
- `include/cosmosim/core/profiling.hpp` (operational run-event report)

## 1) Snapshot schema (GADGET/AREPO interoperability)

Current schema identity:

- `schema_name = gadget_arepo_v3`
- `schema_version = 3`

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

- `name = cosmosim_restart_v5`
- `version = 5`

Restart payload includes:

- `SimulationState`
- `IntegratorState`
- `HierarchicalTimeBinScheduler` persistent state (`TimeBinPersistentState`)
- distributed TreePM continuation metadata (`DistributedRestartState`)
- normalized config text and normalized config hash
- provenance record
- payload integrity hash

Compatibility rule is explicit through `isRestartSchemaCompatible(version)`.

## 2.1) Field ownership table (snapshot vs restart)

| Ownership | Fields |
|---|---|
| Shared metadata contract | normalized config text/hash, provenance payload, schema identity |
| Snapshot-only (interoperable science output) | `/Header` cosmology attrs, `/PartTypeN` particle datasets, read aliases (`Position`, `VEL`, `ID`, etc.) |
| Restart-only (exact continuation state) | full `SimulationState` hot/cold lanes, `StateMetadata`, module sidecars + schema versions, `IntegratorState`, scheduler persistent state (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`), distributed TreePM restart state (`decomposition_epoch`, owning-rank table, PM slab layout, cadence/long-range metadata, restart policy), payload integrity hashes |

Additive softening sidecar persistence:
- Snapshot `/PartTypeN/GravitySofteningComoving` (`float64`, optional; per-particle, comoving units).
- Restart `/state/particle_sidecar/gravity_softening_comoving` (`float64`, optional; per-particle, comoving units).

Snapshot and restart intentionally remain separate contracts: snapshot is analysis/interchange oriented;
restart is execution-resume oriented.

## 3) Provenance payload

`ProvenanceRecord` persists:

- schema tag (`provenance_v3`)
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
- Validated lightweight science: `star_formation_history`, `angular_momentum_budget`, `gas_xy_slice_density`, `gas_xy_projection_density`
- Provisional heavy reference-only: `power_spectrum` (disabled unless `diagnostics_execution_policy = all_including_provisional`)

## 5) Change procedure for schema-affecting work

When changing snapshot/restart/provenance fields:

1. Update the corresponding interface headers and implementation.
2. Update `docs/configuration.md` if config keys or normalized text semantics changed.
3. Update validation expectations in `docs/validation_plan.md`.
4. Add/update tests in `tests/unit` + `tests/integration` + `tests/validation` as applicable.
5. Record rationale in `docs/architecture/decision_log.md`.

## Compatibility notes (2026-04-20)

- Snapshot schema was intentionally bumped to `gadget_arepo_v3` (`schema_version = 3`)
  to add optional per-particle softening sidecar dataset (`GravitySofteningComoving`) per particle group.
- No external `/PartType*` dataset names were changed.
- Restart schema version/name are now `cosmosim_restart_v5`, version `5`, because restart payloads persist optional particle softening sidecar state.
- Restart contract enforcement was tightened: missing continuation-critical metadata
  now fails fast with explicit errors instead of producing weak checkpoints.
- Restart schema is `cosmosim_restart_v5`; distributed TreePM state is persisted under restart-only
  data and covered by restart integrity hashing.
- Diagnostics maturity metadata is additive to diagnostics JSON bundles and does not alter snapshot/restart/provenance schema compatibility.
