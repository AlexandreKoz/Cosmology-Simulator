# Restart and Checkpoint System

## Scope and schema

CosmoSim restart checkpoints are **exact-continuation artifacts** and intentionally richer than analysis snapshots.
The restart schema (`cosmosim_restart_v18`) persists:

- full `SimulationState` hot/cold SoA lanes (through a narrow `RestartPersistentStateView`),
- `StateMetadata` blob,
- hydro geometry ownership lanes: cell centers, gas-cell identity, cell-local velocity, and thermodynamic fields,
  cell `patch_index`, and `PatchSoa` descriptors (`patch_id`, `level`, `first_cell`,
  `cell_count`, `owning_rank`),
- authoritative gas-cell identity records under `/state/gas_cell_identity`, including
  `gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`,
  `local_cell_row`, `local_row_reconstruction_policy=explicit_dense_local_cell_row`,
  and `identity_generation_at_write` as an audit stamp,
- module sidecars (`ModuleSidecarRegistry`) with per-module schema versions,
- `IntegratorState`,
- hierarchical scheduler persistent state (`TimeBinPersistentState`),
- distributed TreePM continuation state (`distributed_gravity_state`, schema-versioned),
- normalized config text/hash and provenance payload,
- output cadence persistence state for deterministic future snapshot/restart naming and cadence decisions,
- root file-kind metadata (`cosmosim_file_kind=restart_checkpoint`),
- compact `/restart_diagnostics` audit metadata,
- payload integrity hash (FNV-1a 64-bit with explicit string/vector length delimiters).

By design, this differs from GADGET/AREPO-style analysis snapshots where scheduler internals and opaque sidecars are not mandatory.

Restart write payloads now carry `RestartPersistentStateView` (`persistent_state.simulation_state`) rather than a broad ad-hoc object pointer. This creates an explicit outer-boundary guard: transient runtime owners such as `TransientStepWorkspace`, `HydroScratchBuffers`, PM work arrays, TreePM traversal scratch, MPI send/recv buffers, and output staging buffers are out-of-scope for restart traversal by type.


## Restart-safe boundary contract

Restart checkpoints may only be written from a completed, globally coherent restart boundary. The runtime predicate `core::evaluateRestartBoundary(...)` is the narrow contract used by workflow output dispatch and HDF5 restart payload validation. It rejects half-step KDK states, local active-bin substeps, non-restart-safe boundary kinds, and PM refresh transitions with an uncommitted long-range refresh event.

Failed restart requests are hard errors, not warnings or silent skips. The diagnostic records the current and last completed boundary kind, `inside_kdk_step`, `last_completed_restart_safe`, local-substep activity, PM refresh legality/commit-pending state, `step_index`, and scheduler tick when available. Intentionally represented half-step or local-substep restart is not implemented in schema v15, so those states must not be serialized as persistent truth.


## Restart diagnostics metadata

Schema v15 retains the compact `/restart_diagnostics` group. It records the schema identity, current and last completed boundary kind, restart-safe decision, scheduler tick/bin/active/pending counts, PM cadence/field-version summary, output cadence summary, and stochastic module count. These fields are audit metadata only: authoritative continuation remains the serialized `SimulationState`, `IntegratorState`, scheduler arrays, output cadence state, stochastic state, distributed TreePM state, normalized config, provenance, and integrity hash.

## File format and compatibility

- Format: HDF5 (`writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`).
- Root file-kind gate: restart readers require `cosmosim_file_kind=restart_checkpoint` and reject ordinary `science_snapshot` files before reading runtime truth.
- Schema version gate: `isRestartSchemaCompatible(file_schema_version)`.
- Current compatibility policy: write current v18; read v18 plus documented legacy v17/v16/v15/v14 paths.
- v14 compatibility materializes `/state/gas_cell_identity` from
  `/state/gas_cells/{gas_cell_id,parent_particle_id}` with `has_parent_particle=true`
  and requires `gas_cell_id == parent_particle_id != 0` for every cell. It does not
  reinterpret `parent_particle_id=0` as absent in old files.
- `v5` and older restart files are intentionally rejected because they can omit exact sidecar truth lanes required by v6, including materialized softening override mask/value datasets.

## Atomic write semantics

Writers always emit to `<final>.tmp` first and only then rename into final path.
Finalization uses a direct rename of the temp artifact (no pre-remove step), so the old
restart path is never explicitly deleted before replacement. This preserves atomic replace
behavior on filesystems where `rename` is atomic.

## Integrity and provenance

- `restartPayloadIntegrityHash` hashes state+integrator+scheduler+config text/hash+provenance and first validates that derived particle `time_bin` mirrors match scheduler `bin_index` authority. For gas cells with a local parent particle, derived cell `time_bin` mirrors are validated through the parent scheduler entry; parentless gas-cell time-bin mirrors remain cell-local diagnostics until a future cell scheduler authority is introduced. Restart validation also requires hydro patch descriptors to cover every gas-cell row exactly once and requires each cell's `patch_index` to agree with the serialized patch range. H1 Cartesian CFL metadata remains derived from persistent cell centers/config and is not serialized as a scratch cache.
- The hash covers particle lanes, sidecars, gas-cell mirror lanes, authoritative gas-cell identity records (`gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`, `local_cell_row`, and reconstruction policy), gas-cell velocity and thermodynamic lanes, species counts, full star/BH/tracer sidecars (including stellar-evolution cumulative lanes), integrator
  time-bin context, scheduler persistent arrays (`bin_index`, `next_activation_tick`,
  `active_flag`, `pending_bin_index`), output cadence persistence fields, and the full serialized provenance payload. State `time_bin` arrays remain hash inputs for corruption detection, but restart continuation imports scheduler state and rebuilds mirrors rather than treating mirrors as fallback authority.
- Because gravity TreePM metadata now lives in `ProvenanceRecord`, restart artifacts carry
  and integrity-protect the exact PM controls (`pm_grid`, assignment, deconvolution,
  `asmth_cells`, `rcut_cells`, cadence), derived scales (`Δmesh`, `r_s`, `r_cut`),
  softening policy/kernel/epsilon, and FFT backend name.
- Distributed continuation metadata is persisted under `/distributed_gravity/state` (text-encoded
  `parallel::DistributedRestartState`) and includes:
  - decomposition epoch,
  - per-particle owning rank mapping (`owning_rank_by_item`),
  - PM grid shape and slab ownership table (`pm_slab_begin_x_by_rank`, `pm_slab_end_x_by_rank`),
  - kick-opportunity cadence state (`gravity_kick_opportunity`, `pm_update_cadence_steps`),
  - long-range field refresh metadata (`long_range_field_version`,
    `last_long_range_refresh_opportunity`, `long_range_field_built_step_index`,
    `long_range_field_built_scale_factor`),
  - restart long-range policy (`long_range_restart_policy`).
  Distributed floating-point restart metadata is text-encoded with enough precision for exact read-back comparisons.
- **Policy:** restart continuation uses `long_range_restart_policy=deterministic_rebuild`.
  Cached PM long-range field arrays are not serialized; on resume, PM long-range state is rebuilt
  deterministically on the next refresh opportunity. This is contractually explicit in restart payload + provenance.
- Tracer restart payload includes host-coupling lanes (`host_cell_index`, `mass_fraction_of_host`,
  `last_host_mass_code`, `cumulative_exchanged_mass_code`) for deterministic continuation.
- Writer stores both integer and hex payload integrity hashes.
- Reader recomputes hash and rejects mismatches.
- Provenance is serialized with the checkpoint for continuation auditing.

## Exact-restart completeness checklist

The checklist exposed by `exactRestartCompletenessChecklist()` and enforced by
restart write/read checks is:

1. `simulation_state_lanes_and_metadata`
2. `particle_identity_softening_and_drift_epoch_lanes`
3. `gas_cell_identity_lanes`
4. `hydro_geometry_patch_state`
5. `species_specific_sidecars`
6. `module_sidecars_with_schema_versions`
7. `integrator_state`
8. `integrator_owned_pm_sync_state`
9. `scheduler_persistent_state`
10. `output_cadence_persistent_state`
11. `stochastic_module_persistent_state`
12. `restart_diagnostics_summary`
13. `distributed_gravity_state`
14. `normalized_config_text_and_hash`
15. `provenance_record`
16. `payload_integrity_hash_and_hex`
17. `amr_pending_flux_register_state`

`validateContinuationMetadata(...)` now requires non-empty normalized config text, a normalized config hash that matches the text, and a provenance config hash that matches the normalized config hash before hashing/writing restart payloads.

## Exactness policy

Restart checkpoints target exact continuation on the same build/feature path. Current tests
expect bitwise equality for persisted arrays and strict scalar equality, with no tolerance
window except tiny floating-point assertion slack in selected test comparisons (`< 1e-15`)
to avoid false negatives from host-side literal conversions.

## Runtime-truth round-trip and compatibility behavior

Restart round-trip invariant coverage now explicitly checks that:

- particle IDs and canonical particle order survive reload,
- species counts/index reconstruction and gas identity-by-particle-ID survive reload,
- scheduler persistent state survives and reconstructs the same active set at resume tick,
- particle `time_bin` mirrors are rejected when stale and rebuilt from scheduler state on successful read,
- particle-bound gas-cell `time_bin` mirrors are rejected when stale against their parent gas particle scheduler entry and are rebuilt from that mapping on successful read,
- hydro patch geometry state is rejected when patch ranges overlap, omit gas cells, extend beyond
  the cell arrays, or disagree with per-cell `patch_index`,
- per-particle softening overrides survive and continue to take precedence over species/global defaults,
- normalized config text/hash and provenance config hash remain aligned after reload.

Compatibility policy for legacy restart payloads is explicit:

- Missing required continuation fields (for example scheduler persistent lanes such as
  `pending_bin_index`) are rejected with clear read errors.
- Stale particle `time_bin` mirrors that disagree with `/scheduler/bin_index` are rejected by hash/write/read validation before exact continuation is accepted. Stale particle-bound gas-cell `time_bin` mirrors are checked against the parent gas particle scheduler entry. This originated as a v6 compatibility check; the current schema is v15 after H2 gas-cell identity promotion added explicit identity records without making diagnostics authoritative truth.
- v15 malformed identity records are rejected on read: duplicate or zero `gas_cell_id`,
  missing `/state/gas_cell_identity/has_parent_particle`, non-binary parent flags,
  `parent_particle_id=0` when `has_parent_particle=true`, nonzero parent values when
  `has_parent_particle=false`, sparse or duplicate `local_cell_row`, and identity
  `owning_patch_id` values that do not match each cell's owning `PatchSoa::patch_id`.
- The v6 reader requires `/state/particle_sidecar/gravity_softening_comoving` and `/state/particle_sidecar/has_gravity_softening_override` datasets to be present. The mask remains authoritative; empty datasets encode no materialized per-particle softening lane, while populated datasets must satisfy `ParticleSidecar` ownership invariants. Missing lanes are rejected instead of being interpreted as implicit defaults.

## Parallel and scale-up note

The current implementation is single-rank write/read but keeps explicit scheduler/rank-related sidecars,
which preserves a direct path to coordinated rank-safe checkpointing in future MPI modes.

## Runtime-state exactness check for reference workflow restart verification

The reference workflow now treats restart verification as exact runtime-state continuation, not a count-only smoke test. The verification path compares all restart-authoritative state families: particle hot SoA lanes, particle metadata/softening override lanes, gas-cell identity and hydro fields, patch state, star/black-hole/tracer sidecars, species count ledgers, opaque module sidecar payloads, scheduler state, provenance/config metadata, and distributed gravity restart metadata.

Module sidecar payload length is explicitly included in the restart payload integrity hash before payload bytes. This preserves hash boundary clarity across adjacent module payloads and makes omitted or reshaped module truth detectable by the hash.


### AMR patch geometry restart lanes (v16)

PatchSoa now persists restart-authoritative AMR patch geometry lanes: `parent_patch_id`, `morton_key`, `origin_x/y/z_comoving`, `extent_x/y/z_comoving`, and `cell_dim_x/y/z`. Production AMR hydro requires these lanes to be explicit; legacy restart inputs without them are accepted only as non-AMR/legacy patch states and do not silently enter the production AMR hydro path.

## AMR hydro restart-equivalence coverage

The H3 hardening pass added `integration_restart_equivalence_amr_hydro`, an HDF5-gated integration test that compares a direct local AMR hydro continuation against a checkpoint/reload/continue path. The test exercises production AMR coverage, explicit PatchSoa geometry lanes, stable `GasCellIdentityMap` records, ghost-fill/reflux safety, scheduler state, integrator state, output cadence state, and stochastic persistent state.

The restart-equivalence harness now compares explicit AMR patch lanes (`parent_patch_id`, `morton_key`, origins, extents, cell dimensions), `GasCellIdentityMap` records, and the gas identity generation value. Restart read restores the serialized gas identity generation with `GasCellIdentityMap::assignWithGeneration(...)` so a checkpoint/reload path does not masquerade as a runtime identity rebuild.

This proves local HDF5 AMR hydro restart equivalence for the exercised synchronized-sweep scenario. It does not prove MPI restart, restart after AMR migration, Berger-Colella subcycling, or replay of persistent pending flux registers.

---

## AMR pending flux-register restart lanes (v17)

Restart schema v17 introduced pending AMR flux-register state; current v18 checkpoints retain it under:

```text
/state/amr_pending_flux_registers
```

The group persists stable identity and coverage metadata for deferred reflux records, including:

- `register_key`
- `coarse_patch_id`
- `coarse_gas_cell_id`
- `coarse_cell_index`
- `level`, `axis`, `orientation`
- expected/coarse/fine area coverage
- interval start/end and coarse timestep
- expected/completed fine substeps and substep coverage mask
- coarse/fine face counts
- gas-cell identity generation and patch geometry generation
- coarse/fine flux-integrated mass, momentum, and total-energy contributions

Old restart files without this group are treated as legacy inputs with an empty pending store only at a synchronized point. A pre-v18 payload that does contain pending deferred-reflux records is rejected by the temporal-history reader, because it cannot prove the matching temporal boundary history needed for safe mid-subcycle continuation. Current v17 restart schema validation required the group and datasets so a new v17 restart could not accidentally omit pending deferred reflux truth.

The restart payload integrity hash includes the pending-register store for v17 payloads. The HDF5-gated `integration_restart_equivalence_amr_flux_registers` test writes an incomplete pending register before restart, reads it back, merges the missing fine contribution, applies completed reflux, and compares direct vs restart continuation.

This proves the exercised single-rank pending-register restart path. It does not prove MPI AMR restart, restart after AMR patch migration, or scheduler-owned Berger-Colella subcycling.

---

## AMR temporal boundary-history restart lanes (v18)

Restart schema `cosmosim_restart_v18` adds `/state/amr_temporal_boundary_history` for open local AMR
coarse intervals. Each history record persists patch ID/level, geometry fingerprint, gas identity generation,
interval start/end, completion state, and stable-ID patch-local conserved start/end records. The restart
payload integrity hash includes this store for v18 payloads.

`integration_restart_equivalence_amr_temporal_ghosts` checkpoints after coarse end-state capture and a
first fine temporal ghost use, with an incomplete pending reflux record present. The restarted branch reads
the history, uses it for the midpoint coarse-to-fine ghost, completes the deferred flux contribution, and
matches the direct continuation. This is single-rank local coverage only.

Legacy schemas before v18 do not contain temporal history. They may be read only for states that do not
require resumption inside an active temporal AMR interval; legacy files are not a supported way to resume
one. The v18 reader conservatively rejects a pre-v18 payload that still contains pending AMR flux-register
records, because that observable deferred-synchronization state cannot prove the required time-aligned
coarse boundary history exists.
