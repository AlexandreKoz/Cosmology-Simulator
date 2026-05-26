# Restart and Checkpoint System

## Scope and schema

CosmoSim restart checkpoints are **exact-continuation artifacts** and intentionally richer than analysis snapshots.
The restart schema (`cosmosim_restart_v11`) persists:

- full `SimulationState` hot/cold SoA lanes (through a narrow `RestartPersistentStateView`),
- `StateMetadata` blob,
- module sidecars (`ModuleSidecarRegistry`) with per-module schema versions,
- `IntegratorState`,
- hierarchical scheduler persistent state (`TimeBinPersistentState`),
- distributed TreePM continuation state (`distributed_gravity_state`, schema-versioned),
- normalized config text/hash and provenance payload,
- payload integrity hash (FNV-1a 64-bit with explicit string/vector length delimiters).

By design, this differs from GADGET/AREPO-style analysis snapshots where scheduler internals and opaque sidecars are not mandatory.

Restart write payloads now carry `RestartPersistentStateView` (`persistent_state.simulation_state`) rather than a broad ad-hoc object pointer. This creates an explicit outer-boundary guard: transient runtime owners such as `TransientStepWorkspace`, `HydroScratchBuffers`, PM work arrays, TreePM traversal scratch, MPI send/recv buffers, and output staging buffers are out-of-scope for restart traversal by type.

## File format and compatibility

- Format: HDF5 (`writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`).
- Schema version gate: `isRestartSchemaCompatible(file_schema_version)`.
- Current compatibility policy: exact version match (`11`).
- `v5` and older restart files are intentionally rejected because they can omit exact sidecar truth lanes required by v6, including materialized softening override mask/value datasets.

## Atomic write semantics

Writers always emit to `<final>.tmp` first and only then rename into final path.
Finalization uses a direct rename of the temp artifact (no pre-remove step), so the old
restart path is never explicitly deleted before replacement. This preserves atomic replace
behavior on filesystems where `rename` is atomic.

## Integrity and provenance

- `restartPayloadIntegrityHash` hashes state+integrator+scheduler+config text/hash+provenance and first validates that derived particle `time_bin` mirrors match scheduler `bin_index` authority. For particle-bound gas cells, derived cell `time_bin` mirrors are validated through each cell's parent gas particle scheduler entry rather than by cell-row count equality.
- The hash covers particle lanes, sidecars, gas-cell identity lanes, species counts, full star/BH/tracer sidecars (including stellar-evolution cumulative lanes), integrator
  time-bin context, scheduler persistent arrays (`bin_index`, `next_activation_tick`,
  `active_flag`, `pending_bin_index`), and the full serialized provenance payload. State `time_bin` arrays remain hash inputs for corruption detection, but restart continuation imports scheduler state and rebuilds mirrors rather than treating mirrors as fallback authority.
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
2. `particle_identity_and_softening_override_lanes`
3. `gas_cell_identity_lanes`
4. `species_specific_sidecars`
5. `module_sidecars_with_schema_versions`
6. `integrator_state`
7. `scheduler_persistent_state`
8. `distributed_gravity_state`
9. `normalized_config_text_and_hash`
10. `provenance_record`
11. `payload_integrity_hash_and_hex`

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
- per-particle softening overrides survive and continue to take precedence over species/global defaults,
- normalized config text/hash and provenance config hash remain aligned after reload.

Compatibility policy for legacy restart payloads is explicit:

- Missing required continuation fields (for example scheduler persistent lanes such as
  `pending_bin_index`) are rejected with clear read errors.
- Stale particle `time_bin` mirrors that disagree with `/scheduler/bin_index` are rejected by hash/write/read validation before exact continuation is accepted. Stale particle-bound gas-cell `time_bin` mirrors are checked against the parent gas particle scheduler entry. This originated as a v6 compatibility check; the current schema is v11 after Stage 6 removed transient hydro reconstruction gradients from restart truth.
- The v6 reader requires `/state/particle_sidecar/gravity_softening_comoving` and `/state/particle_sidecar/has_gravity_softening_override` datasets to be present. The mask remains authoritative; empty datasets encode no materialized per-particle softening lane, while populated datasets must satisfy `ParticleSidecar` ownership invariants. Missing lanes are rejected instead of being interpreted as implicit defaults.

## Parallel and scale-up note

The current implementation is single-rank write/read but keeps explicit scheduler/rank-related sidecars,
which preserves a direct path to coordinated rank-safe checkpointing in future MPI modes.

## Runtime-state exactness check for reference workflow restart verification

The reference workflow now treats restart verification as exact runtime-state continuation, not a count-only smoke test. The verification path compares all restart-authoritative state families: particle hot SoA lanes, particle metadata/softening override lanes, gas-cell identity and hydro fields, patch state, star/black-hole/tracer sidecars, species count ledgers, opaque module sidecar payloads, scheduler state, provenance/config metadata, and distributed gravity restart metadata.

Module sidecar payload length is explicitly included in the restart payload integrity hash before payload bytes. This preserves hash boundary clarity across adjacent module payloads and makes omitted or reshaped module truth detectable by the hash.
