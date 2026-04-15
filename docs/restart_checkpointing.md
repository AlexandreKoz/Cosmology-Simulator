# Restart and Checkpoint System

## Scope and schema

CosmoSim restart checkpoints are **exact-continuation artifacts** and intentionally richer than analysis snapshots.
The restart schema (`cosmosim_restart_v3`) persists:

- full `SimulationState` hot/cold SoA lanes,
- `StateMetadata` blob,
- module sidecars (`ModuleSidecarRegistry`) with per-module schema versions,
- `IntegratorState`,
- hierarchical scheduler persistent state (`TimeBinPersistentState`),
- normalized config text/hash and provenance payload,
- payload integrity hash (FNV-1a 64-bit).

By design, this differs from GADGET/AREPO-style analysis snapshots where scheduler internals and opaque sidecars are not mandatory.

## File format and compatibility

- Format: HDF5 (`writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`).
- Schema version gate: `isRestartSchemaCompatible(file_schema_version)`.
- Current compatibility policy: exact version match (`3`). Older `v2` restart files are intentionally rejected because they omit required stellar-evolution continuation state.

## Atomic write semantics

Writers always emit to `<final>.tmp` first and only then rename into final path.
Finalization uses a direct rename of the temp artifact (no pre-remove step), so the old
restart path is never explicitly deleted before replacement. This preserves atomic replace
behavior on filesystems where `rename` is atomic.

## Integrity and provenance

- `restartPayloadIntegrityHash` hashes state+integrator+scheduler+config text/hash+provenance.
- The hash covers particle lanes, sidecars, species counts, full star/BH/tracer sidecars (including stellar-evolution cumulative lanes), integrator
  time-bin context, and scheduler persistent arrays (`bin_index`, `next_activation_tick`,
  `active_flag`, `pending_bin_index`).
- Tracer restart payload includes host-coupling lanes (`host_cell_index`, `mass_fraction_of_host`,
  `last_host_mass_code`, `cumulative_exchanged_mass_code`) for deterministic continuation.
- Writer stores both integer and hex payload integrity hashes.
- Reader recomputes hash and rejects mismatches.
- Provenance is serialized with the checkpoint for continuation auditing.

## Exact-restart completeness checklist

The checklist exposed by `exactRestartCompletenessChecklist()` and enforced by
restart write/read checks is:

1. `simulation_state_lanes_and_metadata`
2. `module_sidecars_with_schema_versions`
3. `integrator_state`
4. `scheduler_persistent_state`
5. `normalized_config_text_and_hash`
6. `provenance_record`
7. `payload_integrity_hash_and_hex`

`validateContinuationMetadata(...)` now requires non-empty normalized config text, a normalized config hash that matches the text, and a provenance config hash that matches the normalized config hash before hashing/writing restart payloads.

## Exactness policy

Restart checkpoints target exact continuation on the same build/feature path. Current tests
expect bitwise equality for persisted arrays and strict scalar equality, with no tolerance
window except tiny floating-point assertion slack in selected test comparisons (`< 1e-15`)
to avoid false negatives from host-side literal conversions.

## Parallel and scale-up note

The current implementation is single-rank write/read but keeps explicit scheduler/rank-related sidecars,
which preserves a direct path to coordinated rank-safe checkpointing in future MPI modes.
