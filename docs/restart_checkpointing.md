# Restart and Checkpoint System

## Scope and schema

CosmoSim restart checkpoints are **exact-continuation artifacts** and intentionally richer than analysis snapshots.
The restart schema (`cosmosim_restart_v21`) persists:

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
- `IntegratorState`, including the `pm_refresh_enabled` policy bit,
- restart-authoritative `/gravity_force_cache` acceleration triplets used by the next KDK pre-kick,
- hierarchical scheduler persistent state (`TimeBinPersistentState`),
- distributed TreePM continuation state (`distributed_gravity_state`, schema-versioned),
- normalized config text/hash and provenance payload,
- output cadence persistence state for deterministic future snapshot/restart naming and cadence decisions,
- root file-kind metadata (`cosmosim_file_kind=restart_checkpoint`),
- compact `/restart_diagnostics` audit metadata,
- payload integrity hash (FNV-1a 64-bit with explicit string/vector length delimiters).

By design, this differs from GADGET/AREPO-style analysis snapshots where scheduler internals and opaque sidecars are not mandatory.

## Ordered output timeline events (v21)

Schema v21 adds `snapshot_interval_time_code` and `next_snapshot_time_code` to
`/output_cadence`, with the latter also mirrored in `/restart_diagnostics` for audit.
The interval is anchored at the configured begin time. `ReferenceWorkflow` clips a
larger KDK interval before it crosses the next event, writes at the resulting safe
boundary, advances the persisted event, and stores the pre-clip interval as the next
step proposal. This keeps uninterrupted and restart continuation schedules aligned.
Step-modulo cadence remains supported and may be enabled alongside code-time cadence.
The integrity hash covers both new fields. Legacy v20 reads set both fields to zero,
which explicitly disables code-time events rather than inventing a future event.

## H1 fixed/patch workflow restart note (v20)

The v20 addition is deliberately narrow: it persists the committed particle and gas-cell gravity acceleration cache at a restart-safe KDK boundary. `ReferenceWorkflow` restores that cache before its next pre-kick, while continuing to rebuild PM mesh work arrays at their next legal refresh surface. This closes the workflow-level restart gap without serializing transient hydro scratch, ghost descriptors, PM meshes, tree traversal workspaces, or MPI buffers.

Restart write payloads now carry `RestartPersistentStateView` (`persistent_state.simulation_state`) rather than a broad ad-hoc object pointer. This creates an explicit outer-boundary guard: transient runtime owners such as `TransientStepWorkspace`, `HydroScratchBuffers`, PM work arrays, TreePM traversal scratch, MPI send/recv buffers, and output staging buffers are out-of-scope for restart traversal by type.

## Gravity force history, PM cadence, and relative-MAC continuation

The `/gravity_force_cache` is the committed dense particle/gas-cell acceleration
state consumed by the next KDK pre-kick. Stable `particle_id` and `gas_cell_id`
lanes make it possible to remap this cache after a row reorder without treating
dense row as physical identity. On import, all components must be finite and
the stable-ID lanes must be unique and exactly cover the current state.

The same committed particle acceleration is the only production source for the
relative-force tree MAC's previous-acceleration magnitude. The workflow uses it
only when the cache is valid, its dense extents match, and its particle-index
generation matches the current `SimulationState`. Otherwise the traversal uses
the deterministic COM-distance fallback. A restart therefore does not invent a
force scale, and a migration/compaction generation change cannot silently feed
stale row-indexed history to the MAC.

Production config currently requires
`numerics.treepm_update_cadence_steps = 1` and
`numerics.hierarchical_max_rung = 0`. Every integrator-issued,
rank-coordinated production force-refresh surface rebuilds PM. The rung-zero
restriction is fail-closed until per-element kick/drift epochs exist for a
correct mixed-rung KDK update. On resume, both particle and gas-cell scheduler
payloads must have `max_bin=0`, all committed bins zero, and every pending bin
zero or unset. Ranks validate this policy before production stepping. The
integrator/distributed restart state remains authoritative for:

- gravity-kick opportunity and refresh decision;
- long-range field version and last refresh opportunity;
- PM field build step and scale factor;
- same-topology PM slab ownership.

PM real/spectral meshes, periodic unwrapped tree coordinates, tree topology,
multipoles, hierarchy packets, and MPI buffers are transient. They are never
restart truth. The PM mesh follows `deterministic_rebuild` at the next legal
refresh; the force cache preserves the already committed kick input. Tree and
hierarchy state are rebuilt from owner particles. Their actual particle-
decomposition epoch is restored from distributed restart truth, while force and
exchange identity follows the resumed field/exchange sequence.

The in-memory TreePM long-range cache validity signature is transient as well.
It covers force epoch, force-evaluation scale factor, `G_code`, split scale,
rectangular box axes, assignment, boundary, PM decomposition mode, and
deconvolution. It is not serialized as a second owner. A newly constructed or
signature-incompatible coordinator therefore rejects an explicit reuse request
coherently; it does not silently reinterpret reuse as refresh. Divergent
explicit refresh/reuse votes fail before PM collectives. Particle-decomposition
epoch is deliberately excluded because PM ownership remains the fixed FFT slab
map; tree protocol records still require that epoch.

The distributed decomposition epoch is a workflow-owned generation, not a
rank-local dense-row counter. It advances once on every rank only when a
measured rebalance commits an actual particle-ownership change, is persisted in
`/distributed_gravity/state` and provenance, and is restored before TreePM use.
The dense-row acceleration cache is invalidated immediately on such a commit,
so a checkpoint written after migration honestly stores `cache.valid=false`
rather than stale row-indexed force history.

`provenance_v6` adds the normalized tree-opening criterion, theta,
relative-force tolerance, and relative-acceleration floor to snapshot/restart
audit metadata. This is an additive provenance schema change; it does not make
those fields separate live authorities or change the v21 restart ownership contract.


## Restart-safe boundary contract

Restart checkpoints may only be written from a completed, globally coherent restart boundary. The runtime predicate `core::evaluateRestartBoundary(...)` is the narrow contract used by workflow output dispatch and HDF5 restart payload validation. It rejects half-step KDK states, local active-bin substeps, non-restart-safe boundary kinds, and PM refresh transitions with an uncommitted long-range refresh event.

Failed restart requests are hard errors, not warnings or silent skips. The diagnostic records the current and last completed boundary kind, `inside_kdk_step`, `last_completed_restart_safe`, local-substep activity, PM refresh legality/commit-pending state, `step_index`, and scheduler tick when available. Intentionally represented half-step or local-substep restart is not implemented in schema v21, so those states must not be serialized as persistent truth.


## Restart diagnostics metadata

Schema v15 retains the compact `/restart_diagnostics` group. It records the schema identity, current and last completed boundary kind, restart-safe decision, scheduler tick/bin/active/pending counts, PM cadence/field-version summary, output cadence summary, and stochastic module count. These fields are audit metadata only: authoritative continuation remains the serialized `SimulationState`, `IntegratorState`, scheduler arrays, output cadence state, stochastic state, distributed TreePM state, normalized config, provenance, and integrity hash.

## File format and compatibility

- Format: HDF5 (`writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`).
- Root file-kind gate: restart readers require `cosmosim_file_kind=restart_checkpoint` and reject ordinary `science_snapshot` files before reading runtime truth.
- Schema version gate: `isRestartSchemaCompatible(file_schema_version)`.
- Current compatibility policy: write current v21; read v21 plus documented v20/v19/v18/v17/v16/v15/v14 paths. v20 has no code-time output-event fields and materializes that lane as disabled. v19 also has no serialized gravity force cache and therefore follows the explicit safe-bootstrap compatibility path rather than claiming bitwise workflow continuation.
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
- `G_code` is deterministically re-derived from the integrity-protected frozen
  `units.{length_unit,mass_unit,velocity_unit}` through
  `core::newtonGravitationalConstantCode(UnitSystem)`. It is not a second
  persisted authority and requires no restart schema field. PM and short-range
  TreePM both consume that same scale-free value. Scale factor and Hubble rate
  remain integrator truth; collisionless KDK and the gas conservative source
  apply the cosmological response to their respective state.
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
- **Supported distributed restart topology:** same-topology, rank-local restart only. A continuation must use the
  same MPI world size, same runtime rank for each rank-qualified checkpoint file, same normalized config hash,
  same PM grid shape, same PM decomposition mode, same per-rank PM slab table, and coherent TreePM cadence/field
  metadata before mutable runtime state is used. `ReferenceWorkflow` validates these fields on resume and rejects
  mismatches with expected/observed rank, world-size, slab, config, or cadence diagnostics. Rank-qualified filenames
  are a file-mapping guard, not a remapping algorithm.
- **Unsupported topology changes:** rank-count-changing restart, arbitrary particle/cell remap, and remote-row
  assignment from another rank's checkpoint are intentionally rejected. Each rank is authoritative only for its
  local restart payload under the validated topology. Imported PM/tree/hydro ghosts are transient and are rebuilt
  from owner state after restart; they are not persisted as truth.
- **Policy:** restart continuation uses `long_range_restart_policy=deterministic_rebuild` for PM mesh-field storage. The v20 force-cache addition, retained by v21, serializes the acceleration cache that the next KDK pre-kick consumes, so `ReferenceWorkflow` does not discard already-committed force state across a restart. PM mesh arrays themselves remain non-persistent and are rebuilt on the next legal refresh opportunity.
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
9. `gravity_force_cache_at_kdk_boundary`
10. `scheduler_persistent_state`
11. `gas_cell_scheduler_persistent_state_keyed_by_gas_cell_id`
12. `output_cadence_persistent_state`
13. `stochastic_module_persistent_state`
14. `restart_diagnostics_summary`
15. `distributed_gravity_state`
16. `normalized_config_text_and_hash`
17. `provenance_record`
18. `payload_integrity_hash_and_hex`
19. `amr_pending_flux_register_state`
20. `amr_temporal_boundary_history_state`

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

MPI workflow restart uses rank-qualified serial-HDF5 files and supports exact
same-world-size, same-rank, same-PM-layout continuation. It is not parallel
HDF5/MPIO, and it cannot change rank count or remap one rank's checkpoint rows
onto another rank. The np1--np4 DMO gate exercises direct-vs-resumed state under
the supported topology.

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

Restart schema v17 introduced pending AMR flux-register state; current v21 checkpoints retain it under:

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

Restart schema v19 introduced `/state/amr_temporal_boundary_history`; current `cosmosim_restart_v21` retains it for open local AMR
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

---

## Gas-cell scheduler persistence (v19)

Schema v19 introduced a separate gas-cell time-bin state; current v21 (`cosmosim_restart_v21`) persists it in
`/gas_cell_scheduler`. Its identity key is explicitly `gas_cell_id`, not
`parent_particle_id` and not a dense local row. The group contains `gas_cell_id`,
`bin_index`, `next_activation_tick`, `active_flag`, and `pending_bin_index`, together with
`current_tick` and `max_bin` attributes.

Restart validation checks that each scheduler identity agrees with the authoritative
`GasCellIdentityMap` and that `CellSoa::time_bin` is a derived mirror of `bin_index`. The
payload integrity hash and diagnostics include this state. v14--v18 files retain a narrow
compatibility reconstruction route based on persisted cell time-bin mirrors; that route never
looks up a cell through a parent particle.

A state with gas cells but no `PatchSoa` is legal for legacy/non-AMR Cartesian continuation only
when every `cells.patch_index` is the zero sentinel and identity records use
`owning_patch_id = 0`. AMR states continue to require complete patch coverage.

## H1 fixed-patch workflow restart evidence (v20)

For the H1 fixed/Cartesian-patch workflow proof, `/gravity_force_cache` stores stable
`particle_id` and `gas_cell_id` lanes beside cached acceleration components. On import,
`ReferenceWorkflow` remaps those components to the current dense rows by stable ID rather
than assuming write-time row order. A missing, duplicate, zero, or non-covering cache
identity lane is rejected. This is required because dense row is transient storage position,
not physical identity.

The HDF5 workflow-equivalence test covers a checkpoint boundary at which gas-cell rows are
intentionally reordered before resume. It proves the exercised local workflow path, not
multi-rank migration/repartition restart. Validation posture remains: CI guards;
deterministic regression and restart equivalence; multi-resolution reference convergence;
cross-code comparison; then science-readiness benchmarks.
