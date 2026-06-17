# Parallel distributed-memory contracts (repair scope)

This document captures the current **reviewable distributed-memory contract** used by the
Phase 2 TreePM implementation.
Phase 2 TreePM-specific ownership/message terminology is frozen in
`docs/treepm_phase2_distributed_contract.md` and must be used consistently with this file.

## 1) Gravity-aware decomposition contract

`buildMortonSfcDecomposition(...)` uses Morton ordering plus a weighted gravity-cost model.
Rank-range cuts are contiguous in SFC order and now use nearest-target boundary selection when a
prefix crosses each rank load target. This keeps ownership deterministic while reducing clustered
load overshoot compared with a first-crossing-only split.

Per-item weighted load is:

```
weighted_load_i =
  owned_particle_weight * I(item.kind == particle)
  + active_target_weight * active_target_count_recent_i
  + remote_tree_interaction_weight * remote_tree_interactions_recent_i
  + work_weight * max(work_units_i, 0)
  + memory_weight * memory_bytes_i
```

Where:

- `active_target_count_recent` is the recent active-target pressure attributed to the entity.
- `remote_tree_interactions_recent` is recent remote tree export/import interaction pressure.
- `memory_bytes` keeps decomposition memory-aware.

`LoadBalanceMetrics` now also records per-rank decomposition components:

- `owned_particles_by_rank`
- `active_targets_by_rank`
- `remote_tree_interactions_by_rank`
- `memory_bytes_by_rank`
- `weighted_load_by_rank`

This preserves auditable decomposition inputs while still using one deterministic SFC ordering.

## 2) Ownership + transfer lifecycle contract (scaffold-level)

`parallel::LocalGhostDescriptor` encodes local-index residency explicitly:

- `residency = kOwned`: the local index is rank-owned and `owning_rank` must equal `world_rank`.
- `residency = kGhost`: the local index is ghost-state and `owning_rank` must be a remote rank.

`buildGhostExchangePlan(world_rank, span<LocalGhostDescriptor>, bytes_per_ghost)` rejects mismatched ownership combinations and builds send/recv index maps by remote owner rank.

`parallel::GhostTransferDescriptor` now makes lifecycle intent and post-transfer expectation explicit:

- `role = kOutboundSend`: send-side payload descriptor for a neighbor rank.
- `role = kInboundReceive`: receive-side payload descriptor for a neighbor rank.
- `intent`:
  - `kGhostRefreshRequest` (outbound): request side of ghost refresh.
  - `kGhostRefreshReceiveStaging` (inbound): receive staging side of ghost refresh response.
  - `kOwnershipMigrationSend` / `kOwnershipMigrationReceiveStaging`: reserved typed intents for future migration wiring (not produced by current scaffold builder).
- `peer_rank`: remote rank for this payload contract.
- `neighbor_slot`: canonical slot in `neighbor_ranks`.
- `expected_post_transfer_residency`: explicit post-transfer local residency expectation (`kGhost` for current ghost-refresh scaffold).
- `local_indices`: local indices participating in this transfer payload.

`GhostExchangePlan` now publishes both legacy index vectors and typed transfer-role lists:

- `outbound_transfers` (`kOutboundSend` only)
- `inbound_transfers` (`kInboundReceive` only)

`validateGhostExchangePlan()` enforces role/intent/container consistency, peer-rank + slot consistency, and exact descriptor-index equality with canonical per-neighbor vectors.

In this repair scope, the builder only emits ghost-import slots and ghost-refresh intents. `recv_local_indices_by_neighbor` identifies local ghost rows that expect refreshed data, while `send_local_indices_by_neighbor` intentionally remains empty for the descriptor-only planner because peer-owned export rows are not derivable from local ghost placeholders alone. Full migration scheduling and migration ownership commit remain out of scope.

## 3) Migration/pack/unpack invariants

Ghost packet invariants remain:

- SoA lanes must be shape-consistent before pack/unpack.
- Packed local indices must be in range.
- `bytes_per_ghost` must be strictly positive at plan-construction time.
- Decode enforces payload-shape consistency (`encoded_count` must match exact byte payload shape).
- Decode must consume exactly the buffer payload (no trailing bytes).

Current scaffolding behavior keeps only the receive-side ghost staging indices in the descriptor-only plan. Export/send row selection requires a future owned-boundary planner.

For particle ownership migration in `core::SimulationState`:

- `ParticleMigrationRecord` is a schema-like identity-preserving transform contract. A record carries the persistent hot lanes (`position_*`, `velocity_*`, `mass_code`, `time_bin`), common cold metadata (`particle_id`, `sfc_key`, `species_tag`, `particle_flags`, `owning_rank`), materialized gravity-softening values plus the authoritative override mask, and the complete star, black-hole, or tracer sidecar payload required by the record species.
- `packParticleMigrationRecords(local_indices)` packs authoritative hot fields and sidecar metadata for requested local rows, including species-specific sidecar payloads when the species tag requires them. Missing required sidecar rows are hard errors, because silently migrating only hot lanes would corrupt `{particle_id -> sidecar payload}` identity.
- `commitParticleMigration(...)` is the explicit ownership synchronization point:
  - outbound-owned rows are removed,
  - stale ghost/import rows are removed only when their local `owning_rank` is remote from `world_rank`; duplicate stale-ghost indices are rejected,
  - inbound authoritative rows are appended with `owning_rank == world_rank`,
  - final local particle IDs are checked for uniqueness across kept plus inbound rows,
  - species sidecars are rebuilt from the post-commit `species_tag` ledger with remapped particle indices,
  - species counts and `ParticleSpeciesIndex` are rebuilt from `species_tag`,
  - the particle-index generation counter is bumped once for the index-space mutation.
- Gas-cell state is not keyed by old local particle positions during ownership compaction. Particle migration may carry one compatibility gas-cell payload only when a migrating gas particle has a single attached gas-cell row. Parentless cells and split/merged cells with non-unique parent lineage must move through `GasCellMigrationRecord`, keyed by stable nonzero `gas_cell_id`.
- `GasCellMigrationRecord` carries the gas identity record, parent-lineage presence bit, owning patch identity, destination row hint, hydro conserved/persistent sidecar fields, cell timestep mirror, and ghost generation/epoch metadata. `commitGasCellMigration(...)` removes outbound/stale local rows, rejects stale ghosts whose `{gas_cell_id, gas_cell_identity_generation, ghost_hydro_epoch}` do not match the commit boundary, appends inbound authoritative rows, rebuilds `SimulationState::gas_cell_identity`, remaps host-cell sidecars by old row to new `gas_cell_id`, and bumps the cell-index generation once.
- `AmrPatchMigrationRecord` is the AMR ownership-transfer payload for patch-aware hydro. It carries the patch descriptor (`patch_id`, owner rank, level, parent patch, Morton key, origin, extent, cell dimensions) plus every authoritative `GasCellMigrationRecord` in that patch. `commitAmrPatchMigration(...)` removes outbound patch descriptors and their gas-cell rows in the same commit boundary, appends inbound patch descriptors and gas-cell sidecars atomically, rewrites local `cells.patch_index` to the rebuilt dense patch table, rejects stale gas ghosts by generation/epoch, rebuilds `SimulationState::gas_cell_identity`, and bumps the cell-index generation once. A parent particle referenced by a migrated gas cell is lineage metadata unless that particle also migrates through `ParticleMigrationRecord`; patch ownership is the gas-cell authority, so no authoritative gas cell may remain on the old rank after an outbound patch commit.
- Imported gas ghosts are read-only boundary state. They may be discarded by stale-ghost records or restored after solver use, but they must not mutate authoritative hydro truth unless routed through an explicit owner-side conservative correction path.
- Ownership is never committed mid-walk/mid-exchange by this contract; commit must run at a phase boundary so stale active/kernel views are invalidated before further scatter.

## 3) Deterministic reduction contract

`deterministicRankOrderedSum()` provides a deterministic reference aggregation in rank order.
`compareReductionAgreement()` reports:

- deterministic baseline sum (`deterministic_baseline_sum`)
- measured reduction under test (`measured_sum`)
- baseline-vs-measured absolute error
- baseline-vs-measured relative error

`ReductionAgreementPolicy` now includes explicit `mode`:

- `kAbsoluteOnly`
- `kRelativeOnly`
- `kAbsoluteAndRelative`
- `kAbsoluteOrRelative`

`satisfiesReductionAgreement()` uses this explicit mode instead of implicit boolean logic. The pseudo-multi-rank infrastructure tests currently default to `kAbsoluteOrRelative` to preserve scaffold tolerance behavior while making the policy auditable.

This enables explicit reproducibility checks in environment-independent pseudo-multi-rank tests without claiming real MPI collective determinism for all hardware/network stacks.

## 6) Remaining scope limits (post-Phase 2 closeout)

- PM decomposition mode supports:
  - `slab`: baseline x-slab distributed ownership,
  - `pencil`: transposed FFT mode with explicit real-space x-slab ownership and intermediate y-owned transposed spectral ownership.
- Restart continuation of PM long-range field values is policy-limited to `deterministic_rebuild`
  (metadata is serialized; cached field arrays are rebuilt after resume).
- Descriptor-only ghost-plan builder still does not infer owned export rows from ghost-only local
  metadata; explicit owned-boundary export planning remains future work.

## 5) Config-freeze consensus contract across ranks

`evaluateRankConfigConsensus()` compares per-rank digests for:

- normalized config hash identity
- `mpi_ranks_expected` agreement
- `deterministic_reduction` agreement

Mismatch reporting is explicit by rank ID to avoid silent config drift.
Consensus now also records machine-readable mismatch rows (`property`, `baseline_rank`, `rank`, `baseline_value`, `rank_value`) for CI/debug artifacts.

## Migration notes

- Preferred API for ghost planning is now the typed descriptor overload:
  - `buildGhostExchangePlan(int, std::span<const LocalGhostDescriptor>, std::size_t)`.
- The legacy owner-rank vector overload remains, but it is an adapter that infers residency from `owner_rank == world_rank`.
- Callers that already track owned/ghost state explicitly should migrate to `LocalGhostDescriptor` to avoid implicit residency assumptions.
- Transfer-payload role auditing should prefer `GhostExchangePlan::outbound_transfers` / `inbound_transfers` over ad-hoc interpretation of raw send/recv index vectors.

## Phase 2 distributed gravity vocabulary alignment

For upcoming distributed TreePM work:

- `particle owner`: rank that owns and updates a particle's authoritative state.
- `slab owner`: rank that owns a PM x-slab when `treepm_pm_decomposition_mode=slab`.
- `distributed execution topology`: the explicit runtime tuple `(world_size, world_rank, pm_slab, device_assignment)` derived from MPI world state plus validated GPU request state.
- `tree export/import`: bounded payload exchange for short-range tree source support, capped by `treepm_tree_exchange_batch_bytes`.

These names are contractual. They prevent pseudo-distributed ambiguity where all ranks hold replicated state but are described as distributed.

## PM slab layout ownership contract (new explicit type)

`parallel::PmSlabLayout` defines the auditable PM real-space ownership mapping used by gravity storage:

- decomposition axis: x only,
- slab range representation: half-open `[begin_x, end_x)`,
- per-rank slab rule for `global_nx` and `world_size`:
  - `base = global_nx / world_size`,
  - `remainder = global_nx % world_size`,
  - rank `r` owns `base + 1` x-indices if `r < remainder`, else `base`,
  - `begin_x = r * base + min(r, remainder)`.

Contracted helpers:

- `pmOwnedXRangeForRank(...)`: deterministic per-rank slab bounds,
- `pmOwnerRankForGlobalX(...)` / `pmOwnerRankForGlobalCell(...)`: ownership lookup,
- `PmSlabLayout::{localXFromGlobal, globalXFromLocal, localLinearIndex}`:
  centralized global/local conversions with range validation.

Single-rank mapping is the same abstraction, not a special case:

- `world_size = 1`, `world_rank = 0`, owned slab `[0, global_nx)`, and local storage equals full mesh.

This contract is intentionally storage/ownership-centric. Distributed PM communication covers both
owner-routed deposition and reverse interpolation delivery:

- particle owners compute global stencil contributions,
- destination slab owner is resolved by `pmOwnerRankForGlobalX(...)`,
- contributions are batched per destination and exchanged with collective messaging,
- only slab owners accumulate into local PM density arrays after ownership/range validation,
- slab owners compute interpolation contributions for owned stencil cells and return force/potential
  contributions back to particle owners via reverse communication.

## Distributed TreePM restart/debug continuation contract

`parallel::DistributedRestartState` (schema_version `2`) now defines the serialized distributed
continuation contract used by restart checkpoints:

- decomposition ownership: `decomposition_epoch`, `world_size`, `owning_rank_by_item`
- PM layout contract: `pm_grid_nx/ny/nz`, `pm_decomposition_mode`, slab ownership vectors
  `pm_slab_begin_x_by_rank` / `pm_slab_end_x_by_rank`
- cadence/field metadata: `gravity_kick_opportunity`, `pm_update_cadence_steps`,
  `long_range_field_version`, `last_long_range_refresh_opportunity`,
  `long_range_field_built_step_index`, `long_range_field_built_scale_factor`
- long-range continuation policy: `long_range_restart_policy`

Current policy is explicitly fixed to `deterministic_rebuild`:

- cached PM long-range arrays are not serialized;
- restart resumes with deterministic rebuild at the next cadence-triggered refresh;
- this policy is versioned and validated during deserialize + restart hashing.

`evaluateDistributedRestartCompatibility(...)` compares restart metadata against current
`DistributedExecutionTopology` and returns typed booleans + mismatch messages for:

- world size mismatch,
- PM grid-shape mismatch,
- PM decomposition-mode mismatch,
- per-rank local slab ownership mismatch,
- invalid PM cadence metadata (`pm_update_cadence_steps >= 1` required).

This reporting is used for continuation debugging and avoids opaque restart failures when
rank layout/config drifts between write and resume.

## Distributed ownership identity reductions (P1 correctness floor)

Distributed workflow tests must prove partitioned ownership rather than accepting replicated state.
Each rank now emits a local ownership identity summary over the rank-local authoritative particle-ID
lane:

- `local_owned_count`: number of locally owned particle rows after compaction/restart.
- `local_particle_id_sum`: unsigned sum of local authoritative particle IDs.
- `local_particle_id_xor`: bitwise XOR of local authoritative particle IDs.
- `local_particle_ids_unique`: local duplicate-ID rejection signal.

MPI tests reduce these summaries with `SUM` for counts and ID sums and `BXOR` for ID XORs, then
compare the reduced tuple with the generated initial-condition identity tuple. A two-rank test is
not considered distributed-correct if all ranks hold the full generated data set: replicated ranks
would over-count `local_owned_count` and `local_particle_id_sum`, and even-count replication would
zero or otherwise alter the XOR. The local duplicate flag is reduced across ranks before report-level
identity status is marked green.

These reductions are correctness checks only. They do not introduce mature load balancing, pencil FFT,
or performance-tuned migration scheduling. They harden identity and ownership determinism for the
current slab/TreePM path while preserving the long-range restart policy of deterministic rebuild.

## Reproducibility impact of this repair

The repair does not alter solver numerics, force kernels, PM decomposition policy, or restart payload
schema version. It adds auditable report fields and reduced identity checks around the existing
rank-local compacted state:

- restart long-range PM continuation remains `deterministic_rebuild`;
- restart metadata still carries `owning_rank_by_item` and PM slab ownership tables;
- rank-local particle IDs are summarized deterministically using integer count/sum/xor reductions;
- partition-vs-replication detection is a test/report invariant, not a scheduling or load-balancing
  behavior change.

## Ghost refresh payloads are not ownership migration records

`GhostExchangeBuffer` carries the narrow ghost-refresh payload shape only: entity ID plus the small hydro-style ghost fields used by the current exchange path. It is not an ownership-transfer schema. Descriptor-aware `packFrom()` and `unpackAppendTo()` calls therefore reject `kOwnershipMigrationSend` and `kOwnershipMigrationReceiveStaging` intents.

Persistent ownership migration must continue to use `core::ParticleMigrationRecord`, because that record carries authoritative runtime-truth lanes: particle hot state, particle identity metadata, softening value/mask, species tag, owning rank, and species-specific sidecar payloads. Any future distributed ownership packet must either be exactly this schema or a versioned superset with tests proving sidecar, softening, restart, and gas-cell identity preservation.

Partition correctness checks now use count, particle-ID sum, particle-ID square sum, particle-ID XOR, and local uniqueness. The reference workflow compares reduced rank-local identity against the generated pre-partition IC identity. Passing by comparing a rank-reduced partition summary to itself is forbidden because it cannot detect replicated-state success.
