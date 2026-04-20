# Parallel distributed-memory contracts (repair scope)

This document captures the current **reviewable distributed-memory contract** used by the
Phase 2 TreePM implementation.
Phase 2 TreePM-specific ownership/message terminology is frozen in
`docs/treepm_phase2_distributed_contract.md` and must be used consistently with this file.

## 1) Gravity-aware decomposition contract

`buildMortonSfcDecomposition(...)` uses Morton ordering plus a weighted gravity-cost model.

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

- `packParticleMigrationRecords(local_indices)` packs authoritative hot fields and sidecar metadata for requested local rows, including species-specific sidecar payloads when the species tag requires them.
- `commitParticleMigration(...)` is the explicit ownership synchronization point:
  - outbound-owned rows are removed,
  - stale ghost/import rows are removed,
  - inbound authoritative rows are appended with `owning_rank == world_rank`,
  - species sidecars are rebuilt with remapped particle indices,
  - species counts and species index are rebuilt.
- Ownership is never committed mid-walk/mid-exchange by this contract; commit must run at a phase boundary.

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
