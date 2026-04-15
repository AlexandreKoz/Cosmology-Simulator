# Parallel distributed-memory contracts (repair scope)

This document captures the current **reviewable infrastructure contract** for distributed-memory scaffolding.
It does **not** claim full production MPI capability.

## 1) Ownership + transfer lifecycle contract (scaffold-level)

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

## 2) Migration/pack/unpack invariants

Ghost packet invariants remain:

- SoA lanes must be shape-consistent before pack/unpack.
- Packed local indices must be in range.
- `bytes_per_ghost` must be strictly positive at plan-construction time.
- Decode enforces payload-shape consistency (`encoded_count` must match exact byte payload shape).
- Decode must consume exactly the buffer payload (no trailing bytes).

Current scaffolding behavior keeps only the receive-side ghost staging indices in the descriptor-only plan. Export/send row selection requires a future owned-boundary planner.

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

## 5) What remains unproven in this environment

- No production MPI exchange correctness claim is made; this repo path is still pseudo-multi-rank + CPU-only contract scaffolding.
- Migration lifecycle beyond typed intents (actual migration transfer execution + ownership commit across ranks) is not implemented in this pass.
- Descriptor-only ghost planning does not yet compute owned export rows for outbound refresh payloads.

## 4) Config-freeze consensus contract across ranks

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
