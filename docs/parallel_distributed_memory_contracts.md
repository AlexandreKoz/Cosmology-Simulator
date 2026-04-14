# Parallel distributed-memory contracts (repair scope)

This document captures the current **reviewable infrastructure contract** for distributed-memory scaffolding.
It does **not** claim full production MPI capability.

## 1) Ownership + transfer-role contract: owned/ghost/send/receive

`parallel::LocalGhostDescriptor` encodes local-index residency explicitly:

- `residency = kOwned`: the local index is rank-owned and `owning_rank` must equal `world_rank`.
- `residency = kGhost`: the local index is ghost-state and `owning_rank` must be a remote rank.

`buildGhostExchangePlan(world_rank, span<LocalGhostDescriptor>, bytes_per_ghost)` rejects mismatched ownership combinations and builds send/recv index maps by remote owner rank.

`parallel::GhostTransferDescriptor` now makes transfer-payload role explicit:

- `role = kOutboundSend`: send-side payload descriptor for a neighbor rank.
- `role = kInboundReceive`: receive-side payload descriptor for a neighbor rank.
- `peer_rank`: remote rank for this payload contract.
- `local_indices`: local indices participating in this transfer payload.

`GhostExchangePlan` now publishes both legacy index vectors and typed transfer-role lists:

- `outbound_transfers` (`kOutboundSend` only)
- `inbound_transfers` (`kInboundReceive` only)

In this repair scope, request/response sets remain symmetric scaffolding; no full migration scheduler is claimed.

## 2) Migration/pack/unpack invariants

Ghost packet invariants remain:

- SoA lanes must be shape-consistent before pack/unpack.
- Packed local indices must be in range.
- `bytes_per_ghost` must be strictly positive at plan-construction time.
- Decode enforces payload-shape consistency (`encoded_count` must match exact byte payload shape).
- Decode must consume exactly the buffer payload (no trailing bytes).

Current scaffolding behavior keeps symmetric request/response index sets per neighbor while full migration scheduling is still staged.

## 3) Deterministic reduction contract

`deterministicRankOrderedSum()` provides a deterministic reference aggregation in rank order.
`compareReductionAgreement()` reports:

- deterministic baseline sum (`deterministic_baseline_sum`)
- measured reduction under test (`measured_sum`)
- baseline-vs-measured absolute error
- baseline-vs-measured relative error

`satisfiesReductionAgreement()` provides a narrow policy gate over absolute/relative tolerances so callers/tests do not reinterpret raw error fields ad hoc.

This enables explicit reproducibility checks in environment-independent pseudo-multi-rank tests without claiming real MPI collective determinism for all hardware/network stacks.

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
