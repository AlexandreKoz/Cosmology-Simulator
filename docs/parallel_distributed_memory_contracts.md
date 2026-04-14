# Parallel distributed-memory contracts (repair scope)

This document captures the current **reviewable infrastructure contract** for distributed-memory scaffolding.
It does **not** claim full production MPI capability.

## 1) Ownership contract: owned vs ghost local indices

`parallel::LocalGhostDescriptor` encodes local-index residency explicitly:

- `residency = kOwned`: the local index is rank-owned and `owning_rank` must equal `world_rank`.
- `residency = kGhost`: the local index is ghost-state and `owning_rank` must be a remote rank.

`buildGhostExchangePlan(world_rank, span<LocalGhostDescriptor>, bytes_per_ghost)` rejects mismatched ownership combinations and builds send/recv index maps by remote owner rank.

## 2) Migration/pack/unpack invariants

Ghost packet invariants remain:

- SoA lanes must be shape-consistent before pack/unpack.
- Packed local indices must be in range.
- Decode must consume exactly the buffer payload (no trailing bytes).

Current scaffolding behavior keeps symmetric request/response index sets per neighbor while full migration scheduling is still staged.

## 3) Deterministic reduction contract

`deterministicRankOrderedSum()` provides a deterministic reference aggregation in rank order.
`compareReductionAgreement()` reports:

- deterministic reference sum
- measured/reference absolute error
- measured/reference relative error

This enables explicit reproducibility checks in environment-independent pseudo-multi-rank tests without claiming real MPI collective determinism for all hardware/network stacks.

## 4) Config-freeze consensus contract across ranks

`evaluateRankConfigConsensus()` compares per-rank digests for:

- normalized config hash identity
- `mpi_ranks_expected` agreement
- `deterministic_reduction` agreement

Mismatch reporting is explicit by rank ID to avoid silent config drift.

## Migration notes

- Preferred API for ghost planning is now the typed descriptor overload:
  - `buildGhostExchangePlan(int, std::span<const LocalGhostDescriptor>, std::size_t)`.
- The legacy owner-rank vector overload remains, but it is an adapter that infers residency from `owner_rank == world_rank`.
- Callers that already track owned/ghost state explicitly should migrate to `LocalGhostDescriptor` to avoid implicit residency assumptions.
