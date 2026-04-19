# TreePM Phase 2 distributed contract (frozen scaffold)

## Scope and non-goals

This document freezes what **Phase 2** means in this repository for distributed gravity infrastructure.

- In scope: ownership semantics, decomposition terminology, message contract, config/documentation/build surface.
- Out of scope in this stage: introducing new physics behavior beyond the existing TreePM split contract.

Phase 1 remains the numerical baseline. A one-rank run through future distributed code paths must stay numerically consistent with the current Phase 1 single-rank TreePM contract unless an explicit bug-fix ADR states otherwise.

## Phase 2 ownership model

### 1) PM mesh ownership (target contract)

- PM mesh ownership is distributed by **x-slabs**.
- For a global mesh with `N = numerics.treepm_pm_grid` and `R = world_size`, slab owner rank `r` owns a contiguous x-index interval `[x_begin_r, x_end_r)`.
- Slab decomposition mode is configured by `numerics.treepm_pm_decomposition_mode` and currently only allows `slab`.
- Empty slabs are allowed when `R > N`; ownership remains deterministic from rank order.

### 2) Particle owner vs slab owner

- A particle has exactly one **particle-owner rank** at any time.
- PM density assignment and force interpolation require transient communication between particle owners and slab owners.
- Particle ownership is not redefined by PM slab ownership.

### 3) Distributed FFT path (targeted sequence)

For each PM refresh event in Phase 2 runtime:

1. Particle owners export assignment contributions to slab owners.
2. Slab owners assemble local density slabs.
3. Slab owners perform distributed FFT/Poisson/gradient/inverse FFT with periodic zero-mode policy unchanged from Phase 1:
   - `∇²φ = 4π G a² (ρ - ρ̄)`
   - `φ_{k=0}=0`, `a_{k=0}=0`
4. Slab owners send needed force samples/interpolants back to particle owners for active particles.

This sequence is implemented in the current runtime path; the document remains the contract anchor for ownership and ordering.

## Tree export/import contract (implemented)

- Short-range tree work is particle-owner-centric.
- Remote source data are exchanged in explicit batches, bounded by:
  - `numerics.treepm_tree_exchange_batch_bytes` (must be `> 0`)
- Exchange contract:
  - target-owner rank exports compact active-target request batches to each peer, capped by the configured byte limit,
  - peer ranks evaluate each request target against local tree/source data with the current short-range kernel, MAC, and cutoff semantics,
  - peers return partial acceleration responses keyed by batch/request IDs,
  - target-owner rank validates per-peer request/response coverage (duplicate/missing detection) and accumulates returned partials.
- Export/import does not imply ownership transfer.

## Active-set and migration timing

- Active-set restriction remains authoritative: force accumulation only mutates active slots on owner rank.
- Ownership migration is an explicit phase boundary operation (not mid-kick).
- PM refresh cadence semantics from Phase 1 remain unchanged (`treepm_update_cadence_steps`).

## Restart continuity expectations

- No snapshot/restart schema drift is introduced by this contract freeze.
- Continuation artifacts must preserve normalized config/provenance consistency across ranks.
- Future distributed restart extensions require explicit schema versioning, compatibility behavior, docs updates, and tests in the same patch.

## Determinism limits (honest scope)

- Deterministic behavior is required where practical (normalized config, ownership maps, rank-ordered control flow).
- Bitwise identity across heterogeneous MPI networks/collectives is **not** promised by this document.
- Pseudo-multi-rank tests do not prove completion of distributed TreePM physics or communication correctness.

## Why slab first

`slab` is frozen as the first PM decomposition mode because it is the smallest auditable extension from the current single-rank cubic-mesh contract:

- straightforward contiguous ownership mapping,
- clear owner responsibilities per x-interval,
- minimal moving parts while preserving a path to later pencil decomposition and accelerator-aware remaps.
