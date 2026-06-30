# AMR Hydro Subcycling and Pending Flux Registers

Date: 2026-06-17

## Scope

This document describes the post-H3 local AMR hydro subcycling and persistent pending flux-register implementation. The implementation is a credibility step for local AMR hydro. It is not a full distributed Berger-Colella AMR scheduler.

## Sweep modes

The production AMR hydro orchestrator supports two modes:

1. **Synchronized sweep**: preserves the H3 behavior. All active local AMR patches are advanced with the same hydro timestep and complete reflux records are applied immediately.
2. **Local subcycled sweep**: advances explicit AMR levels recursively. For a refinement ratio `r`, a finer level advances with `dt_fine = dt_coarse / r` and is stepped `r` times before the parent level synchronization point applies completed reflux.

The current subcycling path is local and orchestrator-scoped. The global scheduler does not yet persist AMR level timelines, and MPI-distributed subcycling remains out of scope.

## Persistent pending-register owner

Pending flux-register state lives in `core::SimulationState::pending_flux_registers`. This keeps deferred reflux state in restart-authoritative state rather than transient solver scratch.

A pending record includes stable metadata for validation:

- register key;
- coarse patch id;
- coarse gas-cell id;
- coarse patch-local cell index;
- AMR level pair metadata;
- face axis and side;
- expected/coarse/fine area coverage;
- timestep interval and coarse timestep;
- expected/completed fine substeps and coverage mask;
- coarse/fine face counts;
- gas-cell identity generation;
- patch geometry generation;
- flux-integrated mass, momentum, and total-energy contributions.

## Merge and apply policy

Hydro remains geometry-agnostic. Solver-emitted face fluxes are accumulated into AMR flux-register entries. In synchronized mode, complete entries may be applied immediately. In persistent mode, entries are merged into the pending store.

A pending register can mutate state only after it is complete and validated. Apply-time validation resolves the coarse target by stable `coarse_gas_cell_id`, checks patch ownership and patch-local mapping, verifies area and substep coverage, and rejects stale or missing targets.

Incomplete or invalid records are not silently applied.

## Restart path

HDF5 restart schema v17 introduced pending-register serialization; current schema v19 retains and serializes pending registers under `/state/amr_pending_flux_registers`. Legacy restart files without the group load with an empty pending store. Current v19 restart validation requires the group so pending deferred reflux state cannot be omitted from new checkpoints.

The restart-equivalence test `integration_restart_equivalence_amr_flux_registers` writes an incomplete pending record before restart, reloads it, completes the missing fine contribution, applies reflux, and compares final direct-vs-restart state.

## Remaining limitations

- The global production scheduler does not own subcycled AMR timelines.
- Coarse/fine ghost timing is still local and simplified.
- MPI-distributed synchronized AMR patch execution now has a bounded production path: owner ranks exchange
  explicit patch geometry plus gas-cell conserved state, build transient read-only remote patch ghosts only
  for actual adjacency, and route coarse-fine flux-register entries back to the authoritative coarse owner
  before reflux. The old all-gathered patch/cell payload summaries remain validation/debug metadata.
- MPI-distributed AMR subcycling and remote temporal ghost interpolation are not implemented. Distributed
  callers that request local subcycling or temporal coarse-to-fine interpolation fail rather than reusing
  the local temporal-history model across ranks.
- Restart after true MPI AMR migration is limited to the existing migrated patch/gas/scheduler payload and
  owner-local pending flux-register state. Rank-count-changing restart remains out of scope.
- Validation remains CI/regression scale, not cross-code science validation.

## Temporal coarse-to-fine ghost histories (v18)

The local two-level subcycling path now records a restart-authoritative temporal boundary interval in
`core::SimulationState::amr_temporal_boundary_history`. Before the coarse level advances, the
orchestrator captures coarse conserved state at `t_start_code`. After the coarse level advances through
its complete coarse step, it captures the corresponding `t_end_code` state. Fine ghost fills request the
physical start time of the particular fine update. This matches the current hydro staging: ghost values are
consumed before MUSCL-Hancock reconstruction; the predictor itself constructs face-stage states.

For `t_fill` inside `[t_start, t_end]`, the AMR layer linearly interpolates **conserved** density,
momentum, and total energy, then validates and recovers primitives through the hydro EOS conversion.
Endpoint requests use the stored endpoint exactly. Requests outside the interval, non-positive/infinite
intervals, identity generation changes, and geometry fingerprint changes are rejected rather than
silently clamped or copied from the coarse end state.

The history records stable gas-cell IDs and patch-local cell indices, plus a patch geometry fingerprint and
identity generation. Dense-row order is never a temporal-history key. Refinement and derefinement are
prohibited while an active history exists. This is a safe local lifecycle policy, not a migration/remap
implementation.

HDF5 restart schema v19 retains these records under:

```text
/state/amr_temporal_boundary_history
```

`integration_restart_equivalence_amr_temporal_ghosts` checkpoints with both an active temporal history
and an incomplete pending flux register, reloads, consumes the midpoint temporal ghost, completes the
remaining fine contribution, refluxes, and compares direct and restarted continuation.

## Current boundary

This is a **local, two-level** temporal boundary model. It does not implement remote MPI history exchange,
three-or-more nested active temporal intervals, arbitrary scheduler-owned AMR time bins, AMR patch
migration while a history is live, or temporal fine-to-coarse restriction. Fine-to-coarse use at a
non-synchronization time is rejected.
