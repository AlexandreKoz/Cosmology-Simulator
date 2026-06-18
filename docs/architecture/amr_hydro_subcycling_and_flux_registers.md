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

HDF5 restart schema v17 serializes pending registers under `/state/amr_pending_flux_registers`. Legacy restart files without the group load with an empty pending store. Current v17 restart validation requires the group so pending deferred reflux state cannot be omitted from new checkpoints.

The restart-equivalence test `integration_restart_equivalence_amr_flux_registers` writes an incomplete pending record before restart, reloads it, completes the missing fine contribution, applies reflux, and compares final direct-vs-restart state.

## Remaining limitations

- The global production scheduler does not own subcycled AMR timelines.
- Coarse/fine ghost timing is still local and simplified.
- MPI-distributed AMR subcycling is not implemented.
- Restart after true MPI AMR migration is not implemented.
- Validation remains CI/regression scale, not cross-code science validation.
