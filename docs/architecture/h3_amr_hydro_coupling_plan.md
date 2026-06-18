# H3 AMR Hydro Coupling Plan and Boundaries

Date: 2026-06-17

## Current production coupling

The production AMR hydro coupling uses `core::SimulationState` as the authoritative owner of gas-cell identity, hydro primitive/conserved state, patch ownership, and restart truth. AMR patch descriptors provide explicit geometry lanes. Hydro patch geometry is built from `PatchDescriptor` plus stable `GasCellIdentityMap` rows, not from sorted dense row order.

When production AMR coverage exists, the workflow uses the AMR hydro orchestrator path. The old single global Cartesian hydro patch remains a fallback only for non-AMR coverage.

## Restart coupling

AMR patch geometry and gas identity are restart truth. The HDF5 restart path persists PatchSoa geometry lanes and `/state/gas_cell_identity`. The dedicated `integration_restart_equivalence_amr_hydro` test proves local direct vs restart continuation equivalence for an exercised synchronized-sweep AMR hydro case.

Legacy restart states without explicit AMR geometry must not silently enter the production AMR hydro path as if they were complete AMR states.

## Synchronized sweeps, not subcycling

The current AMR hydro execution model is synchronized local sweeps. It does not implement:

- recursive fine-level time subcycling;
- Berger-Colella level scheduling;
- level-local dt authority in the production scheduler;
- temporal interpolation of coarse/fine ghosts;
- deferred multi-step reflux replay.

This is an intentional H3 boundary. Subcycling requires a separate design with scheduler, ghost timing, pending register, restart, and MPI tests.

## Flux-register boundary

Flux registers are generated from real hydro face fluxes and reflux is applied only when records are complete and area-consistent. Incomplete records are skipped/rejected and counted. They are not stored as pending persistent state.

This avoids corrupting conserved state, but it is not a full deferred-register model. A future implementation must persist pending records with patch identity, gas-cell identity, expected patch-local target, area/dt coverage, and generation/epoch metadata, then prove restart replay equivalence.

## MPI/migration boundary

Current H3.8 evidence is local payload/commit contract coverage. It verifies patch geometry lanes, gas sidecars, identity rebuild, row-order changes, stale ghost epoch rejection, and local production AMR coverage after commit.

It does not prove real MPI AMR migration. A future MPI acceptance test must run with at least two ranks and verify patch/gas transfer, receiver identity coverage, stale remote ghost rejection across ranks, local AMR hydro or ghost-fill after receive, and restart after migration if that claim is made.

## Validation boundary

The current AMR shock tube, Sedov, and synchronization stress tests are CI guards. They are useful for regression protection but must not be described as convergence studies or cross-code validation. Scientific validation remains future work.
