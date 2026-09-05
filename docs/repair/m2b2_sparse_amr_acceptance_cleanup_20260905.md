# M2B-2 sparse AMR acceptance cleanup

Date: 2026-09-05

## Scope

This focused cleanup closes the two residual findings from the post-M2B final
audit without redesigning the production AMR solver.

## Changes

- `PatchHierarchy` now retains only the gas-cell IDs retired by the most recent
  derefinement transaction. The retired-ID vector is transaction-local
  invalidation evidence rather than an append-only historical ledger.
- The 12-cycle scaffold regression now asserts that retired-ID cardinality stays
  fixed at one child octet's cells and that each cycle replaces, rather than
  accumulates, the prior IDs.
- A registered MPI+HDF5 integration test now sends an AMR patch from rank 0 to
  rank 1 through `exchangeBoundedAlltoallBytes` and the authoritative AMR
  migration wire. Pending reflux and temporal-boundary history migrate with the
  patch, are remapped by stable identity on destination commit, survive HDF5
  checkpoint/reload, and remain usable by a post-restart migration commit.
- `CURRENT_STATUS.md` and `docs/validation_plan.md` now describe the current
  validation surface without promoting unexecuted MPI evidence.

## Scientific impact

No refinement criterion, prolongation/restriction rule, reflux equation, CFL
rule, reconstruction, Riemann solver, precision, or scientific tolerance is
changed. The patch changes scaffold memory lifetime and validation coverage only.

## Validation contract

CPU validation must at minimum rebuild and run `unit_amr_refinement` plus the
existing local AMR migration/wire regressions. The MPI test is authoritative only
when `mpi-hdf5-fftw-debug` configures and the registered two-rank test executes;
an unavailable MPI dependency is a blocker, not a pass.
