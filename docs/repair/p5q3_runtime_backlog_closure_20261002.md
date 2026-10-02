# P5-Q3 runtime backlog closure — 2026-10-02

## Scope

Repair the remaining source-level backlog exposed by P5-Q2 without reverting P1–P4 or M48 and without starting P5-W.

## Root-cause / repair matrix

| Family | Classification | Repair | Evidence |
|---|---|---|---|
| Periodic-seam top-domain bounds | Production defect | Build/refit compact unwrapped periodic intervals and use minimum-image coverage checks | `unit_parallel_distributed_memory` PASS |
| Uniform rung-zero append after source birth | Production defect | Governed one-time scheduler materialization only when future activation is required | `unit_retained_capacity_transaction` PASS; source-runtime tests PASS |
| Newborn drift epoch | Production defect | Stamp newly created stars at source-evaluation/step-end epoch | star-formation source-runtime PASS |
| Post-kick cache rebuild source epoch | Production defect | Pre-kick uses step-begin; post-kick uses step-end | star-formation source-runtime PASS |
| Final ownership audit after births | Production defect | Compare against evolved expected-ID authority rather than immutable IC digest | star-formation source-runtime PASS through run completion |
| Stage-6 MPI memory fixture | Stale fixture | Assert current shared communication-arena logical accounting | `unit_stage6_final_acceptance` PASS |
| H2 duplicate-ID fixture | Stale fixture | Test duplicates before issuing the generation-scoped validation certificate | `unit_species_state_organization` PASS |
| Diagnostic governor fixture | Stale fixture | Recompute current post-certificate scratch requirement before constructing one-byte-short budget | `unit_analysis_diagnostics` PASS |
| DMO/TreePM restart digest | Already fixed/current source | No production change; exercised exact uninterrupted/resumed paths | `integration_restart_equivalence_dm_only`, `integration_restart_equivalence_treepm` PASS |
| MPI/Parallel-HDF5/P3 distributed qualification | Environment blocked | No semantic workaround | `mpi-hdf5-fftw-debug` configure: missing `MPI_CXX` / `mpi-cxx` |

## Preserved architecture

- No TreePM force parameter or force mathematics changed.
- No permanent all-peer routing introduced.
- No population-scale periodic membership table introduced.
- Homogeneous DMO remains compact until heterogeneous/future-activation scheduling is actually required.
- Scheduler growth is planned through the existing `RetainedCapacityTransaction` / MemoryGovernor path.
- M48-H2 identity certification remains authoritative; tests no longer mutate compatibility storage behind a valid certificate.
- No restart digest tolerance was loosened.
- P5-W was not started.

## Validation executed

### PASS

- `cmake --preset cpu-only-debug`
- focused build of changed CPU targets
- `unit_parallel_distributed_memory`
- `unit_retained_capacity_transaction`
- `unit_stage6_final_acceptance`
- `unit_species_state_organization`
- `unit_analysis_diagnostics`
- `integration_star_formation_source_runtime`
- all seven `integration_star_formation_amr_*` / `integration_effective_ism_amr_*` aliases
- `cmake --preset hdf5-debug`
- `integration_restart_equivalence_dm_only`
- `integration_restart_equivalence_treepm`

### BLOCKED

- `cmake --preset mpi-hdf5-fftw-debug`: host lacks `mpi-cxx` / `MPI_CXX`.

### INCOMPLETE / NOT CLAIMED

A full clean CPU build was started but exceeded individual execution windows while compiling the broad test inventory; all directly affected targets were built and exercised. No full MPI/OpenMP/Parallel-HDF5 qualification is claimed.

## Remaining gates

No known source defect from the scoped backlog remains after the exercised tests. Remaining P5-Q gates are distributed qualification: P3 OFF/1/2/4 plus MPI+OpenMP incoming traversal, event-by-event P2 fallback/recovery telemetry, MPI restart artifact smoke, and Parallel-HDF5 provider-family validation on a dependency-complete host.

## Verdict

**CONDITIONALLY READY** for the next serious 64^3/128^3 DMO rehearsal, conditional on the outstanding dependency-complete MPI/Parallel-HDF5 qualification. This report does not claim 256^3/512^3 qualification.
