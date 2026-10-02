# P5-Q2 / M48-Q runtime and validation closure report

Date: 2026-10-02  
Mode: Repair  
Campaign: P5-Q2 / M48-Q  
Verdict: **NOT READY**

## Scope and authority

The supplied source archive did not contain `.git` metadata, so `git rev-parse HEAD`, `git status --short`, and `git diff --stat` could not be recorded. The current source and repository authority documents were inspected directly. No P1-P4/M48 representation was intentionally reverted.

## Root-cause matrix

| Failure family | Classification | Root cause | Layer repaired | Evidence |
|---|---|---|---|---|
| P2 mixed authoritative/fallback geometry | PRODUCTION DEFECT | fallback was rank-local, so hierarchy exchange could receive mixed geometry frames | `src/gravity/tree_pm_coupling.cpp`, public reason enum | source compiled in CPU build; MPI runtime blocked |
| P2 tied Morton-key refit | PRODUCTION DEFECT | refit selected the first inclusive SFC interval, collapsing tied seed leaves | `src/parallel/distributed_memory.cpp` | `unit_parallel_distributed_memory` PASS |
| periodic-seam top-domain geometry | UNRESOLVED / NEEDS MORE EVIDENCE | raw min/max leaf bounds still lack a compact wrapped representation in this pass | none | remains open |
| H2 owner-query scratch | ALREADY FIXED IN CURRENT SOURCE | prepared bounded stack/owner storage already exists | none | preserved; no reimplementation |
| TreePM split-scale boundary | PRODUCTION DEFECT | strict binary `<` rejected mathematically equal values after rounding | `src/core/config.cpp` | `unit_config_parser` PASS; full CPU build PASS |
| distributed-IC routing collective count | STALE TEST/FIXTURE | communicator-global logical count was summed again across ranks | `tests/integration/test_distributed_ic_reader_mpi.cpp` | source compile; MPI runtime blocked |
| distributed PM routing fixture | STALE TEST/FIXTURE | fixture required remote density traffic on every rank | `tests/integration/test_pm_periodic_mode.cpp` | CPU compile/runtime PASS; MPI branch blocked |
| migration transport limits | STALE TEST/FIXTURE | 8-byte test ceiling is smaller than mandatory wire records/headers | three MPI integration fixtures | source compile; MPI runtime blocked |
| collective HDF5 output roots | STALE TEST/FIXTURE | fixtures supplied different rank-local run roots to collective single-file output | distributed hydro + DMO validation fixtures | source compile; MPI runtime blocked |
| DMO restart digest | UNRESOLVED / NEEDS MORE EVIDENCE | exact first divergence was not reproduced on this host | none | remains open; digest not weakened |
| MPI restart artifact | UNRESOLVED / NEEDS MORE EVIDENCE | dependency-complete runtime unavailable | none | remains open |
| HDF5 provider-family | ENVIRONMENT/DEPENDENCY | MPI C++ toolchain unavailable before provider-family qualification | none | configure BLOCKED |
| P3 current-tree MPI+OpenMP qualification | ENVIRONMENT/DEPENDENCY | MPI toolchain unavailable | none | remains open |

## Production repairs

### Collective TreePM geometry fallback

Each rank still evaluates whether authoritative geometry is safe locally. A communicator `MPI_Allreduce(MAX)` now decides whether any rank requires conservative fallback. When fallback is required, ranks that were otherwise healthy publish their derived local tree-root frame for that force event as well, and report the explicit `collective_peer_fallback` reason. This prevents mixed-frame exchange without permanently disabling authoritative geometry or sparse routing.

### Tied-SFC refit preservation

Refit no longer assigns every source whose key lies in multiple identical/overlapping seed intervals to the first seed. Matching tied intervals are selected deterministically from source ordinal without introducing persistent population-scale membership state. Existing `domain_leaf_id` values remain the seed identity. A regression verifies that two tied seed leaves survive repeated refits and retain stable IDs.

### Split-scale equality

The periodic TreePM split-scale guard now uses an 8-epsilon, scale-aware comparison around the equality boundary. This accepts floating-point-equal values while leaving materially undersized split scales illegal. No TreePM parameter, force split, softening, or tolerance was changed.

## Fixture repairs

- distributed IC logical collective counts use the already-global common value instead of an all-rank sum;
- PM density routing assertions use global traffic, allowing legitimate zero-traffic ranks;
- MPI transport-limit fixtures use 256 bytes rather than the physically impossible 8-byte ceiling while remaining small enough to force bounded fragmentation on these tiny fixtures;
- collective-output workflow fixtures use one shared run root for all ranks.

## Preserved optimizations

The repair does not restore aliased P1 counters, permanent full-peer routing, disabled OpenMP, heavyweight P4 decomposition populations, per-target routing vectors, dense homogeneous-DMO scheduler/cold sidecars, population-sized restart verification copies, or old M48 representations. M48-H2 prepared `QueryScratch` remains intact.

## Commands actually executed

### PASS

```text
bash scripts/ci/check_repo_hygiene.sh
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --test-dir build/cpu-only-debug -R '^(unit_config_parser|unit_parallel_distributed_memory)$' --output-on-failure
ctest --test-dir build/cpu-only-debug -R '^(unit_parallel_distributed_memory|integration_pm_periodic_mode)$' --output-on-failure
```

The complete CPU Debug build finished successfully. Focused tests passed.

### FAIL / current-tree failures observed during attempted full CPU inventory

The full `ctest --preset test-cpu-debug --output-on-failure` run was terminated by the command window before completion, but before termination it exposed failures including:

```text
unit_retained_capacity_transaction
unit_stage6_final_acceptance
unit_species_state_organization
integration_star_formation_source_runtime
integration_star_formation_amr_covered_coarse
integration_star_formation_amr_refine_derefine
integration_star_formation_amr_level_equivalence
integration_star_formation_amr_reflux_ordering
integration_star_formation_amr_patch_reorder
integration_effective_ism_amr_threshold_invariance
integration_effective_ism_amr_eos_restriction
unit_analysis_diagnostics
```

The star-formation/ISM failures progressed past the split-scale boundary and now fail at the separate uniform rung-zero append contract, which requires independent fixture/production classification.

### BLOCKED

```text
cmake --preset mpi-hdf5-fftw-debug
```

CMake fails because `MPI_CXX` / `mpi-cxx` is not available. Therefore no MPI-only repair is claimed as runtime-validated here, and Parallel-HDF5/FFTW qualification was not reached.

### NOT RUN / incomplete campaign gates

- periodic seam compact top-domain routing oracle tests;
- P3 OpenMP OFF/1/2/4 force-vector matrix and MPI+OpenMP incoming-target traversal;
- P2 event-by-event telemetry recovery sequence;
- DMO restart first-divergence instrumentation and exact continuation repair;
- MPI restart publication artifact qualification;
- HDF5 provider-family consumer matrix;
- production-like release performance regression capture.

## Remaining failures and blockers

The periodic seam representation and restart exact-continuation issues remain scientific/runtime blockers under this campaign. The dependency-complete distributed matrix also remains mandatory. The CPU-floor failures listed above must be classified and repaired without weakening scheduler, identity, memory-governor, or M48 contracts.

## Readiness verdict

**NOT READY**

This verdict applies specifically to the P5-Q2 / M48-Q acceptance criteria. It does not erase the previously demonstrated small DMO first-light capability; it means the current repository is not yet qualification-clean for the requested serious 64^3 / 128^3 rehearsal gate under this campaign.
