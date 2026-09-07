# CHUÍ / CosmoSim — M2D-1 Memory Scheduling Closure Repair

Date: 2026-09-06. Mode: focused repair. Base: `Cosmology-Simulator-main(14)(7).zip`.
Base SHA-256: `5257af8039cf1cccac94085e35633afa493595bb586f2ade9de8dab9cc3f5dc5`.

## Verdict

The reproduced P1 partition defects and the identified batch/cadence safety gaps are repaired in source, with focused CPU evidence. The dispatcher-owned reservation facility is implemented and tested, but not every built-in major task has a complete independently certified peak estimator. Full M2D acceptance remains **provisional**, not closed: real MPI+FFTW execution, all owner-specific peak contracts, and broad scientific/topology gates remain outstanding. No numerical method or precision was changed. Do not infer full-physics workstation feasibility from the synthetic byte tests.

## Implemented changes and ownership

### D1 — exact SFC feasibility before work balance

The exact planner groups equal `(Morton key, entity ID)` points, checks aggregate group bytes, constructs right-greedy mandatory suffix boundaries, then cuts for work only when the suffix remains feasible. It preserves contiguous SFC ordering, indivisible groups, optional empty ranks, and the configured transient reserve. Work imbalance is secondary to the hard memory ceiling. The original 100:1:1:1 work counterexample now produces two 8-byte ranks under an 8-byte cap. The existing 260-to-140-byte peak regression remains valid.

### D2 — distributed memory-safe cut repair and progress

Sampled weighted cuts are proposals. When their exact global target memory exceeds the per-rank allowance, the production MPI path uses bounded local prefix blocks, checked two-limb uint64 reductions, fixed-width key/ID searches, and one candidate per rank to refine the cuts. No global entity population is gathered. The same exact feasibility-first helper is exercised by CPU synthetic tests; its MPI callbacks are not claimed to have executed here. All O(P) helper metadata is prepared through collective failure agreement before refinement collectives. Existing bounded migration, collective failure coordination, and transaction reservations remain authoritative.

An existing hard-memory violation bypasses the soft migration-fraction veto. A target that remains unsafe is rejected, and the distributed planner now fails explicitly if allowed migrations cannot provide relief instead of silently reporting an ordinary soft-throttle decision. This does not promise that every feasible final ownership can fit the transient old/new coexistence budget: that transaction remains subject to the governor and may fail safely.

### A/B — explicit task admission ownership

`RuntimeTaskSchedulingProfile` now distinguishes a complete peak from an unknown estimate and dispatcher-owned from owner-managed memory. A dispatcher-owned task acquires one governor reservation for the task callback lifetime and releases it on normal or exceptional exit. Owner-managed tasks may provide a state-dependent preflight callback; its temporary reservation is released before the owner's existing physical lease, preventing duplicate charges. The admission gate uses the existing `FailureCoordinator` before entering task collectives. Unknown peaks, declared dependencies, and non-task-end lifetimes disallow speculative overlap. The registry remains the task-order authority and the production scientific dispatcher remains serial.

The remaining major-task owner-specific estimates are not all complete: gravity, sources, analysis, and output still rely on their existing local admissions and some static graph peaks remain unknown. No unverified constant has been presented as a certified bound. Before enabling concurrent production execution, publish and test complete incremental models from these owners, cover persistent/replacement coexistence, and validate each proposed overlap against actual live reservations and scientific ordering. No new executor, configuration namespace, or memory authority was introduced.

### C — retained-capacity-aware deterministic batching

The central selector accepts multiple workspace widths and already-accounted retained capacities. A workspace that fits its retained capacity costs zero new allocation; growth charges the complete replacement buffer, so old/new coexistence is preserved. The source runtime uses this for the temporary stellar budget, retained feedback-event buffer, and retained contiguous star indices. Full spatial-index rebuild staging remains an explicit fixed reserve. Checked aggregate arithmetic rejects overflow even in unlimited mode. The selected size never replaces actual governor reservations.

### O — explicit coalesced optional cadence

Light/heavy science each retain constant-size first/latest due epoch and checked missed count. Repeated missed products coalesce into one current-state diagnostic. The actual execution epoch is never relabeled as the original due epoch. Events record first/latest due, missed/coalesced counts, actual execution step, and no historical replay. Pending work is cleared only after a successful write. At run-segment termination, remaining pending work is reported as dropped; it is never forced through a hard memory rejection.

Optional cadence is deliberately best-effort and nonpersistent. Existing checkpoint provenance text records the pending summary and policy, but the restart schema is unchanged. On restart a fresh optional cadence begins and the previous checkpoint's recorded summary is emitted as a policy event. Exact historical science products are not guaranteed across persistent pressure, segment boundaries, or restart. Users requiring exact cadence must use an explicitly budgeted exact-output mechanism; this repair does not silently fabricate historical physical states. Required run-health diagnostics and the scientific integrator remain unchanged.

## Quantitative and scientific evidence

- Exact serial SFC: original 4 × 4 B, weights 100:1:1:1, np2/cap8 → [8,8]. Additional cases cover transient reserve, exact capacity, empty ranks, zero-byte items, and indivisible duplicate-key/ID groups.
- Deterministic shared cut-repair helper: 86,016 exhaustive cases (4^6 byte sequences, 2–4 ranks, caps1–7) compared with an independent contiguous-partition feasibility oracle. This is CPU algorithm evidence, not MPI execution.
- Retained workspace test: 8 retained slots, 8+16+4 B widths, 80 B fixed staging, 1,000 B ceiling, 800 B baseline → batch8, incremental64 B, accounted944 B with fixed staging. At baseline500 → batch15, incremental420 B. At baseline920 → no feasible minimum batch. Unlimited and aggregate-overflow cases are covered.
- Runtime task tests verify owner-managed preflight, dispatcher-owned commitment/lifetime, failed admission cleanup, unknown-peak/dependency overlap rejection.
- Reference workflow tests verify transient catch-up, two missed cadence epochs coalescing into a current-state product, and terminal-pressure dropped-product reporting.
- No force criterion, PM assignment, reconstruction, Riemann solver, CFL, refinement tolerance, feedback equation, precision, or physical restart state was changed.

## Commands and outcomes

Fresh base configuration:

```bash
cmake --preset cpu-only-debug
cmake --preset mpi-hdf5-fftw-debug
cmake --preset pm-hdf5-fftw-debug
cmake --preset hdf5-debug
```

CPU and HDF5 configuration PASS (GCC14.2, CMake3.31.6, OpenMP4.5, HDF5 1.14.5). MPI preset blocked at configuration: `Could NOT find MPI (missing: MPI_CXX_FOUND CXX)`; `mpi-cxx` unavailable. PM preset blocked: `COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library was not found.` These were attempted once, not repeatedly retried.

Focused CPU build and tests:

```bash
cmake --build build/cpu-only-debug --target test_unit_memory_governor test_unit_parallel_distributed_memory test_unit_runtime_module_registry test_integration_reference_workflow -j4
ctest --test-dir build/cpu-only-debug -R '^(unit_memory_governor|unit_parallel_distributed_memory|unit_runtime_module_registry|integration_reference_workflow)$' --output-on-failure
```

All four targets compiled and linked; the final changed-source runs of all four tests pass, including the expanded end-to-end and exhaustive SFC regressions. An intermediate test initially failed because a fixture expected batch13 rather than the correct batch15, and a terminal-cadence fixture used four additional steps rather than one; both test expectations were corrected, and the final targeted runs pass. A compile error from passing a governor object instead of its snapshot in the new overflow test was corrected. These were test-harness corrections, not production changes to hide a failure.

The original post-merge audit's 136/136 CPU result is historical. This repair does not claim that full CPU/HDF5 inventories, source-package-completeness, real MPI np2/np3/np4, full TreePM scientific comparisons, or topology benchmarks have been freshly completed. See the final handoff response for any additional validation completed after this source freeze.

## Remaining acceptance gates / Campaign E handoff

1. Run the dependency-complete MPI+HDF5+FFTW matrix including skewed/empty/infeasible ownership, hard-memory relief with low soft migration fraction, np2/np3/np4 agreement, collective fault injection, and actual migration transaction feasibility.
2. Complete owner-specific major-task incremental peak/lifetime contracts before enabling nonserial overlap; retain the existing subsystem reservations and governor authority.
3. Complete full CPU/HDF5 and source-package inventories, then compare scientific results across legal batch sizes and supported scheduling orders. Collect actual rank-max/mean RSS and work imbalance, and Release topology measurements on representative hardware.
4. Keep Campaign E independent work possible, but do not treat full-physics memory scheduling or distributed production as certified until these gates pass. No AIMD, GPU scheduling, mixed precision, shared-memory rewrite, or new physics is part of this repair.

Suggested branch: `campaign-m2d1-memory-scheduling-closure`.
Suggested PR title: `Close memory-constrained SFC scheduling and runtime admission gaps`.

## Final focused HDF5 evidence

After the source freeze, the following HDF5-linked targets compiled and linked successfully:

```bash
cmake --build build/hdf5-debug --target test_unit_memory_governor test_unit_parallel_distributed_memory test_unit_runtime_module_registry test_integration_reference_workflow -j4
ctest --test-dir build/hdf5-debug -R '^(unit_memory_governor|unit_parallel_distributed_memory|unit_runtime_module_registry|integration_reference_workflow)$' --output-on-failure
```

PASS 4/4, 0 failures, 25.13 s. The initial aggregate HDF5 compilation was interrupted after 85/107 steps by the external command limit; incremental continuation completed the remaining 22 steps without source errors. This is focused HDF5 evidence, not the full HDF5 inventory.
