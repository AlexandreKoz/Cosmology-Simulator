# M2D full-physics memory-aware runtime scheduling — 2026-09-06

## Scope and authority

Campaign M2D starts from the accepted M2A–M2C memory architecture and extends
its existing authorities. It does **not** introduce a second configuration
surface, allocator, memory reporter, timestep scheduler, or decomposition
geometry. The process `core::MemoryGovernor` remains the sole runtime memory
admission authority; `RuntimeModuleRegistry` remains the frozen production task
registry; the existing stage dispatcher remains the production numerical-order
authority; and the Morton/SFC decomposition remains the domain-partitioning
geometry.

No scientific precision, tolerance, Riemann solver, reconstruction, CFL rule,
TreePM force criterion, AMR refinement error rule, chemistry mode, feedback
model, or restart schema is changed by this campaign.

## A. Major task/lifetime graph

`RuntimeTaskDeclaration` now carries a bounded scheduling profile for the
major production tasks:

- declared/estimated peak bytes and memory class;
- task- or stage-end release boundary;
- optional/required status;
- conservative compute, DRAM-bandwidth, and communication pressure classes;
- explicit dependencies expressed as fully-qualified task keys.

Registry freeze validates that every declared dependency exists and precedes
its consumer in the already-authoritative execution plan. This makes the
high-water task graph auditable without turning every loop into a task runtime
or changing the accepted stage order.

The built-in reference task chain is explicitly represented from pre-gravity
through drift, gravity refresh, hydro, sources, post-gravity, analysis, and
output. Optional analysis remains identifiable as optional diagnostic work.

## B. Memory dependencies and overlap policy

`runtimeTasksMayOverlap()` is a conservative legality predicate for eligible
future/concurrent execution. It rejects overlap when:

- two task resource grants conflict through a write/read-write dependency;
- both tasks are classified as high DRAM-bandwidth pressure;
- both are classified as high communication pressure;
- both are classified as high compute pressure; or
- two known task peaks cannot fit in the current governor headroom.

The current reference dispatcher is serial, so it already avoids these
pathological overlaps and preserves the existing scientific order. M2D does
not speculate a new asynchronous executor merely to exercise the metadata.

Optional science diagnostics also carry bounded starvation-recovery state. If
a light/heavy science diagnostic is due while the process governor is Red/Trip,
the analysis runtime keeps one pending bit for that diagnostic class rather
than dropping the cadence permanently. The pending work executes on the first
later AnalysisHooks stage whose pressure is below Red, even when the original
cadence is no longer due, and the bit clears only after successful execution.
Sustained unsafe pressure still defers indefinitely: fairness never overrides
the hard memory ceiling. The pending bits are owner-local scheduling metadata,
not restart/scientific truth.

## C. Deterministic headroom-aware batch/chunk sizing

The existing governor now exposes one checked deterministic batch-size policy.
It derives an item count from a governor snapshot, bytes/item, fixed reserve,
minimum item count, alignment, and a deterministic headroom-use fraction.
Unlimited-budget configurations preserve the requested maximum.

Production users now include:

- hydro active batches (4,096 B/cell declared scratch coefficient);
- remote hydro ghost transport, sized from simultaneous send+receive wire
  residency;
- stellar-evolution/feedback source batches, bounded by the existing 4096-star
  correctness/performance cap but reduced deterministically when headroom is
  tighter.

Hydro's automatic configured maximum remains 16,384 cells, therefore the
largest declared automatic hydro batch scratch is
`16,384 * 4,096 = 67,108,864 B = 64 MiB`. This is a task-local upper bound, not
a new persistent reservation.

Feedback batching preserves star/event iteration order and uses reservations
for the evolution-budget rows, feedback-event staging, retained contiguous-star
batch, and feedback spatial-index coexistence. The memory policy changes batch
residency only; source equations and event semantics are unchanged.

## D. Multi-constraint decomposition with hard rank memory

The existing `DecompositionWorkComponents` now distinguishes additional
high-water-relevant costs:

- persistent memory pressure;
- transient memory pressure;
- source/feedback event work;
- communication surface/volume proxy;
- the previously existing particle, gravity, hydro, AMR, gas, and generic work
  terms.

Existing user-visible decomposition weights remain authoritative. The new
components deliberately reuse the current memory-pressure and generic-work
weights rather than create a competing config namespace.

More importantly, rank memory is no longer only a soft weighted term.
`DecompositionConfig` may carry a hard `max_rank_memory_bytes` plus a declared
per-rank transient reserve. SFC construction cuts a range before adding an item
that would exceed the persistent allowance; an item/rank combination that
cannot satisfy the cap fails deterministically. Candidate distributed
rebalances are vetoed before migration if their predicted target peak exceeds
the hard ceiling.

The workflow derives the hard decomposition envelope from the process governor:
local decomposition-owned persistent bytes plus currently available governor
headroom gives the local admissible rank ceiling, and a collective global
minimum makes the most constrained rank authoritative. A conservative major-task
transient reserve is included, covering the maximum of hydro batch scratch,
bidirectional bounded MPI transport, and bounded feedback/source staging.
This envelope is specifically for decomposition-owned persistent state plus the
declared transient reserve; it is not added to whole-process baseline memory as
a second physical authority.

A hard current-rank violation forces a serial rebalance even when the ordinary
load/memory-imbalance thresholds and migration-fraction throttle would otherwise
suppress movement. Safety therefore outranks a soft average target.

### Synthetic clustered acceptance case

The focused unit regression places four 60 B items on rank 0 with a 20 B
transient reserve and a 140 B/rank hard peak limit. The initial peak is
`240 + 20 = 260 B` on rank 0. The memory-aware SFC split produces two items per
rank, `120 + 20 = 140 B` peak/rank, and the current hard violation forces
rebalance even with deliberately non-triggering soft thresholds. A separate
three-item case proves deterministic rejection when no legal two-rank partition
can satisfy the cap.

## E. MPI/OpenMP/node awareness

`MpiContext` now records node-local MPI rank **and node-local MPI size** using
`MPI_Comm_split_type(..., MPI_COMM_TYPE_SHARED, ...)` when MPI is available.
It also exposes a checked global minimum reduction used by the memory-envelope
consensus. The reference workflow records world rank/size, node-local
rank/size, OpenMP compile status, requested threads, and configured threads in
`runtime.topology` profiling telemetry.

No cross-rank shared-memory window is introduced in M2D. The current repository
does not expose a large immutable per-rank resource whose conversion to an MPI
shared window has a clear enough payoff to justify new lifetime/NUMA/thread-
safety complexity in this focused campaign. The new topology telemetry provides
the data needed to tune rank/thread placement without assuming one rank count is
universally optimal.

## F. Bandwidth-aware concurrency

Major task declarations carry compute, memory-bandwidth, and communication
pressure. The conservative overlap predicate prevents two high-pressure tasks
of the same limiting class from being concurrent even when their byte
reservations would fit. Because the current dispatcher remains serial, M2D does
not introduce a bandwidth-contention regression while establishing this
contract. Later task-overlap work can use measured profiler evidence to relax
specific pairs rather than defaulting to optimistic overlap.

## Focused validation completed before artifact freeze

Configuration:

```text
cmake --preset cpu-only-debug
```

Result: PASS. GCC 14.2.0 and OpenMP 4.5 detected; MPI/HDF5/FFTW/CUDA/Python are
disabled in this preset.

Compilation:

```text
cmake --build build/cpu-only-debug --target cosmosim_parallel -j2
cmake --build build/cpu-only-debug --target cosmosim_workflows -j2
cmake --build build/cpu-only-debug --target \
  test_unit_memory_governor \
  test_unit_parallel_distributed_memory \
  test_unit_runtime_module_registry -j2
```

Result: PASS. The first broad build attempts were interrupted by the execution
window while still compiling successfully; the targeted libraries then linked
cleanly. One compile-time test-development defect (the profiling component array
still had its old fixed length) was found and corrected before this freeze.

Focused tests:

```text
ctest --test-dir build/cpu-only-debug \
  -R 'unit_(memory_governor|parallel_distributed_memory|runtime_module_registry)$' \
  --output-on-failure
```

Result: PASS, 3/3 tests, 0 failures.

These focused tests cover finite/unlimited/starved deterministic sizing, task
DAG validation and conservative concurrency guards, hard rank-memory splitting,
forced safety rebalance, expanded cost telemetry, and infeasible-partition
rejection.

## Extended validation after the artifact-first source gate

The complete CPU debug build was resumed incrementally and completed:

```text
cmake --build --preset build-cpu-debug -j8
```

Result: PASS. All production, unit, integration, validation, and benchmark
targets linked. Later starvation-recovery edits touched only
`analysis_runtime.cpp` and the reference-workflow integration test; both were
rebuilt successfully before final acceptance testing.

### Final CPU test evidence

The final source tree passes every executable CPU test except the separately
long-running source-package-completeness test, which is classified below rather
than counted as a failure. Because aggregate commands can be terminated by the
external command window, completion was accumulated across the authoritative
CTest inventory and resumed ranges:

```text
ctest --test-dir build/cpu-only-debug -j8 \
  -E '^integration_source_package_completeness$' --output-on-failure
ctest --test-dir build/cpu-only-debug -R \
  '^(integration_reference_workflow|integration_star_formation_source_runtime|integration_star_formation_amr_covered_coarse|integration_star_formation_amr_refine_derefine|integration_star_formation_amr_level_equivalence|integration_star_formation_amr_reflux_ordering|integration_star_formation_amr_patch_reorder|integration_effective_ism_amr_threshold_invariance|integration_effective_ism_amr_eos_restriction|validation_convergence)$' \
  -j4 --output-on-failure
ctest --test-dir build/cpu-only-debug -I 45,48 --output-on-failure
ctest --test-dir build/cpu-only-debug -I 49,51 --output-on-failure
ctest --test-dir build/cpu-only-debug -R '^integration_reference_workflow$' \
  --output-on-failure
```

Result: **136/136 executable CPU tests PASS, 0 failures**, with
`integration_source_package_completeness` intentionally excluded from that
count. The final `integration_reference_workflow` includes an 8-step tiny
starvation-recovery regression: completed step 4 is forced to Red pressure when
light-science cadence 4 is due, the real production baseline is restored before
OutputCheck, step 5 executes the pending light-science diagnostic outside its
cadence, and later cadence continues normally. The test asserts both
`analysis.memory_pressure_deferral` and `analysis.memory_pressure_catchup`.

### Source-package completeness boundary

The final-tree command:

```text
timeout 25s bash scripts/ci/test_source_package_completeness.sh
```

returned exit code 124 after reaching 57/101 objects in the script's independent
extracted-source build. No source/package assertion failed before the timeout.
Earlier unrestricted attempts also exceeded the external command window. This
regression is therefore **NOT COMPLETED**, not failed. The changed-files M2D
bundle is generated independently and verified with `unzip -t`.

### Dependency/configuration probes

The relevant dependency paths were attempted once and then not retried when the
missing dependency was explicit:

```text
cmake --preset mpi-hdf5-fftw-debug
```

Result: BLOCKED during configure because `mpi-cxx` / `MPI_CXX` is unavailable
(`Could NOT find MPI (missing: MPI_CXX_FOUND CXX)`). No multi-rank runtime,
rank-memory max/mean, or MPI topology performance evidence is claimed.

```text
cmake --preset pm-hdf5-fftw-debug
```

Result: BLOCKED during configure because the FFTW3 serial double-precision
development library is unavailable. The production reference-workflow benchmark
also correctly aborts in the CPU-only build because TreePM production requires a
production FFT backend rather than diagnostic naive DFT.

```text
cmake --preset hdf5-debug
```

Result: PASS configure. HDF5 1.14.5 and OpenMP 4.5 were detected.

The touched HDF5-capable paths and final starvation regression were then rebuilt
and tested:

```text
cmake --build build/hdf5-debug --target \
  test_integration_reference_workflow \
  test_unit_memory_governor \
  test_unit_parallel_distributed_memory \
  test_unit_runtime_module_registry -j8
ctest --test-dir build/hdf5-debug -R \
  '^(unit_memory_governor|unit_parallel_distributed_memory|unit_runtime_module_registry|integration_reference_workflow)$' \
  --output-on-failure
```

Result: **PASS 4/4**.

### Available OpenMP topology comparison

With MPI unavailable, the executable topology comparison is single-rank
OpenMP. The CPU Debug hydro-kernel benchmark was run with the environment and
benchmark metadata both set to 1, 2, and 4 threads:

```text
OMP_NUM_THREADS=1 COSMOSIM_BENCH_THREADS=1 ./build/cpu-only-debug/bench_hydro_kernels
OMP_NUM_THREADS=2 COSMOSIM_BENCH_THREADS=2 ./build/cpu-only-debug/bench_hydro_kernels
OMP_NUM_THREADS=4 COSMOSIM_BENCH_THREADS=4 ./build/cpu-only-debug/bench_hydro_kernels
```

Observed measurement interval and reported face throughput:

| Threads | measurement_ms | face_updates_per_second | effective_bandwidth_gb_s |
| ---: | ---: | ---: | ---: |
| 1 | 4359.710685 | 240515.042342 | 0.030786 |
| 2 | 4343.176647 | 241430.658991 | 0.030903 |
| 4 | 4337.839584 | 241727.703317 | 0.030941 |

The Debug workload is effectively flat across 1–4 threads (4 threads is only
about 0.5% faster in this measurement), so M2D records the evidence rather than
assuming higher thread count is automatically better. This is not a substitute
for Release-mode or MPI rank/thread topology qualification on production
hardware.

## Scientific/reproducibility assessment

M2D changes **when/how much work may be resident simultaneously** and which SFC
partition is admissible under a hard rank-memory ceiling. It does not change the
accepted numerical operator. Legal batching preserves iteration order. The
production task graph validates the pre-existing stage order rather than
inventing alternatives. Domain migration changes ownership, not physical state,
and continues to pass through the existing migration/restart identities and
collective failure-coordination paths. Optional-diagnostic pending bits are
non-serialized runtime scheduling metadata and alter only whether a previously
due optional analysis product is retried after transient pressure; required
run-health diagnostics and solver state are unchanged.

## Handoff

M2D deliberately does not implement AIMD/predictive control, proportional
fairness, network calculus, GPU/device scheduling, or mixed precision. A later
campaign may use the new topology/task-pressure telemetry to benchmark selected
safe overlaps and device residency. The next scheduler campaign should preserve
`MemoryGovernor`, `RuntimeModuleRegistry`, the typed configuration authority,
and the SFC decomposition as the existing authorities rather than forking them.

Suggested branch: `campaign-m2d-memory-aware-runtime`

Suggested PR title: `Schedule full-physics work under explicit memory and load constraints`

## Final handoff verification

Immediately before regenerating the authoritative changed-files bundle from the
final 23-path worktree, the four directly touched M2D/reference-workflow targets
were rebuilt incrementally and rerun:

```text
cmake --build build/cpu-only-debug --target \
  test_unit_memory_governor \
  test_unit_parallel_distributed_memory \
  test_unit_runtime_module_registry \
  test_integration_reference_workflow -j8
ctest --test-dir build/cpu-only-debug -R \
  '^(unit_memory_governor|unit_parallel_distributed_memory|unit_runtime_module_registry|integration_reference_workflow)$' \
  --output-on-failure
```

Result: **PASS 4/4, 0 failures**. Ninja reported no pending compilation work,
confirming the binaries were current with the final starvation-recovery source
state. The changed-files bundle was then regenerated from an untouched extracted
copy of the supplied base ZIP and verified independently for ZIP integrity and
path/hash agreement; build trees, runtime/test outputs, caches, and the supplied
base ZIP are not included.
