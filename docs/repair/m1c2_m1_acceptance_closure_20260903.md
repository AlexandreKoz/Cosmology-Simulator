# M1C-2 M1 memory acceptance closure

Status: source/architecture closure achieved; distributed runtime validation blocked by environment.
Date: 2026-09-03
Campaign: M1C-2 M1 Acceptance Closure

## Scope

M1C-2 is the final acceptance pass over the M1A -> M1B -> M1C-1 memory
architecture. It does not redesign M2 hydro/AMR/full-physics representation.
The repository state audited here already includes the later MPI count,
collective-preparation, sparse-participation, and M1C-1 post-campaign repairs.

The M1 objective remains: minimize wall time inside the configured safe process
memory ceiling while preserving the accepted numerical methods, tolerances,
restart semantics, and reproducibility contracts.

## Campaign lineage outcomes

| Campaign | Outcome in current source | Acceptance interpretation |
| --- | --- | --- |
| M1A structural memory closure | Achieved in current source | PM density routing no longer retains an `O(stencil_size * N_local)` xyz-contribution graph; density and reverse interpolation use bounded x-plane wire records and compact scheduler state. |
| M1B memory-governor foundation | Achieved in current source | One process `MemoryGovernor`, reserve-before-allocate RAII, capacity-based accounting, no double counting of governed physical ranges, retained-capacity lifetime truth, and deterministic finite-budget rejection are present. |
| M1C-1 runtime integration | Achieved in current source | Gravity, migration, output/restart readback, pressure scheduling, RSS/PSS observation, distributed rank-memory telemetry, and collective-safe large-phase admission are integrated; post-campaign migration wire residency and TreePM edge contracts are closed. |
| M1C-2 acceptance | Source/architecture achieved | The matrix below closes every mandatory source/architecture gate. Multi-rank runtime gates remain dependent on an MPI+FFTW-capable environment and are not promoted from source evidence alone. |

## M1 acceptance matrix

| Gate | Verdict | Evidence class | Current evidence |
| --- | --- | --- | --- |
| 1. No `O(stencil_size * N_local)` PM deposition contribution graph | PASS | source-verified + CPU test | `src/gravity/pm_solver.cpp` deposits owner-local x planes immediately and encodes at most one remote record per destination x-plane group. `integration_pm_periodic_mode` additionally rejects regression toward xyz-stencil materialization when MPI is executable. |
| 2. Certified 512^3 / 8-rank / 512^3-mesh / TSC PM deposition communication/scratch <=128 MiB/rank | PASS | analytically verified + unit test | `modelPmRoutingCapacity()` is independent of `N_local`, uses actual rank-metadata capacity, record alignment, and checked arithmetic. The 8-rank interpolation-side worst model is 134,151,808 B; exact derivation below. Runtime routing also checks actual retained wire capacities plus metadata against 134,217,728 B. |
| 3. Force and potential interpolation bounded simultaneous communication workspace | PASS | source-verified + unit/integration test | Force/potential share two reusable request/response wire buffers, compact responses in place, use the same aggregate model, and enforce the same 128 MiB runtime high-water guard. |
| 4. No population-scale dynamic scheduler strings | PASS | source-verified + CPU test | `TimeBinHotMetadata` and persistent scheduler state contain compact integer vectors only. Candidate diagnostic labels are transient; `testSchedulerCandidateLabelsAreTransient` proves a 1 MiB label does not change scheduler logical/capacity/high-water bytes. |
| 5. Large target-scale allocations reject before backing commitment | PASS | source-verified + CPU test | `MemoryGovernor::reserve()` precedes governed allocation; `MonotonicScratchAllocator` reserves before block allocation. Exact-limit/one-byte-over tests pass. |
| 6. Physical retained capacity stays committed until storage release | PASS | source-verified + CPU test | Governed monotonic scratch retains its committed reservation across logical `reset()` and releases only with backing block destruction. |
| 7. No physical byte range double-counted in governor admission | PASS | source-verified + CPU test | `MemoryEntry::governed_commitment` excludes already-governed physical capacity from baseline while ownership reporting still retains the entry. Atomic baseline/commit reconciliation is tested. |
| 8. Impossible finite-budget configurations fail before catastrophic solver staging | PASS | runtime-verified CPU | `integration_reference_workflow` uses a 1-byte process budget and verifies failure at process/gravity admission before large TreePM staging. |
| 9. Rank-dependent admission/allocation failure collective-safe | PASS for source architecture; runtime multi-rank pending environment | source-verified | Gravity, decomposition, migration, snapshot/restart admission failures are reduced through `FailureCoordinator` before subsequent collective phases. Later MPI sparse-participation/count repairs are preserved. |
| 10. Real process residency telemetry and declared/observed separation | PASS | runtime-verified Linux CPU | `/proc/self/status` RSS/HWM and optional `smaps_rollup` PSS are sampled. The profiler separates current declared residency from conservative governor policy demand and records unexplained resident bytes. |
| 11. Repeated equivalent work stabilizes retained capacity | PASS for CPU-owned scheduler/tree/scratch; distributed PM runtime pending environment | runtime-verified CPU + source-verified | Scheduler capacity/high-water, tree node capacity, and monotonic scratch retained capacity have repeated-use tests that require stabilization rather than monotonic growth. PM routing retains bounded reusable exchange buffers; multi-rank repeated routing execution requires MPI. |
| 12. Scientific correctness/restart/rank-count tolerances unchanged and pass where executable | PASS where executable | runtime-verified CPU; dependency paths recorded below | No force criterion, CFL rule, precision, tolerance, restart schema, or rank-count contract is modified by M1C-2. Dependency-enabled results are recorded in the validation section. |

## Critical 128 MiB PM routing audit

The M1 gate excludes persistent PM mesh/FFT plan arrays and opaque MPI/FFTW
internals. Those remain independently reported/budgeted; they are not hidden
inside the communication number.

The source uses:

```text
M_target = 128 * 1024 * 1024 = 134,217,728 B
M_policy_model_limit = M_target - 64 KiB = 134,152,192 B
```

At 8 ranks the interpolation routing workspace has twelve rank-sized `int`
metadata vectors plus one rank-sized `size_t` cursor vector. On the validated
64-bit GNU/Linux ABI used by the unit test:

```text
M_rank_metadata = 12 * 8 * sizeof(int) + 8 * sizeof(size_t)
                = 12 * 8 * 4 + 8 * 8
                = 448 B
remote_peers = 7
wire_record = 96 B
configured_per_peer_max = 16 MiB
```

`modelPmRoutingCapacity()` divides the remaining aggregate payload budget over
the two simultaneously live directions, then rounds down to a whole wire
record:

```text
B_peer_effective = 9,582,240 B
M_send_capacity  = 7 * 9,582,240 = 67,075,680 B
M_recv_capacity  = 7 * 9,582,240 = 67,075,680 B

M_workspace = 448 + 67,075,680 + 67,075,680
            = 134,151,808 B
            = 127.9371337890625 MiB

margin_to_128_MiB = 134,217,728 - 134,151,808
                   = 65,920 B
```

This is not a logical-payload-only claim. The routing preflight obtains rank
metadata from container **capacity**, wire growth uses a fresh exact-size vector
to avoid geometric-growth retention, and every exchange round recomputes the
actual retained send/receive capacities plus metadata. If that actual value
exceeds 128 MiB, production throws before accepting the round as valid.

The 64 KiB policy reserve is therefore not being used to conceal known dynamic
metadata; known rank metadata is in the 448-byte expression. It is fixed
headroom beneath the acceptance ceiling. Record alignment contributes the
remaining 384 B of the 65,920 B analytical margin.

Density routing has fewer rank metadata vectors and a slightly smaller model:
134,151,680 B. Force and potential interpolation use the 134,151,808 B case, so
that is the conservative M1 analytical routing value.

The certified local particle count is 512^3 / 8 = 16,777,216 particles/rank,
but `N_local` is deliberately absent from the capacity model. Increasing the
particle population increases the number of bounded rounds and total traffic,
not the simultaneous routing workspace.

## Governed-allocation coverage

M1 closes the target DMO/general-runtime peak paths without claiming M2 full-
physics compaction:

- canonical `SimulationState` and scheduler capacity: baseline-owned;
- monotonic step scratch: governed physical blocks;
- TreePM phase/cache/tree/PM growth: collective phase admission followed by
  atomic retained-baseline reconciliation;
- bounded PM routing: governed by the enclosing gravity high-water and the
  independent 128 MiB routing contract;
- decomposition planning: governed phase-resident memory;
- migration transaction and bounded fragmented wire packets: governed
  communication memory;
- snapshot/restart readback: governed diagnostic memory consuming predeclared
  overlap allowance without double counting;
- optional analysis: deferred under red/trip pressure while required solver and
  health work remains mandatory;
- opaque MPI/FFTW/HDF5/backend allocations: configured external-runtime reserve
  plus observed-process reconciliation.

Hydro/AMR batching, chemistry fidelity tiers, star/BH/feedback sidecar research,
and device-memory governance are M2/later work and are not required for M1
acceptance.

## Estimator / governor / RSS reconciliation

The three views answer different questions and must not be added blindly:

1. `MemoryReport`: CHUI-owned logical/capacity/high-water ownership truth.
2. `MemoryGovernor`: admission-policy baseline + committed + pending + opaque /
   planned reserves, with governed ownership ranges reconciled out of baseline.
3. process observation: current RSS, peak RSS/HWM, optional PSS, and unexplained
   residency relative to **current declared residency**, not future policy-only
   headroom.

The Linux unit path verifies that real RSS/HWM are available and nonzero and
that PSS, where exposed, is positive. The accounting unit path separately
proves that pending reservations and future output-overlap reserve do not
masquerade as current resident memory.

## Scaling law probes

### Analytical target law

For the certified 8-rank PM routing policy:

```text
M_routing(N_local) = M_rank_metadata + 2 * (ranks - 1) * B_peer_effective
```

There is no `N_local` term. Therefore the 64^3, 96^3, 128^3, and 512^3 global
particle probes all have the same worst-case PM routing capacity at fixed rank
count/configuration; only bounded-round count/total traffic changes.

For eight evenly populated ranks:

| Global population | Local particles/rank | Certified PM routing bound |
| ---: | ---: | ---: |
| 64^3 = 262,144 | 32,768 | 134,151,808 B worst-case interpolation workspace |
| 96^3 = 884,736 | 110,592 | 134,151,808 B worst-case interpolation workspace |
| 128^3 = 2,097,152 | 262,144 | 134,151,808 B worst-case interpolation workspace |
| 512^3 = 134,217,728 | 16,777,216 | 134,151,808 B worst-case interpolation workspace |

This is an analytical capacity proof, not a claim that all four populations were
executed with MPI in this environment.

### Retained-capacity stability

CPU tests cover three independent retained owners:

- scheduler: repeated equivalent candidate/reconciliation work leaves retained
  capacity and retained-capacity high-water unchanged;
- tree: repeated rebuilds of the same topology keep node capacity high-water at
  the first stabilized value;
- monotonic governed scratch: logical reset does not allocate or reserve a
  second physical block when retained capacity suffices.

Distributed PM exchange buffers use the same retained-capacity reuse model, but
multi-rank repeated execution is reported only if the MPI configuration is
available.

## Rank imbalance evidence

The production profiler emits local/sum/max/mean/max-to-mean for governor
accounted demand, RSS, peak RSS, and communication high-water. Rank maximum is
the safety quantity. OS-derived distributed metrics are valid only if all ranks
provide the observation. Multi-rank runtime values are not fabricated when MPI
is unavailable.

## Scientific-correctness invariants

M1C-2 changes only acceptance documentation and strengthens the existing PM
capacity-model regression. It does not alter:

- PM assignment/deconvolution mathematics;
- TreePM split, opening criteria, softening, or force convention;
- time integration or CFL/gravity timestep rules;
- precision or reproducibility policy;
- HDF5/restart schema;
- numerical validation tolerances.

## Validation record

### Environment used for this acceptance pass

```text
OS/kernel: Linux 6.18.35 x86_64
CMake: 3.31.6
C++ compiler: GNU 14.2.0
Visible CPUs: 5
Visible MemTotal: 6,090,748 kB
HDF5: 1.14.5
FFTW3 development package: unavailable
MPI C++ development package: unavailable
```

This is a constrained validation sandbox, not the project's 48 GiB reference
workstation. Target-scale 512^3 memory figures below are therefore analytical
capacity evidence, not a claim of target-scale execution.

### Focused pre-artifact acceptance

```text
bash scripts/ci/check_repo_hygiene.sh
  PASS before generated build trees existed.
  (The ZIP did not preserve executable bits, so the script was invoked through bash.)

cmake --preset cpu-only-debug
  PASS; OpenMP 4.5 enabled, HDF5/FFTW/MPI disabled.

cmake --build --preset build-cpu-debug -j8 --target <focused M1 targets>
  PASS.

ctest --test-dir build/cpu-only-debug --output-on-failure -R <focused M1 regex>
  10/10 PASS.
```

The required changed-files acceptance checkpoint was generated immediately
after this green focused set, before the extended validation below.

### Full CPU debug build and test suite

```text
cmake --build --preset build-cpu-debug -j8
  PASS; all configured CPU debug targets linked successfully.

ctest --preset test-cpu-debug --output-on-failure -j4
  136/138 tests completed and passed before the outer 300 s execution ceiling
  terminated CTest. No test failure was reported.

ctest --test-dir build/cpu-only-debug --output-on-failure \
  -R '^(integration_effective_ism_amr_eos_restriction|integration_source_package_completeness)$' -j2
  integration_effective_ism_amr_eos_restriction: PASS.
  integration_source_package_completeness: execution interrupted by the 120 s outer ceiling.

ctest --test-dir build/cpu-only-debug --output-on-failure \
  -R '^integration_source_package_completeness$'
  repeated with 300 s outer ceiling: no failure output, but did not complete.
  repeated with COSMOSIM_PACKAGE_BUILD_JOBS=16 and 300 s outer ceiling: same result.
```

Therefore 137/138 configured CPU tests were explicitly observed PASS. The only
uncompleted CPU test is the release/source-package completeness integration
test; inspection shows it stages, hashes, ZIPs, extracts, re-hashes, runs
hygiene, then performs a fresh extracted-tree configure/build. Its incomplete
status is an execution-window limitation of this acceptance environment, not a
reported source/test failure and not an M1 scientific-memory gate.

### Measured memory/resource probes

```text
./build/cpu-only-debug/bench_memory_accounting_overhead
  particles=200000 cells=100000
  known_owned_bytes=34,198,583
  baseline_rss_bytes=1,855,488
  allocated_rss_bytes=36,016,128
  rss_delta_bytes=34,160,640
  known_owned_to_rss_delta_ratio=1.00111
  peak_rss_bytes=36,167,680

./build/cpu-only-debug/bench_pm_solver
  backend=naive_dft threads=1
  warmup_solves=1 measured_solves=8
  setup_ms=621.401 steady_ms=4,872.66
  plan_cache_size=1 plan_build_count=1
  (diagnostic CPU backend only; not FFTW or distributed PM evidence.)

./build/cpu-only-debug/bench_core_migration_memory
  initial_cells=50000 outbound_cells=12500 inbound_cells=12500 final_cells=50000
  steady_state_report_bytes=6,250,000
  exchange_payload_capacity_bytes=2,250,000
  baseline_rss_bytes=1,863,680
  steady_rss_bytes=14,721,024
  precommit_rss_bytes=16,953,344
  postcommit_rss_bytes=21,860,352
  peak_rss_bytes=31,444,992
  peak_to_steady_rss_ratio=2.13606
```

The accounting probe is especially relevant to Gate 10: CHUI-owned reported
capacity and the observed incremental RSS agree to about 0.11% in this probe,
while still remaining separate quantities. The migration probe exposes the
transactional coexistence peak rather than treating only steady logical state
as physical truth.

### HDF5 path

```text
cmake --preset hdf5-debug
  PASS; HDF5 1.14.5 found.

cmake --build --preset build-hdf5-debug -j8 --target stage1_runtime_truth_targets
  PASS.

ctest --test-dir build/hdf5-debug --output-on-failure -j4 -R <M1/restart HDF5 regex>
  16/16 PASS.
```

The 16-test set includes snapshot HDF5 roundtrip, restart checkpoint roundtrip,
DM-only/TreePM/hydro/multirate/output/stochastic/AMR restart equivalence,
reference-workflow hydro row-order and restart equivalence, softening ownership,
and provenance roundtrip.

### Dependency-blocked production PM/MPI paths

Each unavailable configuration was attempted once as required:

```text
cmake --preset pm-hdf5-fftw-debug
  BLOCKED: COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library
  was not found; fftw3.pc is unavailable.

cmake --preset mpi-release
  BLOCKED: Package 'mpi-cxx', required by 'virtual:world', not found;
  Could NOT find MPI_CXX / MPI.

cmake --preset mpi-hdf5-fftw-debug
  BLOCKED at MPI discovery for the same missing MPI C++ development support.
```

Consequently this environment cannot execute multi-rank PM routing, rank-count
equivalence, distributed admission-failure convergence, or combined
MPI+FFTW+HDF5 runtime validation. Those results are not inferred from CPU
execution.

## M1 final acceptance verdict

```text
M1 SOURCE/ARCHITECTURE CLOSURE: ACHIEVED
M1 DISTRIBUTED RUNTIME VALIDATION: BLOCKED BY ENVIRONMENT
```

No mandatory source/architecture gate remains open. The remaining distributed
validation debt is precise: rerun the current MPI/FFTW-labelled TreePM, PM,
rank-equivalence, failure-coordination, and repeated-routing tests in an
environment providing MPI C++ and FFTW3 development/runtime support.

## Explicit M2 handoff

M1 acceptance does not claim that a 512^3 uniform full-physics hydro calculation
fits the workstation envelope. M2 owns full-physics memory architecture:
active-batch hydro primitives/gradients/reconstruction/flux scratch, sparse AMR
metadata and regrid coexistence, chemistry fidelity tiers, species-specific
star/BH/feedback sidecars, analysis/halo memory discipline beyond the DMO
acceptance surface, and later device-memory governance. M2 must preserve the M1
single-governor authority and must not reopen a population-scaled PM routing
representation.
