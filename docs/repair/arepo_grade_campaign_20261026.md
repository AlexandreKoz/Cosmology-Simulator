# AREPO-grade gravity / MPI / memory / hierarchical campaign — 2026-10-26

**Status: IMPLEMENTED — UNCOMPILED AND UNTESTED — OWNER VALIDATION REQUIRED.**
Nothing in this campaign was built, run, tested, benchmarked or profiled; the
owner explicitly forbade execution. Status vocabulary below:
*Implemented — uncompiled and untested*, *Partially implemented*, *Existing
functionality retained*, *Not implemented*, *Scientifically unqualified*.

Baselines (verified present on `claude/vigilant-fermat-i8jvqw`): `80cda07`
(campaign round 1) and `a33ff76` (round 2). Each stage below is a separate
commit.

## Stage 1 — correctness hardening

### P0 PM routing admission (defect repaired)

**Defect (confirmed by source inspection).** `PmSolver::planDensityRouting`
and `planForceInterpolationRouting` validated span extents (and resolved box
lengths) *before* `planPmRoutingCollective`; a rank with mismatched spans
threw while its peers entered `MPI_Allreduce` — a deadlock. The force planner
also read `coordinate_source_index` rows without a range check in the measure
loop (silently skipped) and the consumers trusted a plan keyed only by data
pointers, counts and batch bytes: a plan made for mutated coordinates at the
same address, or under a different decomposition/scheme, would be consumed and
could undersize the bounded wire reservations; plan *presence* was also
assumed rank-uniform without being checked (planned ranks skipped the routing
call's reduction, unplanned ranks executed it).

**Repair.**

* Planning is one transaction: `validate_local` (extents, index ranges, box
  lengths), option validation, key construction and the workload measurement
  run inside a single `try`; then exactly one `MPI_Allreduce(MAX)` over
  `[count, M, failure_kind, world_size - rank]`. Any failure throws on every
  rank: the failing rank rethrows its original exception, peers report the
  lowest failing rank and the class (1 input validation, 2 allocation). No
  plan is published and no rank reaches routing data exchange. MPI errors in
  the reduction itself remain unrecoverable (documented).
* Plan validity is semantic, not referential. The plan stores a
  `PmRoutingPlanKey` (world size/rank, mesh shape, owned x range, batch bytes,
  assignment scheme, boundary, decomposition mode, box-length bits) and the
  agreed maxima. Each routing call **always** re-measures its own workload and
  checks, without throwing, that the key matches and current count / per-peer
  maximum / largest local round are bounded by the plan.
* Routing state machine (every rank, fixed order): (1) coordinated workload
  phase (plan withdrawn unconditionally, round model, measurement, plan
  check) → (2) `agreePmRoutingWorkload`: one `MPI_Allreduce(MAX)` over
  `[count, M, plan_present, plan_absent, plan_violated, world_size - rank]` →
  mixed presence or any violation throws identically everywhere → (3)
  coordinated preflight lease from the *current* agreed workload → (4)
  bounded rounds. Reservations are therefore always computed from present
  demand; a valid plan only certifies that the arena was sized before leasing.
* Cost: the routing call keeps one small reduction (as before round 2); the
  extra cost on the planned path is one O(N_local) workload measurement.

### P6 consensus batch (defect repaired)

`RankConsensusBatch` (round 2) is mathematically sound — `~` is an
order-reversing bijection on `uint64`, so `min(~v) = ~max(v)` and equality
holds iff `min == max` — but had two defects:

1. Its MPI element count was `2 × entries`; a rank appending a different
   number of entries produced mismatched collective counts (undefined MPI
   behaviour). The layout now has a fixed capacity (32 entries, 71 lanes) with
   structural lanes for entry count, an FNV-1a name/kind fingerprint and
   capacity overflow; mismatches throw a `logic_error` on every rank.
2. The PM-decision assembly contained rank-local `throw`s (sync-event
   authorization, scale-factor checks) between the first append and the
   batched `verify()`, so one rank could leave while peers waited in the
   collective. The assembly now runs in a `try`; failures are recorded with
   `recordLocalFailure` and reported by `verify()` on every rank (failing rank
   rethrows its own exception, peers name the lowest failing rank).

### Tests written (not run)

`tests/integration/test_pm_routing_admission_mpi.cpp`
(`integration_pm_routing_admission_mpi_{two,four}_rank`, own deadlock
watchdog plus ctest `TIMEOUT 600`): invalid extents / invalid box on one rank,
out-of-range source index on one rank, zero-particle rank, uneven and
zero-width slab ownership, skewed population on one slab, plan reuse after
coordinate mutation, routing-option change, assignment-scheme change and
decomposition change, rank-divergent plan presence, density and interpolation
agreement with the serial reference, and measured reservation ≥ observed wire
high-water.

## Stage 2 — reference-path integrity

### P1 in-place FFT audit (*Implemented — uncompiled and untested*)

Audited every consumer of the aliased padded real view: plan creation (serial
FFTW, FFTW-MPI slab and transposed-out pencil, naive DFT fallback), FFT-backed
density deposition and normalization (`density_z_stride`), compact→FFT density
copy, mean subtraction, spectral gradient (`potential_k` copy, per-axis
inverse into compact lanes), finite-difference gradient (`plane_row`),
naive forward/inverse DFT staging, plan cache move (vector move preserves the
FFTW-planned buffer), isolated/open path (separate zero-padded workspace, no
alias) and memory reports. All physical accesses use the padded stride; FFT
padding is never read as mesh data; `fourier` is never reallocated after
planning. Findings and changes:

* Raw `(local_ix·ny + iy)·real_z_stride` arithmetic replaced by
  `PlanResources::realRowBase()` in the periodic solve paths.
* New plan-creation invariant (all backends): stride `== 2(nz/2+1)` and the
  aliased view covers `local_nx · ny · stride` doubles, else `logic_error`.
* FD gradient: the per-solve `send_planes` vector (an allocation outside
  governed plan storage) now lives in the retained halo buffer (6 planes);
  `estimatePmPlanResourcesMemory` counts 6 planes.
* Alignment: `fourier` keeps `std::allocator` (16-byte on glibc, same as the
  former out-of-place real array). FFTW plans are created on the exact buffer
  they execute on, so alignment affects only SIMD codelet choice, not
  correctness. No change.
* The FFT-backed density token is already consumed before the forward
  transform (in-place r2c destroys the density, by design).

Tests: `tests/integration/test_pm_in_place_fft_layout.cpp` (serial and 2/4-rank
MPI): analytic Fourier-mode response on 16×8×8, 15×9×7, 12×10×6, 8×6×9 for
every axis and two modes; uniform-density zero force; linear normalization;
bitwise repeated solves on the cached plan; isolated point-mass antisymmetry;
in-place memory ownership; distributed versus serial spectral forces on even,
uneven (`4P+1`, `P+1`) and zero-width slabs; distributed versus serial FD.

### P4 remote group walk (defect repaired)

Audit: exactly-once accumulation holds (a group returns before its first emit
when any member is ineligible); the shared MAC uses the group box and the
minimum member acceleration, so acceptance is conservative for every member;
periodic deltas are nearest-image relative to the first member.
**Defect:** the group pair kernel excludes `r2 == 0` (needed to drop the self
pair of owned targets), but a remote target has no self pair in the local
tree, so an exactly coincident distinct source was dropped from the pair count
and `nearest_pair_r2` (forces unaffected). The remote path now uses a kernel
copy with `include_coincident = true`; owned targets are bitwise unchanged.

Test: `tests/integration/test_tree_pm_remote_group_walk_mpi.cpp` (2 and 4
ranks; fails on one rank; requires the grouped remote path and the remote
history fallback to have run): grouped versus per-target remote walk under
slab, 80/20 asymmetric and round-robin ownership, seam clusters, sparse
background, a coincident cross-rank pair, group sizes 0/16/32/64, leaf sizes
4/8/16, two softenings and mixed history fallback; component-wise forces
(1e-9 rel) and exact owned/incoming pair-count agreement in an all-direct
configuration.

### Topology caches (*Existing functionality retained*)

Refit is certified only for an unchanged source-layout generation, so the
cached row→tree-position and leaf maps (keyed by `topologyGeneration()`) are
valid across certified refits; a full build changes the generation and a
failed construction zeroes it.

### P7 checkpoint semantics (*audited; tests added*)

See `docs/restart_checkpointing.md` "Checkpoint request semantics": no
dropped/duplicated requests, PM deferral bounded by `cadence − 1` steps, time
checkpoints do not shorten steps. New cadence tests:
`testCoincidentAndPendingCheckpointEventsCollapse`,
`testTimeCheckpointCrossingSeveralEvents`.
