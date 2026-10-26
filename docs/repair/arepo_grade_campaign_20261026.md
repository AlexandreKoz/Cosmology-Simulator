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

## Stage 3 — fine-tick architecture

### Source-state contract (audited; prediction made fail-closed and cheaper)

Authoritative state and epochs, as implemented:

| Quantity | Owner / epoch | Fine-tick behaviour |
| --- | --- | --- |
| Positions `x` | `SimulationState::particles`, epoch `last_drift_time_code / last_drift_scale_factor` (per particle, or one common epoch for homogeneous DMO; restart-persisted in `/state/particle_sidecar/`) | Hierarchical KDK canonically drifts **every** source each fine tick (`hierarchical all-source drift`), then `updateAllParticleDriftEpoch` (O(1) for homogeneous DMO). |
| Velocity lane | global KDK: canonical momentum `p` (drift `∫dt/a²`); hierarchical KDK: `u = p/a` at the common epoch, Hubble drag in the drift | Kicks touch active rungs only. |
| Predicted sources | gravity scratch only (`rebuildOwnedParticleCompactView`), never written to state | Global KDK cache rebuilds predict inactive sources from their epoch. |

Necessity analysis of the per-fine-tick O(N) work (no measurements; source
reasoning only):

* **Required:** an O(N) source-position update (all sources move between
  force epochs; any lazy scheme still materializes predicted coordinates for
  the tree), O(N) tree refit or rebuild, and LET exchange of moved boundary
  sources. GADGET-2/4 and AREPO likewise drift or predict every particle per
  step. The remaining alternative — node-level analytic prediction
  (GADGET-2 style dynamic tree update with velocity-bounded extents) — is a
  separate tree design and **not implemented**.
* **Avoidable and removed:**
  * Per-source 64-point cosmological drift quadrature in the predictor, now
    memoized exactly per `(source epoch, evaluation epoch)` bit pattern
    (`source_prediction_drift_integrals` telemetry); predicted coordinates
    are bitwise unchanged.
  * The per-solve O(N) uniform-mass scan (round 2), now certified by
    `TreePmOptions::source_mass_generation`. The gravity runtime passes
    `particleIndexGeneration() + 1` for borrowed homogeneous DMO, where masses
    change only by reorder or migration, and 0 (rescan every solve)
    otherwise. Telemetry: `mass_uniformity_scans / reuses`,
    `uniform_source_mass_lane_released`.
  * Topology-only metadata: row→tree-position and leaf maps are cached by
    topology generation (round 2). Visit ordering is O(T) in the active
    target count.
* **Fail-closed contract:** the global-KDK momentum predictor is invalid for
  hierarchical `u = p/a` lanes; it now throws if a hierarchical source would
  need a non-zero predicted drift (cannot happen with the all-source drift,
  and so becomes an enforced invariant instead of an assumption).
* Repeated consensus: one fixed-size batched reduction at gravity entry and
  one for PM decisions (Stage 1).
* Certified refit (`treepm_topology_refit_enabled`) remains opt-in; with it
  off every fine tick rebuilds the tree (O(N log N)). The performance deck
  enables it.

Tests: `tests/integration/test_tree_pm_fine_tick_caches.cpp`. It covers:

* Retained caches versus an eager fresh-coordinator reference across 12 ticks
  with active fractions of 100%, 25%, ~3% and one target.
* Certified refit and topology-cache reuse, plus certified mass reuse.
* Seam-crossing refit rejection, which forces a rebuild.
* A dense-row reorder, which forces a rebuild with reason 8.
* Uniform scalar mass versus the explicit mass lane, and exact mass scaling.
* OpenMP thread-count determinism with and without the group walk.

**Not implemented:** per-particle lazy drift, which would leave inactive
canonical positions at their own epochs. It would need per-bin epoch
derivation in the drift stage, analysis and output. The saving would be only
the canonical write pass, since the tree still reads every predicted
position. The current all-source drift is exact and is kept.

## Stage 4 — hierarchical PM integration (*Implemented — uncompiled and untested; scientifically unqualified*)

The existing hierarchical scheme was already a defined two-level split:
r-RESPA with PM half kicks at both endpoints of every coarse block and tree
kicks on the fine rungs. Round 2's `stale_within_cadence` is a global-KDK
field-reuse approximation and stays restricted to `hierarchical_max_rung = 0`.

The new policy `multirate_outer_kick` adds an outer PM level spanning
`M = treepm_update_cadence_steps` coarse blocks. Equations, the cycle
planning rules, cache validity and restart behaviour are derived in
`docs/time_integration.md` ("Multirate outer PM level"). Changes:

* `core::HierarchicalPmCyclePlan` drives the block driver's PM surfaces: the
  opening kick only in a cycle-opening block, refresh plus closing kick only
  in a cycle-closing block. A disabled plan gives the bit-for-bit reference.
* `StepBoundaryState::multirate_pm_cycle_interior`; interior boundaries are
  not restart- or output-safe.
* `TimeCoordinator` plans each cycle:
  * Δc is frozen for the cycle, and the PM displacement and boundary limits
    apply to the whole cycle.
  * The last block lands exactly on the planned end.
  * The cycle is truncated to the step budget.
  * Interior blocks start from a Tree-only endpoint reuse or bootstrap.
* `GravityRuntime::hierarchicalShortRangeEndpointReusable`.
* `OutputRestartRuntime` holds pending output at interior boundaries.
* The PM sync cadence is 1, so there is one refresh per cycle.
* Config: the enum value `kMultirateOuterKick` and the key value
  `multirate_outer_kick`, with validation and normalized dump.
* Telemetry: event payload `pm_policy` / `pm_cycle_blocks` /
  `pm_cycle_block_index`, and the counters
  `time.hierarchical_multirate_pm_cycles`,
  `time.hierarchical_multirate_tree_bootstrap` and
  `time.hierarchical_multirate_quantum_clamped_blocks`.

Restart: no schema change (checkpoints only at cycle ends).

Tests:

* `test_hierarchical_treepm_workflow` (serial and MPI registrations) gains
  multirate checks:
  * Synchronized endpoints.
  * Refinement convergence toward the fine global reference, with a
    provisional bound of `1e-3`.
  * Fewer PM refreshes than the exact scheme at the same Δc.
  * Bit-identical resume from a cycle-end checkpoint.
  * Step-budget truncation of a cycle.
* `test_config_parser`: validation of `multirate_outer_kick`.

**Not implemented:**

* Analytic two-body and energy-drift gates. The workflow comparison above is
  the trajectory gate.
* An MPI decomposition-consistency run specific to multirate. The existing
  hierarchical MPI registrations run the same test file on 2 ranks.
* Adaptive mid-cycle termination (see the docs for why).

## Stage 5 — memory architecture

### Lossless structural compaction (audit; source-level verdicts)

| Candidate | Verdict |
| --- | --- |
| Species-wide constant mass | **Implemented** (round 2 lane release + Stage 3 certified reuse, no per-tick scan for DMO). |
| Uniform softening | Existing: scalar `uniform_epsilon` path; per-source lanes only when non-uniform. |
| Oversized local indices | Existing: `TreeLocalIndex = uint32` with `checkedIntegralNarrow` at every global→local boundary (48 sites); TreePM rejects > 2³²−1 local sources. |
| Duplicate unwrapped periodic lanes (`treepm.periodic_tree_source_*`, 24 B/source) | Retained by design. They are rewritten every build or refit; releasing them per tick would turn into allocation and page-fault churn on every fine tick. A zero-copy identity unwrap (only when the canonical input needs no shift) would help only seam-free axes. Not implemented. |
| Hierarchical total-force history (`m_particle_accel_*`, 24 B/particle) beside the Tree/PM component lanes | **Not removed.** Recomputing total as Tree+PM differs from the coupling's accumulated total at roundoff. It would change rung decisions and the restart-persisted v24 endpoint bitwise. A candidate for a versioned follow-up. |
| `m_hierarchical_force_generation` (uint64/particle) | Retained: the exact-epoch kick guard, restart persisted. |
| FD gradient halo staging | **Fixed** (Stage 2): the per-solve allocation moved into the governed plan storage. |
| PM real array | **Removed** (round 2 in-place FFT; Stage 2 invariant). |

### Canonical identity (audit)

Authoritative IDs are `uint64` (`particle_sidecar.particle_id`,
`gas_cell_id`, `parent_particle_id`) through migration, HDF5 snapshots
(`ParticleIDs`), restart and species conversion. No ID was narrowed in this
campaign. Compact indices (`TreeLocalIndex`, `uint32` rows) are rank-local, are
never persisted as identity, and are produced only by checked narrowing.

### Experimental precision policies (*Not implemented*; design contract)

No lossy representation was introduced, because none could be made both
useful and safe without owner error budgets. The design contract for the
first candidate, kept for the owner:

* **Candidate: float32 hierarchical total-force history.** Quantity: scale-free
  comoving acceleration, code units, used only for rung criteria and
  relative-MAC tolerances, never for kicks. Dynamic range: about 10⁻⁸–10⁸ of
  the force scale, within float32 normal range. Rounding: round-to-nearest
  narrowing, relative error ≤ 2⁻²⁴. Conversion boundaries: written at the
  scatter after a solve, widened in `criterion()` and in the previous-
  acceleration lane. Accumulation: none, since history is overwritten every
  solve. Forbidden in global KDK, where kicks read this lane. Periodic
  boundaries: not applicable. Restart: persisted widened (exact), so resumed
  runs are bitwise identical to uninterrupted runs. Saving: 12 B/particle.
  Effect: rung decisions near bin thresholds may flip relative to the double
  reference.
* Rejected for now: fixed-point periodic coordinates (they would change every
  force), narrowed LET packets (they would change remote forces), and compact
  PM routing records (they would change deposited density).

### MemoryGovernor coverage of campaign allocations

| Allocation | Accounting |
| --- | --- |
| PM routing reservations | Communication arena, governed. Sized by the routing plan before configure. |
| FD halo planes | PM plan storage. Live report and estimate (6 planes). |
| Topology caches, radix scratch | Live report plus preflight estimate (round 2). |
| Incoming remote group order/keys | Live report; **new** preflight term `gravity.estimate.treepm_incoming_group_order`, bounded by `tree_exchange_batch_bytes`. |
| Drift-factor memo | Fixed 16-entry stack array (no heap). |
| Consensus batch | Fixed-size stack arrays (no heap). |
| Multirate cycle | Scalars only. |
