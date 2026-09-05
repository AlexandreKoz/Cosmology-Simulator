# Runtime memory governance

Status: M1B foundation, 2026-08-29.

This document defines the runtime memory-governance behavior implemented by
Campaign M1B. It complements the pre-run gravity and composed DMO process
memory estimators; it does not replace them.

## Authority and ownership

The reference workflow composition root owns exactly one `core::MemoryGovernor`
for the process. `workflows::RuntimeServices` borrows that governor and selected
runtime components may in turn borrow it. The governor is not a singleton and
`core` does not depend on `workflows`.

The user-visible hard ceiling remains
`parallel.process_memory_budget_bytes`. A value of zero retains the existing
unlimited-budget meaning: runtime accounting and high-water telemetry continue,
but the governor does not reject a request because of a hard ceiling.

The production policy also reuses the existing process-level reserve semantics:

- `parallel.gravity_backend_unknown_reserve_bytes`;
- `parallel.process_mpi_unknown_reserve_bytes`;
- `parallel.process_fftw_unknown_reserve_bytes`;
- `parallel.process_hdf5_unknown_reserve_bytes`;
- `parallel.process_allocator_unknown_reserve_bytes`;
- `parallel.process_output_restart_overlap_bytes`;
- `parallel.process_memory_safety_margin_fraction`.

M1B adds no second configuration namespace and no duplicate memory budget.

## Accounting model

The governor separates four quantities that must not be added blindly to the
ordinary memory report:

1. **baseline owned bytes**: live CHUI-owned capacity observed through the
   existing `MemoryReport` but not represented by a governor commitment;
2. **governed committed bytes**: physical backing storage whose reservation was
   admitted and then committed;
3. **pending reserved bytes**: admitted demand for which backing storage has not
   yet been committed;
4. **opaque/planned reserves**: configured external-runtime reserve and the
   configured output/restart overlap allowance.

`MemoryEntry::governed_commitment` is the reconciliation marker. A physical byte
range represented by a governed entry remains visible in `MemoryReport`, but
`memoryReportBaselineOwnedBytes()` excludes it from the governor baseline because
that same range is already present in `committed_bytes`. Consequently:

```text
raw accounted demand =
    baseline owned
  + governed committed
  + pending reserved
  + external runtime reserve
  + planned overlap reserve
```

The configured process safety margin is then applied with the same
multiplicative/ceiling convention used by process preflight. `MemoryReport`
totals remain a description of owned container capacity; governor fields are a
policy/accounting overlay and are not an additional physical-memory total to be
summed with those report totals.

`MemorySubsystem` and `MemoryClass` are orthogonal. `MemorySubsystem` identifies
who owns a reported allocation; `MemoryClass` identifies its residency/purpose.
M1B memory classes are canonical persistent, persistent cache, phase resident,
scratch arena, communication, diagnostic, and external runtime.

## Reservation lifecycle

`MemoryReservation` is move-only RAII state with the following lifecycle:

```text
reserve -> pending -> backing allocation succeeds -> commit
        -> physical storage lifetime -> release/destruction
```

Pending destruction releases only the pending reservation. Committed
destruction releases committed accounting. Explicit `release()` is idempotent;
moved-from reservations are inert; a second `commit()` is a logic error.

Reservation precedes target physical allocation. If the backing allocator
throws, the pending RAII object unwinds and governor state returns to its prior
value. A reservation is retained for the physical lifetime of its allocation,
not merely for the logical phase that currently uses the bytes.

## Pressure and headroom

For finite budgets the internal M1B policy uses deterministic thresholds at 85%
and 95% of the configured hard ceiling:

- **Green**: below the Amber threshold;
- **Amber**: at or above 85% and below 95%;
- **Red**: at or above 95% while still within the hard ceiling;
- **Trip**: accounted demand is already above the hard ceiling, or a prospective
  reservation is rejected because it would cross the hard ceiling.

Pressure is observational in M1B. Adaptive batch shrinking, optional-task
postponement, cache eviction, and broad phase scheduling are M1C/later work.

## Governed physical path: monotonic scratch

`MonotonicScratchAllocator` can borrow a `MemoryGovernor`. Before each new
physical block it computes the exact intended block capacity, reserves
`MemoryClass::kScratchArena`, allocates backing storage, commits the reservation,
and stores that reservation inside the block.

`reset()` resets logical offsets only. It does not free blocks or release their
committed reservations. Reusing retained blocks therefore creates no second
reservation. Destruction of the block/allocator releases the physical
commitment. The allocator reports current logical use, retained capacity,
historical logical-use high-water, and retained-capacity high-water separately.

The production `TransientStepWorkspace` receives the process governor from the
reference workflow through `RuntimeServices`. After the collective DMO process
preflight succeeds, the production gravity runtime admits the all-particle
`uint32` index lane through this arena. Rank-local reservation/allocation
failures are passed through `FailureCoordinator` before later TreePM
collectives, so a low-headroom rank cannot simply throw while peers continue.
The admitted block is reused by later direct gravity views in the step and its
commitment survives workspace reset as long as the backing block remains
resident. Existing standalone callers may construct an ungoverned scratch
allocator for compatibility and isolated unit work.

Authoritative particle and gas-cell scheduler retained capacities are explicit
persistent `MemoryEntry` owners. Governor baseline reconciliation therefore
comes directly from the merged ownership report rather than adding hidden
scheduler byte totals outside that report.

## Runtime-resource leases remain separate

`workflows::RuntimeResourceLease` remains the stage-view freshness/capability
mechanism. It is not a memory reservation and has not been overloaded or
renamed. A task may eventually hold both a `RuntimeResourceLease` and a
`MemoryReservation` for independent reasons.

## Telemetry

The existing profiler `memory_report` object contains an additive `governor`
object when a governor snapshot is attached. It exposes hard limit, baseline,
external/planned reserves, committed, reserved, raw and safety-adjusted demand,
headroom, pressure, peak committed/reserved/accounted demand, and rejection
count. Human-readable memory reporting includes memory class and the governed
reconciliation marker.

M1B intentionally does not sample RSS/PSS or `/proc`; measured process-memory
reconciliation is an M1C acceptance task.

## M1B boundary / M1C handoff

M1B proves the central authority, accounting model, reservation lifecycle,
production ownership path, profiler integration, and one retained physical
allocation path. It does **not** yet reserve every major allocation. PM routing,
Tree/LET storage and exchange, FFT/PM grids, migration staging, output/restart
staging, halo/analysis buffers, hydro/AMR batches, and device memory remain
broad-integration work for M1C or later campaigns.

## M1C-1 runtime integration addendum (2026-08-31)

Campaign M1C-1 extends the M1B foundation without creating a second memory
authority. Major DMO/general-runtime high-water phases now use the existing
process governor:

- TreePM admits a conservative phase high-water collectively before target-scale
  cache/staging growth. Retained gravity capacity transfers atomically from the
  phase commitment into baseline ownership at phase completion.
- Runtime decomposition planning has a narrow governed lifetime and is destroyed
  before migration staging. Migration then admits transaction/staging coexistence
  collectively, and preserved scheduler state uses compact identity records
  rather than full particle-migration records.
- Snapshot/restart readback uses the existing
  `parallel.process_output_restart_overlap_bytes` allowance as predeclared
  headroom. Actual physical staging commits only the amount beyond that allowance,
  so policy reserve and physical staging are not counted twice.
- Optional science diagnostics are deferred under Red/Trip pressure. Required
  run-health diagnostics and scientific solver semantics are unchanged.

`MemoryReservation::reconcileBaselineOwnedAndRelease()` is the physical-lifetime
handoff primitive for phases whose admitted memory becomes retained capacity. It
updates baseline ownership and releases the committed reservation under one lock,
avoiding either a double-account window or an unaccounted retained-capacity gap.

Production process-memory observation now supplies optional Linux/WSL RSS,
peak-RSS/HWM, and `smaps_rollup` PSS values at workflow reporting boundaries.
The profiler distinguishes conservative governor policy demand from current
declared residency. RSS discrepancy uses baseline-owned retained capacity +
currently committed governed capacity + the configured opaque external-runtime
estimate; inactive future output/restart overlap and uncommitted reservations do
not masquerade as current residency. RSS/PSS remain diagnostic evidence; they do
not become scientific state or replace deterministic allocation admission.

Distributed telemetry reports rank-local/sum/max/mean/imbalance for declared
memory, RSS/peak RSS where all ranks expose it, and communication high-water.
Rank maximum remains the safety quantity.

Post-M1C-1 closure also replaces population-scaled migration wire
materialization with an explicit conditional field codec and reusable fragmented
packets bounded by `parallel::mpiTransportRoundLimitBytes()`. Migration governor
admission is derived from the actual transaction coexistence, native record and
dynamic sidecar capacities, scheduler/index staging, bounded packet buffers, and
maximum single-record assembly requirements.

The detailed coverage inventory and M1C-2 handoff are in
`docs/repair/m1c1_runtime_memory_integration_20260831.md`.

## M1C-2 acceptance closure addendum (2026-09-03)

M1C-2 closes the M1 source/architecture acceptance matrix in
`docs/repair/m1c2_m1_acceptance_closure_20260903.md`. The certified 8-rank PM
routing model now has an explicit regression for the 512^3 target population and
for the exact worst-case interpolation accounting expression:

```text
448 B rank metadata
+ 67,075,680 B retained send capacity
+ 67,075,680 B retained receive capacity
= 134,151,808 B
```

This is 65,920 B below the 128 MiB contract. `N_local` is absent from the
capacity model by design; source population changes round count/traffic, not the
simultaneous routing workspace. Runtime still checks actual retained wire
capacities plus rank metadata against the 128 MiB ceiling.

M1C-2 does not promote unavailable MPI/FFTW execution into runtime evidence and
does not redefine M2 hydro/AMR/full-physics memory scope.

## M2B-1 production AMR regrid transaction addendum (2026-09-05)

Production `SimulationState` refinement and derefinement are explicit governed
transactions. `ProductionAmrHydroOptions::regrid_memory_governor` must point to
the process `MemoryGovernor`; regrid fails closed when that authority is absent.
Before allocating candidate AMR authority or transaction scratch, the regrid
path computes a checked upper bound and reserves it as phase-resident memory.
Only AMR-owned cell, gas-sidecar, patch, and gas-identity state are staged;
unrelated particle/species/module populations are not copied. A reservation or
preparation failure leaves live AMR authority unchanged.

Canonical memory reporting now includes allocated capacity for the gas-cell
stable-identity records and lookup tables, pending flux-register records, and
AMR temporal-boundary histories including nested cell records. These entries
therefore participate in `memoryReportBaselineOwnedBytes()` and the existing
governor baseline reconciliation instead of remaining hidden heap state.

Regrid is permitted only at a synchronization-safe point. Active temporal
boundary history or any unresolved pending flux register causes deterministic
pre-mutation rejection. The pending-register representation does not encode all
fine-contributor topology identities, so this campaign deliberately uses the
strong global-drain contract rather than attempting an unsafe intersection-local
remap.
