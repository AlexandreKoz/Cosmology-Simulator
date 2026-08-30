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
