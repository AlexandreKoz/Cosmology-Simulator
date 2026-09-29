# M48-04 — Unified Gravity Communication Arena and MPI Workspace Closure

Date: 2026-09-28. Mode: focused production memory-architecture repair.

## Previous ownership

The production TreePM solve bounded its major communication protocols but retained
the backing capacities under separate physical owners. `PmSolver` kept an
independent density-routing workspace and an independent plane-interpolation
workspace; `TreePmCoordinator` retained a separate short-range request/response
workspace; and PM slab-halo refresh allocated ordinary temporary send and receive
vectors for each force component before copying three complete component results
into the persistent force-halo cache.

That ownership model charged sequential communication maxima as though they were
simultaneously resident. It also left part of the Tree exchange's structured
request/response validation state outside the retained wire-buffer accounting.

## Shared physical owner

Production TreePM now owns one rank-local `GravityCommunicationArena` in
`TreePmCoordinator`. The arena owns one backing allocation and, when a
`MemoryGovernor` is available, one `MemoryClass::kCommunication` reservation
under the owner label:

```text
treepm.communication_arena
```

The allocation uses `std::pmr::monotonic_buffer_resource` with
`std::pmr::null_memory_resource()` as its upstream. Logical phase containers
therefore cannot overflow onto the ordinary heap. The backing allocation is
retained between gravity calls; phase changes reset only the monotonic cursor.

The arena state distinguishes:

```text
Idle
PmDensity
PmHalo
PmInterpolation
TreeExchange
```

Only one lease can be active. A phase resource can be obtained only while its
matching lease is active. Lease release invalidates all arena-backed views before
the next phase can borrow the same storage.

## Capacity derivation

The physical capacity is source-derived before first distributed use:

```text
capacity = max(
    PM density routing requirement,
    PM interpolation routing requirement,
    PM halo temporary staging requirement,
    short-range Tree communication requirement)
```

PM density/interpolation reuse the existing deterministic PM routing model and
remain individually constrained by the established 128 MiB/rank structural
workspace ceiling. Halo staging is derived from two send planes because receive
payloads now land directly in the retained force-halo cache. Tree capacity uses
`estimateTreePmExchangeMemory(...)` and includes four wire buffers, structured
request batches, exact response masks/counters, remote acceleration lanes,
rank/neighbor metadata, and reusable codec/decoded-record/identity-validation
scratch.

For the certified 512^3 / 8-rank profile with a 16 MiB PM per-peer batch and a
4 MiB Tree exchange batch, the source arithmetic is:

```text
PM density requirement          134,217,216 B   (~127.9995 MiB)
PM interpolation requirement    134,217,344 B   (~127.9996 MiB)
PM halo staging                   4,194,304 B   (4.0000 MiB)
Tree exchange requirement       158,398,648 B   (~151.0607 MiB)
----------------------------------------------------------------
shared physical arena           158,398,648 B   (~151.0607 MiB/rank)
certified internal limit        268,435,456 B   (256 MiB/rank)
```

These are architectural/source-model values, not measured RSS. The certified
profile therefore remains below the intended 256 MiB/rank CHUI-owned gravity
communication contract without weakening the PM-specific 128 MiB ceilings.

The arena is intentionally non-growing after first allocation. A later phase
contract that would require more capacity fails before use rather than silently
allocating a second backing workspace.

## PM density and interpolation leases

`PmSolver` no longer retains independent production density/interpolation wire
owners. Each distributed routing call constructs phase-local `std::pmr::vector`
metadata and wire buffers on the currently attached arena resource. Existing PM
wire records, routing destinations, deterministic request ordering,
count/displacement rules, exchange epochs, payload validation, per-peer limits,
and collective participation are unchanged.

The interpolation path preserves its existing two-buffer request/response reuse:
request and response phases reuse the same send/receive wire objects inside one
`PmInterpolation` lease.

Standalone `PmSolver` use remains supported. If no external TreePM arena is
attached, the solver lazily creates one bounded 128 MiB PM-only fallback arena.
Attaching the production arena destroys any standalone fallback first, so a
full-sized fallback cannot coexist as a hidden duplicate owner.

## PM halo staging and cache publication

PM slab-halo exchange has a production `...Into` path whose two send staging
vectors are PMR-backed by the `PmHalo` lease. Receive payloads are written
directly into caller-provided retained cache spans, and all posted MPI requests
complete before the function returns and the lease can reset.

`PmGridStorage::ForceHaloCache` remains separate physical ownership because its
six final X/Y/Z left/right lanes must survive until interpolation. Refresh is now
transactional:

1. invalidate the cache and prepare the six retained output lanes;
2. exchange X with shared staging and commit directly into the X lanes;
3. reset/reuse staging for Y, then Z;
4. verify common depth/peer metadata;
5. publish `valid=true` and the exchange sequence only after every component
   succeeds.

Any preparation, exchange, or commit failure leaves `valid=false`. The old
simultaneous `force_x_halo` + `force_y_halo` + `force_z_halo` temporary-result
overlap is removed from the production TreePM path.

## Short-range Tree exchange lease

After PM interpolation completes, the same physical backing allocation is
borrowed as `TreeExchange`. The following formerly retained or phase-local large
communication state is arena-backed:

- request/receive and response-send/response-receive wire buffers;
- communicator-wide count/displacement arrays;
- sparse-neighbor count/displacement arrays;
- per-peer structured request batches;
- response expected/seen masks and per-target response counters;
- remote acceleration accumulation lanes;
- reusable request/response codec buffers;
- decoded request/response storage and response-construction storage;
- exact duplicate-detection identity scratch;
- peer-participation/requested-peer metadata.

Duplicate request detection remains exact, but the node-based temporary hash set
is replaced by arena-backed identity scratch plus sorting used only for
validation. This does not change request ordering or force accumulation order.
Response identity, duplicate rejection, finite-value checks, and exact coverage
validation remain intact.

Persistent LET/domain hierarchy state, sparse peer-graph cache/communicator,
tree nodes, PM fields, FFT plans, force-halo cache, and OpenMP DFS stack storage
remain outside the communication arena because their lifetimes are not the same
as transient request/response staging.

## MPI completion and failure boundaries

PM density/interpolation leases are released only after their solver call
returns, which is after collective communication and local decode/validation
using the phase storage. PM halo staging is released only after the blocking
halo helper has completed every posted MPI request. The Tree exchange lease
outlives all arena-backed request/response containers and is destroyed only
when the distributed residual exchange scope completes.

Arena admission and PM/halo preparation use the existing coordinated failure
paths. The arena's null upstream converts any underestimated capacity into a
hard allocation failure rather than an untracked heap workspace. The source
capacity model reserves deterministic structural headroom for alignment and PMR
bookkeeping.

## Physical versus logical accounting

The production memory report exposes the physical owner exactly once as:

```text
treepm.communication_arena
```

Per-phase PM density, PM halo, PM interpolation, and Tree exchange high-water
entries are emitted as logical/non-owning telemetry with zero owned capacity.
`PmSolver` likewise reports logical PM routing high-water and claims physical
capacity only when its standalone fallback arena actually exists.

Gravity preflight now charges:

```text
gravity.estimate.communication_arena = max(PM density, PM halo, PM interpolation, Tree exchange)
```

and charges the six-lane force-halo cache separately. It no longer adds separate
physical PM-routing and sparse-Tree-exchange maxima to the same modeled peak.
The runtime phase reservation excludes the arena because the arena itself is a
first-class MemoryGovernor commitment, preventing double reservation.

## Reproducibility and numerical scope

This campaign changes storage ownership and lifetime scheduling only. It is
intended to preserve PM assignment, PM interpolation, TreePM split mathematics,
tree traversal, softening, wire schemas/versions, request IDs and epochs,
collective ordering, response validation, and force accumulation order. The
identity-scratch sort is validation-only and does not reorder force arithmetic.

No empirical bitwise/rank-equivalence claim is made in this implementation pass.
Runtime and numerical qualification are intentionally deferred.

## Deferred work

This closure does not collapse homogeneous DMO source/target state, change the
uniform-rung scheduler, redesign migration communication, stream restart
verification, compact scalar tree epsilon, alter FFT ownership, or model opaque
MPI/FFTW allocator internals. Those remain later campaign/qualification work.

Build/tests: NOT RUN — intentionally deferred by campaign instruction.
