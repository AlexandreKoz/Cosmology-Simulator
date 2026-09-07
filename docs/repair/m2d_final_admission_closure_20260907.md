# M2D final admission closure — focused repair

Date: 2026-09-07. Base: Cosmology-Simulator-main(15)(8).zip.
Mode: repair. This report records the actual patch, not a declaration of full
M2D acceptance. The source ZIP and its current AGENTS.md remain authoritative.

## Findings addressed

- M2D-2-A: analysis now obtains and holds an owner-managed diagnostic lease
  before materializing population validation scratch, FFT meshes, bin arrays,
  and output data. A checked allocation model covers the real/complex mesh
  coexistence and the built-in FFT worker lines and reductions. The optional
  admission decision is collective; a hard-limit rejection defers/coalesces
  the due product, while allocation, numerical, and I/O errors remain fatal.
  Required health is never silently skipped. Execution failures are coordinated
  before another optional stage begins. The accepted numerical estimator,
  precision, and scientific policy have not changed.
- M2D-2-B: the analysis owner publishes a state-dependent required-health
  estimate to the existing runtime task registry. Its preflight is released
  before the owner acquires the physical lease; optional science is admitted
  separately and cannot turn into a mandatory whole-stage rejection. The
  existing serial stage order and conservative unknown-peak overlap refusal
  remain intact. A complete peak model for every gravity, source, and
  output/restart owner is not established by this repair. Do not enable
  nonserial full-physics execution on the strength of this partial contract.
- M2D-2-C/D: the previous SFC repair is preserved. No new MPI algorithm,
  numerical operator, configuration key, restart dataset, or schema version is
  introduced. Distributed runtime certification and full-process production
  memory measurements remain explicit acceptance requirements.

## Ownership and scale

The diagnostic validation workspace owns one uint64 particle-ID lane, three
uint8 particle markers, and one uint32 cell-owner lane: 11*N_particles +
4*N_cells bytes at requested size. It is phase-local and reused across the
ownership and ID checks; capacity is reported with checked arithmetic. The
bounded diagnostic overload sorts copied IDs; the original public no-argument
hash validator remains unchanged for production callers. No new persistent
scientific state is introduced.

The FFT model includes 24*N_mesh^3 bytes at real/complex coexistence, plus a
separate reduction peak with complex mesh, worker lines, per-worker/global bin
accumulators and result/wrapper arrays. At 256^3 the two meshes require
402,653,184 bytes (384 MiB); at 512^3, 3,221,225,472 bytes (3 GiB).
These are analytical lower bounds on the whole calculation, not RSS results.
The bundle adds the ownership workspace, configured SF-history/quicklook arrays
and an 8 MiB metadata allowance. Opaque FFTW plan/allocator/OS memory remains
under the existing external-runtime reserve; the allowance is not a measured
upper bound for every external allocation. Production capacity reconciliation
and the existing hard ceiling remain authoritative. No claim of a certified
512^3 full-physics workstation envelope follows from these calculations.

The lease lifetime includes generation and local serialization. The public
engine retains its ungoverned compatibility constructor; production workflow
construction injects the existing governor. A caller-supplied enclosing lease
must be committed and sufficient; the caller owns its lifetime. Independent
standalone estimator calls acquire their own lease when the engine is governed.

## Compatibility, numerical and reproducibility impact

Existing runtime_error catches continue to accept MemoryAdmissionError, which
is a more specific subtype for hard-limit refusal only. Other errors are not
converted into optional deferrals. New public APIs are additive; existing
ungoverned calls and simulation-state validators remain source-compatible.
No input schema, normalized configuration key, restart serialization, science
output schema, filename convention, precision, force criterion, CFL rule,
refinement criterion, source equation, or scientific tolerance is changed.
Optional cadence remains coalesced/nonpersistent and products are labeled by
the actual execution epoch. Rejection decisions and diagnostic scheduling can
change with available headroom; the solver's scientific update order is not
changed.

## Validation and unresolved gates

The focused CPU build and four tests are recorded in the bundle acceptance
summary. The analysis unit regression checks mesh-size arithmetic, finite
headroom rejection below Red, exact-boundary recovery, lease release, duplicate
ID detection, ownership-validator equivalence, and numerical equality within
the existing 1e-10 small-mesh reference tolerance. The end-to-end test covers
transient deferral, catch-up, repeated missed cadence, and terminal drops.

The required remaining acceptance is: complete owner-specific peak/lifetime
contracts, full CPU/HDF5 inventories and source-package test, dependency-complete
MPI+FFTW np2/np3/np4 skew/empty/fault-injection/restart matrix, measured rank
max/mean RSS and work imbalance, and representative Release rank/thread
benchmarks. No unavailable runtime evidence is inferred from source review.

Suggested branch: campaign-m2d-final-admission-closure
Suggested PR title: Complete full-physics memory admission and certify M2D scheduling
