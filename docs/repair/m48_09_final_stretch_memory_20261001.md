# M48-09 — final TreePM stretch-memory package

_Date: 2026-10-01_

## Scope

M48-09 removes the three remaining population/mesh-scale convenience owners from the certified homogeneous periodic DMO TreePM path without changing precision, assignment scheme, deconvolution, TreePM splitting, node softening envelopes, or the out-of-place FFT architecture.

The three closures are:

- M48-09A: uniform build-time source softening as one validated scalar;
- M48-09B: periodic physical density deposited directly into the active FFT plan's padded real allocation;
- M48-09C: periodic Poisson multiplication derived from spectral geometry and small per-axis metadata rather than a full scalar spectral lane.

Generic heterogeneous softening, standalone/isolated compact density, real-potential materialization, CUDA staging, and isolated-open PM remain supported by their existing explicit fallback paths.

## M48-09A — uniform source softening

`TreeGravitySolver` now records source softening as either `kUniform` or `kMaterialized` and exposes a `ResolvedSourceSofteningView`. Uniform eligibility is proven structurally: the source per-particle epsilon lane and override mask must be absent, and no source species lane may activate species-dependent softening. The borrowed homogeneous-DMO production path therefore qualifies without an O(N) discovery scan.

Uniform mode validates and retains only `TreeGravityOptions::softening.epsilon_comoving` and releases any capacity left by a previous heterogeneous build. Heterogeneous configurations continue to resolve and store one exact epsilon per source. Tree traversal and TreePM residual evaluation hoist the representation selection out of their pair-interaction loops by selecting a base pointer plus a zero-or-one stride.

Tree node `softening_min_comoving` and `softening_max_comoving` remain authoritative traversal state. Uniform leaves initialize both envelopes directly from the scalar; materialized leaves retain the previous exact min/max aggregation. Legacy callers without a valid `GravitySourceGeneration` still validate current source content and source-softening identity against the build-time representation.

The `tree.source_softening` memory report now has zero vector capacity in uniform mode. The gravity preflight model removes the source-sized `double` term only when uniformity is an explicit runtime/preflight property; the borrowed homogeneous-DMO representation qualifies automatically.

## M48-09B — FFT-resident periodic density

The ordinary periodic CPU TreePM refresh now prepares the authoritative PM FFT plan before deposition and writes density directly into `PlanResources::real`. The real allocation remains the allocation against which FFTW plans were created and is never resized after plan construction.

The density destination uses the backend's explicit `real_z_stride`. Physical mesh cells are addressed as `(local_ix, iy, iz)` within padded rows; FFTW padding is never treated as a physical density cell. The complete real allocation is cleared before deposition, but accumulation, cell-volume normalization, density sums, mean subtraction, and all physical-cell diagnostics iterate physical `iz < nz` locations only.

`PmGridStorage::m_density` is now lazy. Non-const `density()` is the compatibility materialization boundary; const reads do not allocate. Generic `assignDensity()`, isolated/open PM, and the existing CUDA staging path therefore retain their compact-density semantics. The periodic TreePM specialization releases any stale compact-density capacity before direct deposition, so the specialized grid reports zero `pm_mesh.density` capacity.

A small plan/grid generation token records which grid and plan currently own valid physical density. `solvePoissonPeriodic*()` accepts either that exact current FFT-backed generation or an explicitly materialized compact compatibility density. A stale direct-density token fails closed. Immediately before the forward FFT, the token is consumed because `PlanResources::real` ceases to contain authoritative physical density once transform execution begins.

The historical compact-to-padded full-mesh copy is absent from the specialized production path. It remains only as the explicit generic compatibility boundary when a caller deliberately materializes compact density and then requests a periodic solve.

## M48-09C — derived periodic Poisson operator

`PlanResources` no longer owns a population-scale periodic Poisson-kernel vector. It retains only:

- `real`;
- `fourier`;
- `potential_k`;
- FFTW plan/layout metadata;
- O(nx + ny + nz/2) wave-number and assignment-window axis tables.

The axis tables cache signed `kx`, signed `ky`, non-negative r2c `kz`, and the one-dimensional assignment-window factors. For every logical spectral coefficient, the solver derives the same operator used previously:

`(-4*pi*G/k^2) * window_deconvolution * TreePM_Gaussian_long_range_filter`.

The multiplication preserves the prior `wx * wy * wz` transfer-window composition, the `1/max(transfer_window^2, 1e-12)` correction, the exact Gaussian split helper, and an exact zero DC mode. Split scale and gravitational constant are no longer cache-key inputs because they are applied directly on every solve rather than stored in cached full-mesh state.

Poisson multiplication and gradient construction now share one spectral traversal for both slab order `(local_ix, iy, iz)` and FFTW transposed order `(local_iy, ix, iz)`. `PlanResources` records logical and allocated local complex extents separately. Logical modes are processed from authoritative decomposition metadata, and any FFTW over-allocation tail is explicitly zeroed so the previous fail-closed tail behavior is preserved.

`spectral_operator_rebuilds` remains a compatibility telemetry counter, but now counts rebuilds of the small spectral axis/window metadata only.

## Source-derived memory model

For the M48 target of 512^3 particles, a 512^3 periodic PM mesh, and eight MPI ranks, the current source model is:

| Owner | Before M48-09 | After M48-09 specialized path | Source-derived change |
| --- | ---: | ---: | ---: |
| resolved tree source epsilon | 1,073,741,824 B aggregate | one scalar per rank; no population lane | -1.000000 GiB plus negligible scalar ownership |
| compact periodic PM density | 1,073,741,824 B aggregate | 0 B compact capacity; density aliases no new owner and lives in already-accounted `plan.real` | -1.000000 GiB |
| periodic Poisson scalar operator | 538,968,064 B aggregate for 67,371,008 logical r2c coefficients | 163,968 B aggregate for six axis tables across eight ranks | -0.501800 GiB net |

The gross removed population/mesh-scale owners are 2,686,451,712 B = 2.501953 GiB aggregate. After charging the six small axis tables (20,496 B/rank; 163,968 B across eight ranks), the source-derived net architectural reduction is 2,686,287,744 B = 2.501800 GiB aggregate. This is essentially the historical ~2.502-GiB stretch target, but it is recomputed from the current representation widths rather than forced by a subtraction constant.

The Poisson figure above uses the logical 512 x 512 x 257 r2c spectral extent. Backend-specific FFTW over-allocation can make the old full-kernel owner slightly larger on a particular communicator; runtime qualification must use the active backend allocation query and measured process memory. No measured RSS claim is made here.

The dominant periodic plan ownership after M48-09 is therefore the padded real allocation, the Fourier working spectrum, and `potential_k`, plus negligible spectral-axis metadata. The PM grid separately owns the three real-space force fields. Density is charged once, through the FFT real allocation.

## Reproducibility expectations

Tree softening changes representation only: the same fallback/species/override resolver remains authoritative for heterogeneous states, and uniform DMO uses the exact scalar that previously filled every source element. Pair-softening arithmetic and node-envelope values are intended to remain unchanged.

Periodic density changes address calculation, not the mathematical deposition. The same contribution records, routing order, stencil weights, physical-cell accumulation order, cell-volume normalization, and mean subtraction are retained. Padding is explicitly excluded from physical reductions.

The Poisson operator is mathematically identical and preserves the prior expression structure as closely as practical. Axis values and window factors are precomputed once rather than recomputed inside each spectral cell, so M48-Q must still establish the required bitwise/tolerance-level equivalence across supported compilers and FFTW layouts.

## Validation status

Tests: NOT RUN — intentionally deferred to M48-Q qualification.

Build: NOT RUN — intentionally deferred by campaign instruction.

Static source review only is part of this campaign. M48-Q remains responsible for numerical comparison, MPI/rank-count equivalence, FFTW slab/transposed coverage, whole-run high-water telemetry, and the 512^3 / 48-GiB qualification.
