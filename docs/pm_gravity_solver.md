# PM Gravity Solver (Periodic Cosmological Modes)

## Scope

This document defines the **operational contract** for `cosmosim::gravity::PmSolver` in the periodic cosmological path.

- Boundary condition: periodic box only for this stage.
- Assignment/interpolation kernel: runtime-selectable `CIC` or `TSC`, with matched deposit/gather semantics.
- Backend policy:
  - `COSMOSIM_ENABLE_FFTW=ON`: FFTW-backed production path.
  - `COSMOSIM_ENABLE_FFTW=OFF`: fallback `naive_dft` for bring-up/small tests, not production-grade TreePM.

## Mathematical contract

The input mesh field to `solvePoissonPeriodic()` is comoving mass density
`rho_com(x)` in code units: owner particle mass divided by comoving code
volume. It is not a physical-density lane and must not be multiplied by
`a^3` before the solve.

The solver applies explicit mean subtraction:

- `delta_rho_com(x) = rho_com(x) - mean(rho_com)`
- `mean(rho_com) = (1/Ncell) sum_cell rho_com(xcell)`

and solves for the scale-free comoving Newtonian potential kernel `Psi`:

- `nabla_x^2 Psi(x) = 4 pi G_code delta_rho_com(x)`

with Fourier-space relation (for `k != 0`):

- `Psi_k = -4 pi G_code delta_rho_com,k / k^2`

and force/acceleration relation:

- `A_i(k) = -i k_i Psi_k`

There is deliberately no factor of `a^2` in this operator. The production
workflow converts the physical Newton constant into configured code units with
`core::newtonGravitationalConstantCode(UnitSystem)` and supplies that same
`G_code` to PM and the complementary short-range tree. In a cosmological run,
the returned `A` is the common scale-free input to the dynamical update. For
physical peculiar velocity `u = a dx/dt`, both collisionless KDK and the gas
conservative source implement

```text
du/dt + H u = A/a^2.
```

The particle KDK factors and `ComovingGravityExpansionSource` therefore own
the scale-factor response for their respective state, outside the PM kernel.
`PmSolveOptions::scale_factor` remains required positive validity/source-time
metadata, but changing it must not rescale an otherwise identical PM kernel.

The periodic zero mode is pinned by policy:

- `Psi_{k=0} = 0`
- `A_{k=0} = 0`

This enforces the periodic-box gauge convention (potential mean fixed to zero; only differences are physical).

## Discrete Fourier conventions

For real-to-complex (`r2c`) storage with mesh `(Nx, Ny, Nz)`, the reduced-complex shape is `(Nx, Ny, Nz/2 + 1)`.

Mode mapping used by the solver:

- `kx(ix) = 2π/Lx * (ix <= Nx/2 ? ix : ix - Nx)`
- `ky(iy) = 2π/Ly * (iy <= Ny/2 ? iy : iy - Ny)`
- `kz(iz) = 2π/Lz * iz`, `iz ∈ [0, Nz/2]`

The negative `kz` branch is represented by Hermitian conjugates and is not explicitly stored.

## Normalization and backend behavior

- FFTW inverse transforms are unnormalized by FFTW; solver applies `1/Ncell` after each inverse transform.
- Fallback `naive_dft` path already applies inverse normalization internally; no extra factor is applied on that path.
- The contract above (`delta_rho_com`, `Psi_k`, `-ik Psi_k`, zero mode) is
  backend-invariant.

## API-level output guarantees

After `solvePoissonPeriodic(grid, options, ...)` returns successfully:

- `grid.potential()` contains `Psi(x)` with the zero mode pinned.
- `grid.force_x()`, `grid.force_y()`, `grid.force_z()` contain the scale-free
  comoving kernel components from `A_i(k) = -i k_i Psi_k`.

Potential is a supported output, not an incidental side effect.

For particle-space sampling:

- `interpolateForces(...)` gathers mesh acceleration to particles using the same assignment kernel selected for deposition.
- `interpolatePotential(...)` gathers mesh potential to particles using the same geometry and stencil.

Matched deposition + gather is a hard contract for both schemes in this stage.

## Distributed reverse-interpolation integrity contract

For slab-distributed periodic PM force gather, particle rows remain authoritative
only on the particle-owner/origin rank. The DMO force path routes compact
**particle/x-plane group** requests rather than one retained record per mesh
cell. Each request carries:

- `origin_rank` and `destination_rank`;
- an opaque origin-local particle token;
- a deterministic request sequence and monotonic 64-bit exchange epoch;
- up to the three TSC x-plane indices/weights owned by the same destination;
- the wrapped y/z mesh coordinates needed for the destination to expand the
  y/z stencil locally.

The slab owner validates sender/destination identity, epoch, sequence, x-plane
ownership, bounds, and finite numerical fields before expanding the y/z stencil.
It accumulates all mesh contributions represented by that request and returns
one force response. The origin regenerates the deterministic expected request
stream for response validation; it does not retain an O(N) request registry.
Wrong sender/token, duplicate/out-of-order, stale, malformed, or non-finite
responses fail closed before authoritative particle output is mutated.

Communication proceeds in coordinated rounds bounded by
`treepm_pm_exchange_batch_bytes` **per peer**. Local x-plane work is consumed
directly without entering MPI buffers. Send and receive wire buffers are reused
between request and response phases, and the profiler records their actual
high-water marks. Zero-work ranks still enter every collective round.

Potential interpolation retains the older per-cell request/response protocol
for now; it is not used for the first-light DMO force update. Both protocols use
checked byte-count/displacement conversion and coordinated preflight/failure
agreement.

## Assignment/gather kernels (Phase 1)

Let `u_i = x_i / Δi` be the particle coordinate in mesh-cell units on axis `i`, where
`Δx = Lx/Nx`, `Δy = Ly/Ny`, `Δz = Lz/Nz`.

- **CIC** (first-order B-spline):
  - support: 2 cells per axis
  - weights around `i = floor(u)`:
    - `w_i = 1 - (u - floor(u))`
    - `w_{i+1} = u - floor(u)`
  - Fourier window per axis: `W_CIC(k_i) = sinc(k_i Δi / 2)^2`

- **TSC** (second-order B-spline):
  - support: 3 cells per axis
  - with `j = floor(u + 1/2)` and `δ = u - j`:
    - `w_{j-1} = 0.5 * (0.5 - δ)^2`
    - `w_j = 0.75 - δ^2`
    - `w_{j+1} = 0.5 * (0.5 + δ)^2`
  - Fourier window per axis: `W_TSC(k_i) = sinc(k_i Δi / 2)^3`

## Memory and scratch behavior

### PM slab ownership/storage contract (Phase 2 milestone)

PM mesh ownership is now explicit through a centralized slab layout contract:

- global shape is `(Nx, Ny, Nz)`,
- decomposition is by contiguous half-open x-slabs `[x_begin, x_end)`,
- each global real-space cell `(ix, iy, iz)` is owned by exactly one rank: the slab owner of `ix`.

The x-slab authority map intentionally matches FFTW-MPI's default input-block
layout. With `B = ceil(Nx / world_size)`, rank `r` owns

```text
[min(r*B, Nx), min((r+1)*B, Nx))
```

so the final block is truncated and high ranks may own zero planes. This is not
the alternative "give the first `Nx % world_size` ranks one extra plane"
partition. Matching FFTW's layout is required: a different software-side owner
map can configure at np2 yet fail ownership validation at np3 or np4.

Storage is represented by `parallel::PmSlabLayout` and consumed by
`gravity::PmGridStorage`:

- `PmGridStorage(shape)` is the single-rank degenerate case (`world_size=1`, owned slab `[0, Nx)`),
- `PmGridStorage(shape, layout)` allocates only `layout.local_nx * Ny * Nz` cells,
- global/local mapping is contractually centralized in:
  - `local_x = global_x - owned_x.begin` (valid iff `global_x ∈ [owned_x.begin, owned_x.end)`),
  - `global_x = owned_x.begin + local_x` (valid iff `local_x < owned_x.extent`),
  - local linear index: `(local_x * Ny + iy) * Nz + iz`.

`PmSolver::solvePoissonPeriodic` now supports true slab-distributed FFT on MPI ranks when
`COSMOSIM_ENABLE_MPI=ON` and `COSMOSIM_ENABLE_FFTW=ON`:

- The configured `parallel::PmSlabLayout` is validated against FFTW MPI ownership (`local_nx`, `local_0_start`) before plan creation; mismatched user-side slab partitions are rejected instead of silently assuming backend compatibility.

- each rank owns only `layout.local_nx * Ny * Nz` real cells and solves only that slab portion;
- no rank-0 gather path is used in the long-range PM solve;
- FFT plans are created with `fftw_mpi_plan_dft_r2c_3d/c2r_3d` and `MPI_COMM_WORLD`;
- one-rank solves and multi-rank solves share the same k-space Poisson/gradient operator.

For this phase, FFT spectral storage is treated as **non-transposed slab output**:

- local spectral index: `(local_ix, iy, iz)` with flat index `(local_ix * Ny + iy) * (Nz/2+1) + iz`;
- global mode index: `ix = owned_x.begin + local_ix`, `iy`, `iz`;
- mode mapping remains:
  - `kx(ix) = 2π/Lx * (ix <= Nx/2 ? ix : ix - Nx)`
  - `ky(iy) = 2π/Ly * (iy <= Ny/2 ? iy : iy - Ny)`
  - `kz(iz) = 2π/Lz * iz`.

Distributed density assignment is now enabled for slab layouts in `PmSolver::assignDensity`
when MPI is enabled (`COSMOSIM_ENABLE_MPI=ON`) and slab metadata matches
`MPI_COMM_WORLD`.

`PmSolver::interpolateForces` and `PmSolver::interpolatePotential` now support slab-distributed
particle ownership with explicit reverse communication. Particle-owner ranks build global
stencil requests and send them to slab owners (`MPI_Alltoallv`); slab owners validate ownership,
compute weighted contributions from local PM fields, and return per-particle contributions to
the requesting owner ranks (`MPI_Alltoallv`). Final acceleration/potential accumulation occurs
only on the particle-owner rank in local particle order.

Assignment, periodic solve, and force/potential gather use the actual
`MPI_COMM_WORLD` size to enter a coordinated API preflight before any PM
collective or layout-based branch. The preflight fingerprints the API
operation, full-domain-serial versus distributed mode, mesh `nx/ny/nz`, box
axes, scale factor, `G_code`, assignment/deconvolution, execution/residency,
decomposition, boundary, split/workspace controls, and entry-specific control
lanes; disagreement fails on every rank. Validation and payload-decoding
failures are likewise reduced before throwing, so one bad rank cannot leave
peers blocked in the next count or payload exchange. A full-domain serial
reference solve is legal inside an MPI job only when every rank enters the same
PM call and selects that mode. Rank-zero-only PM calls are not legal collective
usage. Distributed callers must construct the slab layout from the active
communicator on every rank. In an MPI-enabled executable that has not called
`MPI_Init`, or has already finalized, the library observes a serial world and
does not call `MPI_Comm_size`/`MPI_Comm_rank`.

### Distributed force-interpolation reverse-message contract

The production DMO force gather uses the bounded plane-group protocol:

1. **Particle-owner rank** computes CIC/TSC x support and groups x planes by
   `pmOwnerRankForGlobalX(Nx, world_size, ix)`.
2. Locally owned x-plane groups are accumulated directly from local PM fields.
3. A remote group is encoded as one fixed-width force-plane request containing
   the x indices/weights plus y/z mesh coordinates. The receiver expands y/z
   weights locally, so TSC does not require 27 retained request records.
4. **Slab-owner rank** validates the record, evaluates all represented local
   stencil cells, and emits one accumulated force response for that request.
5. The origin validates the response against the deterministic request stream
   and accumulates it into the origin-owned target row.

The same deterministic particle/group ordering is used in every round. The
number of rounds is derived collectively from the maximum active/source count,
so empty ranks remain legal participants without dummy particles.

### Distributed PM wire contract

PM routing never transmits native C++ object layouts. The version-1 wire format
is fixed-width, little-endian, and requires IEEE-754 binary64. Every record
starts with the 12-byte header `(magic="PMW1", version=1, record_kind)`. The
first-light DMO routing adds compact plane-group kinds alongside the legacy
per-cell kinds used by potential gather:

- density plane-group request: 96 bytes;
- force plane-group request: 96 bytes;
- force response: 64 bytes;
- legacy density contribution: 56 bytes;
- legacy force request: 56 bytes;
- potential request: 56 bytes;
- potential response: 48 bytes.

Decoders require exact size, magic, version, kind, reserved-lane validity,
identity consistency, and complete record alignment before mutating mesh or
particle fields.

### Distributed density assignment message contract

The density path uses the same x-plane grouping strategy:

- particle owners compute wrapped CIC/TSC x support and group planes by slab
  destination;
- locally owned groups deposit directly into owner-local density storage;
- one compact remote record carries up to three x planes plus y/z coordinates
  and particle mass;
- the destination expands y/z weights and deposits directly into its slab;
- communication is split into deterministic coordinated rounds with at most
  `treepm_pm_exchange_batch_bytes` of payload per destination peer in a round.

For TSC this bounds retained routing state by the configured communication
policy rather than by `27 * N_local`. The path does not simultaneously retain
a nested C++ record population, a flattened C++ copy, a wire copy, and a
decoded receive population. Two reusable physical wire buffers cover send and
receive phases. Profile counters expose routed logical records, participating
peers, wire bytes, measured MPI wait, and actual send/receive buffer high-water
bytes.

Receiver validation rejects out-of-range or non-owned x indices, wrong
sender/origin/destination identity, stale epochs, invalid sequence, or non-finite
coordinates/mass before accumulation. Accepted mass is accumulated only into
owner-local slab storage and is normalized by local cell volume.

The periodic PM solve reuses persistent solver-owned spectral scratch buffers for:

- the copied potential spectrum used for mesh potential reconstruction,
- the temporary gradient spectrum used to recover `a_x`, `a_y`, and `a_z`.

This keeps the PM operator auditable while avoiding repeated per-solve heap allocation churn on the hot periodic solve path.

Plan/scratch caches are keyed by slab layout ownership metadata (`world_size`, `world_rank`,
`owned_x.begin`, `owned_x.end`) and are reused until layout/communicator metadata changes.
Inverse normalization is applied exactly once per inverse field (`φ`, `a_x`, `a_y`, `a_z`).
FFTW-MPI ranks with a legal zero-width slab retain a logical extent of zero but
receive a one-element dummy allocation for backend pointer safety; they still
participate in plan creation, transforms, and all PM collectives.

### PM force-slab halo exchange

TreePM populates a one-dimensional x-halo cache for PM force interpolation.
Logical left/right peers are derived from the owners of the adjacent global x
planes, not from `rank-1`/`rank+1`; this skips zero-width slabs and handles a
truncated final FFTW block. A rank with no owned planes participates in the
call but has no peers, payload, or halo storage. Locally owned periodic sides
are copied without MPI.

Remote sides post all `MPI_Irecv` operations before `MPI_Isend`, then complete
with `MPI_Waitall`. Side- and sequence-qualified tags keep left/right
orientation distinct even when both logical neighbors are the same peer in a
two-rank periodic layout. The exchanged depth is bounded by the smallest
nonempty slab extent, and payload/count/byte arithmetic is checked before MPI.
All preparation, allocation, and local field-shape failures are reduced before
point-to-point traffic. A global fingerprint binds mesh shape, requested halo
depth, periodic mode, communicator size, and exchange sequence, so divergent
rank-local controls fail coherently. The single-rank helper remains legal in an
MPI-enabled library without an active MPI session.

`integration_pm_slab_halo_exchange_mpi_{two,three,four}_rank` covers ordinary
periodic rings, same-peer left/right exchange, `Nx < world_size` zero-width
slabs, and the single-global-plane local-copy case. All three rank counts passed
in the 2026-07-13 MPI+HDF5+FFTW run.

### Transposed-slab distributed PM path and future pencil contract

`numerics.treepm_pm_decomposition_mode = pencil` currently activates FFTW's
transposed spectral ownership over x-slab real-space ownership when FFTW+MPI
distributed plans are available. It is an alternate transposed-slab backend,
not a true 2-D process-grid pencil decomposition:

- real-space ownership remains x-slab (`PmSlabLayout`) for particle-owner deposition/interpolation,
- forward FFT uses `fftw_mpi_plan_dft_r2c_3d(..., FFTW_MPI_TRANSPOSED_OUT)` and switches spectral ownership to y-partitioned transposed storage,
- inverse FFT uses `FFTW_MPI_TRANSPOSED_IN` to return from transposed spectral ownership back to x-slab real-space ownership.

The solver-facing `PmDecompositionDescriptor` is now topology-neutral and can
represent serial, x-slab, transposed-slab spectral ownership, and a future true
2-D pencil process grid without changing TreePM force ownership. No in-house FFT
implementation is fabricated by this refactor.

## PM backend capability tiers

Backend existence is separate from production capability. The non-FFTW direct
DFT is `DiagnosticNaiveDft`; FFTW is `ProductionFftw`. The current CUDA path may
perform device assignment/interpolation but still uses a host FFT, so it does
not claim the reserved `ProductionCudaFft` tier. Production TreePM never silently
falls back to the diagnostic DFT; only the explicit
`treepm_allow_diagnostic_naive_dft=true` test/diagnostic override permits it.
`PmBackendArchitecture` records host/device ownership so a future persistent,
device-resident cuFFT backend can be added without relabeling today's staging
path.

## Long-range field cadence and reuse

Production configuration requires `numerics.treepm_update_cadence_steps = 1`
and `numerics.hierarchical_max_rung = 0`. The latter is a fail-closed
restriction: production `ReferenceWorkflow` does not yet persist the
per-element kick/drift epochs needed for a correct mixed-rung KDK update. Every
integrator-issued, rank-coordinated production force-refresh surface therefore
rebuilds PM. Cadence values greater than one are rejected because no validated
long-range predictor or temporal interpolator exists.

The integrator-owned synchronization event carries kick opportunity, refresh
decision, field version, last refresh opportunity, build step, and build scale
factor. The workflow requires all ranks to agree on those values and on the
refresh vote before entering PM collectives. Empty active sets do not opt a
rank out of this consensus.

The coordinator retains an explicit lower-level reuse API for controlled tests
and future integrator work. Reuse requires an exact transient signature match
for force epoch, force-evaluation scale factor, `G_code`, split scale,
rectangular box axes, assignment, boundary, decomposition mode, and
deconvolution policy on every rank. A missing or incompatible cache makes a
reuse request fail on all ranks; it is not silently converted into an
unrequested solve. Divergent explicit refresh votes likewise fail before
density assignment or FFT collectives.

Particle ownership changes advance the TreePM decomposition epoch used by tree
wire records, but they do not by themselves invalidate a PM mesh field owned by
fixed FFT slabs. The long-range compatibility comparison therefore deliberately
does not include the particle-decomposition epoch. It does include the physical
force epoch and every mesh/operator input above. Dense-row acceleration history
is separate: the workflow invalidates that force cache immediately after a
committed ownership change. Production currently does not exercise local-bin
PM reuse because `hierarchical_max_rung=0`; lower-level reuse evidence must not
be described as mixed-rung production maturity. Operational metadata still
records each solve/reuse decision and its exact field-version/build context.

## Optional modifiers

- `enable_window_deconvolution=true` applies scheme-aware deconvolution to the **combined particle-transfer operator** (`deposit * gather`) in k-space:
  - CIC: divide by `(W_CIC(k_x) W_CIC(k_y) W_CIC(k_z))^2`
  - TSC: divide by `(W_TSC(k_x) W_TSC(k_y) W_TSC(k_z))^2`
  - safeguard floor: denominator is clamped to `>= 1e-12` before division
- `tree_pm_split_scale_comoving > 0` applies TreePM long-range Gaussian filter in k-space.
  - In the reference workflow and TreePM coordinator contract, this value is not an independent UX knob;
    it is derived from normalized config controls as `r_s = asmth_cells * Δmesh`.
  - `Δmesh = cbrt((Lx/Nx)*(Ly/Ny)*(Lz/Nz))` in this stage.
  - Companion TreePM cutoff uses `r_cut = rcut_cells * Δmesh` on the residual tree path.
  - The certified default is `asmth_cells=1.25`, `rcut_cells=6.25`, hence
    `r_cut/r_s=5`. Smaller positive explicit cutoffs remain accepted for
    compatibility/sweeps but are not covered by the default certification.
  - Periodic TreePM additionally requires `r_cut < min(Lx,Ly,Lz)/2` because
    the residual tree evaluates one minimum image per source. Typed config and
    direct coordinator entry fail closed at equality or above; refine the PM
    mesh or lower `rcut_cells` and revalidate the resulting profile.

The TreePM coordinator applies the same physical `G_code`, with no scale-factor
power, to the complementary real-space residual. It rejects unequal PM/tree
`G_code` inputs and derives the effective short-range prefactor without
changing standalone-tree semantics. A non-unit-scale-factor regression checks
that an otherwise identical complete PM-plus-tree kernel is scale invariant;
the particle KDK kick or gas conservative source supplies the time-dependent
`1/a^2` response.

These modifiers do not alter the base sign/normalization contract above.

Default production policy is the independently validated profile:

- `assignment_scheme = TSC`
- `enable_window_deconvolution = true`

This default change is reproducibility-relevant and is captured by normalized
configuration and `provenance_v7`. Existing explicit CIC/deconvolution-off
decks keep their requested behavior.

## Accuracy/cost tradeoffs in this stage

- CIC is cheaper (2-point stencil/axis, 8 points in 3D) and remains available
  for compatibility and diagnostics, but is not the accuracy-certified default
  for the current mesh/split envelope.
- TSC uses 3 points/axis (27 points in 3D) and is the certified CPU/FFTW profile.
- Deconvolution improves resolved-mode transfer matching but can amplify
  high-`k` noise/aliasing near Nyquist. It is enabled for the certified default;
  users changing assignment, deconvolution, mesh, split, or cutoff must treat
  that as a new accuracy profile.
- Raising the cutoff from the historical 4.5-cell profile to 6.25 cells closes
  the observed cutoff-transition accuracy gap, but increases residual tree
  search volume, pair work, and distributed target/response traffic. That
  tradeoff must be included in profiling and scaling evidence.

## Implementation note

- Explicit CUDA PM assignment/gather remains CIC-only in this build; CPU paths
  support both `cic` and `tsc`. The CUDA 12.0 `cuda-debug` build and its
  123/123 test inventory pass, including the PM smoke. Pre-sm_60 targets use a
  CAS-backed double density atomic while newer targets use native
  double-precision `atomicAdd`. This build/runtime closure does not promote
  CIC into the Ewald-certified TSC profile.

## Validation focus for this stage

Validation and tests explicitly cover:

- analytic single-mode potential and force shape,
- uniform-density cancellation (mean subtraction + zero-mode policy),
- potential-force consistency for a simple periodic mode,
- transverse leakage diagnostics in periodic mode integration test,
- FFTW-MPI slab ownership and halo exchange at np2/np3/np4, including empty
  slabs and same-peer periodic neighbors,
- an independent Ewald comparison for the complete periodic TreePM force.

The Ewald gate certifies the current TSC+deconvolution profile with the typed
runtime `opening_theta=0.7` default at `relative_L2 <= 1e-2` and
`p99_normalized <= 5e-2`. The 2026-07-13 run observed worst certified values
`8.059667017e-3` and `8.157551429e-3` on the exact cancellation-dominated
`4^3` DMO initial lattice. CIC is still measured, but
its worst covered uniform and split-transition cases exceed those targets and
remain diagnostic. See `docs/tree_pm_coupling.md` for the complete fixture and
reference contract.

Recommended commands:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug -R "unit_pm_solver|integration_pm_periodic_mode|validation_periodic_ewald_reference|validation_tree_pm_ewald_accuracy"
# when MPI is available:
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure -R "integration_pm_periodic_mode_mpi|integration_pm_slab_halo_exchange_mpi|validation_phase2_mpi_gravity"
```

- Distributed PM source terms are owner-local inputs: ranks contribute only their owned particle subset, while PM deposition and interpolation route remote slab interactions explicitly through MPI.


## Isolated/open PM operator (non-periodic)

For `mode_policy.gravity_boundary = isolated_monopole`, the PM long-range solve
uses a doubled-domain free-space convolution on the mesh. Single-rank execution
solves locally; the bounded multi-rank compatibility route gathers the slab
field to a root solve and scatters the result:

1. embed comoving code density on a padded grid `(2Nx, 2Ny, 2Nz)`,
2. convolve with the free-space kernel in Fourier space,
3. extract the physical `(Nx,Ny,Nz)` block.

This is linear (non-circular) convolution on the physical domain, so isolated mode does not reuse periodic-image PM semantics.

We solve the isolated Poisson problem in a finite domain with far-field gauge
\(\phi(\|x\|\to\infty)=0\):

- unsplit: \(\nabla_x^2\Psi = 4\pi G_{\rm code}\rho_{\rm com}\)
- TreePM long-range split: \(\Psi_{\mathrm{LR}}(x)= -G_{\rm code}\int \rho_{\rm com}(x')\,\frac{\mathrm{erf}(\|x-x'\|/(2r_s))}{\|x-x'\|}\,d^3x'\)

with \(r_s =\) `tree_pm_split_scale_comoving`.

Conventions in this stage:
- **Potential gauge:** \(\Psi(\infty)=0\), approximated by finite-box open convolution.
- **Force recovery:** \(\mathbf A=-\nabla_x\Psi\), using second-order central differences in the interior and second-order one-sided stencils on physical boundaries.
- **Boundary condition model:** open/non-periodic; no periodic image wrapping in the short-range residual path.
- **Assignment/gather boundary stencils:** physical-domain CIC/TSC stencil
  points are clipped, not wrapped to the opposite face. Particles outside the
  open domain contribute/sample zero. Clipped weights are not renormalized into
  the domain, preserving the finite-domain convolution interpretation.
- **Split consistency:** PM supplies the unsoftened Gaussian long-range field. The tree evaluates the residual `F_softened,total - F_Newtonian,long`, so Tree+PM composes to the configured Plummer-softened total force. At zero softening this reduces exactly to the usual complementary Gaussian short-range factor before explicit cutoff.
- **Self term policy:** kernel value at `r=0` is set to zero in the PM convolution.

Current stage limitations and safety policy:
- Typed config validation requires
  `numerics.treepm_enable_window_deconvolution=false` for isolated/open
  gravity. The periodic-profile default is `true`, so an isolated deck that
  omits the key fails during config loading; the solver does not silently
  rewrite frozen configuration. Tracked isolated examples set it explicitly.
- Single-rank isolated PM remains the normal supported path.
- Multi-rank isolated/open PM uses a root-gather/root-scatter implementation only
  as a bounded small-grid compatibility path. It is correct inside the configured
  envelope but is not a scalable distributed isolated solver.
- The root-owned transient workspace is estimated before `MPI_Gatherv`, padded
  convolution workspace allocation, or field `MPI_Scatterv`. The guard names the
  grid shape, rank count, estimated root bytes, configured limit, selected route
  (`root_gather_scatter`), and policy (`bounded_small_grid_only`).
- Configure the envelope with
  `parallel.isolated_pm_root_workspace_limit_bytes`; exceeding it fails on all
  ranks before dangerous root allocation/communication.
- Focused zoom long-range correction currently gathers high-resolution source
  coordinates/masses across ranks for the correction solve. It is similarly
  bounded by `parallel.zoom_high_res_allgather_limit_bytes` and reports the
  gathered bytes/limit in TreePM diagnostics; it must not be described as a
  production-scale source transport.
- The focused correction uses the isolated/open convolution and explicitly
  disables PM window deconvolution even when the global periodic solve uses the
  certified TSC+deconvolution profile. Isolated/open deconvolution is rejected
  and has no accuracy certification in this stage.
