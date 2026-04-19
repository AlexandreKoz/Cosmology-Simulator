# PM Gravity Solver (Periodic Cosmological Modes)

## Scope

This document defines the **operational contract** for `cosmosim::gravity::PmSolver` in the periodic cosmological path.

- Boundary condition: periodic box only for this stage.
- Assignment/interpolation kernel: runtime-selectable `CIC` or `TSC`, with matched deposit/gather semantics.
- Backend policy:
  - `COSMOSIM_ENABLE_FFTW=ON`: FFTW-backed production path.
  - `COSMOSIM_ENABLE_FFTW=OFF`: fallback `naive_dft` for bring-up/small tests, not production-grade TreePM.

## Mathematical contract

Input mesh field to `solvePoissonPeriodic()` is a physical density field `ρ(x)` in solver code units.

The solver applies explicit mean subtraction:

- `δρ(x) = ρ(x) - ρ̄`
- `ρ̄ = (1/Ncell) Σcell ρ(xcell)`

and solves the periodic comoving Poisson equation:

- `∇² φ(x) = 4 π G a² δρ(x)`

with Fourier-space relation (for `k != 0`):

- `φ_k = - 4 π G a² δρ_k / k²`

and force/acceleration relation:

- `a_i(k) = - i k_i φ_k`

The periodic zero mode is pinned by policy:

- `φ_{k=0} = 0`
- `a_{k=0} = 0`

This enforces the periodic-box gauge convention (potential mean fixed to zero; only differences are physical).

## Discrete Fourier conventions

For real-to-complex (`r2c`) storage with mesh `(Nx, Ny, Nz)`, the reduced-complex shape is `(Nx, Ny, Nz/2 + 1)`.

Mode mapping used by the solver:

- `kx(ix) = 2π/L * (ix <= Nx/2 ? ix : ix - Nx)`
- `ky(iy) = 2π/L * (iy <= Ny/2 ? iy : iy - Ny)`
- `kz(iz) = 2π/L * iz`, `iz ∈ [0, Nz/2]`

The negative `kz` branch is represented by Hermitian conjugates and is not explicitly stored.

## Normalization and backend behavior

- FFTW inverse transforms are unnormalized by FFTW; solver applies `1/Ncell` after each inverse transform.
- Fallback `naive_dft` path already applies inverse normalization internally; no extra factor is applied on that path.
- The contract above (`δρ`, `φ_k`, `-ikφ_k`, zero mode) is backend-invariant.

## API-level output guarantees

After `solvePoissonPeriodic(grid, options, ...)` returns successfully:

- `grid.potential()` contains the periodic potential solve `φ(x)` with zero mode pinned.
- `grid.force_x()`, `grid.force_y()`, `grid.force_z()` contain mesh acceleration components from `a_i(k) = -i k_i φ_k`.

Potential is a supported output, not an incidental side effect.

For particle-space sampling:

- `interpolateForces(...)` gathers mesh acceleration to particles using the same assignment kernel selected for deposition.
- `interpolatePotential(...)` gathers mesh potential to particles using the same geometry and stencil.

Matched deposition + gather is a hard contract for both schemes in this stage.

## Assignment/gather kernels (Phase 1)

Let `u = x / Δ` be the particle coordinate in mesh-cell units (per axis).

- **CIC** (first-order B-spline):
  - support: 2 cells per axis
  - weights around `i = floor(u)`:
    - `w_i = 1 - (u - floor(u))`
    - `w_{i+1} = u - floor(u)`
  - Fourier window per axis: `W_CIC(k_i) = sinc(k_i Δ / 2)^2`

- **TSC** (second-order B-spline):
  - support: 3 cells per axis
  - with `j = floor(u + 1/2)` and `δ = u - j`:
    - `w_{j-1} = 0.5 * (0.5 - δ)^2`
    - `w_j = 0.75 - δ^2`
    - `w_{j+1} = 0.5 * (0.5 + δ)^2`
  - Fourier window per axis: `W_TSC(k_i) = sinc(k_i Δ / 2)^3`

## Memory and scratch behavior

### PM slab ownership/storage contract (Phase 2 milestone)

PM mesh ownership is now explicit through a centralized slab layout contract:

- global shape is `(Nx, Ny, Nz)`,
- decomposition is by contiguous half-open x-slabs `[x_begin, x_end)`,
- each global real-space cell `(ix, iy, iz)` is owned by exactly one rank: the slab owner of `ix`.

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

- each rank owns only `layout.local_nx * Ny * Nz` real cells and solves only that slab portion;
- no rank-0 gather path is used in the long-range PM solve;
- FFT plans are created with `fftw_mpi_plan_dft_r2c_3d/c2r_3d` and `MPI_COMM_WORLD`;
- one-rank solves and multi-rank solves share the same k-space Poisson/gradient operator.

For this phase, FFT spectral storage is treated as **non-transposed slab output**:

- local spectral index: `(local_ix, iy, iz)` with flat index `(local_ix * Ny + iy) * (Nz/2+1) + iz`;
- global mode index: `ix = owned_x.begin + local_ix`, `iy`, `iz`;
- mode mapping remains:
  - `kx(ix) = 2π/L * (ix <= Nx/2 ? ix : ix - Nx)`
  - `ky(iy) = 2π/L * (iy <= Ny/2 ? iy : iy - Ny)`
  - `kz(iz) = 2π/L * iz`.

Distributed density assignment is now enabled for slab layouts in `PmSolver::assignDensity`
when MPI is enabled (`COSMOSIM_ENABLE_MPI=ON`) and slab metadata matches
`MPI_COMM_WORLD`.

`PmSolver::interpolateForces` and `PmSolver::interpolatePotential` now support slab-distributed
particle ownership with explicit reverse communication. Particle-owner ranks build global
stencil requests and send them to slab owners (`MPI_Alltoallv`); slab owners validate ownership,
compute weighted contributions from local PM fields, and return per-particle contributions to
the requesting owner ranks (`MPI_Alltoallv`). Final acceleration/potential accumulation occurs
only on the particle-owner rank in local particle order.



### Distributed interpolation reverse-message contract (Phase 2)

Ownership and message flow for both force and potential gather:

1. **Particle-owner rank** computes CIC/TSC stencil nodes in global periodic mesh coordinates.
2. For each stencil node `(ix, iy, iz)`, particle owner routes a request to
   `pmOwnerRankForGlobalX(Nx, world_size, ix)`.
3. Request payload fields are:
   - `particle_index` (owner-local index into the caller spans),
   - `global_ix/global_iy/global_iz`,
   - `weight` (matched deposit/gather kernel weight).
4. **Slab-owner rank** receives requests, validates that the x-index is locally owned and the
   global indices are in range, then computes:
   - force gather: `weight * (ax, ay, az)` from owner-local PM force fields,
   - potential gather: `weight * phi` from owner-local PM potential field.
5. Slab owners send per-request contributions back to the originating particle-owner rank.
6. Particle-owner rank accumulates returned contributions by `particle_index` into output spans.

Ordering/determinism policy:

- Request generation order is deterministic: particle index order, then stencil loop order
  (`x`, `y`, `z`).
- Returned contributions are accumulated on owner rank in MPI receive order; accumulation target is
  owner-local particle order indexed by `particle_index`.
- No rank requires replicated full PM fields for interpolation in distributed mode.

### Distributed density assignment message contract (Phase 2)

Ownership and routing model:

- **particle owner rank** computes assignment stencils in global mesh coordinates from
  wrapped periodic particle positions.
- **slab owner rank** is chosen by x-index ownership:
  - `destination_rank = pmOwnerRankForGlobalX(Nx, world_size, global_ix)`.
- Only slab owners accumulate to local PM density storage.
- No remote direct writes are permitted.

Record format (packed per contribution):

- `global_ix` (`uint32`): global x cell index in `[0, Nx)`.
- `global_iy` (`uint32`): global y cell index in `[0, Ny)`.
- `global_iz` (`uint32`): global z cell index in `[0, Nz)`.
- `mass_contribution` (`double`): already weighted by the matched assignment kernel.

Batching and ordering:

- Each owner rank batches records by destination rank.
- In-batch order is deterministic append order from nested loops:
  particle index order, then stencil axis loop order (`x`, `y`, `z`).
- Batched records are exchanged via `MPI_Alltoallv` as byte payloads.
- PM solver-owned send/recv buffers are reused across solves for stable layout metadata.

Receiver validation before accumulation:

- Reject any received record with out-of-range global indices.
- Reject any record whose `global_ix` is not owned by the receiving slab rank.
- Accepted records are accumulated only into owner-local slab storage and then normalized by local cell volume.

The periodic PM solve reuses persistent solver-owned spectral scratch buffers for:

- the copied potential spectrum used for mesh potential reconstruction,
- the temporary gradient spectrum used to recover `a_x`, `a_y`, and `a_z`.

This keeps the PM operator auditable while avoiding repeated per-solve heap allocation churn on the hot periodic solve path.

Plan/scratch caches are keyed by slab layout ownership metadata (`world_size`, `world_rank`,
`owned_x.begin`, `owned_x.end`) and are reused until layout/communicator metadata changes.
Inverse normalization is applied exactly once per inverse field (`φ`, `a_x`, `a_y`, `a_z`).

## Long-range field cadence and reuse (Reference workflow Phase 1)

The reference TreePM workflow now supports explicit long-range PM reuse controlled by
`numerics.treepm_update_cadence_steps`:

- cadence unit: gravity kick opportunities (`gravity_kick_pre`, `gravity_kick_post`)
- refresh rule: rebuild PM mesh/field every `N` kick opportunities (`N >= 1`)
- reuse rule: between refreshes, reuse the cached PM long-range field and only rebuild/evaluate short-range tree residual forces

Phase 1 scope is deliberately conservative and single-rank. It is an auditable cadence control surface, not a full multirate production design.

Operational metadata records which PM field version each kick used, including the field build step index and build scale factor.

## Optional modifiers

- `enable_window_deconvolution=true` applies scheme-aware deconvolution to the **combined particle-transfer operator** (`deposit * gather`) in k-space:
  - CIC: divide by `(W_CIC(k_x) W_CIC(k_y) W_CIC(k_z))^2`
  - TSC: divide by `(W_TSC(k_x) W_TSC(k_y) W_TSC(k_z))^2`
  - safeguard floor: denominator is clamped to `>= 1e-12` before division
- `tree_pm_split_scale_comoving > 0` applies TreePM long-range Gaussian filter in k-space.
  - In the reference workflow and TreePM coordinator contract, this value is not an independent UX knob;
    it is derived from normalized config controls as `r_s = asmth_cells * Δmesh`.
  - `Δmesh = box_size / PMGRID` (Phase-1 cubic PM UX).
  - Companion TreePM cutoff uses `r_cut = rcut_cells * Δmesh` on the residual tree path.

These modifiers do not alter the base sign/normalization contract above.

Default policy in this phase is conservative:

- `assignment_scheme = CIC`
- `enable_window_deconvolution = false`

## Accuracy/cost tradeoffs in this stage

- CIC is cheaper (2-point stencil/axis, 8 points in 3D) and remains first-class.
- TSC is smoother and generally reduces anisotropy/self-force artifacts, but uses 3 points/axis (27 points in 3D).
- Deconvolution improves transfer amplitude matching for resolved modes, but can amplify high-`k` noise/aliasing near Nyquist; therefore it is opt-in.

## Implementation note

- Explicit CUDA PM assignment/gather remains CIC-only in this build; CPU paths support both `cic` and `tsc`. This guard is intentional and explicit rather than hidden.

## Validation focus for this stage

Validation and tests explicitly cover:

- analytic single-mode potential and force shape,
- uniform-density cancellation (mean subtraction + zero-mode policy),
- potential-force consistency for a simple periodic mode,
- transverse leakage diagnostics in periodic mode integration test.

Recommended commands:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug -R "unit_pm_solver|integration_pm_periodic_mode|validation_integration"
# when MPI is available:
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug -R "integration_pm_periodic_mode_mpi_two_rank"
```
