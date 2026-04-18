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

The periodic PM solve reuses persistent solver-owned spectral scratch buffers for:

- the copied potential spectrum used for mesh potential reconstruction,
- the temporary gradient spectrum used to recover `a_x`, `a_y`, and `a_z`.

This keeps the PM operator auditable while avoiding repeated per-solve heap allocation churn on the hot periodic solve path.

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
```
