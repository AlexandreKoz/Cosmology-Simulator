# PM Gravity Solver (Periodic Cosmological Modes)

## Scope

This document defines the **operational contract** for `cosmosim::gravity::PmSolver` in the periodic cosmological path.

- Boundary condition: periodic box only for this stage.
- Assignment/interpolation kernel: CIC only.
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

- `interpolateForces(...)` gathers mesh acceleration to particles using CIC transpose gather.
- `interpolatePotential(...)` gathers mesh potential to particles using the same CIC geometry/convention.

## Memory and scratch behavior

The periodic PM solve reuses persistent solver-owned spectral scratch buffers for:

- the copied potential spectrum used for mesh potential reconstruction,
- the temporary gradient spectrum used to recover `a_x`, `a_y`, and `a_z`.

This keeps the PM operator auditable while avoiding repeated per-solve heap allocation churn on the hot periodic solve path.

## Optional modifiers

- `enable_window_deconvolution=true` applies CIC-window deconvolution in k-space (`1/W(k)^2`), clamped to avoid division blow-up.
- `tree_pm_split_scale_comoving > 0` applies TreePM long-range Gaussian filter in k-space.

These modifiers do not alter the base sign/normalization contract above.

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
