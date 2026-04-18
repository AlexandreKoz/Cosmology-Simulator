# PM Gravity Solver (Periodic Cosmological Modes)

## Scope and assumptions

- The implementation targets long-range periodic PM gravity and currently supports **CIC** assignment only.
- Runtime selection is config-driven through `numerics.treepm_*` keys (no hard-coded PM mesh/split in the workflow path).
- Force interpolation uses the transpose of the assignment kernel (CIC gather) to preserve consistency.
- FFT backend defaults to FFTW when `COSMOSIM_ENABLE_FFTW=ON`; otherwise a correctness-oriented fallback backend (`naive_dft`) is used only for small-test workflows and not treated as production TreePM support.
- The periodic zero mode is explicitly set to zero (`phi_0 = 0`).
- Units are solver-local and explicit: `box_size_mpc_comoving`, `scale_factor`, and `gravitational_constant_code` are required inputs.

## Solver stages

1. Density assignment (`assignDensity`): particle mass to mesh with CIC.
2. FFT and Poisson solve (`solvePoissonPeriodic`):
   - `rho_k = FFT[delta_rho]`
   - `phi_k = -4π G a^2 rho_k / k^2`, `k != 0`
   - Optional CIC window deconvolution (`1 / W(k)^2`).
3. Gradient in k-space:
   - `g_k = i k phi_k`
4. Inverse FFT to real-space force components.
5. Force interpolation (`interpolateForces`): CIC transpose gather to particles.

## Profiling hooks

`PmProfileEvent` records:

- `assign_ms`
- `fft_forward_ms`
- `poisson_ms`
- `gradient_ms`
- `fft_inverse_ms`
- `interpolate_ms`
- `bytes_moved`

These fields are intended for in-run profiling and benchmark summaries, not for correctness validation.

## FFT backend discovery and validation path

When `COSMOSIM_ENABLE_FFTW=ON`, CMake resolves FFTW in this order:

1. `find_package(FFTW3 CONFIG)`
2. `find_package(FFTW3 MODULE)` (`cmake/FindFFTW3.cmake`)
3. `pkg-config` fallback (`fftw3.pc`)

If all three fail, configuration stops with an actionable fatal error. The build does not silently downgrade a requested FFTW-enabled path.

Recommended PM/TreePM validation workflow:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug
```

## Build gating

TreePM build support can be checked with:

- `treePmSupportedByBuild()`
- `requireTreePmSupportOrThrow(cosmosim::core::GravitySolver::kTreePm)`

If TreePM is requested without FFTW support, the gate emits a clear runtime error message describing the required build flag (`COSMOSIM_ENABLE_FFTW=ON`).

## Extensibility notes

- `PmAssignmentScheme` already reserves `kTsc` for later extension.
- `numerics.treepm_assignment_scheme=tsc` is parsed and normalized now, but runtime intentionally rejects it until TSC deposition/interpolation support is completed.
- FFT backend selection is represented in API (`PmSolver::fftBackendName`) and can be adapted for future cuFFT/pluggable backend wiring.
- Current data ownership keeps PM grids (`PmGridStorage`) independent from particle state arrays.

## Phase-1 runtime config contract

The reference workflow derives TreePM split controls from normalized config using one formula path:

- `N = numerics.treepm_pm_grid`
- `Δmesh = cosmology.box_size / N`
- `r_s = numerics.treepm_asmth_cells * Δmesh`
- `r_cut = numerics.treepm_rcut_cells * Δmesh`

Additional runtime PM controls:

- `numerics.treepm_assignment_scheme`
- `numerics.treepm_enable_window_deconvolution`
- `numerics.treepm_update_cadence_steps`


## Integration stop/go validation ladder

- `tests/integration/test_pm_periodic_mode.cpp` is the periodic PM stop/go test.
- The test now reports backend, cosine similarity, transverse leakage, and build-flag assumptions in failure diagnostics.
- Tolerance policy is explicit by feature path:
  - `COSMOSIM_ENABLE_FFTW=ON`: strict periodic spectral-shape check (`cosine_similarity > 0.9`).
  - `COSMOSIM_ENABLE_FFTW=OFF`: fallback path is checked only for finite/non-zero response and documented as non-production.
- This keeps CPU-only fallback useful for bring-up while preventing it from being mistaken for the FFTW validation grade.
