# Hydrodynamics core solver (finite-volume Godunov scaffold)

## Variable definitions

The hydro core uses comoving conserved variables per cell:

- `mass_density_comoving = rho`
- `momentum_density_*_comoving = rho * u_*_peculiar`
- `total_energy_density_comoving = rho * (e_internal + 0.5 |u_peculiar|^2)`

The explicit update implemented in `HydroCoreSolver::advancePatch` is:

\[
U^{n+1}_i = U^n_i - \frac{\Delta t}{a V_i} \sum_f A_f \hat{F}_f + \Delta t\, S_i.
\]

with `a = scale_factor`, `V_i = cell_volume_comoving`, and `A_f = face.area_comoving`.

## Stage structure and ownership

The update is split into stages to keep AMR reuse and future GPU kernels straightforward:

1. Reconstruction (`HydroReconstruction`): produce left/right primitive states per face.
2. Riemann (`HydroRiemannSolver`): compute face-centered numerical fluxes.
3. Flux accumulation (`HydroCoreSolver`): apply conservative owner/neighbor updates.
4. Source application (`HydroSourceTerm`): apply gravity/expansion and optional future sources.

Implementation is now split along these boundaries:

- `src/hydro/hydro_reconstruction.cpp` + `include/cosmosim/hydro/hydro_reconstruction.hpp`:
  limiter library, piecewise-constant reconstruction, and MUSCL-Hancock predictor/floors.
- `src/hydro/hydro_riemann.cpp` + `include/cosmosim/hydro/hydro_riemann.hpp`:
  HLLE and HLLC flux construction with HLL fallback on HLLC degeneracy.
- `src/hydro/hydro_core_solver.cpp` + `include/cosmosim/hydro/hydro_core_solver.hpp`:
  patch update scaffold (validation, cache fill, conservative flux accumulation, source-term pass, profiling).

Hydro hot data remains SoA in `HydroConservedStateSoa`. Cold patch metadata is isolated in `HydroPatchColdData`.
For hierarchical timestepping, `HydroActiveSetView` allows explicit active cell/face subsets via
`advancePatchActiveSet` without rewriting geometry ownership.
For memory-sensitive loops, `HydroScratchBuffers` and `HydroPrimitiveCacheSoa` can be reused across calls via
`advancePatchWithScratch` / `advancePatchActiveSetWithScratch` to avoid repeated allocations and repeated
primitive reconstruction for piecewise-constant paths.

## Source-term convention

`ComovingGravityExpansionSource` supplies conservative-source hooks using
`gravity_accel_*_peculiar` and `hubble_rate_code`:

- `S_rho = 0`
- `S_m = rho g - H m`
- `S_E = rho (u · g) - 2 H e_internal`

This is intentionally conservative and modular; additional baryonic sources should be implemented through extra `HydroSourceTerm` instances.

## Assumptions

- Face normals are unit vectors.
- Current patch geometry uses a uniform `cell_volume_comoving` per patch.
- `advancePatch` updates the full patch; `advancePatchActiveSet` updates only caller-selected cell/face subsets.
- Boundary handling is left to reconstruction/ghost filling; if no neighbor is present, piecewise-constant reconstruction mirrors owner state.
- The provided integration test uses a periodic 1D Sod-like setup to exercise the full stage path while preserving conservation.
