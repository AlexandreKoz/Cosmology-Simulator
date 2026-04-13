# TreePM coupling

## Scope

This module defines the explicit force split policy used by CosmoSim TreePM:

- PM computes long-range forces with Fourier filter `exp(-k^2 r_s^2)`.
- Tree computes the complementary short-range residual with the Gaussian real-space factor
  `erfc(r/(2 r_s)) + 2/sqrt(pi) * (r/(2 r_s)) * exp(-(r/(2 r_s))^2)`.

The split scale `r_s` is carried in `TreePmSplitPolicy::split_scale_comoving` and is also propagated into
`PmSolveOptions::tree_pm_split_scale_comoving` so the PM filter is auditable.

## Ownership and accumulation order

1. `TreePmCoordinator` zeros a compact active-set force view.
2. PM density assignment, Poisson solve, and interpolation produce long-range acceleration for active particles.
3. Tree build/traversal contributes short-range residual acceleration into the same active-set sidecars.
4. Diagnostics report continuity and composition checks near `r_s`.

## Assumptions

- Periodic wrapping in the residual pass uses minimum-image offsets and the PM box size.
- The residual tree multipole acceptance follows the existing geometric opening criterion.
- `TreeSofteningPolicy` remains authoritative for softened pair behavior.

## Provenance implications

- The split scale and kernel are explicit runtime configuration in `TreePmOptions`.
- Profiling fields expose PM phase timings plus tree short-range and coupling overhead timing.


## Periodic coupling validation ladder

- `tests/integration/test_tree_pm_coupling_periodic.cpp` is the periodic TreePM stop/go test.
- Failure diagnostics now include:
  - Build assumption (`COSMOSIM_ENABLE_FFTW`),
  - Runtime support gate (`treePmSupportedByBuild()`),
  - Force error (`rel_l2`) versus direct periodic minimum-image reference,
  - Split continuity metrics (`composition_error_at_split`, `max_relative_composition_error`).
- Tolerance policy:
  - FFTW-enabled path uses stricter `rel_l2 < 0.75`.
  - Non-FFTW fallback path is allowed a looser envelope (`rel_l2 < 1.8`) and is documented as non-production validation grade.
- Composition continuity remains strict (`< 1e-12`) to protect split-kernel correctness independent of backend performance.
