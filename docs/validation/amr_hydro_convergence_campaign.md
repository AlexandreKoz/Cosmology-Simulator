# AMR Hydro Convergence Campaign

Date: 2026-06-18

## Purpose and separation from CI

`integration_amr_hydro_shock_tube` and `integration_amr_hydro_sedov` remain fast CI guards. They are
not presented as convergence evidence. The reproducible campaign interface is separate:

- input decks: `validation/input_decks/amr_sod_temporal_subcycling_v1.json` and
  `validation/input_decks/amr_sedov_temporal_subcycling_v1.json`;
- campaign template: `validation/campaign/amr_hydro_convergence_campaign_v1.json`;
- analysis tool: `tools/run_amr_hydro_validation.py`;
- profile/metric contract tests: `validation_amr_hydro_convergence`.

The campaign template intentionally contains unresolved output/reference paths. Replace them only with
actual CHUÍ profile outputs and a provenance-tagged exact or vetted reference profile. Generated campaign
results belong under a run/output directory and are not repository artifacts.

## Shock tube

The Sod input deck fixes gamma, left/right primitive states, domain, outflow boundary conditions, final
time, CFL, reconstruction, limiter, Riemann solver, AMR refinement ratio, temporal ghost policy, and a
64/128/256 effective-resolution ladder. A campaign result must include L1 density/velocity/pressure errors,
L2 density error, conservation diagnostics, AMR coverage, temporal ghost/reflux counters, and observed
slopes. Shock and contact discontinuities make globally second-order L1 claims inappropriate; monotonic
error reduction and physically interpretable shock-limited rates are the expected evidence.

## Sedov

The Sedov deck specifies 3D geometry, gamma, ambient density, explosion energy/deposition radius, box,
boundary conditions, final time, hydro method, AMR controls, and 32/64/128 ladder. The campaign must use
a documented Sedov-Taylor self-similar profile or vetted table, not merely blast-radius scaling. Report
radial density/velocity/pressure profile error, blast radius, energy/mass/momentum drift, symmetry metrics,
refinement distribution, and temporal ghost/reflux diagnostics.

## What ran in this repository pass

Only the deterministic C++ metric/contract tests ran. They validate interpolation, common-grid norms,
observed-order calculation, metadata compatibility rejection, and synthetic contract fixtures. No
multi-resolution CHUÍ shock or Sedov production campaign was executed here; therefore this repository does
not claim publication-grade convergence evidence.
