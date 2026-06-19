# AMR Hydro Validation Ladder

Date: 2026-06-17

## Level 0: CI guards

Current AMR hydro shock tube, Sedov, synchronization, reflux, subcycling, and pending-register tests are fast CI/regression guards. They assert production AMR coverage, finite state, positive density/pressure/internal energy where required, stable gas-cell identity coverage, explicit patch geometry, conservation diagnostics, and safe rejection of invalid reflux records.

These tests are necessary but not sufficient for scientific validation.

## Level 1: deterministic regression tests

The restart-equivalence tests compare direct continuation against checkpoint/reload/continue execution for deterministic scenarios. The AMR flux-register restart test exercises pending deferred reflux state across HDF5 restart.

This level protects restart and state-ownership truth. It does not measure physical convergence.

## Level 2: convergence-capable studies

Future AMR hydro validation should run shock tube, Sedov, and synchronization cases at multiple resolutions and collect structured error metrics. Default CI can keep tiny resolutions, but a non-CI validation profile should run enough resolution to test monotonic error behavior and conservation trends.

Suggested metrics:

- shock tube: line-sampled density, velocity, pressure, and L1/L2 errors against an analytic or trusted reference solution;
- Sedov: blast radius, radial profile monotonicity, energy localization, and total-energy conservation;
- synchronization: mass/momentum/energy conservation across coarse/fine interfaces and restart boundaries.

## Level 3: cross-code comparison

No cross-code validation is currently complete. Future comparisons should use published configurations against established AMR or moving-mesh solvers and should record initial conditions, boundary conditions, units, resolution, reconstruction/limiter/solver settings, and tolerances.

## Level 4: production science readiness

Production science readiness requires a documented validation deck, convergence plots, restart/reproducibility evidence, MPI execution evidence if MPI is claimed, and explicit limitations. The current repository has not reached this level for AMR hydro.

## Level 2 update: campaign machinery is implemented, evidence is not yet complete

The repository now contains versioned Sod/Sedov input decks, a common-profile metric implementation,
a command-line analyzer (`tools/run_amr_hydro_validation.py`), and opt-in campaign template. The fast C++
contract tests verify metric algebra and compatibility rejection. They do **not** execute a multi-resolution
production CHUÍ Sod/Sedov campaign and must not be cited as publication-grade convergence results.

## Level 3 update: external data contract exists, no external result exists

A versioned profile/manifest contract and common-grid comparison tool are implemented. The only in-tree
cross-code test is a deliberately synthetic contract fixture. No RAMSES/ENZO/AREPO/Athena++ output was
run or compared. This remains framework readiness, not cross-code validation.
