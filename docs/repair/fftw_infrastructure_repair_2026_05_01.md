# FFTW infrastructure repair — 2026-05-01

## Scope

Dedicated FFTW3 discovery, CI dependency, and PM planner hardening pass. This repair focuses on repository infrastructure that can be validated without a local FFTW3 installation, while preserving the existing PM/TreePM solver contract.

## Repaired issues

- The `mpi-hdf5-fftw-debug` CI row installed `libfftw3-dev` but not `libfftw3-mpi-dev`. That made the distributed PM feature path depend on an unavailable `fftw3_mpi` component even on Ubuntu runners. The row and the infrastructure-gate dependency bundle now install `libfftw3-mpi-dev` explicitly.
- FFTW discovery now distinguishes serial FFTW from MPI FFTW. Serial PM presets require only `fftw3`; distributed MPI+FFTW presets call the finder with `REQUIRE_MPI` and fail at configure time with a precise package/install hint if `fftw3_mpi` is absent.
- FFTW discovery now supports three routes: CMake package config targets, CosmoSim's `FindFFTW3.cmake` module, and pkg-config imported targets for both `fftw3` and `fftw3_mpi`.
- The custom `FindFFTW3.cmake` module now probes `fftw3-mpi.h`, exposes `FFTW3::fftw3_mpi` only when both the MPI library and header are present, and attaches the serial FFTW target as a link dependency.
- FFTW plan creation no longer hard-codes `FFTW_MEASURE`. The build exposes `COSMOSIM_FFTW_PLAN_RIGOR` with validated values `ESTIMATE`, `MEASURE`, `PATIENT`, and `EXHAUSTIVE`. The default is `ESTIMATE` to avoid CI configure/test stalls from expensive planner searches.
- Repository hygiene now includes a static guard that prevents the MPI+FFTW CI path from regressing to a serial-only FFTW dependency bundle.

## Validation performed in this environment

- `bash -n scripts/ci/check_repo_hygiene.sh`
- `bash scripts/ci/check_repo_hygiene.sh`
- `cmake --preset cpu-only-debug`
- `cmake --preset pm-hdf5-fftw-debug` was checked to fail cleanly in this container with an explicit missing `libfftw3-dev` diagnostic, as expected because FFTW3 is not installed here.
- A fake-FFTW configure smoke test was used to verify that `FindFFTW3.cmake` and `cosmosim_find_fftw()` accept an explicit `FFTW3_ROOT` and generate build metadata with `COSMOSIM_ENABLE_FFTW=ON`.
- `src/gravity/pm_solver.cpp` was compiled in CPU/no-FFTW mode.
- `src/gravity/pm_solver.cpp` was compiled in serial FFTW mode against a fake header-only FFTW surface to validate the planner-macro and include/target wiring. This was compile-only and did not execute FFTW kernels.

## Not claimed

This environment still lacks a real FFTW3 installation and real `fftw3_mpi`, so the PM+HDF5+FFTW and MPI+HDF5+FFTW runtime test suites were not executed against the real FFTW libraries here. The repaired CI dependency surface is intended to make those paths executable on Ubuntu runners or local machines with the documented packages installed.
