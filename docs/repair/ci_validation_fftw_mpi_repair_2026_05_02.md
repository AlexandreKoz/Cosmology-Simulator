# CI validation, FFTW, and MPI repair — 2026-05-02

## Scope

This repair addresses CI failures observed in the CPU-only, PM+HDF5+FFTW, and MPI+HDF5+FFTW debug lanes.

## Fixed failures

- CPU-only validation timeouts in `validation_integration` and `validation_phase2_mpi_gravity_single_rank` were caused by validation tests using production FFT-sized PM meshes while the build was intentionally using the naive DFT fallback. The tests now keep the FFTW-backed validation sizes for real FFTW builds and use reduced correctness-sized meshes for no-FFTW fallback builds.
- The FFTW-backed `unit_pm_solver` isolated-open test used a stale box-size assertion. The isolated solver uses zero-padded open-boundary convolution, so a self-similar cell-space point-mass setup should preserve relative discretization error under pure box scaling rather than necessarily improve. The test now checks that invariant.
- The MPI+FFTW build did not propagate MPI's include interface into `cosmosim_gravity`, even though `pm_solver.cpp` includes and uses MPI in distributed FFTW paths. `cosmosim_gravity` now links `MPI::MPI_CXX` directly when MPI+FFTW is enabled, and the local FFTW MPI imported target also propagates MPI when available.
- The source include order now includes `mpi.h` before `fftw3-mpi.h` in MPI+FFTW builds, matching FFTW MPI's header dependency.

## Validation performed in this environment

This container does not provide the real FFTW/MPI development stack used by GitHub Actions, so the FFTW runtime and MPI runtime lanes could not be executed here. The CPU/no-FFTW repaired path was built and the formerly failing tests were run:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug --parallel 2 --target \
  test_validation_integration test_validation_phase2_mpi_gravity test_unit_pm_solver
ctest --test-dir build/cpu-only-debug \
  -R '^(unit_pm_solver|validation_integration|validation_phase2_mpi_gravity_single_rank)$' \
  --output-on-failure --timeout 120
```

Result: all three selected tests passed in the local CPU/no-FFTW build.
