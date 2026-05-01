# CMake and CI hardening pass — 2026-05-01

## Scope

This pass targets CMake configure fragility, CI timeout behavior, stale or overly broad CTest selections, and generated-artifact hygiene. It intentionally avoids numerical rewrites.

## Repairs

- Added a C++20 target helper that prefers `target_compile_features(cxx_std_20)` when the active CMake/compiler pair exposes the feature table, and falls back to per-target `CXX_STANDARD` properties otherwise. This prevents configure-time failure on valid but feature-table-fragile CI toolchains.
- Added project-level CTest timeout defaults:
  - ordinary tests: `COSMOSIM_CTEST_DEFAULT_TIMEOUT_SECONDS`
  - long validation/MPI/I/O/TreePM tests: `COSMOSIM_CTEST_LONG_TIMEOUT_SECONDS`
- Centralized CTest property application after all tests are registered:
  - source-root working directory for every test
  - labels for unit/integration/validation/python/MPI tests
  - processor accounting for two-rank and three-rank MPI tests
  - timeout classification for long tests
- Hardened `scripts/ci/run_preset_pipeline.sh`:
  - uses `cmake --fresh --preset` in CI to avoid stale cache state
  - bounds CTest execution with `--timeout`
  - enables `--output-on-failure`
  - anchors CI matrix test regexes to exact test names
  - makes build parallelism configurable through `COSMOSIM_CI_BUILD_PARALLEL`
- Hardened the optional MPI smoke workflow to use exact test selection and an explicit timeout.
- Extended repository hygiene to syntax-check the CI shell entrypoints.

## Verification performed

- `cmake --preset cpu-only-debug`
- `cmake --preset hdf5-debug`
- `ctest --preset test-cpu-debug -N`
- anchored CTest selection smoke check for exact test-name matching
- `bash -n` on CI shell scripts
- `bash ./scripts/ci/check_repo_hygiene.sh`

`pm-hdf5-fftw-debug` configure was also attempted in this container and correctly failed because FFTW3 development libraries are not installed in the environment. That is an environment dependency failure, not a repository wiring failure.
