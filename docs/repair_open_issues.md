# Repair open issues (P01–P19 freeze ledger)

_Date captured: 2026-04-07 (UTC)_

## Current blocker ledger after dependency-enabled validation

No open P19 blockers remain after the PM/FFTW stabilization repair and full preset revalidation.

## Recently closed items

| ID | Status | Area | Prior symptom | Closure evidence |
|---|---|---|---|---|
| P19-GATE-FFTW-TEST-001 | Closed | PM FFTW analytic mode response | `unit_pm_solver` failed at `tests/unit/test_pm_solver.cpp:82` (`cosine_similarity > 0.98`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"` now passes. |
| P19-GATE-FFTW-TEST-002 | Closed | TreePM periodic PM long-range coupling (FFTW path) | `integration_tree_pm_coupling_periodic` failed with `rel_l2=18129.9` (required `<= 0.75`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` now passes (`36/36`). |
| P19-ARCH-CORE-BOUNDARY-003 | Closed | Core dependency direction | `include/cosmosim/core/reference_workflow.hpp` assembled analysis/I/O/physics workflow concerns inside `core/`. | Workflow assembly moved to `workflows/` and `integration_core_dependency_direction` guard added to fail on forbidden upward includes. |
| P19-CONFIG-CONTRACT-004 | Closed | Frozen configuration contract typing | Policy keys (`gravity_solver`, `hydro_solver`, `coordinate_frame`, boundaries, feedback mode/variant) remained string contracts after parse. | `unit_config_parser`, `unit_simulation_mode`, and `integration_simulation_mode_toy_runs` now validate enum-backed freeze behavior with deterministic normalized text/hash and centralized key+alias registry handling. |

## Verified non-blocking evidence (this run)

- PM HDF5+FFTW path passes:
  - `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"`
  - `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` (`36/36` passed)
- HDF5 path passes:
  - `ctest --preset test-hdf5-debug --output-on-failure` (`36/36` passed)
- CPU-only path passes:
  - `ctest --preset test-cpu-debug --output-on-failure` (`36/36` passed)

## Outcome tags

- **CPU-only preserved**
- **HDF5 path preserved**
- **PM HDF5+FFTW path restored**
- **P19 blocker ledger cleared**
