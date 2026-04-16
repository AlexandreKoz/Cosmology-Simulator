# Initial release validation summary (v0.1.0-initial)

This summary defines the minimum validation evidence required to declare release readiness.

## Validation gate definition

A release-candidate commit must pass:

1. The actual unit and integration suite through `ctest --preset ... --output-on-failure` or an explicitly enumerated equivalent command.
2. Integration smoke through the real config-driven runtime application path.
3. Validation ladder checks:
   - unit tolerance checks,
   - integration tolerance checks,
   - regression checks,
   - convergence checks.

Release reports should prefer the preset-wide or explicitly enumerated command rather than a shorthand subset invocation.

## Required reference tests

- `unit_config_parser`
- `unit_simulation_state`
- `unit_time_integration`
- `unit_pm_solver`
- `integration_reference_workflow`
- `integration_runtime_app_smoke` (when HDF5 is enabled)
- `integration_snapshot_hdf5_roundtrip` (when HDF5 is enabled)
- `integration_restart_checkpoint_roundtrip`
- `validation_unit`
- `validation_integration`
- `validation_regression`
- `validation_convergence`

## Scientific credibility statement

This release gate is based on explicit documented tolerance checks, the real config-driven runtime smoke path, and integration/regression behavior, not on feature count. Experimental modules are documented separately and are not treated as validated production features.

## Limitations and assumptions

- Validation evidence in this document is scoped to repository tests and runnable release smoke configs.
- GPU/CUDA path is not a release-science gate in this initial release.
- MPI build support exists, but the default release gate does not claim mandatory multi-process production validation.
