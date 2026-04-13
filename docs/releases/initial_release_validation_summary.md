# Initial release validation summary (v0.1.0-initial)

This summary defines the minimum validation evidence required to declare release readiness.

## Validation gate definition

A release-candidate commit must pass:

1. Core unit tests (`ctest -L unit` or equivalent selected unit suite).
2. Integration smoke through the reference workflow pipeline.
3. Validation ladder checks:
   - unit tolerance checks,
   - integration tolerance checks,
   - regression checks,
   - convergence checks.

## Required reference tests

- `unit_config_parser`
- `unit_simulation_state`
- `unit_time_integration`
- `integration_reference_workflow`
- `integration_snapshot_hdf5_roundtrip` (when HDF5 is enabled)
- `integration_restart_checkpoint_roundtrip`
- `validation_unit`
- `validation_integration`
- `validation_regression`
- `validation_convergence`

## Scientific credibility statement

This release gate is based on explicit documented tolerance checks and integration/regression behavior, not on feature count. Experimental modules are documented separately and are not treated as validated production features.

## Limitations and assumptions

- Validation evidence in this document is scoped to repository tests and reference decks under `validation/`.
- GPU/CUDA path is not a release-science gate in this initial release.
- MPI path is included via decomposition/restart consistency tests but does not claim large-cluster scaling validation.
