# Reference workflow integration pass

This document defines the integrated reference workflow that wires configuration parsing, generated IC ingest, integrator stage sequencing, diagnostics callbacks, star-formation and black-hole callbacks, and profiler report generation through one auditable entry point (`workflows::ReferenceWorkflowRunner`).

## Layer placement and boundary

- The concrete assembly lives in `include/cosmosim/workflows/reference_workflow.hpp` and `src/workflows/reference_workflow.cpp`.
- `core/` provides only typed state/config/stage contracts (`StepOrchestrator`, `IntegratorState`, etc.) and does not assemble analysis/I/O/physics callbacks.
- Transitional compatibility aliases are exposed in `namespace cosmosim::core` from the workflow header so existing callers can migrate incrementally without pulling workflow assembly back into `core/`.

## Scope and conservative assumptions

- The reference path requires `numerics.gravity_solver=treepm`, `numerics.hydro_solver=godunov_fv`, and `units.coordinate_frame=comoving`.
- Schema compatibility is checked against snapshot schema v1 and restart schema v2.
- By default, integration tests run with `write_outputs=false` so the path remains valid when HDF5 is disabled.
- When HDF5 is enabled, the runner can emit and read back restart and snapshot files in one pass.

## Canonical stage ordering

The runner uses `core::StepOrchestrator` and therefore preserves the canonical order:

1. `gravity_kick_pre`
2. `drift`
3. `hydro_update`
4. `source_terms`
5. `gravity_kick_post`
6. `analysis_hooks`
7. `output_check`

## Profiling and reporting

Each run emits:

- `reference_profile.json`
- `reference_profile.csv`

These files are generated from the built-in `core::ProfilerSession` and provide an end-to-end stage profile suitable for regression gating and mini-run cost comparisons. The workflow also emits a `canonical_stage_order` flag in reports/bench output, backed by the canonical order validator in `core::time_integration`.
