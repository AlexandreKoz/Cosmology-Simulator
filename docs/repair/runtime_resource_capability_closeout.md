# Runtime resource capability closure

Date: 2026-07-17  
Suggested PR: `refactor-enforce-runtime-resource-capabilities-and-complete-descriptor-contracts`

## Scope

This repair closes the public `StepContext` escape found after the Campaign A
runtime decomposition. It does not implement hierarchical KDK, distributed IC
reading, hydro wake-up, elastic restart, asynchronous output, or new physics.

## Implemented contract

- Public stage views no longer contain or expose `core::StepContext`,
  `SimulationState`, `ownerContext()`, or a protected `stageContext()` helper.
- The dispatcher binds each task's descriptor-declared resources to the view
  only for that invocation and clears the grant with RAII even when the task
  throws.
- The source-private `RuntimeStageAccess` bridge verifies the task grant before
  a built-in owner can enter its existing rung-zero implementation.
- Descriptor registration rejects resource keys or access modes that exceed
  the selected typed view. For example, an analysis task cannot declare write
  access to particle positions.
- Negative compile tests prove that an external `AnalysisRuntime` subclass and
  an ordinary stage-view caller cannot recover the removed owner context.
- The repeated gas-cell identity/parent-ownership helpers are consolidated in
  `workflows/internal/gas_cell_ownership.*`.
- Repository hygiene now rejects nested Windows ADS sidecars and generated
  runtime-output directories; leaked `Zone.Identifier` files and validation
  outputs were removed from the source handoff.

## Trust boundary

Built-in gravity, hydro/AMR, source, analysis, drift, and output owners are
trusted repository implementations. Their source-private bridge still enters
existing `StepContext`-based numerical bodies after validating the descriptor
grant. This preserves current rung-zero numerical ordering without exposing a
global mutable state escape to public or third-party module code.

Future lane-by-lane kernel capabilities can replace that internal bridge
incrementally. The current patch does not claim dynamic tracing of every memory
access inside a trusted built-in owner.

## Descriptor scope

`RuntimeModuleDescriptor` is authoritative for production stage construction,
ordering, typed view selection, and task-scoped resource grants. Timestep
criteria, restart payload schemas, and migration payload schemas remain owned by
their current dedicated runtime services. They are not claimed as descriptor
extension points in this patch.

## Validation floor

Executed in the repair environment:

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --test-dir build/cpu-only-debug -R \
  '^(unit_runtime_resources|unit_runtime_module_registry|integration_runtime_descriptor_probe|integration_runtime_decomposition_structure|integration_runtime_resource_capability_compile)$' \
  --output-on-failure
```

Result: 5/5 passed.

```text
ctest --test-dir build/cpu-only-debug -R \
  '^(integration_reference_workflow|integration_reference_workflow_hydro_row_order|integration_restart_equivalence_reference_workflow_hydro)$' \
  --output-on-failure
```

Result: 3/3 passed.

MPI, FFTW, and HDF5 runtime validation were intentionally left for the user's
dependency-complete environment. No pass claim is made for those paths here.
