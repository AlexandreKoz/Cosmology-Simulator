# Stage 6.5-6.8 repair notes

This repair pass hardens the Stage 6 source, diagnostics, transient-state, and memory-reporting contracts without introducing a second particle or gas state owner.

## Scope

- Stage 6.5: source, feedback, AGN, tracer, and diagnostics code now exposes narrow runtime views for hot/update loops.
- Stage 6.6: persistent restart truth remains separated from transient workspaces; reconstruction gradients stay in transient workspace storage from the prior 6.0-6.4 hardening pass.
- Stage 6.7: memory accounting now reports owned capacity for the live `SimulationState`, module sidecars, species lists, patch metadata, and `TransientStepWorkspace` buffers.
- Stage 6.8: profiling JSON reports persistent/transient totals, subsystem entries, high-water bytes where available, unknown bytes, and uncertainty notes.

## Narrow-view contracts added

- `StarFormationRuntimeView` for active-cell star formation updates.
- `StellarFeedbackGeometryView` and `StellarFeedbackDepositionView` for feedback target selection and deposition.
- `BlackHoleAgnAccretionView` for sparse black-hole accretion and feedback updates.
- `StellarEvolutionRuntimeView` for active-star mass return, metals, and feedback bookkeeping.
- `TracerHostMassView` for tracer host-cell mass updates.
- `DiagnosticsStateView`, `ParticleDiagnosticsView`, `GasDiagnosticsView`, and `StarFormationHistoryView` for analysis reductions and quicklook products.

The legacy broad-state entry points remain as compatibility wrappers. They build narrow views and delegate to the view-based implementations.

## Memory-accounting contract

Memory reports distinguish:

- owned bytes vs referenced bytes;
- persistent bytes vs transient bytes;
- high-water marks when the owning subsystem exposes them;
- known owned host allocations vs unknown or external allocations.

Known unmeasured external allocations, such as third-party HDF5/FFTW internals, OS allocator overhead, and future device allocations, are reported as uncertainty notes rather than exact numbers.

## Validation

The repair added focused unit coverage for Stage 6.5-6.8 source/diagnostics view contracts and extended memory-accounting tests for persistent lanes, workspace capacity, required subsystem categories, and pre-run estimates.
