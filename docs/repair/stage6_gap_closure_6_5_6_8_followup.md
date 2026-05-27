# Stage 6 gap closure follow-up: 6.5 through 6.8

This note records the focused repair pass that closed the remaining Stage 6 audit gaps after the final acceptance pass.

## Runtime narrow-view closure

- `reference_workflow.cpp` now builds an explicit `AdaptiveTimeStepCriteriaView` before timestep-candidate evaluation. The criteria loop consumes spans for particle velocities, species tags, softening, force outputs, gas thermodynamic lanes, gas-to-particle mapping, and black-hole sidecar lanes instead of reading broad `SimulationState` inside the loop.
- `halo_workflow` now exposes `HaloParticleView` and `FofHaloFinder::buildCatalogFromView`. The legacy `buildCatalog(SimulationState, ...)` API remains as a compatibility wrapper that builds the narrow read-only view and delegates.
- `DiagnosticsEngine::generateBundle` builds `DiagnosticsStateView` once and uses narrow diagnostics views for run health, star-formation history, angular momentum, quicklook maps, projections, and the provisional power-spectrum path.

## Live memory accounting closure

- Tree nodes, PM mesh storage, TreePM active buffers, and TreePM exchange buffers now expose live capacity accounting hooks.
- `TreePmCoordinator::memoryReport()` reports PM mesh storage, tree storage, active-set work buffers, MPI/exchange buffers, and explicit uncertainty notes for external FFT/allocator/library allocations.
- `mergeMemoryReports()` combines runtime memory reports without turning views into owned bytes.
- The reference workflow now refreshes profiler memory reports with both persistent `SimulationState` memory and transient workspace/TreePM runtime memory at startup and during scheduler steps.

## Tests and validation

The repair is guarded by the final Stage 6 acceptance test, memory accounting tests, source/diagnostics view tests, halo workflow tests, and the Stage 6 benchmark aggregate target.

Validated commands in this environment:

```sh
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure -I 1,49
ctest --preset test-cpu-debug --output-on-failure -I 50,82
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4
ctest --preset test-hdf5-debug --output-on-failure -R 'snapshot_hdf5|restart_checkpoint|provenance_roundtrip'
cmake --build --preset build-cpu-debug --target stage6_benchmarks -j4
```

`pm-hdf5-fftw-debug` remains dependency-blocked in this environment because FFTW3 development libraries are not installed.
