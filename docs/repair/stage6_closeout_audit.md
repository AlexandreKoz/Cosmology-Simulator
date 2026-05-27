# Stage 6 closeout audit

This file records the repository-side closeout gate for Stage 6 after the 6.9/6.10 hardening pass.
It is intentionally evidence-oriented rather than a design essay.

## Final acceptance matrix

| Stage 6 requirement | Status | Evidence |
| --- | --- | --- |
| One primary particle data model | PASS | `SimulationState::particles` is the canonical `ParticleSoa` owner; `ParticleSoaStorage` is a non-restart utility and is rejected by compile-time guards. |
| Explicit owners for required particle lanes | PASS | Hot physics lanes live in `ParticleSoa`; IDs, species tags, rank/ownership, drift epoch, and softening override metadata live in `ParticleSidecar`; acceleration is transient force output only. |
| Other views are derived adapters | PASS | `GravityParticleKernelView`, `HydroCellKernelView`, diagnostics/source views, and active views are span-based derived views with generation guards where they scatter. |
| Drift/kick and active hot kernels use compact/narrow views | PASS | Active particle/cell views are materialized in `TransientStepWorkspace`; stale scatter is rejected by `test_stage6_final_acceptance`. |
| Gravity hot paths use narrow views | PASS | Tree active-set APIs and PM force output views consume position/mass/target/output spans rather than swollen particle objects. |
| Hydro/source/timestep/analysis views are narrow by default | PASS | Source, AGN, tracer, stellar feedback/evolution, and diagnostics modules expose narrow runtime views; compatibility wrappers delegate to those view contracts. |
| Persistent and transient state are separated | PASS | Hydro reconstruction gradients live in `TransientStepWorkspace`, not `GasCellSidecar`; restart payload root is `RestartPersistentStateView`. |
| Restart serializes persistent truth only | PASS | `RestartWritePayload` exposes no transient workspace, hydro scratch, PM workspace, MPI buffer, or output buffer member; hash tests prove transient workspace mutation does not change restart integrity. |
| Memory accounting covers required subsystems | PASS | Runtime and pre-run reports cover particles, gas/hydro, tree, PM mesh, sidecars, active sets, MPI buffers, scratch, and output buffers with capacity-based owned bytes or explicit unknown/estimate notes. |
| Memory reporting overhead is benchmarkable | PASS | `bench_memory_accounting_overhead` and `bench_active_kernel_compact_vs_full` cover accounting overhead; `stage6_benchmarks` builds the Stage 6 benchmark set. |
| Tests/benchmarks enforce Stage 6 behavior | PASS | `test_stage6_layout_contracts`, `test_stage6_source_diagnostics_views`, `test_stage6_final_acceptance`, `test_stage6_active_views`, and memory/restart tests are wired into CTest. |
| Clean repo | PASS | Build directories, generated binaries, caches, logs, and test outputs are excluded from release zip packaging. |

## Stage 6-specific validation targets

Recommended focused commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug --target test_unit_stage6_layout_contracts test_unit_stage6_source_diagnostics_views test_unit_stage6_final_acceptance test_integration_stage6_active_views test_unit_memory_accounting bench_active_kernel_compact_vs_full bench_gravity_kernels bench_hydro_kernels bench_memory_accounting_overhead bench_profiling_overhead -j4
ctest --test-dir build/cpu-only-debug --output-on-failure -R "stage6|memory_accounting|restart_checkpoint_schema|restart_checkpoint_roundtrip"
cmake --build --preset build-cpu-debug --target stage6_benchmarks -j4
```

Feature-dependent validation:

- HDF5 restart/snapshot tests require an HDF5-enabled preset and available HDF5 development libraries.
- FFTW/TreePM PM backend tests require FFTW development libraries.
- MPI/CUDA validation remains feature-gated and must be reported as `NOT RUN` when dependencies are unavailable.

## Remaining risks

- Benchmark numbers are machine- and build-type-dependent; they are regression instruments, not absolute performance claims.
- Tree, PM, MPI, HDF5, and GPU external allocator internals remain explicitly unknown unless subsystem-specific hooks report owned buffers.
- Compatibility wrappers that accept `SimulationState` remain for API stability; their contract is to build narrow views and delegate, not to become new hot-loop authority.
