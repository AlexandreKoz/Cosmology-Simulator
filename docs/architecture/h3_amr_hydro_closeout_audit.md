# H3 AMR Hydro Closeout Audit

Date: 2026-06-16
Prompt range: H3.0-H3.9
Mode: production repair closeout

## Executive verdict

The H3 AMR hydro branch is no longer invalidated by the previous `cosmosim_amr` compile failure. The build blocker in `applyFluxRegistersToSimulationState` was fixed without removing stale-target validation, and the CPU/no-HDF5/no-MPI AMR/hydro test bundle now compiles and passes in this environment.

This is still **not** a full distributed or restart-complete AMR hydro closeout. The accepted claim after this repair is narrower: the local/single-rank production AMR hydro path builds, is selected when production AMR coverage exists, uses `core::SimulationState` as authoritative hydro state, maps gas rows to patch-local cells from explicit patch geometry rather than sorted dense-row order, performs conservative refine/derefine through row-order-independent maps, fills ghost descriptors with auditable status metadata, emits solver-derived flux-register records, and skips unsafe reflux records instead of applying them.

## Environment and configuration used

Configured and tested in the local container with:

- Compiler: GNU C++ 14.2.0 (`c++ (Debian 14.2.0-19) 14.2.0`).
- Build type: `Release`.
- Tests: enabled.
- HDF5: disabled.
- MPI: disabled.
- FFTW/CUDA: disabled by the configured build.

Configuration command:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCOSMOSIM_ENABLE_TESTS=ON \
  -DCOSMOSIM_ENABLE_HDF5=OFF \
  -DCOSMOSIM_ENABLE_MPI=OFF
```

## Commands run and outcomes

Build commands:

```bash
cmake --build build --target cosmosim_amr -j2
cmake --build build --target test_integration_amr_production_hydro_integration -j2
./build/test_integration_amr_production_hydro_integration
cmake --build build --target integration_amr_hydro_shock_tube -j2
cmake --build build --target integration_amr_hydro_sedov -j2
cmake --build build --target integration_amr_synchronization_stress -j2
```

Outcome: all commands above completed successfully.

The full AMR/hydro CTest set was built first because this repository's default `ctest` execution expects test executables to exist. Then the targeted AMR/hydro test regex was run:

```bash
ctest --test-dir build --output-on-failure -R "amr|AMR|hydro"
```

Outcome: **24/24 tests passed**.

Passing tests in that run:

- `unit_hydro_reconstruction`
- `unit_hydro_boundary_conditions`
- `unit_hydro_riemann`
- `integration_amr_patch_migration`
- `unit_amr_refinement`
- `unit_amr_hydro_geometry`
- `unit_amr_ghost_fill`
- `unit_amr_flux_register_generation`
- `unit_hydro_core_solver`
- `integration_hydro_sod_like`
- `integration_hydro_decoupled_gas_cells`
- `integration_hydro_conservation_periodic`
- `integration_hydro_axis_symmetry`
- `integration_amr_static_refinement_sync`
- `integration_amr_patch_boundary_consistency`
- `integration_amr_conservative_refine`
- `integration_amr_conservative_derefine`
- `integration_amr_reflux_conservation`
- `integration_amr_production_hydro_integration`
- `integration_amr_hydro_shock_tube`
- `integration_amr_hydro_sedov`
- `integration_amr_synchronization_stress`
- `integration_restart_equivalence_hydro_toy`
- `validation_hydro_classics`

## H3.0-H3.9 verdicts

### H3.0 -- AMR/hydro coupling design audit

Verdict: **accepted for the local production ownership contract**.

`core::SimulationState` remains the production owner of gas/hydro state. `amr::PatchHierarchy` and `AmrPatch::m_conserved` remain scaffold/topology/test infrastructure. The workflow still enters the AMR-driven hydro path when `amr::hasProductionAmrHydroCoverage(context.state)` is true and does not fall back to the old global Cartesian toy patch in that case.

Limit: this does not close distributed AMR hydro, AMR subcycling, HDF5 AMR hydro restart equivalence, or galaxy-formation physics.

### H3.1 -- AMR patches as first-class hydro geometry providers

Verdict: **accepted for explicit local patch geometry and row-order-independent mapping**.

A shared helper in `include/cosmosim/amr/amr_patch_indexing.hpp` now computes patch-local `(i,j,k)` and linear patch-cell indices from `PatchDescriptor` origin, extent, cell dimensions, and gas-cell center coordinates. It validates positive dimensions/extents, bounds, duplicate cells, missing cells, and identity ownership. Both AMR hydro geometry construction and orchestrator-side validation use the same helper.

`buildAmrHydroPatchGeometry` now builds patch-local row order from explicit geometry and stable identity records. It no longer treats sorted dense rows as patch-local truth.

Limit: production coverage still assumes a patch's dense rows are represented as a complete patch block in the current `SimulationState`/identity model. The physical mapping within that block is row-order independent, but arbitrary interleaving of one patch's rows through another patch's dense block is still outside the current production coverage contract.

### H3.2 -- Conservative prolongation on refine

Verdict: **accepted for piecewise-constant conservative production refinement**.

`refineProductionPatchInSimulationState` now obtains parent rows through the shared geometry map instead of `sortedRowsForPatch(...)[parent_cell]`. The helper validates complete unique parent coverage before mutation. The integration test now includes a regression where parent gas rows are deliberately swapped before refinement; child states still map to the correct physical parent cell and volume-integrated mass, momentum, and total energy are conserved.

Limit: prolongation is piecewise-constant. No slope-limited/high-order conservative prolongation is claimed.

### H3.3 -- Conservative restriction on derefine

Verdict: **accepted for local conservative production derefinement under the explicit octant contract**.

`derefineProductionPatchInSimulationState` now builds row-by-patch-cell maps for every child patch using the shared patch-local indexing helper. Restriction uses geometry-based child cell locations rather than sorted child-row order. The integration test now shuffles rows within child patch blocks before derefine and verifies conservation.

The pre-existing contract checks for child count, child level, extents, octants, duplicates, and patch relationships are preserved.

Limit: the production derefine path remains local/synchronous. It does not claim distributed derefine or subcycled fine/coarse synchronization.

### H3.4 -- Same-level and coarse-fine ghost fill

Verdict: **accepted for local ghost fill metadata and diagnostics**.

`fillAmrHydroGhostCells` remains live in the production AMR hydro path. Ghost descriptors now carry auditable fill statuses beyond the original unfilled states, including filled physical boundary, filled same-level, filled coarse-to-fine, filled fine-to-coarse, stale remote rejection, skipped remote, and missing source states. Diagnostics now include skipped remote/imported ghosts, stale epoch rejections, missing source records, and unresolved ghosts.

Tests inspect per-ghost status, not only aggregate counters.

Limit: full MPI/distributed AMR ghost exchange is not implemented. Remote/imported ghost support is an epoch/identity correctness seam and local contract, not a production multi-rank exchange path.

### H3.5 -- Flux-register generation during hydro sweeps

Verdict: **accepted for local solver-emitted register records**.

`HydroCoreSolver` remains geometry-agnostic and emits `HydroFluxRegisterRecord` values through `HydroFluxRegisterSink` when patch geometry marks coarse-fine faces. The records carry stable coarse gas-cell IDs where reflux targets are required. Existing and strengthened tests verify records are produced by the hydro face loop rather than by a purely synthetic/manual production path.

Limit: no persistent flux-register restart schema was added.

### H3.6 -- Automatic reflux application

Verdict: **accepted for local synchronized sweeps with unsafe-register skipping**.

The compile blocker in `applyFluxRegistersToSimulationState` is fixed by using the shared patch-local indexing helper. Reflux resolves targets by stable `coarse_gas_cell_id`, then validates current patch ownership and patch-local cell position before mutation. Incomplete coarse-only or fine-only records, area mismatches, missing targets, and wrong patch ownership are skipped or rejected and counted. Tests cover complete application, incomplete registers, area mismatch, missing ID, wrong patch ownership, row reorder, and nonzero momentum/energy diagnostics.

Limit: reflux is applied at the end of the current local AMR hydro stage and is safe only for synchronized local sweeps. This is not a complete deferred/subcycling reflux solution.

### H3.7 -- AMR hydro validation tests

Verdict: **accepted as CI-scale guards, not scientific validation decks**.

The requested H3.7 tests compile and run in the CPU/no-HDF5/no-MPI configuration:

- `integration_amr_hydro_shock_tube`
- `integration_amr_hydro_sedov`
- `integration_amr_synchronization_stress`

CMake now also provides build-target aliases with those exact names, while CTest continues to run the same names as tests backed by `test_*` executables.

Limit: these are regression/smoke guards with finite-state, positivity, conservation, and qualitative-behavior checks. They are not publication-grade convergence studies.

### H3.8 -- AMR-aware MPI patch/gas-sidecar migration

Verdict: **partially accepted as local contract coverage only**.

`integration_amr_patch_migration` passes and covers local/pseudo-rank pack/commit behavior, patch descriptors, gas sidecars, identity coverage, and stale ghost epoch rejection.

Limit: no two-rank MPI AMR migration or distributed AMR ghost-exchange smoke test was run because MPI was disabled in this build. H3.8 must not be advertised as production distributed AMR hydro.

### H3.9 -- Closeout audit

Verdict: **accepted as an honest local repair closeout**.

This document records exact commands, configuration, test results, skipped evidence, and remaining limitations. Older aspirational claims that depended only on docs/test presence have been replaced with build/test-backed claims.

## Tests not run / skipped evidence

- **HDF5 AMR hydro restart equivalence:** not run because the requested environment was configured with `COSMOSIM_ENABLE_HDF5=OFF`. The passing `integration_restart_equivalence_hydro_toy` is a toy/H1-style hydro restart guard, not a production AMR hydro run/restart/run equivalence proof.
- **MPI AMR migration/ghost exchange:** not run because the requested environment was configured with `COSMOSIM_ENABLE_MPI=OFF`. No full distributed AMR hydro ghost exchange is claimed.
- **AMR subcycling/deferred flux-register persistence:** not implemented or tested.
- **Publication-grade hydro validation:** not attempted; the H3.7 tests are CI-scale guards.

## Remaining known limitations

- Production AMR hydro is local/single-rank in the tested configuration.
- Full MPI AMR ghost exchange remains future work.
- Full AMR hydro restart equivalence with HDF5 checkpoint/reload remains future work.
- Reflux drains inside the current local hydro stage; there is no persistent deferred-register or subcycling completion model.
- Conservative refine/prolongation is piecewise-constant only.
- No cooling, star formation, feedback, metals, MHD, radiation, chemistry, or galaxy-formation subgrid coupling is closed by H3.

## Allowed claims after this repair

- `cosmosim_amr` compiles in the CPU/no-HDF5/no-MPI test configuration.
- Registered AMR/hydro tests compile and pass in the targeted CTest run.
- The production AMR hydro path remains `SimulationState`-owned and AMR-driven when patch coverage exists.
- Patch-local mapping, refine, derefine, and reflux validation use explicit geometry and stable identity rather than sorted dense-row order.
- Ghost fill has per-ghost fill/skip/reject metadata and expanded diagnostics.
- Unsafe reflux registers are skipped/reported rather than silently applied.

## Forbidden claims after this repair

- Full distributed AMR hydro is complete.
- Production AMR hydro restart equivalence has been proven.
- AMR time subcycling is implemented.
- The CI shock/Sedov/stress guards are publication-grade validation.
- H3 closes galaxy-formation physics or any subgrid model readiness.
