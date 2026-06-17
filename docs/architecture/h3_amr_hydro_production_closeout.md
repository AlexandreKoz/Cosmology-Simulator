# H3 AMR hydro production integration closeout

This note records the narrow production AMR-hydro integration added after the H3 audit.  The ownership
rule is intentionally strict: `core::SimulationState` is the only persistent owner of production gas state.
`amr::PatchHierarchy` and `AmrPatch::conservedView()` remain available for AMR topology/scaffold tests,
but the production hydro path now operates on `SimulationState::{cells, gas_cells, gas_cell_identity,
patches}` and scatters by stable `gas_cell_id`.

## Production path

`src/workflows/reference_workflow.cpp::HydroStageCallback::onStage` now checks for complete AMR patch
coverage with `amr::hasProductionAmrHydroCoverage(state)`.  When coverage is available it calls the
production AMR orchestrator instead of building one global fallback Cartesian patch.  The fallback path is
preserved for toy/uniform states without full patch coverage.

The orchestrator is implemented in `include/cosmosim/amr/amr_hydro_orchestrator.hpp` and
`src/amr/amr_hydro_orchestrator.cpp`.  It derives patch descriptors from `PatchSoa` plus cell centers,
builds patch-local `AmrHydroPatchGeometry` views from `SimulationState`, fills AMR ghosts, runs
`HydroCoreSolver::advancePatchActiveSetWithScratch`, scatters patch-local results back by
`gas_cell_id`, and drains flux registers through a `SimulationState` reflux path before the hydro stage
returns.

## Flux registers and reflux

`populateAmrHydroFluxRegisterFaces` marks coarse-fine ghost faces in `HydroPatchGeometry` so the
existing hydro solver emits real `HydroFluxRegisterRecord` values from the Riemann face loop.  These
records are accumulated by `amr::FluxRegisterAccumulator`.  Reflux is then applied to authoritative
coarse rows in `SimulationState` by resolving the coarse patch/cell metadata through the current
`GasCellIdentityMap` row coverage.  Diagnostics now include corrected mass, momentum x/y/z, total
energy, and an internal-energy refresh delta.

## Production regrid helpers

`refineProductionPatchInSimulationState` and `derefineProductionPatchInSimulationState` provide
conservative, `SimulationState`-owned refinement/derefine helpers for synchronized production tests.
They allocate nonzero child/replacement `gas_cell_id` values, rebuild `GasCellIdentityMap`, update patch
and sidecar lanes in one commit, bump the cell-index generation, and conserve volume-integrated mass,
momentum x/y/z, and total energy in the tested cases.

## Acceptance test

`tests/integration/test_amr_production_hydro_integration.cpp` is the new production acceptance guard.  It
checks:

- AMR patch coverage and geometry derived from `SimulationState`.
- AMR ghost fill and real hydro sweeps over coarse/fine patches.
- Solver-emitted flux-register entries and reflux into `SimulationState` rows.
- Gas identity sidecar consistency after the sweep.
- Conservative production refine and derefine through `SimulationState`, not `AmrPatch::m_conserved`.

## Remaining limitations

This patch does not implement distributed remote AMR ghost exchange.  Remote ghost epochs and
read-only contracts remain local/test-facing scaffolding until the MPI patch migration stage.  Patch
geometry is derived from `PatchSoa` cell ranges and cell centers because `PatchSoa` still lacks explicit
origin/extent/cell-dimension lanes.  That derivation is checked, but an explicit restart-authoritative patch
geometry schema would be cleaner before large production AMR runs.  Reflux state is drained inside the
hydro stage; no new persistent flux-register restart schema was added.

## 2026-06-16 hardening update: explicit patch geometry, stable-ID reflux, and safe incomplete-register handling

This patch tightens the production AMR hydro path without promoting `amr::AmrPatch::m_conserved` to production truth. Persistent hydro state remains owned by `core::SimulationState` lanes: `cells`, `gas_cells`, `gas_cell_identity`, and `patches`.

Implemented hardening:

- `core::PatchSoa` now carries restart-authoritative patch descriptor lanes: parent patch ID, Morton key, origin, extent, and cell dimensions. Production AMR hydro coverage requires these explicit lanes instead of silently deriving patch geometry from sorted cell centers.
- AMR patch-local hydro mapping is row-order robust. `buildAmrHydroPatchGeometry` computes `(i,j,k)` from cell centers and explicit patch geometry, builds a patch-local permutation, rejects duplicate/missing cells, and keeps scatter keyed by stable `gas_cell_id`.
- Production refine/derefine helpers validate caller-provided ID ranges before mutating state and provide seedless overloads that scan for non-colliding patch and gas-cell ID ranges.
- Production derefine now validates the exact octant contract for eight children, including level, extents, origin octants, duplicate octants, and cell dimensions. Restricted parent cells preserve a parent particle ID only when all contributors agree; time bins are remapped conservatively from child minima rather than blindly reset.
- Flux-register records now carry a stable coarse gas-cell correction target. Reflux applies by resolving `coarse_gas_cell_id` through `GasCellIdentityMap`, not by sorted patch rows.
- Reflux skips incomplete or area-mismatched registers. It never applies a missing coarse or fine contribution as an implicit zero side. Diagnostics count complete, incomplete, area-mismatched, and missing-target registers.
- Restart schema is bumped to `cosmosim_restart_v16` so explicit patch geometry lanes are serialized and hashed. Older v14/v15-compatible inputs are still accepted with legacy zero-geometry backfill, but production AMR hydro will not enter unless explicit geometry is present.

Still intentionally deferred:

- True distributed/MPI AMR ghost exchange remains future work. Current production AMR hydro is hardened for local/single-rank patch coverage and keeps remote/imported ghost contracts explicit, but it does not claim multi-rank AMR ghost exchange completeness.
- AMR hydro restart equivalence is not a dedicated standalone test yet. The v16 patch-geometry restart lanes are covered by restart schema/round-trip tests, but a full run/restart/run production AMR hydro comparison should still be added before claiming H3 final closeout.

## 2026-06-16 production repair closeout: build, row-order mapping, ghost status, and reflux acceptance

A subsequent H3 repair pass fixed the compile failure in `applyFluxRegistersToSimulationState` and reran the CPU/no-HDF5/no-MPI AMR/hydro test set. The production AMR hydro path still treats `core::SimulationState` as authoritative state and still runs from the workflow whenever production AMR patch coverage exists.

Additional hardening in this pass:

- Added `include/cosmosim/amr/amr_patch_indexing.hpp` as a shared AMR helper for patch-local `(i,j,k)` and linear cell indexing from explicit `PatchDescriptor` geometry and gas-cell center coordinates.
- Reused that helper in hydro geometry construction, production refine, production derefine, and reflux target validation.
- Removed remaining production assumptions that sorted dense gas rows are equivalent to patch-local physical order.
- Added row-reorder regression coverage for production refine, derefine, and reflux target application.
- Extended AMR ghost descriptors with filled/skipped/rejected/missing-source statuses and expanded ghost-fill diagnostics.
- Preserved unsafe-register behavior: incomplete, area-mismatched, missing-target, and wrong-owner reflux records are skipped or rejected and counted instead of being applied.
- Added CMake build-target aliases for the H3.7 integration guards `integration_amr_hydro_shock_tube`, `integration_amr_hydro_sedov`, and `integration_amr_synchronization_stress`; the backing executables remain `test_*` targets and CTest uses the integration names.

Evidence gathered in this repair pass:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOSMOSIM_ENABLE_TESTS=ON -DCOSMOSIM_ENABLE_HDF5=OFF -DCOSMOSIM_ENABLE_MPI=OFF
cmake --build build --target cosmosim_amr -j2
cmake --build build --target test_integration_amr_production_hydro_integration -j2
./build/test_integration_amr_production_hydro_integration
cmake --build build --target integration_amr_hydro_shock_tube -j2
cmake --build build --target integration_amr_hydro_sedov -j2
cmake --build build --target integration_amr_synchronization_stress -j2
ctest --test-dir build --output-on-failure -R "amr|AMR|hydro"
```

Result: `cosmosim_amr` built successfully and the targeted AMR/hydro CTest run passed 24/24 tests.

Remaining limitations are unchanged in scope: no full MPI AMR ghost exchange, no HDF5 production AMR hydro restart-equivalence proof in this CPU/no-HDF5 build, no AMR time subcycling, and no persistent deferred flux-register restart schema.
