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
