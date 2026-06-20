# H3 AMR Hydro Closeout Audit

Date: 2026-06-17  
Prompt range: H3.0-H3.9  
Mode: hardening pass after build/row-order/reflux repair

## Executive verdict

The H3 branch now has a dedicated AMR hydro restart-equivalence test in addition to the earlier local production AMR hydro, row-order, ghost-fill, reflux, and CI-scale validation guards. The new evidence supports a narrow claim: **single-rank synchronized production AMR hydro can be checkpointed through the HDF5 restart path and continued deterministically for the exercised AMR hydro scenario**.

This still does **not** prove Berger-Colella time subcycling, distributed production AMR hydro, persistent/deferred flux-register restart state, restart after MPI AMR migration, or publication-grade scientific validation. Those remain future work.

## Environment and configuration used

Tested in the local container with:

- Compiler: GNU C++ 14.2.0.
- CPU debug preset: HDF5 disabled, MPI disabled.
- HDF5 debug preset: HDF5 enabled, MPI disabled.
- MPI: not enabled and no real two-rank AMR migration/hydro smoke test was added in this pass.

Configuration commands used:

```bash
cmake --preset cpu-debug
cmake --preset hdf5-debug
```

## Commands run and outcomes

CPU/no-HDF5/no-MPI build and AMR guard tests:

```bash
cmake --build --preset build-cpu-debug --target cosmosim_amr -j2
cmake --build --preset build-cpu-debug --target \
  test_unit_amr_hydro_geometry \
  test_integration_amr_hydro_shock_tube \
  test_integration_amr_hydro_sedov \
  test_integration_amr_synchronization_stress \
  test_integration_amr_patch_migration \
  test_integration_amr_production_hydro_integration \
  test_integration_amr_reflux_conservation \
  test_integration_restart_equivalence_amr_hydro -j2
ctest --preset test-cpu-debug \
  -R "amr_hydro|amr_synchronization|amr_patch_migration|restart_equivalence_amr_hydro" \
  --output-on-failure
```

Outcome: the targeted CPU CTest regex passed after explicitly building the matched executables. In this CMake layout, `ctest` does not build missing test executables automatically.

HDF5 restart-equivalence build and tests:

```bash
cmake --build --preset build-hdf5-debug --target \
  test_integration_restart_equivalence_harness \
  test_integration_restart_equivalence_dm_only \
  test_integration_restart_equivalence_treepm \
  test_integration_restart_equivalence_hydro_toy \
  test_integration_restart_equivalence_multirate_bins \
  test_integration_restart_equivalence_output_enabled \
  test_integration_restart_equivalence_stochastic_sources \
  test_integration_restart_equivalence_amr_hydro -j2
ctest --preset test-hdf5-debug -R "restart_equivalence" --output-on-failure
```

Outcome: **8/8 HDF5 restart-equivalence tests passed**, including `integration_restart_equivalence_amr_hydro`.

## H3.0-H3.9 verdicts

### H3.0 -- AMR/hydro coupling ownership

Verdict: **accepted for the local production ownership contract**.

`core::SimulationState` remains the restart-authoritative owner of gas and hydro state. The production workflow still enters the AMR-driven hydro path when `amr::hasProductionAmrHydroCoverage(context.state)` is true. `amr::PatchHierarchy` and `AmrPatch::m_conserved` remain scaffold/topology/test infrastructure, not production hydro truth.

Limit: this does not imply distributed AMR hydro or subcycling.

### H3.1 -- AMR patches as first-class hydro geometry providers

Verdict: **accepted for explicit local patch geometry and row-order-independent mapping**.

AMR patch descriptors carry explicit restart-authoritative geometry lanes, and patch-local row mapping is built from geometry and stable gas identity rather than sorted dense row order. The restart-equivalence harness now compares explicit PatchSoa geometry lanes and gas identity records/generation.

Limit: production AMR coverage still assumes each active patch has complete local cell coverage in the current single-rank state.

### H3.2 -- Conservative prolongation on refine

Verdict: **accepted for local piecewise-constant conservative refinement**.

Refinement uses row-by-patch-cell maps built from explicit geometry. Existing row-shuffle regressions verify that parent lookup does not depend on sorted dense rows, and conservation checks remain tight.

Limit: no higher-order slope-limited conservative prolongation is claimed.

### H3.3 -- Conservative restriction on derefine

Verdict: **accepted for local conservative derefinement under the explicit child-octant contract**.

Restriction uses geometry-based maps for child patches. Row-order shuffle coverage and negative child-contract tests remain the relevant evidence.

Limit: no distributed derefine or subcycled fine/coarse synchronization is claimed.

### H3.4 -- Same-level and coarse-fine ghost fill

Verdict: **accepted for local ghost fill metadata and diagnostics**.

Ghost descriptors carry auditable fill/skip/reject/missing-source status, and diagnostics count more than filled totals. Local tests inspect individual ghost statuses.

Limit: remote distributed AMR ghost exchange is not production implemented.

### H3.5 -- Flux-register generation during hydro sweeps

Verdict: **accepted for local solver-emitted register records**.

The hydro solver remains geometry-agnostic and emits records through the flux-register sink. Reflux target records carry stable coarse gas-cell IDs where mutation needs identity validation.

Limit: records are not persisted as pending restart state.

### H3.6 -- Automatic reflux application

Verdict: **accepted for synchronized local sweeps with unsafe-register skipping**.

Reflux resolves by stable `coarse_gas_cell_id`, validates patch ownership and patch-local target mapping, applies complete/area-consistent records once, and skips incomplete, area-mismatched, missing-target, or wrong-owner records.

Limit: this is immediate synchronized-sweep reflux. It is not a deferred/subcycling reflux model.

### H3.7 -- AMR hydro validation tests

Verdict: **accepted as CI-scale guards only**.

`integration_amr_hydro_shock_tube`, `integration_amr_hydro_sedov`, and `integration_amr_synchronization_stress` build and run under the CPU preset. The shared validation helper asserts that production AMR coverage exists, then checks finite/positive state and gas identity coverage.

Limit: these tests are not convergence studies and are not cross-code scientific validation.

### H3.8 -- AMR-aware MPI patch/gas-sidecar migration

Verdict: **not production accepted; local contract coverage only**.

`integration_amr_patch_migration` was strengthened to check patch descriptor geometry lanes, gas sidecars, identity rebuild, row-order change across local commit, stale ghost epoch rejection, and production AMR coverage after local commit.

Limit: no real two-rank MPI AMR migration/hydro test was added or run. H3.8 currently has local payload/commit contract coverage only; it does not prove distributed AMR hydro patch migration, remote ghost exchange, AMR hydro after migration, or restart after migration.

### H3.9 -- Closeout audit

Verdict: **accepted as an honest local/HDF5 hardening closeout**.

This document records commands, feature gates, tests passed, and remaining limitations. It deliberately does not claim unavailable MPI, subcycling, deferred-register, or cross-code validation evidence.

## Restart-equivalence evidence added in this pass

New test: `tests/integration/test_restart_equivalence_amr_hydro.cpp`  
CMake/CTest name: `integration_restart_equivalence_amr_hydro`

The test builds an AMR hydro `SimulationState` with:

- explicit PatchSoa geometry lanes for coarse and fine patches;
- stable `gas_cell_id` coverage through `GasCellIdentityMap`;
- nontrivial density, pressure, and velocity fields;
- a coarse/fine interface that exercises ghost-fill and flux-register safety paths;
- direct vs checkpoint/reload/continue comparison through the existing restart-equivalence harness.

The harness now compares:

- gas primitive/conserved lanes through the existing state comparisons;
- cell centers, mass, time bins, and patch indices;
- gas-cell identity sidecar lanes;
- `GasCellIdentityMap` records and generation;
- explicit PatchSoa geometry lanes;
- scheduler, integrator, output cadence, and stochastic persistent state.

## Remaining known limitations

- No Berger-Colella AMR time subcycling.
- No level-by-level AMR dt hierarchy integrated with the production scheduler.
- No persistent/restart-safe pending flux-register model.
- No real two-rank MPI AMR patch/gas migration test.
- No distributed AMR ghost exchange.
- No restart-after-MPI-AMR-migration test.
- H3.7 tests are fast CI guards, not convergence or cross-code validation.

---

## Post-H3 credibility pass: AMR subcycling and pending flux registers

Date: 2026-06-17  
Mode: local AMR hydro credibility extension after H3 closeout

### What changed after the original H3 audit

The repository now contains a first real **local** AMR hydro subcycling path in the AMR hydro orchestrator. This is not a global scheduler-integrated Berger-Colella implementation. It advances explicit AMR levels recursively inside the local orchestrator, uses a validated refinement ratio, advances fine levels with smaller timesteps (`dt_fine = dt_coarse / refinement_ratio`), and delays complete reflux application until the required fine-substep coverage has accumulated.

The repository also now contains a restart-authoritative pending flux-register store in `core::SimulationState`. Pending records are stable-ID based and include coarse patch identity, coarse gas-cell identity, patch-local target index, level pair metadata, face orientation, area coverage, timestep/substep coverage, flux-integrated conserved quantities, and generation metadata. Complete records are resolved by stable `coarse_gas_cell_id` before mutation; incomplete, area-mismatched, stale-generation, wrong-owner, and missing-target records are skipped/rejected with diagnostics rather than silently applied.

### New test evidence

New or strengthened tests:

- `unit_amr_pending_flux_registers`
- `integration_amr_hydro_subcycling`
- `integration_restart_equivalence_amr_flux_registers`

The HDF5 restart-equivalence coverage now includes a case where an incomplete pending register is written to restart state, read back, completed by a later fine contribution, and then applied. The restart-equivalence harness also compares pending flux-register metadata and accumulated flux integrals.

### Updated verdict

The old limitations are partly closed:

- **AMR hydro subcycling:** implemented as local orchestrator-level subcycling for explicit AMR patch coverage. It is not yet scheduler-owned, MPI-distributed, or cosmological timeline-integrated.
- **Persistent pending flux-register state:** implemented in restart-authoritative `SimulationState` and serialized under `/state/amr_pending_flux_registers` in restart schema v17 (retained by current v19 checkpoints).
- **H3.7 validation:** still CI-scale. The new subcycling/pending-register tests are deterministic regression/credibility guards, not scientific convergence or cross-code validation.

### Remaining limitations

- No production scheduler ownership of subcycled AMR level timelines.
- No temporal interpolation model for coarse/fine ghost fill across arbitrary scheduler phases.
- No MPI-distributed AMR hydro subcycling.
- No restart after real MPI AMR patch migration.
- No publication-grade AMR shock/Sedov convergence deck.
- No cross-code comparison against AREPO, RAMSES, ENZO, Athena++, or another reference solver.

---

## Post-H3 temporal-boundary and validation-campaign update (2026-06-18)

The local AMR subcycling scope now includes restart-authoritative coarse temporal boundary histories and
conserved-state coarse-to-fine ghost interpolation at fine-update start time. HDF5 schema v19 retains and serializes
active histories and `integration_restart_equivalence_amr_temporal_ghosts` resumes a midpoint fill after
restart while an incomplete pending flux register is also present. This is accepted only for local two-level
orchestrator-level subcycling; it is not scheduler-owned general AMR subcycling or distributed MPI AMR.

The repository now has convergence/cross-code campaign machinery, but no long CHUÍ Sod/Sedov campaign
or external-code result was run in this pass. H3.7 tests remain CI guards. See
`docs/validation/amr_hydro_convergence_campaign.md` and
`docs/validation/amr_hydro_cross_code_protocol.md`.
