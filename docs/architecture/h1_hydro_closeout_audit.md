# H1 hydro closeout audit

Date: 2026-06-11

Mode: Audit, with repo-local documentation requested by prompt H1.8.

## Verdict

H1 is ready to close from code and registration inspection: the live hydro path no longer depends on a hidden one-dimensional x-chain geometry, and the H1 CPU validation surface is registered. H2 can begin after the requested CMake/CTest command bundle is rerun in an environment where `cmake` is available on `PATH`.

This is not a publication-grade hydro claim. The registered classical cases are CI-scale sentinels for positivity, conservation, symmetry, source coupling, and qualitative trends. Publication-grade Sedov/Evrard convergence and self-gravitating collapse remain later validation work.

## Evidence Basis

Prompt-specified surfaces inspected:

- H1 audit/plan: `docs/architecture/hydro_h1_patch_plan.md`.
- Live workflow: `src/workflows/reference_workflow.cpp`.
- Public hydro interfaces: `include/cosmosim/hydro/`.
- Hydro implementations: `src/hydro/`.
- Hydro unit tests: `tests/unit/test_hydro_*.cpp`.
- Hydro integration tests: `tests/integration/test_hydro_*.cpp`.
- Validation tests: `tests/validation/`.
- Build registration: `CMakeLists.txt`.
- Hydro docs: `docs/hydro_core_solver.md`, `docs/validation_ladder.md`.

Repository-orientation surfaces inspected under the agent contract:

- `README.md`
- `CONTRIBUTING.md`
- `docs/architecture/overview.md`
- `docs/architecture/developer_workflow_contract.md`
- `docs/architecture/runtime_truth_map.md`
- `docs/architecture/adr_runtime_truth_ownership.md`
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`

## 1D Chain Assumption Search

Search command:

```powershell
rg -n "TODO|FIXME|owner_cell \+ 1|neighbor_cell = owner|cell_index \+ 1|normal_x = 1\.0|normal_y = 1\.0|normal_z = 1\.0|dt_over_dx_code|dt_over_cell_width_code|HydroStageCallback::rebuildGeometryIfNeeded|rebuildGeometryIfNeeded" include/cosmosim/hydro src/hydro src/workflows/reference_workflow.cpp
```

Findings:

- No `TODO`/`FIXME` comments were found in production hydro paths by that search.
- `HydroStageCallback::rebuildGeometryIfNeeded(...)` now builds `HydroCartesianPatchSpec`, keys the cache on dimensions, origin, widths, timestep, hydro boundary policy, and patch signature, then calls `makeCartesianPatchGeometry(...)` and `appendCartesianBoundaryGhostFaces(...)` (`src/workflows/reference_workflow.cpp:4269`).
- The reconstruction policy is rebuilt with per-axis `dt_over_cell_width_code = {dt/dx, dt/dy, dt/dz}` while keeping scalar `dt_over_dx_code` only as compatibility fallback (`src/workflows/reference_workflow.cpp:4306`, `include/cosmosim/hydro/hydro_reconstruction.hpp:26`).
- `makeCartesianPatchGeometry(...)` builds explicit x, y, and z internal faces with axis metadata, transverse face areas, owner-minus and neighbor-plus stencil rows, and unit normals (`src/hydro/hydro_cartesian_patch.cpp:84`).
- Remaining `normal_x = 1.0`, `normal_y = 1.0`, and `normal_z = 1.0` hits are axis-specific Cartesian face construction, not a production x-only fallback (`src/hydro/hydro_cartesian_patch.cpp:111`, `src/hydro/hydro_cartesian_patch.cpp:130`, `src/hydro/hydro_cartesian_patch.cpp:149`).

## Production Hydro Path

The live workflow hydro callback uses the production Cartesian geometry path:

- `onStage(...)` calls `rebuildGeometryIfNeeded(context.state, context.integrator_state.dt_time_code)` before active face construction (`src/workflows/reference_workflow.cpp:3920`).
- Active faces are selected from explicit `m_geometry.faces` when either owner or neighbor is active (`src/workflows/reference_workflow.cpp:4317`).
- The solver receives `m_geometry`, `m_reconstruction`, `m_riemann_solver`, reusable scratch/cache, and the active cell/face view through `advancePatchActiveSetWithScratch(...)` (`src/workflows/reference_workflow.cpp:3976`).
- The hydro callback still requires the temporary particle-bound gas-cell contract before packing/unpacking gas rows (`src/workflows/reference_workflow.cpp:3918`). That is acceptable for H1 and is the handoff boundary for H2.

`HydroCoreSolver` validates explicit face ownership and unit normals, rejects non-authoritative real-cell owners, allows explicit ghost rows through metadata, and accumulates conservative owner/neighbor flux deltas using `face.area_comoving / geometry.cell_volume_comoving` (`src/hydro/hydro_core_solver.cpp:62`, `src/hydro/hydro_core_solver.cpp:624`).

`MusclHancockReconstruction` consumes explicit `HydroFace::axis`, per-axis CFL values, owner-minus/neighbor-plus stencil rows, and all three velocity components. It does not reconstruct by raw contiguous row inference (`src/hydro/hydro_reconstruction.cpp:66`, `src/hydro/hydro_reconstruction.cpp:181`).

## Boundary and Ghost Ownership

Boundary behavior is documented in `docs/hydro_core_solver.md` and implemented in `src/hydro/hydro_boundary_conditions.cpp`.

Confirmed code behavior:

- `appendCartesianBoundaryGhostFaces(...)` appends x, y, and z lower/upper ghost faces and rejects interior/imported-MPI as physical-boundary modes (`src/hydro/hydro_boundary_conditions.cpp:170`).
- Physical boundary ghosts are writable scratch rows; imported MPI ghosts are read-only and skipped by `fillHydroBoundaryGhostCells(...)` (`src/hydro/hydro_boundary_conditions.cpp:188`).
- Reflective boundaries reverse only the velocity component normal to the face axis (`src/hydro/hydro_boundary_conditions.cpp:141`).

Confirmed test behavior:

- `unit_hydro_boundary_conditions` checks periodic, open, and reflective ghost rows; verifies x/y/z ghost metadata; verifies owner real cells are not mutated; verifies reconstruction consumes ghost state; and verifies imported ghost metadata remains read-only (`tests/unit/test_hydro_boundary_conditions.cpp:79`).

## Conservation Accounting

`HydroProfileEvent::conservation` records before/after totals, flux deltas, source deltas, floor deltas, floor count, and residuals (`include/cosmosim/hydro/hydro_core_solver.hpp:264`). The solver computes these around active-cell updates (`src/hydro/hydro_core_solver.cpp:569`) and reports residual as:

```text
after - before - flux_delta - source_delta - floor_delta
```

Confirmed tests:

- `unit_hydro_core_solver` checks source-term separation and residuals (`tests/unit/test_hydro_core_solver.cpp:223`).
- `integration_hydro_conservation_periodic` runs a closed 4x3x2 periodic box with x/y/z periodic faces, non-unit transverse areas, no source terms, no floors, and conserved mass/momentum/energy totals over 20 steps (`tests/integration/test_hydro_conservation_periodic.cpp:89`).

## Validation Registration

CTest registration evidence in `CMakeLists.txt`:

- `unit_hydro_reconstruction`, `unit_hydro_boundary_conditions`, `unit_hydro_riemann`, `unit_hydro_core_solver` (`CMakeLists.txt:517`, `CMakeLists.txt:521`, `CMakeLists.txt:525`, `CMakeLists.txt:657`).
- `integration_hydro_sod_like`, `integration_hydro_conservation_periodic`, `integration_hydro_axis_symmetry` (`CMakeLists.txt:667`, `CMakeLists.txt:671`, `CMakeLists.txt:675`).
- `integration_restart_equivalence_hydro_toy` is registered through the restart-equivalence loop and labeled `restart;scheduler;hydro;gas;equivalence;integration` (`CMakeLists.txt:763`, `CMakeLists.txt:774`).
- `validation_hydro_classics` is registered and labeled `validation;integration;hydro` (`CMakeLists.txt:899`).

Classical hydro cases covered by `validation_hydro_classics`:

- Sedov blast (`tests/validation/test_validation_hydro_classics.cpp:53`).
- Noh converging inflow (`tests/validation/test_validation_hydro_classics.cpp:95`).
- Gresho vortex (`tests/validation/test_validation_hydro_classics.cpp:138`).
- Kelvin-Helmholtz (`tests/validation/test_validation_hydro_classics.cpp:213`).
- Evrard-style analytic-gravity collapse toy (`tests/validation/test_validation_hydro_classics.cpp:271`).

Exact CTest regex for H1 closeout registration/evidence:

```powershell
ctest --preset test-cpu-debug --output-on-failure -R "^(unit_hydro_reconstruction|unit_hydro_boundary_conditions|unit_hydro_riemann|unit_hydro_core_solver|integration_hydro_sod_like|integration_hydro_conservation_periodic|integration_hydro_axis_symmetry|validation_hydro_classics|integration_restart_equivalence_hydro_toy)$"
```

Prompt-required broader regex:

```powershell
ctest --preset test-cpu-debug --output-on-failure -R "hydro|validation_hydro|restart_equivalence_hydro_toy"
```

## Restart Equivalence

`integration_restart_equivalence_hydro_toy` exists and is registered. In HDF5-enabled builds, it uses a 4x3x2 Cartesian hydro toy state, evolves with `HydroCoreSolver::advancePatch(...)`, writes/restarts through the shared restart-equivalence harness, and compares direct/restarted gas thermodynamics and velocities (`tests/integration/test_restart_equivalence_hydro_toy.cpp:35`, `tests/integration/test_restart_equivalence_hydro_toy.cpp:79`, `tests/integration/test_restart_equivalence_hydro_toy.cpp:134`).

No H1 closeout evidence indicates a restart schema change in this pass. Hydro scratch, geometry, primitive cache, and profile counters remain transient, while gas thermodynamic/identity lanes are persisted through `SimulationState`.

## Commands and Outcomes

Commands attempted in this audit:

```powershell
cmake --preset cpu-only-debug
```

Outcome: blocked before configure. PowerShell reported:

```text
cmake : O termo 'cmake' não é reconhecido como nome de cmdlet, função, arquivo de script ou programa operável.
```

Additional local tooling checks:

```powershell
Get-Command cmake -ErrorAction SilentlyContinue | Format-List *
Test-Path -LiteralPath 'C:\Program Files\CMake\bin\cmake.exe'
Test-Path -LiteralPath 'C:\Program Files (x86)\CMake\bin\cmake.exe'
```

Outcomes:

- `Get-Command cmake` found no command.
- Both standard Windows CMake install paths returned `False`.
- `codex_app.load_workspace_dependencies` reported no bundled workspace runtime dependencies configured.

Therefore the requested commands below were not executable in this local shell and must be rerun after CMake is installed or placed on PATH:

```powershell
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "hydro|validation_hydro|restart_equivalence_hydro_toy"
```

Docs/interface hygiene command:

```powershell
& 'C:\Program Files\Git\bin\bash.exe' -lc 'export PATH="/usr/bin:/bin:/mingw64/bin:$PATH"; cd "/g/Cosmology Simulator"; ./scripts/ci/check_repo_hygiene.sh'
```

Outcome: pass.

HDF5 restart-schema retest was not required by this audit because no restart schema was changed here.

## Remaining H2 Blockers

H2 may begin, but it should be scoped to gas-cell identity promotion and must close these blockers before declaring H2 stable:

1. `H2-HYDRO-GAS-ID-AUTHORITY`: promote `gas_cell_id` / `parent_particle_id` from particle-bound temporary contract support into the production-authoritative gas-cell identity interface for hydro, with explicit remap behavior for reorder, migration, restart, and active views.
2. `H2-HYDRO-RESTART-SCHEMA`: if identity promotion changes persistent gas identity semantics, bump/version the restart schema, document compatibility behavior, and add negative/roundtrip tests before accepting old ambiguous payloads.
3. `H2-HYDRO-WORKFLOW-DECOUPLING`: remove naked dependence on the particle-bound gas-cell contract from live hydro pack/unpack only after replacement helper APIs and tests prove stable `gas_cell_id` lookup in the workflow path.
4. `H2-HYDRO-MPI-GHOST-ID`: define imported hydro ghost ownership and conservative correction exchange in terms of stable gas-cell identity, not local row position or parent-particle coincidence.

## Remaining H3 Blockers

H3 must not start until an H2 closeout says stable `gas_cell_id` is production-authoritative. Remaining H3 blockers are:

1. `H3-HYDRO-AMR-PATCH-OWNERSHIP`: define AMR patch ownership, prolongation/restriction, synchronization, and reflux/conservation implications against the stable gas-cell identity contract.
2. `H3-HYDRO-MULTIRANK-CONTINUATION`: provide multi-rank restart/migration evidence for gas-cell identity and hydro ghost/correction ownership, including stale ghost rejection and exact continuation tests.
3. `H3-HYDRO-VALIDATION-CAMPAIGN`: upgrade CI-scale Sedov/Evrard guards into documented validation campaigns with reference profiles or convergence targets before using them as publication-grade evidence.
4. `H3-HYDRO-PERFORMANCE-FLOOR`: add benchmark/profiling evidence for the production multidimensional hydro path once identity promotion and AMR-facing ownership are stable.

## Reproducibility Impact

This audit patch is documentation-only. It changes no solver numerics, configuration keys, normalized config output, snapshot/restart schema, HDF5 dataset names, scheduler authority, gas-field persistence, or output naming.
