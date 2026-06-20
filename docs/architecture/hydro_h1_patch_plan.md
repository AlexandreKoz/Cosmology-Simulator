# H1 fixed/patch hydro: production acceptance plan and evidence

**Current revision:** 2026-06-20  
**Scope:** fixed Cartesian patch hydro in the real `ReferenceWorkflow` path. This
is not a claim of publication-grade hydrodynamic validation.

## Production ownership contract

- `core::SimulationState` is restart-authoritative.
- A dense gas-cell row is a local storage position, never a Cartesian identity.
- `gas_cell_id`, `GasCellIdentityMap`, explicit patch geometry, and physical cell
  centers establish identity and topology.
- `HydroCoreSolver` remains geometry-agnostic. `HydroStageCallback` owns workflow
  layout construction, ghost lifecycle, CFL gating, gather/scatter, and diagnostics.
- Complete production AMR coverage takes the AMR hydro path. The fixed Cartesian
  compatibility path is allowed only when a complete, uniform Cartesian cell-center
  layout is proved or validated explicit fixed-patch metadata agrees with it.

## H1.0–H1.8 acceptance ledger

| Item | Status | Production evidence | Remaining boundary |
|---|---|---|---|
| H1.0 baseline audit and patch plan | **Accepted** | This document and the closeout identify the live workflow path, ownership contract, and proof obligations. | The plan does not substitute for science validation. |
| H1.1 3D Cartesian patch geometry | **Accepted, fixed/patch scope** | `workflows/internal/cartesian_gas_cell_layout.*` derives axes from actual centers, proves a complete uniform Cartesian product, maps geometry rows ↔ dense rows, and rejects duplicates, holes, off-lattice and ambiguous centers. `integration_reference_workflow_hydro_row_order` exercises shuffled rows. | Arbitrary non-Cartesian mesh geometry is not an H1 fixed-patch capability. |
| H1.2 ghosts and boundary conditions | **Accepted, fixed/patch scope** | Existing axis-aware periodic/open/reflective coverage remains; reordered periodic workflow coverage uses stable mapping rather than dense topology. | Distributed ghost exchange is covered by existing local contracts, not a multi-rank proof here. |
| H1.3 directional reconstruction | **Accepted, fixed/patch scope** | The workflow uses physical patch widths and geometry-order gather; reconstruction stays face-axis-aware for density, all velocity components and pressure. | Classical tests are guards, not convergence measurements. |
| H1.4 CFL enforcement | **Accepted, fixed/patch scope** | Workflow CFL metadata resolves stable `gas_cell_id`, explicit patch identity and physical directional widths. Invalid geometry is rejected before a fabricated width is used. | No claim of global multi-rank timestep optimality. |
| H1.5 conservation accounting | **Accepted, fixed/patch scope** | Closed periodic workflow equivalence compares mass, three momentum components and total energy by stable identity after reorder. | Does not establish high-resolution convergence of conservation error. |
| H1.6 classical validation ladder | **Partially accepted** | Sod/Sedov/Noh/Gresho/Kelvin–Helmholtz/Evrard guards remain deterministic CI-scale finite/positivity/regression checks. | No documented multi-resolution analytic-reference study, cross-code comparison, or publication-grade error analysis. |
| H1.7 hydro restart round-trip | **Accepted, HDF5-enabled fixed/patch scope** | `integration_restart_equivalence_reference_workflow_hydro` runs real `ReferenceWorkflow`, `HydroStageCallback`, geometry build, ghosts, CFL/schedulers, conservation and a safe-boundary dense-row reorder through checkpoint/reload. v20 persists KDK force cache keyed by stable particle/gas IDs. | Requires HDF5; no multi-rank migration/restart proof. |
| H1.8 closeout audit | **Accepted** | `h1_hydro_closeout_audit.md` records build mode, exact tests and explicit limitations. | Must be revised when scope or evidence changes. |

## Geometry policy

There is no production near-cubic `cell_count` fallback. The workflow either:

1. dispatches to valid AMR production coverage;
2. validates explicit fixed Cartesian patch metadata against live centers and identity;
3. derives and validates a full Cartesian layout from live centers; or
4. fails with an invariant-specific diagnostic.

A dense-row reorder can reuse physical topology only after the row↔geometry mapping is
rebuilt/validated. Cache keys include identity generation, geometry signature and relevant
physical geometry/boundary policy inputs; stale dense mappings are not reusable.

## Validation ladder (intentionally ordered)

1. CI sanity/guard tests: small deterministic finite/positivity and symmetry checks.
2. Deterministic regression and HDF5 restart equivalence.
3. Multi-resolution convergence against analytic/reference solutions.
4. Cross-code comparison on selected benchmarks.
5. Science-readiness benchmarks and documented accuracy budgets.

Only stages 1–2 are established by H1. This document makes no publication-ready claim.
