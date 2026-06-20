# H1 fixed/patch hydro closeout audit

**Revision date:** 2026-06-20  
**Decision:** H1 is accepted only within the stated fixed/Cartesian-patch production
scope, with HDF5 enabled for the workflow restart proof. It is **not** publication-ready
hydrodynamic validation and it does not certify MPI restart/migration behavior.

## Evidence configuration

- CPU configuration: `cpu-only-debug` / `build-cpu-debug` / `test-cpu-debug`.
- HDF5 configuration: `hdf5-debug` / `build-hdf5-debug` / `test-hdf5-debug`.
- HDF5 was enabled and reported by CMake (`HDF5 1.14.5` in the acceptance environment).
- The repository has no `cosmosim_all` target; commands used actual discovered targets.

## Item-by-item closeout

### H1.0 — baseline audit and patch plan: accepted

The plan identifies `ReferenceWorkflow` → `HydroStageCallback` as the production
orchestration path and states the stable-identity, explicit-geometry ownership rules.

### H1.1 — 3D Cartesian hydro patch geometry: accepted in fixed/patch scope

`CartesianGasCellRowLayout` is a single reusable builder. It derives sorted axes from
live centers, checks complete Cartesian occupancy and uniform spacing, builds
`HydroCartesianPatchSpec`, and emits both dense→geometry and geometry→dense maps.
It rejects nonuniform, duplicate, missing, off-lattice and ambiguous layouts. Production
code no longer invents a near-cubic patch from `cell_count`.

Evidence: `integration_reference_workflow_hydro_row_order` and
`integration_restart_equivalence_reference_workflow_hydro`.

### H1.2 — ghost zones and boundary conditions: accepted in fixed/patch scope

The existing periodic, open/outflow and reflective boundary machinery remains intact.
The row-order integration test specifically verifies periodic workflow evolution through
physical mapping. Boundary descriptors remain transient/read-only; real updates scatter
back through stable identity and validated geometry mappings.

### H1.3 — dimensionally correct reconstruction: accepted in fixed/patch scope

The active workflow passes axis-aware geometry and directional widths to reconstruction.
The reorder regression compares the same 3D physical patch after geometry-order gather
and dense-order scatter. No scalar `dt/dx` compatibility lane is used as active workflow
policy.

### H1.4 — CFL enforcement: accepted in fixed/patch scope

CFL metadata derives patch ID, patch-local row, stable gas-cell ID and directional widths
from validated physical geometry/identity. Missing geometry fails loudly. Reordered rows
produce the same limiting physical cell diagnostics when compared by `gas_cell_id`.

### H1.5 — conservation accounting: accepted in fixed/patch scope

The periodic workflow regression compares closed-box mass, momentum x/y/z and total
energy between canonical and shuffled dense storage. Existing solver conservation
accounting remains the authority for flux/source/floor diagnostic breakdowns.

### H1.6 — classical validation ladder: partially accepted

Sedov, Noh, Gresho, Kelvin–Helmholtz and Evrard-related checks are retained as
low-resolution deterministic CI guards. They are not convergence studies, reference-error
studies, cross-code validation, or a publication-grade claim.

### H1.7 — hydro restart round-trip: accepted with HDF5, fixed/patch scope

The HDF5 integration test executes direct two-step workflow evolution against one step,
checkpoint/read, safe-boundary dense-row reorder, and one resumed workflow step. It
exercises workflow geometry reconstruction, ghosts, hydro gather/solve/scatter, directional
CFL scheduling, gas scheduler state, conservation reporting and KDK force-cache import.

Restart schema v20 stores the KDK cache with stable `particle_id` and `gas_cell_id` lanes;
import remaps cache values to current dense rows. This prevents a valid row reorder or
compaction from consuming accelerations by stale dense index.

### H1.8 — closeout audit: accepted

This closeout separates accepted scope from unproven scope and records HDF5 dependence.

## Remaining limitations

- No multi-resolution analytic/reference-solution study was executed in this repair.
- No cross-code benchmark or science-readiness campaign was executed.
- No MPI toolchain was available in the acceptance environment; multi-rank ghost exchange,
  migration and restart-after-migration remain unproven here.
- The accepted geometry scope is complete uniformly spaced Cartesian fixed patches plus
  existing explicit AMR coverage; it is not a general unstructured-mesh claim.
- The KDK cache is remapped by stable IDs on resume. It is not a substitute for distributed
  force-field serialization or a proof of exact continuation across rank repartitioning.
