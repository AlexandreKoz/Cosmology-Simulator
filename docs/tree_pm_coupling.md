# TreePM coupling

## Scope

`TreePmCoordinator` implements a Gaussian TreePM split with explicit runtime-derived split and cutoff controls from mesh-cell UX parameters:

- `r_s = asmth_cells * Δmesh`
- `r_cut = rcut_cells * Δmesh`
- `Δmesh = box_size / PMGRID`

where `asmth_cells` and `rcut_cells` are the authoritative user-facing controls.

## Operational split contract

Long-range PM filter in Fourier space:

- `F_LR(k) = exp(-k^2 r_s^2)`

Short-range real-space residual factor:

- `F_SR(r) = erfc(r / (2 r_s)) + (2 / sqrt(pi)) * (r / (2 r_s)) * exp(-(r / (2 r_s))^2)`

Composition requirement (before explicit cutoff):

- `F_SR(r) + F_LR(r) = 1`

### Softening interaction contract

Residual tree forces use the existing tree softening policy (`TreeSofteningPolicy`) first, then apply the Gaussian short-range factor:

- `a_SR = a_tree_softened * F_SR(r)`

This is the explicit branch contract used by tests and diagnostics. It keeps the TreePM split as a decomposition of the same softened force law prior to truncation.

## Cutoff semantics (`r_cut`)

Tree residual traversal is now explicitly cutoff-bounded:

1. **Node-level pruning**: if the minimum distance from target to node AABB exceeds `r_cut`, the node is skipped.
2. **Node acceptance guard**: accepted internal nodes must be fully enclosed by `r_cut` (maximum target-to-node AABB distance <= `r_cut`), otherwise traversal descends.
3. **Leaf pair culling**: particle-particle residual contributions with `r > r_cut` are skipped.

This ensures `rcut_cells` changes actual traversal behavior (not diagnostics-only metadata).

## Runtime diagnostics and audibility

`TreePmDiagnostics` reports:

- mesh-cell controls (`asmth_cells`, `rcut_cells`) and derived lengths (`Δmesh`, `r_s`, `r_cut`),
- split composition checks (`composition_error_at_split`, `max_relative_composition_error`),
- factors at `r_s` and `r_cut`,
- residual cutoff traversal counters (`residual_pruned_nodes`, `residual_pair_skips_cutoff`, `residual_pair_evaluations`).

## Ownership and accumulation order

1. PM computes long-range acceleration on the compact active set.
2. Tree computes short-range residual with Gaussian factor and explicit `r_cut` truncation.
3. The active-set accumulator stores PM + residual totals.

## Validation notes

- `tests/unit/test_tree_pm_split_kernel.cpp` checks split composition and mesh-cell-to-length derivation.
- `tests/integration/test_tree_pm_coupling_periodic.cpp` checks periodic coupling against a **minimum-image periodic direct reference** (not a full infinite-periodic Ewald sum), plus cutoff-pruning and PM/tree/split consistency checks.
