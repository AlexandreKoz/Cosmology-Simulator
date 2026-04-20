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

Residual tree forces use the existing tree softening policy (`TreeSofteningPolicy`) first, then apply the Gaussian short-range factor.
The effective pair softening is the same as the tree solver:

- source epsilon resolution: per-particle sidecar -> species epsilon table -> scalar fallback,
- target epsilon resolution: per-active-target sidecar -> scalar fallback,
- pair rule: `epsilon_pair = max(epsilon_source, epsilon_target)`.

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

### Distributed short-range export/import contract (Phase 2 active-target model)

The runtime now records and enforces **response coverage** for each active-target batch slot. If a target was exported to `N` peers after periodic-bound pruning, the caller must receive exactly `N` responses before accumulation; missing or duplicate responses are treated as hard runtime errors rather than silently under-accumulating the short-range residual.

When `world_size > 1` and MPI is enabled, short-range residual evaluation is no longer rank-replicated:

1. **Target owner local pass**:
   - owner rank computes local-local short-range residual for its active targets against local sources.
2. **Export requests by peer**:
   - for each remote peer rank, owner exports active target batches capped by `tree_exchange_batch_bytes`.
   - request packet fields:
     - `batch_token` (`uint32`): active-batch start slot on owner rank.
     - `request_id` (`uint32`): index inside batch (`[0, batch_size)`).
    - `target_x_comoving`, `target_y_comoving`, `target_z_comoving` (`double`): target coordinates.
    - `target_softening_epsilon_comoving` (`double`): resolved target-side epsilon for the same pair law.
3. **Remote source evaluation**:
   - remote rank evaluates each request target against its local tree/source particles with the same:
    - opening criterion (selected MAC: geometric `l/r < theta` or COM-distance-aware `(l+δ_com)/r < theta`),
    - multipole truncation order (`kMonopole` or `kQuadrupole`),
    - softening policy (`softenedInvR3`),
    - split factor (`F_SR(r)`),
    - cutoff policy (`r_cut` AABB pruning + leaf pair skip).
4. **Return partial responses**:
   - response packet fields:
     - `batch_token` (`uint32`) and `request_id` (`uint32`) echoed from request,
     - `accel_x_comoving`, `accel_y_comoving`, `accel_z_comoving` (`double`) partial acceleration from that remote rank.
5. **Owner accumulation and validation**:
   - owner sums returned partials across peers into active slots.
   - duplicates or missing responses are detected per peer/per batch by `(batch_token, request_id)` coverage checks; mismatches fail fast with an exception.

Ordering is deterministic by rank order and batch-local request order.

## Validation notes

- `tests/unit/test_tree_pm_split_kernel.cpp` checks split composition and mesh-cell-to-length derivation.
- `tests/integration/test_tree_pm_coupling_periodic.cpp` checks periodic coupling against a **minimum-image periodic direct reference** (not a full infinite-periodic Ewald sum), plus cutoff-pruning and PM/tree/split consistency checks.
- `tests/integration/test_tree_pm_coupling_periodic.cpp` also includes an MPI two-rank distributed short-range export/import check (active-target export, cutoff-boundary peers, one-rank reference agreement).
- `bench/bench_tree_pm_force_error_map.cpp` maps force error against a **periodic spectral + direct short-range proxy reference** across PMGRID/ASMTH/RCUT sweeps and writes `validation/artifacts/tree_pm_force_error_map.csv`.

## Migration notes

- `TreePmOptions` now includes `tree_exchange_batch_bytes` (default `4 MiB`), and workflow wiring sets it from `numerics.treepm_tree_exchange_batch_bytes`.
- Existing callers constructing `TreePmOptions` without this field keep prior behavior via default value.


## Distributed source ownership note

In multi-rank TreePM, the short-range tree is intended to be built from **rank-owned source particles only**. Active targets may be a strict owned subset of that local source set, and export/import validation should compare gathered owner-local results against a single-rank reference rather than relying on replicated full-state digests.


### Boundary-condition coupling

`TreePmCoordinator` now branches by `PmSolveOptions::boundary_condition`:
- `kPeriodic`: periodic FFT PM long-range + periodic minimum-image short-range residual.
- `kIsolatedOpen`: open-boundary doubled-domain free-space PM convolution + non-periodic short-range residual (minimum-image disabled).

Unsupported combinations fail fast:
- isolated/open PM with `world_size > 1` throws.
