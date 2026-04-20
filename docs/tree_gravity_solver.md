# Tree Gravity Solver (Barnes-Hut)

## Scope

This module implements short-range gravity for CosmoSim with an indexed Barnes-Hut octree:

- SoA node storage (`TreeNodeSoa`) for locality and vectorization-friendly traversal.
- Morton-based particle ordering for stable locality-first build inputs.
- Explicit stages: build, multipole accumulation, traversal.
- Monopole+quadrupole force truncation (order `l=2`) with Plummer-equivalent softening.
- Active-set traversal API for compact gravity kernel loops.
- Profile counters for build/traversal runtime and interaction statistics.

## Interfaces

- Header: `include/cosmosim/gravity/tree_gravity.hpp`
- Ordering helper: `include/cosmosim/gravity/tree_ordering.hpp`
- Softening helper: `include/cosmosim/gravity/tree_softening.hpp`

## Opening and force model

For accepted internal nodes, the solver applies a multipole expansion around node COM,
truncated at quadrupole order:

\[
\mathbf{a} = \mathbf{a}_{M} + \mathbf{a}_{Q}
\]

Monopole term:
\[
\mathbf{a}_{M} = G M \frac{\mathbf{r}}{(r^2+\epsilon^2)^{3/2}},\quad \mathbf{r}=\mathbf{x}_{\mathrm{com}}-\mathbf{x}
\]

Traceless quadrupole tensor:
\[
Q_{ij} = \sum_a m_a (3\,\Delta x_{a,i}\Delta x_{a,j} - |\Delta \mathbf{x}_a|^2\delta_{ij})
\]

Quadrupole acceleration correction:
\[
a_{Q,i}=\frac{G}{2}\left[2\frac{(Q\mathbf{r})_i}{(r^2+\epsilon^2)^{5/2}}-5\frac{(\mathbf{r}^TQ\mathbf{r})r_i}{(r^2+\epsilon^2)^{7/2}}\right]
\]

Supported opening criteria (MAC):

1. Geometric BH:
\[
\frac{l}{r} < \theta,\quad l=2h
\]
2. COM-distance aware BH (default):
\[
\frac{l+\delta_{\mathrm{com}}}{r} < \theta,\quad \delta_{\mathrm{com}} = |\mathbf{x}_{\mathrm{center}}-\mathbf{x}_{\mathrm{com}}|
\]

Fallback/accept rules:
- Leaves are always accepted and evaluated P2P.
- Internal nodes that fail MAC descend to children.
- Internal nodes accepted by MAC use multipole evaluation (monopole only or monopole+quadrupole per option).

The legacy geometric opening criterion is:

\[
\frac{l}{r} < \theta
\]

where `l = 2 * half_size` and `r` is target-to-node-COM distance.

## Configuration and provenance implications

Assumptions and implications documented for reproducibility:

1. Current implementation is non-periodic short-range and does not apply PM long-range correction.
2. Multipole order is explicit via `TreeMultipoleOrder` (`kMonopole` / `kQuadrupole`), defaulting to quadrupole.
3. Softening uses `TreeSofteningPolicy::epsilon_comoving` in the same unit frame as particle positions.
4. Ordering is deterministic for identical inputs because Morton keys are stable-sorted.

Integration with normalized config/provenance can map directly to `TreeGravityOptions` fields:

- `opening_theta`
- `opening_criterion`
- `multipole_order`
- `max_leaf_size`
- `gravitational_constant_code`
- `softening.epsilon_comoving`

## Tests and benchmark hooks

- Unit tests: `tests/unit/test_tree_gravity.cpp`
- Integration test: `tests/integration/test_tree_gravity_vs_direct.cpp`
- Benchmark: `bench/bench_tree_gravity.cpp`

## Distributed TreePM reuse contract

`TreePmCoordinator` reuses the same tree traversal kernel for remote short-range requests in distributed runs. Remote request targets are evaluated against rank-local tree/source data with unchanged:

- MAC (`l/r < theta`),
- softening (`softenedInvR3`),
- short-range split factor (`F_SR`),
- cutoff semantics (`r_cut` pruning + pair culling).


## Integration validation tightening

- `tests/integration/test_tree_gravity_vs_direct.cpp` remains the isolated/non-periodic trust anchor.
- The test now checks two particle distributions (quasi-uniform and clustered) and a 3-point MAC schedule (`theta = 0.35, 0.6, 0.9`) against direct summation.
- Stop/go line:
  - Tight opening must satisfy `max_rel_error < 0.02` and `mean_rel_error < 0.01` on both distributions.
  - Error and cost must be monotonic with opening (`theta↑ => error↑ and visited/P2P work↓`).
- Diagnostics include visited/accepted nodes and P2P interactions for auditable error/cost tradeoff checks.
