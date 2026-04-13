# Tree Gravity Solver (Barnes-Hut)

## Scope

This module implements short-range gravity for CosmoSim with an indexed Barnes-Hut octree:

- SoA node storage (`TreeNodeSoa`) for locality and vectorization-friendly traversal.
- Morton-based particle ordering for stable locality-first build inputs.
- Explicit stages: build, multipole accumulation, traversal.
- Monopole force baseline with Plummer-equivalent softening.
- Active-set traversal API for compact gravity kernel loops.
- Profile counters for build/traversal runtime and interaction statistics.

## Interfaces

- Header: `include/cosmosim/gravity/tree_gravity.hpp`
- Ordering helper: `include/cosmosim/gravity/tree_ordering.hpp`
- Softening helper: `include/cosmosim/gravity/tree_softening.hpp`

## Opening and force model

For accepted internal nodes, the solver applies:

\[
\mathbf{a}(\mathbf{x}) = G M \frac{\mathbf{x}_{\mathrm{com}} - \mathbf{x}}{\left(|\mathbf{x}_{\mathrm{com}} - \mathbf{x}|^2 + \epsilon^2\right)^{3/2}}
\]

The geometric opening criterion is:

\[
\frac{l}{r} < \theta
\]

where `l = 2 * half_size` and `r` is target-to-node-COM distance.

## Configuration and provenance implications

Assumptions and implications documented for reproducibility:

1. Current implementation is non-periodic short-range and does not apply PM long-range correction.
2. Multipoles are monopole-only in this revision; node arrays and stage separation leave room for quadrupole sidecars.
3. Softening uses `TreeSofteningPolicy::epsilon_comoving` in the same unit frame as particle positions.
4. Ordering is deterministic for identical inputs because Morton keys are stable-sorted.

Integration with normalized config/provenance can map directly to `TreeGravityOptions` fields:

- `opening_theta`
- `max_leaf_size`
- `gravitational_constant_code`
- `softening.epsilon_comoving`

## Tests and benchmark hooks

- Unit tests: `tests/unit/test_tree_gravity.cpp`
- Integration test: `tests/integration/test_tree_gravity_vs_direct.cpp`
- Benchmark: `bench/bench_tree_gravity.cpp`


## Integration validation tightening

- `tests/integration/test_tree_gravity_vs_direct.cpp` remains the isolated/non-periodic trust anchor.
- The test now checks both **tight** (`theta=0.4`) and **loose** (`theta=0.8`) opening settings against direct summation and reports max/mean relative errors.
- Stop/go line:
  - Tight opening must satisfy `max_rel_error < 0.03` and `mean_rel_error < 0.015`.
  - Loose opening is required to be no more accurate than tight opening, guarding against accidental criterion inversion.
- Diagnostics include visited/accepted node counters for auditability when tolerances fail.
