# M48-03 — Unified Tree-Construction Workspace and Warm-Rebuild Memory Closure

Date: 2026-09-28. Mode: focused production memory-architecture repair.

## Previous ownership

The periodic TreePM path retained three final unwrapped coordinate lanes plus two
population-sized `double` axis workspaces. `TreeGravitySolver` separately retained
the final Morton permutation, Morton keys, resolved source softening, and a
partition vector, while `buildMortonOrdering(...)` allocated a complete returned
ordering plus independent radix key/index scratch. On a warm rebuild, the old
retained ordering could therefore coexist with the new returned ordering and its
radix temporaries before move assignment committed the replacement.

Under the campaign's source-width model this was approximately 64 B/source of
retained periodic/tree particle-scale construction ownership and approximately
88 B/source during the warm Morton rebuild overlap.

## Unified physical workspace

`TreeGravitySolver` now owns and reuses exactly three population-scale
construction lanes:

- `m_ordering.morton_key`: `uint64_t x N` (8 B/source), the primary construction
  key lane and final retained Morton-key lane;
- `m_morton_key_scratch`: `uint64_t x N` (8 B/source), the alternate stable-radix
  key lane;
- `m_construction_index_scratch`: `TreeLocalIndex x N` (4 B/source under the
  checked 32-bit local-index contract), used first as radix-index scratch and
  then as recursive topology partition scratch.

The retained final permutation remains `m_ordering.sorted_particle_index`
(`TreeLocalIndex x N`, 4 B/source), and resolved source epsilon remains
`double x N` (8 B/source). The three final periodic unwrapped coordinates remain
in `TreePmCoordinator` (`double x N` each, 24 B/source). Thus the shared
construction workspace is 20 B/source and the complete particle-scale tree
construction ownership targeted by this campaign is 56 B/source:

```text
three final unwrapped coordinates     24 B/source
final retained permutation             4 B/source
resolved source epsilon                8 B/source
shared construction workspace         20 B/source
-----------------------------------------------
total                                  56 B/source
```

These are source/model-derived ownership figures, not measured process RSS.

## Periodic preprocessing reuse

`TreePmCoordinator` no longer owns dedicated wrapped/ordered axis vectors.
Beginning periodic preprocessing borrows the tree solver's two 64-bit key lanes
and simultaneously establishes the rebuild invalidation boundary.

For each axis, in sequence X -> Y -> Z:

1. every finite coordinate is wrapped into `[0, L_axis)` using the existing
   arithmetic;
2. the wrapped value is written directly into that axis's final output lane;
3. its standards-compliant `std::bit_cast<uint64_t>` representation is written
   into the primary construction key lane;
4. the same stable eight-pass byte-radix ordering alternates between the two
   64-bit construction lanes;
5. the largest circular gap, cyclic final gap, equal-gap minimum-anchor tie rule,
   and post-gap anchor selection are applied to the sorted representation;
6. the final output lane is transformed in place with the existing
   `wrapped < anchor ? wrapped + L : wrapped` rule.

No extra coordinate-sized lane survives or is allocated for the next axis.

## In-place Morton rebuild

The return-by-value `buildMortonOrdering(...)` path is replaced by
`buildMortonOrderingInPlace(...)`. It requires pre-sized retained primary and
shared scratch lanes and performs no population-scale allocation itself.

The tree solver initializes identity indices and the unchanged 21-bit Morton
keys directly into retained authoritative storage, then executes the existing
stable eight 8-bit radix passes between the retained primary arrays and shared
scratch. Because the pass count is even, the authoritative final key and
permutation end in the retained primary arrays without a population-scale final
copy. Static guards document that ownership invariant.

The old ordering is invalidated before any construction lane is borrowed, so a
warm build overwrites retained storage instead of constructing a second complete
ordering beside it.

## Radix-index to partition transition

During Morton sorting, `m_construction_index_scratch` is the alternate index
lane. Once ordering completes, that logical radix lifetime ends. Recursive tree
construction reuses the same physical vector for octant partition staging and
copies each node span back into the retained final permutation exactly as before.
No second N-sized partition vector exists.

## Lifecycle and generation rule

The first construction-workspace borrow marks the prior tree/order invalid and
clears logical node state while retaining useful capacities. New calls to the
node/ordering accessors reject an in-progress or failed rebuild. The next
`TreeBuildGeneration` is published only after Morton ordering, topology build,
and multipole accumulation complete successfully. If construction fails, the
previous generation value is not advanced and the partial rebuild is not exposed
as a valid tree.

Stable-count warm rebuilds resize to the same logical sizes and therefore reuse
retained vector capacities. Population growth may legitimately grow those
owners; no `shrink_to_fit()` phase churn is introduced.

## Memory model and observability

The gravity preflight estimator now charges:

- periodic tree coordinate staging as three doubles/source (24 B/source);
- tree construction ownership as final permutation + resolved epsilon + the
  three shared construction lanes (32 B/source).

The previous two periodic axis-workspace labels are removed from runtime memory
reporting. Runtime reporting now names the three physical construction owners as
`tree.construction.key_primary`, `tree.construction.key_scratch`, and
`tree.construction.index_scratch`, while the final permutation and source
softening remain separately reported. Logical periodic-sort, Morton-radix, and
partition roles are not double-charged.

At 512^3 sources, the model-derived retained reduction from 64 to 56 B/source is
approximately 1 GiB aggregate, and removing the old warm-build overlap from the
88 B/source model to 56 B/source is approximately 4 GiB aggregate. These are
architectural arithmetic results only; runtime RSS qualification is deferred.

## Reproducibility and numerical scope

The patch intentionally preserves the mathematical sequence and numerical
representation of:

- periodic wrapping and finite-coordinate validation;
- axis-specific box lengths;
- stable byte-radix ordering of finite non-negative wrapped doubles;
- largest-circular-gap selection, cyclic gap, equal-gap tie behavior, anchor,
  and unwrapping rule;
- tree bounds and extent floor;
- 21-bit normalized coordinate quantization, clamping, `std::llround`, and
  Morton interleave/bit placement;
- stable duplicate-key ordering;
- recursive octant partitioning;
- node topology, COMs, quadrupoles, softening envelopes, MAC, leaf policy,
  TreePM split, and force equations.

Storage ownership and lifetime change; empirical bitwise/runtime equivalence is
not claimed in this implementation pass.

## Interface migration note

The public gravity ordering helper changed from return-by-value
`buildMortonOrdering(...)` to allocation-free `buildMortonOrderingInPlace(...)`.
Callers must provide pre-sized authoritative ordering arrays plus the key and
index scratch spans. The production tree solver owns those buffers internally;
standalone `TreeGravitySolver::build(...)` callers require no new external
workspace object.

## Deferred work

M48-03 does not scalarize resolved source epsilon, compact tree nodes, unify
communication arenas, redesign migration/restart I/O, collapse homogeneous DMO
state, or alter PM spectral ownership. Tree-node capacity remains
geometry/distribution-sensitive, communication buffers remain independently
owned, and process-level RSS qualification remains pending.

Build/tests: NOT RUN — intentionally deferred by campaign instruction.
