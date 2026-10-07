# P6.1 TreePM hot-path instrumentation cleanup — 2026-10-07

**SOURCE IMPLEMENTED**

**BUILD NOT RUN**

**TESTS NOT RUN**

**BENCHMARKS NOT RUN**

**VALIDATION PENDING**

Mode: repair. Input revision: `72c6ce6da5e4097459c7c3343306f4c1d04acf8b`.
Only the two post-P6 instrumentation issues and their directly required memory
accounting/documentation were changed. The P6 qualification campaign remains
pending; no compilation, simulation, test suite or benchmark was initiated.

## Timing ownership

`TreePmCoordinator::evaluateShortRangeResidual` no longer reads the clock in
`evaluateTargetAgainstLocalTree`. Each existing logical block of at most 64
targets owns one start/stop pair in non-distributed local traversal,
distributed local-overlap traversal and incoming-target traversal. This also
applies to the existing serial fallback.

The executing worker sums block duration into `elapsed_work_ms` and retains
`block_work_ms_max`. After join, local/incoming counter merges sum elapsed
work and retain the largest block duration; existing worker-region summaries
remain available. Block timing includes target preparation, traversal,
force-slot writes and local diagnostic accumulation, excluding scheduling
waits and MPI. There is no retained block table and no per-target timer.
Instrumentation clock calls scale with blocks rather than targets. This is a
structural overhead reduction, not a measured runtime improvement.

## Hot counters and optional spatial history

`TreePmTraversalCounters` retains `alignas(64)` and scalar targets, visits,
opens, direct pairs, accepted internal multipoles/leaves, rejection counters
and work timing. Its two histogram arrays move into a separate private
`alignas(64)` worker scratch type, with independent capacity/high-water reporting.
Scratch is prepared once at the existing outer workspace scope only when
`spatial_work_history_enabled` is true, before any worker region or distributed
request posting. No per-target allocation is introduced.

Owned-target spatial work accumulates across local batches; incoming targets
receive no spatial scratch. A fixed worker/bin-order integer reduction at
solve end precedes the unchanged MPI sums. The fixed 64-bin SFC mapping,
visits + direct pairs work measure, target counts, 0.5 decay, work/target rate
and P5-W weighting mathematics are unchanged. Default-disabled execution
does not allocate, clear, update or reduce histograms. If a direct API caller
previously enabled the feature, retained capacity remains separately accounted.
The existing MemoryGovernor preflight conservatively budgets the separated
optional width alongside the actual hot-counter size; no governor or policy
interface is added.

## Preservation and interface impact

No TreePM force mathematics, strict 0.08 envelope, adaptive/relative MAC,
maximum-angle policy, softening/cutoff/near-node guard, pair accounting, target
or child traversal order, OpenMP scheduling, QueryScratch, tree reuse/refit,
KDK, PM algorithm or MPI protocol changed. All P6 experimental configuration
defaults remain unchanged. Force arithmetic and integer work accounting remain
in the existing order by source inspection; runtime equivalence is untested.

Native counter consumers must rebuild after removal of the raw histogram
members and addition of `block_work_ms_max`; planner history remains exposed
through `TreePmDiagnostics::spatial_work_per_target`. Existing worker timing
fields now summarize block-duration sums. No serialized event fields or
snapshot/restart/provenance schemas change. Diagnostic timing remains
nondeterministic and does not feed force arithmetic or P5-W weighting.

## Files and static evidence

- `include/cosmosim/gravity/tree_pm_coupling.hpp`: aligned scalar counters,
  block maximum and separate optional worker scratch.
- `src/gravity/tree_pm_coupling.cpp`: block clocks, optional scratch lifecycle,
  deterministic spatial reduction and separate capacity reporting.
- `src/gravity/gravity_memory.cpp`: preserve optional scratch coverage in the
  existing conservative worker estimate.
- `CURRENT_STATUS.md`: concise P6.1 status.
- `docs/tree_pm_coupling.md`: timing/scratch ownership and native API migration.
- `docs/profiling.md`: revised work-time interpretation and tail diagnostics.
- `docs/memory_governance.md`: split worker owners and preflight coverage.
- `docs/repair/p6_preproduction_performance_architecture_implementation_20261006.md`:
  narrow P6 follow-up status.
- This note: scope, limitations and changed-file inventory.

Static commands: `git status --short`, `git diff --stat`, `git diff --check`,
targeted `git diff -- <paths>`, `git rev-parse HEAD`, `rg` / `rg --files`,
`cat` / `sed` / `wc` source-document inspection, and Python standard-library
source comparison and ZIP creation/inventory. New-note whitespace is checked
with `git diff --no-index --check /dev/null docs/repair/p6_1_treepm_hotpath_instrumentation_cleanup_20261007.md` (expected
exit 1 for the new-file difference; no whitespace diagnostics).

`git diff --check` returned 0. Static Python comparison against the input Git
blob confirmed the node/pair traversal body is byte-identical, all OpenMP
directives and the spatial EMA update are unchanged, the target evaluator has
zero clock reads, and each of the three block loops has one start/stop pair.
These are source-inspection results only.

Compilation: **NOT RUN**. Tests: **NOT RUN**. Benchmarks: **NOT RUN**.
Scientific validation: **NOT RUN**. Static inspection does not qualify runtime
behavior, numerical equivalence or performance. The existing
[P6 qualification plan](p6_validation_campaign_plan_20261006.md) is unchanged.

Delivery: `chui_p6_1_treepm_hotpath_cleanup_20261007.zip`, changed/new files only,
with repository-relative paths. Unrelated local `opencode.jsonc` and
`session-ses_efe3.md` are preserved and excluded.
