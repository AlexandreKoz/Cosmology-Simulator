# Hierarchical Timestepping Audit (2026-04-06)

## Scope reviewed
- Scheduler implementation and bin legality/synchronization paths:
  - `include/cosmosim/core/time_integration.hpp`
  - `src/core/time_integration.cpp`
- Data layout and ownership contract:
  - `include/cosmosim/core/simulation_state.hpp`
  - `docs/time_integration.md`
  - `docs/state_model_memory_layout.md`
- Representative tests:
  - `tests/unit/test_time_integration.cpp`
  - `tests/integration/test_hierarchical_time_bins.cpp`
  - `tests/integration/test_hierarchical_timestep_regression.cpp`

## 1) Pass/fail summary

**Overall: PARTIAL PASS**

- **PASS:** Bin hierarchy is discrete power-of-two and transitions are guarded by synchronization checks (`tick % period == 0`).
- **PASS:** Gravity/CFL/source/user hooks are explicit and testable through `TimeStepCriteriaRegistry` and `combineTimeStepCriteria`.
- **PASS (qualified):** Active set is not discovered by scanning contiguous full-state SoA arrays; scheduler iterates compact per-bin ownership vectors.
- **FAIL (diagnostics consistency):** `TimeBinDiagnostics` exposes `clipped_to_min_dt` and `clipped_to_max_dt`, but runtime code never increments those counters.
- **FAIL (single-source ownership):** `ParticleSoa`/`CellSoa` retain `time_bin` while scheduler has authoritative `TimeBinHotMetadata`, leaving dual state representations.

## 2) Synchronization or legality failures

### Verified safety mechanisms
- Bin periods are strict powers of two via `1ULL << bin`.
- Pending transitions are only applied for currently active elements.
- Transition is rejected and counted as illegal when destination period does not align with current tick.

### Findings
- **No hard synchronization bug found** in current implementation path.
- **Legality observability gap:** `clampBin()` silently coerces out-of-range requests to `max_bin`; this prevents crashes but masks caller policy/config errors.

## 3) Performance/locality risks from active-set strategy

### Current strengths
- Membership is stored as compact per-bin vectors (`m_elements_by_bin`) with O(1) swap-erase updates.
- Scheduler does not sweep particle/cell SoA arrays to build active sets.

### Risks
- `rebuildActiveSet()` still loops over every bin and every member each substep; this becomes O(total elements) even when few bins are tick-active.
- `std::sort(m_active_elements)` adds O(A log A) work every substep.
- `occupancy_by_bin` is recomputed from scratch each substep, despite occupancy only changing on transition/assignment events.

## 4) Missing diagnostics for timestep pathologies

- **Unwired counters:** `clipped_to_min_dt` / `clipped_to_max_dt` are dead fields from the scheduler’s perspective.
- **Missing coercion visibility:** No count for requests clamped by `max_bin`.
- **Missing backlog visibility:** No metric for age/volume of pending transitions that cannot be applied yet.
- **Missing per-element pathology breadcrumbs:** No compact event stream for repeated illegal transitions or repeated min-dt clipping on the same element.

## 5) Root causes (why these gaps exist)

1. **Split diagnostic ownership without integration seam**
   - `mapDtToTimeBin` is a pure helper returning flags; scheduler diagnostics are mutable state inside `HierarchicalTimeBinScheduler`.
   - There is no shared aggregator API or callback that bridges mapping flags into scheduler counters.

2. **“Rebuild each substep” design choice favored simplicity over asymptotic cost**
   - The scheduler stores enough metadata for incremental accounting, but current implementation intentionally recomputes occupancy and active counts in one pass each substep for determinism and straightforward invariants.

3. **Legacy hot-state field retention during scheduler introduction**
   - `time_bin` remains in particle/cell SoA from earlier layout decisions while a newer scheduler-side hot metadata block was added.
   - Without an explicit deprecation/mirroring policy, both can drift conceptually into dual authority.

4. **Safety-by-clamp policy with no audit trail**
   - Illegal/out-of-range bin requests are normalized (`min(requested, max_bin)`) to keep runtime robust.
   - Because this path does not emit diagnostics, policy mistakes appear as valid behavior downstream.

5. **Transition application intentionally tied to activity**
   - Applying transitions only for active elements preserves synchronization safety.
   - However, this also means there is no built-in visibility into deferred/starved requests, so long-lived pathologies can remain opaque.

## 6) Minimum changes needed for safety/transparency

1. **Wire clipping diagnostics:** introduce a tiny reporting seam (e.g., `recordDtMapping(result)`) at dt-to-bin call sites to update scheduler/global counters.
2. **Add coercion diagnostics:** increment `clamped_transition_requests` whenever `requested_bin > max_bin` before clamping.
3. **Incremental occupancy accounting:** maintain `occupancy_by_bin` in `eraseFromBin`/`insertIntoBin`; keep substep pass focused on tick-active bins only.
4. **Choose one authoritative time-bin owner:** deprecate SoA `time_bin` fields or enforce strict synchronization via explicit invariant checks/tests.
5. **Add pathology-focused tests/metrics:** pending-transition age histogram, repeated-illegal-attempt counters, and clipping hot-element tracking.

## Validation run used for this audit
- Configure/build: `cmake --preset cpu-only-debug && cmake --build --preset build-cpu-debug -j4`
- Tests: `ctest --test-dir build/cpu-only-debug --output-on-failure`
