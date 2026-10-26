# Round 3 — correctness hardening and multirate stabilization (2026-10-26)

**Status: IMPLEMENTED — NOT COMPILED — NOT TESTED — OWNER VALIDATION REQUIRED.**

The owner forbade all execution for this round: no build, test, Python, MPI,
benchmark, lint or CI run. Nothing in this record has passed validation.
Branch `claude/vigilant-fermat-i8jvqw`, on top of the
[AREPO-grade campaign](arepo_grade_campaign_20261026.md).

## Issue A — explicit AREPO ParticleID translation (`tools/compare_dmo_arepo.py`)

**Option.** `--arepo-id-offset N` (default `0`, never inferred) defines
`normalized = original_arepo_id + N`.

**Checked uint64 arithmetic.** Python integers are used for the bounds check.
The tool rejects:

* overflow and underflow;
* negative IDs;
* float, bool and other non-integer IDs or offsets;
* duplicates on either side;
* non-identical normalized sets.

**Matching.** Matching is one-to-one, by sort order of the normalized IDs, and
the original arrays are not mutated. Halo matching uses the same normalized
membership.

**Diagnostics.** On a mismatch the tool reports:

* the counts;
* samples of IDs found on only one side;
* the ID ranges;
* the constant shift that *would* match, as a hint that is never applied.

**Output.** The JSON gains `identity_mapping`, with the convention, offset,
ranges and counts, and `provenance`, with the reader ID fingerprints. The
snapshot reader is unchanged.

**Tests.** `tests/python/test_compare_dmo_arepo.py`, registered as
`python_compare_dmo_arepo_identity` (exit 77 = SKIP without numpy or h5py).

## Issue B — frozen-timestep violation (`src/workflows/time_coordinator.cpp`)

**Defect.** At opening, the old multirate path froze the coarse interval for
the whole PM cycle. Inner blocks reused it even when the particle criteria
tightened, and the only response was a warning counter
(`…_quantum_clamped_blocks`). A block could therefore run with a quantum
coarser than the gravity/encounter criterion.

**Repair.**

* **Inner blocks.** Every block is re-planned by the reference per-block
  rules, so `quantum ≤ min_i dt_i` always holds.
* **Outer interval.** The outer PM interval `T` is fixed at opening as the
  minimum of:
  * the target block count;
  * the PM displacement criterion;
  * the output/end boundary (exact clip);
  * the retry bound.

  The inner blocks tile `T` exactly.
* **Transactional retry.** When the PM displacement criterion, re-checked at
  every interior boundary, rejects `T`, a collective rollback restores the
  governed cycle-opening snapshot. The source generation moves strictly
  forward, and the cycle is retried with `T ← min(T/2, Δt_disp)`. This
  terminates, because a single-block cycle cannot be rejected. The hard
  retry bound is 64, and exceeding it fails closed.
* **Endpoint rounding.** The closing interval is nudged by at most four ulps
  so that cycles close exactly on their planned end.
* **Counters.**
  * coarse and cycle limiters (gravity/particle, expansion, PM displacement,
    explicit, boundary, retry);
  * boundary shortening;
  * rejected/retried;
  * completed;
  * recovery failures;
  * unrepresentable quantum;
  * inexact cycle end;
  * soft-stop extension.

  The old clamp counter was removed, not renamed.

## Issue C — segmentation invariance

**Change.** The `final_step` truncation of the PM cycle is gone. The policy is
now *enforced cycle-boundary-only*:

* **Soft stop.** A budget reached inside a cycle runs on to the cycle end,
  which is reported as extension blocks. An open cycle at segment end is a
  `logic_error`.
* **Interior boundaries.** Checkpoints and output are deferred and analysis is
  skipped. Rebalance and migration are also skipped inside cycles, so
  ownership is frozen.
* **Restart schema.** None is needed. Every checkpoint is at a closed cycle,
  and `T` depends only on cycle-boundary state and physical criteria.

**Tests.** `test_hierarchical_treepm_workflow` (np 1/2/4) compares, at 1e-22
on phase space, with the same schedule and the same summed PM refreshes:

* one run of nominally 8 blocks with M = 4;
* the same run split 4+4, 2+2+2+2 and 1×8;
* budget 3 then the rest;
* budget 5 then the rest.

## Scientific contract

`docs/time_integration.md` (multirate section) states, against the code:

* which variables each operator modifies;
* the executed composition and its factor integrals;
* the time coordinate and scheduling formula;
* the enforced invariants and the recovery protocol;
* the interior PM usage;
* the authoritative epochs;
* the conditions for second order and reversibility, which adaptive steps
  break and which have not been measured;
* the stop/restart policy.

`exact` remains the reference.

## Section 7 dependency audit (source inspection only)

| Dependency | Finding |
| --- | --- |
| PM refresh counts | Refresh only at the cycle close (plus the segment bootstrap when no v24 endpoint exists). Profiler and cadence records include the work of rejected attempts, documented as diagnostic. Block records exclude it. |
| PM generation validity | The opening kick requires the opening field version/generation. A rollback advances the generation past every value the rejected attempt used and re-imports the endpoint under it. |
| Opening/closing kicks | `[a0, a_mid]` and `[a_mid, a_end(actual)]` with `a_mid = a(t0 + T/2)`: no gap or overlap. |
| Endpoint rounding | **Defect fixed**: the cycle-end ulp nudge. The reference path is unchanged. |
| Checkpoint deferral | Interior requests are latched to the cycle end. The soft stop guarantees the end is reached in the segment. |
| Safe-boundary flags | `multirate_pm_cycle_interior` is unchanged. Analysis hooks now also skip interior boundaries (**defect fixed**). |
| Tree-only bootstrap | Unchanged. |
| Momentum vs velocity | The predictor contract comment was wrong (it described a⁻² drift); **corrected** to `x += u ∫dt/a`. The hierarchical path stores `u = p/a`. No code change. |
| Scheduler after recovery | Both schedulers' persistent states are restored, then `syncTimeBinsFromSchedulers`. |
| MPI agreement | Rejection and `T` use globally reduced criteria. The rollback is wrapped in collective failure coordination. The tests assert identical block schedules across ranks. |
| Uninitialized fields | All new records and cycle state have default member initializers. |
| Remaining-step assumptions | Removed from cycle planning. The loop runs `while (step < final || cycle open) && t < t_end`. |

## Not done / limitations

* **No interior-cycle restart persistence.** The cycle-boundary-only policy was
  chosen instead. A future hard stop (wall-clock or signal) inside the
  hierarchical segment must defer to the cycle end in the same way.
* **Rollback snapshot memory.** It holds a full in-memory state copy, governed
  but roughly doubling particle memory during multi-block cycles.
* **Tree reuse/refit.** With these enabled, a rollback or restart can change
  the summation order. Agreement is then at roundoff level, and the tests use
  1e-18 for it.
* **No convergence order measured.** The tests gate monotone convergence and
  print the observed splitting order. PM-only and tree-only isolated
  convergence suites and an analytic two-body problem were **not written**.
* **No fixture forces a genuine displacement rejection.** The forced test seam
  exercises rollback completeness. Genuine rejection is exercised only if the
  dynamics trigger it.
* Gas/source multirate is still unsupported.
