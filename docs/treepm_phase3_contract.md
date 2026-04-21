# TreePM Phase 3 contract and maturity map (audit baseline)

Date: 2026-04-20 (UTC)

## Purpose and scope

This document defines what **Phase 3** means in the current repository, based on an explicit audit of the existing gravity stack. It is a contract for staged implementation prompts, not a claim that Phase 3 numerics are already implemented.

Scope of this document:

- inventory of what is already implemented and evidenced,
- explicit Phase 3 additions still required,
- forbidden completion claims before evidence exists,
- a hard closure gate with command-backed evidence,
- a stage map from requirements to concrete files/subsystems.

Out of scope in this stage:

- implementing new Phase 3 numerics,
- changing existing solver behavior,
- changing snapshot/restart schema semantics beyond current documented contracts.

## Current repository baseline (audited facts)

The following statements are true in the current tree and are treated as baseline constraints for Phase 3 planning:

1. **PM solver production contract is periodic** in the documented path (`solvePoissonPeriodic`, periodic zero-mode pinning, periodic deposition/gather semantics). This is the primary production PM contract described in repository docs and current solver API behavior.
2. **Distributed PM decomposition is slab-only** (`treepm_pm_decomposition_mode=slab`; non-slab modes are rejected in current config/validation paths).
3. **Axis-aware runtime support exists but scalar legacy/reporting and downstream consumers still need periodic cleanup** (workflow now constructs full `(Nx,Ny,Nz)` PM shapes, but legacy scalar aliases and some downstream surfaces remain for compatibility).
4. **Tree short-range path is monopole-only** with geometric MAC `l / r < theta`.
5. **Tree softening remains one-policy / one-epsilon in runtime path** (`comoving_fixed` + Plummer-style softened inverse-r^3 lane via one configured epsilon).
6. **Hierarchical stepping infrastructure exists**, including staged callbacks and time-bin scheduling; however, **integrator-grade PM synchronization/accuracy claims for a future multirate Phase 3 regime are not yet proven** by current evidence.
7. **Validation evidence is strong for Phase 1/2 baseline correctness contracts**, including MPI distributed Phase 2 agreement gates, but this evidence is **not sufficient to claim Phase 3 closure**.
8. **Current scaling artifacts are performance evidence only** and are explicitly not mature strong/weak scaling certification.

## Phase 3 scope (what must be added)

Phase 3 in this repository is defined as a maturity step from Phase 2 distributed correctness toward integrator-grade, evidence-backed production behavior.

Required additions (contract level):

1. **Integrator-grade cadence/synchronization contract**
   - Define and enforce PM long-range refresh/reuse correctness rules under hierarchical stepping beyond current cadence metadata consensus.
   - Include explicit invariants for when reused PM fields are valid relative to active-bin evolution and kick stages.

2. **Evidence-backed multirate accuracy envelope**
   - Add validation cases that detect unacceptable drift from cadence/reuse decisions under hierarchical activity patterns.
   - Provide tolerances and acceptance criteria distinct from Phase 2 distributed equivalence gates.

3. **Distributed maturity evidence expansion**
   - Extend beyond correctness-at-small-rank evidence to staged performance and reproducibility campaigns with explicit pass/fail criteria.
   - Keep the distinction between artifact generation and certification explicit.

4. **Contract-safe restart/provenance continuity for new Phase 3 lanes**
   - Any added cadence/synchronization state must follow existing schema-discipline requirements: versioning, compatibility behavior, tests, docs, migration notes.

5. **Operational claim discipline**
   - Keep public docs and release notes synchronized with actual command-backed evidence.
   - Do not promote Phase 3 language to “complete” until all hard-gate commands and artifacts pass.

## Forbidden claims before evidence exists

The following claims are forbidden until the hard gate in this document is fully satisfied:

- “Phase 3 is complete.”
- “Hierarchical TreePM multirate synchronization is production-proven.”
- “Strong/weak scaling is certified for mature production deployment.”
- “Restart continuation is Phase 3-complete” (unless any new schema/continuation lanes have matching versioning/tests/docs/migration updates).
- “Distributed decomposition is generalized beyond slab” (unless additional modes are implemented with typed config/docs/tests and ADR updates).

## Phase 3 closure hard gate (must all pass)

Phase 3 may be marked closed only when **all** items below are true in the same evidence cycle:

1. **Contract docs and ADR alignment**
   - `docs/treepm_phase3_contract.md` updated to final implemented behavior.
   - Same-patch ADR entry in `docs/architecture/decision_log.md` documenting scope, decision, and consequences.

2. **Validation gate coverage exists and passes**
   - Dedicated validation entries for Phase 3 cadence/synchronization invariants (not only Phase 2 equivalence).
   - Existing Phase 1/2 validation floors remain green.

3. **Restart/provenance contract integrity maintained**
   - If new continuation fields exist, schema/version updates + compatibility tests + docs/migration notes are included.
   - No silent payload drift.

4. **Scaling/reproducibility maturity evidence is explicit**
   - If certification language is used, it must be backed by declared strong/weak scaling campaign criteria and pass evidence.
   - Otherwise release/docs must keep non-certification language.

5. **Command-backed evidence bundle is published**
   - Configure/build/test commands and generated artifacts are listed verbatim in closeout docs.
   - Any blocker is reported with failing command and reason; no implicit closure language.

## Phase 3 staged maturity map (requirement -> files/subsystems)

| Stage | Requirement | Primary files/subsystems expected to change |
|---|---|---|
| P3-S0 (this patch) | Contract freeze and maturity map | `docs/treepm_phase3_contract.md`, `docs/architecture/decision_log.md`, `docs/validation_ladder.md`, `docs/validation_plan.md`, `docs/releases/known_issues.md` |
| P3-S1 | Hierarchical cadence + PM synchronization invariants (contract + instrumentation) | `include/cosmosim/core/time_integration.hpp`, `src/core/time_integration.cpp`, `src/workflows/reference_workflow.cpp`, `include/cosmosim/gravity/tree_pm_coupling.hpp`, `src/gravity/tree_pm_coupling.cpp` |
| P3-S2 | Phase 3 validation gates for multirate/cadence correctness | `tests/validation/test_validation_phase3_*.cpp` (new), `tests/validation/test_validation_phase2_mpi_gravity.cpp` (if shared harness is reused), `docs/validation_ladder.md`, `docs/validation_plan.md`, `validation/reference/*` |
| P3-S3 | Restart/provenance/schema continuity for any new cadence fields | `include/cosmosim/parallel/distributed_memory.hpp`, `src/parallel/distributed_memory.cpp`, restart/provenance I/O layers, `docs/output_schema.md`, migration notes, restart schema tests |
| P3-S4 | Maturity/scaling certification lane (only if claimed) | scaling artifact generators/CI targets, `docs/validation_ladder.md`, `docs/releases/known_issues.md`, closeout docs with explicit certification criteria |

## Reproducibility impact for this patch

- This patch is documentation-only.
- No runtime behavior, config key semantics, solver math, or schema payload bytes are changed.
- Deterministic outputs and provenance behavior are unchanged.

## What remains unimplemented after this stage

- Phase 3 multirate/integrator numerics and synchronization mechanisms.
- Phase 3-specific validation executable/tests and tolerance keys.
- Any expanded restart/provenance schema for additional Phase 3 continuation state.
- Any mature strong/weak scaling certification campaign and closure evidence.
