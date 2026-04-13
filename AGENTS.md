# CosmoSim Agent Guardrails (Infrastructure Repair)

Scope: entire repository.

These rules are mandatory for infrastructure repair sessions and are enforced in code review.

## 1) Session scope: infrastructure repair only
- Allowed: build/test/preset repair, dependency wiring, interface/documentation alignment, schema/config safety hardening.
- Not allowed: new physics models, new numerics/solvers, or new science workflow features.
- Any diff that changes solver behavior must be rejected in repair-only sessions.

## 2) Architecture dependency direction
- `core/` is foundational and must not depend on `analysis/`, `io/`, `physics/`, or workflow/app layers.
- Any exception requires an architecture decision entry in `docs/architecture/decision_log.md` with explicit whitelist rationale before merge.

## 3) Configuration system discipline
- Do not introduce a second configuration system.
- Keep the existing param-style user UX and repair the typed, validated core path (`core::SimulationConfig` and normalization/provenance flow).
- String-literal policy checks are forbidden outside centralized config mapping layers when a typed contract already exists.

## 4) Snapshot/restart schema discipline
- No silent schema changes for snapshot/restart/provenance payloads.
- Any schema change must include: versioning update, compatibility behavior, tests, and docs updates (`docs/output_schema.md` + migration notes).

## 5) Evidence and closure discipline
- Do not claim issue/task closure without command-backed tests, or explicit blocker reporting that names the failing command and reason.
- CPU-only green status is not feature-path closure when HDF5/FFTW paths are in scope.

## 6) Patch hygiene for repair work
- No cosmetic-only refactors presented as infrastructure fixes.
- Keep patches minimal and directly tied to an open issue or documented guardrail.
- All public interface changes under `include/cosmosim/**` require same-patch docs updates and migration notes.

## 7) Required companion docs for contributors/reviewers
- Reviewer checklist: `docs/code_review.md`
- Architecture decisions: `docs/architecture/decision_log.md`
- Contributor workflow: `CONTRIBUTING.md`
## 8) Reproducibility and determinism discipline
- Do not weaken normalized config dumps, provenance completeness, output naming stability, or deterministic test modes.
- Any repair affecting config loading, snapshot/restart writing, scheduling, or rank coordination must state reproducibility impact explicitly.
- If deterministic behavior changes, document why, update tests, and record the tolerance or migration policy.

## 9) Interface and config-key anti-drift rules
- Do not introduce shadow interfaces that duplicate existing config, state, I/O, or scheduler responsibilities.
- No new config keys without typed validation, normalized dump support, docs updates, and tests.
- No renamed or repurposed config keys without compatibility behavior or explicit migration notes.
- New or changed public interfaces under `include/cosmosim/**` must define ownership, invariants, and caller-visible behavior.

## 10) Test-floor for infrastructure repair
- Infrastructure repairs must add targeted tests unless existing tests already cover the repaired invariant and that coverage is named explicitly.
- “Builds successfully” is not sufficient evidence for schema, config, provenance, restart, or interface repairs.
- TODO stubs, disabled assertions, or weakened checks do not count as closure.

## 11) Repair-state maintenance
- After any merged infrastructure repair, update `docs/repair_state_recap.md` and `docs/repair_open_issues.md` if the patch changes current repo assumptions, closes an issue, or changes remaining blockers.
- Architecture exceptions or dependency-boundary changes must be recorded in `docs/architecture/decision_log.md`.
