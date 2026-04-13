# Infrastructure Diff Review Checklist

Use this checklist for repair-scoped PRs.

- [ ] Scope is infrastructure repair only (no new physics/numerics/science features).
- [ ] `core/` dependency direction is preserved (no upward include/use into analysis/io/physics/workflow without logged ADR whitelist).
- [ ] Diff does not introduce shadow abstractions or duplicate ownership across config, state, scheduler, I/O, diagnostics, or parallel layers.
- [ ] No second config system introduced; param-style UX still maps into typed validated config path.
- [ ] No string-literal policy branching added where typed enums/contracts already exist.
- [ ] New or changed config keys are validated, documented, normalized, and tested; renamed keys include compatibility behavior or migration notes.
- [ ] Snapshot/restart/provenance schema behavior unchanged, or change is explicit (versioning + compatibility + docs + migration notes).
- [ ] Reproducibility impact assessed: normalized config dump, provenance completeness, output naming, deterministic test mode, and restart continuity are preserved or explicitly updated.
- [ ] Claimed closure has command-backed test evidence, or blocker is explicitly reported with failing command.
- [ ] Patch adds targeted tests for the repaired invariant, or explicitly names the existing tests that already cover it.
- [ ] Public interface changes under `include/cosmosim/**` include same-patch docs updates and migration notes.
- [ ] Performance-sensitive changes do not introduce hidden full-state scans, cold-field contamination, or unnecessary copies; profiling/benchmark evidence is included when relevant.
- [ ] `docs/repair_state_recap.md`, `docs/repair_open_issues.md`, and `docs/architecture/decision_log.md` are updated when assumptions, blockers, or architecture decisions changed.
- [ ] Diff avoids cosmetic-only refactors unrelated to repair objective.
