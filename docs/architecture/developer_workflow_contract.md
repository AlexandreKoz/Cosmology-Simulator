# Developer workflow contract (Codex + human review)

This contract defines the minimum acceptable quality bar for implementation prompts and PR review in CosmoSim.

## Required response structure for implementation prompts

1. Brief design summary.
2. Exact file list with one-line rationale per file.
3. Explicit assumptions, invariants, and numerical conventions.
4. Full changed-file content.
5. Tests added/updated.
6. Benchmark or profiling hook added/updated.
7. Acceptance checklist.

## Mandatory engineering rules

- No pseudocode or TODO stubs in core logic.
- No unrelated rewrites.
- No hidden interface drift.
- No silent schema/config/provenance changes.
- No ambiguous unit/frame naming where confusion is possible.

## Reviewer checklist

- Are assumptions and limitations explicit?
- Are config/schema/restart/provenance implications documented?
- Do tests cover local invariants and a broader pipeline path?
- Is there at least one benchmark/profiling hook for costly paths?
- Are naming and ownership rules respected?

## Documentation coupling rules

Any PR that changes behavior must update docs in-repo in the same patch:

- user-facing workflow updates (`README.md`, `docs/build_instructions.md`)
- configuration keys or validation semantics (`docs/configuration.md`)
- snapshot/restart/provenance schema behavior (`docs/output_schema.md`)
- validation expectations (`docs/validation_plan.md`)
- profiling/benchmark usage (`docs/profiling.md`)
- architecture decisions (`docs/architecture/decision_log.md`)
