# Contributing to CosmoSim

CosmoSim is a research codebase with strict expectations for numerical integrity, reproducibility, and auditable interfaces.

## 0) Repo-local repair guardrails (read first)

- Mandatory repair-session rules: [`AGENTS.md`](AGENTS.md)
- Infrastructure review checklist: [`docs/code_review.md`](docs/code_review.md)
- Architecture decisions and exceptions: [`docs/architecture/decision_log.md`](docs/architecture/decision_log.md)

If your patch conflicts with these guardrails, resolve the conflict before opening a PR.

## 1) Required developer workflow

1. Choose the smallest preset that exercises your change.
2. Keep changes modular and ownership-respecting.
3. Run relevant unit + integration + validation checks.
4. If behavior or interfaces change, update docs in the same patch.
5. Include exact commands and outcomes in your PR.

See build matrix details in [`docs/build_instructions.md`](docs/build_instructions.md).

## 2) Interface and naming rules (authoritative)

- Files/directories: `lower_snake`.
- Classes/types: `PascalCase`.
- Methods: `camelCase`.
- Variables/parameters: `lower_snake`.
- Members: `m_lower_snake`.
- Use explicit unit/frame suffixes (`_code`, `_si`, `_cgs`, `_phys`, `_comov`) when ambiguous.

Public API ownership remains under `include/cosmosim/<module>/...`; internal helpers remain under `src/<module>/internal/...`.

## 3) Change-control rules

You **must** document these in the same patch when changed:

- config keys or parsing behavior,
- snapshot/restart/provenance schema fields,
- mode-policy behavior (frame/boundary/units conventions),
- benchmark/profiling output contract used by developers.

Required docs for those updates:

- config changes → `docs/configuration.md`
- output schema changes → `docs/output_schema.md`
- validation changes → `docs/validation_plan.md`
- benchmark/profiling workflow changes → `docs/profiling.md`
- architecture-level decisions → `docs/architecture/decision_log.md`

## 4) Codex prompt + review expectations

Every implementation PR should satisfy the workflow contract in [`docs/architecture/developer_workflow_contract.md`](docs/architecture/developer_workflow_contract.md):

- explicit assumptions and invariants,
- no hidden schema/config drift,
- no placeholder core logic,
- exact changed-file listing and rationale,
- test and benchmark/profiling evidence.

## 5) Local checks before opening a PR

```bash
./scripts/ci/check_repo_hygiene.sh
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

If your patch touches HDF5/FFTW/MPI/CUDA/Python paths, also run the matching dependency-enabled preset(s) from `docs/build_instructions.md`.
