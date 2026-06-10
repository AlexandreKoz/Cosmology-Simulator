# CHUI / CosmoSim Agent Contract

Scope: entire repository.

This is the canonical instruction surface for AI coding agents and human reviewers working in this repository. The user-facing project name is **CHUI**. Existing namespaces, file stems, schema names, and compatibility identifiers may still say `cosmosim`/`CosmoSim`; do not rename them opportunistically. Treat naming migration as a dedicated compatibility project, not as drive-by cleanup.

## 1) First principles

- Build a serious computational-astrophysics code, not a demo scaffold.
- Preserve scientific correctness over apparent feature velocity.
- Prefer small, auditable patches with explicit ownership boundaries.
- Keep workstation and small-cluster efficiency as first-class constraints.
- Do not claim completion without command-backed evidence or a named blocker.
- Never hide uncertainty: if a dependency, preset, or environment blocks validation, report the exact command and reason.

## 2) Session-mode classification

Before editing, classify the request into exactly one mode and obey its constraints.

### Audit mode

Use when the request asks for review, readiness, scientific accuracy, maintainability, or implementation status.

- Allowed: code/document inspection, evidence mapping, risk ranking, concrete next-step recommendations.
- Not allowed: modifying code unless the user explicitly asks for a patch.
- Required output: verdict, blockers, evidence references, and a prioritized fix plan.

### Repair mode

Use when the request asks to fix blockers, harden interfaces, repair CI/build/config/restart/schema/state ownership, or close known audit gaps.

- Allowed: minimal code/docs/tests tied directly to the defect.
- Not allowed: new physics models, solver rewrites, speculative abstractions, or unrelated cleanup.
- Any solver-behavior change must be explicitly justified as a bug fix and covered by tests.

### Feature implementation mode

Use when the request explicitly asks for a new capability or stage implementation.

- Allowed: production code, docs, tests, and benchmark/profiling hooks required by the feature.
- Not allowed: pseudocode, TODO-core logic, hollow facades, or claims unsupported by validation.
- New physics/numerics must include assumptions, equations or update laws, invariants, and validation targets.

### Packaging mode

Use when the request asks for a clean zip or handoff artifact.

- Return source-only content unless the user explicitly requests build artifacts.
- Exclude build trees, object files, local caches, generated binaries, temporary output, and user-local presets.
- Keep tracked examples, configs, validation decks, docs, and source tests.

## 3) Mandatory orientation pass

Before a non-trivial patch, inspect the relevant subset of:

- `README.md` for repo map and current maturity boundary.
- `CONTRIBUTING.md` for workflow and review expectations.
- `docs/architecture/overview.md` for module ownership.
- `docs/architecture/developer_workflow_contract.md` for patch-response requirements.
- `docs/architecture/runtime_truth_map.md` for runtime state ownership.
- `docs/architecture/adr_runtime_truth_ownership.md` for single-source-of-truth policy.
- `docs/repair_state_recap.md` and `docs/repair_open_issues.md` for known blockers.
- Module-specific docs before touching the matching module.

Do not edit a subsystem before identifying its owner, public interface, persistence/schema impact, and validation floor.

## 4) Architecture dependency direction

- `core` is foundational and must not depend on `analysis`, `io`, `physics`, workflow/app layers, or module-specific implementation details.
- `utils` must remain narrow support code; do not dump solver or physics logic there.
- Public interfaces belong under `include/cosmosim/<module>/`.
- Implementations belong under `src/<module>/`.
- Internal-only helpers belong under `src/<module>/internal/`.
- Any dependency-boundary exception requires an entry in `docs/architecture/decision_log.md` before merge.

## 5) Runtime truth and ownership discipline

- One live authority per runtime fact.
- Mirrors, caches, sidecars, active views, restart payload fields, and diagnostics must define refresh/invalidation rules.
- A migrated particle or cell must carry every required authoritative and sidecar field, and stale ghosts must never mutate truth.
- Scheduler-owned timestep state must remain authoritative; particle/cell timestep lanes are mirrors unless a documented contract says otherwise.
- PM slabs/pencils, hydro ghost cells, AMR patches, tree pseudo-particles, and imported ghosts must define owner, mutation rights, epoch/version, and rebuild boundary.

## 6) Configuration, schema, and reproducibility gates

- Do not introduce a second configuration system.
- Keep the existing param-style user UX mapped into the typed, validated config path.
- No new config key without typed validation, normalized dump support, docs updates, and tests.
- No renamed or repurposed config key without compatibility behavior or explicit migration notes.
- No silent snapshot/restart/provenance schema changes.
- Any schema change requires versioning, compatibility behavior, tests, and docs updates in `docs/output_schema.md` and/or `docs/restart_checkpointing.md`.
- Any change affecting config loading, snapshot/restart writing, scheduling, or rank coordination must state reproducibility impact.

## 7) Scientific and numerical bar

- Do not add baryonic/subgrid physics as isolated toggles without a coherent model contract.
- Hydrodynamics changes must preserve conservative finite-volume accounting or document exactly where source terms enter.
- AMR changes must define refinement ownership, prolongation/restriction behavior, synchronization, and reflux/conservation implications.
- Gravity/TreePM changes must identify force split, boundary behavior, softening convention, and parallel ownership impact.
- Distributed-memory changes must consider particle count, active fraction, gas/AMR cost, tree/PM load, memory pressure, and future GPU occupancy when load balance is affected.
- Performance changes must protect hot/cold separation, SoA locality, compact active-set paths, and sidecar discipline.

## 8) Naming and file-placement rules

- Files and directories are `lower_snake` except externally mandated files such as `AGENTS.md`, `CMakeLists.txt`, `CMakePresets.json`, and `.github/copilot-instructions.md`.
- C++ types are `PascalCase`.
- Functions and methods are `camelCase`.
- Variables and parameters are `lower_snake`.
- Class/struct data members use `m_lower_snake`.
- Use `_code`, `_si`, `_cgs`, `_phys`, `_comov`, `_x`, `_y`, `_z`, `_begin`, and `_end` when ambiguity exists.
- Preserve canonical GADGET/AREPO-style HDF5 group and dataset names exactly.

## 9) Patch hygiene

- No pseudocode in production paths.
- No TODO stubs in core logic.
- No broad formatting churn unless requested.
- No cosmetic-only refactors presented as fixes.
- No hidden full-state scans, cold-field pollution, unnecessary copies, or accidental hot-path allocations.
- No public interface change under `include/cosmosim/**` without same-patch docs and migration notes.
- Do not weaken assertions, tests, deterministic modes, normalized config dumps, provenance completeness, output naming, or restart validation to make tests pass.

## 10) Required evidence floor

Choose the smallest evidence set that exercises the changed invariant, then report exact commands.

Baseline for docs/interface-only patches:

```bash
./scripts/ci/check_repo_hygiene.sh
```

Baseline for CPU code patches:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

If the patch touches HDF5, FFTW, MPI, CUDA, Python bindings, or validation decks, also run the matching preset(s) from `docs/build_instructions.md`. If unavailable, report the blocked command and missing dependency.

## 11) Documentation coupling

Update documentation in the same patch when behavior changes:

- User workflow: `README.md`, `docs/build_instructions.md`.
- Config keys or semantics: `docs/configuration.md`.
- Snapshot/restart/provenance schema: `docs/output_schema.md`, `docs/restart_checkpointing.md`.
- Validation expectations: `docs/validation_plan.md`.
- Profiling/benchmark workflow: `docs/profiling.md`.
- Architecture decisions or dependency exceptions: `docs/architecture/decision_log.md`.
- Current blocker/assumption state: `docs/repair_state_recap.md`, `docs/repair_open_issues.md`.

## 12) Agent response contract

For implementation or repair work, the final response must include:

1. What changed.
2. Why it changed.
3. Exact files touched.
4. Commands run and outcomes, or exact blockers.
5. Any remaining limitations.
6. Link to the clean artifact when packaging was requested.

For audit work, include the evidence basis and distinguish confirmed defects from risk hypotheses.

## 13) Clean artifact rules

Before creating a zip, remove or exclude:

- `build/`, `build-*`, `cmake-build-*`, `CMakeFiles/`, `CMakeCache.txt`, `CTestTestfile.cmake`, `Testing/`.
- Object/static/shared libraries, executables, generated coverage, and profiler traces.
- Python caches: `__pycache__/`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`, `.ipynb_checkpoints/`.
- Runtime outputs: `output/`, temporary scratch, logs generated during validation.
- User-local presets: `CMakeUserPresets.json`.

Do not remove source files, docs, configs, tests, validation decks, benchmark sources, or examples that are meant to be tracked.
