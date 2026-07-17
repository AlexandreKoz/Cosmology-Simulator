# CHUÍ Runtime Decomposition Goal Specification

**Experiment:** Campaign A — runtime decomposition only  
**Suggested PR title:** `refactor-runtime-workflow-into-owner-services-and-typed-resource-views`  
**Suggested local branch:** `refactor/runtime-owner-services`  
**Primary audit findings:** BRAIN-006, BRAIN-007, BRAIN-015, BRAIN-016, and the orchestration portion of BRAIN-018

---

## 1. Experiment hypothesis

A persistent Codex Goal with one measurable architectural objective will produce a materially better result than a single mega-prompt containing several research-scale projects.

The objective of this experiment is only:

> Transform CHUÍ's rung-zero production workflow from an 8,008-line policy monolith into a small composition/lifecycle root backed by owner services, enforceable stage resource views, and descriptor-driven module registration, while preserving current numerical behavior and the existing passing test matrix.

This campaign deliberately does **not** implement hierarchical KDK, distributed IC ingestion, hydro wake-up, elastic restart, asynchronous output, or new physics. Those are later Goals that should build on the decomposed runtime.

---

## 2. Repository preparation

Work in a dedicated local branch or isolated Codex worktree.

Preserve the current hardening patch. Do not start from the older baseline.

Put the adversarial audit in the repository at:

```text
docs/repair/chui_runtime_brain_prompt_result_adversarial_audit.md
```

Put this experiment specification in the repository at:

```text
docs/repair/runtime_decomposition_goal_spec.md
```

Read `AGENTS.md` before doing anything else.

Do not use destructive Git operations. Preserve all pre-existing local changes.

---

## 3. First command: planning only

Start a fresh Codex thread rooted at the local CHUÍ repository. Select the highest available reasoning level.

Run `/plan` with the following text:

```text
/plan

Read AGENTS.md, docs/repair/chui_runtime_brain_prompt_result_adversarial_audit.md,
and docs/repair/runtime_decomposition_goal_spec.md. Inspect the current repository,
especially src/workflows/reference_workflow.cpp, time integration, runtime services,
module registry, IC startup, gravity, hydro/AMR, migration, output/restart, and their
tests.

Produce a repository-grounded execution plan for Campaign A only: transform the
current 8,008-line production workflow into a small composition/lifecycle root backed
by narrow owner services, enforceable stage resource views, and descriptor-driven
module registration, while preserving all existing rung-zero behavior.

Do not edit files yet. The plan must identify exact extraction seams, ownership and
dependency direction, migration order, temporary compatibility adapters, tests after
each milestone, structural acceptance checks, likely risks, and rollback-safe
checkpoints. It must explicitly avoid implementing hierarchical KDK, distributed IC
reading, hydro wake-up, elastic restart, asynchronous output, or unrelated physics.
```

Review the plan before activating the Goal.

Reject or ask Codex to revise the plan if it:

- proposes one replacement god-file;
- treats line movement alone as decomposition;
- leaves all callbacks with unrestricted mutable `SimulationState`;
- keeps module descriptors as descriptive strings only;
- attempts hierarchical KDK or distributed IC work;
- lacks tests after each extraction milestone;
- does not specify how rung-zero numerical behavior will be compared;
- does not identify a dependency-safe extraction order.

---

## 4. Persistent Goal command

After the plan is acceptable, activate this Goal in the same thread:

```text
/goal Complete the runtime-decomposition objective in
docs/repair/runtime_decomposition_goal_spec.md and the approved plan without
stopping until every mandatory acceptance gate in that specification is satisfied
by source evidence and passing tests. Preserve the current rung-zero numerical
behavior, public configuration semantics, restart compatibility, state ownership,
and all pre-existing user changes. Work milestone by milestone, update the living
execution log, run focused validation after every milestone, repair failures before
continuing, and perform an adversarial diff review before completion. Do not
implement hierarchical KDK, distributed IC ingestion, hydro wake-up, elastic
restart, asynchronous output, new physics, or unrelated cleanup. If no defensible
path remains under the local environment or budget, stop as blocked rather than
complete and report the exact unmet gates, attempted approaches, evidence, and the
minimum input needed to resume.
```

A Goal is not complete because the code compiles or because classes with the desired
names exist. Completion requires every acceptance gate below.

---

## 5. Required target architecture

Use repository-native names where existing conventions provide better terminology.
The responsibilities must nevertheless become equivalent to the following.

### 5.1 Composition and lifecycle root

`ReferenceWorkflowRunner` and `reference_workflow.cpp` may own only:

- top-level run lifecycle;
- construction of immutable/frozen configuration;
- construction of shared runtime services;
- selection and construction of registered modules/services;
- high-level call into the time coordinator;
- final summary and top-level error reporting.

It must not contain numerical gravity, hydro, AMR, IC parsing/conversion, migration
serialization, snapshot field comparison, restart payload verification, or detailed
timestep-criteria implementations.

### 5.2 Owner services

Extract narrow production owners equivalent to:

- `InitialConditionRuntime`
  - canonical startup path;
  - IC manifest invocation;
  - initial state construction and initial decomposition handoff;
  - no new distributed IC implementation in this campaign.

- `TimeCoordinator`
  - rung-zero KDK lifecycle;
  - scheduler and active-set ownership;
  - timestep-criteria plan coordination;
  - endpoint and code-time event clipping;
  - legal output/restart boundaries.

- `GravityRuntime`
  - TreePM/PM cadence orchestration;
  - force-cache lifecycle;
  - gravity task/stage contributions;
  - gravity-specific validation and diagnostics.

- `HydroAmrRuntime`
  - hydro and currently supported AMR orchestration;
  - ghost/patch freshness requirements;
  - hydro/AMR stage contributions;
  - no local-timestep wake-up implementation in this campaign.

- `SourceRuntime`
  - cooling, star formation, black-hole, tracer, and existing source orchestration;
  - source-specific timestep criteria registration;
  - spawned-state handoff through existing supported rung-zero semantics.

- `MigrationBalanceRuntime`
  - decomposition, migration, compaction, scheduler remap, sidecar transport,
    rebalance decisions, and generation invalidation.

- `OutputRestartRuntime`
  - snapshot/checkpoint writing;
  - field-aware verification;
  - restart schema/readback validation;
  - output cadence state and provenance.

- `FailureCoordinator`
  - shared failure-consensus primitive and top-level phase boundaries already
    supportable without redesigning MPI execution.

Services may be divided further where existing repository seams justify it. Do not
combine all extracted code into one `runtime_brain.cpp`.

### 5.3 Runtime services

Expand or preserve `RuntimeServices` as the shared injected dependency bundle. It
must be constructed in one composition root and passed explicitly. It should expose
only dependencies genuinely shared by multiple owners, such as:

- root `MpiContext` and topology references;
- profiler/diagnostics;
- deterministic policy;
- persistent runtime workspaces or memory resources where appropriate;
- logging/failure facilities;
- execution backend metadata already present.

Do not introduce a service locator, hidden global singleton, or new mutable global
state.

### 5.4 Enforceable resource views

Production stage/task callbacks must no longer receive unrestricted mutable access
to all of `SimulationState` merely for convenience.

Introduce typed views or capability objects appropriate to the actual stages, such
as:

- drift particle write view;
- gravity source read view;
- acceleration write view;
- hydro conserved-state view;
- AMR patch/state view;
- source mutation view;
- migration ownership view;
- output/restart read view.

The exact type set should be minimal and repository-driven.

Each view must:

- expose only authorized state lanes;
- distinguish read-only from writable resources;
- carry or validate state generation and scheduler tick/epoch metadata where
  staleness is possible;
- preserve stable-ID and sidecar alignment;
- fail loudly in debug tests when used after invalidation;
- document authority, derivation, invalidation events, and restart policy.

Temporary compatibility adapters are allowed only during staged migration. No
unrestricted production callback path may remain at Goal completion.

### 5.5 Descriptor-driven module composition

Replace the current descriptive string table with typed module descriptors that
actually drive runtime contributions.

A descriptor must be able to declare and register, as applicable:

- module identity and capability prerequisites;
- service/factory contribution;
- stages/tasks;
- required resource views;
- timestep criteria;
- restart payload hooks;
- migration fields/hooks;
- diagnostics;
- configuration validation already present in the central typed config.

The composition root must iterate registered descriptors/factories rather than
manually knowing every callback class.

Add a test-only probe module that is registered through the descriptor mechanism
without editing `ReferenceWorkflowRunner` or a central callback list. The probe must
contribute at least one observable stage/diagnostic action and be verified by a
test. A table of non-empty strings is not sufficient.

---

## 6. Mandatory structural acceptance gates

The Goal is incomplete unless all gates pass.

### GATE-A — workflow reduction

After excluding comments and generated code:

- `src/workflows/reference_workflow.cpp` should target no more than **1,500 lines**;
- an evidence-backed exception may extend the hard limit to **2,000 lines**, but the
  final report must justify every remaining responsibility;
- `ReferenceWorkflowRunner::runImpl()` must be a high-level lifecycle method rather
  than a numerical algorithm;
- no single replacement runtime source may exceed **2,000 lines**;
- no extracted service may simply contain the old callback body unchanged as one
  multi-thousand-line function.

Line counts are guardrails, not the only proof. Passing them by dumping code into a
different monolith fails this gate.

### GATE-B — dependency direction

The component graph must remain acyclic.

- workflow/composition may depend on owner-service interfaces;
- owner services may depend on core and their numerical modules;
- numerical modules must not depend upward on workflow composition;
- `core` must not depend on physics/workflows;
- no owner service may include another owner's private implementation header.

Extend the existing architecture/hygiene checks to enforce the new boundaries.

### GATE-C — production resource restriction

All production callbacks/tasks use typed resource views or narrow owner APIs.

Tests must prove:

- stale generation/tick views are rejected;
- read-only views cannot mutate forbidden lanes at compile time where feasible;
- a stage cannot access unrelated sidecars through its public production interface;
- migration/compaction invalidates prior views.

### GATE-D — real descriptor registration

The test-only probe module must register through the descriptor/factory mechanism
without editing the composition root.

Its stage/diagnostic contribution must execute in a focused integration test.

Removing the probe descriptor from the registry must remove the contribution without
other workflow edits.

### GATE-E — behavior preservation

The refactor must not change accepted rung-zero results.

At minimum:

- all pre-existing CPU-debug tests pass;
- the real reference-workflow integration tests pass;
- focused HDF5 IC, snapshot, and restart tests pass;
- restart direct-vs-resumed behavior remains within its existing exact/tolerance
  policy;
- deterministic ordering/provenance tests remain valid;
- no public configuration key is silently renamed or reinterpreted;
- restart schema changes are avoided unless unavoidable and explicitly justified.

Where practical, add a golden pre/post orchestration equivalence test comparing
stable-ID-keyed particle/state summaries and stage ordering.

### GATE-F — no scope escape

The final diff must not:

- enable nonzero hierarchical rungs;
- claim adaptive production stepping;
- implement or claim multifile/distributed IC ingestion;
- implement hydro wake-up;
- implement rank-remappable restart;
- implement asynchronous output;
- add new physics;
- weaken tolerances, assertions, validation, or timeouts;
- perform unrelated bulk formatting.

The capability report must remain truthful.

### GATE-G — maintainability evidence

The final closeout must contain:

- before/after line counts and largest runtime units;
- service responsibility table;
- dependency graph result;
- list of former workflow responsibilities and their new owners;
- list of remaining workflow responsibilities with justification;
- proof that a module contribution can be added through descriptors;
- proof that production resource access is narrower;
- exact test commands and outcomes;
- unresolved risks.

---

## 7. Milestone sequence

Codex may refine this sequence during `/plan`, but it must preserve dependency safety.

### Milestone 0 — baseline and living log

Create:

```text
docs/repair/runtime_decomposition_exec_plan.md
```

Record:

- current branch and dirty state;
- baseline line counts;
- current service/callback map;
- current test commands and outcomes;
- milestone status;
- design decisions and rejected alternatives;
- exact files changed;
- blockers.

Run the focused baseline and record it before edits.

### Milestone 1 — shared contracts and service skeletons

Introduce stable owner-service interfaces and injected runtime dependencies without
moving large algorithms yet.

Add construction/unit tests.

Preserve behavior.

### Milestone 2 — output/restart extraction

Extract output, snapshot verification, restart verification, cadence/event state,
and provenance orchestration.

This is a lower-risk seam with strong existing tests.

Run focused HDF5/restart tests.

### Milestone 3 — IC startup and migration/balance extraction

Extract IC startup invocation and the ownership/decomposition/migration orchestration
without attempting distributed IC reading.

Run IC and ownership/migration tests.

### Milestone 4 — gravity extraction

Move TreePM/PM cadence orchestration and force-cache lifecycle behind `GravityRuntime`.

Run gravity and workflow tests.

### Milestone 5 — hydro/AMR and source extraction

Move hydro, AMR, and source orchestration into their owners while preserving rung-zero
semantics.

Run hydro, AMR, source, restart, and workflow tests.

### Milestone 6 — time coordination and typed views

Reduce the workflow loop to `TimeCoordinator`, complete typed resource-view migration,
remove compatibility adapters, and prove stale/invalidation behavior.

### Milestone 7 — descriptor-driven composition

Make typed descriptors/factories instantiate and register real contributions. Add the
probe-module test.

### Milestone 8 — full validation and adversarial review

Run the broadest available matrix, inspect the entire diff, and repair:

- replacement god-files;
- cross-owner private dependencies;
- lingering unrestricted state access;
- stale views;
- duplicated runtime authority;
- hidden service construction;
- changed numerical order;
- restart/provenance drift;
- module-specific hardcoding in the composition root.

Do not complete the Goal while any mandatory gate is unmet.

---

## 8. Validation commands

Use repository presets and instructions from `AGENTS.md`. At minimum, attempt:

```bash
bash scripts/ci/check_repo_hygiene.sh

cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure

cmake --preset asan-debug
cmake --build build/asan-debug
ctest --test-dir build/asan-debug --output-on-failure

cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

If MPI dependencies are available, also configure/build the relevant MPI path and run
focused sequential tests. Do not make MPI runtime claims when MPI is unavailable.

After each milestone, run the smallest affected set before proceeding.

Never resolve failures by weakening tolerances, suppressing tests, increasing timeouts
without evidence, or disabling assertions.

---

## 9. Goal status and progress updates

During the run, compact progress reports should state:

- current milestone;
- files/responsibilities moved;
- focused tests run and result;
- current `reference_workflow.cpp` line count;
- remaining mandatory gates;
- blockers, if any.

Do not mark the Goal complete merely because a milestone or budget segment ended.

A budget stop is not completion. The correct state is blocked/budget-limited with the
remaining gates listed.

---

## 10. Final response contract

The final response must contain:

1. completion verdict: complete, blocked, or budget-limited;
2. before/after architecture;
3. before/after workflow and largest-file line counts;
4. owner-service responsibility map;
5. typed resource-view map;
6. module descriptor/factory behavior and probe-module proof;
7. exact files changed;
8. exact tests and outcomes;
9. gate-by-gate result for GATE-A through GATE-G;
10. preserved pre-existing changes;
11. remaining limitations and risks;
12. `git status --short --branch`, `git diff --stat`, and `git diff --check`.

It must not claim that a PR was opened or pushed.

---

## 11. Independent review after the Goal

After Codex reports completion, do not immediately accept it.

In a new Codex thread or with a separate high-reasoning reviewer, run `/review` on
uncommitted changes using this instruction:

```text
Adversarially review the uncommitted CHUÍ runtime-decomposition patch against
docs/repair/runtime_decomposition_goal_spec.md and
docs/repair/chui_runtime_brain_prompt_result_adversarial_audit.md.

Do not assess only compilation or naming. Verify every GATE-A through GATE-G claim
against the actual production path. Look specifically for a replacement god-file,
callbacks that still receive unrestricted mutable SimulationState, descriptors that
remain decorative metadata, changed rung-zero numerical order, stale views after
migration/compaction, hidden MPI-context construction, dependency cycles, restart or
provenance drift, tests that exercise scaffolding rather than production, and line
count gaming.

Return a gate-by-gate pass/fail verdict with exact source anchors and commands. Do not
modify the repository.
```

Accept the experiment only if this independent review confirms all mandatory gates.

---

## 12. Follow-up campaigns

Only after Campaign A passes should separate Goals address:

1. canonical multifile distributed IC ingestion;
2. production hierarchical cosmological KDK with per-element epochs, prediction,
   neighbor limiter/wake-up, and restart;
3. elastic restart;
4. task-based communication overlap and asynchronous output.

Do not merge these objectives back into Campaign A.
