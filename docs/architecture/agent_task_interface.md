# Agent task interface

This document turns the root `AGENTS.md` contract into an executable workflow for Codex, Copilot-style coding agents, Cline/Continue harnesses, and human reviewers. It is intentionally explicit because this repository is a scientific codebase: an agent that only produces plausible-looking C++ can damage numerical validity, restart safety, or future MPI scalability.

## 1) Repository identity

- Project name in current user-facing discussion: **CHUI**.
- Historical/code identity: `CosmoSim`, `cosmosim`, and `include/cosmosim/...` remain compatibility names until a dedicated migration plan exists.
- Primary target: desktop-first and small-cluster-first cosmological simulation with a scale-up path.
- Scientific center of gravity: credible cosmological cubes and zoom-in galaxy formation, with dark-matter-only maturity preceding baryonic complexity.

## 2) Agent entrypoints

Agents should read these files in order for any substantial change:

1. `AGENTS.md`.
2. `README.md`.
3. `CONTRIBUTING.md`.
4. `docs/architecture/overview.md`.
5. `docs/architecture/developer_workflow_contract.md`.
6. `docs/architecture/runtime_truth_map.md`.
7. `docs/repair_state_recap.md`.
8. `docs/repair_open_issues.md`.
9. The module document matching the touched subsystem.

The agent must stop and re-scope if the requested edit conflicts with any of these contracts.

## 3) Task-intent router

| User wording | Agent mode | Default behavior |
|---|---|---|
| "audit", "review", "readiness", "brutal analysis" | Audit | Inspect and report; do not patch unless asked. |
| "fix", "patch", "harden", "repair" | Repair | Make minimal production changes tied to the defect. |
| "implement Stage", "add feature", "one-shot implement" | Feature implementation | Implement code, docs, tests, and evidence floor. |
| "return a clean zip" | Packaging | Produce a source-only archive with build artifacts excluded. |
| "update agent interface" | Documentation/interface repair | Update `AGENTS.md`, agent docs, and any linked contributor surfaces. |

When the request mixes modes, perform the most concrete requested deliverable first. For example, "audit and patch" means inspect, patch confirmed blockers, and report what remains.

## 4) Patch-planning checklist

Before editing, answer internally:

- Which module owns this behavior?
- Is this public API, internal implementation, config, schema, docs, tests, or validation data?
- Is there a single runtime owner, or am I creating a shadow authority?
- Does restart/snapshot/provenance state change?
- Does MPI ownership, ghost validity, sidecar alignment, or scheduler truth change?
- Does this need a benchmark or profiling hook?
- Which smallest command set can prove the invariant?

If those questions cannot be answered from the current files, inspect more of the repo rather than inventing architecture.

## 5) Required implementation shape

A real patch should contain the smallest coherent subset of:

- Header/API change with ownership and caller-visible behavior documented.
- Implementation change without placeholder logic.
- Tests covering the local invariant and at least one integration path when state/schema/parallel behavior changes.
- Docs update for any user-visible, architecture, config, schema, validation, or workflow change.
- Benchmark/profiling hook for hot paths or performance-sensitive data layout changes.

Documentation-only patches still need link consistency and hygiene evidence.

## 6) Scientific red lines

Agents must not:

- Add a hydro, AMR, star-formation, feedback, black-hole, chemistry, or cooling feature as a mere flag without a model contract.
- Claim full galaxy-formation readiness from dark-matter-only evidence.
- Treat single-rank tests as proof of distributed correctness.
- Treat a successful compile as proof of restart/schema/config correctness.
- Replace conservative hydro/AMR accounting with visually plausible but non-conservative logic.
- Ignore unit/frame suffixes in public APIs or persistent data.

## 7) Distributed-memory red lines

Agents touching MPI, decomposition, ghosts, migration, PM slabs/pencils, hydro ghosts, or AMR patch exchange must define:

- Owner rank and mutation rights.
- Ghost/import epoch or version.
- Serialization payload completeness.
- Sidecar migration and rebuild behavior.
- Stale-ghost invalidation behavior.
- Restart continuation behavior.
- Load model impact beyond raw particle count.

If a feature is only single-rank or pseudo-distributed, say so in docs and tests.

## 8) Clean zip procedure

When packaging, create the archive from a clean staging copy and exclude generated artifacts. A suitable exclusion list includes:

```text
.git/
build/
build-*/
cmake-build-*/
CMakeFiles/
CMakeCache.txt
Testing/
output/
__pycache__/
.pytest_cache/
.mypy_cache/
.ruff_cache/
.ipynb_checkpoints/
*.o
*.a
*.so
*.dll
*.dylib
*.exe
*.log
CMakeUserPresets.json
```

Keep `CMakeUserPresets.json.example`, source tests, benchmark sources, validation decks, and documentation.

## 9) Final-report template

Use this shape after a patch:

```text
Changed:
- ...

Evidence:
- command: result

Touched files:
- path: reason

Limitations:
- ...

Artifact:
- clean zip link
```

Do not claim that unrun dependency-enabled paths passed. Report environment blockers explicitly.

## 10) Current practical priority bias

Unless the user says otherwise, bias fixes toward:

1. Config/build/test hygiene and reproducibility.
2. Runtime truth ownership, restart safety, and sidecar exactness.
3. TreePM correctness and first non-baryonic cosmological cube readiness.
4. MPI ownership/decomposition/exchange contracts.
5. AMR and hydro foundations only after their ownership, conservation, and validation interfaces are explicit.
6. Moving Voronoi mesh as an RFC/research branch, not a drive-by core rewrite.
