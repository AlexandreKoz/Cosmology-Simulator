# `src/io` Audit Closure Matrix — 2026-08-09

Source audit lineage: 2026-08-06 `src/io` adversarial audit. The campaign was applied to
the then-current repository rather than blindly using historical line numbers.

Status vocabulary:

- **Fixed** — production code changed and the source-level defect is closed.
- **Superseded** — a broader architectural change removes the original failure mode.
- **Already/current-tree** — the current tree already contained part of the required gate;
  residual drift was corrected during this campaign.
- **Environment-limited validation** — implementation is present but the campaign host
  lacked the runtime/platform needed to execute that branch.

| ID | Status | Closure |
|---|---|---|
| IO-001 | Fixed | Shared parent/transaction handling accepts bare restart filenames. |
| IO-002 | Fixed | Unique same-directory temporary files; empty/path-like suffixes cannot alias the final checkpoint. |
| IO-003 | Superseded | Snapshot HDF5 handles close before transactional publish; no delete-old window; hard-link failures are checked. Windows runtime validation is environment-limited. |
| IO-004 | Superseded | Restart uses the same durable transactional layer with POSIX file+directory sync and a guarded native Windows path. Windows runtime validation is environment-limited. |
| IO-005 | Fixed | Explicit snapshot dialect and caller species mapping; ambiguous populated external PartType2/3 fail closed. |
| IO-006 | Superseded | Central snapshot conversion descriptor handles `a`, `h`, frame and unit semantics in both directions. |
| IO-007 | Fixed | PartType5 science output/readback now preserves current BH science sidecar state using documented CHUÍ extension fields. |
| IO-008 | Fixed | Explicit missing-field policy and availability/evolution-ready reporting replace silent believable zero physics. |
| IO-009 | Fixed | Read budgets, checked arithmetic, bounded strings, shape/type checks and post-import scientific validation precede/guard materialization. |
| IO-010 | Fixed | Restart digests are generated in one logical payload traversal rather than two full state walks. |
| IO-011 | Superseded | v23 adds canonical little-endian SHA-256; legacy FNV is compatibility-only. |
| IO-012 | Superseded | Snapshot datasets are filled through bounded hyperslab staging; scratch memory is policy-bounded instead of species-sized. |
| IO-013 | Fixed | Removed the redundant session-open full-file SHA pass; normal strict lifecycle is inspection hash plus completion verification rather than three complete hashes. |
| IO-014 | Fixed | Successful routing protocol consolidated from 30 to 23 communicator-wide calls per batch without deleting exact reconciliation/audits. MPI runtime validation is environment-limited. |
| IO-015 | Fixed | MPI wrappers request error-return behavior, check operation statuses and abort on communicator failure instead of entering another unsafe consensus. MPI runtime validation is environment-limited. |
| IO-016 | Superseded | Versioned wire v2 uses compact common data plus species-specific payload; DM no longer carries irrelevant sidecar zeros. MPI runtime validation is environment-limited; codec is independently testable. |
| IO-017 | Fixed | Checked narrowing is applied before assignment to bounded integer/count fields. |
| IO-018 | Fixed | Manifest parser now handles Unicode escapes/control characters and enforces JSON file/string/depth/container budgets. |
| IO-019 | Fixed | Canonical metadata reads are bounded; special files and canonical path escapes are rejected under the containment policy. |
| IO-020 | Superseded | Snapshot-set inspection/discovery and aggregation support single and logical multifile inputs with cross-member agreement checks. |
| IO-021 | Fixed | Reader inspects HDF5 dataset creation properties; unknown/mixed storage policy is represented honestly. |
| IO-022 | Fixed | Active snapshot schema/dialect, normalized config/hash, naming-rule versions, member/set identity and provenance are validated/persisted. |
| IO-023 | Fixed by explicit contract | Per-member dense local sidecar indices are explicitly 32-bit and capped; logical multifile totals are 64-bit. The code no longer advertises unsupported >32-bit local materialization. |
| IO-024 | Fixed | Stale `validateRecordScientificState` declaration was aligned with the implementation. |
| IO-025 | Superseded in the high-risk paths | Transaction semantics, snapshot conversion semantics and snapshot-set semantics were centralized. Restart field enumeration remains intentionally explicit for backward-compatible schema review rather than being replaced by risky metaprogramming during the integrity migration. |
| IO-026 | Superseded by responsibility extraction | Transaction, conversion and snapshot-set orchestration were split into dedicated units, substantially reducing the monolithic snapshot/restart coupling. Further cosmetic file splitting is not required for correctness and was deliberately not used as a substitute for production code. |
| IO-027 | Already/current-tree + Fixed residuals | CI already defines strict first-party warning flags. Residual I/O warning findings were fixed; every `src/io` translation unit passed the repository's GCC strict warning set independently. Whole-repo strict build remains blocked by unrelated pre-existing warnings outside I/O. |
| IO-028 | Fixed | Active output/restart docs and release gates were updated to v6/v23 and current dataset/finalization semantics. Historical repair documents remain historical. |
| IO-029 | Fixed | Generated isolated IC helper now uses checked counts/IDs and deterministic positions inside the configured box. |
| IO-030 | Fixed | Exact distributed-ID scratch supports configured run-scoped scratch, collision-resistant job token and RAII cleanup. MPI runtime validation is environment-limited. |
| IO-031 | Superseded | Manifest publication uses the shared transaction primitive; remove-before-rename fallback is gone. |
| IO-032 | Fixed | Internal HDF5 RAII closes general property-list identifiers via `H5Pclose`. |
| IO-033 | Fixed | Canonical-limit diagnostics no longer direct users to a nonexistent general multifile writer. |
| IO-034 | Fixed | IC byte accounting is explicitly named/documented as logical metadata/payload bytes rather than physical filesystem traffic. |

IO-035 was an engineering-process observation, not an audit defect or authorship fact. The
campaign addressed its actionable causes through central contracts, warning cleanliness,
bounded inputs, explicit schemas, and smaller responsibility-specific I/O components.

## Validation note

The final campaign validation record belongs in the delivery handoff. This document does
not convert unavailable MPI/Windows execution into a pass claim. Source implementation
and focused serial/HDF5 evidence are kept separate from environment-limited validation.

## Campaign validation evidence

Executed in the campaign Linux environment (GCC 14.2, HDF5 1.14.5, MPI disabled):

- built `cosmosim_io` and `cosmosim_workflows` successfully;
- 7/7 focused I/O acceptance tests passed:
  `unit_ic_reader`, `unit_ic_record_codec`, `unit_snapshot_hdf5_schema`,
  `unit_restart_checkpoint_schema`, `integration_ic_generated_pipeline`,
  `integration_snapshot_hdf5_roundtrip`, and
  `integration_restart_checkpoint_roundtrip`;
- `integration_docs_scaffold` and `integration_release_readiness_artifacts` passed;
- restart-equivalence cases passed for the harness baseline, DM-only, TreePM,
  multirate bins, output-enabled, AMR hydro, AMR flux registers, and AMR temporal
  ghosts;
- every `src/io` translation unit passed the repository CI GCC warning set
  (`-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Wundef
  -Wformat=2 -Wnull-dereference -Werror`) when compiled independently with its
  configured compile command;
- repository hygiene checks passed.

Two broader unchanged restart-equivalence fixtures currently fail state-validity gates:
`integration_restart_equivalence_hydro_toy` fails its own
`SimulationState::validateOwnershipInvariants()` assertion before checkpoint I/O, and
`integration_restart_equivalence_stochastic_sources` reaches the pre-existing restart
writer state-validity rejection (`cannot checkpoint invalid simulation state`). Neither
fixture file nor `src/core`/`src/physics` was changed by this campaign, and the same
restart state-validity rejection existed in the uploaded original tree. An attempted
clean comparison build of the original repository exceeded the campaign execution
window before those executables were produced, so this note does not claim an executed
original-tree pass/fail comparison. The production restart validator was not weakened to
make invalid fixtures pass.

MPI runtime tests were not executed because the environment lacks MPI development/runtime
support. Windows transactional code was not executed because the environment is Linux.
No MPI/Windows runtime pass is claimed.
