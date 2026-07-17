# CHUÍ Runtime-Brain One-Shot Prompt — Adversarial Implementation Audit

**Audited snapshot:** `Cosmology-Simulator-main(14).zip`  
**Baseline snapshot:** `Cosmology-Simulator-main(13)(1).zip`  
**Prompt assessed:** `chui_runtime_brain_one_shot_codex_prompt.md`  
**Audit date:** 2026-07-16  
**Mode:** code-first adversarial audit; repository source was not modified

---

## 1. Executive verdict

Codex produced a **coherent, testable partial hardening patch**, but it did **not** one-shot the requested runtime-brain transformation.

The patch is not fake work. It introduces several useful changes:

- a truthful machine-readable runtime capability report;
- an explicit IC manifest and field-conversion model;
- stricter IC header/schema validation;
- PartType5 black-hole import;
- synchronization-safe deferred timestep coarsening;
- active-index spans instead of per-step active-list copies;
- a persistent transient workspace;
- exact final-time clipping;
- periodic code-time output events persisted through restart;
- much stronger snapshot roundtrip verification;
- a shared root MPI/runtime-services bundle;
- a narrow collective failure coordinator;
- a candid living execution/closeout document.

The implementation compiles and its tested paths are healthy. A clean CPU-debug build passed all **125/125** configured tests. Six focused HDF5 tests passed, and four focused scheduler/runtime tests passed under ASan/UBSan.

However, the central mission remains unfinished:

- production still rejects every nonzero hierarchical rung;
- the global timestep remains fixed;
- the main workflow was not decomposed and grew to **8,008 lines**;
- every MPI rank still reads and materializes the complete single IC file;
- there is no hydro neighbor timestep limiter or wake-up;
- restart remains same-topology/rank-local;
- stage contracts still do not restrict callback mutation of `SimulationState`;
- module descriptors do not drive composition;
- general rank-local callback exceptions remain capable of diverging collective order;
- asynchronous output and a task/resource execution graph were not implemented.

### Strict prompt-level result

| Classification | Count | Meaning |
|---|---:|---|
| Fully closed | **3 / 19** | The requested production behavior is implemented and meaningfully tested |
| Partially closed | **10 / 19** | Useful infrastructure or a narrower subset exists, but the prompt’s closure condition is unmet |
| Still open | **6 / 19** | The decisive production capability is absent |

The engineering quality of the attempted changes is good. The completion of the requested scope is low. The result is best described as:

> **A respectable Phase-0/Phase-1 hardening pass, not the requested runtime-brain replacement.**

---

## 2. Change profile

Excluding build trees and generated test outputs, the new snapshot contains:

- **34 modified files**
- **9 new files**
- approximately **3,027 inserted lines**
- approximately **257 deleted lines**

Change distribution:

| Area | Files | Insertions | Deletions |
|---|---:|---:|---:|
| `src/` | 10 | 1,656 | 202 |
| `tests/` | 13 | 559 | 22 |
| `docs/` | 9 | 519 | 28 |
| `include/` | 9 | 277 | 4 |
| build/release metadata | 2 | 16 | 1 |

Important new files include:

- `docs/repair/chui_runtime_brain_exec_plan.md`
- `include/cosmosim/workflows/runtime_capabilities.hpp`
- `src/workflows/runtime_capabilities.cpp`
- `include/cosmosim/workflows/runtime_services.hpp`
- `src/workflows/runtime_services.cpp`
- `src/io/ic_manifest.cpp`
- `tests/unit/test_runtime_capabilities.cpp`
- `tests/unit/test_runtime_services.cpp`
- `tests/unit/test_module_registry.cpp`

The dominant source file remains:

```text
src/workflows/reference_workflow.cpp: 8,008 lines
```

The audited baseline was approximately 7,589 lines. The supposed composition root therefore grew by roughly **419 net lines**, with about 545 additions and 126 removals.

---

## 3. Validation performed in this audit

### 3.1 Repository hygiene

Executed:

```bash
bash scripts/ci/check_repo_hygiene.sh
```

Result: passed.

The executable bit was not preserved by the uploaded ZIP, so invoking the script directly produced a permission error. Running it through `bash` passed. This is packaging behavior, not a source failure.

### 3.2 Complete CPU-debug matrix

Executed:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -j4
```

Result:

```text
100% tests passed, 0 tests failed out of 125
```

This includes the real rung-zero workflow, scheduler infrastructure, current restart harnesses, hydro/AMR tests, TreePM validations, and newly added runtime metadata tests.

It does **not** prove production mixed-rung stepping, distributed IC import, elastic restart, wake-up, or multi-rank failure coordination because those paths do not exist.

### 3.3 Focused HDF5 matrix

Built and executed:

```text
integration_reference_workflow
unit_runtime_capabilities
unit_ic_reader
unit_restart_checkpoint_schema
integration_snapshot_hdf5_roundtrip
integration_restart_checkpoint_roundtrip
```

Result:

```text
6/6 passed
```

A complete HDF5 tree build was started but did not finish within the audit execution window. The six directly affected HDF5 paths above were built and run successfully.

### 3.4 Sanitizers

Built and ran under ASan/UBSan:

```text
unit_time_integration
unit_runtime_capabilities
unit_runtime_services
unit_module_registry
```

Result:

```text
4/4 passed
```

### 3.5 MPI limitation

Attempted:

```bash
cmake --preset mpi-hdf5-fftw-debug
```

Configuration failed because the environment has no working MPI C++ installation:

```text
Could NOT find MPI_CXX
Could NOT find MPI
```

No multi-rank runtime test was executed. This audit therefore makes no MPI-pass claim.

---

## 4. What Codex genuinely accomplished

## 4.1 Runtime capability truth — good closure of the safety/reporting milestone

`runtime_capabilities.cpp` introduces a serialized capability report distinguishing:

- fixed global timestep: supported;
- adaptive global timestep: unsupported;
- production hierarchical local timestep: unsupported;
- canonical external IC import: provisional;
- distributed IC import: unsupported;
- rank-remappable restart: unsupported;
- asynchronous output: unsupported.

This is valuable because the repository now describes its real production envelope in machine-readable form rather than relying only on prose.

The associated test verifies the unsupported nonzero-rung request is rejected.

### Limitation

`validateRequestedRuntimeCapabilities()` currently gates only hierarchical timesteps. The report does not yet drive a generalized capability/configuration resolver for all unsupported combinations.

There is also a stale detail:

```text
runtime_capabilities.cpp:81 says restart schema v20
```

The repository now writes restart schema **v21**.

**Verdict:** useful and honest safety infrastructure; not feature closure.

---

## 4.2 Deferred timestep coarsening — genuinely implemented

The prior scheduler threw whenever a larger requested bin was not aligned at the current tick. The new implementation:

- keeps the requested coarsening in `pending_bin_index`;
- continues the element on its current smaller bin;
- commits the larger bin at the next shared synchronization boundary;
- persists the pending transition through restart;
- records deferred waits.

Relevant implementation:

```text
src/core/time_integration.cpp:1489-1498
src/core/time_integration.cpp:1581-1625
```

The new unit tests cover:

- deferral;
- eventual synchronized commit;
- preservation across scheduler export/import.

Because an element active on a coarser old bin is automatically aligned to any finer bin, the immediate-refinement path remains consistent.

**Verdict:** `BRAIN-009` is closed.

---

## 4.3 Endpoint clipping and code-time output events — real improvement

The workflow now:

- computes the remaining interval to `t_code_end`;
- clips the final step rather than overshooting;
- clips a step to the next periodic code-time snapshot event;
- restores the pre-clipped base step afterward;
- persists code-time output cadence and next event in restart schema v21;
- resumes the original base step after a restart written at a clipped event;
- reports endpoint/output-event clipping in runtime diagnostics.

Relevant implementation:

```text
src/workflows/reference_workflow.cpp:7029-7088
src/workflows/reference_workflow.cpp:7713-7755
src/workflows/reference_workflow.cpp:7793-7827
```

The real workflow integration test covers both endpoint clipping and a restart at a code-time event.

**Verdict:**

- `BRAIN-019` is closed.
- `BRAIN-017` is only partially closed because the prompt requested event scheduling in scale factor, redshift, physical/code time, step count, and wall-clock safety. Only periodic code-time plus the existing step cadence were implemented.

---

## 4.4 Stronger snapshot verification — substantial but narrower than claimed

The old workflow treated equal particle count as a successful snapshot roundtrip. The new verification checks:

- snapshot schema and file kind;
- config/provenance identity;
- header time, redshift, box dimensions, and cosmology;
- particle IDs and species;
- coordinates, velocities, masses;
- softening values and override masks;
- tracer sidecars;
- stable-ID keyed equivalence.

This is a meaningful correction and is properly tested.

### Critical limitation

The science snapshot writer itself still writes primarily:

- `Coordinates`
- `Velocities`
- `Masses`
- `ParticleIDs`
- gravity softening fields
- tracer-specific fields

It does **not** serialize the principal gas, star, or black-hole science state such as:

- gas density;
- internal energy;
- pressure/temperature;
- star formation/birth metadata;
- stellar metallicity;
- black-hole subgrid mass and accretion fields.

Consequently, the new verifier is strong for the **current narrow snapshot schema**, but it cannot prove a full hydro/galaxy-formation snapshot roundtrip.

The header also still sets `Flag_StellarAge` and `Flag_Metals` to one even though corresponding stellar-age/metal datasets are not written. That pre-existing inconsistency remains.

**Verdict:** `BRAIN-013` is partially closed, not fully closed under the one-shot prompt’s requirements.

---

## 4.5 IC manifest and header validation — strong infrastructure

The new IC layer adds an explicit model for:

- dialect and version;
- source files and provenance identifiers;
- per-file and total counts;
- high-word counts;
- `NumFilesPerSnapshot`;
- mass table;
- box size, scale factor, redshift, cosmology, and Hubble parameter;
- dataset path/type/rank/dimensions;
- base SI unit;
- `h` exponent;
- scale-factor exponent;
- comoving/physical frame;
- velocity convention;
- extensive/intensive/specific semantics;
- conversion/default/drop/reject disposition;
- particle-type policy.

The HDF5 reader now requires and cross-checks much more header metadata. It validates dataset ranks/classes and positive particle masses. It supports:

- physical/comoving coordinate conversion;
- explicit `h` and `a` scaling;
- physical peculiar, square-root-a-scaled, and comoving-coordinate-rate velocities;
- gas density and internal-energy unit conversion;
- PartType5 black-hole species and sidecar mapping;
- fail-closed PartType2/3 handling unless an explicit policy is supplied;
- local duplicate-ID rejection.

The unit tests exercise manifest conversion math, high-word metadata, PartType5, PartType2 rejection/override, malformed input, and duplicate IDs.

This is the largest useful body of work in the patch.

### Why it is still not a safe general GADGET/AREPO importer

The default production path does not receive a user-supplied manifest. It calls:

```cpp
makeGadgetArepoBridgeV1Manifest(ic_path, schema)
```

That default contract hardcodes:

- kpc;
- solar masses;
- km/s;
- comoving coordinates;
- physical peculiar velocities;
- zero `h` exponent;
- zero scale-factor exponent.

Relevant code:

```text
src/io/ic_reader.cpp:128-215
src/io/ic_reader.cpp:552-557
src/io/ic_reader.cpp:667-676
```

These conventions are not inferred from the HDF5 file and are not universal across GADGET/AREPO IC producers. Calling them a named versioned bridge makes the assumption explicit, which is better than a hidden assumption, but the API still has a generic-looking name:

```cpp
readGadgetArepoHdf5Ic(...)
```

A user can therefore still feed a conventional h-scaled or sqrt(a)-velocity file to the default path and obtain a physically incorrect result.

Additional limitations:

- there is no JSON manifest parser;
- no typed config key selects an external manifest;
- explicit manifests are accessible only to C++ callers;
- no canonical CHUÍ converter is integrated;
- the synthesized manifest reports assumed field types such as `float64` rather than recording the actual HDF5 datatype;
- complete file hashing is performed in the read path and would be repeated by every MPI rank;
- optional gas thermodynamic fields may still default to zero;
- star formation time/metallicity and richer BH fields remain incomplete.

**Verdict:** `BRAIN-003` and `BRAIN-004` are partially closed.

---

## 4.6 Active spans and workspace reuse — real micro-level improvement

The workflow no longer copies scheduler active spans into fresh vectors each step. It now passes read-only spans and retains a `TransientStepWorkspace` outside the loop, clearing and reusing it.

It also evaluates timestep criteria only for active indices after a completed step.

These are good low-level corrections, and profiling counters were added.

### Hidden remaining O(N) work

Before evaluating those active indices, each call still constructs a fresh `AdaptiveTimeStepCriteriaStorage` and calls `buildAdaptiveTimeStepCriteriaView()`, which:

- rebuilds the particle-ID-to-row map;
- copies all gas-cell IDs;
- traverses all gas cells;
- rebuilds all CFL/patch metadata;
- allocates/resizes several full-cell vectors.

Relevant code:

```text
src/workflows/reference_workflow.cpp:494-637
src/workflows/reference_workflow.cpp:852-894
```

Thus the criteria loops are active-only, but their view/materialization path remains full-state and allocation-heavy. In a true fine-rung hierarchy this could still erase much of the expected benefit.

**Verdict:** `BRAIN-008` and `BRAIN-018` are partially closed.

---

## 4.7 Root runtime-services injection — narrow success

The composition root now constructs one `RuntimeServices` bundle and passes the root `MpiContext` into drift, gravity, TreePM construction, and hydro. This removes several workflow-local context constructions.

`BRAIN-016` can reasonably be considered closed for the audited `ReferenceWorkflow` production path.

### Scope caveat

`RuntimeServices` currently contains only:

- MPI context;
- profiler;
- deterministic flag.

It does not yet own the allocator/workspace, topology, output, logging, backend, task executor, or other services requested by the prompt. A fallback `TreePmCoordinator` constructor elsewhere still constructs an `MpiContext`, so the “only the composition root constructs MPI context” claim is not universally true across the repository.

---

## 5. What the prompt did not accomplish

## 5.1 Production hierarchical KDK remains completely disabled

The decisive line still exists:

```text
src/workflows/reference_workflow.cpp:7422-7424
```

```cpp
if (config.numerics.hierarchical_max_rung != 0) {
    throw ...
}
```

Configuration validation also rejects any nonzero rung.

No per-element authoritative state was added for:

- last drift epoch;
- last kick epoch;
- independent begin/end epochs;
- predicted inactive state;
- force/source validity epochs;
- wake-up reason;
- correct elapsed-interval cosmological factors.

The existing hierarchical tests still test scheduler infrastructure, not a nonzero-rung `ReferenceWorkflowRunner`.

**Verdict:** `BRAIN-001` remains open.

---

## 5.2 Production timestep remains fixed

The base global interval is still initialized once from the override or configured segment length. Timestep criteria map candidates into scheduler bins whose minimum interval is the same base `dt_time_code`.

Since production is restricted to rung zero, all candidates collapse to the same fixed global interval. Hydro can reject an unsafe step, but the runtime does not reduce the next global base step.

The new active-only evaluation work is supporting infrastructure, not adaptive global stepping.

**Verdict:** `BRAIN-002` remains open.

---

## 5.3 Multifile/distributed IC ingestion was not implemented

Although `IcManifest` can describe multiple files, the actual reader explicitly rejects them:

```text
src/io/ic_reader.cpp:685-689
```

The workflow still loads the IC before distributed ownership decomposition. Every rank enters the same synchronous single-file read path, hashes the same file, allocates the full global state, and then decomposes it.

There is no:

- reader-group assignment;
- file/chunk distribution;
- direct-to-owner routing;
- per-rank memory proof;
- cross-rank duplicate-ID test;
- np1/np2/np4 distributed IC equivalence test.

The capability report truthfully labels this unsupported.

**Verdict:** `BRAIN-005` remains open.

---

## 5.4 The maintainability refactor failed

The one-shot prompt explicitly required `reference_workflow.cpp` to become a small composition/lifecycle root and rejected merely moving the monolith.

Instead:

- the file grew from roughly 7,589 to **8,008 lines**;
- gravity, hydro, AMR, IC startup, migration, timestep criteria, output, restart verification, scheduler wiring, and the main loop remain inside it;
- no `InitialConditionService`, `TimeCoordinator`, `GravityRuntime`, `HydroAmrRuntime`, `MigrationBalanceService`, or `OutputRestartService` was extracted;
- `runtime_services.cpp` is only 34 lines;
- `runtime_capabilities.cpp` is reporting logic, not runtime decomposition.

This is the clearest failure against the maintainability objective.

**Verdict:** `BRAIN-006` remains open and the concentration metric worsened.

---

## 5.5 Stage contracts remain non-enforceable

Production callbacks still receive a mutable `StepContext` containing full mutable `SimulationState` and `IntegratorState`. Nothing in the type system prevents a hydro callback from mutating unrelated particle, scheduler-mirror, BH, or cache state.

The new `AdaptiveTimeStepCriteriaView` is a useful narrow read view for one helper path. It is not a replacement for production stage/task resource views.

No runtime hazard table or typed stage-specific read/write capability system was implemented.

**Verdict:** `BRAIN-007` is only partially addressed.

---

## 5.6 Hydro limiter/wake-up remains absent

No production mechanism was added to:

- bound neighboring hydro timestep ratios;
- wake inactive local neighbors;
- propagate wake-ups across MPI neighbors;
- reactivate cells affected by shocks or feedback;
- prove no stale inactive flux/source state is consumed.

The repository remains rung zero, so this omission is not currently activated as a mixed-rung correctness bug, but it blocks enabling local hydro timesteps.

**Verdict:** `BRAIN-010` remains open.

---

## 5.7 Spawned-element registration remains non-atomic

The patch ensures scheduler coverage is refreshed before criteria are indexed after source/AMR stages. This is useful defensive cleanup.

It does not create the requested single ownership transaction coupling:

- stable identity;
- sidecars/species;
- scheduler epochs;
- active-list membership;
- mirrors;
- force/source caches;
- migration;
- restart payload.

**Verdict:** `BRAIN-011` is partially closed.

---

## 5.8 Elastic restart remains absent

Restart schema v21 adds code-time output cadence state. It does not add:

- logical distributed checkpoint manifests;
- rank-count remapping;
- stable-ID redistribution into a new topology;
- PM topology reconstruction across changed rank counts;
- 1→2, 2→1, or 2→4 continuation.

The same-world-size/same-rank topology restriction remains.

**Verdict:** `BRAIN-012` remains open.

---

## 5.9 General MPI failure coordination remains unsafe

`FailureCoordinator` performs an allreduce over a local-failure flag and preserves the local exception cause. That is a useful primitive.

In the production workflow it is used narrowly around restart-topology validation. General stage callbacks are still surrounded only by local C++ exception propagation and the top-level local catch.

A rank throwing before peers enter a later collective can therefore still strand peers. No failpoint-driven multi-rank tests were added.

The unit test runs only a synthetic single-rank context and cannot establish collective safety.

**Verdict:** `BRAIN-014` is partially closed.

---

## 5.10 Module descriptors do not drive composition

The new `ModuleDescriptor` table is a compile-time array of semicolon-delimited strings such as:

```text
"gravity_kick_pre;force_refresh;gravity_kick_post"
```

The unit test checks that every string is non-empty and retrievable.

The table does not:

- instantiate modules;
- contribute typed tasks;
- register timestep criteria;
- register migration fields;
- register restart serializers;
- validate resource hazards;
- eliminate workflow edits when adding a module.

The workflow still manually constructs and registers every callback.

This is useful architecture metadata, but it is not a module-composition system.

**Verdict:** `BRAIN-015` is partially closed.

---

## 5.11 No task/resource graph, communication overlap, or asynchronous output

The fixed callback sequence remains. `StageScheduler` was not transformed into a dynamic active-resource task planner. Communication was not represented as dependency tasks, and no output staging pipeline was added.

The capability report correctly says asynchronous output is unsupported.

These omissions correspond to most of the prompt’s Milestone 5.

---

## 6. Strict BRAIN-001 through BRAIN-019 closure matrix

| ID | Strict status | Audit evidence |
|---|---|---|
| **BRAIN-001** | **Open** | Real workflow and config still reject `hierarchical_max_rung != 0`; no per-element KDK epochs |
| **BRAIN-002** | **Open** | Base global `dt_time_code` remains fixed; criteria do not resize it |
| **BRAIN-003** | **Partial** | Manifest/conversion model exists, but default bridge still hardcodes a non-universal convention and no config manifest parser exists |
| **BRAIN-004** | **Partial** | PartType5 fixed and Type2/3 fail closed; star/BH/auxiliary mapping remains incomplete |
| **BRAIN-005** | **Open** | Reader explicitly rejects multifile manifests; every rank still reads full single-file IC |
| **BRAIN-006** | **Open** | `reference_workflow.cpp` grew to 8,008 lines; no owner-service decomposition |
| **BRAIN-007** | **Partial** | Narrow criteria view added, but production callbacks retain full mutable state access |
| **BRAIN-008** | **Partial** | Active indices drive candidate loops, but full-state criteria metadata is rebuilt each step |
| **BRAIN-009** | **Closed** | Coarsening is deferred to a legal synchronization point and persists through restart |
| **BRAIN-010** | **Open** | No hydro neighbor timestep limiter or wake-up path |
| **BRAIN-011** | **Partial** | Scheduler coverage repair exists after spawning, but no atomic scheduler/state/migration/restart transaction |
| **BRAIN-012** | **Open** | Restart remains same-topology/rank-local |
| **BRAIN-013** | **Partial** | Strong verification for current schema; gas/star/BH science fields are absent from snapshot schema |
| **BRAIN-014** | **Partial** | Failure coordinator primitive exists but covers only a narrow phase; no general callback/failpoint MPI proof |
| **BRAIN-015** | **Partial** | Descriptive string table exists; it does not drive runtime composition |
| **BRAIN-016** | **Closed, workflow scope** | Root MPI context is injected through the audited production workflow |
| **BRAIN-017** | **Partial** | Periodic code-time events added; redshift/scale-factor/physical-time/wall-clock event system absent |
| **BRAIN-018** | **Partial** | Active-list copies removed and workspace reused; criteria-view rebuilding remains O(N) and allocation-heavy |
| **BRAIN-019** | **Closed** | Endpoint and ordered code-time event clipping implemented and tested |

### Aggregate

```text
Closed:   3
Partial: 10
Open:     6
```

---

## 7. New or residual issues found during this audit

### NEW-001 — capability report has stale restart schema version

`runtime_capabilities.cpp` reports schema v20 while the new restart schema is v21.

**Severity:** low, but undermines the purpose of a machine-readable truth report.

### NEW-002 — default “GADGET/AREPO” bridge remains scientifically hazardous

The new manifest makes the assumptions explicit but still applies them by default to a generic-looking importer. Many legitimate external ICs use different h and velocity scalings.

**Severity:** high for unsupervised external IC use.

### NEW-003 — synthesized manifest can misreport physical file datatype

The recognized manifest populates fields with `scalar_type = "float64"` from expected conventions. The HDF5 reader only verifies a floating-point class and can read float32. The audit manifest may therefore describe the converted target expectation rather than the actual source datatype.

**Severity:** medium provenance/schema-integrity issue.

### NEW-004 — per-rank full-file hashing compounds distributed startup cost

The source file is read completely to compute an FNV-1a provenance identifier inside the same path every rank enters.

**Severity:** medium now; high at realistic multifile/large-IC scale.

### NEW-005 — snapshot flags overstate written science fields

`Flag_StellarAge` and `Flag_Metals` are written as enabled even though corresponding star-age/metal datasets are not part of the current snapshot payload.

**Severity:** medium interoperability/provenance issue.

### NEW-006 — “closed” module registry test verifies strings, not behavior

`test_module_registry` proves that non-empty strings exist. It does not prove a module can be added without editing the workflow or that descriptors register tasks/restart/migration behavior.

**Severity:** medium risk of overstating architectural completion.

### NEW-007 — general collective failure risk remains

The new failure primitive is technically reasonable, but its limited deployment can create false confidence because most callback/collective boundaries remain outside it.

**Severity:** high-confidence MPI architectural risk; not runtime reproduced here.

### NEW-008 — criteria view is still rebuilt globally

Active-loop counters can show only active criteria evaluations while omitting the full-state cost of materializing the view. Current profiling can therefore make local-timestep readiness look better than it is.

**Severity:** medium performance-accounting issue.

---

## 8. Prompt milestone assessment

| Prompt milestone | Result |
|---|---|
| Milestone 0 — capability truth/safety | **Mostly achieved** |
| Milestone 1 — runtime decomposition/resource views | **Failed, with small supporting pieces** |
| Milestone 2 — canonical distributed IC v2 | **Partially achieved at metadata/conversion layer; distributed/canonical tooling absent** |
| Milestone 3 — production hierarchical KDK | **Not achieved** |
| Milestone 4 — restart/output/failure continuity | **Partially achieved** |
| Milestone 5 — task graph/performance architecture | **Small micro-optimizations only** |
| Required closure matrix/reporting | **Achieved candidly** |
| CPU validation discipline | **Strong** |
| HDF5 affected-path validation | **Strong** |
| MPI validation | **Environment-blocked and feature-incomplete** |

---

## 9. Maintainability verdict after the patch

### Before

CHUÍ had modular numerical organs connected through a 7,589-line workflow brain.

### After

CHUÍ has:

- a better capability report;
- a better IC-conversion abstraction;
- a small runtime-services bundle;
- a small failure primitive;
- descriptive module metadata;
- several improved local algorithms;

but the same central workflow now contains **8,008 lines**.

Therefore the maintainability result is:

> **Local maintainability improved; systemic maintainability did not.**

The patch reduces ambiguity in a few subsystems but increases the size of the central policy concentration. It does not change the fact that a future change to timestepping, IC startup, hydro, output, migration, or restart requires understanding a massive shared compilation unit.

Codex’s own closeout correctly admits this. That candor is a strength, but it does not count as implementation.

---

## 10. Final judgment of the Codex experiment

### What this experiment proves

At highest reasoning, Codex can:

- digest a very large architecture prompt;
- identify a safe subset;
- implement several cross-cutting changes coherently;
- add tests and documentation;
- preserve a large existing CPU suite;
- refrain from falsely enabling unsafe hierarchical stepping;
- write a candid closure matrix instead of claiming complete success.

That is valuable. In fact, refusing to “enable” nonzero rungs without per-element epochs was the scientifically correct choice.

### What it does not prove

It does not show that a single giant prompt can reliably deliver three major research-code projects at once:

1. a runtime-service/task-graph refactor;
2. a distributed, canonical cosmological IC system;
3. a physically correct hierarchical cosmological KDK with hydro wake-up and elastic restart.

The prompt exceeded a practical one-shot change horizon. Codex spent its budget on safe scaffolding, validation, and lower-risk closures, then left the hardest coupled features open.

### Experiment score

| Dimension | Grade |
|---|---:|
| Honesty of closeout | **A** |
| Code quality of implemented subset | **B+** |
| Regression discipline | **A-** |
| Scientific safety | **B+** |
| Maintainability improvement | **D+** |
| Completion of requested production capabilities | **D** |
| Overall one-shot outcome | **C-** |

The outcome is not wasted effort. It is a useful preparatory patch. It simply should not be mistaken for completion of the runtime-brain campaign.

---

## 11. Recommended next execution strategy

Do not issue another prompt containing all remaining findings. Split the work into three ordered campaigns.

### Campaign A — runtime decomposition only

Target:

```text
refactor-runtime-workflow-into-owner-services-and-typed-resource-views
```

Acceptance must be structural:

- `reference_workflow.cpp` becomes a small lifecycle/composition root;
- gravity, hydro/AMR, IC, migration, output/restart, and time coordination move into distinct owners;
- no replacement god-file;
- typed resource views restrict production mutation;
- module descriptors register real contributions;
- rung-zero behavior remains equivalent.

This should happen first because adding hierarchical KDK to the current 8,008-line file would fossilize the monolith.

### Campaign B — canonical and distributed IC ingestion

Target:

```text
feat-canonical-distributed-ic-manifest-conversion-and-owner-routing
```

Acceptance:

- manifest parser and typed config integration;
- canonical converter;
- actual source datatype discovery;
- multifile file-set discovery;
- reader-rank assignment;
- direct/bounded routing to final owners;
- no full global state on every rank;
- np1/np2/np4 equivalence and global duplicate-ID tests.

### Campaign C — production hierarchical KDK

Target:

```text
feat-production-hierarchical-kdk-epochs-prediction-wakeup-and-restart
```

Acceptance:

- per-element epochs;
- correct elapsed-interval cosmological factors;
- inactive prediction;
- immediate refinement/deferred coarsening;
- hydro neighbor limiter/wake-up;
- spawned-element atomic registration;
- active-fraction work proof;
- mixed-rung direct/reference/restart tests;
- MPI failure coordination around every new collective phase.

Elastic restart and asynchronous output can then be separate follow-ups rather than hidden inside the KDK campaign.

---

## 12. Bottom line

Codex successfully implemented **the safe perimeter around the brain**, but not the brain itself.

The patch should be retained because it improves:

- truth reporting;
- IC metadata and conversion plumbing;
- scheduler transition semantics;
- endpoint/output-event correctness;
- snapshot verification;
- runtime service injection;
- active-list/workspace behavior;
- tests and documentation.

But the repository is still:

- rung-zero;
- fixed-global-step;
- single-file/full-replication IC;
- same-topology restart;
- callback-sequential;
- non-wake-up hydro;
- centrally orchestrated by an 8,008-line workflow file.

The correct acceptance decision is:

> **Accept as a partial hardening patch, reject as completion of the one-shot prompt.**
