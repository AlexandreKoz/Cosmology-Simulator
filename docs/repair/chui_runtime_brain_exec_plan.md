# CHUI runtime-brain repair execution plan and closeout

Date started: 2026-07-16  
Mode: repair + feature implementation  
Status: in progress

This document is the command-backed execution record for the runtime-brain
repair. It is updated as implementation and validation progress; a checked box
means the named production path and its focused test evidence exist, not merely
that an interface was added.

## Baseline and environment

- Workspace: `/home/xandy/dev/chui`, branch `main...origin/main`.
- The worktree was already heavily modified before this task: 108 tracked files
  changed, 15,403 insertions, 1,931 deletions, plus new MPI and validation test
  sources. All pre-existing changes are treated as user work and preserved.
- Toolchain: CMake 3.28.3, Ninja 1.11.1, GCC/G++ 13.3.0.
- Dependencies discovered: HDF5 1.10.10, FFTW 3.3.10, OpenMPI/MPI 3.1.
- OpenMPI wrapper discovery emits a sandbox socket warning and CMake initially
  omitted `/usr/lib/x86_64-linux-gnu/openmpi/include`. The feature preset
  configures when supplied
  `-DMPI_CXX_HEADER_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include`.
- MPI launch requires execution outside the desktop filesystem/network sandbox
  because PMIx must create local listener sockets.

### Commands executed before edits

| Command | Outcome |
| --- | --- |
| `pwd` | PASS; `/home/xandy/dev/chui` |
| `git status --short --branch` | PASS; recorded the pre-existing dirty tree |
| `git diff --stat` | PASS; 108 tracked files, 15,403 insertions, 1,931 deletions |
| `find .. \\( -name AGENTS.md -o -name AGENTS.override.md \\) -print` | PASS; only repository-root `AGENTS.md` found |
| `cmake --version` | PASS; 3.28.3 |
| `ninja --version` | PASS; 1.11.1 |
| `c++ --version` | PASS; GCC 13.3.0 |
| `mpicxx --version` | Toolchain present; OpenMPI emitted sandbox `socket() failed with errno=1` |
| `cmake --preset cpu-only-debug` | PASS |
| `cmake --build --preset build-cpu-debug -j4` | PASS; 388 build steps |
| focused CPU scheduler/workflow/IC CTest selection | PASS; 8/8 |
| `cmake --preset hdf5-debug` | PASS; HDF5 1.10.10 |
| `cmake --preset pm-hdf5-fftw-debug` | PASS; FFTW 3.3.10 |
| `cmake --preset mpi-hdf5-fftw-debug` | FAIL at MPI C++ header discovery; `mpi.h` omitted by wrapper probe |
| `cmake --preset mpi-hdf5-fftw-debug -DMPI_CXX_HEADER_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include` | PASS |
| `cmake --build --preset build-mpi-hdf5-fftw-debug -j4` | PASS; 398 build steps |
| focused MPI-feature non-launch tests | PASS; 8/8 |
| focused np2 TreePM and hydro tests inside sandbox | ENVIRONMENT BLOCKED; PMIx listener socket creation denied |
| same focused np2 tests outside sandbox | PASS; 2/2 |

## Repository facts discovered during orientation

- `core::HierarchicalTimeBinScheduler` owns bins, activation ticks, pending
  transitions, and active sets. Particle/cell `time_bin` lanes are mirrors.
- Particle and gas-cell schedulers are already separate and gas-cell scheduler
  identity is keyed by stable `gas_cell_id`.
- `StepOrchestrator` dispatches typed stage buckets but every production
  callback still receives a broad mutable `StepContext::state` reference.
- `StageScheduler::schedule` returns a fixed KDK table; it does not build a
  resource/freshness/dependency plan.
- The workflow predicts inactive particle sources to TreePM evaluation epochs
  from persisted `last_drift_*` sidecars, but active drift/kick factors still
  use the current global base-step interval.
- Scheduler `appendElements` exists, but source-created particle registration is
  post-hoc coverage growth rather than one atomic identity/scheduler transaction.
- Coarsening currently throws when the requested larger bin is not aligned;
  the pending transition is not safely deferred.
- Adaptive criteria are represented by a narrow view, but workflow update
  functions still traverse broad state after every substep and rebuild hook
  registries inside loops.
- `TransientStepWorkspace` and active-index vectors are created/copied inside
  the production substep loop.
- `src/workflows/reference_workflow.cpp` is 7,589 lines and owns IC selection,
  MPI migration protocols, gravity, hydro/AMR, output/restart, and lifecycle.
- `core::module_registry` is currently only a static module-name list; it does
  not describe runtime contributions or capability prerequisites.
- External IC import is single-file and every rank calls it independently.
  The reader assumes `kpc`, `Msun`, `km/s`, comoving coordinates, and peculiar
  velocities; PartType2/3/5 silently map to dark matter; header/cosmology,
  multifile, high-word counts, field shapes/types, duplicate IDs, and
  distributed materialization are not adequately validated.
- Restart schema v20 is exact only for the same MPI topology and rank-local
  files. Rank-count changes are explicitly rejected.
- Output cadence is step-modulo only. Snapshot verification is materially
  stronger than count-only in current dirty work, but remains embedded in the
  workflow rather than an output service.
- The top-level workflow catches rank-local exceptions only after many possible
  collectives; no general phase-boundary failure coordinator prevents peers
  from entering a later collective.
- The main loop checks `current_time_code < t_code_end` but does not clip
  `dt_time_code`, so the last step can overshoot.

## Ownership and numerical invariants

1. Scheduler persistent state remains the sole bin/activation authority.
2. Per-element integration epochs must be keyed by stable identity, migrated
   with the element, persisted in restart, and mirrored nowhere as authority.
3. Particle positions remain comoving; particle velocities remain physical
   peculiar velocity `u = a dx/dt`; TreePM returns scale-free `A`.
4. Each active element must use drift/kick/Hubble-drag integrals over its own
   stored epoch interval. Inactive sources are predicted, not committed, to a
   force-evaluation epoch.
5. Refinement takes effect at the current safe boundary; coarsening remains
   pending until its larger period is aligned.
6. Gas wake-up may mutate only owner scheduler truth. Imported ghosts are
   read-only carriers; owner-directed wake requests are reconciled before the
   next stage/collective boundary.
7. IC conversion is driven by a versioned manifest. No unit, `h`, scale-factor,
   frame, velocity, thermodynamic, or species convention may be inferred from
   a dialect name alone.
8. Distributed startup materializes only assigned chunks/records on a rank;
   global count/species/mass/ID checks are reductions over local truth.
9. Restart/output/failure transactions either commit coherently or remain
   uncommitted on every rank.
10. Endpoint and output events clip the base interval before KDK preparation;
    a required event is never crossed.

## Milestones

- [x] M0: establish command-backed CPU and MPI-feature baseline.
- [x] M0: add a machine-readable runtime capability report and fail-closed
  compatibility gates.
- [ ] M1: decompose workflow policy into injected runtime services and a
  descriptor-driven composition root while preserving rung-zero behavior.
- [ ] M1: enforce stage/resource access with generation/epoch-stamped narrow
  views and a concrete task/resource plan.
- [ ] M2: implement versioned canonical IC manifest, explicit conversions,
  species policy, validation, multifile discovery, and distributed loading.
- [ ] M3: implement production per-element hierarchical cosmological KDK,
  active-only criteria plans, deferred coarsening, hydro limiter/wake-up,
  atomic scheduler registration, and endpoint/event clipping.
- [ ] M4: persist new temporal/event/wake state, implement supported
  rank-remappable restart, field-aware snapshot verification, event cadence,
  and rank-consistent failure coordination.
- [ ] M5: remove/bound hot-path allocation and copies, add work/communication
  counters, and demonstrate active-fraction and distributed-startup behavior.
- [ ] Final: run focused tests per milestone, broad preset matrix, hygiene,
  `git diff --check`, and adversarial self-review.

## Files changed per milestone

- M0: `docs/repair/chui_runtime_brain_exec_plan.md`,
  `include/cosmosim/workflows/runtime_capabilities.hpp`,
  `src/workflows/runtime_capabilities.cpp`,
  `include/cosmosim/workflows/reference_workflow.hpp`,
  `src/workflows/reference_workflow.cpp`, `tests/unit/test_runtime_capabilities.cpp`,
  and `CMakeLists.txt`.
- M3 scheduler/endpoint slice: `include/cosmosim/core/time_integration.hpp`,
  `src/core/time_integration.cpp`, `tests/unit/test_time_integration.cpp`,
  `docs/time_integration.md`, `include/cosmosim/workflows/reference_workflow.hpp`,
  `src/workflows/reference_workflow.cpp`,
  `tests/integration/test_reference_workflow_end_to_end.cpp`, and
  `docs/reference_workflow.md`.
- M1 services/registry slice: `include/cosmosim/workflows/runtime_services.hpp`,
  `src/workflows/runtime_services.cpp`, `tests/unit/test_runtime_services.cpp`,
  `include/cosmosim/core/module_registry.hpp`, `src/core/module_registry.cpp`,
  `tests/unit/test_module_registry.cpp`, `docs/architecture/overview.md`,
  `src/workflows/reference_workflow.cpp`, and `CMakeLists.txt`.
- M2 manifest/validation slice: `include/cosmosim/io/ic_reader.hpp`,
  `src/io/ic_manifest.cpp`, `src/io/ic_reader.cpp`,
  `tests/unit/test_ic_reader.cpp`, `docs/ic_reader_compatibility.md`,
  `src/workflows/runtime_capabilities.cpp`, and `CMakeLists.txt`.

## Unresolved blockers

- No external `CHUI_adversarial_brain_audit*.md` file was present in or adjacent
  to the repository; the checked-out code and requested prompt are the evidence
  map.
- MPI tests require an unsandboxed PMIx listener environment. This is an
  execution-environment constraint; the focused np2 baseline passes there.
- Optional source/AMR local-timestep combinations must remain fail-closed until
  their conservation, wake-up, and restart tests pass; DMO and the current
  supported fixed-grid hydro path are the minimum feature target.

## BRAIN-001 through BRAIN-019 closure matrix

Initial classifications below are baseline facts and will be updated with
source/test evidence as milestones close.

| ID | Status | Baseline evidence / required closure |
| --- | --- | --- |
| BRAIN-001 | still open | Production throws for `hierarchical_max_rung != 0`; global factors are applied to active rows. |
| BRAIN-002 | partially closed | Criteria submit scheduler candidates, but production is rung-zero and uses a fixed base interval. |
| BRAIN-003 | partially closed | A validated versioned `IcManifest` now owns field units, h/a exponents, frame, velocity convention, semantics, header/cosmology, and dispositions. The direct bridge uses a named v1 contract, but typed-config manifest-file selection and canonical converter tooling remain open. |
| BRAIN-004 | partially closed | PartType2/3 now reject by default and require distinct explicit mapping policies; PartType5 maps to BH and fills required sidecars. Complete PartType4 stellar field import and auxiliary family preservation remain open. |
| BRAIN-005 | still open | The direct path now rejects a multifile manifest instead of replicating ambiguous reads, but coordinated file/chunk assignment and distributed materialization are not implemented. |
| BRAIN-006 | still open | `reference_workflow.cpp` is a 7,589-line policy owner. |
| BRAIN-007 | partially closed | Contracts and compact views exist, but callbacks retain mutable full-state access. |
| BRAIN-008 | partially closed | Post-step criteria now evaluate scheduler-active rows and no longer construct per-element hook registries; criteria-view metadata construction still performs broad work and needs a persistent compiled plan. |
| BRAIN-009 | closed | Unaligned coarsening remains pending, survives scheduler restart state, and commits only at the next shared destination-period boundary; focused scheduler tests cover both paths. |
| BRAIN-010 | still open | No production neighbor timestep limiter/wake-up owner protocol. |
| BRAIN-011 | partially closed | Scheduler append exists; creation is not an atomic owner transaction. |
| BRAIN-012 | still open | v21 restart explicitly rejects rank-count changes. |
| BRAIN-013 | closed for the current snapshot schema | Production verification is stable-ID keyed and checks schema/file kind, normalized config and provenance identity, time/redshift/box/cosmology headers, every exported phase-space/mass/species/softening lane, and tracer fields; mismatch fails the output step. Fields not represented by `gadget_arepo_v4` remain restart-only rather than being falsely claimed as snapshot coverage. |
| BRAIN-014 | partially closed | A shared `FailureCoordinator` now gates restart-topology validation and preserves the local cause; remaining collective phases must migrate to it. |
| BRAIN-015 | partially closed | The compile-time registry now describes config, state/sidecars, tasks, criteria, restart, migration, diagnostics, prerequisites, and incompatibilities; composition-root activation is not yet descriptor-driven. |
| BRAIN-016 | closed | One process `RuntimeServices` bundle injects the root MPI context into drift, gravity, TreePM construction, and hydro; only the composition root constructs `MpiContext`. |
| BRAIN-017 | closed for periodic code-time events | Typed config now supports step-modulo and code-time cadence. Production clips at ordered code-time events, persists/restores the next event and pre-clip next-step proposal in restart v21, hashes the new state, and reads v20 with the lane explicitly disabled. Arbitrary explicit redshift/scale-factor lists remain outside this bounded closure. |
| BRAIN-018 | partially closed | The workflow reuses one capacity-preserving step workspace and consumes scheduler active spans without index copies, with profiler counters; remaining solver/internal staging allocation bounds need evidence. |
| BRAIN-019 | closed | Production clips the KDK interval before preparation, reports the committed endpoint, and emits `time.endpoint_clip`; integration coverage verifies no overshoot. |

## Milestone command log

Focused baseline CTest selection:

```text
unit_time_integration
integration_time_integration_loop
integration_hierarchical_time_bins
integration_hierarchical_timestep_regression
integration_reference_workflow
unit_ic_reader
integration_stage6_active_views
integration_restart_equivalence_multirate_bins
```

Result: 8/8 passed on `cpu-only-debug`. On the MPI feature build, the eight
non-launch focused tests passed and the np2 TreePM/hydro workflow tests passed
outside the PMIx socket-restricted sandbox.

Capability-report slice:

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_runtime_capabilities
ctest --test-dir build/cpu-only-debug -R '^unit_runtime_capabilities$' --output-on-failure
```

Result: configuration/build passed and `unit_runtime_capabilities` passed. Each
workflow run now writes `runtime_capabilities.json`; an explicitly requested
unsupported hierarchical runtime is rejected through the same capability
contract rather than being silently advertised.

Ordered output-event and restart-v21 slice:

```text
cmake --build --preset build-cpu-debug -j4 --target test_unit_restart_checkpoint_schema test_integration_reference_workflow
TMPDIR=/tmp ctest --test-dir build/cpu-only-debug -R '(unit_restart_checkpoint_schema|integration_reference_workflow)$' --output-on-failure
cmake --build --preset build-hdf5-debug -j4 --target test_integration_reference_workflow test_integration_restart_checkpoint_roundtrip test_unit_restart_checkpoint_schema
TMPDIR=/tmp ctest --test-dir build/hdf5-debug -R '(integration_reference_workflow|integration_restart_checkpoint_roundtrip|unit_restart_checkpoint_schema)$' --output-on-failure
```

Result: CPU 2/2 and HDF5 3/3 passed. The HDF5 production-path test clips a
`0.0002` proposal onto `0.00015` code-time events, verifies checkpointed next
event and restored `0.0002` next-step proposal, and resumes exactly to the
configured endpoint.

Final validation floor (2026-07-16):

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
TMPDIR=/tmp ctest --preset test-cpu-debug --output-on-failure
# 125/125 passed

cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4
TMPDIR=/tmp ctest --test-dir build/hdf5-debug --output-on-failure
# 127/127 passed

cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug -j4
TMPDIR=/tmp ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
# 127/127 passed

cmake --preset mpi-hdf5-fftw-debug -DMPI_CXX_HEADER_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include
cmake --build --preset build-mpi-hdf5-fftw-debug -j4
TMPDIR=/tmp ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure
# 153/153 passed outside the PMIx socket-restricted sandbox

cmake --preset asan-debug
cmake --build build/asan-debug -j4
TMPDIR=/tmp ctest --test-dir build/asan-debug --output-on-failure
# LeakSanitizer environment failure: "LeakSanitizer does not work under ptrace"
ASAN_OPTIONS=detect_leaks=0 TMPDIR=/tmp ctest --test-dir build/asan-debug --output-on-failure
# 125/125 passed with AddressSanitizer and UndefinedBehaviorSanitizer active;
# leak detection remains environment-blocked by ptrace

./scripts/ci/check_repo_hygiene.sh
git diff --check
# passed
```

The first combined CPU configure/build/test invocation hit its 60-second command
wrapper after tests 1--85 had passed. The complete CTest preset was immediately
rerun with a larger command budget and passed 125/125; the timeout is not counted
as a test failure.

Scheduler transition and endpoint slice:

```text
cmake --build --preset build-cpu-debug -j4 --target test_unit_time_integration test_integration_hierarchical_time_bins test_integration_hierarchical_timestep_regression
ctest --test-dir build/cpu-only-debug -R '^(unit_time_integration|integration_hierarchical_time_bins|integration_hierarchical_timestep_regression)$' --output-on-failure
cmake --build --preset build-cpu-debug -j4 --target test_integration_reference_workflow
TMPDIR=/tmp ctest --test-dir build/cpu-only-debug -R '^integration_reference_workflow$' --output-on-failure
```

Result: 3/3 scheduler tests and the reference-workflow integration test passed.
The first workflow-test invocation without `TMPDIR=/tmp` was environment-blocked
because `std::filesystem::temp_directory_path()` resolved to the read-only host
mount `/mnt/g/Temp`; no code failure occurred.

Runtime-services, registry, active-criteria, and IC-manifest slices:

```text
cmake --build --preset build-cpu-debug -j4 --target test_unit_module_registry test_unit_runtime_services test_integration_reference_workflow
TMPDIR=/tmp ctest --test-dir build/cpu-only-debug -R '^(unit_module_registry|unit_runtime_services|integration_reference_workflow)$' --output-on-failure
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4 --target test_unit_ic_reader
TMPDIR=/tmp ctest --test-dir build/hdf5-debug -R '^unit_ic_reader$' --output-on-failure
```

Result: the two service/registry units and production workflow integration
passed; the HDF5 IC unit passed with manifest conversion, PartType2 explicit
policy, PartType5 sidecar, duplicate-ID, malformed-shape, high-word, and
cosmology-mismatch coverage.
