# Campaign B final acceptance execution plan

## Mode and baseline

Mode: repair plus final Campaign B acceptance.

Baseline snapshot: `Cosmology-Simulator-main(28)(1).zip`.

Baseline findings reproduced before editing:

- the recursive hygiene check rejected two encoded Windows ADS sidecars;
- supplied-manifest inspection still resolved missing fields from the live configuration before comparing the manifest;
- distributed batch readers rotated by global batch index;
- source-session hashing was path-based and only performed before payload reads;
- manifest digest and several post-exchange accounting operations could throw outside rank-consistent phases;
- canonical single-file output truncated `NumPart_ThisFile` above `UINT32_MAX`;
- `src/io/ic_reader_file_set.cpp` contained 4,136 lines.

## Scientific and ownership invariants

- A schema-v4 supplied manifest is the only scientific authority for source order, hashes, field aliases and conversion contracts, species policy, missing-field policies, and replacement values.
- Runtime configuration may name a source in manifest mode only as a compatibility consistency check; it cannot override manifest science.
- Every rank reaches the same failure vote before a later MPI collective.
- Each source file has one deterministic reader rank in the current exact-audit implementation.
- Full-file hashing is bounded by source count, not routing-batch count.
- Payload reads remain bounded by configured chunk and staging capacities.
- Exact count, mass, owner, route-loss/duplication, and duplicate-ID audits remain enabled.
- A canonical one-file artifact is never finalized when a local family count cannot be represented by `NumPart_ThisFile`.
- MPI acceptance remains provisional unless actual np1/np2/np4 tests execute successfully.

## Milestones

1. Make manifest missing-field contracts authoritative and cover the production startup path.
2. Guard digest, serialization, exchange, and source-to-final accounting phases and extend fault-injection seams.
3. Assign each source to a stable reader rank, reuse one persistent session, and expose hash/open/reader/exchange counters.
4. Bind path identity to the open HDF5 session with size, timestamp, native device/inode where available, and start/end SHA-256 validation.
5. Add fail-closed canonical single-file count validation at `UINT32_MAX + 1` before HDF5 creation.
6. Split the IC subsystem into durable translation units with no source file above 2,000 lines.
7. Correct capability and user documentation to match the implemented guarantees and environmental limits.
8. Run supported CPU/HDF5 validation, record blocked MPI commands, remove generated artifacts, and package a recursively clean source ZIP.

## Validation log

Commands and exact outcomes are appended during the repair. Build trees and test outputs are excluded from the final source handoff.

## Closure matrix

The final section will classify every F0-F8 item as `closed`, `partially closed`, `environment-blocked`, or `still open` using command-backed evidence.

## Implemented repair record

### F0 — supplied-manifest authority

- `manifest_v1` resolves the source file set, order, hashes, selected aliases,
  conversion contracts, species policy, and missing-field contracts from the
  supplied schema-v4 manifest.
- Live configuration missing-field defaults are not re-applied while verifying
  an authoritative manifest.
- `mode.ic_file` is optional in manifest mode. When present, the production
  startup path treats it only as a source-set consistency check and rejects a
  conflict.
- Production-path integration coverage starts `InitialConditionRuntime` with no
  `mode.ic_file`, uses a manifest `use_config_value` contract while live config
  remains `reject`, compares intentionally equivalent direct/manifest state,
  and rejects an optional conflicting path.

### F1 — collective failure consensus

The distributed ingestion implementation now guards manifest digest,
serialization accounting, post-exchange accounting, exact source-to-final
accounting, duplicate-audit accounting, reader-session creation, and source
identity completion validation. Rank-selective fault seams are registered for
those phases with finite CTest timeouts. Static MPI-enabled syntax validation
was performed because no MPI C++ toolchain is installed here.

### F2/F3 — stable readers, bounded hashing, and source identity

- Reader ownership is deterministic by source-file index; batches for a file do
  not rotate among ranks.
- One persistent payload session is reused per nonempty source file.
- Complete-file SHA-256 work is bounded by three passes per source file:
  inspection, payload-session start, and payload-session completion.
- Sessions capture path size and modification time, plus POSIX device/inode
  where available, before/after open and hashing. Completion revalidates
  identity, size, and SHA-256.
- Documentation does not claim same-file-descriptor hashing; it states the
  narrower portable guarantee implemented by the session abstraction.

### F4 — distributed cost observability

The import report exposes routing batches, chunks, reader batches and records,
reader-record imbalance, main exchanges, exact-audit exchanges, bytes,
collective phases, routing collective phases, peak staging, and ingestion wall
clock nanoseconds. The MPI fixture requires main exchanges to scale with
routing batches rather than source chunks and caps routing logical phases at
32 per batch. Exact source/final and duplicate-ID audits remain mandatory; no
weaker default mode was introduced.

### F5 — canonical single-file count safety

Canonical conversion validates all six per-family local counts before creating
an HDF5 file. `UINT32_MAX` is accepted; `UINT32_MAX + 1` and mixed-family
boundary violations fail with a format-limit diagnostic and cannot create a
finalized artifact.

### F6 — IC subsystem decomposition

The former 4,136-line `ic_reader_file_set.cpp` was removed. Durable units now
separate file-set/manifest/schema inspection, serial ingestion, streaming
record production, distributed failure protocol, distributed audits, and
routing ingestion. Structural
CI rejects any `src/io/ic_*.cpp` above 2,000 lines, a replacement `ic_utils.cpp`,
or leakage of the narrow private common header. At this repair point the largest
units are `ic_file_set_manifest.cpp` (1,700 lines),
`ic_distributed_ingestion.cpp` (1,134 lines),
`ic_distributed_audit.cpp` (585 lines), and
`ic_failure_protocol.cpp` (313 lines).

### F7/F8 — truth and handoff

Capability text remains provisional for distributed IC ingestion until real MPI
execution. Configuration, schema, profiling, workflow, and output documentation
now state the manifest authority, reader/hash bounds, source-identity limits,
collective counters, and canonical large-count behavior. The two encoded ADS
sidecars were deleted, recursive source packaging was added, and the final ZIP
will be validated with the same denylist as the worktree.

## Commands and evidence captured so far

```text
bash scripts/ci/check_repo_hygiene.sh
  baseline: failed on two encoded Zone.Identifier sidecars

cmake --preset hdf5-debug
  passed

focused HDF5 build and tests
  unit_ic_reader: passed
  unit_runtime_capabilities: passed
  integration_initial_condition_manifest_runtime: passed
  integration_ic_decomposition_structure: passed

MPI-enabled static syntax checks for:
  src/io/ic_distributed_ingestion.cpp: passed
  tests/integration/test_distributed_ic_reader_mpi.cpp: passed

cmake --preset mpi-hdf5-fftw-debug
  environment-blocked: MPI_CXX was not found
```

The final command log, complete closure matrix, archive hash, and blocked MPI
commands are appended after the supported CPU/HDF5 gates and source packaging.

## Final supported validation evidence

The final source state was configured, built, and tested with the dependency
sets available in this environment.

```text
cmake --preset cpu-only-debug
  passed; GNU 14.2.0, MPI/HDF5/FFTW disabled

cmake --build --preset build-cpu-debug
  passed; final incremental build reported no work after the complete rebuild

ctest --preset test-cpu-debug --output-on-failure -j1
  passed: 132/132 tests, 0 failures, 153.35 s

cmake --preset hdf5-debug
  passed; HDF5 1.14.5 detected

cmake --build --preset build-hdf5-debug
  passed; final incremental build reported no work

ctest --preset test-hdf5-debug --output-on-failure -j1
  passed: 136/136 tests, 0 failures, 431.03 s

focused HDF5 Campaign B selection
  passed: 7/7 tests, including InitialConditionRuntime manifest startup,
  reader/source-identity tests, converter tests, record codec, capability,
  configuration, and IC decomposition structure

MPI-enabled declaration-only syntax checks
  passed for ic_failure_protocol.cpp, ic_distributed_audit.cpp,
  ic_distributed_ingestion.cpp, and test_distributed_ic_reader_mpi.cpp

cmake --preset mpi-hdf5-fftw-debug
  environment-blocked during configuration: MPI_CXX was not found
```

A diagnostic `-j4` HDF5 run initially produced subprocess aborts in the two
transform-fuzz tests because those fixtures ran concurrently. Both tests passed
when isolated, and the required serial preset then passed all 136 tests. No
product assertion or Campaign B test failed in the final serial gate.

No MPI executable, FFTW/FFTW-MPI distributed path, parallel-HDF5 path, or
np1/np2/np4 distributed runtime test was executed. The declaration-only syntax
check is compilation evidence only and is not represented as MPI runtime
acceptance.

## Exact source changes

Modified:

```text
CMakeLists.txt
docs/campaign_b_acceptance_repair_report.md
docs/configuration.md
docs/ic_reader_compatibility.md
docs/output_schema.md
docs/profiling.md
docs/reference_workflow.md
docs/repair/campaign_b_completion_exec_plan.md
docs/repair_open_issues.md
docs/repair_state_recap.md
include/cosmosim/io/ic_reader.hpp
scripts/ci/check_repo_hygiene.sh
src/io/ic_reader.cpp
src/io/ic_reader_session.cpp
src/io/internal/ic_reader_session.hpp
src/workflows/initial_condition_runtime.cpp
src/workflows/runtime_capabilities.cpp
tests/integration/test_distributed_ic_reader_mpi.cpp
tests/unit/test_ic_reader.cpp
tools/convert_ic.cpp
```

Added:

```text
docs/repair/campaign_b_final_acceptance_exec_plan.md
scripts/ci/package_source_zip.sh
src/io/ic_distributed_audit.cpp
src/io/ic_distributed_ingestion.cpp
src/io/ic_failure_protocol.cpp
src/io/ic_file_set_manifest.cpp
src/io/ic_serial_ingestion.cpp
src/io/ic_stream_ingestion.cpp
src/io/internal/ic_canonical_limits.hpp
src/io/internal/ic_distributed_audit.hpp
src/io/internal/ic_failure_protocol.hpp
src/io/internal/ic_file_set_common.hpp
tests/integration/test_ic_decomposition_structure.cmake.in
tests/integration/test_initial_condition_manifest_runtime.cpp
```

Removed:

```text
src/io/ic_reader_file_set.cpp
docs/repair/CHUI_runtime_brain_prompt_result_adversarial_audit.md#Uf03aZone.Identifier
docs/repair/chui_runtime_decomposition_goal_experiment.md#Uf03aZone.Identifier
```

## Tests added or materially extended

- `integration_initial_condition_manifest_runtime`: real startup-path coverage
  for manifest-only configuration, authoritative missing-field values,
  direct/manifest equivalence, provenance authority, and conflicting optional
  `mode.ic_file` rejection.
- `unit_ic_reader`: source-path/hash/schema authority, dialect/config-independent
  missing-field behavior, source replacement/size/in-place mutation detection,
  unchanged-session success, manifest-hash mismatch, and canonical count
  boundaries.
- `integration_distributed_ic_reader_mpi`: stable reader/hash bounds, exchange
  and collective counters, reader imbalance, np1/np2/np4 modes, duplicate and
  route mutation cases, plus finite-time rank-selective failure phases.
- `integration_ic_decomposition_structure`: no old/replacement monolith, no IC
  source over 2,000 lines, required responsibility units, and private-header
  dependency limits.
- `integration_convert_ic`: existing production converter coverage continues to
  exercise bounded streaming and canonical bundle behavior.

## Final closure matrix

| ID | Status | Evidence and remaining condition |
|---|---|---|
| F0-1 | closed | Production startup accepts `manifest_v1` without `mode.ic_file`. |
| F0-2 | closed | Manifest source set/order/paths are used; optional config path is only a conflict check. |
| F0-3 | closed | Manifest missing-field policy/value contracts bypass live config defaults. |
| F0-4 | closed | Conflicting optional `mode.ic_file` fails in the real runtime path. |
| F1-1 | environment-blocked | Digest phase is guarded and fault test is registered/static-compiled; real MPI bounded termination remains unexecuted. |
| F1-2 | environment-blocked | Serialization/post-exchange accounting is guarded and fault tests compile; real MPI remains unexecuted. |
| F1-3 | environment-blocked | Source/final and global-ID audit accounting is guarded and fault tests compile; real MPI remains unexecuted. |
| F2-1 | environment-blocked | Stable per-file reader ownership is implemented; multi-rank runtime evidence remains unexecuted. |
| F2-2 | environment-blocked | Three full-file hash passes per nonempty source is encoded/counted/tested in the MPI fixture; real multi-rank bound remains unexecuted. |
| F2-3 | environment-blocked | Files map deterministically across ranks rather than rank zero; real reader-distribution evidence remains unexecuted. |
| F3-1 | closed | Persistent sessions capture and revalidate size/time/native identity and start/end SHA-256. |
| F3-2 | closed | Replacement, size change, identity change, detectable mutation, hash mismatch, and unchanged success are tested. |
| F4-1 | environment-blocked | Logical collective phases and exchange counters are implemented/static-compiled; runtime accuracy remains unexecuted. |
| F4-2 | environment-blocked | Fixture requires main exchanges to equal batches and batches to be fewer than chunks; MPI execution remains blocked. |
| F4-3 | environment-blocked | Strict count/mass/ownership/route/duplicate audits remain mandatory; distributed test execution remains blocked. |
| F5-1 | closed | Single-file local counts above `UINT32_MAX` fail before HDF5 creation; boundary tests pass. |
| F6-1 | closed | Former 4,136-line monolith removed; schema/manifest, serial, stream, transport, failure, and audit duties are split. |
| F6-2 | closed | Largest IC source is 1,700 lines; structural and hygiene gates enforce 2,000. |
| F7-1 | closed | Current and historical Campaign B documents identify authority, bounds, limitations, and supersession correctly. |
| F7-2 | closed | Distributed capability remains explicitly provisional pending real MPI. |
| F8-1 | closed | Source-staging hygiene uses the recursive worktree policy and rejects ADS/build/output/cache/binary content. |
| F8-2 | closed | Source-only packaging script validates ZIP integrity and archive denylist; final artifact evidence is reported with the handoff. |
| F8-3 | closed | Missing MPI toolchain and all unexecuted distributed claims are stated explicitly. |
| F8-4 | closed | Exact full and focused user-side MPI commands are supplied below. |

## User-side MPI acceptance commands

Run on a dependency-complete MPI/HDF5/FFTW installation:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1
```

Focused np1/np2/np4 authority/equivalence matrix:

```bash
timeout 300s ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader(_(canonical|manifest|alternate|type_policy))?_mpi_(1|2|4)_rank$'
```

Scaling, reconciliation, duplicate, and failure-consensus gates:

```bash
timeout 300s ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_scaling_mpi_(2|4)_rank$'

timeout 300s ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_(duplicate|route_loss|route_duplicate)_mpi_two_rank$'
timeout 600s ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_fault_.*_mpi_two_rank$'
```

## Final verdict

The Campaign B source repair, serial/HDF5 acceptance, maintainability split, and
clean-handoff implementation are closed. Final campaign acceptance is
**environment-blocked**, not source-blocked: the remaining evidence gate is the
registered real-MPI np1/np2/np4 correctness, fault-termination, and cost matrix.
