# Campaign B distributed final gates execution plan

Mode: repair and final Campaign B distributed acceptance.

Suggested PR title:

```text
fix-campaign-b-distributed-manifest-verification-consensus-and-handoff
```

Suggested local branch:

```text
fix/campaign-b-distributed-final-gates
```

## Baseline

The audited snapshot still:

- inferred supplied-manifest authority from a nonempty `source_sha256` vector;
- selected the distributed inspection dialect from runtime configuration;
- allowed exact-audit capacity/path accounting to throw outside consensus;
- reported logical failure-consensus phases as though they were actual MPI
  communicator calls;
- lacked a real supplied-canonical distributed fixture and actual-call formula;
- contained two encoded Windows ADS sidecars in the delivered ZIP.

The accepted serial scientific conversion, canonical converter, source identity,
stable file ownership, large-count guard, and IC decomposition were preserved.

## Scientific, ownership, and protocol invariants

- A non-null supplied manifest is validated and authoritative or rejected.
- The supplied manifest dialect controls distributed source inspection.
- Canonical supplied manifests execute embedded/sidecar/digest/marker verification.
- Runtime configuration never silently replaces supplied manifest science policy.
- Every potentially throwing local operation before a later collective is inside
  a rank-consistent phase; capacity arithmetic is checked and cannot wrap.
- Actual MPI calls and logical consensus phases are measured separately.
- Actual-call counters increment saturatingly and `noexcept` in one wrapper unit.
- Count, mass, ownership, route, source-to-final, and duplicate-ID audits remain exact.
- File ownership remains stable by source-file index and full hashes remain bounded.
- No IC implementation source exceeds 2,000 lines.
- The user-facing handoff is produced only by the verified packaging script.

## Implemented milestones

### G0 — Explicit manifest authority

- Serial and streaming entry points validate every non-null manifest.
- Distributed validation occurs inside `IC distributed configuration validation`
  before discovery or scientific ingestion.
- The old `manifest != nullptr && !source_sha256.empty()` heuristic is removed.
- Invalid/incomplete manifests fail instead of falling back to runtime authority.
- Direct API tests cover empty hashes, incomplete source metadata, unsupported
  dialect version, and mismatched per-file vector lengths.

### G1 — Distributed canonical manifest verification

- Distributed inspection uses `options.manifest->dialect` whenever a manifest is
  supplied; runtime convention is used only without a manifest.
- Canonical sources require the canonical schema marker and complete embedded,
  HDF5-digest, sidecar, and completion-marker verification.
- Bridge/canonical mismatches fail explicitly.
- Registered distributed cases cover:
  - valid supplied canonical manifests at np1/np2/np4;
  - direct canonical versus supplied canonical exact rank-local state equivalence;
  - tampered sidecar, missing marker, bad HDF5 digest;
  - source-path and source-hash mismatches;
  - canonical-manifest/bridge-source and bridge-manifest/canonical-source mismatch;
  - invalid supplied manifest.

### G2 — Failure consensus

The reviewed production phase order is:

1. collective configuration validation;
2. distributed file discovery and path broadcast/decode;
3. assigned-file inspection and hashing;
4. manifest-fragment gather, assembly, broadcast, and decode;
5. runtime-cosmology, serialization, digest, and routing configuration;
6. per routing batch: read/convert, coverage, owner serialization, accounting,
   main exchange, decode, exact source-to-final audit, append, staging accounting;
7. per-file session completion identity validation;
8. reader imbalance, final state, exact global duplicate-ID audit, scientific
   totals, and final routed-record accounting.

Exact source/final and duplicate-ID capacity work is guarded. Added finite-time
fault seams:

```text
exact_audit_capacity_accounting
duplicate_audit_capacity_accounting
audit_path_storage_accounting
audit_peak_staging_accounting
```

Capacity multiplication uses checked arithmetic; it cannot silently wrap.

### G3 — Truthful actual MPI-call accounting

Production IC collectives now pass through:

```text
src/io/ic_mpi_collectives.cpp
src/io/internal/ic_mpi_collectives.hpp
```

The wrapper counts actual calls by operation:

```text
Allreduce
Bcast
Gather
Gatherv
Alltoall
Alltoallv
```

It also tracks total, routing, and non-routing calls. Logical phase counters are
retained under explicit logical names; compatibility aliases are documented as
logical only. Structural CI rejects raw production IC collective calls outside
the wrapper.

Successful routing protocol v1 is exact:

```text
routing_mpi_collective_call_count
  = 30 * routing_batch_count
```

The 30 calls are 22 failure-consensus votes, three coverage reductions, two
`Alltoall`/`Alltoallv` pairs, and one source/final reconciliation reduction.

Successful non-routing protocol v1 is also exact:

```text
nonrouting_mpi_collective_call_count
  = 40
  + (validate_runtime_cosmology ? 1 : 0)
  + source_file_count
  + 10 * distributed_id_audit_round_count
  + mpi_bcast_call_count
```

The Bcast term remains measured because length-prefixed metadata may span more
than one 64 MiB payload chunk. Tests require:

- both formulas;
- total = routing + non-routing;
- per-operation sum = total;
- min/max agreement for total, routing/non-routing, each operation, logical
  phases, and duplicate-audit rounds across ranks;
- correct `collectives_per_million_records`, including zero-record behavior.

### G4 — Scientific integrity retained

No integrity mode was weakened. The path still enforces:

- global and per-species counts;
- global mass totals;
- ownership completeness and exclusivity;
- route-loss and route-duplication detection;
- exact source-to-final ID multiset reconciliation;
- global duplicate-ID rejection;
- start/end source identity and SHA-256 validation;
- stable per-file readers and the three-full-hash-pass source ceiling.

### G5 — Documentation and capability truth

Configuration, IC compatibility, output schema, profiling, workflow, repair-state,
and open-issue documentation now distinguish:

- supplied manifest authority;
- manifest-selected distributed dialect;
- canonical bundle verification;
- logical phases versus actual calls;
- per-operation counters and versioned cost formulas;
- duplicate-audit rounds and normalized collective cost;
- provisional MPI status where execution is unavailable.

Distributed capability remains `provisional`; no MPI runtime pass is claimed.

### G6 — Handoff hygiene

- Removed both encoded ADS sidecars.
- Strengthened `scripts/ci/package_source_zip.sh` to stage, recursively audit,
  archive, run `unzip -t`, enforce one repository root, apply the archive
  denylist, count entries, and print SHA-256.
- Strengthened repository hygiene to require the MPI wrapper and reject raw
  production collectives.

## Files changed

### Added

```text
docs/repair/campaign_b_distributed_final_gates_exec_plan.md
src/io/ic_mpi_collectives.cpp
src/io/internal/ic_mpi_collectives.hpp
```

### Modified

```text
CMakeLists.txt
docs/configuration.md
docs/ic_reader_compatibility.md
docs/output_schema.md
docs/profiling.md
docs/reference_workflow.md
docs/repair_open_issues.md
docs/repair_state_recap.md
include/cosmosim/io/ic_reader.hpp
scripts/ci/check_repo_hygiene.sh
scripts/ci/package_source_zip.sh
src/io/ic_distributed_audit.cpp
src/io/ic_distributed_ingestion.cpp
src/io/ic_failure_protocol.cpp
src/io/ic_file_set_manifest.cpp
src/io/ic_serial_ingestion.cpp
src/workflows/initial_condition_runtime.cpp
tests/integration/test_distributed_ic_reader_mpi.cpp
tests/integration/test_ic_decomposition_structure.cmake.in
tests/unit/test_ic_reader.cpp
```

### Removed

```text
docs/repair/CHUI_runtime_brain_prompt_result_adversarial_audit.md#Uf03aZone.Identifier
docs/repair/chui_runtime_decomposition_goal_experiment.md#Uf03aZone.Identifier
```

## Validation evidence

### Hygiene and static review

- No remaining `Zone.Identifier`/encoded ADS files.
- No old manifest-authority hash heuristic.
- No raw production `MPI_Allreduce`, `MPI_Bcast`, `MPI_Gather`,
  `MPI_Gatherv`, `MPI_Alltoall`, or `MPI_Alltoallv` outside
  `ic_mpi_collectives.cpp`.
- All `src/io/ic_*.cpp` units remain below 2,000 lines; largest is
  `ic_file_set_manifest.cpp` at 1,707 lines.
- Packaging and hygiene scripts pass `bash -n`.

### MPI-enabled declaration compilation

Using a local declaration-only MPI header because no system MPI C++ toolchain is
installed, the following five translation units passed `-fsyntax-only` with
MPI/HDF5 enabled:

```text
src/io/ic_mpi_collectives.cpp
src/io/ic_failure_protocol.cpp
src/io/ic_distributed_audit.cpp
src/io/ic_distributed_ingestion.cpp
tests/integration/test_distributed_ic_reader_mpi.cpp
```

This is compilation evidence only, not runtime MPI acceptance.

### CPU-only

```text
cmake --preset cpu-only-debug
PASS

cmake --build --preset build-cpu-debug
PASS

CPU test inventory
132/132 passed across bounded serial invocations.
The long TreePM-Ewald test was executed directly with its existing tolerance and
returned status 0 after the all-suite command window expired.
```

After the final counter/report additions, the complete CPU tree rebuilt and the
focused capability/docs/structure set passed 3/3.

### HDF5

```text
cmake --preset hdf5-debug
PASS; HDF5 1.14.5

cmake --build --preset build-hdf5-debug
PASS

focused Campaign B/capability/docs/structure set
7/7 passed
```

The broad HDF5 inventory produced 135/136 successful results with zero observed
assertion/test failures. `integration_runtime_app_smoke` exceeded its existing
300-second CTest timeout again; it is unchanged and was already identified as a
long-running incomplete test in the input audit. No timeout or tolerance was
weakened.

### MPI preset

```text
cmake --preset mpi-hdf5-fftw-debug
BLOCKED
```

Exact configure failure:

```text
Package 'mpi-cxx', required by 'virtual:world', not found
Could NOT find MPI_CXX
Could NOT find MPI (missing: MPI_CXX_FOUND CXX)
```

No np1/np2/np4, malformed-input, fault-injection, counter, or scaling runtime
pass is claimed.

## Adversarial self-review

- Supplied manifest authority no longer depends on hash-vector contents.
- Distributed canonical inspection cannot select bridge semantics from runtime
  convention.
- Invalid manifests cannot silently become absent manifests.
- Canonical supplied fixtures use actual canonical embedded/sidecar/marker
  artifacts.
- Source path/hash, dialect, sidecar, marker, and digest mismatches are registered.
- Capacity/path/peak calculations between collectives are guarded and checked.
- Test-harness MPI calls occur after report capture and do not contaminate
  production counters.
- Actual-call formulas are protocol-derived, not arbitrary ceilings.
- Scientific count/mass/ownership/route/ID checks remain in production.
- Stable reader assignment and hash bounds were not changed.
- No replacement IC monolith or private-header leak was introduced.
- Distributed capability remains provisional.
- The final archive must be the exact packaging-script output.

## Final closure matrix

| ID | Status | Evidence |
|---|---|---|
| G0-1 | closed | Every non-null manifest validated at serial/stream/distributed entry |
| G0-2 | closed | Invalid/incomplete manifests reject; no fallback heuristic |
| G1-1 | closed | Distributed inspection selects supplied `manifest.dialect` |
| G1-2 | closed | Canonical supplied path executes canonical bundle verifier |
| G1-3 | closed | Bridge/canonical mismatches registered and fail closed |
| G1-4 | environment-blocked | np1/np2/np4 tests exist; runtime unavailable |
| G2-1 | closed | Source/final capacity calculations guarded and checked |
| G2-2 | closed | Duplicate-audit path/storage/peak calculations guarded |
| G2-3 | environment-blocked | 30-second fault tests exist; MPI launch unavailable |
| G3-1 | closed | Logical counters explicitly named; aliases documented |
| G3-2 | closed | Central actual-call wrapper covers production IC collectives |
| G3-3 | environment-blocked | Per-kind assertions exist; MPI runtime unavailable |
| G3-4 | environment-blocked | Per-counter rank agreement assertions exist |
| G3-5 | environment-blocked | Exact routing/non-routing formulas registered |
| G3-6 | closed | Metric implemented with zero-record behavior |
| G4-1 | partially closed | Source contract retained; distributed runtime evidence blocked |
| G4-2 | partially closed | Stable reader/hash assertions registered; runtime blocked |
| G4-3 | environment-blocked | Scaling test exists; MPI runtime unavailable |
| G5-1 | closed | Docs match source and limitations |
| G5-2 | closed | Distributed capability remains provisional |
| G6-1 | closed | Source-tree hygiene prepared; final rerun occurs before packaging |
| G6-2 | closed | Packaging script is the only handoff path |
| G6-3 | closed | Final clean archive produced after generated-tree cleanup |
| G6-4 | closed | Script performs `unzip -t`, root, and denylist checks |
| G6-5 | closed | Packaging output reports SHA-256 |
| G7-1 | closed | CPU 132/132 passed across bounded invocations |
| G7-2 | partially closed | 135/136 HDF5; runtime smoke hit unchanged 300-second timeout |
| G7-3 | environment-blocked | `MPI_CXX` unavailable; syntax compilation passed |
| G7-4 | closed | Exact commands listed below |

## User-side MPI acceptance commands

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1
```

Focused supplied-canonical and general distributed matrix:

```bash
timeout 300s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_.*_mpi_(1|2|4)_rank$'
```

Fault-consensus matrix:

```bash
timeout 600s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_fault_.*_mpi_two_rank$'
```

Malformed manifest, canonical bundle, route, and duplicate cases:

```bash
timeout 300s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_(canonical_manifest_.*|bridge_manifest_.*|invalid_manifest|duplicate|route_loss|route_duplicate)_mpi_two_rank$'
```

Scaling and counter formula cases:

```bash
timeout 300s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_(canonical_manifest|scaling)_mpi_(1|2|4)_rank$'
```

## Final verdict

Source-side Campaign B final-gate defects are closed. Distributed runtime
acceptance is environment-blocked until the registered MPI matrix executes on a
host with a working MPI C++ toolchain. The capability must remain provisional
until that evidence exists.
