# Campaign B completion execution plan

## Task classification

Mode: **REPAIR + CAMPAIGN B COMPLETION + PACKAGING**.

The repository snapshot is authoritative. Existing Campaign A/B architecture is preserved. The repair prioritizes scientific correctness and fail-closed provenance over apparent feature completeness.

## Initial evidence

- `AGENTS.md` read and applied.
- Repository was supplied as a source archive without Git metadata; no destructive Git operation is possible or used.
- Baseline hygiene failed on two encoded Windows ADS sidecars under `docs/repair/`.
- The audited IC implementation is concentrated in `src/io/ic_reader_file_set.cpp` (4,032 lines).
- Confirmed baseline defects:
  - unsigned-only reads for `NumPart_ThisFile` and `NumFilesPerSnapshot`;
  - velocity-convention powers applied to `InternalEnergy` and `BH_Mdot`;
  - missing gas/star/BH fields silently materialized from implicit zero/current-state defaults;
  - canonical digest syntax checked without recomputing a manifest digest;
  - converter finalizes HDF5 before its manifest and materializes a complete `SimulationState`;
  - repeated file/group/dataset opens per source chunk;
  - per-chunk distributed exchange and failure-vote phases;
  - manifest mode redundantly requires `mode.ic_file`;
  - submitted archive contains unsafe sidecar names and generated validation outputs.

## Milestones

### B0 — Header compatibility

- Add checked signed-or-unsigned readers for AREPO count attributes.
- Keep total low/high word attributes unsigned.
- Reject negative, overflow, non-integer, malformed-rank, and malformed-extent inputs.
- Preserve actual HDF5 datatype metadata in manifest fields for count/header attributes.

### B1 — Field-specific conversion semantics

- Replace semantics-only special handling with explicit per-field contracts.
- Apply velocity storage convention only to `Velocities`.
- Define independent `InternalEnergy` and `BH_Mdot` units/exponents.
- Bump manifest semantic schema and provide compatibility diagnostics.
- Add independent golden conversion tests at nontrivial `a`, `h`, and unit scales.

### B2 — Missing-field policies

- Add typed missing-field policies to normalized configuration and manifests.
- Fail closed by default for gas thermodynamics and star/BH sidecars.
- Permit defaults only through explicit recorded policy/value pairs.
- Reject unsupported reconstruction/preservation modes before payload reads.
- Reject unsupported converter tracer output during CLI validation.
- Diagnose negative stellar formation times as unsupported wind particles.

### B3 — Canonical bundle integrity

- Bind canonical HDF5 to normalized manifest JSON and SHA-256.
- Recompute and verify the embedded/sidecar manifest digest at read time.
- Finalize HDF5, manifest, and completion marker as a recoverable bundle.
- Clean orphaned partial/final files on injected or ordinary failure.

### B4 — Streaming converter

- Introduce a source-chunk callback/session API so the converter writes extensible datasets directly.
- Keep authoritative converter staging bounded by configured chunk size.
- Maintain bounded count/mass/statistics/ID audits and deterministic source order.
- Retain full-state import only as the runtime/explicit small-input path.

### B5 — Persistent sessions and batched routing

- Cache file/group/dataset handles per reader session.
- Batch several chunks within the configured staging bound before one routing exchange.
- Add open/exchange/collective counters and deterministic batch-size tests.
- Keep failure-prone operations within coordinated phases.

### B6 — Manifest authority and source identity

- Let manifest mode derive source files without requiring a redundant `mode.ic_file`.
- Reject a conflicting explicitly supplied source path.
- Revalidate source identity (size/mtime/hash) at reader-session start and payload completion.
- Document deterministic x-slab ingestion ownership and required downstream rebalance.

### B7 — IC subsystem decomposition

- Split schema/conversion/session/codec/serial/distributed/audit responsibilities into narrow internal units.
- Add line-count/dependency guardrails.

### B8 — Documentation and packaging

- Update configuration, schema, compatibility, distributed-ingestion, capability, repair-state, and validation documentation.
- Strengthen recursive hygiene checks.
- Remove generated output/ADS/cache/build artifacts and create a source-only ZIP.

## Validation order

1. Focused HDF5 unit and converter tests.
2. Config/parser and capability tests.
3. HDF5 preset build/test.
4. CPU preset build/test.
5. MPI preset configure/build/test only when wrappers and dependencies are available.
6. Repository hygiene, archive listing audit, `unzip -t`, and SHA-256.

## Implementation result

### Files and responsibility changes

The repair preserved the existing owner-service and distributed-ingestion design while adding narrow IC units for:

- `src/io/ic_conversion_catalog.cpp`: field-specific dimensional/cosmological contracts;
- `src/io/ic_canonical_bundle.cpp`: embedded/sidecar manifest binding and completion-marker verification;
- `src/io/ic_reader_session.cpp`: persistent source-file and dataset handles plus source-identity revalidation;
- `src/io/ic_byte_codec.cpp`: reusable endian-safe scalar coding;
- `src/io/ic_record_codec.cpp`: independently testable 168-byte distributed particle record;
- `src/io/ic_sha256.cpp`: source and manifest SHA-256 implementation.

`src/io/ic_reader_file_set.cpp` decreased in responsibility but remains 4,136 lines. The hygiene guard rejects a replacement `ic_utils` dumping ground, limits the remaining file to 4,200 lines, and limits split IC units to 1,200 lines. Full separation of schema inspection, serial ingestion, distributed ingestion, and distributed audit remains a maintainability follow-up.

### Validation commands and exact outcomes

Focused HDF5 build and Campaign B tests:

```bash
cmake --build --preset build-hdf5-debug \
  --target test_unit_ic_reader test_unit_ic_record_codec \
           test_unit_config_parser test_unit_runtime_capabilities \
           cosmosim_convert_ic -j5
ctest --test-dir build/hdf5-debug --output-on-failure \
  -R '^(unit_config_parser|unit_runtime_capabilities|unit_ic_reader|unit_ic_record_codec|integration_convert_ic)$'
```

Result: **5/5 passed**. This includes real signed AREPO-style headers, malformed/negative/overflow integer rejection, independent conversion golden values, explicit missing-field policies, verified canonical bundles, converter streaming bounds, duplicate-ID rejection, relative manifest authority, tracer early rejection, and six bundle-failure rollback seams.

Focused CPU build and tests:

```bash
cmake --build --preset build-cpu-debug \
  --target test_unit_config_parser test_unit_runtime_capabilities -j5
ctest --test-dir build/cpu-only-debug --output-on-failure \
  -R '^(unit_config_parser|unit_runtime_capabilities)$'
```

Result: **2/2 passed**.

Full CPU build attempt:

```bash
cmake --build --preset build-cpu-debug
```

Result: progressed without a compiler diagnostic to **37/334** build actions, then the command was stopped by the 20-minute execution limit. The full CPU CTest matrix was therefore not run.

Full HDF5 build attempts:

```bash
cmake --build --preset build-hdf5-debug -j2
cmake --build --preset build-hdf5-debug -j8
```

Result: both progressed without a compiler diagnostic but were stopped by 20-minute execution limits; the second attempt reached **39/320** remaining actions. The complete HDF5 CTest matrix was therefore not run.

MPI configuration attempt:

```bash
cmake --preset mpi-hdf5-fftw-debug
```

Result: **environment-blocked** because no `mpicxx`/`MPI_CXX` wrapper or `mpirun` was available. The MPI build and runtime tests were not claimed. As an additional source-only check, `src/io/ic_reader_file_set.cpp` and `tests/integration/test_distributed_ic_reader_mpi.cpp` both passed `c++ -fsyntax-only` with `COSMOSIM_ENABLE_MPI=1`, HDF5 enabled, and a temporary declaration-only MPI API header. This is not runtime evidence.

Repository hygiene after deleting validation build trees and Python caches:

```bash
bash scripts/ci/check_repo_hygiene.sh
```

Result: **passed**. Archive integrity and listing checks are recorded in the final handoff response.

## Campaign B closure matrix

| ID | Status | Closure evidence / remaining boundary |
|---|---|---|
| B0-1 | closed | Signed and unsigned AREPO per-file/file-count attributes use checked wide reads; signed-int32 and mixed-multifile fixtures pass. |
| B0-2 | closed | Negative, overflow, floating, malformed-rank, and malformed-extent fixtures fail closed. |
| B1-1 | closed | Manifest v4 validates velocity convention power from velocity semantics only. |
| B1-2 | closed | `InternalEnergy` uses an independent specific-energy contract and is invariant under velocity-storage convention changes. |
| B1-3 | closed | `BH_Mdot` uses an independent mass/time contract and is invariant under velocity-storage convention changes. |
| B1-4 | closed | Golden expectations are hand-calculated at nontrivial `a`, `h`, and source units and do not call production conversion helpers. |
| B2-1 | closed | Missing gas `InternalEnergy` and `Density` reject by default; explicit configured values are recorded and tested. |
| B2-2 | closed | Star and BH defaults require typed normalized policy/value selection and are serialized into manifest/provenance. |
| B2-3 | closed | Canonical tracer mapping is rejected during CLI validation; negative stellar birth times receive an AREPO wind diagnostic. |
| B3-1 | closed | Embedded and sidecar manifest SHA-256 values are recomputed and compared; verified status is recorded. |
| B3-2 | closed | HDF5, sidecar, and completion marker use partial files, validation, ordered rename, and rollback; six injected failures are tested. |
| B4-1 | closed | Converter appends source chunks directly to canonical hyperslabs, uses bounded external duplicate-ID runs, and records `FullStateMaterialized=0` plus peak staging counters. |
| B5-1 | environment-blocked | Persistent source/dataset handles and open counters are implemented and HDF5-tested; multi-rank runtime evidence is unavailable. |
| B5-2 | environment-blocked | Source batches multiple chunks within `ic_staging_particle_count`, reports routing batches/collective phases, and MPI syntax-checks; collective reduction requires real np1/np2/np4 execution. |
| B5-3 | partially closed | New session/read/route/decode/append work is inside rank-consistent phases and fault seams remain; real MPI fault execution and further audit of allocation-capable accounting between collectives remain. |
| B6-1 | closed | Manifest mode derives the source set, resolves relative paths from the manifest directory, and rejects conflicting explicit source authority. |
| B6-2 | closed | Reader sessions revalidate size and SHA-256 and retain the verified file handle through payload reads. |
| B7-1 | partially closed | Six durable responsibilities were split with a codec unit test and structural guardrails; the remaining 4,136-line ingestion unit still contains schema, serial, distributed, and audit responsibilities. |
| B8-1 | closed | Config, schema, compatibility, distributed status, capability reporting, repair ledgers, and output docs state the supported behavior and provisional MPI status. |
| B8-2 | closed at packaging | Recursive hygiene and archive listing checks reject build/output/cache/binary/ADS/partial artifacts; final evidence is emitted with the ZIP. |
| B8-3 | closed | No MPI/FFTW/parallel-HDF5 pass is claimed; exact blocked and unexecuted commands are recorded. |

No P0 scientific-correctness item remains open. Campaign B is nevertheless classified as **partially complete** until the real MPI acceptance matrix runs, the automatic pre-production imbalance/rebalance trigger is implemented, and the remaining IC monolith is split further.
