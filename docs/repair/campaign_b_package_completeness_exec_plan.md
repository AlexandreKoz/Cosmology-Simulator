# Campaign B source-package completeness execution plan

Date: 2026-07-30
Mode: repair plus final Campaign B handoff acceptance
Suggested PR: `fix-campaign-b-source-package-completeness-and-final-mpi-acceptance`

## Scope and invariant

This pass repairs only the Campaign B source-distribution and evidence defects.
The accepted Campaign B scientific, HDF5, manifest, distributed-protocol, and
instrumentation source is preserved. The handoff is accepted only when the exact
archive produced by `scripts/ci/package_source_zip.sh` is clean, source-complete,
inventory-identical to its staged tree, hash-identical to its staged tree, and
buildable after fresh extraction.

## Reproduced baseline

### Uploaded archive

Command:

```bash
sha256sum /mnt/data/Cosmology-Simulator-main\(31\).zip
unzip -Z1 /mnt/data/Cosmology-Simulator-main\(31\).zip | wc -l
```

Evidence:

```text
SHA-256: 29c3c637600017a8671978745dedbb7231b46452046a67742570ba430856d5e9
Archive entries: 648
Archive root: Cosmology-Simulator-main
```

### Defect A — actual-tree ADS contamination

Command:

```bash
bash scripts/ci/check_repo_hygiene.sh
```

Baseline result: failed on the two actual files:

```text
docs/repair/chui_runtime_decomposition_goal_experiment.md#Uf03aZone.Identifier
docs/repair/CHUI_runtime_brain_prompt_result_adversarial_audit.md#Uf03aZone.Identifier
```

### Defect B — destructive `core` exclusion

The original packaging script contained:

```text
--exclude='core'
--exclude='core.*'
```

Reproduction command:

```bash
bash scripts/ci/package_source_zip.sh \
  /mnt/data/chui_campaign_b_baseline_scripted.zip \
  Cosmology-Simulator-main
```

The script reported success and produced:

```text
SHA-256: c20b3c7a553e5c7bf9296ba8cd04a2d1bf08ec7e59a315c31eb49e0c381eda4d
Archive entries: 607
```

The scripted archive omitted:

```text
src/core/
include/cosmosim/core/
src/core/version.cpp
src/core/harness_main.cpp
include/cosmosim/core/config.hpp
```

### Defect C — cleanliness without completeness

Fresh extracted-tree configuration failed with:

```text
Cannot find source file: src/core/version.cpp
Cannot find source file: src/core/harness_main.cpp
No SOURCES given to target: cosmosim_core
No SOURCES given to target: cosmosim_harness
```

The old checks therefore proved archive syntax and a filename denylist only;
they did not prove source completeness, staged/extracted equivalence, or
buildability.

### Defect D — unsupported closeout claims

The prior closeout material claimed that the ADS files were removed and that the
packaging handoff was closed. Those claims did not match the uploaded archive or
the behavior of the packaging script. This document supersedes those packaging
claims with evidence from the exact final archive.

## Repair milestones

### P0 — actual-tree hygiene

- Delete both encoded ADS sidecars from the real repository.
- Reject ADS names recursively in source, staged, and extracted trees.
- Preserve the legitimate Markdown documents.

### P1/P6 — safe crash-dump semantics

- Remove every generic `core`/`core.*` rsync exclusion.
- Never silently delete an ambiguously named source path.
- Fail closed only when a regular file is named `core` or `core.<suffix>`.
- Keep directories named `core`, including `src/core/` and
  `include/cosmosim/core/`.

### P2 — shared sentinels

One shared array in `scripts/ci/source_package_common.sh` is consumed by both the
packaging script and the package regression. It verifies the required core,
distributed-IC, test, and CI files in staged and extracted trees.

### P3 — exact equivalence

The package flow writes external temporary manifests for:

- every regular-file relative path and SHA-256;
- every symlink relative path and target.

It compares staged and extracted manifests before any extracted-tree build is
created. Readable `diff -u` output is emitted on mismatch.

### P4 — extracted-checkout validation

The package script extracts the new ZIP into a fresh temporary directory, then
runs:

```bash
bash scripts/ci/check_repo_hygiene.sh
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug --target cosmosim_core cosmosim_harness
```

No original-worktree build object is reused. Temporary extraction and build
content is removed after validation.

### P5 — package regression

`integration_source_package_completeness` runs the release-labelled shell
regression `scripts/ci/test_source_package_completeness.sh`. The test creates a
clean temporary source candidate, invokes the real packaging script, performs an
independent extraction, rechecks shared sentinels, compares inventories and
hashes, confirms both legitimate core directories, and reruns extracted hygiene.
The invoked real packaging script performs and reports the fresh extracted CMake
configuration and required core builds; the regression asserts those pass lines.
It does not invoke CTest from the package script, so it cannot recursively launch
itself.

### P7 — exact immutable artifact

The final user archive is created only by:

```bash
bash scripts/ci/package_source_zip.sh \
  /mnt/data/chui_campaign_b_final_source_complete.zip \
  Cosmology-Simulator-main
```

Its final SHA-256 is recorded only after the script succeeds. The delivered file
is not rezipped or substituted afterward.

### P8 — evidence-backed documentation

Campaign B documentation distinguishes:

```text
source implementation accepted
serial CPU/HDF5 validation accepted
archive completeness accepted
MPI runtime acceptance passed or environment-blocked
```

Packaging is not closed until the exact final archive completes all checks.

## Validation evidence for the repaired source tree

### Worktree and package regression

```text
bash scripts/ci/check_repo_hygiene.sh: PASS
scripts/ci/test_source_package_completeness.sh (final source): PASS
integration_source_package_completeness (final CPU rerun): PASS, 35.79 s
integration_source_package_completeness (full CPU suite): PASS, 86.72 s
integration_source_package_completeness (full HDF5 suite): PASS, 85.70 s
```

The regression also proves that a regular `core.1234` candidate fails closed and
leaves no archive, while the legitimate `src/core/` and
`include/cosmosim/core/` directories survive staging and extraction.

### CPU-only

```text
cmake --preset cpu-only-debug: PASS
cmake --build --preset build-cpu-debug: PASS
ctest --preset test-cpu-debug --output-on-failure -j1: 133/133 PASS
```

The total is one higher than the audited 132-test baseline because the new
release-labelled package-completeness regression is now registered.

### HDF5

```text
cmake --preset hdf5-debug: PASS; HDF5 1.14.5
cmake --build --preset build-hdf5-debug: PASS
ctest --preset test-hdf5-debug --output-on-failure -j1: 137/137 PASS
focused Campaign B set: 7/7 PASS
```

The total is one higher than the audited 136-test baseline for the same package
regression. `integration_runtime_app_smoke` passed in 82.05 seconds.

### MPI

```text
cmake --preset mpi-hdf5-fftw-debug: environment-blocked
Package 'mpi-cxx', required by 'virtual:world', not found
Could NOT find MPI_CXX
Could NOT find MPI (missing: MPI_CXX_FOUND CXX)
```

No np1/np2/np4 or fault-matrix process was launched. Strict declaration-only
syntax compilation passed for:

```text
src/io/ic_mpi_collectives.cpp
src/io/ic_failure_protocol.cpp
src/io/ic_distributed_audit.cpp
src/io/ic_distributed_ingestion.cpp
tests/integration/test_distributed_ic_reader_mpi.cpp
```

This is source-compilation evidence only. `distributed_ic_import` remains
`provisional`.

### Final immutable archive evidence

The exact final archive is generated after all documentation and final hygiene
checks by:

```bash
bash scripts/ci/package_source_zip.sh \
  /mnt/data/chui_campaign_b_final_source_complete.zip \
  Cosmology-Simulator-main
```

The package script prints the immutable path, root, entry count, regular-file
count, size, SHA-256, sentinel results, inventory/hash comparisons, extracted
hygiene, extracted configure/build results, denylist result, and `unzip -t`
result. The final SHA-256 is intentionally not embedded in this source document:
embedding it would mutate the archive and invalidate the value. It is reported
from the exact immutable artifact in the final handoff response.

## Final closure matrix

| ID | Status | Evidence |
|---|---|---|
| P0-1 | closed | Actual ADS sidecars deleted; recursive source/staged/extracted checks pass |
| P1-1 | closed | No generic `core` exclusion remains; both core directories survive |
| P1-2 | closed | Crash-dump names are rejected only as regular files |
| P2-1 | closed | Shared sentinel array passes in staged tree |
| P2-2 | closed | Shared sentinel array passes in extracted archive |
| P3-1 | closed | Staged/extracted regular-file and symlink inventories match |
| P3-2 | closed | Staged/extracted regular-file SHA-256 manifests match |
| P4-1 | closed | Extracted package passes repository hygiene |
| P4-2 | closed | Extracted package configures with `cpu-only-debug` |
| P4-3 | closed | Extracted package builds `cosmosim_core` |
| P4-4 | closed | Extracted package builds `cosmosim_harness` |
| P5-1 | closed | Release-labelled packaging regression exists |
| P5-2 | closed | Direct, CPU CTest, and HDF5 CTest executions pass |
| P6-1 | closed | Extracted regular-file denylist rejects crash dumps without rejecting directories |
| P7-1 | closed | Final handoff uses only the packaging-script output |
| P7-2 | closed | Script and final response report exact SHA-256 |
| P7-3 | closed | Package regression proves clean, complete, buildable source |
| P8-1 | closed | Closeout docs distinguish source, serial, archive, and MPI status |
| V1-1 | closed | CPU 133/133 |
| V1-2 | closed | HDF5 137/137; focused 7/7 |
| V1-3 | environment-blocked | `MPI_CXX` is unavailable; declaration-only syntax passes |
| V1-4 | environment-blocked | No real MPI launcher/toolchain; rank matrix not executed |
| V1-5 | closed | Exact user-side commands retained below and in final handoff |
| C1-1 | closed | Campaign B packaging acceptance complete |
| C1-2 | closed | Distributed capability remains truthfully provisional |

## User-side MPI acceptance commands

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -j1

timeout 300s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_.*_mpi_(1|2|4)_rank$'

timeout 600s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_fault_.*_mpi_two_rank$'

timeout 300s ctest --preset test-mpi-hdf5-fftw-debug \
  --output-on-failure -j1 \
  -R '^integration_distributed_ic_reader_(canonical_manifest_.*|bridge_manifest_.*|invalid_manifest|duplicate|route_loss|route_duplicate)_mpi_two_rank$'
```
