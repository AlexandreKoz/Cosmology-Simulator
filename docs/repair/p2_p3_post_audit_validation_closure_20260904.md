# P2/P3 post-audit validation and runtime-hygiene closure

Status: **PARTIALLY ACHIEVED**
Date: 2026-09-04
Campaign: P2/P3 Post-Adversarial-Audit Validation, Runtime-Hygiene, and Residual-Defect Closure

## Scope and provenance

This campaign starts from the supplied post-P0/P1 authoritative source archive. The archive did not contain upstream `.git` metadata, so the upstream starting commit is unavailable. For exact local diff generation only, the extracted source was committed as synthetic baseline `8ec6883`. That hash is not upstream project provenance and must not be reported as the repository's original commit.

The inherited P0/P1 state is documented in `docs/repair/p0_p1_post_audit_runtime_closure_20260903.md`. This campaign preserves those TreePM/PM consensus, generic DMO ghost-lane, gas parent-identity, TreePM complementarity, distributed-IC, failure-coordination, bounded-communication, HDF5/restart, and memory-governor repairs.

No solver tolerance or scientific acceptance threshold was weakened.

## Reproduction and root-cause summary

### Runtime application smoke artifact contract

The serial and MPI runtime smoke scripts still asserted obsolete snapshot paths. Production output now writes one logical snapshot set under `snapdir_###`, with member files plus exactly one `*.complete` completion manifest. The smoke scripts were corrected to assert the authoritative snapshot-set layout and stable completion-manifest fields, and to require the production operational report to record successful snapshot write completion and scientific readback/validation.

The production smoke remains an FFT-backed TreePM acceptance gate. It is not registered when FFTW is disabled.

### Runtime capability / manifest-v4 authority

The post-P0/P1 source already advertised `manifest v4`, so the original prose mismatch no longer reproduced. A residual duplication remained: schema name/version/converter-version literals were repeated across the manifest record, validator, runtime capability prose, and test.

The IC reader contract now owns canonical constexpr schema/version identifiers. Manifest defaults, validation, runtime capability text, and the capability unit test derive from that authority. No duplicate schema-version authority was introduced.

### ASan allocation-failure test infrastructure

`unit_scratch_memory_governor` used an actual `SIZE_MAX` array allocation to force backing-allocation failure. GCC ASan aborts on that impossible request before the C++ failure path can prove reservation rollback.

The unit-test executable now provides a test-local deterministic one-shot `operator new[]` failure. The allocator still successfully reserves through the `MemoryGovernor` first, the backing allocation then throws `std::bad_alloc`, and the test proves both reserved and committed bytes return to zero. Production allocator/governor semantics are unchanged; sanitizer behavior is not globally suppressed.

### HDF5-without-FFTW production TreePM validation gating

`validation_phase2_mpi_gravity_single_rank` and its MPI rank variants were registered without requiring a production FFT backend. This made the HDF5-only preset present a diagnostic naive-DFT path as if it were production-gravity acceptance.

The complete Phase-2 production gravity validation executable/matrix is now registered only when `COSMOSIM_ENABLE_FFTW=ON`. Its labels explicitly include `fftw`, and the runtime-truth CTest guard asserts absence when FFTW is disabled. The diagnostic naive-DFT backend remains available to focused diagnostic fixtures; it is not promoted to production TreePM validation.

The validation plan now documents the FFTW requirement and the existing eight-rank Phase-2 entry.

### Star-formation MPI and periodic TreePM MPI residuals

No speculative star-formation or TreePM source edits were made. The required multi-rank reruns could not be performed because the environment has no MPI C++ development support, and FFTW development files are also unavailable. These residual audit findings therefore remain runtime-unclosed rather than being guessed away.

## Files and symbols changed

- `CMakeLists.txt`: gate Phase-2 production gravity validation on FFTW and label it truthfully.
- `include/cosmosim/io/ic_reader.hpp`: define the single IC audit-manifest schema/version/converter authority.
- `src/io/ic_manifest.cpp`: validate manifests against the canonical constants.
- `src/workflows/runtime_capabilities.cpp`: derive manifest-version capability detail from the same authority.
- `tests/unit/test_runtime_capabilities.cpp`: assert the canonical manifest version without duplicating literal `4`.
- `tests/unit/test_scratch_memory_governor.cpp`: deterministic test-local backing-allocation failure after reservation.
- `tests/integration/test_runtime_app_smoke.cmake.in`: validate serial snapshot-set output/completion/readback contract.
- `tests/integration/test_runtime_app_mpi_treepm_smoke.cmake.in`: validate MPI snapshot-set output/completion/readback contract.
- `tests/integration/test_runtime_truth_ctest_labels.cmake.in`: assert FFTW-dependent smoke/Phase-2 registration truth, including eight-rank coverage.
- `docs/validation_plan.md`: document FFTW gating and the existing eight-rank Phase-2 matrix.
- `docs/repair_open_issues.md`: supersede the stale flat-snapshot smoke diagnosis with the corrected snapshot-set root cause.

No snapshot/restart schema, runtime physics state, configuration key, numerical force split, or MPI wire protocol was changed.

## Validation evidence

### Focused CPU Debug

Commands used during the campaign included:

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --test-dir build/cpu-only-debug -R 'unit_scratch_memory_governor|unit_runtime_capabilities|integration_runtime_truth_ctest_labels' --output-on-failure
```

Final focused result: **3/3 PASS**.

The broad CPU Debug run was executed in bounded ranges because the command execution window is finite. All registered ordinary CPU tests completed successfully; combined campaign evidence yielded **137/137 PASS** with no unexpected CPU Debug failure.

### HDF5 Debug, FFTW disabled

The HDF5-only preset configured with:

```text
COSMOSIM_ENABLE_HDF5=ON
COSMOSIM_ENABLE_FFTW=OFF
COSMOSIM_ENABLE_MPI=OFF
```

The post-patch CTest inventory no longer registers either `integration_runtime_app_smoke` or `validation_phase2_mpi_gravity_single_rank` in this unsupported production-TreePM configuration. The broad HDF5 campaign sweep reached test 133 with all completed tests passing before the command window ended; therefore a complete HDF5-suite green total is not claimed here.

### ASan + UBSan

The sanitizer preset is GCC Debug with:

```text
-fsanitize=address,undefined -fno-omit-frame-pointer
```

Focused command:

```text
ctest --test-dir build/asan-debug -R '^unit_scratch_memory_governor$' --output-on-failure
```

Final result: **1/1 PASS**. The deterministic failure still proves governor reservation rollback and does not use a global ASan suppression.

### FFTW and MPI dependency-complete matrices

The current environment blocks the required production runtime matrices.

```text
cmake --preset pm-hdf5-fftw-debug
```

fails because the FFTW3 serial double-precision development library is unavailable (`libfftw3-dev` / equivalent not installed).

```text
cmake --preset mpi-hdf5-fftw-debug
```

fails during MPI discovery because `mpi-cxx` / `MPI_CXX` development support is unavailable.

Consequently the following were not runtime-executed on this environment after the patch:

- `integration_runtime_app_smoke` in HDF5+FFTW;
- `integration_runtime_app_mpi_treepm_smoke_two_rank`;
- Phase-2 gravity np1/np2/np3/np4/np8 with production FFTW;
- the registered MPI star-formation family;
- periodic TreePM MPI rank-count variants.

These remain validation blockers, not evidence of a new production defect.

### Release assertion smoke

Not rerun after the final source edits in this environment. No release-specific source path was modified.

### Repository hygiene

After removing generated `build/` and `test_outputs/` trees, the final source tree passed:

```text
bash ./scripts/ci/check_repo_hygiene.sh
```

`git diff --check` also passed.

## Final finding classification

| Finding | Final state | Root cause / evidence |
|---|---|---|
| CD-001 | REPAIRED | Manifest-v4 runtime capability already matched post-P0/P1; remaining duplicated version authority was consolidated into the IC contract constants. |
| CD-006 | TEST-CONTRACT CORRECTED | Runtime smoke expected obsolete flat snapshot output instead of the authoritative `snapdir_###` logical snapshot set and completion manifest. |
| PD-002 | REMAINS OPEN | No independent star-formation source defect demonstrated; required MPI rerun is blocked by unavailable MPI C++ development support. |
| TID-001 | TEST-CONTRACT CORRECTED | Sanitizer-hostile `SIZE_MAX` allocation replaced by deterministic bounded test-local backing-allocation failure while preserving reservation rollback coverage. |
| TID-002 | TEST-CONTRACT CORRECTED | Capability test now derives the manifest-version assertion from the canonical IC schema authority rather than a duplicated literal. |
| TID-003 | TEST-CONTRACT CORRECTED | Runtime-app smoke scripts now validate the production snapshot-set artifact and completion/readback contract. |
| DRT-001 | DOCUMENTATION ALIGNED | Validation/runtime-truth documentation now states the canonical manifest-v4 authority and FFTW production-gate semantics without creating another schema authority. |
| INC-001 | REMAINS OPEN | Periodic TreePM MPI residual cannot be runtime reclassified without MPI+FFTW execution; no speculative TreePM edit made. |
| HDF5-without-FFTW feature-gating issue | FEATURE-GATING CORRECTED | Phase-2 production gravity validation and real runtime smoke are absent when FFTW is unavailable; diagnostic naive DFT is not labeled as production acceptance. |

## Remaining coverage gaps

The audit's broader OpenMP determinism, adversarial configuration boundaries, MPI AMR reflux conservation, install/export consumer validation, large-volume DMO, CUDA parity, rank-count-changing restart rejection, and schema-boundary expansion items remain coverage work unless independently demonstrated as defects. They were intentionally not expanded into this focused campaign.

## Closure verdict

```text
P2/P3 POST-AUDIT VALIDATION CLOSURE: PARTIALLY ACHIEVED
```

The concrete runtime-smoke contract, manifest-version authority, sanitizer allocation-failure infrastructure, and FFTW feature-gating defects are repaired at source/test-contract level. Full closure is withheld because dependency-complete FFTW and MPI runtime acceptance, including the residual star-formation and periodic TreePM MPI families, could not be executed in this environment, and the broad HDF5 suite did not finish its final tail before the campaign execution window ended.
