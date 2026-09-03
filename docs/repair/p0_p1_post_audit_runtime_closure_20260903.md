# P0/P1 post-audit runtime closure

Status: source repairs implemented; dependency-complete distributed runtime acceptance blocked by environment.
Date: 2026-09-03
Campaign: P0/P1 Post-Adversarial-Audit Runtime Closure

## Scope and provenance

This closure addresses the concrete P0/P1 failures exposed by the broader
post-campaign build/test audit. The supplied authoritative source ZIP did not
contain `.git`, so no upstream commit hash can be recovered from the artifact.
For exact local diff generation only, the campaign created a synthetic baseline
commit `39d868e` from the untouched extracted ZIP. That hash is not upstream
project provenance and must not be reported as the repository's original
commit.

The repair preserves the existing runtime-truth, memory-governor,
FailureCoordinator, distributed-IC routing/reconciliation, TreePM, HDF5, and
bounded-communication authorities. No numerical tolerance was loosened to make
a runtime defect pass.

## Root causes and repairs

### Distributed TreePM / PM rank consensus

Two distinct contracts had been conflated.

- `particle_index_generation` and `cell_index_generation` are rank-local dense
  indexing/view epochs. Migration or local reorder can legitimately advance
  them differently on different ranks. They no longer implicitly advance the
  communicator-global `gravity_source_generation`.
- Physical gravity-source mutation continues to advance the explicit global
  gravity-source generation at the authoritative source-update boundaries.
- PM collective-entry fingerprints no longer require target coordinate/output
  layout enums to be communicator-identical. Those layouts describe rank-local
  target representations; communicator-global PM solve semantics remain in the
  fingerprint and remain fail-closed.

This keeps real collective semantic disagreement detectable while removing the
false invariant that produced the observed rank-local generation/layout
failures.

### Generic DMO ghost payload optional-lane contamination

The generic ghost wire/commit path erased optional-lane absence by resizing
all destination lanes, which could make a gravity-only DMO payload advertise
hydro state after commit.

The wire now carries an explicit optional-lane presence mask. Unpack and commit
preserve absence, validate incompatible non-empty peer schemas, and materialize
only lanes that are actually present. Empty payloads no longer classify empty
hydro vectors as hydro data. Hydro-specific state remains owned by the
`gas_cell_id`-keyed hydro exchange.

### Gas optional parent identity / HDF5 invariant

`GasCellIdentityMap` previously allowed the contradictory in-memory state
`parent_particle_id = optional{0}`, while snapshot/restart readers correctly
rejected a present zero parent ID. The authoritative assignment/rebuild
boundary now rejects present-zero parent identity. Parentless cells use an
absent optional; valid parents remain nonzero values. Star-formation and restart
fixtures were aligned with that canonical invariant.

### Periodic TreePM complementarity contract

The single-rank periodic integration test compared the TreePM result against a
minimum-image direct reference and additionally required the split error to be
no worse than PM-only. That monotonic comparison is not a valid periodic image
accuracy contract because the PM field contains periodic-image forces that the
minimum-image reference omits.

The invalid incidental comparison was removed without changing force
parameters or broad error gates. The independent periodic Ewald validation
remains the authoritative periodic accuracy/complementarity evidence and passes
in the executable matrix.

### Distributed IC acceptance closure

The distributed reader architecture was preserved. Three narrow acceptance
issues were repaired:

- routing collective-phase telemetry is recorded on every participating rank
  because it describes communicator-global protocol work, rather than being
  populated on root only;
- malformed-manifest fixtures now first match frozen runtime cosmology,
  `BoxSize`, and start epoch so each fixture reaches the corruption it claims
  to test;
- the single-rank distributed entry point reports `already_partitioned=true`
  consistently with its API boundary: the returned state is final rank-local
  ownership even when implementation delegates to the serial reader.

No duplicate-ID, route-integrity, failure-coordination, bounded-staging,
source/final reconciliation, or provenance check was weakened.

### Secondary HDF5 acceptance fixture/metadata closure

The broad HDF5 sweep exposed three stale acceptance details outside the
original defect list. They were repaired without weakening production guards:

- runtime capability detail now advertises the already-authoritative manifest-v4
  strict audit and source-chunk-to-canonical streaming path;
- `convert_ic`'s internal typed-unit configuration fixes the fixture box size in
  physical Mpc so selecting an output unit cannot accidentally change a
  dimensionless default and violate the certified TreePM softening profile;
- the HDF5 restart-contract gravity fixture explicitly uses the diagnostic naive
  DFT only when no fast FFT backend exists and reduces that diagnostic mesh to
  4^3. With a fast FFT backend it retains the production 16^3 profile. The
  production fail-closed backend policy is unchanged.

## Validation

### CPU Debug

Build command:

```text
cmake --build --preset build-cpu-debug -j4
```

The functional CPU matrix passed all 137 ordinary tests. The suite was
executed in bounded CTest ranges because the campaign execution harness limits
individual command windows. After the final HDF5-related source adjustments,
the affected CPU `unit_runtime_capabilities` and
`validation_phase2_mpi_gravity_single_rank` tests were rebuilt and passed.

The separately executed final-tree `integration_source_package_completeness`
script also passed, yielding combined evidence for all 138 registered CPU Debug
tests.

### HDF5 Debug

The HDF5 build completed and all 153 ordinary tests passed. The separately
executed final-tree source-package completeness script also passed, yielding
combined evidence for all 154 registered HDF5 Debug tests. This includes the
repaired snapshot/restart gas identity paths, star-formation/AMR roundtrips,
runtime capabilities, `convert_ic`, periodic TreePM integration, and phase-2
gravity single-rank fixture.

### FFTW and MPI + HDF5 + FFTW

The dependency-enabled presets could not be configured in this campaign
environment:

```text
cmake --preset pm-hdf5-fftw-debug
```

is blocked because the FFTW3 development package is unavailable, and

```text
cmake --preset mpi-hdf5-fftw-debug
```

is blocked because `mpi-cxx` / `MPI_CXX` development support is unavailable.

Therefore the required np2/np3/np4/np8 TreePM, DMO, PM, and distributed-IC
runtime matrix was not executed. Source-level fixes must not be promoted to
multi-rank runtime acceptance without that evidence.

### ASan + UBSan

The `asan-debug` configuration and complete build succeeded. The executable
suite showed no production memory-corruption finding. All observed runnable
coverage outside the source-package test and the deliberate maximum-allocation
case passed, including star-formation/AMR, validation integration,
phase-2 gravity single-rank, periodic Ewald, and DMO empty-rank single-rank.

`unit_scratch_memory_governor` deliberately requests `SIZE_MAX` bytes to test
allocation rejection. Under this GCC ASan runtime, the sanitizer aborts on the
impossible allocation before C++ can surface `bad_alloc`/`length_error`; this is
the previously identified test-infrastructure behavior, not a production
memory-corruption finding. The ordinary CPU/HDF5 versions of that test pass.

### Scientific/runtime checks executed

- periodic Ewald reference validation: PASS;
- periodic TreePM integration: PASS in executable backend configuration;
- phase-2 gravity single-rank: PASS;
- gas parentless/valid-parent identity unit and restart/snapshot paths: PASS;
- star-formation + AMR identity/roundtrip paths: PASS;
- generic ghost optional-lane presence/commit regression: PASS;
- DMO empty-rank smoke single-rank: PASS;
- distributed IC multi-rank scientific/accounting checks: environment-blocked.

### Source-package and hygiene

The final source tree passed `scripts/ci/check_repo_hygiene.sh` after generated
`build/` and `test_outputs/` validation artifacts were moved outside the source
tree. The final-tree source-package completeness script passed its staged and
extracted hygiene checks, staged/extracted inventory and hash identity checks,
CMake configure, extracted `cosmosim_core`/`cosmosim_harness` build, archive
denylist, and ZIP integrity test. No validation build/runtime artifact is part
of the campaign source patch.

## Acceptance verdict

`P0/P1 POST-AUDIT RUNTIME CLOSURE: PARTIALLY ACHIEVED`

The confirmed source defects have focused repairs and executable single-rank /
HDF5 / sanitizer evidence. Full closure cannot be declared because the
mandatory FFTW and MPI multi-rank matrices are blocked by missing development
dependencies in the campaign environment. In particular, communicator-wide
TreePM/PM consensus, DMO np2/3/4/8, and distributed-IC fault/acceptance behavior
must still be executed on an MPI+FFTW-capable system.

Distributed IC remains provisional at exactly the repository-authoritative
boundary; this campaign does not promote it from source evidence alone.

## Final package/hygiene result

PASS. Final repository hygiene and source-package completeness both pass on the
campaign tree. The patch bundle is generated separately from the exact final
diff and contains only campaign-modified/new repository files plus the bundle
manifest and unified patch artifact.
