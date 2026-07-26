# Campaign B acceptance repair report

> **Historical report — superseded 2026-07-26.** This document records the
> earlier Campaign B acceptance pass and must not be used as the current closure
> verdict. The authoritative final repair record is
> `docs/repair/campaign_b_final_acceptance_exec_plan.md`. Real MPI np1/np2/np4
> runtime acceptance remains environment-blocked until executed on a
> dependency-complete installation.

Date: 2026-07-24

## Verdict

The Campaign B source-level acceptance repair is complete for the inspected
repository snapshot. CPU-only, HDF5-focused, converter, configuration, reader,
and sanitizer evidence is available. Distributed runtime capability remains
**provisional** because this environment cannot configure an MPI C++ toolchain;
the complete MPI acceptance matrix is implemented and registered but was not
executed here.

The complete HDF5 preset is not reported as fully green because the unrelated
pre-existing `integration_runtime_app_smoke` generated-IC harness process did
not return within repeated bounded command windows. All Campaign B HDF5 tests
passed.

## Audit finding closure

| Audit finding | Repair | Evidence status |
|---|---|---|
| P0-1 direct bridge guessed units/frame/velocity | Added complete typed direct-bridge SI scales, frame, velocity convention, and h/a exponents; incomplete bridges fail before HDF5 access | CPU/HDF5 configuration tests pass |
| P0-2 converter/runtime conversion divergence | Added one dimension-aware C++ conversion engine; compiled converter links the same library; Python tools are launchers only | Gas counterexample, BH, star, direct/manifest/canonical equivalence tests pass |
| P0-3 rank-local exceptions could strand peers | Added phase guards and separated local preparation from later collectives; registered one-rank phase fault tests | Strict MPI-conditioned syntax compile passes; runtime execution pending |
| P0-4 unsafe header attribute reads | Attribute datatype and dataspace are validated before destination allocation/read | Malformed 5/7-element and scalar/vector tests pass under ASan/UBSan |
| P1 semantic datatype gaps | IDs require supported unsigned integer storage; physical fields require supported floating types; lossy classes fail | Floating/signed ID and wrong-type tests pass |
| P1 inaccurate counters | Split metadata/hash/payload/converted/serialized/sent/received/manifest counters and capacity-based transient-memory high-water tracking | Serial/HDF5 tests pass; distributed runtime evidence pending |
| P1 root serial full-file hash | File schema inspection and SHA-256 hashing are assigned deterministically across ranks; rank zero assembles bounded metadata only | Static protocol review and strict syntax compile pass; runtime pending |
| P1 weak completeness proof | Added exact file/chunk coverage, per-chunk source-to-final ID balance, all-axis checks, and bounded external-merge global duplicate audit | Serial reconciliation tests pass; route loss/duplicate MPI tests registered |
| P1 capability contradiction | Production and tests now report distributed import provisional only for MPI+HDF5 builds and unsupported otherwise | CPU/HDF5 capability tests pass |
| P1 narrow MPI evidence | Expanded matrix to direct, canonical, manifest, alternate conversion, type policies, mixed species, scaling, route mutation, duplicate IDs, and phase faults on one/two/four ranks | Registered; strict syntax compile passes; runtime pending |
| P2 unknown fields/families | Enumerate particle groups/datasets; every known-group dataset receives converted/dropped disposition; unknown populated families fail closed | HDF5 malformed tests pass |
| P2 documentation overclaims | Updated configuration, schema, profiling, workflow, validation, and architecture documents to actual behavior and provisional status | Documentation scaffold/release tests pass |
| P2 packaging hygiene | Removed ADS metadata/generated output/build artifacts and restored executable script modes; final archive is verified after extraction | Final packaging verification recorded separately |

## Scientific conversion contract

Manifest schema v3 records explicit physical dimensions and the complete source
scientific convention. The shared engine covers:

- coordinates: length;
- velocities: length/time plus explicit velocity convention;
- particle, initial stellar, and black-hole mass: mass;
- density: mass/length cubed;
- internal energy: length squared/time squared;
- smoothing length: length;
- black-hole accretion rate: mass/time;
- tracer mass history: mass;
- dimensionless fractions and metallicities;
- identity/index fields without numerical unit conversion.

The permanent audit-counterexample regression uses `a=0.5`, `h=0.5`, physical
coordinates, `sqrt_a_scaled_peculiar` velocities, length/mass h exponents `-1`,
input internal energy `4`, and input density `2`. Direct bridge,
manifest-driven ingestion, and canonical conversion agree on internal energy
`8` and density `0.0625` in the configured target code units.

## Distributed collective review

The phase-by-phase collective-order and failure-vote table is maintained in
`docs/architecture/distributed_ic_ingestion.md`. In summary:

- every local source inspection and hash phase votes before metadata exchange;
- conversion and owner/serialization vote before send-layout collectives;
- count/displacement validation and allocations vote before payload exchange;
- wire decode, sidecar append, finalization, source/final reconciliation,
  duplicate audit, and totals validation vote before the next collective;
- empty ranks and zero-count peers remain in identical collective order.

The wire format is field-by-field, fixed-width, and little-endian; no padded C++
object representation is transmitted.

## Validation evidence

### CPU-only

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

- Configure: passed.
- Build: passed.
- Tests: 130/130 passed. The command window ended while the heavy
  `validation_tree_pm_ewald_accuracy` test was running; that test and the final
  test were then run separately and passed without changing tolerances.

### HDF5

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

- Configure: passed.
- Build: passed.
- Campaign B reader, config, converter, manifest, malformed-input, docs, and
  capability tests: passed.
- Overall registered HDF5 suite: 132/133 passed or completed successfully.
- `integration_runtime_app_smoke`: did not return within repeated bounded
  command windows. It uses generated ICs and is outside the external Campaign B
  ingestion path; its timeout was not raised and assertions were not weakened.

### Sanitizers

An external untracked debug build used:

```text
-fsanitize=address,undefined -fno-omit-frame-pointer
```

The following passed with AddressSanitizer and UndefinedBehaviorSanitizer:

- `unit_config_parser`;
- `unit_ic_reader`, including malformed attributes/datatypes and tracer remap;
- `integration_convert_ic` using the compiled converter.

### MPI/FFTW-MPI

One configuration attempt was made:

```bash
cmake --preset mpi-hdf5-fftw-debug
```

It was blocked by the environment:

```text
Could NOT find MPI_CXX
Could NOT find MPI (missing: MPI_CXX_FOUND CXX)
```

MPI-conditioned production and test sources were syntax-compiled with strict
warnings-as-errors against a controlled declaration-only MPI header:

```text
-std=c++20 -Wall -Wextra -Wpedantic -Werror -fsyntax-only
```

This is compile-structure evidence only, not MPI runtime acceptance.

## MPI acceptance handoff

Run on the dependency-complete user system:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure
ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R 'integration_distributed_ic_reader'
```

The focused regex includes:

- normal mixed-species ingestion at one, two, and four ranks;
- canonical, manifest, alternate-convention, and type-policy modes at one, two,
  and four ranks;
- scaling modes at two and four ranks;
- duplicate-ID, route-loss, and route-duplicate tests;
- owner/serialization, send-layout, payload-validation, deserialization,
  sidecar-append, final-state, and source/final-reconciliation fault tests.

Capability must remain provisional until this matrix passes without timeout.

## Remaining boundary

No known source-level requirement from the Campaign B acceptance-repair prompt is
intentionally deferred. The remaining Campaign B boundary is genuine MPI
runtime evidence. The unrelated generated-IC runtime smoke timeout remains a
repository-wide validation issue and is not presented as a Campaign B success.
