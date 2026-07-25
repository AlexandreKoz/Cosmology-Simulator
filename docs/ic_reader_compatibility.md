# Initial-condition ingestion, conversion, and compatibility

The external-IC path is an explicit scientific contract shared by typed
configuration, the runtime reader, the compiled converter, and
`chui_ic_audit_manifest` schema version 4. A GADGET/AREPO-looking filename never
supplies implicit units, coordinate frame, velocity convention, Hubble scaling,
scale-factor powers, species mapping, or missing-field defaults.

Authoritative interfaces:

- `include/cosmosim/io/ic_reader.hpp`
- `src/io/ic_manifest.cpp`
- `src/io/ic_conversion_catalog.cpp`
- `src/io/ic_canonical_bundle.cpp`
- `src/io/ic_reader_session.cpp`
- `src/io/ic_reader_file_set.cpp`
- `src/io/ic_byte_codec.cpp`
- `src/io/ic_record_codec.cpp`
- `src/io/ic_sha256.cpp`
- `tools/convert_ic.cpp`
- `src/workflows/initial_condition_runtime.cpp`

## Typed convention and source authority

`mode.ic_convention` is one of `generated`, `chui_canonical_v1`,
`gadget_arepo_bridge_v1`, or `manifest_v1`.

A direct bridge requires explicit source length, mass, and velocity SI-per-code
factors, source coordinate frame, source velocity convention, all relevant
Hubble and scale-factor exponents, and PartType2/3 policies. Zero or unspecified
scientific values are fail-closed sentinels.

Manifest mode requires only `mode.ic_manifest_file`; source files are derived
from the manifest. A separately supplied non-generated `mode.ic_file` is
accepted only when it resolves to the same first source member, otherwise the
configuration is rejected as conflicting provenance. In the converter,
manifest mode also derives species and missing-field policies/values from the
manifest; policy CLI overrides are rejected rather than silently ignored.

## Real AREPO header integer compatibility

`NumPart_ThisFile` and `NumFilesPerSnapshot` accept signed or unsigned HDF5
integer storage up to 64 bits. Reads use a wide intermediate and reject:

- negative signed values;
- overflow beyond the destination range;
- floating-point or other non-integer types;
- malformed rank or extent.

The manifest preserves the actual source scalar class, byte width, signedness,
byte order, rank, and dimensions. `NumPart_Total` and
`NumPart_Total_HighWord` retain strict unsigned low/high-word handling. Mixed
signed and unsigned per-file count attributes are allowed only when their values
and all other file-set contracts agree.

## Deterministic multifile discovery and source identity

`NumFilesPerSnapshot` determines an exact ordered sequence such as
`snapshot.0.hdf5` through `snapshot.(N-1).hdf5`; directory iteration is never
used as authority. Every member must agree on cosmology, box, epoch,
`MassTable`, file count, reconstructed 64-bit totals, populated families, and
field schema. Per-file counts must match dataset extents and sum to declared
totals.

Discovery records path, size, header metadata, and SHA-256. Payload ingestion
then opens a persistent reader session per assigned file, revalidates file size,
rehashes the source, and keeps the HDF5 file and opened dataset handles stable
across multiple chunks. A changed source fails before its payload is accepted.
This guarantees that the validated payload comes from a file with the same size
and SHA-256 as the inspected source while the session handle remains open.

## Audit manifest schema v4

Schema v4 records:

- ordered source files, sizes, SHA-256 values, original headers, and source
  integer storage metadata;
- actual dataset datatype, byte order, rank, dimensions, and selected alias;
- a field-specific conversion contract with dimensional powers, Hubble and
  scale-factor exponents, frame contribution, and velocity-convention power;
- explicit species policies and typed missing-field contracts;
- converted, preserved, defaulted, dropped, or rejected disposition;
- canonical-bundle verification status and verified manifest digest;
- converter version, warnings, diagnostics, and normalized provenance.

Version 3 is rejected because its prior semantics permitted the velocity-storage
convention to leak into non-velocity fields and did not encode the new
missing-field and verified-bundle contracts.

## Field-specific conversion catalog

Runtime ingestion and `cosmosim_convert_ic` use one C++ conversion equation:

```text
target = stored
       * base_unit_to_si
       * h^hubble_exponent
       * a^(scale_factor_exponent + frame_scale_factor_exponent)
       * velocity_convention_multiplier^velocity_convention_power
       / target_si_per_code(L^length_power M^mass_power T^time_power)
```

The catalog is field-specific:

- `Coordinates`: length units, explicit Hubble/scale-factor behavior, and an
  explicit physical-to-comoving frame transform;
- `Velocities`: velocity units plus the selected AREPO snapshot velocity-storage
  convention;
- `Masses`: mass units and configured Hubble/scale-factor exponents only;
- `Density`: independent `M/L^3` dimensional and cosmological exponents;
- `InternalEnergy`: independent `L^2/T^2` units with velocity-convention power
  zero;
- `BH_Mdot`: independent `M/T` units with velocity-convention power zero;
- IDs and indices: exact unsigned identity values with no physical conversion.

Changing the selected velocity convention therefore changes `Velocities` only;
it does not change `InternalEnergy` or `BH_Mdot`. Golden tests use independent
hand-calculated values at nontrivial `a`, `h`, and source units and compare
direct, manifest-driven, canonical, and converter paths against that truth.

## Explicit missing-field policy

Supported policy values are `reject`, `reconstruct`, `use_config_value`,
`dialect_defined_default`, and `preserve_unavailable`. This build implements
`reject`, validated `use_config_value`, and narrowly defined dialect defaults;
it rejects `reconstruct` and `preserve_unavailable` until their scientific
storage/reconstruction paths exist.

Production defaults are fail-closed:

- missing gas `InternalEnergy`: reject;
- missing gas `Density`: reject;
- missing stellar formation time, initial mass, or metallicity: reject;
- missing `BH_Mdot`: reject unless an explicit recorded default/value is chosen.

Every non-reject resolution is recorded in normalized configuration, manifest,
conversion audit, runtime provenance, and diagnostics. Numerical floors are not
used to disguise absent source physics. Negative stellar formation times are
rejected with an AREPO wind-particle diagnostic. The standalone canonical
converter rejects tracer mapping during CLI validation because canonical schema
v2 has no tracer output family. Direct converter mode exposes the same typed
policy/value pairs through `--gas-internal-energy-policy`,
`--gas-density-policy`, `--star-formation-time-policy`,
`--star-initial-mass-policy`, `--star-metallicity-policy`, and
`--bh-mdot-policy`, together with their matching `--*-value` arguments for
`use_config_value`; it does not have a separate implicit-default path.

## Canonical bundle verification and finalization

Canonical HDF5 uses header schema name `chui_canonical_v1`, schema version 2.
The normalized schema-v4 manifest JSON is embedded in
`/Provenance/ConversionManifestJson`, written as a sidecar, and bound by
`ConversionManifestSha256`. Runtime import:

1. reads the embedded manifest;
2. recomputes its SHA-256;
3. reads and hashes the named sidecar;
4. requires embedded and sidecar JSON to agree exactly;
5. validates the completion marker binding canonical filename, sidecar filename,
   and digest;
6. records `manifest_verified=true` and the verified digest.

A syntactically valid 64-character string is not sufficient. Absence, mismatch,
or tampering fails closed.

The converter writes HDF5, manifest, and completion-marker `.part` files,
validates their binding, and finalizes the three-member bundle with rollback of
partial or already-renamed members on any failure. Fault-injection tests cover
HDF5 write, manifest write, digest mismatch, and every rename/finalization seam.

## Genuinely streaming canonical conversion

`cosmosim_convert_ic` follows:

```text
source chunk -> inspect/validate -> convert -> append canonical hyperslab
             -> bounded ID/count/statistics audit
```

It pre-creates final datasets from validated manifest counts, keeps output
handles open, and appends each converted chunk directly. It does not construct a
complete converted `SimulationState`. Exact duplicate-ID validation uses bounded
sorted runs and a fixed-fan-in external merge. Header counters record chunk
capacity, peak batch records, reader staging bytes, flow identity, and
`ChuiConverterFullStateMaterialized=0`.

## Distributed runtime status

The distributed reader keeps per-file HDF5 sessions and dataset handles open
across chunks. Several deterministic chunks are accumulated into bounded
per-peer buffers up to `mode.ic_staging_particle_count`, then exchanged once per
routing batch. Exact route, source/final identity, duplicate, count, mass, and
ownership audits remain enabled. Initial x-slab ownership is ingestion
ownership; it is not a claim of final work balance.

The source and MPI test matrix are present, but this capability remains
`provisional` because the Campaign B completion environment had no `MPI_CXX`
wrapper and could not execute the one-, two-, and four-rank acceptance matrix.
No MPI, FFTW, FFTW-MPI, or parallel-HDF5 pass is claimed.
