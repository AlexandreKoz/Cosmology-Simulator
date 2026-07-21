# Initial-condition ingestion, conversion, and compatibility

The external-IC path is an explicit scientific contract shared by typed
configuration, the runtime reader, the standalone converter, and the audit
manifest. The runtime never infers unit, frame, velocity, h-scaling, or species
semantics merely from a GADGET/AREPO-looking filename.

Authoritative interfaces:

- `include/cosmosim/io/ic_reader.hpp`
- `src/io/ic_reader_file_set.cpp`
- `src/io/ic_manifest.cpp`
- `tools/convert_ic.py`
- `src/workflows/initial_condition_runtime.cpp`

## Typed convention selection

`mode.ic_convention` must be one of:

- `generated`
- `chui_canonical_v1`
- `gadget_arepo_bridge_v1`
- `manifest_v1`

An external `mode.ic_file` without an explicit convention fails configuration
validation. `manifest_v1` additionally requires `mode.ic_manifest_file`.
PartType2 and PartType3 each have an explicit typed policy (`reject`,
`dark_matter`, `star`, `black_hole`, or `tracer`).

## Source container and deterministic file-set discovery

The bridge accepts HDF5 files containing `/Header` and populated
`/PartType0` through `/PartType5` groups. Required header attributes are:

- `NumPart_ThisFile`
- `NumPart_Total`
- optional `NumPart_Total_HighWord`
- `NumFilesPerSnapshot`
- `MassTable`
- `BoxSize`
- `Time`, `Redshift`, `Omega0`, `OmegaLambda`, and `HubbleParam`

For a multifile set, a member such as `snapshot.0.hdf5` determines the exact
ordered sequence `snapshot.0.hdf5` ... `snapshot.(N-1).hdf5`. Discovery does
not depend on directory iteration order. Missing members, duplicate canonical
paths, conflicting `NumFilesPerSnapshot`, inconsistent cosmology/box/epoch,
inconsistent mass tables, conflicting declared totals, or a sum of per-file
counts that differs from the reconstructed 64-bit total are hard failures.
Low/high count words are reconstructed with checked 64-bit arithmetic before
authoritative state allocation.

## Actual HDF5 schema inspection

Every selected dataset is inspected before authoritative state allocation. The
manifest records the real:

- HDF5 class (integer or floating point),
- byte width,
- signedness where applicable,
- byte order,
- rank,
- dimensions,
- selected source path and alias.

A source `float32` field remains recorded as `float32`; conversion into a C++
`double` does not rewrite source provenance. Scalar attributes retain rank zero
and an empty dimension vector. Across a multifile family, path/alias, class,
width, signedness, byte order, rank, and non-record dimensions must agree.

Accepted aliases include:

- Coordinates: `Coordinates`, `Position`, `POS`
- Velocities: `Velocities`, `Velocity`, `VEL`
- Particle IDs: `ParticleIDs`, `ParticleID`, `ID`
- Masses: `Masses`, `Mass`
- Gas internal energy: `InternalEnergy`, `U`, `Internal_Energy`
- Gas density: `Density`, `Rho`
- Gas metallicity: `Metallicity`, `GFM_Metallicity`
- Gas smoothing length: `SmoothingLength`, `Hsml`, `Smoothing_Length`
- Star formation time: `StellarFormationTime`, `GFM_StellarFormationTime`, `BirthTime`
- Star metallicity: `Metallicity`, `GFM_Metallicity`
- Black-hole mass/accretion: `BH_Mass`, `BlackHoleMass`; `BH_Mdot`, `BlackHoleAccretionRate`

The selected alias is part of the audit manifest.

## Strict audit manifest v2

`IcManifest` uses schema name `chui_ic_audit_manifest` and strict schema version
`2`. Serialization and deserialization are round-tripped and validated; an
unsupported or incomplete version is not silently upgraded.

The manifest records:

- ordered source paths, byte sizes, SHA-256 hashes, and provenance IDs;
- original relevant header values for every source member;
- per-file counts, low/high words, and reconstructed 64-bit totals;
- actual dataset schema and aliases;
- source/target units, h and scale-factor exponents, coordinate frame, velocity
  convention, and explicit conversion equations;
- the policy for every particle family;
- converted, defaulted, dropped, rejected, preserved auxiliary, and warning lists;
- converter/tool provenance.

The reference workflow writes the validated manifest to
`<run_directory>/ic_manifest.json` on rank zero. Distributed import compares a
canonical SHA-256 digest of the manifest on every rank before materialization.

## Canonical CHUÍ conversion tool

`tools/convert_ic.py` is the standalone streaming converter. Typical use:

```bash
python3 tools/convert_ic.py \
  --input path/to/snapshot.0.hdf5 \
  --output path/to/chui_ic.hdf5 \
  --manifest path/to/chui_ic.ic_manifest.json \
  --source-convention gadget_arepo_bridge_v1 \
  --coordinate-frame comoving \
  --velocity-convention physical_peculiar
```

A previously validated audit manifest can be used as the authoritative source
contract instead of repeating the convention and conversion flags:

```bash
python3 tools/convert_ic.py \
  --source-manifest path/to/source.ic_manifest.json \
  --output path/to/chui_ic.hdf5 \
  --manifest path/to/chui_ic.ic_manifest.json
```

In manifest mode the converter resolves the ordered source list relative to the
manifest, re-hashes every member, checks byte sizes, counts, species policies,
and the complete real HDF5 datatype/dimension signature before conversion. A
stale or edited source set fails closed.

Use `--help` for the complete set of unit, h-exponent, scale-factor-exponent,
velocity, chunk, and species-policy controls. The converter:

- validates the complete file set before output;
- hashes source members with streaming SHA-256;
- streams source datasets in bounded chunks;
- performs exact duplicate-ID detection over the full unsigned 64-bit domain;
- writes canonical HDF5 names and CHUÍ v1 schema/unit/frame attributes;
- writes the complete audit manifest;
- embeds the final manifest digest in the canonical HDF5 header;
- finalizes both artifacts atomically through temporary files;
- returns nonzero on ambiguity, malformed schema, inconsistent files, duplicate
  IDs, or unsupported populated species.

`tools/convert_ic_manifest.py` is retained only as a compatibility entry point
to the same converter implementation.

## Canonical conversion equations

The source contract stores a base conversion and explicit powers of h and scale
factor. Runtime canonical positions are comoving and runtime velocities are
physical peculiar velocities. Supported velocity conventions include physical
peculiar velocity, sqrt(a)-scaled peculiar velocity, and comoving coordinate
rate. The converter and runtime reader share the same declared equations and
manifest semantics so scientific conversion does not diverge between tools.

Density and specific internal energy are converted through their physical
SI dimensions rather than copied as unlabelled scalars.

## Species and sidecars

- PartType0 maps to gas particles and materializes matching gas-cell and gas
  identity state. Internal energy and density are imported; missing optional
  values use explicit zero defaults recorded in the audit.
- PartType1 maps to dark matter.
- PartType4 maps to stars and materializes one star sidecar per imported star,
  including formation time, birth mass, and metallicity where present or an
  explicitly recorded default where permitted.
- PartType5 maps to black holes and materializes black-hole mass and accretion
  sidecars. Required black-hole mass is validated; optional accretion rate may
  default to zero with an audit record.
- PartType2 and PartType3 obey their separately configured policies. Mapping to
  star, black hole, or tracer is accepted only when all required sidecar fields
  can be constructed.
- Tracer mapping requires the complete versioned tracer field set; an incomplete
  tracer family fails closed.
- Populated unsupported or unknown families are never silently discarded.

Duplicate IDs, nonfinite required values, nonpositive masses, invalid positions,
and sidecar/species-ledger inconsistencies are hard failures.

## Distributed ingestion and ownership

MPI ownership is workflow-owned. The generic file-set inspector does not create
or initialize MPI. For `world_size > 1` the workflow passes its borrowed
`RuntimeServices::mpi_context` into the distributed reader.

The distributed path:

1. inspects and validates the file set once and broadcasts the strict contract;
2. assigns deterministic bounded chunks by global chunk index;
3. reads and converts each source chunk exactly once globally;
4. packs generic particle data and required sidecars into an explicit
   little-endian versioned wire record;
5. routes each record directly to its deterministic x-slab initial owner in
   bounded rounds;
6. constructs only the final rank-local `SimulationState`;
7. performs exact distributed duplicate-ID validation by hash partition;
8. collectively validates species counts, species mass totals, overall totals,
   coverage, ownership exclusivity/completeness, manifest agreement, finite
   state, and sidecar invariants.

There is no root full-universe gather, all-rank full-state broadcast, or
per-rank authoritative allocation sized to the global particle count. The
startup result is marked `already_partitioned`, so the legacy replicated-state
ownership compactor is skipped.

## Import counters

Each rank reports files/chunks assigned, bytes and records read, records
converted/routed, bytes sent/received, measured peak staging bytes, and final
local particle/gas/sidecar counts. The reference workflow emits these through
`io.ic_ingestion.summary`; see `profiling.md`.

## Current narrow limitations

- Imported gas metallicity and smoothing length are inspected and audited, but
  remain unsupported authoritative lanes until `SimulationState` owns them.
- The canonical converter currently writes one canonical output HDF5 file;
  input discovery and runtime ingestion are multifile-capable.
- Distributed tests require an MPI implementation and launcher. Builds without
  MPI retain the complete serial multifile/converter path and report the
  distributed capability as unavailable.
