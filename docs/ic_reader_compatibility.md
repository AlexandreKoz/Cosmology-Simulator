# Initial-condition ingestion, conversion, and compatibility

The external-IC path is an explicit scientific contract shared by typed
configuration, the runtime reader, the compiled converter, and audit manifest
schema v3. No GADGET/AREPO-looking filename supplies implicit units, coordinate
frame, velocity convention, h scaling, scale-factor powers, or species policy.

Authoritative interfaces:

- `include/cosmosim/io/ic_reader.hpp`
- `src/io/ic_manifest.cpp`
- `src/io/ic_reader_file_set.cpp`
- `tools/convert_ic.cpp`
- `tools/convert_ic.py` (launcher only)
- `src/workflows/initial_condition_runtime.cpp`

## Typed convention selection

`mode.ic_convention` is one of `generated`, `chui_canonical_v1`,
`gadget_arepo_bridge_v1`, or `manifest_v1`.

A generic direct bridge is valid only when configuration supplies the complete
source contract:

- source length, mass, and velocity SI-per-code factors;
- `comoving` or `physical` source coordinates;
- `physical_peculiar`, `sqrt_a_scaled_peculiar`, or
  `comoving_coordinate_rate` source velocity;
- finite h and scale-factor exponents for length, mass, and velocity;
- explicit PartType2 and PartType3 policy.

Zero or unspecified values are fail-closed sentinels. `manifest_v1` requires
`mode.ic_manifest_file` and obtains the contract from the strict manifest.
Canonical input validates the CHUÍ canonical header. Generated input remains
independent of external conventions.

## Deterministic multifile discovery

The bridge accepts `/Header` plus populated `/PartType0` through `/PartType5`.
`NumFilesPerSnapshot` determines an exact ordered sequence such as
`snapshot.0.hdf5` through `snapshot.(N-1).hdf5`; arbitrary directory iteration is
not used.

Every member must agree on cosmology, box, epoch, `MassTable`,
`NumFilesPerSnapshot`, reconstructed 64-bit totals, and schema. Per-file counts
must equal dataset dimensions and sum to declared totals. Missing, repeated, or
inconsistent members fail before authoritative state allocation.

In MPI mode, file index modulo rank count assigns schema inspection and SHA-256
work. Every file is opened and hashed exactly once globally; rank zero assembles
metadata fragments but does not reread the complete set.

## Safe HDF5 schema inspection

Header attributes are inspected before `H5Aread`:

- count and `MassTable` attributes require exactly six elements;
- scalar cosmology, box, epoch, file count, and canonical provenance attributes
  require scalar dataspaces;
- numeric class, width, sign, and representable range are validated.

Selected physical datasets must be float32 or float64. Particle and tracer
identity/index datasets must be unsigned integers no wider than uint64. Floating
IDs, signed IDs, unsupported compound/string/reference data, wrong ranks, and
wrong vector component counts fail closed rather than relying on implicit HDF5
conversion.

The manifest records actual class, width, signedness, byte order, rank,
dimensions, record count, and selected alias. Every dataset in a populated
supported PartType has a converted or explicitly dropped disposition. Unknown
populated PartType families fail closed. Ambiguous simultaneous aliases fail.

Accepted black-hole aliases are implemented consistently by runtime and
converter:

- mass: `BH_Mass`, `BlackHoleMass`;
- accretion rate: `BH_Mdot`, `BlackHoleAccretionRate`.

## Audit manifest schema v3

`IcManifest` uses `chui_ic_audit_manifest` version 3. Version 2 is rejected
because it cannot fully represent density, specific energy, black-hole
accretion, frame transforms, and velocity-convention powers.

For each field, v3 records:

- source path and alias;
- real HDF5 schema;
- SI base factor;
- h and scale-factor exponent;
- physical dimension `L^p M^q T^r`;
- additional physical-to-comoving frame exponent;
- velocity-convention multiplier power;
- coordinate frame and velocity convention;
- converted/preserved/defaulted/dropped/rejected disposition;
- source/target unit descriptions and canonical equation.

Manifest validation enforces the physical dimensions of coordinates,
velocities, masses, density, internal energy, BH accretion, identifiers, and
dimensionless fields. Source files, byte sizes, SHA-256 hashes, original header
records, family policies, warnings, defaults, dropped fields, and converter
provenance are also retained.

A supplied manifest is authoritative only after ordered files, hashes, sizes,
headers, counts, policies, and actual HDF5 schema are revalidated. All ranks
agree on the serialized manifest digest before materialization.

## One authoritative conversion engine

Runtime ingestion and `cosmosim_convert_ic` call the same C++ functions:

```text
target = stored
       * base_unit_to_si
       * h^hubble_exponent
       * a^(scale_factor_exponent + frame_scale_factor_exponent)
       * velocity_convention_multiplier^velocity_convention_power
       / target_si_per_code(L^length_power M^mass_power T^time_power)
```

The dimensional contract includes:

- coordinates `L`;
- velocities `L/T`;
- masses and stellar/BH/tracer mass fields `M`;
- density `M/L^3`;
- internal energy `L^2/T^2`;
- BH accretion `M/T`;
- IDs and indices as unitless integer identity;
- formation epoch and metallicities as explicit dimensionless values.

The audit counterexample (`a=0.5`, `h=0.5`, physical coordinates,
`sqrt(a)` velocity, length/mass h exponent `-1`, `U=4`, `rho=2`) is a permanent
regression. Direct bridge, authoritative manifest, and canonical converter must
produce the same gas and black-hole values.

## Species and sidecars

- PartType0 creates gas particles, cells, identity, velocity, density, and
  internal-energy lanes.
- PartType1 maps to dark matter.
- PartType4 creates complete star sidecars including formation epoch, initial
  mass, and metallicity.
- PartType5 creates complete black-hole sidecars including subgrid mass and
  accretion rate.
- PartType2/3 follow explicit `reject`, `dark_matter`, `star`, `black_hole`, or
  `tracer` policy. A mapping is accepted only when required fields and sidecars
  are complete.

Unsupported source datasets are retained in the audit as explicit dropped
fields; unsupported populated families fail.

## Compiled canonical converter

Build target: `cosmosim_convert_ic`.

`tools/convert_ic.py` and `tools/convert_ic_manifest.py` are compatibility
launchers; they do not implement scientific equations.

Direct example:

```bash
tools/convert_ic.py \
  --input snapshot.0.hdf5 \
  --output chui_ic.hdf5 \
  --manifest chui_ic.ic_manifest.json \
  --source-convention gadget_arepo_bridge_v1 \
  --source-length-unit-to-si 3.0856775814913673e22 \
  --source-mass-unit-to-si 1.98847e30 \
  --source-velocity-unit-to-si 1000 \
  --coordinate-frame comoving \
  --velocity-convention physical_peculiar \
  --length-h-exponent -1 --length-a-exponent 0 \
  --mass-h-exponent -1 --mass-a-exponent 0 \
  --velocity-h-exponent 0 --velocity-a-exponent 0
```

Manifest example:

```bash
tools/convert_ic.py \
  --source-manifest source.ic_manifest.json \
  --output chui_ic.hdf5 \
  --manifest chui_ic.ic_manifest.json
```

The converter validates and hashes the complete source set, performs exact ID
checks, reads in bounded chunks, writes canonical GADGET-compatible names,
embeds the manifest digest, and atomically finalizes the HDF5 and JSON outputs.
The canonical artifact is interchange IC state, not a restart checkpoint.

## Distributed runtime status

### Tracer host identity after routing

A source `HostCellIndex` is a file-local row number and cannot remain authoritative
after chunking and MPI routing. The reader records that source field as dropped
for authoritative-state construction. Tracer `ParentParticleIDs` are preserved,
and finalization resolves each tracer host to the owner-local gas row with the
matching stable particle ID. Import fails if the parent gas cell is absent from
the final rank or if the mapping is ambiguous. This keeps tracer identity stable
without retaining source-row assumptions.

The distributed reader uses the workflow-borrowed MPI context and never
initializes MPI. It assigns file inspection and chunks across ranks, routes
explicit records directly to final x-slab owners, performs exact file/chunk and
source/final identity reconciliation, validates all coordinates and sidecars,
and reports truthful I/O and staging counters.

The implementation and registered tests exist, but capability is reported as
`provisional` until the full MPI+HDF5+FFTW acceptance preset runs successfully.
See `docs/architecture/distributed_ic_ingestion.md` for collective phase order
and fault-injection coverage.
