# Initial-condition reader compatibility and assumptions

This document records the explicit schema and conversion assumptions used by `cosmosim::io::readGadgetArepoHdf5Ic`.

## Supported container and groups

- File format: HDF5 with `/Header` and `/PartType[0..5]` groups.
- Required `/Header` attributes:
  - `NumPart_ThisFile` (6 entries)
  - `NumPart_Total` and optional `NumPart_Total_HighWord` (6 entries each;
    normalized to canonical 64-bit totals)
  - `NumFilesPerSnapshot`
  - `MassTable` (6 entries)
  - `BoxSize`
  - `Time` (scale factor, positive)
  - `Redshift`, `Omega0`, `OmegaLambda`, and `HubbleParam`

Header counts, cosmology, cubic box size, and start epoch must agree with the
validated runtime configuration before state allocation is committed.

## Dataset aliases

Aliases are accepted to improve interoperability with common GADGET/AREPO derivatives:

- Coordinates: `Coordinates`, `Position`, `POS`
- Velocities: `Velocities`, `Velocity`, `VEL`
- Particle IDs: `ParticleIDs`, `ParticleID`, `ID`
- Masses: `Masses`, `Mass`
- Gas internal energy: `InternalEnergy`, `U`, `Internal_Energy`
- Gas density: `Density`, `Rho`
- Gas metallicity: `Metallicity`, `GFM_Metallicity`
- Gas smoothing length: `SmoothingLength`, `Hsml`, `Smoothing_Length`

The exact alias selected is stored in `IcImportReport::present_aliases` for provenance and audits.

## Versioned manifest and conversions

Every external read owns an `IcManifest` in `IcImportReport`. The manifest
records deterministic source ordering/provenance IDs, dialect/version, per-file
and 64-bit total counts, high words, header cosmology, field path/type/rank/
dimensions/count, base-SI conversion, `h` and scale-factor exponents, frame,
velocity convention, extensive/intensive/specific semantics, disposition, and
the six particle-type policies.

The reference workflow writes the validated import contract as
`<run_directory>/ic_manifest.json` before integration begins. Generated ICs do
not claim an external import manifest.

With no caller manifest, direct import selects the named
`gadget_arepo_bridge_v1` contract: kpc, Msun, km/s, zero `h`/`a` exponents,
comoving coordinates, and physical peculiar velocities. This is a versioned
CHUI convention, not an inference from the word “GADGET.” Producers using
`h`-scaled fields, physical coordinates, sqrt(a)-scaled peculiar velocities,
or comoving coordinate rates must supply an explicit validated manifest.
Runtime particle positions remain comoving and velocities remain physical
peculiar values. Density and specific internal energy are converted through
their derived SI dimensions, not copied as raw numbers.

`IcManifest` validation rejects inconsistent per-file totals/high words,
nondeterministic or duplicate sources, invalid field shapes/conversions,
inconsistent Time/Redshift, nonphysical cosmology/box data, and populated
particle types whose species policy remains `reject`.

## Missing-field behavior

- Missing required datasets raise an exception.
- Missing optional fields are recorded in `IcImportReport::missing_optional_fields`.
- Defaulted values are recorded in `IcImportReport::defaulted_fields`, including:
  - generated particle IDs
  - velocity zero-fill
  - mass fallback from `MassTable`
  - gas `InternalEnergy` zero-fill when missing
  - gas `Density` zero-fill when missing

## Gas thermodynamic mapping behavior

- `PartType0/InternalEnergy` is ingested into `SimulationState::gas_cells.internal_energy_code`.
- `PartType0/Density` is ingested into `SimulationState::gas_cells.density_code`.
- Gas particles are mirrored into the gas-cell skeleton for current hydro ownership:
  - cell centers from imported particle coordinates,
  - cell mass from imported particle mass.
- `GasCellSidecar` fields without direct IC equivalents in this reader revision
  (`pressure_code`, `temperature_code`, `sound_speed_code`, reconstruction gradients) are initialized to zero.

## Species policy

- PartType0 maps to gas and PartType1 to dark matter.
- PartType2 and PartType3 are rejected by default. A caller may explicitly map
  each family to dark matter using its distinct manifest policy; this makes the
  loss of external family identity an auditable choice.
- PartType4 maps to stars. Full stellar formation/metallicity import remains a
  limitation below.
- PartType5 maps to black holes and requires canonical `BH_Mass`; optional
  `BH_Mdot` defaults to zero with an audit record. Both are converted into the
  restart-authoritative CHUI black-hole sidecar.
- Duplicate IDs and nonpositive/nonfinite particle or black-hole masses are
  hard failures.

Floating datasets must have the declared rank and `[count]` or `[count,3]`
shape; particle IDs must be a rank-one integer dataset. Header/group counts and
all selected hyperslabs are range checked.

## Current limitations (explicit and narrow)

- `PartType0/Metallicity` is recognized but currently reported as unsupported because the current
  `SimulationState` schema has no gas-metallicity ownership lane.
- `PartType0/SmoothingLength` is recognized but currently reported as unsupported because the current
  `SimulationState` schema has no gas smoothing-length ownership lane.
- This reader currently imports particle-centric IC payloads only.
- Direct runtime import is deliberately single-file. A manifest declaring more
  than one source fails with guidance instead of reading the whole global IC on
  every rank. Multifile discovery/routing and distributed materialization
  remain unsupported and are reported as such in `runtime_capabilities.json`.
- Canonical CHUI conversion tooling and typed-config manifest-file selection
  are not yet implemented; direct external import therefore remains
  provisional.
- PartType4 formation-time and metallicity sidecar import is not yet complete.
