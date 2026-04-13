# Initial-condition reader compatibility and assumptions

This document records the explicit schema and conversion assumptions used by `cosmosim::io::readGadgetArepoHdf5Ic`.

## Supported container and groups

- File format: HDF5 with `/Header` and `/PartType[0..5]` groups.
- Required `/Header` attributes:
  - `NumPart_ThisFile` (6 entries)
  - `Time` (scale factor, positive)
- Optional `/Header` attribute:
  - `MassTable` (6 entries, used when `Masses` is absent and fallback is enabled)

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

## Conversion assumptions (conservative and explicit)

- Imported coordinates are interpreted as comoving coordinates and then converted to target frame
  (`comoving` or `physical`) according to `config.units.coordinate_frame`.
- Imported base units are assumed to be `kpc`, `msun`, and `km_s`.
- Conversions use `core::makeUnitSystem`, `comovingToPhysicalLength`, and
  `physicalToComovingLength`.

If a producer uses non-standard unit metadata, pre-convert externally or extend the reader with
explicit unit extraction before ingestion.

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

## Current limitations (explicit and narrow)

- `PartType0/Metallicity` is recognized but currently reported as unsupported because the current
  `SimulationState` schema has no gas-metallicity ownership lane.
- `PartType0/SmoothingLength` is recognized but currently reported as unsupported because the current
  `SimulationState` schema has no gas smoothing-length ownership lane.
- This reader currently imports particle-centric IC payloads only.
