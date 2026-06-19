# AMR Hydro Cross-Code Comparison Protocol

Date: 2026-06-18

## Data contract

The canonical manifest schema is `cosmosim_amr_profile_v1`, documented in
`validation/reference/amr_hydro_canonical_profile_schema_v1.json`. It requires source code and version,
case, geometry, gamma, units, final time, boundary conditions, enabled gravity/cooling/cosmology flags,
and provenance. Shock profiles use `x_code,rho_code,vel_x_code,pressure_code`; Sedov profiles use
`radius_code,rho_code,radial_velocity_code,pressure_code`.

`tools/run_amr_hydro_validation.py compare` validates compatible metadata, resamples profiles onto a
common physical coordinate grid, computes L1/L2 metrics, and reports the strongest density-gradient
coordinate. It rejects mismatched EOS, units, time, boundary conditions, geometry, or physics flags.
It never compares AMR rows cell-for-cell.

## External execution policy

External solver execution is opt-in. A user must provide an actual external executable/container or already
produced external profile, an exact input deck, a converter into the canonical profile contract, version and
provenance information, and the expected output path. No solver is bundled or inferred. A synthetic fixture
tests parser/metric behavior only and is explicitly named `synthetic_contract_fixture`; it is not AREPO,
RAMSES, ENZO, Athena++, or an external validation result.

## Required future comparison protocol

For any real comparison retain: exact initial conditions, EOS and units, boundary conditions, resolution
and AMR policy, final time, hydro solver/limiter/reconstruction, sampling coordinates, metric definitions,
reference checksums, compiler/hardware notes, and an interpretation of method-dependent differences.
Moving-mesh and AMR/static-mesh methods are not expected to match cell by cell. A disagreement is failed,
inconclusive, or interpretable only after this metadata is checked.

## Status

Comparison infrastructure and synthetic contract tests exist. No actual RAMSES, ENZO, AREPO, Athena++,
or other external-code output was supplied or run in this pass. No cross-code agreement claim is made.
