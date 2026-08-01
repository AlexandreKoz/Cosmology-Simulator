# RFC 0001: adaptive bound–Jeans star formation

Status: implemented reference model

## Decision

Adopt `adaptive_bound_jeans` as the explicitly selected reference star-particle birth model while preserving `legacy_schmidt_threshold` for compatibility. The model requires local self-gravity, thermal Jeans instability, and optionally converging flow; uses exact finite-step depletion; and performs conservative, batched gas-to-star transfer keyed by stable gas-cell identity.

## Rationale

A fixed physical density threshold imposes a resolution-dependent preferred scale. The selected criteria instead test support and unresolved collapse at the actual AMR cell scale. The implementation remains deliberately narrower than a sink-particle system and does not combine unvalidated turbulence, magnetic, chemistry, or molecular-fraction closures.

## Data ownership

- `GasCellIdentityMap` is the stochastic parent identity authority.
- patch geometry and hydro state are production inputs;
- `gas_cells.metal_mass_code` is the conserved metallicity authority;
- generic particle SoA plus star sidecars own newborn state;
- scheduler and distributed ownership ledgers are updated at the legal source boundary;
- canonical snapshot and versioned restart schemas own persistence.

No module-owned zero-filled divergence or metallicity vector is accepted on the production runtime path.

## Determinism and distribution

The RNG key excludes rank and dense row. Plans are normalized by stable birth identity. Exact global ownership validation catches cross-rank ID collisions, duplicate births, missing IDs, and extra IDs. Restart continuation preserves the global tick and birth identity fields.

## Conservation interpretation

Mass, momentum, and metal mass are transferred to roundoff. Gas internal energy removed with converted mass is recorded as an explicit source sink. Collisionless stars carry kinetic energy through their particle state but no hydrodynamic internal-energy variable.

## Consequences

The source stage requires fresh synchronized patch geometry, gradients, metals, ownership, and scheduler truth. Gas metal mass becomes a cross-cutting conserved scalar in hydro, AMR, MPI, cooling, feedback, I/O, and restart. The implementation adds state and test complexity but removes conflicting metallicity authorities and decomposition-dependent stochastic behavior.

## Non-goals

Sink particles, accretion, merging, individual stars, molecular chemistry, MHD, radiation transport, cosmic rays, and full galaxy calibration remain separate campaigns.
