# Star formation

CHUÍ provides three explicitly named stellar-population birth models:

- `legacy_schmidt_threshold`: the historical density/temperature/convergence prescription retained for old configurations, restarts, and controlled comparisons.
- `adaptive_bound_jeans`: the resolved/semi-resolved collapse model for high-resolution zoom-in and isolated-disk work.
- `effective_multiphase_tng_like`: an SH03/TNG-inspired unresolved-ISM model for coarse cosmological baryonic volumes and population studies.

Star particles represent coeval stellar populations. They are not individual stars, sink particles, accreting clusters, or molecular-cloud control volumes.

## Adaptive eligibility

Only active, locally owned, non-ghost leaf gas cells with positive finite mass, density, and volume are considered. Production inputs come from authoritative patch geometry, gas state, velocity fields, `GasCellIdentityMap`, and conserved gas metal mass.

For stored comoving density and geometry,

\[
\rho_{\rm phys}=\rho_{\rm comov}/a^3,\qquad
\Delta x_{\rm phys}=a\,V_{\rm comov}^{1/3}.
\]

For proper storage, the scale-factor conversions are omitted. The effective cell length is the cube root of the actual cell volume. The production velocity-gradient tensor is evaluated from peculiar velocities using physical patch spacing, so its entries have inverse code-time units. The model records this as `physical_inverse_time_patch_finite_difference`.

The local support parameter is

\[
\alpha_{\rm vir}=
\frac{\lVert\nabla\otimes\mathbf v\rVert_F^2+2(c_s/\Delta x_{\rm phys})^2}
{8\pi G\rho_{\rm phys}},
\]

and eligibility requires `alpha_vir < sf_bound_alpha_vir_max`. No magnetic or subgrid-turbulent support is fabricated when those validated state variables do not exist.

The thermal Jeans mass convention is

\[
M_J=\frac{\pi^{5/2}}{6}\frac{c_s^3}{G^{3/2}\rho_{\rm phys}^{1/2}},
\]

with eligibility requiring

\[
M_J < \max(M_{J,\rm floor},M_{\rm cell}).
\]

When `sf_require_converging_flow=true`, the peculiar-flow divergence must be negative. The compression time is `-1/div(v)`.

## Rate and timestep contract

The free-fall time is

\[
t_{\rm ff}=\sqrt{\frac{3\pi}{32G\rho_{\rm phys}}}.
\]

`sf_collapse_timescale` selects either `free_fall` or `minimum_free_fall_or_compression`. The exact finite-step expected converted fraction is

\[
f_*=1-\exp(-\epsilon_{\rm ff}\Delta t/t_{\rm sf}).
\]

The implementation uses `expm1` for small arguments. The scheduler registers the corresponding source limit

\[
\Delta t_{\max}=-\frac{t_{\rm sf}}{\epsilon_{\rm ff}}\ln(1-f_{\max}),
\]

including the compression time when that collapse-time mode is selected. The mutation path still checks bounds defensively, but clipping is not the primary timestep policy.

## Stochastic spawning and identity

The target mass is fixed or a fraction of the parent gas resolution. The sampler uses integer-plus-Bernoulli spawning with an adjusted feasible target when needed so both possible outcomes obey gas-mass, minimum-particle-mass, maximum-particle-mass, and maximum-count constraints without biased post-draw clipping.

The counter-based draw is keyed by:

- global star-formation seed;
- stable `gas_cell_id`;
- global integration tick;
- birth-attempt ordinal;
- RNG key-schema version.

It excludes dense row, MPI rank, thread number, pointer value, and iteration order. Birth keys additionally include the model schema version. Numeric particle IDs are deterministic hashes of birth keys, use a reserved generated-ID domain, are collision-checked locally before mutation, and are checked by the exact distributed ownership ledger before migration and at workflow completion.

## Conservative transfer

For total birth mass `dm`, the newborn population inherits the parent centroid, peculiar velocity, and metallicity. Gas mass, density, pressure, and conserved metal mass are reduced consistently while gas specific internal energy and velocity remain unchanged. The ledgers record:

- equal gas mass removed and star mass created;
- equal gas momentum removed and stellar momentum created;
- equal gas metal mass removed and stellar metal mass created;
- removed gas internal energy as `star_formation_internal_energy_sink_code`;
- stellar kinetic energy represented by newborn particle mass and velocity.

The code does not claim conservation of hydrodynamic internal energy across conversion: that removed thermal energy is an explicit sink ledger, not a collisionless particle degree of freedom.

After a production birth batch, the source runtime refreshes generic gas-particle compatibility/gravity mirrors from the authoritative owned gas cells. This prevents the next gravity boundary from seeing both pre-birth gas mass and newborn stellar mass; multiple gas cells sharing one lineage particle are aggregated deterministically by stable gas-cell ID.

Births are planned, sorted by stable gas identity, validated, allocated in one batch, appended in one particle resize and one stellar-sidecar resize, and followed by one species-index rebuild. New stars enter scheduler bin zero at the next legal tick, receive no retroactive force kick, and become visible at the next force boundary.

## Metals, AMR, MPI, restart, and output

`metal_mass_code` is the authoritative gas metallicity scalar. Hydro fluxes, AMR prolongation/restriction/refluxing, temporal boundary state, MPI ghost exchange/migration, cooling, stellar feedback, snapshots, and restarts use that same conserved authority. Derived metallicity is `metal_mass_code / gas_mass_code` with zero-mass protection.

Only owned leaf cells spawn. Covered coarse cells and ghosts are rejected. Stable gas-cell identity survives row reorder and migration. Source-created particle IDs are appended to a distinct ownership-origin ledger before any rebalance, so migration is not mistaken for particle creation.

Restart schema 22 persists gas metal state, AMR metal registers/history, star birth identity, scheduler state, RNG metadata, and source sidecars. Schema 21 remains explicitly loadable with zero/default initialization for fields that did not exist there.

Snapshots use canonical GADGET-style names. Gas writes `PartType0/Metallicity`; stars write `PartType4/Metallicity`, `StellarFormationTime`, `BirthMass`, and the CHUÍ birth-identity datasets documented in `output_schema.md`. Science-light diagnostics write `sfr_history.csv` from actual stellar birth mass and formation scale factor.

## Configuration

Mass-valued `*_code` parameters accept either a raw code-unit number or an explicit `kg`, `g`, or `msun` suffix and are normalized to the configured code mass unit. See `configs/adaptive_bound_jeans_isolated_galaxy.param.txt`.

The adaptive model does not use a fixed density threshold as its physical trigger. `sf_temperature_safety_ceiling_k=0` disables the optional numerical ceiling. Historical threshold keys retain their meaning only in legacy mode.

## Scope and limitations

This is a validated reference star-particle birth model, not a complete calibrated galaxy-formation model. It does not implement sink accretion/merging, individual-star IMF sampling, molecular chemistry, radiation hydrodynamics, MHD support, cosmic rays, or a new feedback calibration. Galaxy-scale credibility still depends on cooling, stellar evolution, feedback, enrichment, hydrodynamic response, resolution studies, and comparison against external reference calculations.

## Particle-ID collision contract

Numeric star-particle IDs are deterministic hashes of the immutable birth key. Each source stage sorts and checks the complete local newborn-ID batch before mutation. The workflow then performs exact global duplicate, missing, and unexpected ID validation through the distributed ownership acceptance ledger. This avoids a full scan of the existing particle array and avoids per-birth allocator nodes while retaining an explicit collision-failure contract.


## Shared birth pipeline and exact ID precommit

All production models share one plan/precommit/apply backend. Local physics determines eligibility and expected mass; the backend owns deterministic sampling, exact collision-checked IDs, conservation, batch append, scheduler and stellar sidecars, gravity invalidation, and diagnostics. Particle IDs are checked against existing and same-batch IDs before any gas mutation. MPI precommit gathers immutable birth keys once after planning and resolves collisions deterministically without changing physical birth decisions.

## Effective multiphase model

The equations, pressure closure, cooling interaction, parameter provenance, feedback policy, and limitations are documented in `docs/effective_multiphase_ism.md` and RFC 0002. The effective pressure and derivative signal speed enter the hydro reconstruction/Riemann/CFL path; they are not merely eligibility diagnostics.

## Profile policy

Every shipped configuration with `enable_star_formation = true` explicitly declares `star_formation_model`. High-resolution zoom and resolved-disk examples select `adaptive_bound_jeans`; coarse cubes and population-oriented references select `effective_multiphase_tng_like`; only clearly named compatibility fixtures select the legacy model.
