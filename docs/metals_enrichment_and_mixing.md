# Metals, delayed enrichment, and turbulent mixing

## Status and runtime ownership

CHUÍ uses one authoritative total-metal state in the default `total_only` mode:

- gas mass: `CellSoa::mass_code`;
- gas metal mass: `GasCellSidecar::metal_mass_code`;
- hydro-stage metal density: `metal_mass_density_comoving`, reconstructed from the authoritative metal mass and geometry;
- gas metallicity: the derived ratio `metal_mass_code / mass_code`;
- stellar birth metallicity: immutable `StarParticleSidecar::metallicity_mass_fraction`;
- stellar evolution and unresolved deposition: persistent star-sidecar cumulative and carry lanes.

No independently mutable gas-metallicity field was added. The existing finite-volume passive-scalar, AMR reflux, migration, restart, IC, and scalar `Metallicity` snapshot paths remain the authority for advection and serialization.

The implementation in this release is a **scalar-metals production integration with bounded scope**. `core_elements` is parsed but rejected because a complete element sidecar, yield schema, cooling table, AMR transport, and MPI exchange contract are not yet present. It is intentionally not a decorative mode.

## Governing equation and coordinate treatment

In physical coordinates the modeled scalar obeys

\[
\frac{\partial(\rho Z)}{\partial t}+\nabla\cdot(\rho Z\mathbf v)
=S_Z+\nabla\cdot(\rho\kappa_Z\nabla Z).
\]

Advection is handled by the existing conservative hydro scalar. Source terms transfer metal **mass**, not metallicity. Diffusion is evaluated using physical cell lengths, physical volumes, peculiar velocities, and physical velocity gradients. For a comoving patch at scale factor \(a\), CHUÍ converts stored lengths by \(\ell_{\rm phys}=a\ell_{\rm comov}\) and volumes by \(V_{\rm phys}=a^3V_{\rm comov}\) before constructing the operator. This avoids inserting guessed scale-factor powers into a code-unit diffusion equation.

Each internal diffusion face transfers one amount with equal and opposite signs. The operator changes neither gas mass, momentum, nor energy. A symmetric donor/receiver limiter enforces

\[
0\le M_Z\le M_{\rm gas}
\]

without independent per-cell clipping. Reflective boundaries have zero normal diffusive flux. Open-boundary diffusive escape is not enabled in this implementation; the local operator treats an open exterior face as zero flux rather than claiming a measured loss.

## Gas-to-star transfer

The star-formation path retains the existing conservative rule

\[
\Delta M_{Z,\star}=Z_{\rm gas}\Delta M_\star.
\]

The newborn stellar population records birth mass, formation scale factor, birth metallicity, parent gas-cell ID, deterministic birth key/tick/ordinal, and zero-initialized enrichment/carry/deposition bookkeeping. IC import initializes the same complete state.

## Delayed SSP model

`StellarEvolutionTable` now represents a tensor-product grid in SSP age and birth metallicity. Cumulative fields are stored for each of three channels:

1. stellar winds and AGB mass loss;
2. core-collapse supernovae;
3. Type Ia supernovae.

The v2 table distinguishes:

- cumulative returned total mass;
- cumulative total ejected metals, including returned birth-composition metals;
- cumulative newly synthesized metals;
- cumulative event count;
- cumulative feedback energy.

An interval budget is the difference of cumulative values at its endpoints. The loader validates finite ordered grids, monotonic cumulative quantities, returned fraction no larger than one, metals no larger than returned mass, and channel sums equal to totals.

Age interpolation is logarithmic between positive-age nodes, with explicit handling of age zero. Birth metallicity uses bounded linear interpolation. Values beyond the final age node remain at the final cumulative value; birth metallicity is clamped to the tabulated range.

### Time authority

Cosmological stellar ages use the FLRW proper-time integral

\[
\Delta t=\int_{a_0}^{a_1}\frac{da}{aH(a)}
\]

through `LambdaCdmBackground::cosmicTimeIntervalSi`. Isolated simulations use elapsed physical simulation time from the configured unit system. The old logarithmic scale-factor/Hubble-time proxy is retained only as a deprecated configuration field for normalized-config compatibility; production `SourceRuntime` does not use it.

### Transactional bookkeeping

Production source execution follows this order:

1. evaluate per-star interval budgets without mutating stellar mass;
2. deposit or persist every returned mass, metal, energy, and momentum budget;
3. only after the budget is durably represented, decrease the star mass and advance cumulative SSP counters.

This prevents a restart or no-neighbor condition from losing material after stellar mass has already been reduced. Pending budgets migrate and restart with the star.

The built-in table is deliberately a zero-return, non-calibrated compatibility table. Production users must select and audit a licensed table and may set `stellar_evolution_require_production_table=true` to reject the built-in fallback. CHUÍ does not ship invented “realistic-looking” yields.

## Enrichment deposition

Returned gas mass and returned metal mass use one deterministic target list and one normalized inverse-distance weight set. Ties are broken by distance, stable gas-cell ID, and local row. Only caller-designated owned leaf cells are legal targets.

Deposition is dimensionally explicit:

- cell density is recomputed as new mass divided by cell volume;
- specific internal energy is updated from old total internal energy plus deposited total energy converted from erg to code energy;
- the named kinetic energy budget is conservatively thermalized until a momentum/total-energy injection operator exists;
- momentum is never written into internal energy and remains in persistent carry.

Target absence, stochastic deferral, delayed-cooling energy, and unsupported momentum deposition therefore produce durable carry rather than loss.

Current limitation: the production target search is local and scans its owned-leaf candidate set. Cross-rank enrichment kernels and a bounded spatial-index exchange protocol remain unimplemented; this release does not claim MPI-boundary enrichment completeness.

## Turbulent diffusion

The implemented model is constant-coefficient Smagorinsky:

\[
\kappa_Z=C_{\rm mix}\Delta^2\left|S^\ast\right|,
\]

where \(\Delta=V^{1/3}\), and \(S^\ast\) is the trace-free symmetric peculiar-velocity gradient. Solid-body rotation produces zero strain. The coefficient is a calibration parameter, not a universal constant.

For an internal face,

\[
\dot M_{Z,L\rightarrow R}=\left(\rho\kappa_Z\right)_f
\frac{A_f}{d_{LR}}(Z_L-Z_R),
\]

using a harmonic mean for \((\rho\kappa_Z)_f\). The explicit stability estimate is assembled from each cell's face conductance. The timestep candidate participates in the production scheduler as `gas_cell_metal_diffusion_parabolic`.

Two integrators are available:

- `explicit_subcycling`: second-order SSPRK2 subcycles with the conservative limiter;
- `rkl2`: second-order Runge-Kutta-Legendre super-time-stepping, with a stage-capacity check and fail-closed bounds check.

RKL2 is validated against explicit subcycling for the linear periodic operator. Because high-stage super-time-stepping can be less robust near sharp boundaries, the implementation rejects an unbounded stage result instead of clipping it.

Current production geometry covers owned, leaf, patch-internal Cartesian faces. Same-level inter-patch, coarse-fine, and rank-boundary diffusion flux exchange/reflux are not yet implemented. The configuration consequently fails closed for mixed time rungs, but it does not claim full AMR/MPI diffusion completeness.

## Cooling coupling

The existing cooling path continues to derive scalar metallicity from authoritative gas metal mass. SourceRuntime updates gas authority and compatibility mirrors before later legal source/cooling boundaries. Element-resolved cooling is not silently enabled from scalar metallicity. `core_elements` is rejected until a compatible abundance and cooling-table schema exists.

## Persistent schema and diagnostics

Restart files persist:

- cumulative returned mass and total metals;
- cumulative newly synthesized metals;
- cumulative feedback energy and channel counters;
- unresolved mass, metals, energy, and momentum;
- cumulative deposited mass, metals, and energy.

Snapshots preserve canonical gas and star `Metallicity`. HDF5 headers also record metal species mode, diffusion model, diffusion integrator, diffusion coefficient, and stellar-evolution table path. The normalized configuration and provenance hash remain the reproducibility authority.

`computeMetalBudgetDiagnostics` reports gas metal mass, stellar birth and locked metal estimates, newly synthesized/returned/deposited/carried metals, deposition residual, metallicity extrema and mass-weighted mean, metal-free gas fraction, and invalid state count. `globalMetalAuditResidualCode` evaluates the user-supplied initial/escaped global audit baseline.

## Table provenance and validation data

`resources/stellar_evolution/test_synthetic_v2.txt` is deterministic **test-only** data generated by `tools/generate_stellar_evolution_test_table.py`.

- table ID: `chui_synthetic_ssp_validation`;
- version: `v2_test_2026_08`;
- license: CC0-1.0;
- IMF/mass range/solar scale: not applicable (synthetic);
- tracked channels: winds/AGB, CCSN, SNIa;
- payload SHA-256: `2ac5f1330e423f3e5f9735bcaacc8c1bd31c9832f41c4e4909fabf0482e6f513`;
- full-file SHA-256: `d98478d4ab360dff7a08b90100f346f89dfdcc088ff7a4ff5e50133ce81a5177`.

It must never be used as a production calibration.

No third-party production yield table is redistributed in this release. This avoids silently importing a data license or unsupported transformations. A future production table must record paper, repository/release, license, retrieval date, original hash, IMF, mass/metallicity grids, abundance scale, transformations, interpolation, and extrapolation.

## Scientific rationale and primary references

The design follows the delayed-channel structure used by modern galaxy-formation models rather than instantaneous recycling:

- Pillepich et al., *Simulating galaxy formation with the IllustrisTNG model*, MNRAS 473, 4077 (2018): discrete AGB, core-collapse SN, and SNIa enrichment.
- Wiersma et al., *Chemical enrichment in cosmological, smoothed particle hydrodynamics simulations*, MNRAS 399, 574 (2009), and Wiersma, Schaye & Smith, MNRAS 393, 99 (2009): delayed element production and abundance-dependent cooling.
- Hopkins et al., *FIRE-3: Updated Stellar Evolution Models, Yields, & Microphysics*, MNRAS 519, 3154 (2023), arXiv:2203.00040: updated stellar mass loss, SN rates, yields, and fitting functions.
- Hough et al., *SIMBA-C*, MNRAS 525, 1061 (2023), arXiv:2308.03436: removal of instantaneous recycling in favor of time-dependent enrichment.
- Correa et al., *A subgrid model for chemical enrichment in cosmological simulations*, MNRAS 548, stag645 (2026), arXiv:2604.00980: updated COLIBRE enrichment and turbulent diffusion.
- Rennehan et al., MNRAS 483, 3810 (2019), arXiv:1807.11509: dynamic localized Smagorinsky behavior and the risk of over-diffusion in laminar shear. CHUÍ implements only the lower-cost constant model and documents this calibration limitation.
- Meyer, Balsara & Aslam, JCP 257, 594 (2014): stabilized RKL1/RKL2 super-time-stepping for parabolic operators.

## Known limitations

This release does **not** claim completion of the entire metals campaign. Remaining work includes:

- licensed, calibrated production SSP/yield data;
- optional individual-element storage, transport, yields, snapshot/restart schema, and cooling;
- cross-rank enrichment target discovery and owner-rank deposition;
- same-level patch, coarse-fine reflux, and MPI diffusion faces;
- open-boundary diffusive metal-loss measurement;
- spatial indexing to eliminate the local candidate scan;
- conservative momentum injection rather than persistent carry;
- full workflow restart-equivalence after simultaneous enrichment and diffusion;
- MPI decomposition-equivalence tests and calibrated galaxy-scale validation.
