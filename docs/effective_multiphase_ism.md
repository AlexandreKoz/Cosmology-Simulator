# Effective multiphase ISM and TNG-like star formation

## Scope

`effective_multiphase_tng_like` is CHUÍ's unresolved-interstellar-medium option for coarse cosmological hydrodynamic volumes, population studies, and runs that do not resolve molecular-cloud structure or trustworthy cloud-scale velocity gradients. It is inspired by the Springel--Hernquist (2003) equilibrium multiphase model and the softened effective equation of state used by the IllustrisTNG model family. It is **not** an exact reproduction of IllustrisTNG: CHUÍ does not claim identical winds, cooling calibration, IMF, feedback coupling, or parameter tuning.

A star-forming gas cell represents an unresolved mixture of cold clouds and a supernova-heated ambient phase. Cold fraction, effective internal energy, pressure, closure signal speed, and star-formation timescale are derived from one immutable density table. Star particles remain coeval stellar populations, not individual stars or sink particles.

## Physical threshold and frames

The configured threshold is a physical hydrogen number density:

```text
physics.sf_effective_n_h_threshold = 0.13 cm^-3
```

It is converted once using

\[
\rho_{\rm th,phys}=n_{\rm H,th}m_p/X_H.
\]

When hydro density is comoving, the runtime converts it to physical density before lookup. The optional `sf_effective_min_baryon_overdensity` guard is explicit and is not used in isolated modes. Tests cover scale-factor invariance of the physical threshold and AMR-level independence.

## Equilibrium closure

At density ratio \(r=\rho_{\rm phys}/\rho_{\rm th,phys}\),

\[
t_\star=t_{\star,0}r^{-1/2},\qquad
A=A_0r^{-4/5},
\]

\[
u_{\rm hot}=u_{\rm cold}+\frac{u_{\rm SN}}{1+A}.
\]

The cooling time is obtained from CHUÍ's authoritative cooling provider under the deterministic reference environment recorded by the table metadata. With the configured massive-star fraction \(\beta\), the equilibrium helper is

\[
y=\frac{t_\star}{t_{\rm cool}}
\frac{u_{\rm hot}}{\beta u_{\rm SN}-(1-\beta)u_{\rm cold}}.
\]

The cold-cloud fraction is evaluated using the stable equivalent form

\[
x_{\rm cold}=\frac{\sqrt{1+4y}-1}{\sqrt{1+4y}+1},
\]

which is algebraically equivalent to the usual GADGET-style equilibrium expression and avoids cancellation at extreme \(y\). Invalid or non-positive equilibrium inputs fail closed.

The full mixture energy is

\[
u_{\rm eff,full}=(1-x_{\rm cold})u_{\rm hot}+x_{\rm cold}u_{\rm cold}.
\]

CHUÍ applies the softened closure in specific internal energy:

\[
u_{\rm eff}=q_{\rm EOS}u_{\rm eff,full}+(1-q_{\rm EOS})u_{\rm iso}.
\]

The reference `q_EOS=0.3` is TNG-like. `q_EOS=1` selects the full multiphase branch and `q_EOS=0` the configured isothermal branch. The table is normalized for pressure continuity at the threshold.

Pressure and the closure signal speed are

\[
P_{\rm eff}=(\gamma-1)\rho_{\rm phys}u_{\rm eff},\qquad
c_{\rm eff}^2=\frac{\partial P_{\rm eff}}{\partial\rho}.
\]

The derivative is tabulated from the normalized pressure curve. It reaches primitive-state construction, reconstruction, Riemann wave speeds, and CFL estimation; it is not a diagnostic-only pressure floor.

## Table construction and provenance

`EffectiveMultiphaseEosTable` is built once after units and cooling configuration are frozen. It uses logarithmic density spacing, monotonic linear interpolation, positive derivative checks, no lookup allocation, and a deterministic FNV-based content hash. Provenance records:

- EOS and star-formation schema versions;
- parameter-set name (`chui_sh03_tng_like_v1` by default);
- density range and bin count;
- EOS table hash;
- cooling reference description;
- composition and softening parameters.

Local metallicity continues through ordinary cooling, enrichment, and stellar evolution. It does not silently alter the one-dimensional reference EOS table cell by cell.

## Hot gas, cooling, and energy accounting

Dense shock-heated gas above the configured tolerance relative to the equilibrium branch is not immediately classified as star forming. The closure preserves resolved excess thermal energy and suppresses birth until the cell approaches the branch. Gas on or below the branch is relaxed upward to the effective floor; finite-timescale relaxation is available as a typed option.

Effective-EOS adjustments are explicit source terms and use separate ledgers:

```text
effective_ism_energy_added_code
effective_ism_energy_removed_code
effective_ism_net_energy_adjustment_code
```

They are not merged with the star-birth internal-energy sink, ordinary radiative cooling loss, or resolved stellar-feedback injection. Latest and cumulative values are restart persistent in the `effective_multiphase_ism` module sidecar.

## Star-formation rate and birth mass

The unresolved cold mass is \(M_{\rm cold}=x_{\rm cold}M_g\). For the default `initial_ssp_mass` convention,

\[
\dot M_{\star,\rm birth}=\frac{x_{\rm cold}M_g}{t_\star}.
\]

The existing stellar-evolution module then owns prompt and delayed mass return. The alternative `long_lived_mass` convention includes the factor \((1-\beta)\) and must only be used with a compatible mass-return setup. Configuration validation prevents known double-counting combinations.

Finite-step expectation is

\[
f_{\rm expected}=1-\exp\left[-\frac{\dot M_{\star,\rm birth}}{M_g}\Delta t\right].
\]

Both adaptive and effective models use the same deterministic sampling, exact particle-ID precommit, conservative mutation, scheduler insertion, gravity invalidation, metal transfer, and output path.

## Feedback compatibility

`external_feedback_calibrated` retains CHUÍ's explicit feedback together with effective pressure support and treats them as a coupled calibration problem. It is not TNG wind feedback. `effective_eos_only` rejects incompatible explicit feedback at configuration time while retaining documented stellar mass/metal return semantics.

## Recommended regimes

| Regime | Recommended model |
| --- | --- |
| High-resolution zoom or resolved disk | `adaptive_bound_jeans` |
| Coarse cosmological baryonic cube | `effective_multiphase_tng_like` |
| Population survey / parameter sweep | `effective_multiphase_tng_like` |
| Historical comparison | `legacy_schmidt_threshold` |
| Dark-matter-only | Star formation disabled |

Neither model is universally superior. The adaptive model responds to resolved dynamics but is gradient- and resolution-sensitive. The effective model is stable and inexpensive at coarse resolution but suppresses unresolved structure and uses pseudo-thermodynamic phase averages.

## Limitations

This implementation does not resolve molecular clouds, cold/hot interfaces, individual stars, sinks, magnetic/turbulent support, radiation, TNG kinetic winds, or observational calibration. Effective temperatures are pressure-support diagnostics, not resolved single-phase temperatures. The retained validation products are analytic and controlled-model evidence, not a full galaxy-population calibration.
