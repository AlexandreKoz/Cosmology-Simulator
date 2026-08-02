# Hydrodynamics and thermodynamic closures

CHUÍ's finite-volume hydro state keeps conserved mass, momentum, total energy, and passive metal mass authoritative. Primitive pressure and signal speed are derived through one thermodynamic-closure interface.

The ordinary path uses the ideal-gas closure. When `effective_multiphase_tng_like` is selected and physical density exceeds the configured threshold, `EffectiveIsmThermodynamicClosure` supplies effective pressure, target specific internal energy, and the barotropic derivative signal speed. These values are consumed consistently by primitive caching, MUSCL reconstruction, Riemann wave-speed estimates, and CFL diagnostics.

Resolved hot gas above the effective branch keeps its ideal-gas pressure until the explicit relaxation/cooling policy permits approach to equilibrium. The closure does not overwrite shock energy during primitive conversion.

AMR prolongation, restriction, reflux, and MPI exchange move conserved state. Effective pressure, cold fraction, and signal speed are derived again from synchronized conserved state; they are not transported as independent conserved scalars.
