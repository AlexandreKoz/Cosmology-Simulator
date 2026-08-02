# RFC 0002: Effective multiphase star formation

**Status:** accepted for the CHUÍ reference implementation  
**Model schema:** star formation v3; effective EOS v1  
**Parameter set:** `chui_sh03_tng_like_v1`

## Decision

CHUÍ shall expose three typed star-formation models: explicit legacy compatibility, adaptive bound--Jeans collapse for resolved or semi-resolved ISM structure, and an SH03/TNG-inspired effective multiphase closure for unresolved-ISM simulations. The latter is named `effective_multiphase_tng_like` to communicate lineage without claiming exact IllustrisTNG equivalence.

Eligibility/rate physics and thermodynamic closure are separated from the shared star-birth pipeline. Both production models use the same stable RNG key, exact global particle-ID precommit, conservative gas-to-star transfer, scheduler/stellar sidecars, gravity invalidation, diagnostics, snapshot, and restart contracts.

## Rationale

A single local-collapse model is not scientifically optimal across CHUÍ's intended dynamic range. At high resolution the adaptive model can use resolved gradients and binding. At coarse resolution those quantities are noisy and the unresolved ISM fragments artificially. A pressure-supporting equilibrium closure is a defensible, cheap alternative for population-volume work.

## Thermodynamic contract

The effective table is immutable, hashed, and constructed from frozen units, composition, cooling policy, and typed parameters. The selected pressure and \(\partial P/\partial\rho\) propagate through primitive conversion, reconstruction, Riemann fluxes, and CFL. Conserved density/momentum/energy remain authoritative; cold fraction and effective pressure are derived and are recomputed after AMR restriction/prolongation or restart.

Hot gas above the equilibrium branch retains resolved shock energy and cannot form stars prematurely. Relaxation toward the branch is an explicit source with a separate energy ledger.

## Identity and atomicity

Particle IDs are precommitted before mutation. The registry validates all existing IDs, all local candidates, and—under MPI—all rank candidates. Deterministic collision ordinals rehash only identity, never the stochastic birth decision. If global precommit fails, gas is unchanged.

The registry uses sorted contiguous vectors rather than per-ID node allocation. Source evaluation has no collective; distributed precommit is one batched protocol after local planning.

## AMR and MPI

Only owned active leaf cells spawn. Covered coarse cells are rejected, source application follows hydro/reflux synchronization, and derived EOS state is reconstructed from conserved quantities. Stable gas-cell identity—not dense row or rank—keys physics and stochastic decisions.

The repository includes production-runtime two-rank tests for owner-only spawning, an empty eligible rank, distributed ID registry, restart, AMR ownership/migration path, one-versus-two-rank equivalence contract, and EOS hash equality. These require an MPI+HDF5 environment to execute.

## Compatibility

Historical configs without a selector retain the documented legacy interpretation, but every shipped star-forming production profile now names its model explicitly. Old restart files retain their encoded model; no restart silently changes closure.

## Consequences and non-goals

The effective closure increases disk pressure support and can change fragmentation and scale height. EOS and explicit feedback must be calibrated together. This RFC does not implement exact TNG winds, individual-star IMF sampling, sink accretion, molecular chemistry, MHD, cosmic rays, radiation hydrodynamics, or observational calibration.
