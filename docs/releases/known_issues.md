# Known issues (initial release)

## Gravity hardening current-status note (2026-07-13)

This file preserves initial-release and dated Phase 3 history below. The
current gravity readiness authority is `docs/gravity_production_readiness.md`.
The repaired CPU/MPI TreePM path now has focused correctness evidence through
four ranks, including empty ranks, exact same-topology restart continuation,
an independent periodic Ewald comparison, and a deterministic binned 3D power
spectrum artifact from the production-workflow Zel'dovich gate.

That spectrum proof uses a `12^3` CIC mesh with matched assignment-window
deconvolution and 32 requested bins. It accounts for all 1727 non-DC mesh
modes, retains four empty bins explicitly, reports rather than subtracts the
`V/N = 0.015625` Poisson level, and emits a versioned JSON artifact. The
  initial/evolved fundamental powers are `2.1315894140899176e-6` and
  `2.1827523417862502e-6`; their square-root growth is
  `1.011929959672577` versus the expected `1.0119399775244009`. After
  normalizing the intentionally rank-specific `world_size`, np1--np4 produced
  SHA-256
  `53ee518b63d21cd28790415b8f076614497536bb6a62c95276b8ed6c14984bf3`.
  A true non-MPI serial/MPI1 comparison produced byte-identical JSON SHA-256
  `8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`.

The controlled gate now compares displacement and square-root-power growth
increments against `a_final/a_initial-1` within 7.5 percent and against each
other within relative `1e-3`. Its exact initial `4^3` lattice is separately
force-certified against periodic Ewald; response/ballistic RMS and growing-mode
projection must each lie in `[0.01,0.03]` with alignment at least `0.99`.
This remains a short controlled integration regression, not an analytic LCDM
or cross-code proof.

The remaining gravity release limitations are not erased by those repairs:

- production PM and short-range exchanges still use communicator-wide
  collectives, and the remote tree path is not a locally essential tree;
- correctness/runtime evidence stops at four ranks and does not constitute a
  strong- or weak-scaling campaign;
- rank-count-changing restart, scalable distributed isolated/open PM, and
  distributed zoom PM without high-resolution-source all-gather are absent;
- production hierarchical KDK is restricted to rung zero until per-element
  kick/drift epochs are implemented and restart-validated;
- focused zoom correction disables unsupported isolated/open window
  deconvolution and is outside the periodic TSC accuracy certification;
- isolated/open user decks must explicitly set
  `numerics.treepm_enable_window_deconvolution=false`; omission inherits the
  periodic default `true` and is rejected during typed config loading rather
  than silently rewritten at runtime;
- periodic decks require
  `r_cut=treepm_rcut_cells*cbrt((Lx/Nx)(Ly/Ny)(Lz/Nz))` to be strictly below
  half the shortest box axis. Older coarse decks at or above the bound now fail
  closed; increase PM resolution or reduce the cutoff and revalidate the
  changed force profile;
- the validated DMO case is a 64-particle, two-step linear-mode check, not a
  large-volume convergence, halo-statistics, cross-code, or publication-grade
  cosmology campaign;
- GPU PM remains CIC-only and is outside the certified TSC plus deconvolution
  force-accuracy profile.

The force convention is now explicit: PM and the complementary short-range
tree use the physical Newton constant converted from configured units and return
a scale-free comoving kernel. Cosmological scale-factor response is owned by
collisionless KDK or the gas conservative source. This repairs the former
hard-coded-`G=1`/extra-`a^2` behavior. The gas source additionally uses the
post-drift timeline's code-unit `a,H` and the general-gamma expansion energy
term `-H(2K+3P)`. These are reproducibility-visible numerical corrections for
physical-unit cosmological runs. Gravity timestep proposals likewise convert
`A` to comoving-coordinate acceleration `A/a^3` before combining it with
comoving softening.

Accordingly, the current path is suitable for controlled engineering
first-light and small-cluster validation, with limitations. It is not evidence
for publishable DMO production or AREPO/GADGET-class gravity parity.

## GPU PM acceleration is experimental

- Scope: CUDA PM path.
- Impact: not release-validated for production science runs.
- Mitigation: treat CUDA runs as exploratory; use CPU/MPI production path for published results.

## Black-hole AGN model calibration is incomplete

- Scope: black-hole AGN toy model.
- Impact: physically plausible exploratory behavior but not yet a validated production calibration.
- Mitigation: disable AGN for baseline release reference runs unless explicitly studying the model's exploratory behavior.

## HDF5 interoperability requires feature-enabled build

- Scope: snapshot/restart I/O workflows.
- Impact: users may assume HDF5 works in all presets, but it is gated by `COSMOSIM_ENABLE_HDF5`.
- Mitigation: use HDF5-enabled presets and follow `docs/build_instructions.md` dependency requirements.

## Large-cluster scaling evidence not included in this initial release

- Scope: high-rank production deployment.
- Impact: architecture supports scale-up, but this release package does not claim full strong/weak scaling certification on large HPC systems.
- Mitigation: use desktop/small-cluster scale for this release; add dedicated scaling campaign in a follow-up release.


## TreePM Phase 3 maturity is not yet implemented

- Scope: integrator-grade hierarchical cadence/PM synchronization maturity and Phase 3 closure evidence.
- Impact: current release evidence supports Phase 1/2 baseline correctness contracts, but not Phase 3 completion claims.
- Mitigation: follow `docs/treepm_phase3_contract.md` hard-gate requirements before using Phase 3-complete language.


## Phase 3 scaling evidence remains baseline-only in this stage

- Scope: strong/weak scaling certification.
- Impact: campaign now exports explicit baseline scaling artifacts and provenance, but certification-grade large-rank sweeps are still pending.
- Mitigation: treat `validation/artifacts/research_grade/phase3/scaling/phase2_baseline_scaling_summary.json` as baseline evidence only and complete the larger-rank follow-on campaign before certification claims.

## Phase 3 closeout hard gate remains open (2026-04-21 audit)

- Scope: final Phase 3 integration closeout.
- Impact: current audit cycle did not reach a full hard-gate pass; periodic/zoom integration and reference-workflow/schema lanes include failing tests, and distributed MPI+FFTW command bundle is environment-blocked.
- Mitigation: resolve blockers listed in `docs/treepm_phase3_closeout.md` before claiming Phase 3 completion.
