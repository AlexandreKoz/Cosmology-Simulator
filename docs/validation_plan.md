# Validation plan

CosmoSim uses a four-level validation ladder: unit, integration, regression, convergence.

## Validation structure in-repo

- Unit: `tests/unit/`
- Integration: `tests/integration/`
- Validation ladder: `tests/validation/`
- Frozen tolerances: `validation/reference/validation_tolerances_v1.txt`

`validation/input_decks/` is currently documentation/reference material only. It is **not** the authoritative executable runtime path in this repository build. The real runnable path is the config-driven application entry point (`cosmosim_harness <config.param.txt>`) plus the tested release/runtime smoke configs.

## Required command sequence

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

For HDF5 snapshot/restart flows and the real runtime smoke path:

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

For PM/TreePM FFTW-specific validation when FFTW is available:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

## Gravity Phase 1 evidence classes (explicit references)

- **Open-boundary small-`N` tree check**
  - Method: direct summation, softened with the same `epsilon_comoving` as tree solver.
  - Claim: local tree operator correctness for non-periodic pair accumulation.
- **Periodic PM plane-wave mode check**
  - Method: analytic periodic spectral solution for selected Fourier mode(s).
  - Claim: PM operator/sign/phase correctness in periodic domain.
- **Periodic PM uniform-density cancellation**
  - Method: exact lattice plus deterministic parity-balanced micro-jitter (glass-like proxy) particle placement, periodic Poisson solve with mean subtraction.
  - Claim: near-zero force under uniform density.
- **Periodic TreePM consistency check (integration test)**
  - Method: minimum-image periodic direct reference.
  - Claim: seam, cutoff, split, and distributed protocol regression only.
  - Explicit limitation: the minimum-image reference is not Ewald exact and is
    not the periodic force-accuracy certification.
- **Independent periodic Ewald reference and TreePM accuracy gate**
  - Method: test-only double-precision real-plus-reciprocal Ewald sum with
    independent alpha/truncation convergence, homogeneous-background `k=0`
    convention, central-self removal, and rectangular boxes.
  - Claim: small-N periodic force accuracy for the explicitly named production
    profile. Current certification is FFTW-backed TSC with matched window
    deconvolution, `asmth_cells=1.25`, `rcut_cells=6.25`, quadrupoles, and
    COM-distance opening at the typed/public runtime default `theta=0.7`.
    Lower positive cutoffs remain compatibility inputs, not certified defaults.
  - Gate: `relative_L2 <= 1e-2`, `p99_normalized <= 5e-2`, net-force fraction
    `<= 1e-12`, and seam translation drift `<= 1e-11` on the covered fixtures.
  - Limitation: CIC is measured but diagnostic for this envelope; it is not
    silently accepted when it exceeds the TSC gate.
- **Zoom long-range correction + contamination diagnostics check (integration test)**
  - Method: synthetic high-resolution membership mask with focused PM correction enabled; compare high-res vs low-res target force deltas against no-zoom baseline and assert contamination counters.
  - Claim: explicit zoom force decomposition is wired, high-res membership is not implicit, and contamination failure modes remain visible.
- **Periodic TreePM force-error mapping (benchmark artifact)**
  - Method: fine spectral PM plus exact pairwise short-range residual periodic proxy reference.
  - Explicit limitation: not Ewald exact.
- **Split-kernel matching checks**
  - Method: Gaussian split SR/LR factor complementarity and TreePM PM-only/tree-only/split comparisons.
  - Claim: split coupling semantics are wired and monotone in expected direction.
- **Two-body orbit controlled run + energy drift monitor**
  - Method: leapfrog integration with tree accelerations and softened two-body energy accounting.
  - Claim: bounded long-run drift in a controlled setup.
- **Static halo radial-force profile**
  - Method: tree force on shell probes versus softened direct-summation reference.
  - Claim: radial-force profile agreement for a static mass distribution.
- **Regression-only checks (AMR reflux / SF mass budget)**
  - Claim: implementation stability against historical behavior, not standalone validation proof.

## Force-error mapping artifact requirement

Run:

```bash
./build/pm-hdf5-fftw-debug/bench_tree_pm_force_error_map
```

Expected deterministic artifact:

- `validation/artifacts/tree_pm_force_error_map.csv` (generated artifact; do not rely on a source-controlled copy)

The CSV is the required Phase 1 force-error map over PMGRID/ASMTH/RCUT and must be attached/referenced in review evidence.

## Validation invariants

- No silent comoving/proper or code/SI unit mixing.
- Mode and boundary policy must remain explicit and validated.
- Config unit coverage must reject an otherwise valid isolated/open deck when
  `numerics.treepm_enable_window_deconvolution` is `true` (including inherited
  periodic default `true`) and must accept the same deck when the key is
  explicitly `false`; no silent runtime override is allowed.
- The runtime application must load a real config and emit a real run directory outcome.
- Restart roundtrip must preserve state consistency and provenance.
- Snapshot schema compatibility and alias-read behavior remain tested.
- Deterministic/reproducible config normalization and hashes remain stable for unchanged inputs.

## Expected tests for feature work

- **Unit:** parser/math/layout or local invariants.
- **Integration:** real pipeline path through owning modules.
- **Regression:** preserve prior output metadata/behavior with tolerances.
- **Convergence:** when numerical order or resolution behavior is affected.

## Hydro classical CI guards

The default CPU validation ladder now includes:

- Executable: `test_validation_hydro_classics`
- CTest entry: `validation_hydro_classics`
- Labels: `validation;integration;hydro`

This executable runs small fixed-grid Cartesian versions of Sedov, Noh, Gresho
vortex, Kelvin-Helmholtz, and an Evrard-style collapse toy. The assertions are
deliberately CI-scale: every case checks finite state, positive density,
positive pressure/internal energy, and at least one physically meaningful
invariant or trend. These tests are not paper-grade reference-profile or
convergence campaigns.

The Evrard-style case uses a fixed analytic inward gravity source through the
hydro source-term interface. It is a bounded source-coupling and collapse-trend
guard until a self-gravitating hydro workflow can own a restartable Evrard
reference path with documented gravity/hydro ownership.

## AMR hydro CI guards

The default CPU integration ladder includes three AMR hydro validation guards:

- `integration_amr_hydro_shock_tube`: refines the patch adjacent to a Sod-like
  discontinuity, advances through production AMR hydro, and checks finite
  positive state, coarse-fine ghost fill, flux-register/reflux completion,
  conservation bounds, and shock/contact trend.
- `integration_amr_hydro_sedov`: refines a compact Sedov-style blast fixture,
  verifies child states are nonzero after refinement, advances the refined
  hierarchy, and checks finite positive state, approximate symmetry, outward
  shell motion, and total-energy bound.
- `integration_amr_synchronization_stress`: evolves, refines, evolves across a
  coarse-fine boundary with flux registers and reflux, derefines through
  conservative restriction, evolves again, and checks global conserved totals.

These tests are CI-scale synchronization and conservation sentinels for H3 AMR
hydro. They intentionally do not replace future high-resolution AMR convergence
decks, restart equivalence decks, or cross-code shock/blast comparisons.

Restart comparison blocker: these AMR hydro fixtures currently validate the
in-memory production patch/regrid/ghost/reflux path only. A short AMR
run/restart/run comparison remains a follow-up requirement once the HDF5-enabled
restart gate owns an AMR patch-hierarchy fixture with explicit patch, gas-cell
identity, and reflux-synchronization expectations.

## Validation reporting in PRs

Every PR should include:

- build preset(s),
- exact test commands,
- pass/fail status,
- any intentionally deferred validation and rationale.

## Runnable smoke/config examples

Use one of these config-driven paths for honest runtime smoke checks:

- `configs/release/release_smoke_zoom_in.param.txt`
- `configs/release/release_smoke_isolated_galaxy.param.txt`
- `configs/release/release_smoke_cosmo_cube.param.txt`

Or run the built-in integration smoke gate:

- `integration_runtime_app_smoke` (HDF5-enabled builds)
- `integration_runtime_app_mpi_treepm_smoke_two_rank` (MPI+HDF5+FFTW-enabled builds; launches the real `cosmosim_harness` through `mpiexec` with `parallel.mpi_ranks_expected = 2`)



## Phase 3 validation status

- Phase 3 validation closure criteria are defined in `docs/treepm_phase3_contract.md`.
- This repository stage does not yet provide a dedicated `test_validation_phase3_*` executable/gate.
- Current Phase 2 MPI gravity validation is necessary baseline evidence, but not sufficient evidence for Phase 3 closure.

## Phase 2 distributed TreePM validation gate (MPI gravity gate)

- Distributed workflow snapshot/restart paths are rank-qualified (`..._rank###.hdf5`) so MPI restart validation does not rely on unsafe multi-rank writes to a single file.

The Phase 2 gate is now a dedicated MPI validation suite:

- Executable: `test_validation_phase2_mpi_gravity`
- CTest entries:
  - `validation_phase2_mpi_gravity_single_rank`
  - `validation_phase2_mpi_gravity_two_rank`
  - `validation_phase2_mpi_gravity_three_rank`
  - `validation_phase2_mpi_gravity_four_rank`

### Numerical contracts enforced

- Distributed PM equivalence vs one-rank reference: `rel_L2 <= 1e-10`.
- Distributed full TreePM equivalence vs one-rank reference: `rel_L2 <= 5e-6` and `max_rel <= 5e-5`.
- Force-convention contract: the PM and complementary tree return one
  scale-free comoving kernel using `G_code` derived from configured units. The
  non-unit-scale regression requires force ratio `1` between otherwise
  identical `a=0.5` and `a=1` solves; cosmological response belongs to the
  particle/gas update layer rather than an extra kernel `a^2` factor.
- Gravity-timestep coordinate contract: for comoving softening, production
  converts `|A|` to comoving-coordinate acceleration `|A|/a^3` after step
  commit and evaluates `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`.
- Cosmological gas-source contract: at non-unit `a`, unit coverage requires
  `rho A/a^2-Hm` momentum and
  `rho(u dot A)/a^2-H(2K+3P)` energy. A zero-gravity homogeneous-expansion
  case isolates the energy term. Fixed-grid and AMR workflow paths must both
  use the post-drift timeline's step-end `a,H_code`.
- Comoving gravity-timestep public API contract: unit coverage at `a=0.5`
  checks `computeComovingGravityTimeStep` against
  `eta sqrt(a^3 epsilon_com/|A|)` and rejects `a=0`. The production particle
  and gas criteria must use this helper; `computeGravityTimeStep` remains the
  coordinate-neutral primitive.
- Communication stress path: tiny tree exchange batches
  (`tree_exchange_batch_bytes=64`) plus explicit coordinator refresh/local-reuse
  calls, checked against the same TreePM thresholds. This does not enable
  production config cadence greater than one.
- Coordinator cache-safety contract: reuse requires matching
  force epoch, force-evaluation scale, `G_code`, split and rectangular-box
  geometry, assignment, boundary, PM decomposition mode, and deconvolution.
  Particle decomposition is deliberately not a fixed-slab PM-field invalidator.
  Missing/incompatible cache state on one rank must make explicit reuse fail
  coherently; divergent explicit refresh votes must fail before PM collectives.
- Observability contract: focused periodic/workflow tests inspect local
  source/target/tree occupancy, global empty-rank counts, hierarchy/peer
  traffic, PM solve/reuse, halo/local-slab dimensions, and tree
  build/multipole/opened-node counters. Counters are deterministic structural
  evidence, not a substitute for timing-based strong/weak scaling artifacts.
- PM-field event coherence: a refresh event emitted before cadence commit must
  carry the current `PmRefreshDirective` opportunity/version/build step/build
  scale factor, including the initial bootstrap, rather than the previous
  committed `PmSynchronizationState` values.
- Operational precision: end-to-end report coverage requires the first PM
  refresh to carry opportunity/version `1` and rejects `0.000000` formatting
  for nonzero `gravitational_constant_code` and
  `tree_relative_force_acceleration_floor`. Gravity diagnostic doubles use
  scientific `max_digits10` strings.
- Restart continuation contract in MPI mode: reference workflow restart write/read roundtrip must report `restart_roundtrip_ok=true`, preserve rank-qualified restart naming (`..._rank###.hdf5`) for multi-rank runs, preserve distributed ownership/slab metadata exactly across write/read, and pass a production direct-vs-resumed continuation test under the same topology. Rank-count-changing restart is an explicit negative test and must fail as unsupported.
- Production time-integration honesty contract: normalized config requires
  `hierarchical_max_rung=0`; resume rejects nonzero scheduler max/committed bins
  or nonzero pending bins. Mixed-rung KDK is not accepted until per-element
  kick/drift epochs exist.
- Gravity invalid-state trap contract:
  - NaN/Inf force or PM diagnostic norms must raise fatal gravity-state runtime errors.
  - Illegal gravity sync-state transitions (cadence/version/opportunity regressions) must raise fatal gravity-state runtime errors.
  - Distributed restart compatibility must reject gravity-specific metadata mismatches (cadence state and field-version/refresh coherence).
- Periodic TreePM geometry contract: typed config and direct coordinator tests
  must reject `r_cut >= min(Lx,Ly,Lz)/2`, while the certified DMO/Ewald fixture
  keeps the cutoff strictly inside that one-minimum-image bound.
- Distributed workflow honesty floor: the reference workflow MPI test must prove that the final runtime state is partitioned rather than silently replicated by checking (a) reduced local particle/cell counts equal the reported global counts and (b) the reduced particle-ID sum/xor matches the deterministic generated-IC set.
- Required MPI CI lane: `mpi-hdf5-fftw-debug` must select the registered two-rank periodic PM/TreePM workflow tests, uneven PM slab/potential routing coverage, two/three/four-rank gravity validation, distributed gas migration, hydro interface conservation, AMR boundary/reflux, and distributed restart continuation. Missing MPI/HDF5/FFTW-MPI dependencies or insufficient runner rank capacity are blockers, not silent skips.

### Phase 2 scaling artifacts

Run in MPI-enabled builds:

```bash
cmake --build --preset build-mpi-hdf5-fftw-debug --target generate_mpi_gravity_scaling_artifacts
```

Outputs (artifact files):

- `validation/artifacts/pm_only_scaling_np1.csv`
- `validation/artifacts/pm_only_scaling_np2.csv`
- `validation/artifacts/tree_only_scaling_np1.csv`
- `validation/artifacts/tree_only_scaling_np2.csv`

For final integration hard-gates and exact command bundle, see:

- `docs/treepm_phase2_closeout.md`

## 2026-07-13 gravity hardening gates

The current gravity acceptance bundle adds these deterministic gates:

- `validation_periodic_ewald_reference`: Ewald convergence, alpha independence,
  integer-box translation invariance, zero self force, symmetric cancellation,
  sign, `G` scaling, and invalid-input rejection.
- `validation_tree_pm_ewald_accuracy`: quasi-uniform, clustered, rectangular
  seam, Gaussian-transition, and cutoff-separation classifications with full
  error distributions and tree-work diagnostics. The release gate applies only
  to the named default TSC/deconvolution/quadrupole/COM profile; diagnostic
  records cross CIC/TSC and deconvolution off/on, monopole/quadrupole and all
  three MACs, and axis/diagonal separation classes without relaxing the gate.
  An asserted independent-target MAC matrix calibrates
  `theta=sqrt(alpha)=0.07071067811865475` for `alpha=0.005`, requires internal
  accept/open work with zero particle pairs, and keeps all three MAC errors
  within 25 percent while each meets fixed Ewald gates.
- `integration_tree_pm_coupling_periodic_mpi_{two,three,four}_rank`: periodic
  hierarchy geometry, empty/source-only/target-only ranks, independent targets,
  all-empty rounds, migration-to-empty epochs, repeated exchange epochs, and
  distributed agreement. Divergent layout/context metadata must be rejected on
  every rank at coordinated solve entry; periodic unwrap/local-tree and
  PM halo-cache commit, compact active/PM/zoom target preparation, worst-case
  traversal-stack reservation, and exchange-workspace/request preparation must
  remain failure-voted before the next collective.
- `integration_pm_slab_halo_exchange_mpi_{two,three,four}_rank`: FFTW-compatible
  fixed-block slab ownership, same-peer left/right orientation, zero-width
  slabs, and nonblocking halo completion.
- the periodic PM `Nx=2,np4` case: zero-width high ranks enter the complete
  FFTW-MPI plan/transform/collective solve with backend-safe dummy allocation.
- `validation_phase2_mpi_gravity_{single,two,three,four}_rank`: identical
  physical fixtures across rank counts and the unchanged strict TreePM
  thresholds `relative_L2 <= 5e-6`, `max_relative <= 5e-5`.
- `validation_dmo_zeldovich_workflow_single_rank` plus MPI np2/np3/np4 and
  `validation_dmo_zeldovich_workflow_rank_equivalence`: real config loader,
  production KDK/TreePM, HDF5 restart, stable ownership, periodic linear-mode
  amplitude/phase/coherence, mass/momentum/COM checks, direct-vs-resumed
  equality, and stable-ID np1--np4 comparison. A manual cross-preset comparator
  separately checks a true non-MPI serial binary against MPI world size one.

The DMO displacement-growth gate is increment-based. With
`Delta D_expected=a_final/a_initial-1`, both the direct-mode growth and
`sqrt(P_final/P_initial)` must grow and satisfy
`abs((D_measured-1)-Delta D_expected)/Delta D_expected <= 0.075`. Their two
measured increments must agree within relative `1e-3`. The final velocity is
also compared with the exact zero-force Hubble-drag trajectory
`u_ballistic=u_initial a_initial/a_final`: response/ballistic RMS must lie in
`[0.01,0.03]`, its growing-mode projection must lie in `[0.01,0.03]`, and
alignment must be at least `0.99`. The exact initial `4^3` lattice is a named
periodic-Ewald force fixture, whose short-step kick predicts response near
`0.019` of ballistic RMS. This is an independently anchored integration
regression, not an analytic LCDM growth proof.

Command-backed results from the MPI+HDF5+FFTW build on 2026-07-13:

- the focused ownership/halo/periodic/Phase-2 np1--np4 matrix passed 12/12;
- the DMO producer/comparator matrix passed 5/5;
- the certified Ewald profile's worst observed `relative_L2` was
  `8.059667017e-3` and worst `p99` was `8.157551429e-3`, both on the exact
  cancellation-dominated `4^3` DMO initial lattice. Its net-force fraction was
  `5.985844233e-17`, translation drift `1.523007584e-13`, and
  Ewald/PM/tree/total force L2 was
  `10.16780756` / `8.296819318` / `1.789203432` / `10.08586460`. The separate
  rectangular seam fixture retained translation drift `3.031142721e-14`;
  default-setting diagnostic cutoff-classified relative L2 values
  were `6.209e-3` (axis),
  `8.072e-3` (diagonal), `6.152e-3` (immediately beyond cutoff on-axis), and
  `2.397e-3` (1.1 cutoff on the diagonal);
- DMO direct-mode amplitude grew from `0.0049999895833398901` to
  `0.0050596390031745634`, or `1.0119299088208957` versus expected
  `1.0119399775244009`. Maximum Fourier phase drift was
  `7.2112604976761041e-14` radians with coherence `1`; uninterrupted/resumed
  phase/state drift was zero. Response/ballistic RMS was
  `0.018731640935373602`, growing-mode projection `0.01873163714589127`, and
  alignment `0.99999979769618985`.

The DMO gate retains the independent direct particle-mode diagnostic and now
also runs a deterministic binned three-dimensional `P(k)` estimator. The
estimator deposits on a `12^3` CIC mesh, deconvolves the matched assignment
window, uses 32 linear-k bins to the three-dimensional mesh corner, excludes
DC, accounts for all 1727 remaining discrete modes, retains four empty bins,
and reports without subtracting the `V/N = 0.015625` Poisson level. Its
versioned JSON records normalization, code units, assignment/window/shot-noise
policies, bin edges, centers, powers, mode counts, and explicit empty-bin state.

The np1--np4 gate reported fundamental powers
`2.1315894140899176e-6` initially and `2.1827523417862502e-6` after evolution,
giving `sqrt(P_evolved/P_initial) = 1.011929959672577` versus expected
scale-factor growth `1.0119399775244009`. After normalizing the intentionally
rank-specific `world_size`, all four JSON artifacts had SHA-256
`53ee518b63d21cd28790415b8f076614497536bb6a62c95276b8ed6c14984bf3`.
The cross-preset serial/MPI1 comparison covered 64 stable IDs with periodic
position error `0`, maximum velocity error
`0`, and relative mass error `0`; their JSON artifacts
were byte-identical with SHA-256
`8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`.
This closes the missing-artifact blocker for the controlled gate, but the
direct-DFT implementation and tiny linear fixture are not a performant
large-mesh estimator or publishable/high-dynamic-range DMO campaign. The
estimator's normalization, window, bin, unit, and empty-bin contract is
maintained in `docs/power_spectrum_diagnostics.md`.

## Phase 3 research-grade evidence campaign

For Phase 3 gravity maturity evidence generation and separation of correctness/force-accuracy/scaling classes, use:

- Observable contract: `validation/campaign/phase3_gravity_observables.json`
- Campaign command bundle: `docs/validation_phase3_campaign.md`
- Provenance collector: `scripts/validation/collect_phase3_evidence.py`

The campaign requires explicit cosmological and zoom runtime runs plus force/time-integration/scaling artifacts. Existing `np1/np2` scaling outputs remain baseline performance evidence only and are not standalone certification.


## MPI AMR directed-exchange validation update

The MPI feature lane must include the real MPI-launched AMR directed-exchange
registration `integration_reference_workflow_distributed_amr_mpi_two_rank` and
must not count the serial `integration_amr_distributed_remote_patch_boundary`
fixture as MPI coverage. The serial test remains useful as an AMR data-structure
and orchestrator guard only.

The current two-rank AMR MPI test exercises the new directed patch/cell payload
exchange and owner-routed flux-register exchange with asymmetric ranks and
nonzero mass/momentum/energy payload fields. It is not yet full AMR workflow
acceptance: a later validation pass still needs to prove production
`ReferenceWorkflowRunner` AMR restart continuation, named migration/regrid across
ranks, and all five conserved totals through direct and restart-resumed paths.
