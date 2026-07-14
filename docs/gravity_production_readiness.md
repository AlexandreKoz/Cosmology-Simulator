# Gravity production readiness

_Current-state assessment: 2026-07-13. Historical Phase 2/Phase 3 closeout
documents remain historical records; this file is the current gravity status._

## Direct verdict

| Claim | Status | Basis and boundary |
| --- | --- | --- |
| Engineering first-light ready | PASS WITH LIMITATIONS | Core PM/tree/TreePM/restart paths build and have CPU plus MPI+HDF5+FFTW correctness evidence. Scientific and scaling gates below still limit the claim. |
| DMO validation ready | PASS WITH LIMITATIONS | The production-workflow Zel'dovich gate passes at np1--np4 with exact restart equivalence; it is a controlled 64-particle, two-step linear-mode validation, not a production cosmology campaign. |
| Small-cluster MPI validation ready | PASS WITH LIMITATIONS | Focused PM ownership/halo, periodic TreePM, strict Phase-2, and DMO gates run through four ranks. The transport remains communicator-wide and has no beyond-np4 scaling certification. |
| Publishable DMO production ready | FAIL | The controlled workflow now emits a validated binned 3D power-spectrum artifact, but no large-volume convergence ladder, halo statistics, cross-code comparison, or production-scale campaign exists. |
| AREPO/GADGET-class gravity comparable | FAIL | CHUI does not yet have mature LET scaling, large-rank evidence, production cosmological validation, GPU assignment parity, or the validation history required for that comparison. |

## Readiness matrix

The status vocabulary is intentionally restricted to the project-approved
values.

| Surface | Status | Evidence and limitation |
| --- | --- | --- |
| PM serial | PASS WITH LIMITATIONS | Periodic FFTW operator, CIC/TSC matched assignment/gather, scale-free comoving kernel, physical `G_code` unit conversion, and potential/force checks pass; the TSC PM branch participates in the Ewald-certified composite TreePM force. Certification covers small deterministic fixtures, not a large DMO campaign. |
| PM single-rank MPI | PASS WITH LIMITATIONS | MPI build/runtime path executes and matches the same owner/layout contracts. One rank does not exercise communication scaling. |
| PM multi-rank MPI | PASS WITH LIMITATIONS | FFTW-compatible fixed-block slabs and the transposed-pencil spectral path are implemented with np1--np4 correctness coverage. Explicit PM wire records, routed density/interpolation, exact response coverage, and halo np2/np3/np4 pass. A dedicated `Nx=2,np4` case proves zero-width ranks enter the full FFTW solve. Sparse records still move through communicator-wide `MPI_Alltoallv`. |
| Standalone tree | PASS WITH LIMITATIONS | Direct-sum tests cover clustered/uniform distributions, three MACs, corrected quadrupoles, softened second moments, and work/error trends. It is deliberately nonperiodic and its regression tolerances are not a universal science envelope. |
| Periodic TreePM serial | PASS WITH LIMITATIONS | Seam-safe largest-gap tree geometry, corrected screened multipoles, cutoff, relative MAC, and TSC-vs-Ewald gate pass. CIC is diagnostic for the certified envelope. |
| Periodic TreePM single-rank MPI | PASS WITH LIMITATIONS | Same production coordinator executes in the MPI feature build with Ewald and workflow validation. |
| Periodic TreePM multi-rank MPI | PASS WITH LIMITATIONS | np2/np3/np4 periodic and strict Phase-2 equivalence pass, including empty and target-only ranks. Hierarchy all-gather plus target all-to-all is small-cluster only, not LET. |
| Empty-rank MPI | PASS | Explicit zero-source hierarchy sentinel, zero-count collective participation, source-only/target-only/all-empty rounds, zero-width PM slabs, and migration-to-empty epochs are exercised through four ranks. |
| Restart continuation | PASS WITH LIMITATIONS | Production DMO direct-vs-resumed state/digest is exact at np1--np4. The actual decomposition epoch is restored, PM mesh/tree scratch rebuilds, and committed force history is stable-ID remapped or invalidated after migration. Only rung-zero, same-rank-count/topology restart is supported. |
| Ewald force validation | PASS | Independent rectangular, unsoftened, double-precision Ewald reference converges and checks alpha independence, self force, symmetry, translation, sign, and `G` scaling. It is test-only and small-N. |
| Periodic force accuracy | PASS WITH LIMITATIONS | Certified TSC+deconvolution, `asmth=1.25`, `rcut=6.25`, runtime-default `theta=0.7` profile passes `relative_L2 <= 1e-2` and `p99 <= 5e-2`; worst values are `8.059667017e-3` and `8.157551429e-3` on the exact cancellation-dominated `4^3` DMO initial lattice. CIC and the former 4.5-cell cutoff are diagnostic. |
| Small DMO validation | PASS WITH LIMITATIONS | 64-particle periodic Zel'dovich mode uses production config/KDK/TreePM/HDF5 restart. Its exact cancellation-dominated initial lattice is a named Ewald force fixture, so the trajectory gate is no longer accepted on seeded-velocity growth alone. Increment-relative direct displacement-mode and power gates reject zero growth; a narrow Ewald-anchored short-step response check covers force magnitude and direction. It also checks Fourier phase/coherence, stable-ID physical state, exact resume, true non-MPI serial versus MPI1, and MPI np1--np4 equivalence. |
| Power-spectrum validation | PASS WITH LIMITATIONS | The DMO gate uses deterministic CIC deposition on a 12-cubed mesh, matched window deconvolution, explicit Fourier normalization/units/shot-noise policy, 32 requested bins including empty bins, and a versioned JSON artifact. It is a small controlled linear-mode validation, not a high-dynamic-range survey estimator. |
| Large DMO production | NOT IMPLEMENTED | No large-volume resolution/convergence/scaling/science campaign is certified. |
| Small-cluster production | EXPERIMENTAL | Correctness evidence reaches four ranks, but performance/scaling, fault behavior, and publishable science validation are not closed. |
| Zoom TreePM | EXPERIMENTAL | Focused high-resolution correction exists with contamination diagnostics and a hard all-gather byte limit. Its isolated/open correction explicitly disables unsupported window deconvolution even when global PM uses it. It lacks Ewald/cosmological certification and scalable source transport. |
| Isolated/open PM | PASS WITH LIMITATIONS | Single-rank doubled-domain convolution is supported. Typed config rejects window deconvolution for this nonperiodic operator, and isolated decks must explicitly set it false rather than relying on a silent override. Distributed open PM is only a bounded root gather/solve/scatter compatibility path and is not scalable. |
| GPU PM | EXPERIMENTAL | CUDA 12.0 builds and passes 123/123, including the PM smoke; `sm_52` uses a checked CAS-compatible double density atomic. CUDA PM remains CIC-only, with no TSC parity or Ewald-certified GPU profile. |
| Publishable DMO production | FAIL | Controlled first-light validation is not a publication-grade cosmological result. |
| State-of-the-art gravity parity | FAIL | LET, high-rank scaling, large DMO observables, adaptive-softening validation, and mature GPU parity remain absent. |

## Implemented correctness contracts

### Periodic geometry and multipoles

- Periodic TreePM builds a transient per-axis largest-gap unwrapped source
  frame. Tree topology, root size, COM, quadrupole, raw second moment, MAC, node
  bounds, and remote pruning agree in that frame; particle truth remains
  wrapped and caller-owned.
- Axis-aware minimum images preserve rectangular boxes. Integer-image
  translations, x/y/z seams, edges/corners, interior/seam target-source
  combinations, leaf/internal interactions, monopoles/quadrupoles, and
  geometric/COM MACs have deterministic regression coverage.
- Periodic config and direct coordinator entry reject
  `r_cut >= min(Lx,Ly,Lz)/2`; the residual tree evaluates one minimum image per
  source. This is a validity bound, not an accuracy-tolerance relaxation.
- The Newtonian quadrupole uses the corrected sign for
  `d=x_com-x_target`. Softened and Gaussian-screened accepted nodes use a true
  second-order raw-moment expansion with derivatives of the full radial kernel.
- The relative MAC uses
  `G M l^2 <= alpha max(|a_previous|,a_floor) r^4`, opens a target-containing
  node, and falls back to COM-distance when compatible history is absent.

### MPI ownership and protocol

- FFTW-MPI's fixed ceil-block x ownership is the PM authority map; a truncated
  tail and zero-width high ranks are legal.
- Halo exchange derives logical peers from global plane owners and uses
  receive-first nonblocking point-to-point completion with sequence-qualified
  side tags.
- Hierarchy, request, and response packets are explicit version-1 little-endian
  wire records. They carry rank/record identity and decomposition, force, and
  exchange epochs; payload shape, finite values, bounds, peer identity,
  duplicates, missing responses, and byte/displacement overflow are validated.
- PM density, force-request/response, and potential-request/response traffic
  likewise uses explicit fixed-width version-1 little-endian records. The
  actual MPI world votes operation kind and serial-versus-distributed layout
  mode before any layout branch or collective.
- Every rank contributes one hierarchy root. Empty trees contribute an explicit
  zero-source sentinel. Independent target coordinates allow a source-empty
  rank to own work without dummy mass.
- All ranks enter globally coordinated request and response batches even when
  their byte counts are zero.
- Periodic unwrap/local-tree construction failures are reduced across the
  actual MPI world before residual collectives. Layout/context mismatch is
  rejected by solve-entry consensus rather than a rank-local constructor
  throw.
- PM halo-cache commit, compact active/PM/zoom target preparation, and
  worst-case residual traversal-stack reservation are also actual-world
  failure-voted before subsequent collectives.
- Reusable TreePM exchange workspace owns counts, displacements, payloads,
  remote acceleration scratch, and expected/received response coverage.
  Throwing capacity growth occurs inside coordinated protocol phases and is
  exposed as transient MPI-buffer capacity in the coordinator memory report.
- Zero-width FFTW-MPI ranks retain logical extent zero with backend-safe dummy
  allocation and still participate in plan creation, transforms, and every PM
  collective.

This protocol is fail-fast, not fault tolerant. A process failure or arbitrary
rank-local exception outside the validated collective entry contract is not
recovered; MPI job abort remains the failure boundary.

### Time, restart, and reproducibility

- Production config accepts only PM cadence one and
  `hierarchical_max_rung=0`. The latter blocks mixed-rung KDK until
  per-element kick/drift epochs exist. Every integrator-issued,
  rank-coordinated production force-refresh surface rebuilds PM.
- Explicit lower-level coordinator reuse requires an exact transient signature
  match for force epoch, field-build scale factor, `G_code`, split scale,
  rectangular box axes, assignment, boundary, PM decomposition mode, and
  deconvolution. Missing/incompatible cache on any rank makes reuse fail rather
  than silently triggering an unrequested solve; mixed explicit votes fail
  before PM collectives.
- The pre-commit `gravity.pm_long_range_field` event reports the pending
  refresh directive's opportunity/version/build step/build scale factor. It no
  longer pairs a refresh message with stale previous-field metadata; cadence
  authority and restart semantics are unchanged.
- Gravity operational-event doubles use scientific `max_digits10` formatting
  for physical `G_code`, MAC floors, derived length scales, force norms, and
  scale factors. Audit artifacts no longer round small nonzero evidence to
  `0.000000`.
- The actual particle-decomposition epoch advances only after a committed
  ownership transition and is restored from restart. It protects tree wire
  records but is intentionally not a PM mesh invalidator because FFT slab
  ownership is fixed. Dense-row acceleration history is invalidated immediately
  after the transition.
- PM and the short-range residual return one scale-free comoving kernel with
  the physical Newton constant converted from configured units. Neither inserts
  `a^2`. Collisionless KDK and the gas conservative source apply the
  cosmological response to their respective state.
- Post-drift fixed-grid and AMR gas updates use
  `timeline_step.scale_factor_end/hubble_end_code` and apply
  `S_m=rho A/a^2-Hm`,
  `S_E=rho(u dot A)/a^2-H(2K+3P)`. The still-step-begin committed
  `IntegratorState` and independently recomputed SI Hubble rates are not source
  truth at that stage.
- Post-step adaptive gravity criteria use the committed scale factor and
  pass scale-free `|A|`, `epsilon_com`, and `a` through the public
  `computeComovingGravityTimeStep` helper. It validates finite positive `a` and
  converts to comoving-coordinate acceleration `|A|/a^3`:
  `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`.
- PM mesh arrays, tree geometry/topology, hierarchy records, and exchange
  buffers are transient. Restart persists integrator/distributed PM metadata
  and the committed stable-ID-keyed gravity-force cache, then deterministically
  rebuilds transient fields.
- TSC plus matched window deconvolution and `rcut_cells=6.25` are the normalized
  defaults. With `asmth_cells=1.25`, this gives `r_cut/r_s=5`. The former
  4.5-cell cutoff showed about nine-percent error at the cutoff transition.
  These default changes are reproducibility-relevant; explicit legacy
  CIC/deconvolution-off/lower-cutoff decks retain their values but are not the
  certified profile.
- Opening criterion, theta, relative tolerance, and acceleration floor are
  typed, validated, normalized, and recorded in `provenance_v6`.

## Command-backed evidence snapshot

The 2026-07-13 feature environment provided MPI, HDF5, FFTW, and FFTW-MPI. The
focused gravity run executed real MPI launchers through four ranks.

The final fresh/incrementally rebuilt preset matrix was:

- CPU-only debug: 122/122 passed;
- HDF5 debug (HDF5 1.10.10): 124/124 passed;
- PM+HDF5+FFTW debug (FFTW 3.3.10): 124/124 passed;
- MPI+HDF5+FFTW debug (Open MPI 3.1, np1--np4 registrations): 150/150 passed;
- ASan+UBSan+LSan debug: 122/122 passed with leak detection enabled.
- CUDA debug (CUDA 12.0.140, configured `sm_52` compatibility): 123/123
  passed, including the GPU PM smoke.

Repository hygiene and `git diff --check` also passed. The final focused
gravity/MPI selection passed 36/36, including SFC np2/np3/np4,
halo np2/np3/np4, periodic TreePM np2/np3/np4, strict Phase-2 np1--np4, and
DMO serial/MPI rank equivalence. MPI-enabled serial executables are explicitly
covered without calling communicator APIs before `MPI_Init`.

DMO workflow matrix:

```bash
TMPDIR=/tmp ctest --test-dir build/mpi-hdf5-fftw-debug \
  -R '^validation_dmo_zeldovich_workflow' --output-on-failure
```

Result: 5/5 passed: the MPI-world-one producer, MPI np2/np3/np4 producers, and
their stable-ID rank-equivalence comparator. The physical mode evidence was:

- initial direct fundamental amplitude: `0.0049999895833398901`;
- final direct fundamental amplitude: `0.0050596390031745634`;
- measured direct-mode growth: `1.0119299088208957`;
- expected scale-factor growth: `1.0119399775244009`;
- maximum Fourier phase drift: `7.2112604976761041e-14` radians;
- initial/final mode coherence: `1`;
- maximum np1--np4 periodic position error: `0`;
- maximum np1--np4 velocity error: `2.7755575615628914e-17`;
- uninterrupted/resumed phase and restart-authoritative state drift: `0`.
- response/ballistic RMS: `0.018731640935373602`;
- growing-mode response projection: `0.01873163714589127`;
- response alignment: `0.99999979769618985`.

The acceptance rule compares displacement/power growth increments, not just
amplitudes:
`abs((D_measured-1)-(a_final/a_initial-1))/(a_final/a_initial-1) <= 0.075`
for both the direct mode and square-root power growth. Their increments agree
within relative `1e-3`. A separate zero-force Hubble-drag baseline requires
response/ballistic RMS and growing-mode projection in `[0.01,0.03]`, with
alignment at least `0.99`. The exact initial `4^3` cancellation-dominated
lattice is separately certified against periodic Ewald; its short-step kick
predicts a response near `0.019` of the ballistic RMS. These bounds are an
Ewald-anchored integration regression, not an analytic LCDM growth proof.

The detailed spectrum estimator reported the same values at np1--np4:

- mesh and requested bins: `12^3` CIC and `32` linear-k bins;
- non-DC discrete Fourier modes accounted for: `1727`, with `4` requested
  bins explicitly empty;
- CIC Poisson shot-noise level `V/N`: `0.015625` code-volume units;
- initial recovered density amplitude: `0.0050575757996375115`;
- initial/evolved fundamental power: `2.1315894140899176e-6` /
  `2.1827523417862502e-6`;
- power-derived growth `sqrt(P_evolved/P_initial)`:
  `1.011929959672577`;
- expected scale-factor growth: `1.0119399775244009`.

Each run wrote
`/tmp/cosmosim_dmo_zeldovich_workflow_mpi_npN/dmo_zeldovich_power_spectrum.json`.
After normalizing only the intentionally rank-count-specific `world_size`
field, the MPI np1--np4 artifacts had the identical SHA-256 digest
`53ee518b63d21cd28790415b8f076614497536bb6a62c95276b8ed6c14984bf3`.

True non-MPI serial versus MPI-world-one evidence was run across the HDF5 and
MPI+HDF5+FFTW build trees:

```bash
cmake --build --preset build-hdf5-debug \
  --target test_validation_dmo_zeldovich_workflow --parallel 4
cmake --build --preset build-mpi-hdf5-fftw-debug \
  --target test_validation_dmo_zeldovich_workflow --parallel 4
TMPDIR=/tmp ctest --test-dir build/mpi-hdf5-fftw-debug \
  -R '^validation_dmo_zeldovich_workflow_rank_equivalence$' \
  -j4 --output-on-failure
TMPDIR=/tmp build/hdf5-debug/test_validation_dmo_zeldovich_workflow \
  --compare-serial-mpi-artifacts
```

The cross-preset comparator covered 64 stable IDs: periodic position error was
`0`, maximum velocity error was `0`, and relative mass
error was `0`. The serial and MPI1 power JSON files were byte-identical with
SHA-256
`8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`.

Focused Ewald/periodic TreePM matrix: 6/6 passed, comprising the two Ewald
validation executables plus serial and MPI np2/np3/np4 periodic TreePM. Focused
PM ownership/halo/periodic/strict-Phase-2 matrix: 12/12 passed across np1--np4.
The strict distributed TreePM thresholds remain
`relative_L2 <= 5e-6` and `max_relative <= 5e-5`; they were not relaxed.
The follow-up focused periodic/workflow matrix after cache-validity and
diagnostics hardening passed 6/6.

Ewald metrics (aggregate gate plus default-setting cutoff diagnostics):

- worst TSC relative L2: `8.059667017e-3` on the exact DMO lattice;
- worst TSC p99 normalized error: `8.157551429e-3` on that lattice;
- DMO-lattice net-force fraction and translation drift: `5.985844233e-17` /
  `1.523007584e-13`;
- DMO-lattice Ewald/PM/tree/total force L2:
  `10.16780756` / `8.296819318` / `1.789203432` / `10.08586460`;
- rectangular seam fixture relative L2/p99/net/translation drift:
  `2.003945987e-4` / `2.383921190e-3` / `4.861775315e-17` /
  `3.031142721e-14`;
- diagnostic cutoff axis/diagonal relative L2: `6.209e-3` / `8.072e-3`;
- immediately beyond cutoff axis and 1.1-cutoff diagonal relative L2:
  `6.152e-3` / `2.397e-3`;
- CIC diagnostic uniform fixture: `relative_L2=1.0368e-2`,
  `p99=7.5286e-2`;
- CIC diagnostic split-transition fixture: `relative_L2/p99=1.2138e-2`.

The certified profile reads the public/runtime default
`TreeGravityOptions::opening_theta=0.7`. A separate independent-target,
quadrupole MAC comparison calibrated `theta=sqrt(alpha)=0.07071067811865475`
for `alpha=0.005` below the `l/r<0.08` safety envelope. Geometric,
COM-distance, and relative-force results respectively reported relative L2
`4.541113485e-4`, `4.541439201e-4`, and `4.541347197e-4`; all three p99 values
were `6.721724831e-4`. Their visited/accepted/opened counts were
`48/44/4`, `56/51/5`, and `64/58/6`, with zero particle-pair work. The gate
requires fixed Ewald limits, internal accept/open behavior, less than all-pairs
work, and errors within 25 percent across MACs.

Communication packet/byte/peer counters are implemented alongside local
source/target/tree counts, global empty-rank counts, PM solve/reuse,
halo/local-slab dimensions, and tree build/multipole/opened-node counters. This
run is not a strong/weak scaling artifact, and no counter is elevated to a
scalability claim here.

## Remaining limitations

### Correctness and interface limits

- Rank-count-changing restart and arbitrary topology remap are rejected.
- Production hierarchical KDK is restricted to rung zero; per-element
  kick/drift epochs and exact mixed-rung restart are not implemented.
- Distributed open PM remains a bounded root-workspace path.
- The current MPI failure model is fail-fast; no process-failure recovery is
  implemented.

### Scientific-validation blockers

- The validated power-spectrum path is deliberately small and direct-DFT
  based; it has no interlacing, PCS deposition, survey-window treatment, or
  large-mesh performance/aliasing certification.
- Ewald certification is unsoftened and small-N; adaptive/heterogeneous
  softening requires a separate accuracy envelope.
- CIC is not certified for the current default mesh/split target.
- Isolated/open PM window deconvolution is unsupported. Focused zoom correction
  forces it off and is outside the global periodic TSC certification.
- No large DMO convergence, halo mass function/profile, long-horizon growth,
  cross-code, or publication-grade cosmology campaign has passed.

### Scalability blockers

- Tree hierarchy summaries use `MPI_Allgatherv`; target and response transport
  uses communicator-wide `MPI_Alltoallv`. This is not LET.
- PM routed records use communicator-wide count/data collectives rather than a
  validated sparse-neighbor plan.
- Zoom correction all-gathers high-resolution sources behind a hard byte cap.
- Evidence stops at four ranks.

### Performance blockers

- No current strong/weak scaling campaign quantifies tree, PM transpose,
  request imbalance, memory pressure, or communication volume beyond np4.
- TSC raises deposition/gather work from 8 to 27 cells per particle.
- The certified 6.25-cell cutoff increases short-range traversal volume, pair
  work, and distributed target/response traffic relative to the former
  4.5-cell default.
- CUDA lacks TSC parity and a certified GPU accuracy profile.

### Future contender-grade work

- Implement and validate a true locally essential tree or selective remote-node
  import protocol with bounded summaries, epoch-safe branch requests, and
  communication-volume equivalence evidence.
- Add sparse-peer PM transport while retaining the collective path as a
  correctness reference until equivalence is proven.
- Extend the deterministic 3D `P(k)` reference path with an FFT-backed
  production-scale implementation, optional interlacing/PCS, and explicit
  equivalence tests against the current direct-DFT correctness reference.
- Run large DMO resolution/scaling and science-validation campaigns before any
  publishable or external-code-parity claim.
