# Campaign M2C — Full-Physics Sidecar and Event-Memory Architecture

Date: 2026-09-05
Status: source/architecture closure implemented; focused CPU validation passed

## Scope and authority

This campaign starts from the accepted M2A/M2B hydro/AMR memory architecture and
closes the remaining full-physics memory defects that were visible in current
source. It does not change the configuration authority, restart schema, HDF5
schema, force/hydro methods, source-term efficiencies, feedback coupling law,
chemistry network, numerical precision, or accepted tolerances.

The current state model already has the correct primary species split:
`ParticleSoa` contains the universal particle skeleton; gas thermodynamics live
in `GasCellSidecar`; stellar evolution/feedback truth lives in
`StarParticleSidecar`; black-hole state lives in `BlackHoleParticleSidecar`; and
tracer state lives in `TracerParticleSidecar`. M2C preserves those authorities
rather than adding another sidecar registry.

## A. Species/module sidecar audit

Logical persistent widths below are sums of the declared SoA lane widths. They
are species-local costs, not universal-particle costs. Allocator alignment and
capacity slack remain visible through the existing central memory report.

| State | Applicability | Logical bytes / active row | M2C classification |
| --- | --- | ---: | --- |
| `CellSoa` + `GasCellSidecar` | gas cells only | 125 B/gas cell | conserved/restart gas truth plus the previously accepted pressure/temperature/sound-speed compatibility caches |
| `StarParticleSidecar` | stars only | 224 B/star | birth identity, SSP integration watermark, cumulative enrichment/feedback ledgers, unresolved carry authority; no runtime carry mirror after M2C |
| `BlackHoleParticleSidecar` | BHs only | 72 B/BH | compact BH model/restart/observable state; gas and non-BH particles do not allocate it |
| `TracerParticleSidecar` | tracers only | 48 B/tracer | attachment and conservative mass-exchange truth only |
| dust sidecar | not implemented | 0 B | no placeholder allocation is introduced |

At the 512^3 reference population (134,217,728 rows), the hypothetical cost of
materializing each species sidecar over the *entire* population would be 15.625
GiB for gas cell+sidecar state, 28 GiB for star sidecar state, 9 GiB for BH
sidecar state, and 6 GiB for tracer sidecar state. Production cost instead
scales with the actual count of that species. A single new 8-byte full-population
lane at this reference scale is exactly 1 GiB and therefore still requires an
explicit review justification.

Focused memory-accounting coverage now asserts that a 4096-particle DMO state
owns zero gas/star/BH/tracer cold-sidecar capacity, and a 2048-cell gas state
owns zero star/BH/tracer cold-sidecar capacity.

## B. Shared read-only microphysics tables

The effective multiphase ISM EOS table was previously constructed separately by
SourceRuntime, HydroAmrRuntime, and TimeCoordinator. The reference runtime now
constructs one immutable `shared_ptr<const EffectiveMultiphaseEosTable>` and
passes the same table to all three consumers. Standalone factories retain a
local fallback so focused tests and non-reference callers remain self-contained.

`EffectiveIsmThermodynamicClosure` now accepts shared immutable table ownership;
its value-taking compatibility constructor wraps one private shared instance.
A unit regression asserts object identity for the shared constructor. Stellar
evolution remains one table per SourceRuntime model, not per star/thread, and
no per-element cooling table was found in the current production path.

## C. Chemistry/enrichment fidelity tiers

M2C uses the existing typed scientific configuration rather than creating a
competing fidelity selector:

- `physics.cooling_model=primordial` is the released primordial cooling tier;
- `physics.cooling_model=primordial_metal_line` explicitly activates the
  configured metal-line table;
- `physics.metal_species_mode=total_only` is the released scalar-enrichment
  fidelity and stores authoritative conserved total metal mass on gas cells;
- `physics.metal_species_mode=core_elements` remains deliberately fail-closed in
  schema v2 because element-resolved sidecars, yields, transport, cooling-table
  schema, AMR exchange, and restart contracts are not yet complete;
- no detailed chemistry network is currently implemented, so M2C does not
  fabricate a decorative tier or silently change a running network under memory
  pressure.

This preserves the existing scientific/provenance choice: memory pressure may
change batch size but never the configured cooling/enrichment fidelity.

## D. Compact persistent star/BH truth

The major current star-state violation was a second set of six runtime
`double` arrays in `StellarFeedbackModuleState` mirroring unresolved carry that
was already restart/migration-authoritative in `StarParticleSidecar`. M2C makes
that runtime object allocation-free. This removes 48 B/star of persistent
runtime duplication: 6 GiB at a hypothetical 512^3 star population, or about
4.47 GiB at 100 million stars.

`stellar_age_years_last` is retained because it is the persistent SSP integration
watermark used to integrate only the unprocessed age interval; it is not merely
a display-age cache. BH state is already isolated to BH rows and no permanent
per-BH neighbor/search/scratch structure was found. Schema-visible BH rate and
duty-cycle lanes remain unchanged in this campaign.

## E. Event-driven feedback and tracer-local scratch

Stellar feedback now accepts compact `StellarFeedbackEvent` records and
SourceRuntime processes stellar-evolution evaluation, feedback deposition, and
stellar-ledger commit in batches capped at 4096 stars. The former all-active
`StellarEvolutionStepReport::budgets` materialization and three full-star delta
arrays (`returned_mass`, `returned_metals`, `feedback_energy`) therefore no
longer create population-scale source-step staging. The retained event batch is
bounded to about 128 KiB at the current ABI (`4096 * 32 B`), while the matching
stellar-evolution report is likewise bounded to one 4096-star batch.

Feedback target selection no longer allocates an O(N_gas) distance vector for
each event. A phase-resident x-sorted `uint32_t` gas-cell index is rebuilt from
owned leaf cells. Exact Euclidean nearest-k selection expands deterministically
around the star in x and terminates only when the remaining x-distance lower
bound cannot beat the current worst true neighbor. Per-event nearest-distance
and deposition-target vectors are therefore O(k), where k is the configured
neighbor count. With the default k=8, retained vector payload is approximately
704 B (`8 * (48 B distance record + 40 B target)`) rather than O(N_gas).
The index itself is 4 B per indexed owned leaf cell and keeps capacity high-water
truth. Production index rebuild and event-buffer growth are admitted through the
existing process `MemoryGovernor`; rank-local reservation failure is passed
through the existing `FailureCoordinator` before allocation so no rank can enter
later source-stage collectives after a peer failed memory admission.

Before M2C, the production feedback path also retained/constructed a candidate
index, ownership mask, and full gas-cell volume vector totaling 13 B/gas row for
feedback geometry. Feedback now uses only the 4 B/cell spatial index and an
on-demand exact AMR volume provider. At 512^3 gas cells that changes the feedback
geometry payload from 1.625 GiB to 0.5 GiB, a 1.125 GiB reduction. Metal
diffusion may still build its own volume/mask phase data when that independent
module is enabled.

Tracer source updates no longer retain a full `0..N_cell-1` cache (4 B/gas
cell) for the all-active case. The all-active contract is represented by an
empty active span. Subset updates sort/unique only the active cell indices and
perform membership lookup against that compact list, replacing a 1 B/N_gas
mask with O(N_active) scratch. Invalid active rows now fail deterministically.

No processed feedback-event history is retained. Unresolved scientifically
required budget survives only in the star sidecar carry lanes.

## F. Memory accounting and scale law

The existing central capacity-based report remains the sole persistent-state
memory authority and already reports each gas/star/BH/tracer lane independently.
M2C adds no new 8-byte full-population field. The new feedback index exposes
owned capacity and historical high-water directly; the feedback runtime carry
state reports zero owned capacity.

Key structural deltas:

| Removed/replaced residency | Before | After | Scale |
| --- | ---: | ---: | --- |
| feedback carry runtime mirror | 48 B/star | 0 | N_star |
| feedback/evolution source staging | population-scale evolution report + 24 B/star delta lanes | one <=4096-star evolution report + <=~128 KiB event batch | active stars, hard batch cap 4096 |
| feedback geometry support | 13 B/gas row | 4 B/indexed owned leaf | N_owned_leaf_gas |
| per-event neighbor-distance staging | ~48 B/candidate gas row | ~48 B * k | configured neighbor count |
| tracer all-active cell index cache | 4 B/gas row retained | 0 | eliminated |
| tracer subset activity mask | 1 B/gas row | 4 B/active unique cell | N_active |
| effective-ISM EOS table owners in reference runtime | 3 instances | 1 immutable shared instance | subsystem-scale |

## Scientific-correctness contract and focused evidence

M2C does not change the feedback budget calculation, thermal/kinetic/momentum
mode selection, stochastic decision function, target weights, deposition
conservation equations, AMR volume semantics, stellar carry semantics, tracer
mass-exchange law, cooling table contents, or EOS evaluation.

Focused regressions establish:

- indexed and direct exact nearest-target selection return identical cell order,
  weights, and displacement vectors on deterministic irregular geometry;
- event-batch deposition matches sequential event application for gas mass,
  density, internal energy, metal mass, stellar cumulative deposition, and
  aggregate deposited mass/metals/thermal energy;
- the existing stellar-feedback box integration conservation test passes;
- unresolved feedback budget remains in persistent star carry lanes when no gas
  target exists;
- repeated event processing does not grow module-state memory and does not grow
  the spatial index unless geometry capacity itself grows;
- production feedback index/event-buffer growth is governor-admitted before
  allocation and reservation failure is collectively coordinated;
- shared EOS closure points to the exact shared table object;
- disabled species allocate no cold sidecar capacity in the central memory report.

## Validation executed before artifact-first packaging

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target \
  test_unit_stellar_feedback test_unit_tracer_support \
  test_unit_effective_multiphase_ism
ctest --test-dir build/cpu-only-debug \
  -R 'unit_(stellar_feedback|tracer_support|effective_multiphase_ism)$' \
  --output-on-failure

cmake --build --preset build-cpu-debug -j4 --target \
  test_unit_memory_accounting test_unit_stellar_feedback \
  test_unit_tracer_support test_unit_effective_multiphase_ism \
  test_integration_stellar_feedback_box
ctest --test-dir build/cpu-only-debug \
  -R '(unit_memory_accounting|unit_stellar_feedback|unit_tracer_support|unit_effective_multiphase_ism|integration_stellar_feedback_box)$' \
  --output-on-failure
```

Result at the artifact-first gate: all five focused tests pass (4 unit + 1
integration). Broader validation is attempted only after the changed-files
bundle is generated and is reported separately in the patch manifest/final
handoff.

## Continuation validation after artifact-first checkpoint

The artifact-first checkpoint was intentionally generated before extended
validation. Tool execution later resumed from the same worktree, including the
post-checkpoint 4096-star source batching and governor/FaultCoordinator
preflight changes. The final handoff therefore uses the following additional
evidence.

```text
cmake --build --preset build-cpu-debug -j4
  PASS: full CPU debug build completed, 175/175 Ninja steps.

ctest --preset test-cpu-debug --output-on-failure
ctest --test-dir build/cpu-only-debug -I 76,137 --output-on-failure
ctest --test-dir build/cpu-only-debug -I 126,137 --output-on-failure
  PASS: tests 1-124 and 126-137; 136/137 CPU tests passed.
  INCOMPLETE: integration_source_package_completeness (#125) exceeded the
  command window while its independent extracted-source build was compiling;
  no assertion/compiler failure was observed before termination.

cmake --preset hdf5-debug
  PASS: HDF5 1.14.5 detected.
cmake --build build/hdf5-debug -j4 --target \
  test_integration_restart_equivalence_stochastic_sources \
  test_integration_restart_checkpoint_roundtrip \
  test_integration_snapshot_hdf5_roundtrip \
  test_integration_stellar_feedback_box \
  test_unit_stellar_feedback test_unit_memory_accounting
  PASS.
ctest --test-dir build/hdf5-debug \
  -R '^(unit_(stellar_feedback|memory_accounting)|integration_(stellar_feedback_box|snapshot_hdf5_roundtrip|restart_checkpoint_roundtrip|restart_equivalence_stochastic_sources))$' \
  --output-on-failure
  PASS: 6/6 selected HDF5/M2C tests.

cmake --preset pm-hdf5-fftw-debug
  BLOCKED: FFTW3 serial double-precision development library not found.
cmake --preset mpi-hdf5-fftw-debug
  BLOCKED: MPI C++ support not found (`mpi-cxx` missing / MPI_CXX not usable).

git diff --cached --check
  PASS.
```

Restart evidence includes exact HDF5 roundtrip of the authoritative stellar
carry/deposition/return ledgers in `StarParticleSidecar` and deterministic
restart equivalence for stochastic source evolution. Snapshot HDF5 roundtrip
also passes. This closes the available-environment restart contract for the
state moved/rematerialized by M2C; no alternate chemistry network or numerical
fidelity is selected under memory pressure.

Final available-environment verdict: M2C source/architecture and focused
runtime acceptance are achieved. MPI/FFTW-enabled distributed validation is
environment-blocked and is not claimed.

## Deferred / next-campaign handoff

- `core_elements` remains intentionally unavailable until its complete
  element-resolved scientific/restart/AMR/MPI contract exists; M2C does not
  pretend that a sparse sidecar alone constitutes a chemistry network.
- Previously accepted gas pressure/temperature/sound-speed compatibility caches
  remain a separate restart-schema migration question from M2A; M2C does not
  silently delete them.
- Global task-DAG/resource co-scheduling, memory-aware overlap of independent
  physics/analysis phases, GPU/device physics residency, and broad analysis
  streaming are out of M2C scope and belong to the next memory campaign.
