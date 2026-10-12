# CHUÍ TreePM Attempt 1 source closure — 2026-10-12

**SOURCE IMPLEMENTATION; QUALIFICATION PENDING.** Session classification:
feature implementation mode, with mandatory defect repairs first. The authoritative checkout
started clean on `main`, HEAD `ce0721268a065ccff9a5a9ba242427a994fac192`.
No archive files replaced accepted checkout source. Attempt 1 facts below come
from the supplied request and were traced to the current source; the named
report/archive and `snap_1698` data were not present as repository files.

Compilation: **NOT RUN — deferred by user**. Scientific tests: **NOT RUN —
deferred by user**. MPI tests: **NOT RUN — deferred by user**. Performance
benchmarks: **NOT RUN — deferred by user**. Long-run convergence: **NOT RUN —
deferred by user**. Commands below are instructions, not execution evidence.
No speedup, force-error qualification, completed integration or 48 GiB fit is
claimed. Static delivery checks are recorded at the end.

## Campaign disposition

| Campaign item | Exact status | Source result / remainder |
|---|---|---|
| P0 drift abort | IMPLEMENTED IN SOURCE — NOT VALIDATED | Drift view now grants velocity read/write; task declaration and validator stay strict. |
| Composed runtime regression | IMPLEMENTED IN SOURCE — NOT VALIDATED | CPU reference workflow inspects the installed built-in drift descriptor and executes two steps plus bootstrap/endpoint forces; invalid registry grants remain rejected. CLI smoke adds three steps; distributed reference workflow checks grants collectively. |
| A bounded periodic arithmetic | IMPLEMENTED IN SOURCE — NOT VALIDATED | Actual source/node bounds select one-correction arithmetic; general/tie/rounding fallbacks remain. |
| B mandatory early descent | IMPLEMENTED IN SOURCE — NOT VALIDATED | Prune before leaves; fast containing/intersecting nodes descend before COM/MAC; explicit skipped geometry counters. |
| C node-invariant geometry | IMPLEMENTED IN SOURCE — NOT VALIDATED | Offset evaluation is lazy; geometric MAC omits it. No per-node cache added without measured benefit. |
| D homogeneous full residual | IMPLEMENTED IN SOURCE — NOT VALIDATED | Bounded error-gated Hermite table reaches production direct pairs; analytic fallback and physical governor charge. |
| D multipole radial acceleration | IMPLEMENTED IN SOURCE — NOT VALIDATED | Optional shared analytic value/derivatives, full moment trace, original analytic reference retained. |
| D adaptive proxy / shared tables | PARTIAL — REMAINDER IDENTIFIED | Existing lower-bound rejection preserved. MAC remains analytic; no derivative interpolation or cross-governor global table pool. |
| E bounded leaf SIMD | IMPLEMENTED IN SOURCE — NOT VALIDATED | Runtime AVX2 distance/contribution kernels, scalar fallback, 16-source staging per worker. Radial/table evaluation remains scalar per valid lane. |
| F reuse/refit | ALREADY PRESENT — SOURCE VERIFIED | Existing source/layout/frame/membership authority retained. No competing cache or reuse rule. |
| F preprocessing and telemetry | IMPLEMENTED IN SOURCE — NOT VALIDATED | Proven canonical lanes skip floor normalization; explicit refit attempted/accepted/rejected and geometry preparation timing. |
| F broader build elimination | PARTIAL — REMAINDER IDENTIFIED | Existing preflight/source scans and unknown-generation fingerprints remain; removing them would need a stronger caller certificate. |
| G production scheduling | IMPLEMENTED IN SOURCE — NOT VALIDATED | Typed block sizes 16/32/64/128/256 and leaf size [1,1024], default 64/16; existing dynamic target scheduling and bounded stacks retained. |
| G target work-history scheduling | NOT IMPLEMENTED — REASON | No measured imbalance evidence or new ordering contract; existing aggregate spatial feedback retained without an O(N) weight array. |
| H split tuning support | IMPLEMENTED IN SOURCE — NOT VALIDATED | Bounded staged sweep axes, candidate-specific oracle, total-reference identity checks and actual split/build/memory telemetry. No default mesh change. |
| H PM workspace redesign | NOT IMPLEMENTED — REASON | No source-confirmed redundant FFT/communication buffer removed; established ownership and overlap retained. |
| I toolchain/provenance | IMPLEMENTED IN SOURCE — NOT VALIDATED | Native/LTO variants, fresh capability/provider checks, real Git plus separate archive manifest identity. |
| Launcher permissions | IMPLEMENTED IN SOURCE — NOT VALIDATED | Tracked launcher mode corrected from 100644 to 100755; archive-safe Python invocation documented. |
| Snapshot-as-IC workflow | PARTIAL — REMAINDER IDENTIFIED | Conditional explicit-unit replay deck and bridge conversion regressions; actual snapshot/header/manifest and softening must be supplied and checked. |
| J hierarchical safeguards | ALREADY PRESENT — SOURCE VERIFIED | Source/PM epochs, bootstrap and closed-boundary restart ownership retained; no source-confirmed new integrator defect found. |
| J optimized hierarchy regression | IMPLEMENTED IN SOURCE — NOT VALIDATED | Existing production dispatcher refinement/restart matrix extended to optional kernels/reuse/refit at registered np1/np2/np4. Rung-zero default unchanged. |
| Scientific/memory qualification | PARTIAL — REMAINDER IDENTIFIED | Meaningful automated cases written; all execution, large-population force references, convergence and capacity qualification deferred. |
| Historical force oracle / total reference | BLOCKED — SPECIFIC CONTRACT | Original target list/data and independently qualified total periodic reference are external evidence, absent here; they cannot be fabricated. |

## Runtime repair and coverage deficiency

`src/workflows/runtime_module_registry.cpp` granted only read access to
`particle_velocity` in `RuntimeStageViewKind::kDriftParticles`, while the actual
`time.drift_particles` descriptor in `reference_runtime_composition.cpp`
requests read/write. The hierarchical drift keeps canonical momentum p=a*u fixed by multiplying
peculiar velocities by the timeline Hubble-drag factor. The view now matches that live contract:
position read/write, velocity read/write, migration read-only. No task grant was
downgraded and no validator check was removed.

The old registry unit fixture admitted synthetic declarations, not the actual
composed built-ins. The real CLI TreePM smoke is correctly registered only with
HDF5+FFTW; CPU-only testing cannot exercise that backend. This establishes a
coverage deficiency, not evidence of which Attempt 1 tests were actually run.
The unconditional `integration_reference_workflow` now inspects built-in
descriptors via the additive read-only `RuntimeModuleRegistry::descriptors()`
view and uses always-active checks for two completed steps and at least three
gravity refreshes. `run()` uses the real module admission/dispatcher validation
path. The span is invalidated by registry mutation/destruction; it confers no
mutation authority. Negative registration tests reject migration writes,
unauthorized scheduler access and invalid diagnostic reads. Existing official
test-only `-UNDEBUG` remains intact. MPI checks the actual drift declaration on
all participating ranks; CLI smoke requires three completed steps and a gravity
event beyond IC readiness.

## Geometry, traversal and arithmetic contracts

Production owner: `TreePmCoordinator` in `src/gravity/tree_pm_coupling.cpp`.
Public options/diagnostics: `include/cosmosim/gravity/tree_pm_coupling.hpp`.
Internal arithmetic: `src/gravity/internal/tree_pm_geometry.hpp` and
`tree_pm_acceptance.hpp`. Tree build/refit/moments remain owned by the existing
`TreeGravitySolver`; PM ownership and MPI protocols are unchanged.

Previously repeated nearest-image arithmetic was
`delta - L*nearbyint(delta/L)`. Preflight now folds original-source min/max into
its existing finite-state scan. After a build/refit, a small geometry pass finds
actual center AND COM bounds, including zero-mass nodes. For each target/axis,
`min-target >= -L && max-target <= L` certifies at most one correction. Source
and node frames receive separate certificates. Nodes may be unwrapped;
canonical particle coordinates alone never certify node deltas. Noncanonical
targets/sources and multi-box displacements fall back. Rectangular axes use
their own lengths. A narrow half-box rounding band uses original nearbyint,
preserving half ties and adjacent rounded ratios. Non-nearest rounding also
uses the general helper. Equivalence is nearest-image/force arithmetic, not a
promise about signed-zero bit patterns. Imported targets are certified against
actual receiving-node bounds; other MPI packet/domain geometry remains general.

Pruning uses conservative node AABB minimum distance first. Leaves proceed
directly to pairs. Internal target-containing nodes or nodes not wholly inside
cutoff descend before COM only in fast diagnostics. Remaining guards compute
COM distance as needed; center-COM offset is a lazy norm in the common unwrapped
frame, omitted if an earlier guard settles the decision or the selected MAC
does not require it. No node-offset vector is retained. Full-forensic mode
retains overlapping evaluated failures. Fast mode reports skipped geometry and
MAC explicitly; a skipped test is never a pass or a rejection. Inclusive pair
cutoff and deterministic octant/source order remain unchanged.

## Numerical kernel design

Let `x=r^2`, effective pair epsilon follow the existing combination rule, and
`a_s` denote Gaussian split length, distinct from cosmological scale factor:

```text
f(x) = (x+epsilon^2)^(-3/2) - L(r,a_s)/r^3
L = erf(q) - 2*q*exp(-q^2)/sqrt(pi), q=r/(2*a_s)
H(t)=L(sqrt(t))/t^(3/2), t=x/(4*a_s^2)
Gaussian coefficient = H(t)/(8*a_s^3)
```

The new optional shared analytic kernel simultaneously evaluates H, H' and H''
using the established small-t polynomial or one erf/exp pair, yielding f,
`2*df/dx` and `4*d²f/dx²`. Accepted multipoles keep all original contractions,
including `second_moment_trace`. A direct coefficient never substitutes for
these derivatives. The original radial implementation remains selectable and
handles nonfinite shared coefficients. The adaptive residual proxy retains
the original analytic expression and its safe Newtonian lower-bound rejection;
no interpolated acceptance estimate can make opening easier.

The full-residual table uses up to 16384 geometric cells in
`s=x+epsilon^2`, with 512 mantissa bins per binary exponent. Bit indexing avoids
per-pair log/sqrt for accepted table evaluations. Cubic Hermite data use analytic
value and derivative endpoints. The fourth-derivative envelope combines
`(945/16)*s^(-11/2)` for Plummer and the Gaussian integral
`[1/(2*sqrt(pi)*a_s^3)]*integral_0^1 u^2 exp(-x*u^2/(4*a_s^2)) du`.
The latter's fourth derivative is bounded by its u^10 integral, or a tighter
unscreened bound for x>0. Each cell retains an inflated absolute error envelope
including floating roundoff. Evaluation is accepted only when that envelope
is <= `5e-11*abs(interpolated coefficient)` and finite. This is a qualification
target pending tests, not a measured force guarantee. Near cancellation,
unsupported scales, zero softening, initial uncovered core and out-of-domain
values fall back analytically. Exact coincidence contributes zero vector and
still counts the physical pair. The cutoff test precedes coefficient use.

Uniform resolved SOURCE softening is mandatory; target epsilon must equal that
source epsilon for table selection. Heterogeneous resolved lanes never use the
homogeneous table. SIMD staging also requires uniform source epsilon, combines
the target epsilon correctly, preserves identity-based self exclusion, pads
tails with finite values, applies exact cutoff masks and accumulates valid
contributions in source order. AVX2 performs real four-wide distance and
contribution arithmetic with exact sqrt/div where Plummer is needed; radial
transcendentals/table error gates remain scalar. Generic builds runtime-dispatch
AVX2 only on supported x86 GCC/Clang hosts and retain portable scalar execution.
No approximate reciprocal, explicit FMA or fast-math law is introduced.

## Memory ownership and live set

| Added storage | Owner/lifetime and growth | Physical charge / overlap | Refresh and fallback |
|---|---|---|---|
| Full residual table | Coordinator, retained optional cache; allocate for eligible key | `sizeof(Table)`: 655360 cell bytes + metadata, per rank; reserve before allocation; old released before replacement | Exact epsilon/split/cutoff key; disabled/heterogeneous/key change invalidates; error gate or denied admission uses analytic law |
| Leaf scratch | Coordinator, optional worker cache; grows with planned worker count | `capacity_workers*sizeof(TreePmLeafScratch)`, including alignment; eleven 16-double lanes plus validity, normally 1472 bytes/worker; reserve first, release before growth | Refilled every block; source state is never retained; disabled/heterogeneous/denied admission uses indexed scalar path; no-AVX host uses scalar staged arithmetic |
| Source/node bounds | Inline coordinator metadata, entire coordinator lifetime | 96 bytes total, included in fixed-feedback accounting; no allocation replacement | Source preflight and node build/refit; retained on identical reuse; failed certification uses general nearest image |
| Existing Gaussian table | Existing coordinator cache | Existing `sizeof` charge, separately reported; unchanged | Existing dimensionless table contract retained |
| Per-node offset / packed source lanes | None added | 0 bytes | Lazy calculation/existing indexed source lanes |

Memory reports identify both optional caches as governed commitments. Production
serial memory estimation can evict them before mandatory task admission via
`trimOptionalAccelerationCaches`; no solve may run concurrently with eviction.
Coordinator destruction releases physical storage before its reservation.
No global table pool duplicates or obscures rank-local governor charges.
Cross-instance table sharing is deliberately deferred until lifetime/accounting
can be represented by an established owner. Optional caches do not inflate
mandatory phase estimates; remaining headroom controls admission, and their
physical capacities remain in reports. Existing source staging, tree SoA,
communication arena, PM and output owners retain their own accounting.

Projected size arithmetic (capacity=size only; not a full-process fit claim):

| Particle population | Existing 24-byte periodic XYZ staging | Any additional 8 / 16 / 32 bytes per particle |
|---|---:|---:|
| 128³ | 0.046875 GiB | 0.015625 / 0.03125 / 0.0625 GiB |
| 256³ | 0.375 GiB | 0.125 / 0.25 / 0.5 GiB |
| 512³ | 3 GiB | 1 / 2 / 4 GiB |

Measure tree vector capacities (including moment/softening/order lanes), owned
particle/sidecar capacities and source staging separately. Add each rank's
actual PM real/complex/padded fields, FFTW allocation/workspace, slab halos,
interpolation buffers, communication arenas and rank replication. Measure
simultaneous build/refit/replacement live sets and snapshot/restart staging;
do not sum unrelated phase peaks or scale one small RSS measurement by N.
Report both governor capacity and process/node high-water evidence. **48 GiB:
PROJECTED COMPONENT COSTS ONLY — TOTAL FIT UNPROVEN**, for every large candidate.

## Production configuration, scheduling and reuse

Five new typed/normalized keys default to false/false/false/64/16. Production
passes them to `TreePmOptions` and gravity estimates; no benchmark-only block
setting remains. Dynamic target blocks, fixed-order block diagnostic reduction,
bounded worker stacks and MPI owner authority are preserved. Block changes can
change diagnostic rounding grouping; per-target pair summation order is kept.
No persistent per-particle scheduling weights were added.

Identical-source reuse still requires authoritative source generation/frame/
options; moving positions must advance generation. Refit still requires dense
layout identity, compatible source/frame/options, uniform softening and every
source inside its original leaf cell. Failed certification rebuilds. The new
canonical-coordinate shortcut removes floor/div only with whole-axis bounds
proof. Node bounds refresh after accepted refit. Existing source fingerprints
and validation scans are retained where caller authority is insufficient.
Profile fields distinguish attempted/accepted/rejected refit, reason codes,
source preprocessing and geometry preparation. Historical `reuse_accepted`
includes refit; inspect explicit refit fields and `rebuild_reason` (0 disabled; 1 no valid tree; 2 source identity; 3 frame/owner;
4 build options/softening; 5 exact reuse; 6 refit; 7 leaf escape; 8 layout). Do not
interpret an unchanged benchmark snapshot as evidence of reuse in evolution.

The typed generator writes `reference`, `optimized_matched`, isolated
`matched_gaussian/full_residual/simd/radial/reuse/refit`, `adaptive_025` and
`adaptive_050_experimental` decks. Physics comes from the supplied validated
base, acceptance is matched except for explicitly adaptive decks, and global
KDK is enforced. Runtime's emitted `normalized_config.param.txt` and hash are
the actual effective input; preserve them beside event/provenance logs.
Requested switches and eligible target counters are distinct: lookup target
counts do not claim per-pair hits; SIMD target counts do not measure vector
occupancy. Worker/block telemetry has no per-interaction clock or new atomics.

## Build, source identity and checkpoint compatibility

`CMakeLists.txt` retains existing C++20, MPI, HDF5, FFTW and IPO owners, refreshes
capability checks, verifies parallel headers/MPIO linkage, diagnoses known
Conda/system HDF5 mixtures and checks native header/runtime version agreement.
Cross-compiling explicitly defers execution. Optional native tuning and six
preset variants preserve portable Release. Compiler/version/target/LTO/provider
metadata reach generated build information and runtime provenance. See
[Linux commands](../build_instructions.md).

`cmake/cosmosim_source_identity.cmake` hashes a sorted versioned manifest of
selected repository source/config/doc/test/build/launcher paths. Real Git HEAD
and clean/modified state remain separate. Archives get an explicitly labeled
manifest identity with unavailable Git SHA; a file hash never becomes a commit.
Reconfigure after source/Git changes; build/source identity does not replace
the normalized scientific configuration identity. The launcher's tracked mode
was actually 100644 despite its POSIX shebang; repaired to 100755. Python
invocation remains usable after archive extraction loses Unix mode bits.

Optional fallback and CPU dispatch can change floating roundoff within their
qualification budgets; eligible-mode counters expose which path was selected.
Rank reduction policy and per-target source order are retained. A normalized
flag alone does not promise bitwise equality across memory/CPU variants.

New normalized keys change hashes. The existing verified-v23/v7 closed-global-
KDK migration accepts only exact predecessor default groups and authoritative
roundtrip equality, records `treepm_normalized_defaults_v2`, and retains v1/v2
receipts on exact resumes. Nondefault settings do not gain compatibility by
omitting a key. Snapshot/restart/provenance schemas remain v6/v23/v7; source
identity/provider additions use existing derived provenance text.

## Evolved snapshot as external IC

`configs/treepm_evolved_snapshot_replay.param.txt` is a conditional four-step
first-light template at a=0.76366757582395572, z=0.30947028741013982. It assumes
AREPO-format Mpc/h coordinates, Msun/h masses, km/s/sqrt(a) stored velocities,
uniform 15 kpc comoving Plummer softening and the stated cosmology/box. Those
assumptions must match the actual snapshot header, Units and complete-set
manifest. The dataset was unavailable here; ingestion count alone proves none
of these physical assumptions. CHUI-native snapshots require their actual
native convention, not this interoperable bridge declaration. Softening
sidecars are not reconstructed by declaring a uniform epsilon.

`ic_file_set_manifest.cpp`, `ic_conversion_catalog.cpp`, `ic_manifest.cpp` and
`ic_stream_ingestion.cpp` own conversion. Per field:
`stored * base_SI * h^h_exp * a^a_exp * frame_factor * velocity_factor / target_SI`.
With the template's Mpc/Msun/km_s target units, x_code=x_stored/h,
m_code=m_stored/h and v_peculiar=v_stored*sqrt(a). A monofonIC mass declaration
of 10^10 Msun/h would introduce a factor 10^10 and is not interchangeable.
G_code derives from G_SI and target length/mass/time units; force kernels return
the established unscaled comoving Tree+PM quantity, whose cosmological factors
are applied by the existing integration contract. Config epsilon uses target
comoving length units and the existing pair-combination law, not a guessed
snapshot convention. Unit tests now assert mass/coordinate h conversions and
velocity conversion at the evolved epoch against the snapshot conversion path.

External IC admission establishes new synchronized integration/force state.
The legitimate initial force bootstrap remains; unavailable previous total
force history triggers existing geometric MAC fallback. A same-snapshot strict
bootstrap is separately labeled benchmark history, not invented historical
acceleration. Exact restart additionally restores validated integrator,
scheduler and force/PM cache epochs at a closed boundary. Snapshot replay,
frozen force evaluation and checkpoint continuation are different experiments.

## Tests written, not run

| Source | Independent reference / invariant |
|---|---|
| `tests/unit/test_runtime_module_registry.cpp` | Unauthorized grants rejected by actual registration, including drift migration writes. |
| `tests/integration/test_reference_workflow_end_to_end.cpp` | Actual built-in descriptor + real admission/dispatch, two completed steps and bootstrap/endpoint force refresh count. |
| `tests/integration/test_reference_workflow_distributed_treepm_mpi.cpp` | Same actual drift permissions on all ranks and existing force/workflow invariants. |
| `tests/integration/test_runtime_app_smoke.cmake.in` | Real CLI reaches gravity and completes three steps. |
| `tests/unit/test_tree_pm_split_kernel.cpp` | Long-double independent Gaussian derivative integral, Plummer derivatives, residual table envelope, scalar/AVX2 lanes, half ties/multi-wrap/rectangular frame certification. |
| `tests/integration/test_tree_pm_coupling_periodic.cpp` | Analytic reference vs staged/table kernels, coincident IDs, tails, seam/quadrupole shared derivatives, fast/full accepted-work equality, governed denial/admission/eviction, generation/layout refit invalidation, MPI force comparison. |
| `tests/unit/test_tree_gravity.cpp` | Production block memory estimate independently checks ceil(targets/block) capacity for partial blocks. |
| `tests/unit/test_units_cosmology_provenance.cpp` | Explicit caller Git identity retained and generated Git/manifest identity reaches emitted provenance. |
| `tests/unit/test_config_parser.cpp` | All production blocks/flags/leaf roundtrip, invalid inputs and normalized identity. |
| `tests/integration/test_treepm_experiment_configs.cmake.in` | Ten typed decks, separate adaptive angles, matched diagnostics and overwrite refusal. |
| `tests/unit/test_ic_reader.cpp` | Evolved velocity conversion, Msun/h vs 10^10 Msun/h and Mpc/h normalization. |
| `tests/integration/test_hierarchical_treepm_workflow.cpp` | Production hierarchical dispatcher optional-kernel phase/restart equivalence; existing timestep refinement/rank matrix retained. |

No existing tolerance was weakened. Large-state force/timestep convergence,
extreme-scale multipole qualification, independent nonlinear periodic reference,
collective allocation fault injection for the new optional paths and exhaustive
toolchain/provider variants remain required external qualification.

## Deferred manual commands and acceptance gates

### A — First light

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R 'unit_runtime_module_registry|integration_reference_workflow$|integration_treepm_experiment_configs'
cmake --preset pm-hdf5-fftw-debug -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R 'runtime_app_smoke|integration_reference_workflow$'
```

Require actual force refresh and multiple completed steps, not only
`[CHUI][IC] status=ready`. Keep exact stdout/stderr and normalized config.

### B — Coefficients, force and MPI invariants

```bash
ctest --preset test-cpu-debug --output-on-failure -R 'unit_tree_pm_split_kernel|unit_tree_gravity|unit_units_cosmology_provenance|unit_config_parser|unit_ic_reader|integration_tree_pm_coupling_periodic'
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R 'validation_tree_pm_recovery_ewald_accuracy|hierarchical_treepm_workflow|unit_restart_checkpoint_schema|restart_equivalence'
cmake --preset mpi-hdf5-fftw-debug -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R 'reference_workflow_distributed_treepm|tree_pm_coupling_periodic|hierarchical_treepm_workflow|runtime_app_mpi_treepm|distributed_ic_reader'
ctest --preset test-cpu-debug --output-on-failure
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure
```

Use the coherent providers from the build instructions. Inventory np1/2/3/4/8
force tests and np1/2/4 hierarchical tests; unavailable backends are blockers,
not passing/skipped qualification. Repeat relevant force cases with 1/2/4
OpenMP workers and compare exact discrete work plus bounded vector error.

### C — Frozen-state timing

After verifying snapshot identity/units/softening, set CHUI_SNAPSHOT to the
complete snapshot or marker, CHUI_ORIGINAL_IDS to the original 4096-ID oracle
list, CHUI_BENCH to the qualified Release benchmark, and CHUI_G_CODE,
CHUI_BOX_X/Y/Z to exact solver-unit provenance values. Set CHUI_EXPERIMENT to a
new external output directory. These are required local evidence, not supplied
data. Modern v6 input needs its completion marker; external legacy input needs
explicit `--legacy-dialect` and complete unit declarations supported by export.

```bash
mkdir -p "$CHUI_EXPERIMENT"
python3 tools/treepm_sweep.py export "$CHUI_SNAPSHOT" --group PartType1 \
  --output "$CHUI_EXPERIMENT/particles.csv" --target-ids "$CHUI_EXPERIMENT/targets.txt" \
  --fixed-target-ids "$CHUI_ORIGINAL_IDS" --targets 4096 --epsilon 0.015
wc -l "$CHUI_EXPERIMENT/particles.csv" "$CHUI_EXPERIMENT/targets.txt"
python3 tools/treepm_sweep.py run --exe "$CHUI_BENCH" --output "$CHUI_EXPERIMENT/kernels" \
  --particles "$CHUI_EXPERIMENT/particles.csv" --target-ids "$CHUI_EXPERIMENT/targets.txt" \
  --box-x "$CHUI_BOX_X" --box-y "$CHUI_BOX_Y" --box-z "$CHUI_BOX_Z" --g "$CHUI_G_CODE" \
  --scale-factor 0.76366757582395572 --epsilon 0.015 --meshes 64 --leaves 16 --blocks 64 \
  --policies strict --kernels analytic,lookup,full_residual,simd,simd_full_residual,radial \
  --accounting full --threads 4 --repeats 5 --warmups 1 --oracle
```

Require 262145 CSV lines (header + all sources) and 4096 target-ID lines.
Verify metadata JSON before launch; fixed target sampling applies to oracle
evaluation only, while all sources/targets are timed. If epsilon or units differ
from the conditional template, use verified values consistently everywhere.
Repeat analytic alone with fast accounting, then repeat accepted kernel(s) with
fast accounting. Run adaptive separately with `--policies adaptive --history
bootstrap --maximum-angle 0.25`, and a new output directory for 0.5. Do not merge
the angle comparisons. `--reuse on` measures identical frozen-source reuse;
`--refit on --reuse off` measures only the frozen legal-refit experiment.

### D — Split/leaf/block tuning

Reuse the exact Stage C source/target/box/G/epoch/epsilon arguments in each
command below (shown as a Bash array to avoid silently falling back to synthetic
sources). Select the qualified acceptance/kernel first.

```bash
chui_frozen_args=(--particles "$CHUI_EXPERIMENT/particles.csv" --target-ids "$CHUI_EXPERIMENT/targets.txt" --box-x "$CHUI_BOX_X" --box-y "$CHUI_BOX_Y" --box-z "$CHUI_BOX_Z" --g "$CHUI_G_CODE" --scale-factor 0.76366757582395572 --epsilon 0.015 --threads 4 --oracle)
python3 tools/treepm_sweep.py run --exe "$CHUI_BENCH" --output "$CHUI_EXPERIMENT/meshes" "${chui_frozen_args[@]}" --meshes 64,96,128 --leaves 16 --blocks 64
python3 tools/treepm_sweep.py run --exe "$CHUI_BENCH" --output "$CHUI_EXPERIMENT/leaves" "${chui_frozen_args[@]}" --meshes 64 --leaves 4,8,16,32 --blocks 64
python3 tools/treepm_sweep.py run --exe "$CHUI_BENCH" --output "$CHUI_EXPERIMENT/blocks" "${chui_frozen_args[@]}" --meshes 64 --leaves 16 --blocks 16,32,64,128,256
```

Substitute the winning qualified mesh in later stages. Optional 192 must first
pass memory admission/projection. `--asmths`/`--rcuts` support small justified
split studies; default candidate cap is 40, rather than a blind Cartesian sweep.
At fixed cell-valued split, changing mesh changes physical split/cutoff. Each
candidate regenerates its own matching direct short oracle. For total force,
add `--total-reference PATH --total-reference-manifest PATH` with independently
converged periodic evidence bound to the same frozen state; short-force error
alone cannot qualify PM+Tree. Small uniform fixtures can use `--ewald-level`
and require convergence with increasing Ewald level. Record wall/PM/traversal/
build/refit/preprocessing, visits/pairs/multipoles, actual split/epsilon,
worker imbalance, capacities/RSS, compiler/CPU and normalized effective settings.

### E — Actual evolved-state replay

```bash
build/mpi-hdf5-fftw-release/cosmosim_treepm_experiment_configs \
  configs/treepm_evolved_snapshot_replay.param.txt /tmp/chui_treepm_decks_20261012
python3 ./chui run /tmp/chui_treepm_decks_20261012/reference.param.txt --preset mpi-hdf5-fftw-release
python3 ./chui run /tmp/chui_treepm_decks_20261012/optimized_matched.param.txt --preset mpi-hdf5-fftw-release
python3 ./chui run /tmp/chui_treepm_decks_20261012/adaptive_025.param.txt --preset mpi-hdf5-fftw-release
```

Fix the base's actual path/header/cosmology/unit/softening declarations before
generation. The generator refuses an existing directory and makes input paths
absolute. Inspect emitted normalized config/hashes, completed steps, finite
phase space/accelerations, established health/conservation diagnostics, source
epochs, reuse/refit reasons, active optimization counters and memory high-water.
Compare kernel variants separately. Qualify timestep refinement and trajectory
error before increasing the four-step cap or running to z=0. Adaptive 0.5 is a
separate experimental run after its own force/trajectory qualification. Exact
checkpoint restart uses the existing restart path and its equivalence tests.

### F — Memory and release acceptance

Collect actual `MemoryReport` capacities and process/node high-water snapshots
at source preparation, tree build/refit, PM/FFTW, exchange and output/restart
boundaries for 128³, 256³ and 512³. Apply the live-set accounting above. Exercise
optional admission denial, cache eviction and mandatory collective-failure
paths. Label each fit as measured, projected or exceeded. A successful small
run or table allocation does not prove 48 GiB production fit.

## Remaining uncertainty and delivery

All new accelerators remain opt-in; full diagnostics, strict acceptance,
block64, leaf16, mesh64 and global KDK defaults are retained. Geometry/default
code changes still need compilation and equivalence validation. Table envelope
roundoff, extreme scales, accepted multipole trace, SIMD CPU/compiler variants,
MPI failure coordination, restart migration and nonlinear convergence need the
listed gates. No table-hit sampling, cross-instance table pooling, persistent
packed particle layout, PM redesign or new integrator physics is claimed.

Working-tree files are authoritative. The source-only patch ZIP contains only
the exact changed/new paths below, preserving relative paths and launcher mode;
no build trees, data, binaries, caches, local presets or runtime outputs.
Suggested branch: `perf/treepm-hotpath-memory-bounded-recovery` (the original
branch was preserved). Suggested title: **CHUÍ TreePM Deep Performance Recovery:
Runtime Contract Repair, Memory-Bounded SIMD, Fast Geometry and Scientific
Qualification**.

## Static delivery evidence

- `git diff --check HEAD`: passed (whitespace/static patch check).
- `./scripts/ci/check_repo_hygiene.sh`: could not execute, exit 126 because the existing script lacks its executable bit.
- `bash ./scripts/ci/check_repo_hygiene.sh`: passed. This performs static filename/artifact/config-alias/preset/shell checks; its CMake script/list-presets calls do not configure or compile.
- Inline Python `ast.parse(tools/treepm_sweep.py)` and `json.loads(CMakePresets.json)`: passed grammar/JSON inspection only; no module or test execution.
- Source-only ZIP uses the exact manifest below and is reviewed by member bytes and Unix modes, without executing its contents.

No CMake configuration, build, CTest, Python test, MPI, simulation, benchmark,
dependency installation or hardware profiling was executed.

## Exact modified/new file manifest

<!-- CHUI_ATTEMPT1_MANIFEST -->

56 files: 50 tracked modifications (including launcher mode), 6 new files.

| Repository-relative path | Kind | Change |
|---|---|---|
| [CMakeLists.txt](../../CMakeLists.txt) | Modified | Capability/provider checks, source identity, native/IPO target scope, generator and CTest registration. |
| [CMakePresets.json](../../CMakePresets.json) | Modified | Portable, native, LTO and combined Release variants. |
| [CURRENT_STATUS.md](../../CURRENT_STATUS.md) | Modified | Current source-only status and qualification boundary. |
| [README.md](../../README.md) | Modified | User-facing recovery, experiment and launcher workflow. |
| [bench/bench_tree_pm_sweep.cpp](../../bench/bench_tree_pm_sweep.cpp) | Modified | Production kernel/reuse/refit controls and effective timing/build telemetry. |
| [chui](../../chui) | Modified | Tracked executable mode only: 100644 to 100755. |
| [cmake/cosmosim_feature_summary.cmake](../../cmake/cosmosim_feature_summary.cmake) | Modified | Generated source/provider/native/LTO metadata. |
| [cmake/cosmosim_source_identity.cmake](../../cmake/cosmosim_source_identity.cmake) | New | Real Git identity and independent sorted archive manifest hashes. |
| [cmake/templates/build_config.hpp.in](../../cmake/templates/build_config.hpp.in) | Modified | Compile-time source/compiler/provider/native metadata. |
| [configs/treepm_evolved_snapshot_replay.param.txt](../../configs/treepm_evolved_snapshot_replay.param.txt) | New | Conditional explicit-unit four-step evolved external-IC deck. |
| [docs/architecture/runtime_truth_map.md](../../docs/architecture/runtime_truth_map.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/build_instructions.md](../../docs/build_instructions.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/configuration.md](../../docs/configuration.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/memory_governance.md](../../docs/memory_governance.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/output_schema.md](../../docs/output_schema.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/profiling.md](../../docs/profiling.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/repair/treepm_deep_hotpath_attempt1_closure_20261012.md](../../docs/repair/treepm_deep_hotpath_attempt1_closure_20261012.md) | New | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/repair_open_issues.md](../../docs/repair_open_issues.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/repair_state_recap.md](../../docs/repair_state_recap.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/restart_checkpointing.md](../../docs/restart_checkpointing.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/tree_pm_coupling.md](../../docs/tree_pm_coupling.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [docs/validation_plan.md](../../docs/validation_plan.md) | Modified | Document matching behavior, ownership, qualification and manual workflow. |
| [include/cosmosim/core/config.hpp](../../include/cosmosim/core/config.hpp) | Modified | Typed conservative arithmetic/block/leaf options. |
| [include/cosmosim/gravity/gravity_memory.hpp](../../include/cosmosim/gravity/gravity_memory.hpp) | Modified | Authoritative block size for memory estimates. |
| [include/cosmosim/gravity/tree_gravity.hpp](../../include/cosmosim/gravity/tree_gravity.hpp) | Modified | Explicit refit and geometry preparation profile fields. |
| [include/cosmosim/gravity/tree_pm_coupling.hpp](../../include/cosmosim/gravity/tree_pm_coupling.hpp) | Modified | Options, bounded scratch, derived-cache ownership and diagnostics. |
| [include/cosmosim/gravity/tree_pm_split_kernel.hpp](../../include/cosmosim/gravity/tree_pm_split_kernel.hpp) | Modified | Shared analytic radial derivatives and error-gated full residual table. |
| [include/cosmosim/io/restart_checkpoint.hpp](../../include/cosmosim/io/restart_checkpoint.hpp) | Modified | Document additive v2 compatibility receipt behavior. |
| [include/cosmosim/workflows/gravity_runtime.hpp](../../include/cosmosim/workflows/gravity_runtime.hpp) | Modified | Serial memory-estimation cache eviction contract. |
| [include/cosmosim/workflows/runtime_module_registry.hpp](../../include/cosmosim/workflows/runtime_module_registry.hpp) | Modified | Read-only installed descriptor inspection interface. |
| [src/core/config.cpp](../../src/core/config.cpp) | Modified | Register, validate, parse and normalize all new production keys. |
| [src/core/provenance.cpp](../../src/core/provenance.cpp) | Modified | Real configured Git identity and native/LTO feature provenance. |
| [src/core/version.cpp](../../src/core/version.cpp) | Modified | Build identity includes source/compiler/provider metadata. |
| [src/gravity/gravity_memory.cpp](../../src/gravity/gravity_memory.cpp) | Modified | Actual block diagnostic estimate and inline geometry metadata charge. |
| [src/gravity/internal/tree_pm_acceptance.hpp](../../src/gravity/internal/tree_pm_acceptance.hpp) | Modified | Lazy COM-offset evaluation with analytic MAC retained. |
| [src/gravity/internal/tree_pm_geometry.hpp](../../src/gravity/internal/tree_pm_geometry.hpp) | Modified | Bounds-certified one-correction periodic metric with general/tie fallback. |
| [src/gravity/internal/tree_pm_leaf_kernel.hpp](../../src/gravity/internal/tree_pm_leaf_kernel.hpp) | New | Explicit runtime-selected AVX2 kernels and portable scalar arithmetic. |
| [src/gravity/tree_pm_coupling.cpp](../../src/gravity/tree_pm_coupling.cpp) | Modified | Production geometry, traversal, residual table, bounded staging and cache lifecycle. |
| [src/io/restart_checkpoint.cpp](../../src/io/restart_checkpoint.cpp) | Modified | Verified exact default-group migration and retained v1/v2 receipts. |
| [src/workflows/gravity_runtime.cpp](../../src/workflows/gravity_runtime.cpp) | Modified | Wire production flags/block/leaf, evict optional caches and emit telemetry. |
| [src/workflows/output_restart_runtime.cpp](../../src/workflows/output_restart_runtime.cpp) | Modified | Attach source-manifest/provider identity to existing provenance. |
| [src/workflows/runtime_module_registry.cpp](../../src/workflows/runtime_module_registry.cpp) | Modified | Repair drift velocity grant and expose immutable descriptors. |
| [tests/integration/test_hierarchical_treepm_workflow.cpp](../../tests/integration/test_hierarchical_treepm_workflow.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/integration/test_reference_workflow_distributed_treepm_mpi.cpp](../../tests/integration/test_reference_workflow_distributed_treepm_mpi.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/integration/test_reference_workflow_end_to_end.cpp](../../tests/integration/test_reference_workflow_end_to_end.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/integration/test_runtime_app_smoke.cmake.in](../../tests/integration/test_runtime_app_smoke.cmake.in) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/integration/test_tree_pm_coupling_periodic.cpp](../../tests/integration/test_tree_pm_coupling_periodic.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/integration/test_treepm_experiment_configs.cmake.in](../../tests/integration/test_treepm_experiment_configs.cmake.in) | New | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/unit/test_config_parser.cpp](../../tests/unit/test_config_parser.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/unit/test_ic_reader.cpp](../../tests/unit/test_ic_reader.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/unit/test_runtime_module_registry.cpp](../../tests/unit/test_runtime_module_registry.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/unit/test_tree_gravity.cpp](../../tests/unit/test_tree_gravity.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/unit/test_tree_pm_split_kernel.cpp](../../tests/unit/test_tree_pm_split_kernel.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tests/unit/test_units_cosmology_provenance.cpp](../../tests/unit/test_units_cosmology_provenance.cpp) | Modified | Add source regression coverage for the corresponding production invariant; not run. |
| [tools/treepm_experiment_configs.cpp](../../tools/treepm_experiment_configs.cpp) | New | Typed normalized isolated experiment generator. |
| [tools/treepm_sweep.py](../../tools/treepm_sweep.py) | Modified | Staged split-aware experiments and effective metadata. |
