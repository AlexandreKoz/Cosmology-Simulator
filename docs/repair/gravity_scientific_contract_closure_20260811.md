# Gravity and scientific-workflow contract closure — 2026-08-11

This note records the implementation closure performed against the 2026-08-02
`src/gravity` adversarial audit and the 2026-08-08 gravity/workflow scientific
accuracy audit. The current repository, not the older audited snapshot, was the
source of truth. Findings already repaired or superseded before this campaign
were preserved rather than reimplemented.

## Production/scientific changes

### Authoritative gas gravity

- Owned authoritative **leaf gas cells** are now TreePM sources at
  `CellSoa.center_*_comoving` with `CellSoa.mass_code`.
- Generic gas-tagged particles are compatibility/lineage mirrors only and are
  excluded from the gravity source population, preventing gas-mass double
  counting.
- Gas cells are explicit gravity targets. Parentless cells no longer receive
  zero acceleration merely because `parent_particle_id == 0`.
- Multiple leaf cells may share one lineage parent and remain spatially distinct
  gravity sources/targets.
- Coarse cells covered by child patches are excluded from the authoritative
  leaf source set.
- Generic gas compatibility particles are no longer advanced by the generic
  collisionless drift/kick path.
- Force-cache validity now tracks particle-row, cell-row, and stable gas-cell
  identity generations, so refinement/coarsening/migration invalidates stale
  gas acceleration rows.
- The existing fixed/species comoving Plummer-softening policy is retained for
  gas. The workflow does not introduce continuously adaptive cell-size
  softening without the missing conservation terms.

### Tree and TreePM correctness

- Standalone tree build/traversal now has an explicit source-generation
  contract. Production callers can use a nonzero generation; legacy zero-token
  callers get an exact source-content plus resolved-softening fingerprint check.
- Source count, build-time topology/multipole/softening options, and 32-bit
  index limits are validated before traversal/build.
- Every MAC opens an internal node that contains the target before considering
  multipole acceptance. The invariant is enforced in both standalone tree and
  TreePM residual traversal.
- Tree construction reuses one partition workspace rather than allocating one
  range-sized scratch vector at every recursive node.
- Reused profile averages now divide cumulative interactions by cumulative
  traversed targets.
- Dead quadrupole algebra was removed from both tree implementations.
- Tree ordering rejects non-finite coordinates and guards the 32-bit particle
  index contract.
- Pairwise softened inverse-r3 uses `x * sqrt(x)` instead of a general
  `pow(x, 1.5)` call.

### PM hardening

- PM assignment selection is exhaustive: invalid enum values no longer fall
  through to TSC.
- Public PM option validation is exhaustive for assignment, boundary,
  residency, execution, and decomposition enums.
- Ordinary gravity density assignment rejects negative or non-finite source
  mass, including routed distributed contributions.
- `solveForParticles` now performs a backend-independent source/target/grid
  preflight outside the MPI-only branch, closing the single-rank CUDA validation
  hole.
- Checked cell-count/byte/launch arithmetic and zero-work CUDA launch guards
  were already present in the current repository before this campaign and were
  preserved; overflow regression coverage was added.

### TreePM cache and memory behavior

- Long-range cache validity now includes an explicit `source_generation`, kept
  separate from force/decomposition epochs.
- Periodic unwrapping reuses axis workspaces instead of repeatedly allocating
  full-population wrapped/ordered vectors.
- Traversal stack preallocation is bounded to a small initial frontier rather
  than the complete node count.
- Zoom correction no longer keeps coarse and focused PM grids alive
  simultaneously and reuses the persistent zoom-correction acceleration lanes
  for the focused solve, materially reducing peak temporary memory.
- The three periodic seam sorts remain; replacing the algorithm itself is
  performance work, not a correctness prerequisite.

### Cosmological source epoch

`SourceTerms` acts on post-drift/post-hydro state. The source runtime therefore
uses `timeline_step.scale_factor_end` explicitly for cosmological physical
conversions rather than the still-uncommitted begin-step
`IntegratorState.current_scale_factor`.

This applies to:

- star-formation density/length/threshold/free-fall inputs;
- newborn star formation-scale metadata;
- metal-diffusion physical volume, density, and mixing/filter length;
- BH/AGN physical gas-density conversion.

### Analysis epoch

`AnalysisHooks` still executes before `commitStep()`, but diagnostics are now
labelled/scheduled using the state actually being measured:

- completed step = `integrator_state.step_index + 1`;
- scale factor = `timeline_step.scale_factor_end`.

The already-correct post-commit snapshot/checkpoint metadata path was not
changed.

### BH/AGN mass, frame, and energy ledger

The BH baseline now has explicit semantics instead of a one-line Bondi patch:

- stored comoving density is converted to physical density with `a^-3` when the
  coordinate frame is cosmological;
- Bondi relative velocity is the actual BH/gas peculiar-velocity difference;
- requested accretion is capped by available host gas;
- gas mass is removed by exactly the mass transferred to the BH;
- gas density is reduced consistently at fixed cell volume;
- swallowed metal mass is removed proportionally and accounted in counters;
- BH dynamical/subgrid mass grows by the transferred mass;
- BH peculiar velocity is mass-weighted with swallowed gas, conserving the
  local linear-momentum transfer represented by this host-cell model;
- feedback energy is computed from the **actual** transferred mass/rate after
  the gas availability cap;
- `internal_energy_code` is treated as specific internal energy, so deposited
  total feedback energy is divided by remaining gas mass before addition;
- seed BH mass is also removed from authoritative host gas rather than created
  as duplicate gravitating mass;
- seeds inherit host-gas velocity;
- production seeding uses the existing particle-ID precommit registry and
  registers/schedules newborn BH particles through the normal source-runtime
  birth path.

The chosen baseline is an explicit Newtonian tracked-mass ledger: `mdot_acc`
is gas rest mass transferred into BH tracked mass; `epsilon_r` controls the
unresolved radiative-energy yield and is not an implicit deletion of tracked
mass. This is documented in `docs/black_hole_agn.md`.

## TreePM softening/split accuracy envelope

The current long-range PM filter is Newtonian while the short-range tree force
is Plummer softened. Instead of introducing a new unvalidated softened Fourier
operator, production workflow source assembly now enforces

`epsilon / r_s <= 0.20`

for the effective fixed/species/particle softening. For the implemented
Gaussian long-range plus Plummer short-range combination, the analytic
split/softening mismatch is below about 0.65% over separation at the envelope
limit; normal reference decks are far inside this bound. Continuously adaptive
gas softening remains unsupported pending a conservation-aware formulation and
validation campaign.

## Finding closure map

| Finding | Status | Current disposition |
|---|---|---|
| GRAV-001 | **Repaired** | Tree source generation + legacy content/softening identity checks bind traversal to build state. |
| GRAV-002 | **Repaired** | Target-containing internal nodes are opened before every standalone MAC. |
| GRAV-003 | **Repaired** | Common backend-independent PM preflight now covers single-rank CUDA dispatch. |
| GRAV-004 | **Already repaired in current baseline; statically preserved** | Current CUDA code already had zero-work guards and checked launch/byte arithmetic. CUDA runtime not available here. |
| GRAV-005 | **Already repaired in current baseline + regression** | `PmGridShape::cellCount()` uses checked products; overflow test added. |
| GRAV-006 | **Repaired** | Invalid assignment/options enums fail explicitly. |
| GRAV-007 | **Repaired** | Tree source/node/Morton ordering 32-bit boundaries are guarded before narrowing. |
| GRAV-008 | **Repaired** | Recursive range-sized allocation replaced by reusable partition workspace. |
| GRAV-009 | **Repaired** | Ordinary PM gravity rejects negative mass. |
| GRAV-010 | **Repaired** | Long-range validity/fingerprint includes explicit source generation. |
| GRAV-011 | **Repaired** | Standalone tree validates opening criterion, multipole order, softening, and numerical bounds. |
| GRAV-012 | **Repaired** | Cumulative target count fixes reused profile average in tree and TreePM. |
| GRAV-013 | **Repaired** | TreePM stack starts with a bounded small reserve, not O(nodes). |
| GRAV-014 | **Partially closed** | Full-size periodic scratch allocations are reused; three O(N log N) seam sorts remain. |
| GRAV-015 | **Repaired materially** | Coarse/focused PM grids no longer overlap in lifetime; focused acceleration lanes are reused. |
| GRAV-016 | **Clarified** | `child_base` is explicitly diagnostic; traversal uses `child_index` and does not assume contiguous child roots. |
| GRAV-017 | **Partially closed** | Critical MAC/self-force invariants are consistent in both paths; full standalone/TreePM numerical-kernel deduplication remains maintainability work. |
| GRAV-018 | **Repaired** | Calculated-and-discarded quadrupole terms removed from both implementations. |
| GRAV-019 | **Already repaired in current baseline** | OpenMP active-target tree traversal was present before this campaign. |
| GRAV-020 | **Remaining bounded architecture limitation** | CUDA path still stages assignment/interpolation around a host PM FFT; no speculative full-GPU rewrite was attempted without CUDA tooling. |
| GRAV-021 | **Remaining diagnostic limitation** | Non-FFTW naive DFT fallback remains tiny-grid diagnostic infrastructure, not production certification. |
| GRAV-022 | **Partially superseded** | Invariants were centralized/clarified without a cosmetic file-splitting campaign; large translation units remain. |
| GRAV-023 | **Partially closed** | Feature-path cleanup was performed where touched; a complete strict-warning matrix was not run in this environment. |
| GRAV-024 | **Superseded by current repository policy** | Current naming governance uses `kPascalCase`, `k_lower_snake`, `_comoving`, and `.cpp`; no dangerous legacy-rule mass rename performed. |
| GRAV-025 | **Already clarified in current baseline** | Per-particle softening mask is authoritative; unmasked value lane is not an implicit override. |
| GRAV-026 | **Repaired** | Hot `pow(x,1.5)` replaced by multiply/sqrt form. |
| GRAV-027 | **Bounded internal contract** | Coordinator/public PM and accumulator shape preflights remain the ownership boundary; no unsafe public production path was introduced. |
| GRAV-028 | **Obsolete in current repository** | The audited ceremonial gravity-module boundary files are absent from the current source tree. |

### Scientific workflow findings

| Finding | Status | Current disposition |
|---|---|---|
| SCI-WF-001 | **Implementation repaired** | Authoritative owned leaf gas cells source/receive gravity directly; compatibility particles no longer own gas gravity. Empirical self-gravitating hydro/AMR certification still remains. |
| SCI-WF-002 | **Repaired** | Source-stage cosmological conversions use the end-of-step evaluation epoch. |
| SCI-WF-003 | **Repaired** | Bondi density receives explicit physical/comoving frame conversion. |
| SCI-WF-004 | **Repaired at implementation level** | Accretion and seeding conservatively transfer tracked mass from host gas to BH; feedback specific-energy units corrected. |
| SCI-WF-005 | **Repaired** | Analysis diagnostics label completed state with completed step/end scale factor. |
| SCI-WF-006 | **Repaired** | Target-containing node rule is universal in standalone tree and TreePM. |
| SCI-WF-007 | **Closed by certified envelope** | Workflow enforces `epsilon/r_s <= 0.20` rather than inventing an unvalidated softened PM kernel. |
| SCI-WF-008 | **Closed by policy** | Gas uses fixed/species comoving softening; continuously adaptive cell-size softening is not silently enabled. |

## Verification actually executed

Environment capability check was performed once:

- CMake/Ninja/G++: available;
- OpenMP: available;
- HDF5 development package: available;
- MPI compiler/runtime: unavailable;
- FFTW development headers/library: unavailable;
- CUDA/nvcc/compute-sanitizer: unavailable.

CPU Debug configure succeeded with `cmake --preset cpu-only-debug`.
Affected production libraries (`cosmosim_gravity`, `cosmosim_physics`, and
`cosmosim_workflows`) built successfully after the final repairs.

The following rebuilt focused tests passed:

- `unit_tree_gravity`
- `unit_pm_solver`
- `unit_black_hole_agn`
- `unit_stage6_source_diagnostics_views`
- `integration_tree_gravity_vs_direct`
- `integration_tree_pm_coupling_periodic`
- `integration_black_hole_agn_toy`
- `integration_star_formation_source_runtime`
- `integration_reference_workflow_hydro_row_order`
- `integration_restart_equivalence_reference_workflow_hydro`

The star-formation/source-runtime integration fixture additionally runs legal
parentless-gas and shared-lineage-gas states through the repaired production
workflow.

The periodic TreePM normalization regression continues to report zero scaling
error in its tested normalization case.

A broader selected run had 10/11 tests pass. The single failure,
`integration_hydro_decoupled_gas_cells`, fails in its own fixture construction
at `validateOwnershipInvariants()` before any repaired gravity/source path is
executed. Its test source is byte-for-byte identical to the supplied baseline,
so this campaign did not mutate unrelated production ownership rules merely to
force that pre-existing fixture green.

ASan/UBSan configuration succeeded and the sanitizer build compiled through the
repaired gravity sources/libraries, but two bounded build invocations exhausted
the execution window while compiling the repository dependency graph before the
requested sanitizer test executables were available. No sanitizer test result
is therefore claimed.

The complete supported CPU Debug build subsequently finished successfully.
All configured production libraries, test/validation executables, and benchmark
targets linked without a compiler or linker error.

## Verification not executed

No claim is made for:

- multi-rank MPI runtime behavior;
- FFTW/Ewald production certification;
- CUDA execution or compute-sanitizer;
- large nonlinear DMO convergence;
- self-gravitating Evrard/Jeans/hydrostatic convergence;
- publication-grade baryonic calibration.

Those require toolchains or scientific campaigns not available in this repair
environment.

## Remaining scientific/engineering limitations

- The authoritative gas-gravity coupling is implemented, but self-gravitating
  hydro/AMR still needs empirical convergence and multi-rank validation before
  being promoted to scientific production.
- BH accretion remains a host-cell-local closure; there is no kernel-neighborhood
  accretion, BH repositioning, or BH merger model in this closure.
- Swallowed metal mass is removed/accounted but is not yet a dedicated persistent
  BH metal reservoir.
- TreePM distributed residual exchange is not a mature LET implementation.
- CUDA is not a full device-resident PM solver.
- The naive DFT fallback is not a production large-mesh backend.
- Large DMO force/observable convergence and cross-code evidence remain release
  validation tasks, not implementation claims.
