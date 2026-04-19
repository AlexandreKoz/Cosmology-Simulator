# TreePM Phase 1 integration closeout (feature branch)

## Scope of this closeout

This closeout is limited to **TreePM Phase 1 integration coherence** for the current gravity-upgrade feature branch.
It does not claim Phase 2 features, multi-rank TreePM production readiness, or Ewald-exact periodic references.

## Final numerical contract (implemented)

- Periodic comoving Poisson solve in PM:
  - `∇² φ(x) = 4 π G a² δρ(x)` with `δρ = ρ - ρ̄`
  - Fourier solve for non-zero mode: `φ_k = - 4 π G a² δρ_k / k²`
  - force relation: `a(k) = - i k φ_k`
  - zero-mode policy: `φ_{k=0}=0`, `a_{k=0}=0`
- TreePM Gaussian split:
  - long-range PM filter: `F_LR(k) = exp(-k² r_s²)`
  - short-range residual factor: `F_SR(r) = erfc(r/(2r_s)) + (2/sqrt(pi))*(r/(2r_s))*exp(-(r/(2r_s))²)`
  - composition target: `F_SR + F_LR = 1` (before explicit cutoff truncation)
- Split/cutoff mapping from config (authoritative in Phase 1):
  - `Δmesh = box_size / PMGRID`
  - `r_s = asmth_cells * Δmesh`
  - `r_cut = rcut_cells * Δmesh`
- Softening contract:
  - residual tree acceleration is existing softened tree law multiplied by Gaussian short-range factor (`a_SR = a_tree_softened * F_SR`)
- PM cadence semantics:
  - config control: `treepm_update_cadence_steps >= 1`
  - cadence unit: gravity kick opportunities (`gravity_kick_pre`, `gravity_kick_post`)
  - PM long-range field refreshes every `N` opportunities and is reused between refreshes

## Phase 1 hard-gate checklist

- **No hard-coded PM mesh in workflow path:** runtime PM grid comes from typed config `numerics.treepm_pm_grid`; workflow emits PMGRID in operational metadata/provenance.
- **No hidden split assumptions:** `asmth_cells` and `rcut_cells` are typed controls, normalized, mapped to `r_s`/`r_cut`, and consumed in runtime TreePM options.
- **Reproducibility evidence:** deterministic repeated-run check is wired in `test_reference_workflow_end_to_end` by executing the same config twice and requiring identical `final_state_digest` plus cadence record equality.
- **Force-error map evidence:** `bench_tree_pm_force_error_map` sweeps `(PMGRID, ASMTH, RCUT)` and writes `validation/artifacts/tree_pm_force_error_map.csv`.

## Validation commands for this stage

Required preset path for this stage:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
./build/pm-hdf5-fftw-debug/bench_tree_pm_force_error_map
```

Repeated-run reproducibility check command used in this stage:

```bash
ctest --preset test-pm-hdf5-fftw-debug -R integration_reference_workflow_end_to_end --output-on-failure
```

## Generated evidence artifact

- `validation/artifacts/tree_pm_force_error_map.csv`
  - columns: `reference_method,pm_grid,asmth_cells,rcut_cells,relative_l2_error`
  - required sweep coordinates: PMGRID `{16,24,32}`, ASMTH `{0.8,1.25,2.0}`, RCUT `{3.0,4.5,6.0}`

## Known limits retained in Phase 1

- Single-rank TreePM integration path only.
- Periodic references are explicit approximations (minimum-image direct or spectral+direct proxy), not Ewald exact.
- CUDA PM path remains CIC-only in this stage and rejects unsupported assignment settings rather than silently diverging.

## Source-tree blocker status after Stage 8 follow-up

- The earlier `gravity_pm_glass_like` validation blocker was caused by an over-aggressive
  deterministic lattice jitter that injected large aliasing into the PM solve.
- The validation now uses a tiny parity-balanced micro-jitter proxy instead of the earlier
  large sinusoidal displacement, so the source tree no longer intentionally carries that known
  failing validation case.
- Full Phase 1 closure still requires the command-backed preset path above to pass on the actual
  PM + HDF5 + FFTW environment used for review/CI.
