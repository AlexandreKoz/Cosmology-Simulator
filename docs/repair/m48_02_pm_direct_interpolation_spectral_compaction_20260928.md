# M48-02 — Direct Indexed PM Interpolation and Spectral Workspace Compaction

## Closure

M48-02 removes avoidable population- and spectral-volume ownership from the periodic PM/TreePM path without changing the PM mesh, precision, assignment scheme, deconvolution, TreePM split, or force convention.

### Direct indexed interpolation

`PmSolver::interpolateForces(...)` no longer gathers indexed target coordinates into three compact XYZ vectors. Indexed targets retain their validated `coordinate_source_index` and borrow `pos_x/pos_y/pos_z` directly. The coordinate-layout choice is invariant for one interpolation call; each target resolves its source row before stencil construction, while compact-coordinate callers continue to use their already-contiguous lanes. Indexed-global acceleration output staging remains unchanged because it has a separate scatter-ownership role.

### Spectral ownership

Periodic `PlanResources` now retains only the large plan-owned arrays required by the out-of-place FFT contract:

- padded `real` FFT storage;
- `fourier`, used for the forward result and reused as destructive inverse input;
- `potential_k`, the preserved master potential spectrum;
- `poisson_kernel`, the cached Poisson/window-deconvolution/TreePM-split scalar operator.

The former `working_k` and full `grad_kx/grad_ky/grad_kz` arrays are removed. After Poisson multiplication, `potential_k` preserves the master spectrum. Before each force inverse, the solver derives the signed wave number from the current spectral layout and writes `(-i*k_axis) * potential_k` directly into `fourier`.

The derivative construction preserves both supported spectral ownership layouts: ordinary local-X slab ordering and FFTW-MPI transposed local-Y ordering. It uses the same signed-mode arithmetic as the previous precomputed gradient arrays and preserves the per-axis zero multiplier for an even-length Nyquist mode. The Poisson zero mode remains zero.

### Demand-driven real potential

`PmGridStorage` no longer allocates its real-space potential lane in the constructor. The default `solvePoissonPeriodic(...)` remains potential-producing for source compatibility. `ensurePotentialStorage()` and non-const `potential()` are explicit materialization boundaries; const potential access does not allocate, and `interpolatePotential(...)` requires materialized storage.

Production periodic TreePM and its periodic coarse zoom correction use `solvePoissonPeriodicForcesOnly(...)`, which performs the same force inverses without a potential inverse or potential-lane materialization. Isolated/open PM remains potential-producing because its current force construction consumes that representation.

### Memory model and observability

`estimatePmPlanResourcesMemory(...)` now models one real FFT array, two complex spectral arrays, and one scalar spectral array. For the usual distributed padded-real geometry this is conceptually `56 * alloc_local` bytes, while the backend-reported allocation extent remains authoritative.

The gravity preflight no longer charges an indexed-target XYZ gather. The periodic production TreePM grid estimate now charges density plus three force fields; potential is not charged to the force-only path. Runtime memory reporting contains no removed spectral/gather owners, and the PM potential entry reports its actual vector capacity, which is zero until materialized.

`spectral_operator_rebuilds` is retained as a compatibility telemetry field and now denotes reconstruction of cached Poisson/deconvolution/split operator state rather than reconstruction of full gradient arrays.

## Reproducibility and numerical scope

The PM equations, double precision, mesh dimensions, CIC/TSC stencil order, periodic wrapping, Poisson prefactor, assignment-window deconvolution, Gaussian TreePM split, zero mode, derivative sign, inverse normalization, and MPI request/response ordering are unchanged. Gradient multipliers are now recomputed from spectral geometry rather than loaded from retained arrays. The expressions and mode definitions intentionally match the prior construction, but bitwise equivalence has not been empirically demonstrated and must be checked by the later qualification campaign.

## Validation status

Build/tests: NOT RUN — intentionally deferred by campaign instruction.

Static source review and patch-hygiene inspection are the only validation activities permitted in this campaign. Runtime RSS savings are therefore not claimed here; all memory reductions described above are source-/ownership-derived architectural changes.
