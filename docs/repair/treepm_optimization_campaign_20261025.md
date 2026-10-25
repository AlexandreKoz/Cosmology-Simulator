# TreePM optimization and memory campaign — 2026-10-25

**Status: IMPLEMENTED — NOT COMPILED — NOT TESTED — OWNER VALIDATION REQUIRED.**

The project owner explicitly deferred all compilation, execution, testing and
benchmarking for this session. Nothing below is measured; byte counts are
derived from allocation sites in the source. Every faster policy is opt-in; all
defaults reproduce the previous behavior.

## Implementation matrix

| Campaign | Files | What was done | Status |
| --- | --- | --- | --- |
| A PM memory | `pm_solver.hpp/.cpp`, `gravity_memory.hpp/.cpp`, `config.*`, `gravity_runtime.cpp` | `potential_k` (pristine `phi_k` copy) is now allocated on the first spectral solve instead of at plan creation. New opt-in `PmForceGradientBackend::kFiniteDifference4` (`numerics.treepm_pm_force_gradient = finite_difference_4`): one inverse FFT, 4th-order centered difference, two-plane MPI slab halo exchange, never allocates `potential_k`. Preflight `estimatePmPlanResourcesMemory` and live memory report account backend-specific arrays and FD halo planes. In-place real/complex FFT aliasing was **not** implemented (see risks). | Partially implemented — untested |
| B communication arena | `gravity_communication_arena.hpp`, `tree_pm_coupling.*`, `gravity_runtime.cpp` | `configure()` grows at idle phase boundaries with a governed, strong-exception-safe swap (old capacity intact on refusal), never shrinks implicitly. Telemetry: per-phase high-water, growth events, unused capacity. The estimate already combined phases with `max()` (the audit's "sum" premise does not hold here); its dominant term is `2 (P-1) treepm_pm_exchange_batch_bytes`, a configured bound, so no fixed size was substituted. | Implemented — untested |
| C group-walk observability | `tree_pm_coupling.*`, `gravity_runtime.cpp`, `bench_tree_pm_sweep.cpp`, `configs/performance/…` | `TreePmGroupWalkStatus` (first failing runtime condition), per-group fallback reasons (no leaf / softening / history / geometry), groups, grouped targets, max group size, list high-water. **Repair:** the group walk was silently never used for owned targets whenever the short-range path was distributed (any MPI run with peers); owned targets now go through one shared evaluator in both paths. Performance profile deck added. | Implemented — untested |
| D leaf culling / grouping | `internal/tree_pm_group_walk.hpp`, `tree_pm_coupling.*`, `config.*` | Interaction lists keep one entry per source leaf with minimum-image cell geometry; each target culls leaves whose cell is entirely beyond `r_cut` (conservative `1+1e-9` slack, no square roots), surviving contiguous ranges are re-coalesced for the vectorized kernel. `numerics.treepm_group_walk_target_size` decouples target groups (2..block size, extent ≤ 0.5 r_cut) from the source leaf size. Counters: candidate/culled leaves, candidate particles. | Implemented — untested |
| E2 tree refit | — | `refitWithinOriginalLeafCells` already certifies geometric leaf membership, layout generation and softening, and falls back to rebuild; unchanged. | Already present in source |
| E3 LET graph reuse | `tree_pm_coupling.*`, `gravity_runtime.cpp` | Cached symmetric sparse-peer graph is reused when the request set is a subset of its peers (previously exact match only), still decided collectively; `graph_cache_superset_hits` counter. Source payloads are still exchanged every solve. | Implemented — untested |
| E4/E5 fine-tick cost, PM cadence | — | Analysis/output stages already return early on non-synchronization fine ticks; `treepm_update_cadence_steps = 1` restriction retained (no PM reuse qualification implemented). | Not implemented beyond E3 |
| F MPI consensus | `tree_pm_coupling.cpp`, `time_coordinator.cpp` | Two diagnostic empty-rank reductions per solve fused into one; four population-count reductions per global step fused into one (exact integer sums). Floating-point timestep reductions were deliberately not fused (reduction-order reproducibility). No correctness collective removed. | Partially implemented — untested |
| G output cadence | `config.*`, `output_restart_runtime.*`, `snapshot_hdf5.hpp`, `snapshot_set.cpp` | `output.restart_interval_steps` (independent checkpoints; derived next index, no schema change) and `output.snapshot_verification = structural` (header-only member check via new `inspectSnapshotMemberHdf5`; SHA-256 sidecar, completion marker and set validation retained). Counters for snapshot/restart writes and verification time. | Implemented — untested |
| H gravity eta | `config.*`, `time_coordinator.cpp` | `numerics.gravity_timestep_eta` replaces the three hardcoded `0.2` uses (particle, gas cell, hierarchical acceleration + encounter criteria). | Implemented — untested |
| I analysis cadence | `configs/performance/…` | Cadence already honored and skipped on fine ticks; performance deck sets `run_health_interval_steps = 8`. | Already present in source |
| J governor reconciliation | — | `observeProcessMemory()` (RSS, HWM, PSS via `smaps_rollup`) is already attached to the step memory report with known/unexplained reconciliation. Arena/PM accounting updated as above. | Already present in source |
| K particle compaction | — | No precision or ID narrowing (requires empirical qualification). | Not implemented |
| L hybrid MPI/OpenMP | — | New scratch is per worker (group lists) or per plan (FD halo); FD halo exchange is collective-uniform and requires ≥ 2 x planes per rank (checked collectively). | Partially implemented — untested |

## PM buffer residency (from allocation sites)

Let `N = nx ny nz` (local `N_l`), `A` = FFTW complex cells (`nx ny (nz/2+1)`
serial, backend-queried per rank when distributed).

| Path | Before | After (spectral, default) | After (finite_difference_4) |
| --- | --- | --- | --- |
| Periodic serial host | real `N` dbl + fourier `A` cplx + potential_k `A` cplx (at plan creation) + 3 force `N` dbl (+ optional potential `N`) | same live set, potential_k allocated on first spectral solve | real + fourier + 3 force; **no potential_k** |
| Periodic distributed | real `2A` dbl + fourier `A` + potential_k `A` + 3 force `N_l` + 6 halo planes | same | potential_k removed; + 4 retained `ny nz` halo planes |
| Isolated / CUDA | unchanged | unchanged | rejected (validation) |

Compact density is already released in favor of FFT-backed density (M48-09)
and is not counted. For a 256³ periodic serial mesh: `A = 8.45e6`, so
`potential_k` is 135 MB; the FD backend removes it (~17% of the 806 MB PM live
set) and two of three inverse FFTs. The audit's eight-grid assumption does not
match this checkout.

## Algorithm notes

* **Leaf culling.** For leaf cell center `c` (half width `h`) and target `p`,
  the per-axis distance `max(0, |minimg(c-p)| - h)` lower-bounds every pair
  separation in the leaf (tree cells contain their particles after build and
  after certified refit). Culling skips only pairs the kernel masks; forces
  differ from unculled lists only by summation order.
* **Group fallback.** A group falls back before any accumulation; fallback
  targets use the unchanged per-target walk, so forces of fallback targets are
  bitwise those of the per-target path.
* **Distributed group walk.** The distributed path now evaluates owned
  targets in tree visit order over the batch's visit positions (the batches
  still partition the active set); remote (incoming) targets keep the
  per-target walk.
* **LET subset reuse.** Graph communicators are symmetric; when each rank's
  request set is contained in its cached neighbours, every needed edge exists
  on both endpoints. Reuse remains a collective MIN decision.

## Remaining correctness risks

* FFTW/MPI: in-place r2c/c2r aliasing was not implemented — the serial real
  array is unpadded and several paths index it directly. The FD halo exchange
  assumes x-slab real-space layout (true for slab and transposed-output pencil
  plans) and contiguous rank ordering of x planes.
* Force accuracy: `finite_difference_4` attenuates mesh-scale forces by `T(kh)`
  and its interaction with window deconvolution and the Gaussian split has not
  been measured. Group-target sizes > leaf size change MAC decisions (larger
  group AABB → more conservative acceptance), not the MAC formula.
* LET subset reuse keeps stale extra neighbours until the next miss; they
  exchange zero bytes but count in `let_candidate_peer_count`.
* Arena growth briefly holds old+new capacity under the governor.
* Restart: no schema change. New normalized keys make older restarts fail
  closed on normalized-text mismatch; the v23 migration projection does not
  include them.
* Structural snapshot verification does not compare particle payloads; the
  SHA-256 sidecar proves file integrity, not semantic equality with state.

## Owner validation checklist

1. Build presets from `docs/build_instructions.md` (cpu-only debug, HDF5+FFTW,
   MPI) — nothing was compiled.
2. `unit_config_parser`, `unit_output_cadence` (new),
   `unit_gravity_communication_arena` (new), `integration_config_examples`.
3. `integration_tree_pm_coupling_periodic` single rank and MPI 2/4 ranks:
   new `testGroupWalkLeafCullingMatchesPerTargetWalk`,
   `testPmFiniteDifferenceGradientBackend`, plus existing distributed tests
   (distributed group walk path changed).
4. TreePM accuracy validation (`test_tree_pm_ewald_accuracy`) with group walk
   target size 0/16/32 and with `finite_difference_4`.
5. Restart/output: snapshot+checkpoint cadence combinations, resume
   equivalence with `restart_interval_steps`, structural verification on
   sharded multi-rank output.
6. Hierarchical workflow tests (eta wiring, graph reuse on fine ticks).
7. Benchmarks: `bench_tree_pm_sweep --group-walk on --group-target-size 32`
   and the performance deck versus a strict reference run and AREPO.
