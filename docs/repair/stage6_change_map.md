# Stage 6 Change Map (Current-State Reconnaissance)

Date: 2026-05-24
Scope: repository reconnaissance only (no production refactor in this patch).

## 1) Current particle containers
- **Canonical particle-state owner candidate:** `core::SimulationState::particles` (`ParticleSoa`) + `core::SimulationState::particle_sidecar` (`ParticleSidecar`) together hold core particle truth lanes today.
- Additional persistent state owners under `SimulationState`:
  - `cells` (`CellSoa`), `gas_cells` (`GasCellSidecar`), `patches` (`PatchSoa`)
  - species sidecars: `star_particles`, `black_holes`, `tracers`
  - auxiliary persistent stores: `gas_identity`, `species_counts`, `sidecars`.

## 2) SoA / AoS / adapter patterns
- **SoA dominant lanes:** `ParticleSoa`, `ParticleSidecar`, `CellSoa`, `GasCellSidecar`, species sidecars (`include/cosmosim/core/simulation_state.hpp`).
- **Generic SoA adapter:** `ParticleSoaStorage` and field enum adapters (`include/cosmosim/core/soa_storage.hpp`) used only by tests/benches as a packed SoA utility façade; it exposes `k_owns_persistent_particle_truth=false` and `k_is_restart_serializable=false`.
- **AoS in hydro math kernels:** `HydroPrimitiveState`, `HydroConservedState`, `HydroFaceFlux` are per-item structs; vectorized storage exists via hydro SoA/cache wrappers in `hydro_core_solver`.
- **Adapters/views:** active-set builders in `src/core/simulation_state_active_views.cpp` map sparse global rows to compact workspace spans and read-only/read-write view structs.

## 3) Current active-set handling
- Scheduler authority: `HierarchicalTimeBinScheduler` + `IntegratorState` + `ActiveSetDescriptor` (`include/cosmosim/core/time_integration.hpp`).
- Per-step active indices stored transiently on `SimulationState` (`active_particle_indices`, `active_cell_indices`) and/or passed via `StepExecutionContext`.
- Active compact mirrors built via `build*ActiveView` helpers using `TransientStepWorkspace`.
- Hydro has dedicated active-set view path (`HydroActiveSetView`) and active-only patch advancement APIs.

## 4) Persistent fields (restart-authoritative intent)
- Particle hot/canonical lanes: positions, velocities, mass, scheduler mirror `time_bin` (mirror only), ids/species/flags/owning rank, drift-epoch lanes, softening-override lanes.
- Gas/cell/piecewise persistent lanes: cell centers/mass/patch/time-bin mirror, gas identity/thermodynamic lanes; reconstruction gradients are transient scratch, patch topology ranges.
- Species sidecars: star/BH/tracer persistent bookkeeping lanes.
- Integrator + scheduler persisted state and provenance/config metadata are explicit restart payload members.

## 5) Transient/scratch fields
- `TransientStepWorkspace` explicit compact buffers for particle/cell/hydro active kernels.
- Hydro reconstruction/riemann scratch/cache buffers (`HydroScratchBuffers`, primitive caches).
- PM/tree/hydro kernel temporary vectors and staging buffers in gravity/hydro source implementations.
- Active descriptor metadata freshness fields (`source_generation`, etc.) are derived runtime-only descriptors.

## 6) Restart serialization fields (current)
- Restart payload and schema are centralized in `io/restart_checkpoint.*` with explicit schema/version contract (`cosmosim_restart_v11`).
- Serialized families include simulation state lanes, species sidecars, module sidecars, integrator state, scheduler persistent state, distributed gravity restart state, normalized config, provenance, integrity hash.
- Docs explicitly state cached PM long-range field arrays are **not** serialized; policy is deterministic rebuild on resume.

## 7) Memory/accounting/profiling utilities
- Alignment/memory lane primitives: `AlignedAllocator`, `AlignedVector`, `SoaFieldArray` in `soa_storage.hpp`.
- Runtime profiling session + scopes in `core/profiling.hpp/.cpp` and integration hooks in time integration and kernel modules.
- Bench coverage for layout and overhead: `bench_layout_smoke`, `bench_species_hot_cold_access`, `bench_active_kernel_compact_vs_full`, `bench_profiling_overhead`, `bench_io_restart_kernel`, etc.

## 8) Kernels still operating on broad full-state objects (Stage 6 targets)
- Time-integration stage handlers receive `SimulationState&` through callback context, and many stage implementations continue to traverse broad state surfaces.
- Tree gravity / TreePM coupling entry points still accept broad state-oriented inputs and can pull mixed hot/cold lanes unless narrowed per-kernel.
- Restart read/write hash/validation paths walk many full-state lanes by design; keep broad in I/O but isolate from hot stepping kernels.
- Physics modules (stellar, BH, tracers) often accept global state + active indices; hot-path access narrowing is incomplete.

## 9) Duplicated or ambiguous ownership risks
- `ParticleSoaStorage` is explicitly quarantined as utility/test storage; current call sites are tests/benchmarks only.
- `time_bin` mirrors exist on particle/cell SoA while scheduler `bin_index` is authority; drift risk exists if mirror refresh discipline regresses.
- Gas identity has both `gas_cells` parent-particle lanes and `GasCellIdentityMap`; comments indicate seam not fully wired everywhere (intentional migration seam).
- Optional/module sidecars can duplicate derived counters/metadata unless invariants are consistently documented and tested.

## 10) Existing tests/benchmarks covering Stage 6-like contracts
- Unit: `unit_simulation_state`, `unit_time_integration`, `unit_hot_cold_sidecar_layout`, `unit_profiling`, `unit_restart_checkpoint_schema`, `unit_hydro_core_solver`.
- Integration: `integration_reorder_compaction_sidecars`, `integration_species_migration_invariants`, `integration_softening_ownership_invariants`, `integration_hierarchical_time_bins`, `integration_hierarchical_timestep_regression`, `integration_restart_checkpoint_roundtrip`, `integration_profiling_mini_run`, `integration_time_integration_loop`, `integration_tree_pm_coupling_periodic`.
- Benchmarks: `bench_active_kernel_compact_vs_full`, `bench_species_hot_cold_access`, `bench_soa_storage`, `bench_tree_gravity`, `bench_tree_pm_coupling`, `bench_hydro_kernels`, `bench_io_restart_kernel`, `bench_restart_checkpoint_latency`, `bench_profiling_overhead`.

## 11) File-touch map for Stage 6.x prompts
- **6.1 Canonical particle ownership + adapter demotion:** `include/src/core/simulation_state*`, `include/src/core/soa_storage*`, targeted gravity/hydro/physics call sites + docs/state model/update notes.
- **6.2 Hot/cold lane hardening:** `include/src/core/simulation_state*`, species sidecars, invariants tests under `tests/unit` + `tests/integration`.
- **6.3 Active view enforcement:** `src/core/simulation_state_active_views.cpp`, `include/core/time_integration.hpp`, gravity/hydro consumers.
- **6.4 Gravity kernel surface narrowing:** `include/src/gravity/tree_gravity*`, `tree_pm_coupling*`, `pm_solver*`, relevant benches.
- **6.5 Hydro kernel surface narrowing:** `include/src/hydro/hydro_core_solver*`, hydro integration bridges.
- **6.6 Transient scratch non-persistence audit:** `include/src/core/simulation_state*`, `src/io/restart_checkpoint*`, restart schema tests/docs.
- **6.7 Memory accounting + profiling lanes:** `include/src/core/profiling*`, potentially per-module counters and benches.
- **6.8 Ownership ambiguity cleanup:** `core adapters + module interfaces` with migration glue in tests.
- **6.9 Targeted invariant tests:** add/update focused unit+integration tests for active views, ownership, non-serialization of transient state.
- **6.10 Perf/regression benchmark gate:** benchmark updates under `bench/`, baseline notes, CMake preset/test label integration.

## Data-layout classification snapshot
- **Persistent hot:** particle position/velocity/mass; gravity-facing cell geometry/mass; active force-update needed lanes.
- **Persistent cold:** ids/species/flags/rank, stellar/BH/tracer rich bookkeeping, patch metadata, slow-changing identity/provenance metadata.
- **Optional sidecar:** module sidecars + species-specific sidecars not required by all kernels.
- **Derived view:** active compact mirrors + read-only descriptors + scheduler summaries.
- **Transient scratch:** `TransientStepWorkspace`, hydro reconstruction scratch, PM/tree temporary work arrays/stacks, MPI staging buffers.
- **External I/O schema:** restart payload schema families + HDF5 path contracts.
- **Diagnostics/provenance:** profiler events, module metadata blocks, provenance/config hashes, runtime counters.

## Restart/provenance risk flags
- Must remain non-authoritative: active index lists, compact active mirrors, hydro reconstruction scratch, PM working arrays, tree traversal stacks, MPI exchange buffers, output staging scratch.
- At-risk surfaces for accidental serialization creep: module sidecar payload expansion, distributed-gravity temporary caches, any new workspace lanes added near restart payload builders.

## Unknowns requiring explicit verification in later prompts
- ParticleSoaStorage call sites were enumerated during the 6.0-6.4 repair pass: tests and benchmarks only, no production runtime truth use.
- Some gravity/hydro hot loops still need line-level audit to quantify broad-state dereference frequency.
- GPU/CUDA parity for active-view narrowing requires dedicated review path (`pm_cuda_kernels.cu` + CUDA benches/tests).
