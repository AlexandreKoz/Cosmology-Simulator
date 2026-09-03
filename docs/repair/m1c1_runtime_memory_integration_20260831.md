# M1C-1 runtime memory integration closure

Status: source integration complete; focused CPU validation complete; broader dependency validation recorded separately.
Date: 2026-08-31
Campaign: M1C-1 Runtime Memory Integration

## Scope

M1C-1 integrates the M1B process-level `core::MemoryGovernor` with the remaining
important DMO/general-runtime high-water boundaries. It preserves M1A bounded PM
routing and later MPI collective/count/topology repairs. This document is not an
M1 acceptance declaration; M1C-2 owns final acceptance and scaling closure.

## Governed-allocation coverage inventory

| Runtime owner / path | M1C-1 classification | Physical/accounting contract |
| --- | --- | --- |
| Canonical `SimulationState` particle/cell/species/sidecar state | baseline-owned | Capacity-based `MemoryReport` ownership; never reserved a second time. |
| Particle/gas hierarchical schedulers | baseline-owned | Persistent retained capacity is included in governor baseline reconciliation. |
| `TransientStepWorkspace` monotonic scratch | governed | Existing M1B block reservation survives logical reset for the full backing-block lifetime. |
| CHUI-owned PM/FFT mesh/solver arrays | baseline-owned after phase | TreePM-owned retained capacities remain visible once in the gravity memory report; M1C-1 transfers newly retained capacity atomically from the phase commitment into baseline accounting. |
| PM density/force/potential routing | bounded + governed phase | M1A peer-bounded coordinated-round routing is unchanged; the enclosing gravity high-water is admitted before materialization. Total traffic is not used as live-memory demand. |
| Opaque FFTW/MPI/HDF5/backend allocations | external-runtime reserve | Remain outside CHUI-owned byte claims and are represented by the existing configured unknown-reserve policy. |
| Tree nodes/order/source/target/periodic staging/LET exchange | governed phase -> baseline-owned retained capacity | Conservative TreePM high-water is admitted collectively before O(N) cache/vector growth. Known source/target vector cardinalities reserve explicitly to avoid allocator-geometric surprise capacity. Retained capacity is atomically reconciled into baseline before the commitment is released. |
| Runtime decomposition planning items | governed phase resident | Exact explicit vector footprint is admitted collectively and scoped so decomposition items are destroyed before migration staging begins. |
| Migration transaction old/new coexistence | governed communication | A conservative transaction/staging reservation is admitted collectively before maps, migration records, wire payloads, decode, scheduler remap, and commit. Retained canonical/scheduler capacity is atomically transferred to baseline after commit. |
| Migration preserved scheduler state | compact derived state | Preserved particles use `TimeBinSchedulerIdentityRecord`; the former full `ParticleMigrationRecord` reconstruction copy is removed. |
| Migration wire transport | bounded packet + governed | Post-campaign closure uses an explicit conditional field codec and reusable fragmented packets capped by the authoritative MPI transport-round budget. Total logical traffic may grow, but simultaneous application wire staging is bounded; see `docs/repair/m1c1_postcampaign_closure_20260902.md`. |
| Snapshot scientific readback | governed diagnostic | Actual readback staging consumes the existing planned output/restart overlap allowance first; only bytes beyond the predeclared allowance become an additional commitment. Reservation outlives the readback object. |
| Restart write/readback verification | governed diagnostic | State/scheduler/readback/force-cache coexistence is admitted collectively. Required restart verification is unchanged. The reservation lifetime covers the physical readback/export objects. |
| Optional science diagnostics | pressure scheduled | Run-health diagnostics remain required. Science-light/heavy products are deferred under Red/Trip pressure and record `analysis.memory_pressure_deferral`; numerical solver semantics are unchanged. |
| Hydro/AMR canonical-state redesign, chemistry, star/BH memory tiers, feedback memory | deferred to M2 | M1C-1 preserves compatibility but does not redesign full-physics memory architecture. |
| GPU/device memory governor research | deferred | No device-memory policy redesign in M1C-1. |

## Process-memory observation and reconciliation

Linux/WSL production observation uses low-frequency workflow-boundary sampling:

- current RSS and RSS high-water from `/proc/self/status` (`VmRSS`, `VmHWM`);
- PSS from `/proc/self/smaps_rollup` when available;
- unavailable or malformed observations remain `std::nullopt`, never fake zero.

The process observation is diagnostic evidence, not allocation authority. The
hard admission path remains deterministic governor accounting plus configured
external-runtime reserves.

The profiler reports:

```text
process_memory.current_declared_residency_bytes
process_memory.known_accounted_bytes  # compatibility alias
process_memory.observed_rss_bytes
process_memory.observed_peak_rss_bytes
process_memory.observed_pss_bytes
process_memory.unexplained_resident_bytes
process_memory.observed_to_known_ratio
```

with

```text
unexplained_resident_bytes = max(0, observed_rss_bytes - current_declared_residency_bytes)
```

when RSS exists. Post-campaign closure separates admission policy from current
residency: `current_declared_residency_bytes` includes baseline-owned retained
capacity, committed governed allocations, and the configured opaque external
runtime estimate, while planned future output overlap and uncommitted
reservations remain policy-only demand. `known_accounted_bytes` is retained as a
compatibility alias for the corrected current-residency quantity.

## Distributed reconciliation

`attachDistributedMemoryTelemetry()` performs existing-authoritative MPI
reductions for:

- known/accounted demand;
- current RSS;
- peak RSS;
- communication high-water.

Each metric records local, global sum, rank maximum, rank mean, and
`rank_max / rank_mean`. Optional OS observations are marked valid only when every
rank provides the metric. Rank maximum is the safety quantity; averages are
telemetry only.

## Physical-lifetime transfer

M1C-1 adds `MemoryReservation::reconcileBaselineOwnedAndRelease()`. A committed
phase can replace the measured baseline and remove its own commitment under the
same governor lock. This prevents both failure modes at phase closure:

1. temporarily double-counting retained bytes as both baseline and commitment;
2. releasing a commitment before physically retained capacity becomes baseline-owned.

Gravity and migration use this operation after successful phase completion.

## PM routing acceptance continuity

M1C-1 does not alter M1A PM algorithms or batch policy. The current engineering
model remains:

```text
M_workspace = M_fixed + 2 * (world_size - 1) * B_peer_effective
M_workspace < 128 MiB/rank
```

for the certified 512^3 / 8-rank TSC routing profile. Current project PM
documentation records the deposition-routing model at 127.937012 MiB/rank; M1C-2
must re-audit the source formula and final acceptance evidence, not change it to
make a different number.

## Focused M1C-1 validation

The CPU-only focused campaign build compiled the modified core and workflow
objects, then passed:

```text
unit_memory_governor
unit_scratch_memory_governor
unit_process_memory
unit_memory_accounting
unit_profiling
integration_reference_workflow
integration_species_migration_invariants
integration_gas_cell_identity_migration
unit_pm_solver
integration_pm_periodic_mode
```

`integration_reference_workflow` includes:

- finite process-governor execution;
- a one-byte hard-limit rejection before gravity high-water materialization;
- Red-pressure optional-analysis deferral with the physics step still completing;
- process/distributed-memory profiler schema checks.

Dependency-enabled HDF5/FFTW/MPI validation belongs to the broader post-artifact
validation record and M1C-2 final acceptance evidence.

## M1C-2 handoff

M1C-2 should perform the final hostile acceptance audit and scaling evidence,
including dependency-enabled distributed execution where the environment
permits it. It should not reopen M2 full-physics representation research unless
a defect directly invalidates M1 acceptance.
