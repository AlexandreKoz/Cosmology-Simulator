# M1C-1 post-campaign closure repair and TreePM edge-contract hardening

Status: source closure achieved; focused CPU and HDF5 validation passed; FFTW/MPI runtime validation blocked by environment.
Date: 2026-09-02
Campaign: M1C-1 Post-Campaign Closure Repair and TreePM Edge-Contract Hardening

## Scope and verdict

This repair closes the source-level residuals identified after M1C-1 without
reopening M1C-2 final acceptance or M2 full-physics memory architecture.
Existing M1A bounded PM routing, M1B process memory governance, M1C-1 gravity
admission/lifetime integration, later MPI collective-safety repairs, and the
validated TreePM force split remain authoritative.

Campaign verdict:

```text
CLOSURE REPAIR ACHIEVED
DISTRIBUTED/FFTW VALIDATION: BLOCKED BY ENVIRONMENT
```

This is not an `M1 ACHIEVED` declaration. M1C-2 owns final M1 acceptance and
production-scale evidence.

## Finding closure matrix

| Finding | Status | Closure |
| --- | --- | --- |
| `M1C1-CLOSURE-001` migration wire residency | CLOSED | Full-population send/receive wire materialization is replaced by reusable bounded per-rank fragment packets. Logical migration traffic may grow, but simultaneous packet staging is capped by the existing MPI transport-round budget. |
| `M1C1-CLOSURE-002` migration peak admission | CLOSED | The old `2 * baseline` heuristic is replaced by a checked source-grounded plan covering old/new commit coexistence, native record capacities, dynamic sidecar heap, scheduler remap, index maps, bounded packet buffers, and maximum single-record fragment assembly. |
| `M1C1-CLOSURE-003` snapshot readback admission | CLOSED | Snapshot write failure is coordinated first; readback reservation is then collectively agreed before any rank enters HDF5; readback/verification failures are coordinated afterward. |
| `M1C1-CLOSURE-004` policy versus current residency | CLOSED | Current RSS reconciliation uses current declared residency rather than future policy headroom. Planned overlap and uncommitted reservations remain governor policy demand, not current residency. |
| `TREEPM-CLOSURE-001` even-grid Nyquist derivative | CLOSED | First-derivative multipliers are zero on an even-axis Nyquist mode while the scalar Poisson mode is preserved; near-Nyquist modes remain differentiated. |
| `TREEPM-CLOSURE-002` softening/cutoff contract | CLOSED | Periodic Plummer TreePM is restricted to the existing source-tested `epsilon/r_s <= 0.20` envelope. Static configured softenings and runtime per-particle overrides are both guarded. |
| `TREEPM-CLOSURE-003` aggregate classic-MPI displacement ceiling | CLOSED BY SOURCE/PLANNER | Sparse TreePM target batches are further capped from peer degree and classic-MPI aggregate byte capacity, with a global minimum agreed across ranks. Large logical traffic therefore spans more physical rounds rather than one aggregate `INT_MAX` displacement span. Real MPI execution is environment-blocked in this validation environment. |
| `TREEPM-CLOSURE-004` rectangular PM scientific contract | CLOSED | Cubic/default behavior is unchanged. Periodic rectangular PM is accepted only when the spherical split scale based on the current representative spacing is at least the coarsest PM cell spacing; unsupported strong anisotropy fails with a precise configuration diagnostic. |

## Migration wire format

Migration now uses an explicit versioned little-endian field codec rather than
native struct-layout/padding as wire authority. Common particle fields are
encoded once, and gas/star/BH/tracer state is encoded only when the matching
presence bit is set. Module sidecar payloads remain explicit payload entries.

Measured source-level reference widths from the repaired codec are:

| Record content | Wire bytes |
| --- | ---: |
| DMO/common particle with scheduler state, no per-particle softening lane | 125 |
| Optional per-particle softening lane | +8 |
| Gas fields | +154 |
| Star fields | +220 |
| Black-hole fields | +68 |
| Tracer fields | +44 |

The post-M1C-1 audit measured the previous fixed particle wire at approximately
660 bytes before module sidecars, including inactive baryonic structs. The
repaired common DMO record is therefore about 81% smaller before optional
softening/sidecar data:

```text
old audited fixed DMO wire ~= 660 B/particle
new common DMO wire          = 125 B/particle
new DMO + softening lane     = 133 B/particle
```

## Bounded migration packet contract

`parallel::mpiTransportRoundLimitBytes()` remains the transport authority. Its
production default is 16 MiB. `migration_wire::planPacketCapacity()` divides
that aggregate physical budget across ranks and reserves a fixed fragment
header per active peer.

For the M1 architectural reference topology of eight ranks:

```text
transport round limit                  = 16,777,216 B (16 MiB)
per-peer packet capacity               =  2,097,152 B (2 MiB)
fragment header                        =         32 B
per-peer fragment payload capacity     =  2,097,120 B
aggregate application packet capacity  = 16,777,216 B (16 MiB)
```

The application reuses these packet vectors across rounds. One oversized record
may itself be fragmented, so one record does not force a larger packet bound.
The MPI bounded-exchange implementation may maintain its own round-local send
and receive buffers, so the migration governor conservatively budgets both the
application packet vectors and the MPI internal round buffers.

Consequently, at fixed rank count and transport limit:

```text
M_packet_peak does not scale with total migrated record count.
```

Total wire traffic still scales with migration volume, as it should.

## Migration physical peak admission

Before target-scale migration record/wire materialization, the runtime performs
checked local sizing followed by bounded control exchanges for inbound demand.
The resulting admission plan includes:

- required old/new canonical commit coexistence;
- outbound and inbound native `ParticleMigrationRecord` / AMR record capacity;
- dynamic module-sidecar/AMR heap upper bounds;
- scheduler identity/remap and destination-ID staging;
- material index-map capacity;
- bounded send/receive application packet capacity;
- bounded send/receive MPI internal round capacity;
- maximum simultaneous outbound/inbound single-record fragment assembly;
- AMR patch ownership changes that may induce parent gas-particle migration.

All arithmetic uses checked size/count operations. The governor reservation is
attempted before migration records, wire packets, maps, and commit staging are
materialized. Where transaction ownership permits it, inbound record vectors are
moved into commit state rather than copied wholesale.

The runtime profiler records the repaired contract through migration telemetry,
including transaction reservation, record capacity, packet send/receive
capacity, scheduler remap, index-map, and communication high-water fields.

## Snapshot collective-admission contract

The snapshot path now follows:

```text
local snapshot write
-> collective write-failure agreement
-> local readback reservation attempt
-> collective admission agreement
-> local HDF5 readback + validation
-> collective readback/validation failure agreement
```

This matches the stronger restart pattern: no rank starts governed readback
because only its local reservation succeeded. No new global barrier was added;
the existing failure coordinator remains the authority.

## Current residency versus governor policy

Governor `accounted_bytes` remains the conservative scheduling/admission value
and may include baseline ownership, opaque external reserve, planned future
output/restart overlap, committed allocations, and reservations.

Process residency reconciliation now exposes
`current_declared_residency_bytes`, defined from current baseline-owned bytes,
current committed governed bytes, and the configured opaque external-runtime
estimate. It excludes inactive future output/restart overlap and excludes
uncommitted reservations. `known_accounted_bytes` remains as a compatibility
alias for this corrected current-residency value.

When RSS exists:

```text
unexplained_resident_bytes =
    max(0, observed_rss_bytes - current_declared_residency_bytes)
```

RSS/PSS remain diagnostic evidence and do not become allocation authority.

## PM Nyquist derivative contract

The periodic PM solver continues to solve the same scalar Poisson problem and
uses the same Gaussian TreePM split, assignment/deconvolution, normalization,
and physical wave numbers. Only the odd first-derivative edge convention is
changed.

For each differentiated axis independently:

```text
if N_axis is even and mode_index == N_axis / 2:
    first-derivative multiplier = 0
else:
    first-derivative multiplier = physical k_axis
```

This preserves the Hermitian contract of a real inverse FFT at the self-conjugate
Nyquist mode. The scalar potential coefficient is not zeroed solely because one
derivative component is Nyquist. Tests cover pure Nyquist, mixed
Nyquist/ordinary-axis, and near-Nyquist modes.

## TreePM scientific support contracts

### Plummer softening versus split scale

The existing softened TreePM residual formulation is preserved. The current
source-controlled composition evidence covers `epsilon/r_s <= 0.20`; periodic
TreePM now fails explicitly beyond that envelope instead of silently truncating
an uncontrolled non-compact Plummer residual at fixed `r_cut`. Static configured
softenings are checked in configuration validation and per-particle overrides are
checked during gravity-source materialization.

### Rectangular PM

The existing spherical split based on the representative geometric-mean PM cell
spacing remains unchanged for validated configurations. Periodic rectangular PM
is rejected when that split scale is smaller than the coarsest PM cell spacing.
Cubic configurations are unaffected. This is an explicit support boundary, not
a claim that arbitrary anisotropic PM meshes have been scientifically certified.

## TreePM aggregate classic-MPI rounds

The existing sparse peer graph and per-peer batch configuration are preserved.
The repair adds a planner that computes the maximum targets per peer that fit the
aggregate classic-MPI byte/displacement limit for the current sparse peer degree.
Ranks collectively agree on the minimum safe target count before executing the
exchange so all ranks preserve collective ordering.

Planner tests cover zero peers, multi-peer aggregate limits, deterministic
multi-round coverage, and a case where even one target per peer cannot fit the
artificial classic-MPI byte limit. Real MPI execution remains pending because
this environment has no MPI C++ development/runtime support.

## Validation evidence

### Repository hygiene

Passed after validation build directories were moved outside the source tree:

```bash
bash ./scripts/ci/check_repo_hygiene.sh
```

### Focused CPU validation

The repaired CPU Debug tree passes the campaign-specific tests. An early focused
combined invocation completed the first 14 tests before the command execution cap
stopped while starting `integration_tree_pm_coupling_periodic`; the final two were
then run separately and passed.

Passed focused tests:

```text
integration_reference_workflow
unit_config_parser
unit_profiling
unit_memory_accounting
unit_process_memory
unit_memory_governor
unit_scratch_memory_governor
integration_species_migration_invariants
unit_migration_wire
integration_gas_cell_identity_migration
unit_pm_solver
unit_tree_gravity
integration_pm_periodic_mode
unit_tree_pm_split_kernel
integration_tree_pm_coupling_periodic
integration_tree_gravity_vs_direct
```

No force tolerance was relaxed.

### Full CPU build and suite

The complete CPU Debug build subsequently finished successfully:

```bash
cmake --build --preset build-cpu-debug -j2
```

Because the execution harness caps one shell invocation, the full 138-test CTest
preset was completed in two pieces. The first `ctest --preset test-cpu-debug
--output-on-failure -j2` invocation passed tests 1 through 127 before the command
cap interrupted while tests 128/129 were starting. The explicit continuation for
tests 128 through 138 then passed all 11 remaining tests. Across the two
invocations:

```text
CPU Debug tests: 138 / 138 PASS
failures: 0
```

The continuation includes `validation_periodic_ewald_reference`,
`validation_phase2_mpi_gravity_single_rank`, and
`integration_dmo_empty_rank_smoke_single_rank`, providing additional local
gravity/periodic/empty-rank evidence without pretending to be multi-rank MPI
runtime validation.

### HDF5 validation

`cmake --preset hdf5-debug` configured successfully with HDF5 1.14.5. The
following repaired-tree tests passed:

```text
integration_reference_workflow
integration_snapshot_hdf5_roundtrip
integration_restart_checkpoint_roundtrip
integration_restart_equivalence_output_enabled
```

This directly exercises the changed output/restart implementation with HDF5
enabled.

### FFTW blocker

Attempted once:

```bash
cmake --preset pm-hdf5-fftw-debug
```

Configuration failed because the FFTW3 serial double-precision development
library is unavailable:

```text
COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library was not found.
Debian/Ubuntu package: libfftw3-dev.
```

No FFTW runtime pass is claimed.

### MPI blocker

Attempted once:

```bash
cmake --preset mpi-hdf5-fftw-debug
```

Configuration failed because MPI C++ support is unavailable:

```text
Package 'mpi-cxx', required by 'virtual:world', not found
Could NOT find MPI (missing: MPI_CXX_FOUND CXX)
```

No MPI or distributed FFTW runtime pass is claimed.

## M1C-2 handoff

M1C-2 should now concentrate on final acceptance evidence rather than reopening
these source contracts. In a dependency-complete environment it should run the
real MPI+FFTW gravity/memory matrix, including migration many-round cases,
collective snapshot admission failure, TreePM aggregate-round execution,
Nyquist FFTW reference cases, empty/uneven ranks, repeated-step memory plateau,
rank-max RSS/governor high-water, and the planned 512^3/8-rank analytical and
runtime acceptance evidence.

M2 remains responsible for hydro/AMR/chemistry/star/BH/full-physics memory
architecture.
