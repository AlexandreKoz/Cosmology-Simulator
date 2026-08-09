# I/O Post-Campaign Acceptance Closure Matrix — 2026-08-09

This matrix maps the acceptance re-audit findings to the focused repair. "Closed in source /
environment-limited runtime" means the implementation and durable CI acceptance gate are
present, but the local campaign host did not provide that platform/runtime.

## PC-01 through PC-13

| Finding | Status | Closure evidence |
|---|---|---|
| PC-01 multifile scientific agreement | **Closed** | Full box/cosmology/dialect/unit/config/governance comparison; generation/stem-aware discovery; exact contiguous indices; v2 strong set manifest. |
| PC-02 logical-set budgets | **Closed** | Cumulative preflight/merge budgets and actual merged-state memory accounting; member-count/gas/sidecar limits. |
| PC-03 star/tracer missing fields | **Closed** | Gas/star/tracer/BH fields use explicit reject/unavailable/documented-default policy; defaults/reconstructions reported. |
| PC-04 `evolution_ready` overclaim | **Closed** | Dedicated readiness logic proves persistent IDs, ownership and gas patch-or-parent resolution; analysis-only states remain explicit. |
| PC-05 zero persistent IDs | **Closed** | New persistence-boundary nonzero+unique validator is used by snapshot writing/readiness. |
| PC-06 Windows durability | **Closed in source / environment-limited runtime** | Unsupported `ReplaceFileW` flag removed; supported replace/move + explicit `FlushFileBuffers`; Windows CI transactional test added. |
| PC-07 independent validator | **Closed** | Direct bounded HDF5 validator checks set/header identity, dataset shape/type/counts, finite values, masses, bounds and global IDs without normal state import. |
| PC-08 hydro IC startup abort | **Closed** | Root patch is fully materialized before identity synchronization; `integration_initial_condition_manifest_runtime` passes locally. |
| PC-09 IC integrity modes | **Closed** | `kVerifiedIdentity` default uses one content hash plus stable identity; explicit `kStrictFullRehash` adds completion content rehash. |
| PC-10 collective density | **Closed at source contract / environment-limited MPI runtime** | Routing protocol v3 reduces successful per-batch calls 23 -> 20 while retaining exact reconciliation; MPI tests/CI carry runtime evidence. |
| PC-11 maintenance/schema drift | **Substantially closed for acceptance scope** | Validator/readiness are separate components; required native full-physics field contract is centralized and shared by writer/validator. Remaining large legacy units are tracked maintenance debt, not a known acceptance defect. |
| PC-12 strict-warning governance | **Closed** | `COSMOSIM_STRICT_IO_WARNINGS`; GCC+Clang I/O CI gate; local GCC strict build passes. |
| PC-13 MPI snapshot acceptance | **Closed in source / environment-limited local runtime** | Dedicated 2/4-rank snapshot-set tests registered and selected by MPI HDF5 CI; local host lacks MPI. |

## Previously partial IO findings

| Original finding | Status after closure | Evidence |
|---|---|---|
| IO-003 snapshot finalization | **Closed at supported source contract** | Shared transactional publication; POSIX local path validated; corrected Windows native path + Windows CI. |
| IO-004 restart finalization | **Closed at supported source contract** | Same transaction layer, durable POSIX file+directory sync, corrected Windows flush/replace path. |
| IO-008 silent missing scientific fields | **Closed** | Missing-field policy now covers full-physics star/tracer/BH/gas state instead of plausible silent defaults. |
| IO-009 hostile HDF5 budgets/validation | **Closed** | Single-file and cumulative logical-set budgets; direct structural/scientific validator; checked materialization accounting. |
| IO-014 collective storm | **Substantially closed / MPI evidence via CI** | 30 historical -> 23 -> 20 successful routing calls/batch without deleting exact audits. |
| IO-020 multifile snapshot import | **Closed** | Strict set discovery, scientific agreement, contiguous members, strong completion manifest, cumulative merge and direct validation. |
| IO-025 duplicated schema knowledge | **Substantially closed for active snapshot acceptance path** | Native required full-physics field contract centralized across creation/independent validation; set/readiness/conversion contracts isolated. |
| IO-026 oversized implementation units | **Substantially closed for acceptance risk** | New responsibilities remain extracted rather than returned to the monolith. Further file splitting is non-blocking maintenance debt. |
| IO-027 warning governance | **Closed** | Target-scoped strict I/O warnings and GCC/Clang CI; local GCC strict build passes; Windows native transaction CI added separately. |

## Acceptance interpretation

All concrete correctness gaps PC-01 through PC-09 are locally source/test closed where their
dependencies exist. PC-10 and PC-13 require MPI runtime evidence and PC-06 requires Windows
runtime evidence; the repository now contains dedicated CI gates rather than pretending those
platforms were executed on the Linux campaign host.
