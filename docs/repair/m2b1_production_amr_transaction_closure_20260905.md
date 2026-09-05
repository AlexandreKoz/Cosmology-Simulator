# M2B-1 Production AMR Transaction Closure

Date: 2026-09-05

This focused closure repairs the post-M2B acceptance findings in the production
`core::SimulationState` AMR regrid path.

- Refinement/derefinement stage only AMR-owned authority, not a complete simulation-state copy.
- The complete candidate/scratch high-water is admitted through the existing process `MemoryGovernor` before allocation or mutation.
- Pending reflux or temporal-boundary history blocks topology-changing regrid until synchronization is drained.
- Gas identity, pending reflux, and temporal-history allocated capacity are included in canonical memory accounting.
- Focused tests cover atomic governor rejection, synchronization-state rejection, memory-report participation, conservation, and population-independent regrid reservation scaling.
- The misleading star-formation MPI test registration that claimed AMR migration without executing it is removed. A genuine cross-rank AMR migration/restart continuation remains an MPI-capable follow-up acceptance item.

No hydro equation, prolongation/restriction method, reflux formula, CFL rule, refinement criterion, scientific precision, or tolerance is changed.
