# M2D owner-lifetime continuation — 9 September 2026

Status: partial coding closure. This document records source changes made against the supplied post-merge archive, not historical claims of completion.

## Implemented

The AMR temporal-history owner admits the exact outer and nested record capacities before construction, accounts for the maximum temporary geometry and its construction scratch, and commits a candidate through a no-throw replacement. The old history remains intact on failure. End capture updates fixed-size records without allocating. The existing topology-quiescence barrier remains authoritative.

Pending flux-register merging constructs and validates a bounded replacement before replacing the retained store. Its upper bound is `(old_count + incoming_count) * sizeof(PendingFluxRegisterRecord)`, and actual vector capacity is reconciled with the existing MemoryGovernor. The original record merge order, coverage rules and flux arithmetic remain unchanged. This safe replacement may cost an additional full pending-store copy; it is not advertised as a memory-saving algorithm.

Local subcycling no longer duplicates two full-population density arrays. The identical density conversion factor is applied to the already admitted patch-local source fields. The public source context remains a borrowed view, and no physical model or integration tolerance changes.

The restart integrity traversal now sorts non-owning temporal-history indexes
in a checked, governed arena. The output writer passes its existing governor
through the backward-compatible write policy. No serialized field, version or
hash algorithm changes. The exact incremental bound is P pointer widths plus
the largest cell-index array plus 256 bytes of arena alignment. This removes
an avoidable copy of the full conserved temporal history, but does not certify
the remaining writer or reader allocations.

## Evidence

The focused CPU build and four AMR tests pass. The CPU restart-schema test also passes, including governed integrity-index rejection and retry. New tests exercise tight-headroom rejection before mutation, retained baseline reconciliation and deterministic retry. Existing temporal interpolation, pending-register behavior, subcycling and reflux conservation remain passing. The CPU restart-schema test passes, including governed integrity-index rejection and retry. Repository hygiene passes in the clean patch-verified extraction. No actual RSS reduction or representative 48 GiB envelope is claimed.

## Remaining exact boundaries

The distributed remote flux destination and MPI exchange coexistence, complete reflux-materialization scratch, active-level selection, regrid candidate/communication coexistence, source event/report retained reconciliation and complete restart readback ownership require further implementation. Full MPI/FFTW and whole-process acceptance remain unavailable in this environment. No overlap or complete major-task certificate is enabled by this patch.
