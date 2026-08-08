# CHUI / CosmoSim C++ and schema naming conventions

**Version 2.0 — 2026-08-08**

This document is the repository-authoritative naming policy. It supersedes the
2025 naming document where that document conflicts with the mature C++ API or
the adopted `.param.txt` configuration system. This is a governance update,
not an ABI-renaming campaign.

## C++

- Namespaces: `lower_snake`.
- Public types: `PascalCase`.
- Functions and methods: `camelCase`.
- Local variables and parameters: `lower_snake`.
- Encapsulated private/protected data members: `m_lower_snake`.
- Public aggregate/record fields intentionally remain unprefixed
  `lower_snake`; adding `m_` to POD-style data contracts would obscure their
  record semantics and cause needless ABI/source churn.
- Scoped enum values use the established `kPascalCase` form.
- `constexpr`/compile-time C++ constants use `k_lower_snake` unless an external
  schema requires another spelling. Preprocessor macros remain
  `UPPER_SNAKE` with a `COSMOSIM_` prefix where project-owned.
- Template parameters use `T`, `TIndex`, `TPolicy`, or another concise
  PascalCase `T...` spelling.

## Units and coordinate frames

Raw numerical values must expose units/frame when ambiguity is possible.

- `_code`: code units.
- `_si`, `_cgs`: explicit physical unit systems.
- `_phys`: proper/physical frame.
- `_comoving`: comoving frame. This is the canonical repository spelling;
  legacy documentation that prescribed `_comov` is superseded.
- Unit-bearing domain names may combine both, for example
  `box_size_mpc_comoving` and `softening_kpc_comoving`.
- Cartesian SoA lanes use `_x`, `_y`, `_z`.
- Half-open ranges use `_begin`, `_end` or the established semantic equivalent.

Strong types are preferred when a unit error would be scientifically dangerous
and the type does not add hot-loop overhead.

## Configuration

`.param.txt` is the only authoritative simulation-input format. Keys and enum
values use `lower_snake`, except literature-standard cosmology identifiers that
are explicitly preserved by schema. Canonical booleans are only `true` and
`false`; legacy `yes/no/on/off/1/0` spellings are rejected rather than silently
normalized.

## HDF5 and external schemas

External compatibility names are not restyled. GADGET/AREPO-style groups,
datasets, and attributes keep their canonical spellings (`/Header`,
`Coordinates`, `ParticleIDs`, `NumPart_Total`, etc.).

## Concurrency and device state

Use semantic prefixes when they materially communicate ownership/lifetime:
`tls_` for thread-local state, `atomic_` for atomics, `mtx_` for mutexes, and
`d_`/`h_`/`p_` for device/host/pinned buffers where those distinctions are not
already encoded by the type.

## Governance

New code follows this document. Existing public names are changed only for a
correctness/unit/ownership reason or as part of an explicit compatibility
migration. Reviewers must not use naming cleanup as justification for broad
unrelated source churn.
