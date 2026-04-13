# Units, cosmology, constants, and provenance

## Scope

This module centralizes foundational helpers shared by gravity, hydro, and I/O:

- Physical constants in one authoritative header.
- Unit-system conversions between code units and SI/CGS.
- LambdaCDM background helpers (`H(a)`, `rho_crit(a)`).
- Comoving/proper conversion utilities.
- Provenance serialization for run outputs.

## Assumptions

- `hubble_param` follows the standard `H0 = h * 100 km/s/Mpc` convention.
- Provenance records are authored by rank zero (`author_rank = 0`) and non-zero ranks skip file writes.
- Current hardware summary uses `std::thread::hardware_concurrency` as a portable baseline.
- Config hashes are deterministic FNV-1a over normalized config text.

## File format

Provenance is written to `provenance.meta.txt` in `key=value` lines with stable field order.
This format is intentionally lightweight and dependency-free for early pipeline bring-up.
