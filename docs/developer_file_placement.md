# Developer note: file placement and naming

## Naming conventions (authoritative)

- Files/directories: `lower_snake`
- Classes/types: `PascalCase`
- Methods: `camelCase`
- Variables: `lower_snake`
- Members: `m_lower_snake`
- Add explicit unit/frame suffixes when ambiguity exists.

## Placement checklist

For any new feature:

1. Pick exactly one owning module (`gravity`, `hydro`, `amr`, `physics`, `io`, `analysis`, `parallel`, `core`, `utils`).
2. Put public interfaces in `include/cosmosim/<owner>/`.
3. Put implementation in `src/<owner>/`.
4. Put internal-only headers in `src/<owner>/internal/`.
5. Add at least one unit test in `tests/unit/` and one integration/regression check in `tests/integration/`.
6. Add or extend a benchmark hook under `bench/` for hot paths.
7. Document config/schema/provenance impact under `docs/` and/or `configs/`.

## Anti-patterns to avoid

- Monolithic libraries with no ownership boundaries.
- Dumping physics/solver logic into `utils`.
- Public headers that include private internal headers.
- Placeholder code that pretends numerical completeness.
