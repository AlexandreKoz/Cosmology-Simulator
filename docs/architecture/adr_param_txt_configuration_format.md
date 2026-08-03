# ADR: `.param.txt` is the authoritative simulation configuration format

**Status:** Accepted  
**Date:** 2026-08-02  
**Scope:** simulation input, examples, release profiles, normalized copies, and test fixtures

## Decision

CosmoSim formally uses the GADGET/AREPO-style plain-text parameter language
with the filename suffix `.param.txt` as its only authoritative simulation
configuration format. Filenames use lower snake case.

The runtime accepts the typed `key=value`/sectioned grammar documented in
`docs/configuration.md`. The normalized copy written into each run remains
`normalized_config.param.txt`.

YAML and JSON are not alternative simulation-input formats. Introducing either
would require a separate, explicitly approved compatibility-layer decision; the
core parser must not guess a format from a filename.

This decision does not prohibit JSON where it is already the correct generated
machine-readable representation, including provenance, capability summaries,
IC audit manifests, operational events, counters, and profiles.

## Rationale

The repository, examples, fixtures, parser, normalization path, and runtime
artifacts already converge on `.param.txt`. Formalizing that implementation
removes the stale design/naming-document claim that YAML or JSON is the intended
simulation language, avoids multiple partially equivalent parsers, and keeps
configuration hashing and provenance deterministic.

## Consequences

- New simulation decks and fixtures must end in `.param.txt`.
- Configuration diagnostics refer to parameter files, not YAML/JSON profiles.
- `configs/` contains canonical profiles plus explicitly checked compatibility
  entry points; see `configs/README.md`.
- JSON output schemas remain unchanged unless their own contract changes.
- Historical documents that describe YAML/JSON configuration are superseded by
  this ADR for the configuration-format question only.
