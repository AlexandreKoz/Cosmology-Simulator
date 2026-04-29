# Release checklist (v0.1.0-initial)

## Scope and claims

- [ ] Supported operating modes are listed (`zoom_in`, `cosmo_cube`, `isolated_galaxy`, `isolated_cluster`).
- [ ] Experimental modules are explicitly marked as non-production.
- [ ] Known issues document is present and linked.

## Build and packaging

- [ ] Preset-based builds are documented and reproducible.
- [ ] Project version and release manifest metadata are consistent.
- [ ] `release/release_manifest_v1.json` includes schema, naming-rule, and provenance policy versions.

## Configuration and examples

- [ ] Release example configs exist in `configs/release/`.
- [ ] Example configs parse through the standard config pipeline.
- [ ] Output naming stems in example configs are stable and explicit.

## Validation evidence

- [ ] Core unit tests pass.
- [ ] Integration smoke/reference workflow passes.
- [ ] Validation ladder (`unit`, `integration`, `regression`, `convergence`) passes.

## Benchmark evidence

- [ ] Release benchmark summary exists with environment disclosure guidance.
- [ ] Required benchmark hooks are buildable and runnable.
- [ ] Benchmark interpretations avoid correctness claims.

## Reproducibility and provenance

- [ ] Provenance schema (`provenance_v4`) is documented.
- [ ] Snapshot and restart schema versions are documented.
- [ ] Config normalization/hash policy is documented for reference runs.

## Final release gate

- [ ] `RELEASE_NOTES.md` links all release artifacts.
- [ ] All new release docs and config artifacts are covered by integration scaffold checks.
