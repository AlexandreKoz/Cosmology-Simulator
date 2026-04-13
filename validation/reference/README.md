# Validation Reference Data Policy

- **Versioning:** golden/reference policy files are versioned by explicit suffix (`*_v1`, `*_v2`, ...).
- **Stability:** do not edit an existing version in place; create a new version and update consumers deliberately.
- **Scope:** small CI-scale validation artifacts live in-repo; larger datasets should be stored externally with checksums and retrieval metadata.
- **Tolerance intent:** tolerance envelopes are for correctness validation, not performance assessment.
- **Schema discipline:** any key rename or semantics change requires a schema/version bump and docs update.
