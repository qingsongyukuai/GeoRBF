# GeoRBF

GeoRBF is a Rust library for fitting implicit geological scalar fields from
geological observations.

Version 0.1.0 establishes the production equality/KKT execution spine and its
capacity evidence. The domain-facing fit API remains private until a later
milestone can expose an end-to-end supported capability without placeholders.

## Verification

The crate is pinned to Rust 1.85.0 and uses the checked-in lockfile:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
python3 scripts/audit.py
cargo package --locked
```

Implementation evidence for issue #16 is recorded in
[`docs/implementation/16-production-equality-spine.md`](docs/implementation/16-production-equality-spine.md).
