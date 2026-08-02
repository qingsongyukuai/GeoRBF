# GeoRBF

GeoRBF is a Rust library for fitting implicit geological scalar fields from
geological observations.

Version 0.1.0 establishes the product-internal Cubic generalized-functional
Equality core, its versioned numerical policy, bounded KKT attempt path, and
physical Recover and Verify boundary. The domain-facing fit API remains private
until a later milestone can expose an end-to-end supported capability without
placeholders.

## Verification

The crate is pinned to Rust 1.85.0 and uses the checked-in lockfile:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
python3 scripts/audit.py
cargo package --locked
```

Implementation evidence is recorded for
[#16](docs/implementation/16-production-equality-spine.md) and
[#17](docs/implementation/17-cubic-equality-core.md), with the closed numerical
and recovery path recorded for
[#18](docs/implementation/18-cubic-equality-numerical-policy.md).
