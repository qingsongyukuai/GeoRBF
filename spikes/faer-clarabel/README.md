# faer + Clarabel production-envelope probe

This throwaway crate is the replayable T15 probe for [GeoRBF issue #13](https://github.com/qingsongyukuai/GeoRBF/issues/13). It is deliberately outside the product crate tree and publishes nothing.

The manifest fixes Rust 2024/MSRV 1.85, `faer = 0.24.4` with only `linalg,std`, and `clarabel = 0.11.1` with only `serde`. Clarabel's `serde` feature is explicit because the published 0.11.1 source does not compile with every feature disabled. `Cargo.lock` freezes the complete transitive closure.

Replay on a native supported target:

```text
cargo test --locked --all-targets
cargo run --locked
python3 scripts/audit.py
cargo tree --locked --target all -e features
```

`scripts/audit.py` fails if the active target selects native-link metadata, known native build-tool packages, BLAS/LAPACK/PARDISO/SDP features, an unexpected backend feature graph, or a license outside the audited permissive set. It also prints every selected build script and license expression for the target.

The root workflow `.github/workflows/risk-spike-13.yml` replays the same commands on the five native target families required by NUM-002. A configured workflow is not itself platform proof: update the fact report only from completed run logs.
