# Risk spike 13: faer + Clarabel production dependency envelope

Issue: [#13](https://github.com/qingsongyukuai/GeoRBF/issues/13)

Evidence seam: T15

Requirements: PAPI-002, PAPI-012, NUM-001, NUM-002, NUM-012, NUM-013, VAL-013

Probe: [`spikes/faer-clarabel`](../../spikes/faer-clarabel/README.md)

## Verdict

**Production admission remains gated.** The Linux x86-64 probe proves that the pinned pure-Rust routes can supply the required dense KKT, rank/factorization, QP, SOCP, candidate, residual, certificate, settings, and exact-one-thread evidence. The other four required native platforms are still **ambiguous** until the committed native matrix has completed. Independently, `faer` factor-workspace capacity remains **ambiguous** until its workspace behavior is measured against the 8 GiB plan. Under #13, either ambiguity is unproven and cannot release the product implementation gate.

No product dependency has been added. The probe is a standalone, unpublished Rust crate outside the future `georbf` product tree.

## Fixed envelope

| Field | Fixed value | State |
| --- | --- | --- |
| Rust edition | 2024 | proven |
| MSRV/toolchain | Rust 1.85.0 | proven on Linux x86-64 |
| `faer` | `=0.24.4`, defaults off, `linalg,std` | proven on Linux x86-64 |
| `clarabel` | `=0.11.1`, defaults off, `serde` | proven on Linux x86-64 |
| Clarabel linear solver | `qdldl` | proven; one actual thread |
| Lockfile | Cargo lock format 4; SHA-256 `16250c7d3102c11f9613555623049b48dea7dcda02e533fc0389d4ed148ebe0c` | proven |
| Native/BLAS/LAPACK/PARDISO/SDP | no selected package, link, or feature on Linux x86-64 | proven |
| License policy | MIT/Apache-2.0-compatible permissive closure | proven on Linux x86-64 |

Clarabel 0.11.1 is **unavailable** with a truly empty feature set: its published source references `serde_json` from an unguarded error variant. The envelope therefore explicitly selects only `serde`; this adds serialization crates but no native dependency.

## Capability facts

| Capability | State | Evidence and boundary |
| --- | --- | --- |
| Symmetric indefinite KKT solve | proven | `faer` LBLT solved the manufactured 3×3 KKT candidate `[0.5, -0.5, 1.0]`; normalized backward error was zero on the observed run. |
| Inertia | proven with adapter work | `faer` exposes LBLT 1×1/2×2 block diagonals, from which the probe independently derived inertia `2/1/0`. There is no first-class inertia verdict; GeoRBF must own thresholds and classification. |
| Column-pivoted QR/rank evidence | proven with adapter policy | The rank-one 3×2 matrix reported one accepted diagonal under the probe threshold. `faer` supplies factors/pivots, not GeoRBF rank policy. |
| SVD/rank evidence | proven with adapter policy | Singular values were exposed and independently classified as rank one. Threshold/gray-zone policy remains GeoRBF-owned. |
| Cholesky/error behavior | proven | SPD input succeeded and symmetric-indefinite input returned an LLT error. Cholesky failure alone is not a diagnosis. |
| Capacity behavior | ambiguous | Checked arithmetic rejects both a representable 32,769² `f64` plan just above 8 GiB and an overflowing plan before allocation. The probe does not exercise `faer` allocation/factor-workspace behavior, so GeoRBF must still measure and budget all factor workspaces before production admission. |
| Convex QP | proven | Clarabel returned primal, dual, slack, `Solved`, residuals, gap, iterations, settings, and linear-solver information for the manufactured QP. |
| SOCP/SOC block | proven | Clarabel returned the expected `[1, 1]` candidate for one three-dimensional SOC block with residual/gap evidence. |
| Primal infeasibility | proven | Clarabel returned `PrimalInfeasible` and a dual-cone ray; the probe independently checked `Aᵀz≈0` and `−bᵀz>0`. |
| Unboundedness | proven | Clarabel returned `DualInfeasible` and a recession direction; the probe independently checked the homogeneous constraint and descent margin. |
| Configurable settings and scaling | proven | The probe fixes and records feasibility tolerance, equilibration, regularization/refinement profile, direct solver, iteration limit, and max threads; every Clarabel attempt also emits the actual variable, constraint, inverse, and objective equilibration factors. |
| `ThreadBudget::Exact(1)` | proven for this envelope | `faer` is built without Rayon and reports sequential execution; Clarabel `qdldl` reports one actual thread for `max_threads=1`. No global thread setter or environment variable is used. |
| Automatic/multithreaded route | unavailable in this minimal envelope | Deliberately not selected. A future Rayon or Clarabel faer-sparse envelope needs a new probe of per-call thread control and its larger closure. |

## Dependency and license audit

On native Linux x86-64, Cargo selected 80 packages including the probe. The active graph had no package with Cargo `links` metadata and none of the forbidden native/toolchain packages or BLAS, LAPACK, MKL, Netlib, PARDISO, OpenBLAS, or SDP features. The exact versions and checksums are in `Cargo.lock`; `cargo tree --locked --target all -e features` emits the complete feature graph.

Twenty selected crates have build scripts: `clarabel`, `libc`, `libm`, `nano-gemm-{c32,c64,f32,f64}`, `num-traits`, `paste`, `private-gemm-x86`, `proc-macro2`, `pulp`, `quote`, `serde`, `serde_core`, `serde_json`, `syn 1`, `thiserror`, `zerocopy`, and `zmij`. Inspection found Rust cfg/version detection or Rust code generation only. Clarabel's script contains dormant BLAS/SDP configuration branches, but the feature audit proves none is selected. No `cc`, CMake, bindgen, pkg-config, vcpkg, BLAS/LAPACK provider, or native `links` package is active.

The active non-product license expressions are MIT, Apache-2.0, MIT/Apache alternatives, BSD-2-Clause, BSD-3-Clause, Zlib, Unlicense, and Unicode-3.0 in conjunction with MIT/Apache. They are permissive and compatible with GeoRBF's MIT/Apache-2.0 baseline. `scripts/audit.py` freezes the accepted expressions and fails closed on drift.

## Native platform matrix

| Required native platform | Runner/target | State | Evidence |
| --- | --- | --- | --- |
| Linux x86-64 | local `x86_64-unknown-linux-gnu` | proven | Rust 1.85 check, focused tests, full probe, and audit completed. |
| Linux AArch64 | `ubuntu-24.04-arm` / `aarch64-unknown-linux-gnu` | ambiguous | Native job committed; no completed run attached yet. |
| macOS x86-64 | `macos-15-intel` / `x86_64-apple-darwin` | ambiguous | Native job committed; no completed run attached yet. |
| macOS Apple Silicon | `macos-15` / `aarch64-apple-darwin` | ambiguous | Native job committed; no completed run attached yet. |
| Windows MSVC x86-64 | `windows-2025` / `x86_64-pc-windows-msvc` | ambiguous | Native job committed; no completed run attached yet. |

The labels follow GitHub's current hosted-runner matrix. Workflow configuration is replay machinery, not substitute evidence for a completed native run.

## BackendFingerprint fields

Every replay emits the following implementation evidence, separate from NumericalPolicy semantics:

- schema version; `rustc` version; native target triple; lockfile identity;
- backend crate name/version and selected features;
- Clarabel direct linear solver and cone families;
- resolved settings: tolerances, equilibration, regularization/refinement, iteration limit, direct solver, requested threads;
- actual adapter/backend scaling factors, actual threads, termination, primal/dual/slack or certificate candidate, residuals, gaps, iterations, and optional failure reason;
- dependency audit result, active package count, build scripts, native links/features, and license expressions.

Setup, factorization, non-finite-candidate, and unexpected-termination failures use structured `ProbeError` evidence containing backend, version, problem class, typed failure reason, and detail; they do not collapse to an unclassified error string.

## Replay

From `spikes/faer-clarabel`:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo run --locked
python3 scripts/audit.py
cargo tree --locked --target all -e features
```

The workflow repeats these commands under Rust 1.85.0 on every required native runner. Completed jobs must replace each ambiguous platform row with a run URL and observed fingerprint. Production admission additionally requires measured `faer` factor-workspace capacity evidence; completing the platform matrix alone is insufficient.

## Primary sources

- [`faer 0.24.4` documentation](https://docs.rs/faer/0.24.4/faer/)
- [`Clarabel 0.11.1` settings](https://docs.rs/clarabel/0.11.1/clarabel/solver/implementations/default/struct.DefaultSettings.html)
- [`Clarabel 0.11.1` solver information](https://docs.rs/clarabel/0.11.1/clarabel/solver/implementations/default/struct.DefaultInfo.html)
- [GitHub-hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners)
