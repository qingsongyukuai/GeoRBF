# Issue 25: Equality Spine acceptance and v0.1.0 release

Issue: [#25](https://github.com/qingsongyukuai/GeoRBF/issues/25)

Evidence seams: T01–T09, T11–T17 where applicable to the Equality Spine

Requirements: VAL-007, VAL-008, VAL-013–VAL-015

## Accepted product boundary

v0.1.0 is one complete, deliberately narrow public tracer. A GeoRBF user can
build hard absolute field-value and complete-gradient observations, unoriented
tangent directions, shared level sets or geological horizons, and an explicit
additive field gauge. The immutable snapshot lowers these inputs into the one
physical Canonical Problem IR, fits with the default Cubic contract and full
`Pi1` through the sequential faer augmented KKT, recovers and independently
verifies the candidate, and returns an immutable model plus typed report.

`examples/equality_spine.rs` is both the public runnable example and the body of
the T01 release integration test. Its manufactured planar stratigraphic field
combines every v0.1.0 input kind, recovers the horizon's semantic shared value,
checks all source and semantic-role provenance, and proves ordered batch values
and gradients equal the single-query path. Keeping one executable example as
the test body prevents the documented workflow from drifting into a second
product meaning.

## Cumulative corpus

The release workflow runs all existing product and contract tests rather than
selecting a new release-only subset. The accumulated evidence covers:

| Boundary | Release evidence |
| --- | --- |
| Public happy path | absolute affine field, combined planar geology, frame/metric covariance, shared latent and gauge invariance, tangent semantics |
| Builder and snapshot | checked leaves, atomic add, forward references, aggregate failure and repair, stable permutations, concurrent immutable snapshots |
| Diagnostics | direct and graph conflicts, uninformative groups, unidentified gauge and field modes, rank deficiency, gray zone, capacity priority, backend contract and recovery verification |
| Cubic/CPD/recovery | independent general/origin/generalized-functional fixtures, complete `Pi1`, reduced positivity, affine reproduction, reversible scaling and physical acceptance |
| Query | empty/repeated/large/invalid batches, first invalid index, chunk boundaries, concurrent clones, and a 100,000-location reduced-center smoke |
| Repository/release | rustdoc and compile tests, runnable example, dependency/license/native-link closure, traceability, package contents, and placeholder scan |

The source corpus contains no `todo!`, `unimplemented!`, string-form
not-implemented path, hidden Cargo feature, or fallback API for later
milestones. `scripts/release_audit.py` enforces that statement for product
source and verifies that the adopted oracle declarations, fixtures, and source
manifest remain byte-identical to the independent generator outputs.

## Property and platform profiles

Proptest is pinned at 1.11.0, whose MSRV is Rust 1.85 and whose fixed RNG seed
is configurable. The ordinary Linux x86-64 gate sets 256 cases for every
property family; each other supported native target sets 32. Every job uses a
checked-in fixed seed value. The nightly/release T12 job uses four fixed seed
ranges with 2,500 cases each, totaling 10,000 cases per family while retaining
independent shards for failure localization and promotion into the frozen
corpus.

The dev-only upgrade moves the test RNG stack to `rand` 0.9 and selects
`getrandom` 0.3.4. Its build script only reads rustc/target sanitizer and legacy
Windows cfg facts and emits Rust cfg values; it invokes no native compiler or
linker. The reviewed lockfile identity and this build-script version are pinned
by the same fail-closed audit used on every native target.

The native matrix is Linux x86-64, Linux AArch64, macOS x86-64, macOS Apple
Silicon, and Windows MSVC x86-64. Each native product job runs typechecking,
the cumulative tests, rustdoc/compile tests, strict Clippy, the public example,
repository audits, the feature graph, and package verification. The release
build is compiled with the optimized Cargo profile on every target. The release
profile additionally replays the dependency probe on all five native targets,
the Cubic KKT/QP/SOCP recovery probe, and both interpreter- and OCI-pinned
oracle pipelines. Ordinary CI timings remain smoke evidence, not a formal
performance promise.

## Traceability audit

`validation/v0.1.0/traceability.json` is the single scoped release map. It
lists each of the 80 requirement IDs reached by issues #16–#25 exactly once and
connects every requirement record to:

1. one named public API or Canonical IR behavior;
2. one or more repository test/oracle evidence markers; and
3. one or more user or implementation documentation markers.

The release audit owns the expected v0.1.0 requirement set independently. It
fails on missing, duplicate, or out-of-scope IDs; duplicate behaviors or API/IR
paths; empty API/IR paths; missing files; or stale evidence/documentation markers.
This gives T16 a mechanical dangling/duplication check while the issue #12
specification remains the only normative source of each requirement's wording.

## Risk, capacity, and package gates

Issues #13–#15 remain disposable probes and issues #16–#24 remain their product
adoption evidence. The release profile reruns those facts from the accepted
commit: native pure-Rust dependency closure, KKT/QP/SOCP recovery agreement,
and independent 120-digit oracle regeneration and consumption. The product
continues to enforce its 8 GiB fit peak before allocation and its 256 MiB query
scratch bound; the release adds no hidden normalization, regularization, or
timing waiver.

`cargo package --locked` must include the public example, release notes,
implementation evidence, adopted oracle corpus, and release traceability while
excluding throwaway product capabilities. The package is a source artifact for
real Rust integration feedback; this ticket does not authorize a crates.io
publication.

## Release procedure

1. Run the local replay below and review the complete diff from the issue #25
   starting commit.
2. Commit and push the accepted implementation to `main`.
3. Dispatch `.github/workflows/product-v0.1.yml` at that exact commit with
   `profile=release`; require every native, property, dependency, recovery,
   oracle, documentation, audit, and package job to pass with no waiver.
4. Create and push the annotated `v0.1.0` tag at that accepted commit. The tag
   triggers the same release profile because GitHub does not apply path filters
   to tag pushes.
5. After the tag-triggered run is green, create the GitHub release with
   `RELEASE_NOTES.md` and verify the release target equals the accepted commit.

Local replay under Rust 1.85.0:

```text
cargo fmt --check
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo test --locked --doc
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo run --locked --example equality_spine
python3 -m unittest discover -s scripts -p "test_*.py" -v
python3 scripts/release_audit.py
python3 scripts/audit.py
python3 spikes/oracle-fixtures/generate.py --check
python3 spikes/oracle-fixtures/verify.py
cargo package --locked
```
