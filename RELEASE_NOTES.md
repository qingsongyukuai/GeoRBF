# GeoRBF v0.1.0 — Equality Spine

GeoRBF v0.1.0 is the first end-to-end preview for Rust integration feedback.
It turns supported geological observations into one immutable implicit scalar
field through a fully verified Cubic Equality path.

## Supported scope

- One in-process Rust 2024 `georbf` crate for three-dimensional `f64` data,
  with Rust 1.85 as the exact preview MSRV.
- Explicit orthogonal input frame, handedness, length-unit label, field-unit
  label, stable `SourceId`, and stable `GroupId` inputs.
- Hard field-value observations, complete gradient observations, unoriented
  tangent directions, shared level sets, geological horizons, and explicit
  additive field gauges.
- Default and only fitting route: the scale-free Cubic kernel with complete
  `Pi1`, a sequential faer augmented KKT, reversible numerical transforms, and
  physical Recover and Verify.
- Immutable `SolvedModel` values with simultaneous field value and complete
  gradient queries, including ordered atomic logical batches and recovered
  shared values by `GroupId`.
- A typed `FitReport` on success or failure with stable provenance, problem
  sizes, numerical evidence, backend attempts, physical residuals, and
  canonical acceptance where reached.

The runnable [Equality Spine example](examples/equality_spine.rs) combines an
absolute field value, complete gradient, tangent direction, planar horizon,
explicit horizon gauge, shared-value recovery, provenance inspection, and
single/batch queries through the public API.

## Compatibility boundary

This is a `0.x` preview. The supported public surface is intentionally narrow:
one crate, one 3D `f64` field, no public Cargo features, Cubic only, all-hard
Equality only, and sequential fitting. The crate owns its public geometry and
diagnostic types and does not expose matrices, coefficients, polynomial bases,
backend variables, solver traits, or extension traits.

Patch releases will not silently change the mathematical Cubic contract,
diagnosis meanings, stable identity semantics, or accepted numerical policy.
As a pre-1.0 release, broader source compatibility is not yet promised.

## Diagnostic semantics

Fit success means a backend candidate passed both the adapter contract and
independent recovery in physical canonical semantics. A backend termination is
only evidence; it is never itself a `ProblemDiagnosis`. Failures return no
partial, candidate, or best-effort model.

The Equality preview distinguishes invalid builder input, uninformative shared
level sets, direct hard conflicts, an unidentified additive gauge or field
mode, checked capacity exhaustion, backend contract violation, recovery
verification failure, numerical decision gray zones, and other numerical or
limit failures. General certificate-backed convex infeasibility is not claimed
in this release.

## Out of scope

v0.1.0 does not expose soft relations, statistical covariance, affine bounds
or ordering, field-separation intervals, point-side relations, directed or
axial normals, normal cones, QP/SOCP product routes, Gaussian, inverse
multiquadric, or Wendland kernels, custom functionals, CRS conversion,
Hessians, mesh extraction, persistence formats, uncertainty estimation, or a
formal reference-runner performance SLA. Those capabilities remain reserved
for later milestones and have no placeholder API in this release.

## Verification

The accepted commit must pass the `release` profile of
`.github/workflows/product-v0.1.yml` before tagging and again from the
`v0.1.0` tag. That profile runs:

- the cumulative public, canonical, failure, capacity, query, concurrency, and
  oracle consumer corpus on Linux x86-64/AArch64, macOS x86-64/Apple Silicon,
  and Windows MSVC under Rust 1.85.0;
- four fixed property seed ranges totaling 10,000 cases per T12 family;
- all three locked risk-spike replays, the 120-digit oracle regeneration twice,
  adopted-fixture byte identity, and the digest-pinned OCI replay;
- rustdoc/compile tests, the public example, strict Clippy, dependency/feature/
  native-link/license audits, traceability, `cargo package --locked`, and the
  optimized release build plus the 100,000-query reduced-center smoke.

The scoped requirement-to-behavior/evidence/documentation map is
`validation/v0.1.0/traceability.json`. Formal v1 ValidationProfile workloads
and performance claims are deliberately not frozen by this preview.
