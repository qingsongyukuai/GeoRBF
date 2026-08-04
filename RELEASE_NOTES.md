# GeoRBF v0.2.0 — Convex Relations

GeoRBF v0.2.0 is an end-to-end Convex Relations preview for real Rust
integration feedback. It keeps the complete v0.1.0 Equality Spine and adds
physical soft objectives and Cubic affine-inequality fitting through one
canonical KKT/QP pipeline.

## Supported scope

- One in-process Rust 2024 `georbf` crate for three-dimensional `f64` data,
  with Rust 1.85 as the exact preview MSRV and no public Cargo features.
- Hard and soft Field Value, complete Gradient, and Tangent observations;
  positive quadratic penalties and standard deviations; checked SPD
  covariance for vector residuals and named same-dimension Covariance Groups.
- Hard and soft Field Value Bounds, Directional Derivative Intervals,
  Younger/Older Horizon relations, non-strict Field Level Order, signed Field
  Separation Intervals, and finite Point-to-Level-Set side relations.
- Directed Normal and Axial Normal observations with independent direction and
  positive Minimum Normal Slope channels. Axial inputs require a separately
  identified explicit Polarity Resolution before fitting.
- Default and only kernel: scale-free Cubic with complete `Pi1`. Algebraic
  capability selects the shared faer augmented KKT or Clarabel Convex QP form;
  both use reversible transforms and physical Recover and Verify.
- Immutable `SolvedModel` values with simultaneous field value and complete
  gradient queries, ordered atomic logical batches, and recovered shared values
  by `GroupId`.
- Typed success/failure `FitReport` evidence for relations, residuals, slacks,
  losses, covariance group contribution, FieldEnergy, total objective,
  attempts, certificates, recovery, provenance, and canonical acceptance.

The runnable [Convex Relations example](examples/convex_relations.rs) combines
all of these public relation families in one manufactured affine field. It
audits the successful report and shared values, then checks ordered batch
queries against their single-query equivalents. Its companion release test
also runs a 512-scalar-constraint QP and 10,000-query quantity smoke.

## Compatibility boundary

This remains a `0.x` preview. A future `0.x` minor may adjust source APIs before
v1 freezes them. Patch releases will not silently change the published Cubic
mathematics, hard/soft relation meaning, diagnosis categories, stable identity
semantics, or the `georbf-v1` numerical policy; a mathematical contract fix
will be called out explicitly in release notes.

The crate owns its public geometry, relation, configuration, model, and
diagnostic types. It does not expose matrices, coefficients, polynomial bases,
Clarabel/faer enums, raw primal/dual vectors, solver settings, or implementation
extension traits. Backend termination remains evidence, not a domain diagnosis.

## Diagnostic semantics

Solved and almost-solved backend terminations produce candidates only. A model
is returned only after the same unrelaxed backend-standard-form and physical
canonical acceptance checks pass. Validated Farkas certificates and recession
rays are required before `InfeasibleProblem` or `Unbounded`; unverified status,
limit, insufficient progress, backend-contract failure, recovery corruption,
or inconsistent attempts produce a typed failure with no partial model.

Primary diagnosis follows the fixed semantic order: unresolved input semantics;
direct input conflict; unidentified gauge/mode or interpretable rank deficiency;
capacity; validated infeasible/unbounded; backend contract; recovery
verification; then other numerical/limit failure. Secondary proof and
provenance evidence is retained.

## Out of scope

v0.2.0 does not expose Directed Normal Cone, Angle, Cone Violation, vector
linear-norm epigraphs, SOC/SOCP product routes, Gaussian, Inverse Multiquadric,
Wendland C2, advanced functional inputs, physical-thickness or continuous-area
guarantees, CRS/unit conversion, Hessians, mesh extraction, persistence,
uncertainty estimation, or any v0.3.0+ placeholder.

The 512-constraint/10,000-query check and ordinary CI timings are smoke
evidence only. This release does not claim a formal reference-runner SLA,
freeze v1 workloads or `ValidationProfileId`, or publish to crates.io.

## Verification

The accepted commit must pass the zero-waiver `release` profile of
`.github/workflows/product-v0.2.yml` before tagging and again from the annotated
`v0.2.0` tag. That profile runs:

- the cumulative v0.1.0 and v0.2.0 public, mandatory, pairwise, higher-order,
  failure, certificate, recovery, capacity, query, concurrency, and oracle
  corpus on Linux x86-64/AArch64, macOS x86-64/Apple Silicon, and Windows MSVC;
- Linux property families at 256 cases, other platforms at 32 cases, and four
  fixed release seed ranges totaling 10,000 cases per property family;
- all locked dependency, Cubic CPD recovery, and independent-oracle risk
  replays, including repeated regeneration, adopted-fixture LF byte identity,
  and the digest-pinned OCI replay;
- Rust 1.85 typechecking, rustdoc, strict Clippy, optimized build, both public
  examples, dependency/feature/native-link/license audit, placeholder scan,
  traceability audit, and `cargo package --locked` verification.

The one-to-one scoped requirement map is
`validation/v0.2.0/traceability.json`. Release design and publication procedure
are recorded in
[`docs/implementation/36-convex-relations-release.md`](docs/implementation/36-convex-relations-release.md).
