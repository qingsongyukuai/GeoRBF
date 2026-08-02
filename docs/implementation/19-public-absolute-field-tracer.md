# Issue 19: public absolute field tracer

Issue: [#19](https://github.com/qingsongyukuai/GeoRBF/issues/19)

Evidence seams: T01, T03, T04, T11

Requirements: PAPI-001, PAPI-003–PAPI-008, PAPI-011–PAPI-014, PAPI-017;
DOM-001, DOM-005, DOM-017, DOM-021, DOM-022; IR-005, IR-008; VAL-015

## Public boundary

The `georbf` crate now exposes its first complete domain-facing vertical slice.
Every problem declares an ordered input coordinate frame, handedness, length
unit, and field unit. Crate-owned checked `Point3` and `Vector3` values prevent
non-finite spatial input from entering observations. Opaque owning `SourceId`
and `GroupId` values are independent of insertion order and backend layout.

`FieldValueObservation` represents one hard absolute scalar equality.
`GradientObservation` represents one hard complete-vector equality and lowers
to three scalar derivative functionals while retaining a common SourceId and
stable component role paths. `ProblemBuilder::add` atomically rejects duplicate
SourceIds; `build` owns and deterministically sorts the inputs, resolves the
default Cubic kernel, all-hard FieldEnergy normalization `1`, and the explicit
identity metric when no anisotropy was supplied. A failed build retains its
builder, while a successful `ProblemSnapshot` is immutable and cheaply cloned.

## Fit, recovery, and model

`ProblemSnapshot::fit` lowers only the supported hard observations into the
physical Canonical Problem IR and synchronously traverses the issue #18
Cubic/complete-Π₁/faer-KKT/Recover-and-Verify path. A `FitSuccess` owns both the
immutable `SolvedModel` and a typed `FitReport`. A `FitFailure` owns a typed
`ProblemDiagnosis` and the same report type; no public failure contains a
candidate or partial model.

The success report restores one assessment per scalar hard component with its
physical target, recovered value, residual, tolerance, SourceId, dimension,
and semantic role. It also records resolved kernel and normalization,
NumericalPolicyId, FieldEnergy, total objective, backend fingerprint, and
bounded attempt terminations. Backend termination remains distinct from the
problem diagnosis.

`SolvedModel::evaluate` returns one coherent `FieldSample` containing scalar
value and the complete gradient in the caller's declared input frame. The
model owns its snapshot, hides representer and polynomial coefficients, and is
cheaply cloneable and `Send + Sync` for repeated concurrent reads.

## Validation evidence

- T01 manufactures a known affine field from five absolute values and three
  complete gradients, fits it through the public crate, and checks queried
  value/gradient, hard residuals, zero Cubic quotient FieldEnergy, SourceIds,
  and component role provenance.
- T03 rejects non-finite coordinates, values, and gradient components; invalid
  frames; asymmetric, non-SPD, and non-det-one metrics; duplicate SourceIds;
  and empty builds. It proves rejected adds are atomic and failed builds are
  repairable.
- T04 is covered through public lowering and the recovered per-component
  provenance asserted by T01.
- T11 fits a non-affine hard problem under a translated, uniformly scaled,
  reflected orthogonal frame and the covariantly transformed determinant-one
  anisotropy metric. Recovered value and gradient obey the specified chain
  rule in both frames.

Logical batch queries, Shared Level Sets/Horizons, additive gauges, tangent or
normal observations, soft relations, and non-Cubic kernels remain outside this
tracer and have no premature public constructors.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
