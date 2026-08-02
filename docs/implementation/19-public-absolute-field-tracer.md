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
The admitted backend is sequential: `ThreadBudget::Exact(1)` is honored, while
larger exact requests are rejected during build instead of being silently
reported without changing execution.

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
independent input, scalar-relation, latent, auxiliary, cone, primal, equality,
and KKT dimensions. Exact duplicate scalar facts are merged before
representation so they cannot create artificial multiplier nullity, while an
assessment for every original SourceId remains in the report. Every bounded
attempt retains its algorithm, complete settings, scaling summary, refinement
count, residual evidence, certificate presence, rejection reason, and full
backend fingerprint. The report also retains Cubic polynomial/reduced-pairing
condition evidence, backend rank and inertia, and the full physical canonical
acceptance envelope. Backend termination remains distinct from the problem
diagnosis. Failed backend and recovery paths retain the same available attempt,
rank/inertia/capacity, per-source relation, side-condition, and physical
rejection evidence. Directly contradictory exact values or gradient components
at one point are diagnosed before backend execution with both stable SourceIds,
the semantic component, and incompatible targets.

Lowering derives stable internal RelationId and ResidualId values from each
SourceId and semantic role. Equality assembly adds stable derived block, row,
and representer-column identities alongside backend indices, and Recover and
Verify checks that the complete source → relation → residual → derived-artifact
association survives assembly before a model can exist. Cone identity is not
fabricated because this ticket admits no conic relation.

`SolvedModel::evaluate` returns one coherent `FieldSample` containing scalar
value and the complete gradient in the caller's declared input frame. The
model owns its snapshot, hides representer and polynomial coefficients, and is
cheaply cloneable and `Send + Sync` for repeated concurrent reads. Its custom
`Debug` representation is intentionally redacted so formatting cannot expose
the hidden coefficients or polynomial.

## Validation evidence

- T01 manufactures a known affine field from five absolute values and three
  complete gradients, fits it through the public crate, and checks queried
  value/gradient, hard residuals, zero Cubic quotient FieldEnergy, SourceIds,
  and component role provenance.
- T03 rejects non-finite coordinates, values, and gradient components; invalid
  frames; asymmetric, non-SPD, and non-det-one metrics; duplicate SourceIds;
  empty builds; and unsupported exact thread counts. It proves rejected adds
  are atomic, failed builds are repairable, and contradictory hard inputs
  return typed source evidence without starting a solver attempt.
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
