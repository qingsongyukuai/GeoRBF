# Issue 30: Field Value Bounds

Issue: [#30](https://github.com/qingsongyukuai/GeoRBF/issues/30)

Primary evidence seams: PAPI-008, PAPI-010, PAPI-013; DOM-004,
DOM-017–DOM-020; IR-005, IR-008–IR-010, IR-014; NUM-009–NUM-014;
DIA-002, DIA-004–DIA-009; VAL-004–VAL-005, VAL-008, VAL-015

## Checked public relation

`FieldValueBound` is an atomic problem input with checked lower, upper, and
closed-interval constructors. Endpoints must be finite, accepted signed zero is
canonicalized to positive zero, and an interval whose lower endpoint exceeds
its upper endpoint is rejected as `FieldValueBoundError::EmptyInterval`.
Hard constructors have no weight, user slack, or tolerance parameter.

Soft constructors create one nonnegative violation channel per side. A channel
uses either a checked positive `QuadraticPenalty` or checked positive
`LinearViolationPenalty`; `FieldValueViolationPenalty` lets the two sides of an
interval choose those loss families independently. Each side retains a stable
residual identity derived from its `SourceId` and `field-value-bound/lower` or
`field-value-bound/upper` semantic role.

The builder enforces global `SourceId` uniqueness atomically. Immutable
snapshots retain bounds separately from observations, expose their count, and
require an explicit `FieldEnergyNormalization` whenever any side is soft.

## Canonical lowering and QP realization

Bounds lower into the one physical `CubicCanonicalProblem` as affine
inequalities. Lower sides are normalized to upper-form backend rows by a sign
change; the canonical relation itself retains its physical sense and bound.
Identical hard sides share one canonical/backend relation while all source
relations remain attached to that relation through complete source-to-canonical
and canonical-to-backend provenance edges. Soft sides are never merged, so
duplicate evidence contributes independent violation variables and loss terms;
their augmented rows and nonnegativity rows each retain their source, derived
block/row/column, backend row/column, and cone role.

The capability-driven executor from issue 29 selects Clarabel QP whenever the
canonical problem contains a bound. Every soft side adds exactly one explicit
violation variable, one augmented bound row, and one nonnegativity row.
Quadratic losses enter the Hessian diagonal and linear losses enter the linear
objective. Checked QP capacity includes those variables and rows before form
allocation or backend entry.

## Recovery, reporting, and diagnoses

Recovery reverses deterministic QP scaling and the Cubic reduction, evaluates
each relation in physical field units, and independently checks hard
feasibility, soft augmented-row feasibility, violation nonnegativity, backend
slack equations, provenance, side conditions, and objective round trips.
`FitReport::field_value_bounds` is sorted by `SourceId` and semantic role and
reports the physical bound, recovered value, satisfaction slack, violation,
activity, configured loss, and per-side loss contribution. The shared report
also retains `FieldEnergy` and total objective.

Same-support hard lower/upper contradictions and exact field-value/bound
contradictions are proven during canonical preflight and reported as
`DirectInputConflict` without backend attempts. Soft contradictions are not
classified as infeasible.

General infeasibility is reported only after a Clarabel primal-infeasibility
candidate passes an independent normalized Farkas check: all quantities are
finite, the ray has unit infinity norm, `A^T z` and dual-cone violations are
within the fixed numerical-policy limit, and normalized `-b^T z` exceeds the
fixed strict-separation margin. Invalid or absent rays remain numerical
failures. Valid evidence is exposed as
`FitReport::infeasibility_certificate`; no failed fit can publish a model.
Every legal Field Value Bound objective is a sum of nonnegative FieldEnergy and
positive-weight violation losses, so it cannot have a recession direction with
negative objective. Consequently a backend dual-infeasibility candidate is not
promoted to an `UnboundedProblem` diagnosis for this capability; it remains an
unverified numerical failure.

## Conformance evidence

Public tests cover checked constructors, hard lower/upper/interval fits, both
soft loss families, mixed interval losses, duplicate semantics, stable reports,
direct and certified general infeasibility, frame and field-unit covariance,
and input permutation. Contract tests inject QP recovery corruption and verify
structured rejection with no model. The issue 29 capacity and corruption tests,
issue 27 objective covariance tests, and the cumulative workspace suite remain
part of the acceptance boundary.
