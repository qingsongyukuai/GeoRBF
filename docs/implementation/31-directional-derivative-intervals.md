# Issue 31: Directional Derivative Intervals

Issue: [#31](https://github.com/qingsongyukuai/GeoRBF/issues/31)

Primary evidence seams: PAPI-010, PAPI-013; DOM-010, DOM-017–DOM-018,
DOM-020–DOM-022; IR-002–IR-003, IR-005, IR-007–IR-010, IR-014;
NUM-009–NUM-010; DIA-002, DIA-006–DIA-009; VAL-004–VAL-005,
VAL-008, VAL-015

## Checked public relation

`DirectionalDerivativeInterval` is an atomic problem input with checked lower,
upper, and closed-interval constructors. It stores one explicit oriented
direction, normalized with the physical input frame's Euclidean metric by a
scaled algorithm that remains defined for finite subnormal and near-maximum
components. The zero vector, non-finite endpoints, and an interval whose lower
endpoint exceeds its upper endpoint are rejected. Accepted signed-zero
endpoints and direction components are canonicalized to positive zero.

Bounds and recovered quantities use field-value-per-length units. The relation
does not introduce an angle, a Tangent tolerance, a normal cone, a complete
gradient magnitude, or a user numerical tolerance. A zero-width hard interval
at zero and a hard `TangentDirectionObservation` have equivalent field
observables, while remaining distinct public domain types with distinct report
roles.

Hard constructors carry no loss or slack configuration. Soft lower and upper
sides each own an independent nonnegative violation channel and select either a
checked positive `QuadraticPenalty` or `LinearViolationPenalty` through the
relation-specific `DirectionalDerivativeViolationPenalty` entry point. A soft
relation requires explicit `FieldEnergyNormalization` at snapshot build time.

## Canonical lowering and Convex QP execution

Each side lowers into the sole physical `CubicCanonicalProblem` as a
field-value-per-length affine inequality. Its field functional is the
contraction of the stored unit direction with the gradient at the finite
support. Stable relation and residual identities derive from `SourceId` and the
`directional-derivative-interval/lower` or
`directional-derivative-interval/upper` semantic role.

Exact hard identity uses the complete normalized functional, sense, and bound,
not support alone. Positively scaled equivalent directions therefore share one
canonical/backend hard row while retaining every source provenance; soft
duplicates remain independent objective evidence. Exact hard lower/upper or
Tangent/interval contradictions are diagnosed during preflight. Other affine
infeasibility retains issue 30's independently verified Farkas-certificate
boundary, and an unverified ray remains a numerical failure with no model.

The algebraic capability selects the issue 29 Clarabel QP path. Every soft side
adds one explicit violation variable, one augmented affine row, and one
nonnegativity row. Derivative quadratic and linear penalty reference scales are
converted to field-value scale with the physical solve length before they enter
the gauge-invariant characteristic scale. Backend stopping targets retain
issue 30's established Standard/Robust policy; the unchanged public convex
residual envelope verifies both sides before a zero-width interval is accepted.

## Recovery, reporting, and covariance

Recover and Verify independently evaluates each directional functional from
the recovered physical field, restores derivative-unit slack and violation,
checks hard feasibility, violation nonnegativity, backend slack equations,
objective and transformation round trips, and verifies the complete provenance
graph before publishing a model.

`FitReport::directional_derivative_intervals` is sorted by `SourceId` and
semantic role. Every side reports its normalized physical direction, bound,
recovered directional derivative, satisfaction slack, violation, tolerance,
active state, configured loss, and optional objective contribution. Problem
sizes include every source side, deduplicated canonical hard row, soft
auxiliary variable, affine/nonnegative row, and quadratic or linear objective
term.

For a frame similarity `x' = sQx + b`, the direction becomes `Qd`, gradients
and derivative quantities become the field-unit scale divided by `s` times
their rotated values, quadratic penalty weights transform by the inverse square
of that derivative scale, and linear weights by its inverse. Tolerances and
diagnostic quantities follow the same chain rule. The determinant-one Global
Anisotropy Metric covaries only in kernel distance and is never used to
normalize the physical direction or residual.

## Conformance evidence

`tests/public_directional_derivative_intervals.rs` covers checked construction,
typed soft configuration, atomic builder insertion, hard and soft public QP
fits, manufactured affine and quadratic data, independent side losses,
Tangent equivalence, exact direct conflicts, derivative-only gauge provenance,
validated general infeasibility, positive direction scaling and input
permutation, stable report ordering, and rotation/reflection/uniform-scale
covariance with non-identity anisotropy.

Narrow contract tests inject QP recovery-map corruption, a checked capacity
rejection, and an invalid derivative-problem Farkas ray through public
snapshots. They verify structured failures, retained evidence, and no model.
Issue 29's backend/corruption corpus and the complete pre-existing public,
canonical, audit, and package suites remain part of this capability's release
boundary.
