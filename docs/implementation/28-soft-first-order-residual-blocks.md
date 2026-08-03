# Issue 28: Soft first-order residual blocks

Issue: [#28](https://github.com/qingsongyukuai/GeoRBF/issues/28)

Primary evidence seams: T01, T03–T05, T07–T11

Requirements: PAPI-008–PAPI-010; DOM-005, DOM-010, DOM-017–DOM-018; IR-003,
IR-005–IR-008, IR-014; NUM-005, NUM-009–NUM-010; DIA-006, DIA-008–DIA-009;
VAL-004–VAL-005, VAL-008, VAL-015

## Public residual boundaries

`GradientObservation` retains its hard constructor and adds only configurations
legal for the complete ordered three-component residual: an isotropic
`QuadraticPenalty`, an isotropic `StandardDeviation`, or an explicit checked
three-dimensional `CovarianceMatrix`. `TangentDirectionObservation` adds scalar
quadratic-penalty and standard-deviation constructors, while preserving exactly
one zero directional-derivative residual. Zero gradient continues to satisfy a
Tangent; neither polarity nor minimum-slope semantics is introduced.

`CovarianceMatrix` owns a finite, exactly symmetric, strictly positive-definite
matrix. Its fixed-array and dynamic-row constructors reject empty, non-square,
non-finite, asymmetric, and non-positive-definite inputs. A Gradient rejects a
matrix whose dimension is not three.

`CovarianceGroupBuilder` owns its members in explicit insertion order and
accepts Field Value, complete Gradient, or Tangent members. The first member
fixes the physical residual dimension; later members of another dimension and
duplicate member `SourceId` values are rejected without mutating the draft. A
completed group rejects an empty member set or a covariance dimension different
from the flattened residual dimension. `ProblemBuilder::add` then atomically
checks the complete group's `GroupId` and every member `SourceId` against all
other problem inputs. Because only the group builder can create group members,
hard relations and independently configured soft relations cannot be mixed into
the group or duplicated outside it.

## Canonical objective and numerical path

Canonical soft equality channels remain scalar physical functionals with stable
source, group, relation, residual, and role provenance. The canonical objective
now references ordered residual blocks rather than assuming one independent
loss per scalar row. A block stores one explicit loss and either a diagonal
isotropic precision or the full inverse covariance. Covariance cross terms enter
the symmetric faer KKT as `A^T P A` and `A^T P target`; the route therefore
remains the form-driven quadratic-objective-plus-affine-equalities path.

Covariance Cholesky whitening, its inverse, precision formation, coordinate
normalization, and Ruiz scaling remain private derived transformations.
Whitening does not define physical residuals, tolerance, or diagnosis, and
numerical scaling does not change covariance or objective meaning. The
gauge-invariant physical scale includes the largest covariance marginal standard
deviation or the equivalent penalty/statistical scale.

## Recover and Verify

Recovery independently evaluates every scalar functional in original units,
reassembles each ordered observation/member block, applies the derived whitening
map, and verifies that its inverse recovers the physical residual under the
fixed `1e-11` round-trip envelope. A damaged whitening map produces the public
`WhiteningRoundTripViolation` recovery reason and no model. Physical and
standard-form objectives are independently recomputed as

```text
1/2 FieldEnergy + 1/2 sum(block residual^T block precision block residual).
```

Independent Gradient and Tangent reports retain typed configuration and their
one block-level loss. A named covariance report retains `GroupId`, ordered
member `SourceId` and original scalar/vector residuals, covariance, whitened
residual, round-trip error, and exactly one group-level objective contribution.
It has no member-level objective accessor because covariance cross terms make
such a split non-identifiable.

## Evidence

`tests/public_soft_residual_blocks.rs` covers checked covariance and constructor
boundaries, atomic group construction and problem insertion, Euclidean vector
penalties under a non-identity anisotropy metric, single-observation covariance
cross terms, scalar Tangent semantics including zero gradient, named
cross-member covariance, contradictory soft evidence, public model queries,
and rotation/reflection plus uniform length/field-unit covariance.

`src/cubic_equality.rs` includes a narrow recovery corruption case proving that
a damaged whitening inverse is rejected without a model. The issue #27 scalar
soft corpus and the cumulative v0.1.0 public, canonical, failure, query, audit,
and packaging gates remain unchanged and passing.
