# Issue 33: Field Separation and Point-to-Level-Set Relations

Issue: [#33](https://github.com/qingsongyukuai/GeoRBF/issues/33)

Primary evidence seams: PAPI-008–PAPI-010; DOM-003, DOM-014–DOM-020;
IR-004–IR-005, IR-008–IR-010, IR-014; NUM-009–NUM-010; DIA-002,
DIA-006–DIA-009; VAL-004, VAL-008, VAL-015

## Checked public semantics

`FieldSeparationInterval` names ordered `reference` and `target` Shared Level
Sets and constrains the signed quantity `target - reference` to one finite,
closed interval. Its constructor rejects self-reference, NaN, infinity, and an
empty interval while canonicalizing signed zero. Reversing the group roles
preserves meaning only with the explicit interval transform
`[lower, upper] -> [-upper, -lower]`; no age or spatial orientation is added.

`MinimumFieldOffset` accepts only finite, strictly positive field-value
quantities. `PointToLevelSetRelation` combines it with one finite `Point3`, one
referenced Shared Level Set, and an explicit `PointToLevelSetSide::Increasing`
or `Decreasing`. Side is never inferred from coordinates, gradients,
stratigraphic direction, names, IDs, or insertion order.

Both inputs retain caller-owned `SourceId` and `GroupId` provenance, allow
forward references, and enter the sealed `ProblemBuilder` boundary atomically.
Build reports every dangling reference in deterministic source/group order and
returns the builder for repair. Hard constructors carry no penalty; soft
constructors expose only checked quadratic or linear violation penalties and
require an explicit `FieldEnergyNormalization`.

## Canonical lowering, preflight, and recovery

Each Field Separation Interval side lowers to the physical semantic-latent
expression `target - reference` with its own lower or upper affine geometry.
Each Point to Level Set Relation lowers `point - level` as a lower bound for
Increasing Side or an upper bound for Decreasing Side. The latter form shares
the same field functional instead of inventing a sign-dependent representer
basis. All coefficients sum to zero across field and latent terms, preserving
the additive-gauge semantics.

Hard difference edges participate in deterministic graph preflight alongside
Shared Level Set member equalities and existing level relations. The proof is
derived from Canonical Problem IR and sums signed bounds as exact dyadic values
of the checked finite `f64` inputs, including absolute value and gauge anchors;
it never uses rounded or overflow-prone floating-point path accumulation. Such
positive cycles, including a strict point side at one of its own level-set
members, produce complete source/group/role evidence before backend execution.
Soft relations never enter those conflict proofs. Complete affine
infeasibility outside this exact graph boundary still requires a validated
backend certificate.

Any affine side selects the capability-driven Cubic QP route. Hard duplicates
may share canonical mathematics while retaining every source; soft sides keep
independent violation variables and objective terms. Recover and Verify
restores shared values, sampled point values, signed separations or oriented
offsets, physical slack/violation, tolerance, activity, configured penalty,
loss, FieldEnergy, total objective, and full provenance. Reports sort Field
Separation sides by `SourceId`/semantic role and Point to Level Set assessments
by `SourceId`.

## Units and guarantee boundary

Field Separation Interval bounds and Minimum Field Offset are field-value
quantities. A coordinate-frame or length-unit change does not alter them;
field-value rescaling multiplies their bounds, offsets, recovered values,
slacks, violations, and corresponding penalty semantics according to the
field-unit covariance rules.

GeoRBF verifies these relations only at their finite input support. They do not
claim a continuous-domain inequality, spatial above/below relation, distance,
mesh or isosurface result, stratigraphic age inference, or Physical Thickness.
No report field labels or converts a separation/offset as a length.

## Conformance evidence

`tests/public_level_set_affine_relations.rs` covers checked finite inputs,
self-reference, explicit side and positive offset, legal hard/soft APIs,
forward and dangling references, source atomicity, hard and soft Cubic QP
manufactured cases, independent quadratic/linear violations, active state,
objective and provenance, recovered shared values, exact graph conflicts,
signed role swapping, input permutation, coordinate-frame invariance, and
field-value rescaling.

The cumulative Field Value Bound, Directional Derivative Interval,
stratigraphic relation, QP recovery/certificate/capacity, audit, and package
suites remain part of this capability's release boundary.
