# Issue 20: shared level sets, horizons, and additive field gauge

Issue: [#20](https://github.com/qingsongyukuai/GeoRBF/issues/20)

Evidence seams: T01, T03, T04, T06, T08

Requirements: PAPI-009, PAPI-015; DOM-002, DOM-003, DOM-019; IR-004,
IR-005, IR-008

## Public groups and gauges

`SharedLevelSetBuilder` and `HorizonBuilder` atomically own a stable `GroupId`
and complete members, each with its own `SourceId` and checked `Point3`.
Finishing an empty builder fails, duplicate member sources are rejected without
mutating the draft, and the immutable completed group is the only value that
can enter `ProblemBuilder`. A group cannot receive members through the problem
after insertion. General shared level sets and geological horizons are distinct
public types but lower through the same mathematical shared-level semantics.

`ProblemBuilder::add` atomically checks group and source identities. Gauge and
future relation inputs may forward-reference a `GroupId`; `build` aggregates
unresolved references deterministically while retaining the builder for repair.
The owning snapshot keeps the original group kind, member provenance, gauge
reference, caller frame, and field-unit label.

`AdditiveFieldGauge::at_point` and `at_level_set` are checked, hard conventions
with independent `SourceId`s. They are not field-value observations and retain
distinct semantic role paths in recovery reports. Non-finite representatives
are rejected at construction. No automatic first-member anchor, constant-term
zeroing, or backend minimum-norm representative is used.

## Canonical latent and derived equality form

Lowering materializes one stable semantic latent record per group. The record
owns its `GroupId`, field unit, and every member `SourceId` and support. Thus the
canonical shared value remains present independently of input ordering and of
the algebraic form selected below.

The current all-hard Cubic path derives a deterministic equality realization:

- an explicitly gauged level set lowers each member to the selected recovered
  representative;
- an ungauged multi-member level set lowers to a deterministic full-rank set of
  Helmert-style mean contrasts in stable SourceId order;
- a point gauge lowers as a distinct absolute convention functional.

No member is assigned a value or selected as a semantic reference, so the
contrast realization introduces no first-point anchor. Exact equalities are
normalized, merged, and retain all original public relation assessments.
Opposite-orientation duplicate contrasts share one backend row. This form is a
recovery-mapped optimization only: after solving,
GeoRBF evaluates every original member and recovers the latent as the stable
mean of those equal values. The report independently verifies each canonical
member residual against that recovered latent and each level-set gauge residual
against its declared representative.

The recovered `SharedLevelValue` retains `GroupId`, field unit, and all member
sources. `SolvedModel::shared_level_value` exposes the same immutable value by
`GroupId`. `HardRelationAssessment::group_id` completes source/group provenance
for member and gauge relations without exposing a latent column, coefficient,
matrix, or backend variable.

## Preflight and physical verification

A group whose distinct support contains only one point and whose latent is not
referenced produces typed `UninformativeSharedLevelSet` evidence before backend
execution. A level-set gauge is a real latent reference, so a gauged singleton
is not misreported. If every relation is invariant to `f -> f + c` and there is
neither a true field-value observation nor an explicit gauge, fit returns typed
`UnidentifiedAdditiveGauge` evidence with stably ordered source and group IDs,
again before backend execution.

Manufactured affine horizons pass through the default Cubic/full-`Pi1`/faer KKT
path. Recovery verifies member equalities, the gauge, full gradient, field value,
zero affine Cubic quotient energy, and complete provenance. Changing only the
gauge shifts queried and shared field values by the same constant while leaving
the full gradient and FieldEnergy unchanged. Reordering member construction and
top-level inputs leaves the recovered shared-value map and all canonical
observables unchanged.

The manufactured affine case also closes a tolerance-policy edge: because an
affine field has zero Cubic seminorm, NUM-009's characteristic field scale now
also includes the explicit gradient magnitude times the solve-coordinate
length. The scale remains gauge-invariant while preventing a physical nonzero
affine field from receiving a near-zero value tolerance.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
