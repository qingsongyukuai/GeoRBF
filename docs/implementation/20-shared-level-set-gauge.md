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

Lowering materializes one stable semantic latent variable per group. The
variable owns its `GroupId`, field unit, and every member `SourceId`. It remains
an explicit primal variable from Canonical Problem IR through the augmented KKT
and physical recovery; it is not reconstructed from the field afterward.

The current all-hard Cubic path derives a deterministic equality realization:

- every member lowers to `f(x_member) - h_group = 0` with its own source and
  group provenance;
- a level-set gauge lowers to `h_group = value`;
- a point gauge lowers to `f(x_gauge) = value` as a distinct convention
  relation.

No member is assigned a known value or selected as a semantic reference. Exact
duplicate relations share a canonical row while retaining every caller-owned
assessment. Member equations, true field-value observations, and the active
gauge enter the KKT through a deterministic spanning forest of the
group/point/absolute-reference incidence graph; cycle-closing equations stay in
Canonical IR as verification-only relations. This avoids redundant KKT rows
for coincident, partially overlapping, or absolutely observed groups without
merging `GroupId` identity.

If no true field-value observation exists, the first gauge in stable `SourceId`
order supplies the one solver constraint that selects the additive constant.
All additional gauges remain canonical verification relations. If a true
field-value observation already selects the absolute field, every explicit
gauge is verification-only. Consequently additional conventions can accept or
reject the chosen representative but cannot alter gradient, value differences,
geometry, or FieldEnergy.

Recover and Verify reads each semantic latent directly from the backend
candidate, independently evaluates every member and gauge equation, applies a
relation-specific physical tolerance, and rejects the candidate if any
verification-only convention is incompatible. There is no post-hoc averaging
or tolerance borrowed from another relation.

The recovered `SharedLevelValue` retains `GroupId`, field unit, and all member
sources. `SolvedModel::shared_level_value` exposes the same immutable value by
`GroupId`. `HardRelationAssessment::group_id` completes source/group provenance
for member and gauge relations without exposing a backend column, coefficient,
or matrix.

## Preflight and physical verification

A group with exactly one member and whose latent is not referenced produces
typed `UninformativeSharedLevelSet` evidence before backend execution. A
level-set gauge is a real latent reference, so a gauged singleton is not
misreported. Repeated locations in a multi-member group do not redefine the
one-member semantic rule. If every relation is invariant to the simultaneous
shift `f -> f + c`, `h_group -> h_group + c` and there is neither a true
field-value observation nor an explicit gauge, fit returns typed
`UnidentifiedAdditiveGauge` evidence with stably ordered source and group IDs,
again before backend execution.

Manufactured affine horizons pass through the default Cubic/full-`Pi1`/faer KKT
path. Recovery verifies member equalities, the gauge, full gradient, field value,
zero affine Cubic quotient energy, and complete provenance. Changing only the
gauge shifts queried and shared field values by the same constant while leaving
the full gradient and FieldEnergy unchanged. Reordering member construction and
top-level inputs leaves the recovered shared-value map and all canonical
observables unchanged. Compatible secondary gauges leave geometry and energy
unchanged; an incompatible secondary gauge is retained in the failed recovery
report with its own provenance, residual, and physical tolerance.

Capacity preflight counts every scalar source relation and semantic latent
before Canonical IR allocation. Its conservative dense plan assumes no
duplicate or verification-only row can be removed, so canonical and report
storage remain bounded even when the eventual independent KKT is smaller.

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
