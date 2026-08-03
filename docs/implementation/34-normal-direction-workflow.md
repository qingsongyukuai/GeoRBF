# Issue 34: Directed/Axial Normal and Polarity Resolution

Issue: [#34](https://github.com/qingsongyukuai/GeoRBF/issues/34)

Primary evidence seams: T01, T03–T05, T07–T11, T14, T16

Requirements: PAPI-008–PAPI-010, PAPI-017; DOM-006–DOM-007, DOM-009,
DOM-017–DOM-018, DOM-021–DOM-022; IR-003, IR-005, IR-007–IR-010,
IR-014; NUM-009–NUM-010; DIA-002, DIA-006–DIA-009; VAL-004,
VAL-008, VAL-015

## Checked public semantics

`DirectedNormalObservation` and `AxialNormalObservation` accept only finite,
nonzero vectors and normalize them with the physical Euclidean frame's scaled
algorithm. Directed inputs preserve polarity. Axial inputs retain the normalized
caller orientation for later resolution and a sign-canonical identity shared by
positive and negative scaling. `MinimumNormalSlope` accepts only finite, strictly
positive field-value-per-length quantities.

The direction and minimum-slope channels use separate family-specific
enforcement types. Direction permits hard projection, a Euclidean quadratic
penalty, or isotropic statistical standard deviation. Minimum slope permits a
hard lower bound or an independently penalized nonnegative violation with a
quadratic or linear loss. There is no generic Enforcement/Loss input, default
slope, direction-component L1 loss, cone, angle, or SOCP surface.

`PolarityResolution` owns a separate `SourceId`, may forward-reference an Axial
Normal `SourceId`, and explicitly chooses `AlongInputAxis` or
`AgainstInputAxis`. Repeated, contradictory, dangling, and wrong-kind decisions
are retained through insertion and rejected together as deterministic cross-record
build evidence. The immutable snapshot retains the original Axial input and the
resolution record rather than rewriting either one.

## Preflight and canonical lowering

Every unresolved Axial Normal produces stable `UnresolvedAxialNormalEvidence`.
Unresolved semantics has the highest diagnosis priority, stops before backend
execution, returns no model, and coexists with any lower-priority conflict
evidence found in the same preflight.

A Directed or resolved Axial Normal lowers into two independent physical
channels. Direction uses the complete rotation-invariant projection
`(I - n n^T) grad(f)`. The canonical residual retains every nonzero coordinate
of that physical vector. The representation plan uses the complete physical
gradient-component span at each Normal support and reconstructs every projection
row from it; the derived hard solver form selects an independent row set while
keeping the omitted dependent component as a verification-only canonical
relation. This prevents multiple Normal directions at one support from creating
a dependent representer span without choosing a caller-visible tangent basis.
Hard recovery retains those physical coordinates as one ordered block and judges
its Euclidean residual norm against one channel tolerance, so finite-tolerance
acceptance is rotation invariant.
Minimum slope lowers independently as `n^T grad(f) >= s_min` and therefore selects
the form-driven Convex QP route. Hard slope facts also participate in exact
directional conflict preflight.

Soft direction uses one residual block and one isotropic precision, so its loss
is invariant under rotations and reflections. Soft slope owns a separate
nonnegative violation variable, stable residual identity, and objective term.
The determinant-one Global Anisotropy Metric enters only Cubic kernel distance;
normalization, projection, slope, and diagnostics remain physical Euclidean
quantities.

## Recover, verify, and evidence

Recover and Verify checks the projection equalities or soft block, the affine
slope relation and violation equation, provenance, transformation round trips,
FieldEnergy, and the total objective in physical units before exposing a model.
`DirectedNormalAssessment` reports source and channel roles, resolved direction
and complete gradient, projection vector and Euclidean norm, directed slope,
minimum slope, slack, violation, tolerance, active state, configured losses,
separate objective contributions, and optional Axial input/resolution provenance.

`tests/public_normal_directions.rs` covers checked extreme-magnitude
normalization, Directed polarity, Axial identity and retained input orientation,
independent enforcement types, forward resolution, aggregated repeated/conflicting
decisions, dangling references, unresolved preflight, diagnosis priority, hard
and soft manufactured affine QP cases, literal split objectives, zero-gradient
failure, resolved-Axial/Directed query equivalence, and rotation, reflection,
uniform scaling, and non-identity anisotropy. The cumulative #28, #30, public,
canonical, corruption, capacity, audit, and package suites remain part of this
capability's release boundary.
