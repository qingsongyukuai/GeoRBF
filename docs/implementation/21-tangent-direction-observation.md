# Issue 21: tangent direction observations

Issue: [#21](https://github.com/qingsongyukuai/GeoRBF/issues/21)

Evidence seams: T03, T04, T11

Requirements: DOM-010, DOM-017; IR-002, IR-003, IR-005, IR-008

## Public tangent semantics

`TangentDirectionObservation::try_new` accepts a caller-owned `SourceId`, a
checked physical-frame point, and a checked finite vector. It rejects the zero
vector and uses scaled normalization so finite subnormal and near-maximum
components cannot underflow or overflow during normalization. The first
nonzero component selects one deterministic sign representative, so opposite
and nonzero scaled inputs have identical axial semantics. The stored direction
is unit length in the declared physical input coordinates; its original length
does not become a weight, confidence, gradient magnitude, or polarity.

The observation means only `t^T grad(f) = 0`. A zero gradient satisfies it. It
does not assert a regular level set, a full gradient, a directed or axial
normal, a minimum slope, or an angular tolerance. Those distinctions are also
part of the public rustdoc; normal observations remain outside this milestone.

## Canonical Equality path

Each tangent lowers to one first-order scalar functional with
field-value-per-length dimension, target zero, and semantic role
`tangent-direction-observation/directional-derivative`. Stable RelationId and
ResidualId values derive from the caller's SourceId and that role, independently
of input ordering and backend row layout. The normalized direction is the
functional's gradient coefficient, so `t` and `-t` intern as the same
generalized functional and exact duplicate hard relations can share one
canonical equality without losing per-source assessment.

The functional participates in the existing Cubic generalized-functional
pairing, complete-Pi1 representation analysis, faer augmented KKT, and physical
recovery. Recovery evaluates the accepted field at the original support and
recomputes the direction/gradient contraction; it does not accept an internal
row, coefficient, or backend residual as proof. Provenance verification covers
source, relation, residual, and every derived equality artifact before a model
can be returned.

Zero-target tangent relations also close a tolerance-scale edge. The
gauge-invariant characteristic field scale now includes variation from explicit
field-value semantics, alongside FieldEnergy and explicit derivative targets.
This gives derivative relations a physical field-per-length acceptance scale
without using an absolute gauge value, latent value, backend row norm, scaled
right-hand side, or candidate.

## Validation evidence

- T03 checks finite/nonzero construction, physical normalization, canonical
  sign, canonical positive zero, and atomic duplicate-SourceId rejection.
- T04 fits a manufactured affine field through the public Cubic Equality path,
  independently recomputes `t^T grad(f)` from `SolvedModel::evaluate`, and
  checks derivative units, target, stable role path, SourceId, hard residual,
  FieldEnergy, and recovered provenance. Three independent tangent axes at a
  stationary support prove that zero gradient is admitted without a slope
  assumption.
- T04 also proves that opposite/scaled directions and top-level input
  permutation produce identical samples, relation assessments, sizes, energy,
  objective, and provenance acceptance.
- T11 applies both a proper rotation and a reflection together with translation
  and uniform scale. Recovered value/gradient obey the frame chain rule, while
  tangent recovered values, characteristic scales, and tolerances covary in
  field-value-per-length units.

## Replay

Under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
python3 scripts/audit.py
cargo package --locked
```
