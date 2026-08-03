# GeoRBF

GeoRBF is a Rust library for fitting implicit geological scalar fields from
geological observations.

Version 0.1.0 exposes complete public Cubic Equality tracers. The current v0.2
development surface adds the first complete soft tracer: a Field
Value Observation may carry either a checked positive `QuadraticPenalty` or a
checked positive statistical `StandardDeviation`, while the problem supplies
an explicit checked `FieldEnergyNormalization`. Callers otherwise declare
an input coordinate frame and units; add hard absolute field-value and complete
gradient observations, unoriented tangent directions, or atomically built
shared level sets and geological horizons; choose an explicit additive field
gauge where the relations are translation-invariant; build an immutable
snapshot; and query field value and complete gradient through single points or
ordered, atomic logical batches from an immutable model. Recovered shared values
remain addressable by `GroupId`. A tangent direction constrains only
`t^T grad(f) = 0`: it allows zero gradient and carries neither the magnitude of
a complete gradient observation nor the polarity and nonzero-slope semantics of
a normal direction.

The complete release workflow is available as a runnable
[Equality Spine example](examples/equality_spine.rs). It combines a planar
horizon, explicit additive gauge, absolute field value, complete gradient,
tangent direction, typed fit evidence, shared-value recovery, and ordered
single/batch queries. The supported scope and compatibility boundary are
recorded in the [v0.1.0 release notes](RELEASE_NOTES.md).

```rust,no_run
use georbf::geometry::{
    FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel,
};
use georbf::observation::{
    FieldValueObservation, GradientObservation, TangentDirectionObservation,
};
use georbf::{Point3, ProblemBuilder, SourceId, Vector3};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let frame = InputCoordinateFrame::try_new(
    ["east", "north", "elevation"],
    Handedness::Right,
    LengthUnitLabel::new("m"),
)?;
let mut problem = ProblemBuilder::new(frame, FieldUnitLabel::new("field-unit"));

let supports = [
    ([-1.0, -1.0, -1.0], 2.0),
    ([1.0, -1.0, -1.0], 3.0),
    ([-1.0, 1.0, -1.0], -0.5),
    ([-1.0, -1.0, 1.0], 3.5),
    ([1.0, 1.0, 0.5], 1.625),
];
for (index, (coordinates, value)) in supports.into_iter().enumerate() {
    let [x, y, z] = coordinates;
    problem.add(FieldValueObservation::try_new(
        SourceId::new(format!("value-{index}")),
        Point3::try_new(x, y, z)?,
        value,
    )?)?;
}
problem.add(GradientObservation::new(
    SourceId::new("gradient"),
    Point3::try_new(0.25, -0.5, 0.75)?,
    Vector3::try_new(0.5, -1.25, 0.75)?,
))?;
problem.add(TangentDirectionObservation::try_new(
    SourceId::new("tangent"),
    Point3::try_new(0.25, -0.5, 0.75)?,
    Vector3::try_new(2.5, 1.0, 0.0)?,
)?)?;

let snapshot = problem.build()?;
let fit = snapshot.fit()?;
let sample = fit.model().evaluate(Point3::try_new(0.2, -0.3, 0.4)?)?;
assert_eq!(sample.gradient().components().len(), 3);
let samples = fit.model().evaluate_batch(&[
    Point3::try_new(0.2, -0.3, 0.4)?,
    Point3::try_new(-0.5, 0.25, 0.75)?,
])?;
assert_eq!(samples.len(), 2);
# Ok(())
# }
```

The input boundary is sealed: downstream crates cannot add custom observation,
solver, matrix, or backend inputs to `ProblemBuilder`.

## Soft Field Value objective

Soft Field Value configuration uses named constructors rather than a generic
enforcement or loss enum. A soft problem must explicitly set the physical
FieldEnergy scale before it can become an immutable snapshot:

```rust,no_run
use georbf::geometry::{
    FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel,
};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{FieldValueObservation, QuadraticPenalty};
use georbf::{Point3, ProblemBuilder, SourceId};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let frame = InputCoordinateFrame::try_new(
    ["east", "north", "elevation"],
    Handedness::Right,
    LengthUnitLabel::new("m"),
)?;
let mut problem = ProblemBuilder::new(frame, FieldUnitLabel::new("field-unit"));
problem.add(FieldValueObservation::try_with_quadratic_penalty(
    SourceId::new("soft-value"),
    Point3::try_new(0.0, 0.0, 0.0)?,
    1.25,
    QuadraticPenalty::try_new(2.0)?,
)?)?;
problem.set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0)?)?;
let snapshot = problem.build()?;
# let _ = snapshot;
# Ok(())
# }
```

`FitReport::soft_field_values` returns original-unit targets, recovered values
and residuals, typed penalty/statistical configuration, and each independent
loss contribution. `field_energy` and `total_objective` are independently
recomputed during recovery; a failed objective or provenance round trip returns
a structured fit failure and never a model.

```compile_fail
struct CustomInput;

impl georbf::problem::ProblemInput for CustomInput {
    fn add_to(
        self,
        _builder: &mut georbf::ProblemBuilder,
    ) -> Result<(), georbf::problem::AddError> {
        Ok(())
    }
}
```

## Verification

The crate is pinned to Rust 1.85.0 and uses the checked-in lockfile:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
python3 scripts/audit.py
cargo package --locked
```

Implementation evidence is recorded for
[#16](docs/implementation/16-production-equality-spine.md),
[#17](docs/implementation/17-cubic-equality-core.md),
[#18](docs/implementation/18-cubic-equality-numerical-policy.md), and the public
tracers in [#19](docs/implementation/19-public-absolute-field-tracer.md) and
[#20](docs/implementation/20-shared-level-set-gauge.md), plus
[#21](docs/implementation/21-tangent-direction-observation.md) and the immutable
query model in [#24](docs/implementation/24-immutable-model-batch-query.md).
The cumulative release evidence and traceability audit are recorded in
[#25](docs/implementation/25-equality-spine-release.md).
The first v0.2 soft-objective tracer and its requirement mapping are recorded in
[#27](docs/implementation/27-soft-field-value-objective.md).
