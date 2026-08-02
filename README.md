# GeoRBF

GeoRBF is a Rust library for fitting implicit geological scalar fields from
geological observations.

Version 0.1.0 exposes the first complete public tracer: callers declare an
input coordinate frame and units, add hard absolute field-value and complete
gradient observations with stable `SourceId`s, build an immutable snapshot,
fit through the Cubic Equality path, and query field value and gradient
together from an immutable model.

```rust,no_run
use georbf::geometry::{
    FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel,
};
use georbf::observation::{FieldValueObservation, GradientObservation};
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

let snapshot = problem.build()?;
let fit = snapshot.fit()?;
let sample = fit.model().evaluate(Point3::try_new(0.2, -0.3, 0.4)?)?;
assert_eq!(sample.gradient().components().len(), 3);
# Ok(())
# }
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
tracer in [#19](docs/implementation/19-public-absolute-field-tracer.md).
