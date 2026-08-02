use georbf::geometry::{
    FieldUnitLabel, GlobalAnisotropyMetric, Handedness, InputCoordinateFrame, LengthUnitLabel,
};
use georbf::observation::{FieldValueObservation, GradientObservation};
use georbf::{Point3, ProblemBuilder, SourceId, Vector3};

const SCALE: f64 = 2.5;
const TRANSLATION: [f64; 3] = [10.0, -3.0, 4.0];

fn point(components: [f64; 3]) -> Point3 {
    Point3::try_new(components[0], components[1], components[2]).unwrap()
}

fn vector(components: [f64; 3]) -> Vector3 {
    Vector3::try_new(components[0], components[1], components[2]).unwrap()
}

fn transform_point([x, y, z]: [f64; 3]) -> [f64; 3] {
    [
        SCALE * y + TRANSLATION[0],
        SCALE * x + TRANSLATION[1],
        SCALE * z + TRANSLATION[2],
    ]
}

fn transform_gradient([gx, gy, gz]: [f64; 3]) -> [f64; 3] {
    [gy / SCALE, gx / SCALE, gz / SCALE]
}

fn truth([x, y, z]: [f64; 3]) -> (f64, [f64; 3]) {
    (
        2.0 + 0.5 * x - 1.25 * y + 0.75 * z + 0.2 * x * x - 0.15 * x * y + 0.1 * y * z,
        [
            0.5 + 0.4 * x - 0.15 * y,
            -1.25 - 0.15 * x + 0.1 * z,
            0.75 + 0.1 * y,
        ],
    )
}

fn build_problem(transformed: bool) -> georbf::ProblemSnapshot {
    let frame = if transformed {
        InputCoordinateFrame::try_new(
            ["north-prime", "east-prime", "elevation-prime"],
            Handedness::Left,
            LengthUnitLabel::new("scaled-length"),
        )
        .unwrap()
    } else {
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("length"),
        )
        .unwrap()
    };
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("field"));
    let metric = if transformed {
        [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.5]]
    } else {
        [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]]
    };
    builder
        .set_global_anisotropy_metric(GlobalAnisotropyMetric::try_from_matrix(metric).unwrap())
        .unwrap();

    for (index, support) in [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, 1.0, 0.5],
    ]
    .into_iter()
    .enumerate()
    {
        let (value, _) = truth(support);
        let location = if transformed {
            transform_point(support)
        } else {
            support
        };
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("value-{index}")),
                    point(location),
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    for (index, support) in [[0.25, -0.5, 0.75], [-0.75, 0.25, 0.5], [0.5, 0.75, -0.25]]
        .into_iter()
        .enumerate()
    {
        let (_, gradient) = truth(support);
        let (location, gradient) = if transformed {
            (transform_point(support), transform_gradient(gradient))
        } else {
            (support, gradient)
        };
        builder
            .add(GradientObservation::new(
                SourceId::new(format!("gradient-{index}")),
                point(location),
                vector(gradient),
            ))
            .unwrap();
    }
    builder.build().unwrap()
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 2.0e-8 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn field_observables_are_covariant_under_a_legal_frame_transform() {
    let original = build_problem(false).fit().unwrap();
    let transformed = build_problem(true).fit().unwrap();
    let query = [0.2, -0.3, 0.4];
    let original_sample = original.model().evaluate(point(query)).unwrap();
    let transformed_sample = transformed
        .model()
        .evaluate(point(transform_point(query)))
        .unwrap();

    assert_close(transformed_sample.value(), original_sample.value());
    for (actual, expected) in transformed_sample
        .gradient()
        .components()
        .into_iter()
        .zip(transform_gradient(original_sample.gradient().components()))
    {
        assert_close(actual, expected);
    }
    assert!(
        original
            .report()
            .hard_relations()
            .iter()
            .all(|relation| relation.residual().abs() <= relation.tolerance())
    );
    assert!(
        transformed
            .report()
            .hard_relations()
            .iter()
            .all(|relation| relation.residual().abs() <= relation.tolerance())
    );
}
