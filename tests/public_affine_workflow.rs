use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::KernelKind;
use georbf::observation::{FieldValueObservation, GradientObservation};
use georbf::{Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the manufactured point is finite")
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).expect("the manufactured vector is finite")
}

fn affine_value(point: Point3) -> f64 {
    let [x, y, z] = point.components();
    2.0 + 0.5 * x - 1.25 * y + 0.75 * z
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-9 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn user_can_fit_and_sample_an_absolute_affine_field() {
    let frame = InputCoordinateFrame::try_new(
        ["east", "north", "elevation"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )
    .expect("the frame has three distinct axis labels");
    let field_unit = FieldUnitLabel::new("stratigraphic-coordinate");
    let mut builder = ProblemBuilder::new(frame.clone(), field_unit.clone());

    for (index, location) in [
        point(-1.0, -1.0, -1.0),
        point(1.0, -1.0, -1.0),
        point(-1.0, 1.0, -1.0),
        point(-1.0, -1.0, 1.0),
        point(1.0, 1.0, 0.5),
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("value-{index}")),
                    location,
                    affine_value(location),
                )
                .expect("the manufactured field value is finite"),
            )
            .expect("the SourceId is unique");
    }

    for (index, location) in [
        point(0.25, -0.5, 0.75),
        point(-0.75, 0.25, 0.5),
        point(0.5, 0.75, -0.25),
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(GradientObservation::new(
                SourceId::new(format!("gradient-{index}")),
                location,
                vector(0.5, -1.25, 0.75),
            ))
            .expect("the SourceId is unique");
    }

    let snapshot = builder.build().expect("the hard problem is valid");
    assert_eq!(snapshot.input_coordinate_frame(), &frame);
    assert_eq!(snapshot.field_unit(), &field_unit);
    assert_eq!(snapshot.resolved_kernel().kind(), KernelKind::Cubic);
    assert_eq!(snapshot.field_energy_normalization().factor(), 1.0);

    let success = snapshot.fit().expect("the affine field is identifiable");
    let query = point(0.2, -0.3, 0.4);
    let sample = success
        .model()
        .evaluate(query)
        .expect("a finite query has a finite result");
    assert_close(sample.value(), affine_value(query));
    for (actual, expected) in sample
        .gradient()
        .components()
        .into_iter()
        .zip([0.5, -1.25, 0.75])
    {
        assert_close(actual, expected);
    }

    let report = success.report();
    assert_eq!(report.resolved_kernel().kind(), KernelKind::Cubic);
    assert_eq!(report.field_energy_normalization().factor(), 1.0);
    assert_close(report.field_energy().expect("fit succeeded"), 0.0);
    assert_eq!(report.hard_relations().len(), 14);
    assert!(report.hard_relations().iter().all(|relation| {
        relation.residual().abs() <= relation.tolerance()
            && !relation.source_id().as_str().is_empty()
            && !relation.semantic_role().is_empty()
    }));
    assert_eq!(
        report
            .hard_relations()
            .iter()
            .filter(|relation| relation.source_id().as_str() == "gradient-0")
            .map(|relation| relation.semantic_role())
            .collect::<Vec<_>>(),
        [
            "gradient-observation/component/0",
            "gradient-observation/component/1",
            "gradient-observation/component/2",
        ]
    );

    fn assert_send_sync_clone<T: Send + Sync + Clone>() {}
    assert_send_sync_clone::<georbf::SolvedModel>();
}
