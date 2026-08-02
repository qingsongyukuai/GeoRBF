use georbf::diagnostics::ProblemDiagnosis;
use georbf::geometry::{
    FieldUnitLabel, GeometryError, GlobalAnisotropyMetric, GlobalAnisotropyMetricError, Handedness,
    InputCoordinateFrame, LengthUnitLabel,
};
use georbf::observation::{FieldValueObservation, ObservationError};
use georbf::problem::{AddError, BuildError};
use georbf::{Point3, ProblemBuilder, SourceId, Vector3};

fn frame() -> InputCoordinateFrame {
    InputCoordinateFrame::try_new(
        ["axis-0", "axis-1", "axis-2"],
        Handedness::Right,
        LengthUnitLabel::new("metre"),
    )
    .expect("the fixture frame is valid")
}

#[test]
fn checked_leaf_constructors_reject_invalid_values() {
    assert!(matches!(
        Point3::try_new(0.0, f64::NAN, 0.0),
        Err(GeometryError::NonFinitePoint { axis: 1 })
    ));
    assert!(matches!(
        Vector3::try_new(0.0, 0.0, f64::INFINITY),
        Err(GeometryError::NonFiniteVector { axis: 2 })
    ));
    assert!(matches!(
        InputCoordinateFrame::try_new(
            ["east", "east", "up"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        ),
        Err(GeometryError::DuplicateAxisLabel {
            first: 0,
            second: 1
        })
    ));
    assert!(matches!(
        FieldValueObservation::try_new(
            SourceId::new("bad-value"),
            Point3::try_new(0.0, 0.0, 0.0).unwrap(),
            f64::NEG_INFINITY,
        ),
        Err(ObservationError::NonFiniteFieldValue)
    ));

    assert!(matches!(
        GlobalAnisotropyMetric::try_from_matrix([
            [1.0, 1.0e-15, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]),
        Err(GeometryError::InvalidGlobalAnisotropyMetric(
            GlobalAnisotropyMetricError::NotSymmetric { .. }
        ))
    ));
    assert!(matches!(
        GlobalAnisotropyMetric::try_from_matrix([
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
        ]),
        Err(GeometryError::InvalidGlobalAnisotropyMetric(
            GlobalAnisotropyMetricError::NotPositiveDefinite
        ))
    ));
    assert!(matches!(
        GlobalAnisotropyMetric::try_from_matrix([
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]),
        Err(GeometryError::InvalidGlobalAnisotropyMetric(
            GlobalAnisotropyMetricError::DeterminantNotOne { .. }
        ))
    ));
}

#[test]
fn builder_rejections_are_atomic_and_failed_builds_are_repairable() {
    let location = Point3::try_new(0.0, 0.0, 0.0).unwrap();
    let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field-unit"));
    builder
        .add(FieldValueObservation::try_new(SourceId::new("same"), location, 1.0).unwrap())
        .unwrap();
    assert!(matches!(
        builder.add(FieldValueObservation::try_new(SourceId::new("same"), location, 2.0).unwrap()),
        Err(AddError::DuplicateSourceId { .. })
    ));
    builder
        .add(FieldValueObservation::try_new(SourceId::new("different"), location, 3.0).unwrap())
        .expect("a rejected duplicate did not mutate the builder");
    assert_eq!(builder.build().unwrap().observation_count(), 2);

    let empty = ProblemBuilder::new(frame(), FieldUnitLabel::new("field-unit"));
    let failure = empty.build().expect_err("an empty snapshot is invalid");
    assert_eq!(failure.errors(), &[BuildError::NoObservations]);
    let mut repaired = failure.into_builder();
    repaired
        .add(FieldValueObservation::try_new(SourceId::new("repair"), location, 0.0).unwrap())
        .unwrap();
    assert_eq!(repaired.build().unwrap().observation_count(), 1);
}

#[test]
fn an_unidentified_field_returns_a_typed_failure_without_model_data() {
    let location = Point3::try_new(0.0, 0.0, 0.0).unwrap();
    let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field-unit"));
    builder
        .add(FieldValueObservation::try_new(SourceId::new("only-value"), location, 1.0).unwrap())
        .unwrap();
    let snapshot = builder.build().unwrap();

    let failure = snapshot.fit().expect_err("one value cannot identify Π₁");
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::UnidentifiedFieldMode);
    assert_eq!(failure.report().problem_size().observations(), 1);
    assert_eq!(failure.report().problem_size().scalar_hard_relations(), 1);
    assert_eq!(failure.report().field_energy(), None);
    assert_eq!(failure.report().total_objective(), None);
}
