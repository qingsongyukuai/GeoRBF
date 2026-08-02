use georbf::diagnostics::{
    AnalysisFailureEvidence, ProblemDiagnosis, RankDecision, RankEvidenceDomain,
    SolveCoordinateFailureReason,
};
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
    assert_eq!(failure.report().problem_size().input_observations(), 1);
    assert_eq!(failure.report().problem_size().scalar_hard_relations(), 1);
    assert_eq!(failure.report().field_energy(), None);
    assert_eq!(failure.report().total_objective(), None);
    let rank = failure
        .report()
        .rank_evidence()
        .expect("the field-mode diagnosis retains its rank proof");
    assert_eq!(rank.domain(), RankEvidenceDomain::CubicPolynomialPairing);
    assert_eq!(rank.decision(), RankDecision::RankDeficient);
    assert_eq!(rank.rank(), Some(1));
    assert!(failure.report().backend_rank().is_none());
}

#[test]
fn contradictory_hard_values_return_typed_source_evidence_without_solver_attempts() {
    let location = Point3::try_new(0.0, 0.0, 0.0).unwrap();
    let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field-unit"));
    builder
        .add(FieldValueObservation::try_new(SourceId::new("source-b"), location, 2.0).unwrap())
        .unwrap();
    builder
        .add(FieldValueObservation::try_new(SourceId::new("source-a"), location, 1.0).unwrap())
        .unwrap();

    let failure = builder
        .build()
        .unwrap()
        .fit()
        .expect_err("incompatible exact values at one location are infeasible");
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
    let conflict = failure
        .report()
        .direct_input_conflict()
        .expect("the diagnosis carries stable source and target evidence");
    assert_eq!(conflict.first_source().as_str(), "source-a");
    assert_eq!(conflict.second_source().as_str(), "source-b");
    assert_eq!(
        conflict.semantic_role().as_str(),
        "field-value-observation/value"
    );
    assert_eq!(conflict.first_target(), 1.0);
    assert_eq!(conflict.second_target(), 2.0);
}

#[test]
fn finite_but_unscalable_coordinates_retain_their_typed_failure_reason() {
    let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field-unit"));
    for (source, x, target) in [("left", -1.0e103, -1.0), ("right", 1.0e103, 1.0)] {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(source),
                    Point3::try_new(x, 0.0, 0.0).unwrap(),
                    target,
                )
                .unwrap(),
            )
            .unwrap();
    }

    let failure = builder
        .build()
        .unwrap()
        .fit()
        .expect_err("cubing the characteristic solve length must fail closed");
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::NumericalFailure);
    assert!(matches!(
        failure.report().analysis_failure(),
        Some(AnalysisFailureEvidence::InvalidSolveCoordinateTransform {
            reason: SolveCoordinateFailureReason::FieldRecoveryScaleNotInvertible,
            backend_invoked: false,
        })
    ));
}

#[test]
fn unsupported_exact_thread_requests_fail_before_snapshot_creation() {
    use std::num::NonZeroUsize;

    use georbf::problem::{FitConfiguration, ThreadBudget};

    let location = Point3::try_new(0.0, 0.0, 0.0).unwrap();
    let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field-unit"));
    builder
        .add(FieldValueObservation::try_new(SourceId::new("value"), location, 1.0).unwrap())
        .unwrap();
    builder.set_fit_configuration(
        FitConfiguration::default()
            .with_thread_budget(ThreadBudget::Exact(NonZeroUsize::new(2).unwrap())),
    );

    let failure = builder
        .build()
        .expect_err("the admitted Cubic Equality adapter is sequential");
    assert_eq!(
        failure.errors(),
        &[BuildError::UnsupportedThreadBudget { requested: 2 }]
    );
}
