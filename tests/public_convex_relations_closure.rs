use georbf::diagnostics::SolveAttemptTermination;
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{
    AxialNormalObservation, CovarianceGroupBuilder, CovarianceMatrix, DirectedNormalObservation,
    FieldValueObservation, MinimumNormalSlope,
};
use georbf::relation::{
    DirectionalDerivativeInterval, FieldLevelOrder, FieldSeparationInterval, FieldValueBound,
    MinimumFieldOffset, MinimumFieldSeparation, PointToLevelSetRelation, PointToLevelSetSide,
    PolarityResolution, PolaritySelection, SharedLevelSetBuilder, StratigraphicFieldDirection,
    YoungerThan,
};
use georbf::{GroupId, Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).unwrap()
}

#[test]
fn every_public_v02_convex_relation_shares_one_canonical_acceptance_report() {
    // Manufactured truth: f(x, y, z) = x + 2y + 3z.
    let mut problem = ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field"),
    );
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    for (source, support, value) in [
        ("value-origin", [0.0, 0.0, 0.0], 0.0),
        ("value-east", [1.0, 0.0, 0.0], 1.0),
        ("value-north", [0.0, 1.0, 0.0], 2.0),
        ("value-up", [0.0, 0.0, 1.0], 3.0),
        ("value-diagonal", [1.0, 1.0, 0.0], 3.0),
    ] {
        problem
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(source),
                    point(support[0], support[1], support[2]),
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }

    let mut reference = SharedLevelSetBuilder::new(GroupId::new("reference"));
    reference
        .add_member(SourceId::new("reference/origin"), point(0.0, 0.0, 0.0))
        .unwrap();
    reference
        .add_member(SourceId::new("reference/second"), point(2.0, -1.0, 0.0))
        .unwrap();
    problem.add(reference.build().unwrap()).unwrap();

    let mut target = SharedLevelSetBuilder::new(GroupId::new("target"));
    target
        .add_member(SourceId::new("target/diagonal"), point(1.0, 1.0, 0.0))
        .unwrap();
    target
        .add_member(SourceId::new("target/up"), point(0.0, 0.0, 1.0))
        .unwrap();
    problem.add(target.build().unwrap()).unwrap();

    problem
        .add(
            FieldValueBound::try_interval(
                SourceId::new("field-bound"),
                point(0.0, 2.0, 0.0),
                3.5,
                4.5,
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(
            DirectionalDerivativeInterval::try_interval(
                SourceId::new("derivative-interval"),
                point(0.25, -0.5, 0.75),
                vector(1.0, 0.0, 0.0),
                0.9,
                1.1,
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(YoungerThan::hard(
            SourceId::new("younger-than"),
            GroupId::new("target"),
            GroupId::new("reference"),
            MinimumFieldSeparation::try_new(2.0).unwrap(),
        ))
        .unwrap();
    problem
        .add(FieldLevelOrder::hard(
            SourceId::new("field-order"),
            GroupId::new("reference"),
            GroupId::new("target"),
        ))
        .unwrap();
    problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("field-separation"),
                GroupId::new("reference"),
                GroupId::new("target"),
                2.5,
                3.5,
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(PointToLevelSetRelation::hard(
            SourceId::new("point-side"),
            point(1.0, 0.0, 1.0),
            GroupId::new("target"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(0.5).unwrap(),
        ))
        .unwrap();

    let normal = vector(1.0, 2.0, 3.0);
    problem
        .add(
            DirectedNormalObservation::try_new(
                SourceId::new("directed-normal"),
                point(0.2, -0.3, 0.4),
                normal,
                MinimumNormalSlope::try_new(3.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(PolarityResolution::new(
            SourceId::new("axial-resolution"),
            SourceId::new("axial-normal"),
            PolaritySelection::AlongInputAxis,
        ))
        .unwrap();
    problem
        .add(
            AxialNormalObservation::try_new(
                SourceId::new("axial-normal"),
                point(-0.4, 0.2, 0.1),
                normal,
                MinimumNormalSlope::try_new(3.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();

    let mut covariance = CovarianceGroupBuilder::new(GroupId::new("field-covariance"));
    covariance
        .add_field_value_member(SourceId::new("soft-east"), point(1.0, 0.0, 0.0), 1.0)
        .unwrap();
    covariance
        .add_field_value_member(SourceId::new("soft-north"), point(0.0, 1.0, 0.0), 2.0)
        .unwrap();
    problem
        .add(
            covariance
                .build(CovarianceMatrix::try_new([[1.0, 0.25], [0.25, 2.0]]).unwrap())
                .unwrap(),
        )
        .unwrap();
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let success = problem.build().unwrap().fit().unwrap();
    let report = success.report();
    assert!(report.infeasibility_certificate().is_none());
    assert!(report.recession_ray().is_none());
    assert_eq!(report.field_value_bounds().len(), 2);
    assert_eq!(report.directional_derivative_intervals().len(), 2);
    assert_eq!(report.shared_level_set_relations().len(), 2);
    assert_eq!(report.field_separation_intervals().len(), 2);
    assert_eq!(report.point_to_level_set_relations().len(), 1);
    assert_eq!(report.directed_normals().len(), 2);
    assert_eq!(report.covariance_groups().len(), 1);
    assert_eq!(report.shared_level_values().len(), 2);
    assert!(report.field_energy().unwrap().abs() <= 1.0e-8);
    assert!(report.total_objective().unwrap().abs() <= 1.0e-8);

    let acceptance = report.canonical_acceptance().unwrap();
    assert!(acceptance.accepted());
    assert!(acceptance.recovery_finite());
    assert!(acceptance.provenance_verified());
    assert!(acceptance.objective_verified());
    assert!(acceptance.hard_affine_inequality_violation_max().unwrap() <= 1.0e-8);
    assert!(acceptance.scaling_round_trip_error().unwrap() <= 1.0e-11);
    assert!(acceptance.reduction_round_trip_error().unwrap() <= 1.0e-11);
    assert!(
        report
            .cubic_analysis()
            .unwrap()
            .reduced_condition_estimate()
            .is_some()
    );

    let attempt = report.attempts().last().unwrap();
    fn assert_eq_contract<T: Eq>() {}
    assert_eq_contract::<georbf::diagnostics::ScalingSummary>();
    assert!(attempt.scaling_round_trip_error().unwrap() <= 1.0e-11);
    assert!(matches!(
        attempt.termination(),
        SolveAttemptTermination::CandidateProduced
            | SolveAttemptTermination::ReducedAccuracyCandidateProduced
    ));
    let residual = attempt.convex_residual().unwrap();
    assert!(residual.primal() <= 1.0e-8);
    assert!(residual.dual() <= 1.0e-8);
    assert!(residual.stationarity() <= 1.0e-8);
    assert!(residual.complementarity() <= 1.0e-8);
    assert!(residual.relative_gap() <= 1.0e-8);
    let physical_residual = acceptance.physical_convex_residual().unwrap();
    assert!(physical_residual.primal() <= 1.0e-8);
    assert!(physical_residual.dual() <= 1.0e-8);
    assert!(physical_residual.stationarity() <= 1.0e-8);
    assert!(physical_residual.complementarity() <= 1.0e-8);
    assert!(physical_residual.relative_gap() <= 1.0e-8);

    assert_eq!(
        report
            .covariance_groups()
            .iter()
            .map(|group| group.objective_contribution())
            .sum::<f64>(),
        report.covariance_groups()[0].objective_contribution()
    );
}
