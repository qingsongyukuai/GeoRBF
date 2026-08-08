use georbf::geometry::{
    FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel, Point3, Vector3,
};
use georbf::observation::{
    DirectedNormalObservation, FieldValueObservation, GradientObservation, MinimumNormalSlope,
};
use georbf::relation::FieldValueBound;
use georbf::{ProblemBuilder, SourceId};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn affine_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    1.0 + x + 2.0 * y - 0.5 * z
}

fn fit_with_same_gradient_at_close_supports(use_qp: bool) -> georbf::fit::FitSuccess {
    let mut builder = ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field"),
    );
    for (index, location) in [
        point(0.0, 0.0, 0.0),
        point(1.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0),
        point(0.0, 0.0, 1.0),
        point(1.0, 1.0, 1.0),
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
                .unwrap(),
            )
            .unwrap();
    }
    for (source, location) in [
        ("close-gradient-a", point(0.25, 0.5, 0.75)),
        ("close-gradient-b", point(0.255, 0.5, 0.75)),
    ] {
        builder
            .add(GradientObservation::new(
                SourceId::new(source),
                location,
                Vector3::try_new(1.0, 2.0, -0.5).unwrap(),
            ))
            .unwrap();
    }
    if use_qp {
        builder
            .add(
                FieldValueBound::try_lower(
                    SourceId::new("loose-bound"),
                    point(0.0, 0.0, 0.0),
                    -100.0,
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

fn fit_with_consistent_hard_duplicate(use_qp: bool) -> georbf::fit::FitSuccess {
    let mut builder = ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field"),
    );
    let locations = [
        point(0.0, 0.0, 0.0),
        point(1.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0),
        point(0.0, 0.0, 1.0),
        point(1.0, 1.0, 1.0),
    ];
    for (index, location) in locations.into_iter().enumerate() {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("value-{index}")),
                    location,
                    affine_value(location),
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder
        .add(
            FieldValueObservation::try_new(
                SourceId::new("value-duplicate"),
                locations[0],
                affine_value(locations[0]),
            )
            .unwrap(),
        )
        .unwrap();
    if use_qp {
        builder
            .add(
                FieldValueBound::try_lower(SourceId::new("loose-bound"), locations[0], -100.0)
                    .unwrap(),
            )
            .unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

#[test]
fn consistent_hard_duplicates_keep_both_sources_after_solver_row_compression() {
    for fit in [
        fit_with_consistent_hard_duplicate(false),
        fit_with_consistent_hard_duplicate(true),
    ] {
        let report = fit.report();
        let ledger = report.all_source_recovery().unwrap();
        assert!(ledger.verified());
        assert_eq!(ledger.recovered_sources(), ledger.participating_sources());
        assert_eq!(
            ledger.recovery_edge_count(),
            report.problem_size().scalar_hard_relations()
        );
        assert!(ledger.recovery_edge_count() > ledger.canonical_hard_relation_count());
        for source in ["value-0", "value-duplicate"] {
            let recovered = report
                .hard_relations()
                .iter()
                .find(|relation| relation.source_id().as_str() == source)
                .unwrap();
            assert!(recovered.residual().abs() <= recovered.tolerance());
        }
    }
}

#[test]
fn close_supports_with_the_same_gradient_remain_distinct_in_kkt_and_qp() {
    let equality = fit_with_same_gradient_at_close_supports(false);
    let qp = fit_with_same_gradient_at_close_supports(true);

    assert_eq!(
        equality.report().problem_size().equality_constraints(),
        Some(15)
    );
    assert_eq!(qp.report().problem_size().equality_constraints(), Some(11));
    for success in [&equality, &qp] {
        let ledger = success
            .report()
            .all_source_recovery()
            .expect("an accepted model retains an all-source recovery ledger");
        assert!(ledger.verified());
        assert_eq!(ledger.representer_count(), 11);
        assert_eq!(
            ledger.canonical_hard_relation_count(),
            success.report().problem_size().scalar_hard_relations()
        );
        assert_eq!(ledger.canonical_soft_relation_count(), 0);
        assert_eq!(
            ledger.participating_sources().len(),
            if std::ptr::eq(success, &equality) {
                7
            } else {
                8
            }
        );
        assert_eq!(ledger.recovered_sources(), ledger.participating_sources());
        assert_eq!(
            ledger.recovery_edge_count(),
            success.report().problem_size().scalar_hard_relations()
        );
        assert_eq!(
            ledger.solver_relation_row_count(),
            success
                .report()
                .problem_size()
                .canonical_hard_equalities()
                .unwrap()
                + success
                    .report()
                    .problem_size()
                    .affine_inequality_constraints()
        );
        for source in ["close-gradient-a", "close-gradient-b"] {
            let assessments = success
                .report()
                .hard_relations()
                .iter()
                .filter(|relation| relation.source_id().as_str() == source)
                .collect::<Vec<_>>();
            assert_eq!(assessments.len(), 3);
            assert!(
                assessments
                    .iter()
                    .all(|relation| relation.residual().abs() <= relation.tolerance())
            );
        }
    }
}

#[test]
fn same_normal_direction_at_different_supports_retains_both_sources() {
    let mut builder = ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field"),
    );
    for (index, location) in [
        point(0.0, 0.0, 0.0),
        point(1.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0),
        point(0.0, 0.0, 1.0),
        point(1.0, 1.0, 1.0),
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("normal-value-{index}")),
                    location,
                    affine_value(location),
                )
                .unwrap(),
            )
            .unwrap();
    }
    for (source, location) in [
        ("close-normal-a", point(0.25, 0.5, 0.75)),
        ("close-normal-b", point(0.255, 0.5, 0.75)),
    ] {
        builder
            .add(
                DirectedNormalObservation::try_new(
                    SourceId::new(source),
                    location,
                    Vector3::try_new(1.0, 2.0, -0.5).unwrap(),
                    MinimumNormalSlope::try_new(1.0).unwrap(),
                )
                .unwrap(),
            )
            .unwrap();
    }

    let success = builder.build().unwrap().fit().unwrap();
    let report = success.report();

    assert_eq!(report.problem_size().equality_constraints(), Some(9));
    assert_eq!(report.problem_size().affine_inequality_constraints(), 2);
    assert_eq!(report.directed_normals().len(), 2);
    assert_eq!(
        report.directed_normals()[0].source_id().as_str(),
        "close-normal-a"
    );
    assert_eq!(
        report.directed_normals()[1].source_id().as_str(),
        "close-normal-b"
    );
}
