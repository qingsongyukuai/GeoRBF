use georbf::diagnostics::ProblemDiagnosis;
use georbf::fit::{BoundActiveState, BoundSide};
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::FieldValueObservation;
use georbf::observation::QuadraticPenalty;
use georbf::problem::BuildError;
use georbf::relation::SharedLevelSetBuilder;
use georbf::relation::{
    AdditiveFieldGauge, FieldSeparationInterval, FieldSeparationIntervalError,
    FieldSeparationViolationPenalty, LinearViolationPenalty, MinimumFieldOffset,
    MinimumFieldOffsetError, PointToLevelSetRelation, PointToLevelSetSide,
};
use georbf::{GroupId, Point3, ProblemBuilder, SourceId};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field"),
    )
}

fn level(group: &str, source: &str, location: Point3) -> georbf::relation::SharedLevelSet {
    let mut level = SharedLevelSetBuilder::new(GroupId::new(group));
    level.add_member(SourceId::new(source), location).unwrap();
    level.build().unwrap()
}

#[test]
fn field_separation_interval_checks_ordered_signed_bounds() {
    let reference = GroupId::new("reference");
    let target = GroupId::new("target");

    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            FieldSeparationInterval::try_hard(
                SourceId::new("non-finite-lower"),
                reference.clone(),
                target.clone(),
                value,
                1.0,
            ),
            Err(FieldSeparationIntervalError::NonFiniteBound)
        );
        assert_eq!(
            FieldSeparationInterval::try_hard(
                SourceId::new("non-finite-upper"),
                reference.clone(),
                target.clone(),
                -1.0,
                value,
            ),
            Err(FieldSeparationIntervalError::NonFiniteBound)
        );
    }
    assert_eq!(
        FieldSeparationInterval::try_hard(
            SourceId::new("empty"),
            reference.clone(),
            target.clone(),
            2.0,
            1.0,
        ),
        Err(FieldSeparationIntervalError::EmptyInterval {
            lower: 2.0,
            upper: 1.0,
        })
    );
    assert_eq!(
        FieldSeparationInterval::try_hard(
            SourceId::new("self"),
            reference.clone(),
            reference.clone(),
            -1.0,
            1.0,
        ),
        Err(FieldSeparationIntervalError::SelfReference {
            group_id: reference.clone(),
        })
    );

    let interval = FieldSeparationInterval::try_hard(
        SourceId::new("signed"),
        reference.clone(),
        target.clone(),
        -0.0,
        2.5,
    )
    .unwrap();
    assert_eq!(interval.reference_group_id(), &reference);
    assert_eq!(interval.target_group_id(), &target);
    assert_eq!(interval.lower_bound().to_bits(), 0.0_f64.to_bits());
    assert_eq!(interval.upper_bound(), 2.5);
    assert!(!interval.is_soft());
}

#[test]
fn point_to_level_set_relation_requires_an_explicit_side_and_positive_offset() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            MinimumFieldOffset::try_new(value),
            Err(MinimumFieldOffsetError::NotFinite)
        );
    }
    for value in [-1.0, -0.0, 0.0] {
        assert_eq!(
            MinimumFieldOffset::try_new(value),
            Err(MinimumFieldOffsetError::NotPositive)
        );
    }

    let offset = MinimumFieldOffset::try_new(1.25).unwrap();
    let group_id = GroupId::new("level");
    let increasing = PointToLevelSetRelation::hard(
        SourceId::new("increasing"),
        point(1.0, 2.0, 3.0),
        group_id.clone(),
        PointToLevelSetSide::Increasing,
        offset,
    );
    assert_eq!(increasing.group_id(), &group_id);
    assert_eq!(increasing.location(), point(1.0, 2.0, 3.0));
    assert_eq!(increasing.side(), PointToLevelSetSide::Increasing);
    assert_eq!(increasing.minimum_offset(), offset);
    assert!(!increasing.is_soft());

    let decreasing = PointToLevelSetRelation::hard(
        SourceId::new("decreasing"),
        point(-1.0, -2.0, -3.0),
        group_id,
        PointToLevelSetSide::Decreasing,
        offset,
    );
    assert_eq!(decreasing.side(), PointToLevelSetSide::Decreasing);
}

#[test]
fn affine_level_set_soft_relations_expose_only_positive_legal_losses() {
    let quadratic = QuadraticPenalty::try_new(2.0).unwrap();
    let linear = LinearViolationPenalty::try_new(3.0).unwrap();
    let interval = FieldSeparationInterval::try_with_violation_penalties(
        SourceId::new("soft-interval"),
        GroupId::new("reference"),
        GroupId::new("target"),
        -1.0,
        FieldSeparationViolationPenalty::Quadratic(quadratic),
        2.0,
        FieldSeparationViolationPenalty::Linear(linear),
    )
    .unwrap();
    assert!(interval.is_soft());
    assert!(
        FieldSeparationInterval::try_with_quadratic_penalties(
            SourceId::new("quadratic-interval"),
            GroupId::new("reference"),
            GroupId::new("target"),
            -1.0,
            quadratic,
            2.0,
            quadratic,
        )
        .unwrap()
        .is_soft()
    );
    assert!(
        FieldSeparationInterval::try_with_linear_violation_penalties(
            SourceId::new("linear-interval"),
            GroupId::new("reference"),
            GroupId::new("target"),
            -1.0,
            linear,
            2.0,
            linear,
        )
        .unwrap()
        .is_soft()
    );

    let offset = MinimumFieldOffset::try_new(0.5).unwrap();
    let quadratic_relation = PointToLevelSetRelation::with_quadratic_penalty(
        SourceId::new("soft-point-quadratic"),
        point(0.0, 0.0, 0.0),
        GroupId::new("level"),
        PointToLevelSetSide::Increasing,
        offset,
        quadratic,
    );
    assert!(quadratic_relation.is_soft());
    let linear_relation = PointToLevelSetRelation::with_linear_violation_penalty(
        SourceId::new("soft-point-linear"),
        point(0.0, 0.0, 0.0),
        GroupId::new("level"),
        PointToLevelSetSide::Decreasing,
        offset,
        linear,
    );
    assert!(linear_relation.is_soft());
}

#[test]
fn affine_level_set_relations_allow_forward_references_and_report_dangling_groups() {
    let offset = MinimumFieldOffset::try_new(0.5).unwrap();
    let mut problem = builder();
    problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("b-separation"),
                GroupId::new("reference"),
                GroupId::new("missing-target"),
                -1.0,
                2.0,
            )
            .unwrap(),
        )
        .unwrap();
    assert!(
        problem
            .add(PointToLevelSetRelation::hard(
                SourceId::new("b-separation"),
                point(2.0, 0.0, 0.0),
                GroupId::new("not-inserted"),
                PointToLevelSetSide::Increasing,
                offset,
            ))
            .is_err()
    );
    problem
        .add(PointToLevelSetRelation::hard(
            SourceId::new("a-point-side"),
            point(1.0, 0.0, 0.0),
            GroupId::new("missing-point-level"),
            PointToLevelSetSide::Increasing,
            offset,
        ))
        .unwrap();
    problem
        .add(level("reference", "reference/member", point(0.0, 0.0, 0.0)))
        .unwrap();

    let failure = problem.build().unwrap_err();
    assert_eq!(
        failure.errors(),
        &[
            BuildError::UnknownGroupReference {
                source_id: SourceId::new("a-point-side"),
                group_id: GroupId::new("missing-point-level"),
            },
            BuildError::UnknownGroupReference {
                source_id: SourceId::new("b-separation"),
                group_id: GroupId::new("missing-target"),
            },
        ]
    );

    let mut problem = failure.into_builder();
    problem
        .add(level(
            "missing-point-level",
            "point-level/member",
            point(0.0, 1.0, 0.0),
        ))
        .unwrap();
    problem
        .add(level(
            "missing-target",
            "target/member",
            point(0.0, 0.0, 1.0),
        ))
        .unwrap();
    let snapshot = problem.build().unwrap();
    assert_eq!(snapshot.field_separation_interval_count(), 1);
    assert_eq!(snapshot.point_to_level_set_relation_count(), 1);
}

#[test]
fn hard_affine_level_set_relations_fit_and_report_original_unit_observables() {
    let mut problem = builder();
    for (source, location, value) in [
        ("origin", point(0.0, 0.0, 0.0), 0.0),
        ("east", point(1.0, 0.0, 0.0), 1.0),
        ("north", point(0.0, 1.0, 0.0), 2.0),
        ("up", point(0.0, 0.0, 1.0), 3.0),
        ("target-value", point(1.0, 1.0, 0.0), 3.0),
    ] {
        problem
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }
    problem
        .add(level("reference", "reference/member", point(0.0, 0.0, 0.0)))
        .unwrap();
    problem
        .add(level("target", "target/member", point(1.0, 1.0, 0.0)))
        .unwrap();
    problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("separation"),
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
            point(0.0, 0.0, 1.0),
            GroupId::new("reference"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(2.5).unwrap(),
        ))
        .unwrap();

    let success = problem.build().unwrap().fit().unwrap();
    assert!(
        (success
            .model()
            .shared_level_value(&GroupId::new("reference"))
            .unwrap()
            - 0.0)
            .abs()
            <= 1.0e-7
    );
    assert!(
        (success
            .model()
            .shared_level_value(&GroupId::new("target"))
            .unwrap()
            - 3.0)
            .abs()
            <= 1.0e-7
    );

    let separation = success.report().field_separation_intervals();
    assert_eq!(separation.len(), 2);
    for assessment in separation {
        assert_eq!(assessment.source_id(), &SourceId::new("separation"));
        assert_eq!(assessment.reference_group_id(), &GroupId::new("reference"));
        assert_eq!(assessment.target_group_id(), &GroupId::new("target"));
        assert!((assessment.recovered_reference_value() - 0.0).abs() <= 1.0e-7);
        assert!((assessment.recovered_target_value() - 3.0).abs() <= 1.0e-7);
        assert!((assessment.recovered_field_separation() - 3.0).abs() <= 1.0e-7);
        assert_eq!(assessment.violation(), 0.0);
        assert_eq!(assessment.active_state(), BoundActiveState::Inactive);
        assert!(assessment.loss().is_none());
    }
    assert_eq!(separation[0].side(), BoundSide::Lower);
    assert_eq!(separation[0].bound(), 2.5);
    assert_eq!(separation[1].side(), BoundSide::Upper);
    assert_eq!(separation[1].bound(), 3.5);

    let point_side = &success.report().point_to_level_set_relations()[0];
    assert_eq!(point_side.source_id(), &SourceId::new("point-side"));
    assert_eq!(point_side.group_id(), &GroupId::new("reference"));
    assert_eq!(point_side.side(), PointToLevelSetSide::Increasing);
    assert_eq!(point_side.minimum_offset().value(), 2.5);
    assert!((point_side.recovered_point_value() - 3.0).abs() <= 1.0e-7);
    assert!((point_side.recovered_level_value() - 0.0).abs() <= 1.0e-7);
    assert!((point_side.recovered_field_offset() - 3.0).abs() <= 1.0e-7);
    assert_eq!(point_side.violation(), 0.0);
    assert_eq!(point_side.active_state(), BoundActiveState::Inactive);
    assert!(point_side.loss().is_none());
}

#[test]
fn soft_affine_level_set_relations_recover_independent_violations_and_losses() {
    let mut problem = builder();
    for (source, location, value) in [
        ("origin", point(0.0, 0.0, 0.0), 0.0),
        ("east", point(1.0, 0.0, 0.0), 1.0),
        ("north", point(0.0, 1.0, 0.0), 2.0),
        ("up", point(0.0, 0.0, 1.0), 3.0),
        ("target-value", point(1.0, 1.0, 0.0), 3.0),
    ] {
        problem
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }
    problem
        .add(level("reference", "reference/member", point(0.0, 0.0, 0.0)))
        .unwrap();
    problem
        .add(level("target", "target/member", point(1.0, 1.0, 0.0)))
        .unwrap();
    for (source, lower, upper) in [("high", 4.0, 5.0), ("low", -5.0, -4.0)] {
        problem
            .add(
                FieldSeparationInterval::try_with_violation_penalties(
                    SourceId::new(source),
                    GroupId::new("reference"),
                    GroupId::new("target"),
                    lower,
                    FieldSeparationViolationPenalty::Quadratic(
                        QuadraticPenalty::try_new(2.0).unwrap(),
                    ),
                    upper,
                    FieldSeparationViolationPenalty::Linear(
                        LinearViolationPenalty::try_new(3.0).unwrap(),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
    }
    problem
        .add(PointToLevelSetRelation::with_quadratic_penalty(
            SourceId::new("point-quadratic"),
            point(0.0, 0.0, 1.0),
            GroupId::new("reference"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(4.0).unwrap(),
            QuadraticPenalty::try_new(2.0).unwrap(),
        ))
        .unwrap();
    problem
        .add(PointToLevelSetRelation::with_linear_violation_penalty(
            SourceId::new("point-linear"),
            point(0.0, 0.0, 1.0),
            GroupId::new("reference"),
            PointToLevelSetSide::Decreasing,
            MinimumFieldOffset::try_new(1.0).unwrap(),
            LinearViolationPenalty::try_new(3.0).unwrap(),
        ))
        .unwrap();
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let success = problem.build().unwrap().fit().unwrap();
    let report = success.report();
    assert_eq!(report.field_separation_intervals().len(), 4);
    let interval_side = |source: &str, side: BoundSide| {
        report
            .field_separation_intervals()
            .iter()
            .find(|assessment| {
                assessment.source_id() == &SourceId::new(source) && assessment.side() == side
            })
            .unwrap()
    };
    let high_lower = interval_side("high", BoundSide::Lower);
    assert!((high_lower.violation() - 1.0).abs() <= 5.0e-4);
    assert!((high_lower.loss().unwrap() - 1.0).abs() <= 5.0e-4);
    assert!(interval_side("high", BoundSide::Upper).violation() <= 5.0e-6);
    let low_upper = interval_side("low", BoundSide::Upper);
    assert!((low_upper.violation() - 7.0).abs() <= 5.0e-4);
    assert!((low_upper.loss().unwrap() - 21.0).abs() <= 5.0e-4);

    assert_eq!(report.point_to_level_set_relations().len(), 2);
    let point_relation = |source: &str| {
        report
            .point_to_level_set_relations()
            .iter()
            .find(|assessment| assessment.source_id() == &SourceId::new(source))
            .unwrap()
    };
    let quadratic = point_relation("point-quadratic");
    assert!((quadratic.violation() - 1.0).abs() <= 5.0e-4);
    assert!((quadratic.loss().unwrap() - 1.0).abs() <= 5.0e-4);
    let linear = point_relation("point-linear");
    assert!((linear.violation() - 4.0).abs() <= 5.0e-4);
    assert!((linear.loss().unwrap() - 12.0).abs() <= 5.0e-4);
    assert!(report.field_energy().unwrap().is_finite());
    assert!(report.total_objective().unwrap() >= 34.99);
}

#[test]
fn provable_affine_level_set_graph_conflicts_stop_before_the_backend() {
    let mut interval_problem = builder();
    interval_problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("a-to-b"),
                GroupId::new("a"),
                GroupId::new("b"),
                1.0,
                2.0,
            )
            .unwrap(),
        )
        .unwrap();
    interval_problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("b-to-a"),
                GroupId::new("b"),
                GroupId::new("a"),
                0.0,
                1.0,
            )
            .unwrap(),
        )
        .unwrap();
    interval_problem
        .add(level("a", "a/member", point(0.0, 0.0, 0.0)))
        .unwrap();
    interval_problem
        .add(level("b", "b/member", point(1.0, 0.0, 0.0)))
        .unwrap();
    let failure = interval_problem.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
    let conflict = failure
        .report()
        .shared_level_set_relation_conflict()
        .unwrap();
    assert_eq!(
        conflict.source_ids(),
        &[SourceId::new("a-to-b"), SourceId::new("b-to-a")]
    );

    let mut point_problem = builder();
    point_problem
        .add(PointToLevelSetRelation::hard(
            SourceId::new("impossible-side"),
            point(0.0, 0.0, 0.0),
            GroupId::new("level"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(1.0).unwrap(),
        ))
        .unwrap();
    point_problem
        .add(level("level", "level/member", point(0.0, 0.0, 0.0)))
        .unwrap();
    let failure = point_problem.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
    let conflict = failure
        .report()
        .shared_level_set_relation_conflict()
        .unwrap();
    assert_eq!(
        conflict.source_ids(),
        &[
            SourceId::new("impossible-side"),
            SourceId::new("level/member")
        ]
    );
}

fn covariant_affine_level_set_fit(length_scale: f64, field_scale: f64) -> georbf::fit::FitSuccess {
    let scale_point =
        |[x, y, z]: [f64; 3]| point(length_scale * x, length_scale * y, length_scale * z);
    let mut problem = ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new(if length_scale == 1.0 { "m" } else { "scaled-m" }),
        )
        .unwrap(),
        FieldUnitLabel::new(if field_scale == 1.0 {
            "field"
        } else {
            "scaled-field"
        }),
    );
    for (source, coordinates, value) in [
        ("origin", [0.0, 0.0, 0.0], 0.0),
        ("east", [1.0, 0.0, 0.0], 1.0),
        ("north", [0.0, 1.0, 0.0], 2.0),
        ("up", [0.0, 0.0, 1.0], 3.0),
        ("target-value", [1.0, 1.0, 0.0], 3.0),
    ] {
        problem
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(source),
                    scale_point(coordinates),
                    field_scale * value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    problem
        .add(level(
            "reference",
            "reference/member",
            scale_point([0.0, 0.0, 0.0]),
        ))
        .unwrap();
    problem
        .add(level(
            "target",
            "target/member",
            scale_point([1.0, 1.0, 0.0]),
        ))
        .unwrap();
    problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("separation"),
                GroupId::new("reference"),
                GroupId::new("target"),
                field_scale * 2.5,
                field_scale * 3.5,
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(PointToLevelSetRelation::hard(
            SourceId::new("point-side"),
            scale_point([0.0, 0.0, 1.0]),
            GroupId::new("reference"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(field_scale * 2.5).unwrap(),
        ))
        .unwrap();
    problem.build().unwrap().fit().unwrap()
}

#[test]
fn frame_changes_leave_field_offsets_unchanged_and_field_rescaling_is_covariant() {
    let original = covariant_affine_level_set_fit(1.0, 1.0);
    let length_changed = covariant_affine_level_set_fit(2.5, 1.0);
    let field_changed = covariant_affine_level_set_fit(1.0, 4.0);

    for (transformed, field_scale) in [(&length_changed, 1.0), (&field_changed, 4.0)] {
        for (left, right) in original
            .report()
            .field_separation_intervals()
            .iter()
            .zip(transformed.report().field_separation_intervals())
        {
            assert_eq!(left.side(), right.side());
            assert!((right.bound() - field_scale * left.bound()).abs() <= 1.0e-8);
            assert!(
                (right.recovered_field_separation()
                    - field_scale * left.recovered_field_separation())
                .abs()
                    <= 1.0e-7
            );
        }
        let left = &original.report().point_to_level_set_relations()[0];
        let right = &transformed.report().point_to_level_set_relations()[0];
        assert_eq!(left.side(), right.side());
        assert!(
            (right.minimum_offset().value() - field_scale * left.minimum_offset().value()).abs()
                <= 1.0e-8
        );
        assert!(
            (right.recovered_field_offset() - field_scale * left.recovered_field_offset()).abs()
                <= 1.0e-7
        );
    }
}

#[test]
fn swapping_interval_roles_is_equivalent_only_with_the_signed_interval_transform() {
    let original = FieldSeparationInterval::try_hard(
        SourceId::new("original"),
        GroupId::new("reference"),
        GroupId::new("target"),
        -2.0,
        5.0,
    )
    .unwrap();
    let swapped = FieldSeparationInterval::try_hard(
        SourceId::new("swapped"),
        GroupId::new("target"),
        GroupId::new("reference"),
        -5.0,
        2.0,
    )
    .unwrap();
    assert_eq!(swapped.lower_bound(), -original.upper_bound());
    assert_eq!(swapped.upper_bound(), -original.lower_bound());
    assert_ne!(swapped.lower_bound(), original.lower_bound());
}

fn permuted_affine_level_set_fit(reverse: bool) -> georbf::fit::FitSuccess {
    let mut problem = builder();
    let mut observations = vec![
        ("origin", point(0.0, 0.0, 0.0), 0.0),
        ("east", point(1.0, 0.0, 0.0), 1.0),
        ("north", point(0.0, 1.0, 0.0), 2.0),
        ("up", point(0.0, 0.0, 1.0), 3.0),
        ("target-value", point(1.0, 1.0, 0.0), 3.0),
    ];
    let mut levels = vec![
        level("reference", "reference/member", point(0.0, 0.0, 0.0)),
        level("target", "target/member", point(1.0, 1.0, 0.0)),
    ];
    let mut intervals = vec![
        FieldSeparationInterval::try_hard(
            SourceId::new("separation-a"),
            GroupId::new("reference"),
            GroupId::new("target"),
            2.5,
            3.5,
        )
        .unwrap(),
        FieldSeparationInterval::try_hard(
            SourceId::new("separation-b"),
            GroupId::new("reference"),
            GroupId::new("target"),
            2.0,
            4.0,
        )
        .unwrap(),
    ];
    let offset = MinimumFieldOffset::try_new(2.0).unwrap();
    let mut point_relations = vec![
        PointToLevelSetRelation::hard(
            SourceId::new("point-a"),
            point(0.0, 0.0, 1.0),
            GroupId::new("reference"),
            PointToLevelSetSide::Increasing,
            offset,
        ),
        PointToLevelSetRelation::hard(
            SourceId::new("point-b"),
            point(0.0, 0.0, 1.0),
            GroupId::new("reference"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(2.5).unwrap(),
        ),
    ];
    if reverse {
        observations.reverse();
        levels.reverse();
        intervals.reverse();
        point_relations.reverse();
    }
    for (source, location, value) in observations {
        problem
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }
    for level in levels {
        problem.add(level).unwrap();
    }
    for interval in intervals {
        problem.add(interval).unwrap();
    }
    for relation in point_relations {
        problem.add(relation).unwrap();
    }
    problem.build().unwrap().fit().unwrap()
}

#[test]
fn input_permutation_preserves_affine_level_set_results_and_report_order() {
    let forward = permuted_affine_level_set_fit(false);
    let reverse = permuted_affine_level_set_fit(true);
    assert_eq!(
        forward
            .report()
            .field_separation_intervals()
            .iter()
            .map(|assessment| (assessment.source_id().clone(), assessment.side()))
            .collect::<Vec<_>>(),
        reverse
            .report()
            .field_separation_intervals()
            .iter()
            .map(|assessment| (assessment.source_id().clone(), assessment.side()))
            .collect::<Vec<_>>()
    );
    assert_eq!(
        forward
            .report()
            .point_to_level_set_relations()
            .iter()
            .map(|assessment| assessment.source_id().clone())
            .collect::<Vec<_>>(),
        reverse
            .report()
            .point_to_level_set_relations()
            .iter()
            .map(|assessment| assessment.source_id().clone())
            .collect::<Vec<_>>()
    );
    let query = point(0.25, -0.5, 0.75);
    let left = forward.model().evaluate(query).unwrap();
    let right = reverse.model().evaluate(query).unwrap();
    assert!((left.value() - right.value()).abs() <= 1.0e-10);
    for (left, right) in left
        .gradient()
        .components()
        .into_iter()
        .zip(right.gradient().components())
    {
        assert!((left - right).abs() <= 1.0e-10);
    }
}

#[test]
fn an_explicit_level_set_gauge_closes_a_difference_only_affine_problem() {
    let mut reference = SharedLevelSetBuilder::new(GroupId::new("reference"));
    for (source, location) in [
        ("reference/origin", point(0.0, 0.0, 0.0)),
        ("reference/east", point(1.0, 0.0, 0.0)),
        ("reference/north", point(0.0, 1.0, 0.0)),
        ("reference/up", point(0.0, 0.0, 1.0)),
    ] {
        reference
            .add_member(SourceId::new(source), location)
            .unwrap();
    }
    let mut problem = builder();
    problem.add(reference.build().unwrap()).unwrap();
    problem
        .add(level("target", "target/member", point(1.0, 1.0, 1.0)))
        .unwrap();
    problem
        .add(
            FieldSeparationInterval::try_hard(
                SourceId::new("exact-separation"),
                GroupId::new("reference"),
                GroupId::new("target"),
                1.0,
                1.0,
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(PointToLevelSetRelation::hard(
            SourceId::new("target-side"),
            point(1.0, 1.0, 1.0),
            GroupId::new("reference"),
            PointToLevelSetSide::Increasing,
            MinimumFieldOffset::try_new(0.5).unwrap(),
        ))
        .unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(
                SourceId::new("reference-gauge"),
                GroupId::new("reference"),
                0.0,
            )
            .unwrap(),
        )
        .unwrap();

    let success = problem.build().unwrap().fit().unwrap();
    assert!(
        success
            .model()
            .shared_level_value(&GroupId::new("reference"))
            .unwrap()
            .abs()
            <= 1.0e-7
    );
    assert!(
        (success
            .model()
            .shared_level_value(&GroupId::new("target"))
            .unwrap()
            - 1.0)
            .abs()
            <= 1.0e-7
    );
}
