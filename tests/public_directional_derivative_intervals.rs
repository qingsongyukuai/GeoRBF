use georbf::diagnostics::ProblemDiagnosis;
use georbf::fit::{BoundActiveState, BoundSide};
use georbf::geometry::GlobalAnisotropyMetric;
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{
    FieldValueObservation, GradientObservation, QuadraticPenalty, TangentDirectionObservation,
};
use georbf::relation::{
    DirectionalDerivativeInterval, DirectionalDerivativeIntervalError,
    DirectionalDerivativeViolationPenalty, LinearViolationPenalty,
};
use georbf::{Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).unwrap()
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

#[test]
fn checked_directional_derivative_constructors_normalize_and_reject_invalid_inputs() {
    let location = point(1.0, 2.0, 3.0);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            DirectionalDerivativeInterval::try_lower(
                SourceId::new("lower"),
                location,
                vector(2.0, -3.0, 6.0),
                value,
            ),
            Err(DirectionalDerivativeIntervalError::NonFiniteBound)
        );
        assert_eq!(
            DirectionalDerivativeInterval::try_upper(
                SourceId::new("upper"),
                location,
                vector(2.0, -3.0, 6.0),
                value,
            ),
            Err(DirectionalDerivativeIntervalError::NonFiniteBound)
        );
        assert_eq!(
            DirectionalDerivativeInterval::try_interval(
                SourceId::new("interval-lower"),
                location,
                vector(2.0, -3.0, 6.0),
                value,
                1.0,
            ),
            Err(DirectionalDerivativeIntervalError::NonFiniteBound)
        );
        assert_eq!(
            DirectionalDerivativeInterval::try_interval(
                SourceId::new("interval-upper"),
                location,
                vector(2.0, -3.0, 6.0),
                -1.0,
                value,
            ),
            Err(DirectionalDerivativeIntervalError::NonFiniteBound)
        );
    }
    assert_eq!(
        DirectionalDerivativeInterval::try_interval(
            SourceId::new("empty"),
            location,
            vector(1.0, 0.0, 0.0),
            2.0,
            1.0,
        ),
        Err(DirectionalDerivativeIntervalError::EmptyInterval {
            lower: 2.0,
            upper: 1.0,
        })
    );
    assert_eq!(
        DirectionalDerivativeInterval::try_lower(
            SourceId::new("zero-direction"),
            location,
            vector(0.0, -0.0, 0.0),
            0.0,
        ),
        Err(DirectionalDerivativeIntervalError::ZeroDirection)
    );

    let lower = DirectionalDerivativeInterval::try_lower(
        SourceId::new("normalized"),
        location,
        vector(2.0, -3.0, 6.0),
        -0.0,
    )
    .unwrap();
    assert_eq!(lower.location(), location);
    assert_eq!(
        lower.direction().components(),
        [2.0 / 7.0, -3.0 / 7.0, 6.0 / 7.0]
    );
    assert_eq!(lower.lower_bound().unwrap().to_bits(), 0.0_f64.to_bits());
    assert_eq!(lower.upper_bound(), None);

    let scaled = DirectionalDerivativeInterval::try_lower(
        SourceId::new("scaled"),
        location,
        vector(20.0, -30.0, 60.0),
        0.0,
    )
    .unwrap();
    assert_eq!(scaled.direction(), lower.direction());

    let opposite = DirectionalDerivativeInterval::try_lower(
        SourceId::new("opposite"),
        location,
        vector(-2.0, 3.0, -6.0),
        0.0,
    )
    .unwrap();
    assert_eq!(
        opposite.direction().components(),
        [-2.0 / 7.0, 3.0 / 7.0, -6.0 / 7.0]
    );

    for direction in [
        vector(f64::MAX, -f64::MAX, 0.0),
        vector(f64::MIN_POSITIVE, f64::MIN_POSITIVE, 0.0),
        vector(f64::from_bits(1), 0.0, 0.0),
    ] {
        let unit = DirectionalDerivativeInterval::try_upper(
            SourceId::new("extreme"),
            location,
            direction,
            1.0,
        )
        .expect("every finite nonzero magnitude normalizes without overflow or underflow")
        .direction()
        .components();
        assert!(unit.iter().all(|component| component.is_finite()));
        let norm = unit
            .into_iter()
            .map(|component| component * component)
            .sum::<f64>()
            .sqrt();
        assert!((norm - 1.0).abs() <= 1.0e-15);
    }
}

#[test]
fn typed_soft_sides_and_builder_insertion_remain_independent_and_atomic() {
    let location = point(0.0, 0.0, 0.0);
    let direction = vector(1.0, 2.0, 2.0);
    let interval = DirectionalDerivativeInterval::try_interval_with_violation_penalties(
        SourceId::new("soft-interval"),
        location,
        direction,
        -1.0,
        DirectionalDerivativeViolationPenalty::Quadratic(QuadraticPenalty::try_new(2.0).unwrap()),
        3.0,
        DirectionalDerivativeViolationPenalty::Linear(
            LinearViolationPenalty::try_new(4.0).unwrap(),
        ),
    )
    .unwrap();
    assert!(interval.is_soft());
    assert_eq!(
        interval.direction().components(),
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]
    );
    assert_eq!(interval.lower_bound(), Some(-1.0));
    assert_eq!(interval.upper_bound(), Some(3.0));

    let mut builder = builder();
    builder
        .add(FieldValueObservation::try_new(SourceId::new("value"), location, 0.0).unwrap())
        .unwrap();
    builder
        .add(
            DirectionalDerivativeInterval::try_lower(
                SourceId::new("derivative"),
                location,
                direction,
                -1.0,
            )
            .unwrap(),
        )
        .unwrap();
    assert!(
        builder
            .add(
                DirectionalDerivativeInterval::try_upper(
                    SourceId::new("derivative"),
                    location,
                    direction,
                    1.0,
                )
                .unwrap(),
            )
            .is_err()
    );
    let snapshot = builder.build().unwrap();
    assert_eq!(snapshot.observation_count(), 1);
    assert_eq!(snapshot.directional_derivative_interval_count(), 1);
    assert_eq!(snapshot.source_count(), 2);
}

#[test]
fn hard_lower_upper_and_interval_fit_through_the_public_qp_path() {
    let mut builder = builder();
    for (source, location, value) in [
        ("origin", point(0.0, 0.0, 0.0), 0.0),
        ("east", point(1.0, 0.0, 0.0), 1.0),
        ("north", point(0.0, 1.0, 0.0), 2.0),
        ("up", point(0.0, 0.0, 1.0), 3.0),
    ] {
        builder
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }
    for relation in [
        DirectionalDerivativeInterval::try_lower(
            SourceId::new("hard-lower"),
            point(0.2, 0.3, 0.4),
            vector(1.0, 0.0, 0.0),
            0.9,
        )
        .unwrap(),
        DirectionalDerivativeInterval::try_upper(
            SourceId::new("hard-upper"),
            point(-0.2, 0.1, 0.6),
            vector(0.0, 5.0, 0.0),
            2.1,
        )
        .unwrap(),
        DirectionalDerivativeInterval::try_interval(
            SourceId::new("hard-interval"),
            point(0.7, -0.4, 0.3),
            vector(0.0, 0.0, 2.0),
            2.9,
            3.0,
        )
        .unwrap(),
    ] {
        builder.add(relation).unwrap();
    }

    let success = builder.build().unwrap().fit().unwrap();
    let sample = success.model().evaluate(point(0.25, -0.5, 0.75)).unwrap();
    assert!(
        (sample.value() - 1.5).abs() <= 5.0e-4,
        "recovered value {}",
        sample.value()
    );
    for (actual, expected) in sample
        .gradient()
        .components()
        .into_iter()
        .zip([1.0, 2.0, 3.0])
    {
        assert!((actual - expected).abs() <= 5.0e-4);
    }

    let assessments = success.report().directional_derivative_intervals();
    assert_eq!(assessments.len(), 4);
    let side = |source: &str, expected_side: BoundSide| {
        assessments
            .iter()
            .find(|assessment| {
                assessment.source_id().as_str() == source && assessment.side() == expected_side
            })
            .unwrap()
    };
    let lower = side("hard-lower", BoundSide::Lower);
    assert_eq!(
        lower.semantic_role().as_str(),
        "directional-derivative-interval/lower"
    );
    assert_eq!(lower.direction(), vector(1.0, 0.0, 0.0));
    assert!((lower.recovered_directional_derivative() - 1.0).abs() <= 5.0e-4);
    assert!((lower.slack() - 0.1).abs() <= 5.0e-4);
    assert_eq!(lower.violation(), 0.0);
    assert_eq!(lower.active_state(), BoundActiveState::Inactive);
    assert!(lower.loss().is_none());

    let upper = side("hard-upper", BoundSide::Upper);
    assert!((upper.recovered_directional_derivative() - 2.0).abs() <= 5.0e-4);
    assert!((upper.slack() - 0.1).abs() <= 5.0e-4);
    let interval_upper = side("hard-interval", BoundSide::Upper);
    assert!((interval_upper.recovered_directional_derivative() - 3.0).abs() <= 5.0e-4);
    assert_eq!(
        interval_upper.active_state(),
        if interval_upper.slack() <= interval_upper.tolerance() {
            BoundActiveState::Active
        } else {
            BoundActiveState::Inactive
        }
    );
}

#[test]
fn soft_directional_derivative_sides_recover_independent_physical_losses() {
    let mut missing_normalization = builder();
    missing_normalization
        .add(
            DirectionalDerivativeInterval::try_lower_with_quadratic_penalty(
                SourceId::new("soft"),
                point(0.0, 0.0, 0.0),
                vector(1.0, 0.0, 0.0),
                1.0,
                QuadraticPenalty::try_new(2.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    assert!(missing_normalization.build().is_err());

    let mut builder = builder();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();
    for (source, location) in [
        ("origin", point(0.0, 0.0, 0.0)),
        ("east", point(1.0, 0.0, 0.0)),
        ("north", point(0.0, 1.0, 0.0)),
        ("up", point(0.0, 0.0, 1.0)),
    ] {
        builder
            .add(FieldValueObservation::try_new(SourceId::new(source), location, 0.0).unwrap())
            .unwrap();
    }
    let location = point(0.2, -0.3, 0.4);
    let direction = vector(2.0, 0.0, 0.0);
    builder
        .add(GradientObservation::new(
            SourceId::new("gradient-anchor"),
            location,
            vector(0.0, 0.0, 0.0),
        ))
        .unwrap();
    for relation in [
        DirectionalDerivativeInterval::try_lower_with_quadratic_penalty(
            SourceId::new("soft-lower"),
            location,
            direction,
            1.0,
            QuadraticPenalty::try_new(2.0).unwrap(),
        )
        .unwrap(),
        DirectionalDerivativeInterval::try_lower_with_quadratic_penalty(
            SourceId::new("soft-lower-duplicate"),
            location,
            vector(20.0, 0.0, 0.0),
            1.0,
            QuadraticPenalty::try_new(2.0).unwrap(),
        )
        .unwrap(),
        DirectionalDerivativeInterval::try_upper_with_linear_violation_penalty(
            SourceId::new("soft-upper"),
            location,
            direction,
            -2.0,
            LinearViolationPenalty::try_new(3.0).unwrap(),
        )
        .unwrap(),
        DirectionalDerivativeInterval::try_interval_with_quadratic_penalties(
            SourceId::new("soft-interval"),
            location,
            direction,
            1.0,
            QuadraticPenalty::try_new(4.0).unwrap(),
            2.0,
            QuadraticPenalty::try_new(5.0).unwrap(),
        )
        .unwrap(),
    ] {
        builder.add(relation).unwrap();
    }

    let success = builder.build().unwrap().fit().unwrap();
    let report = success.report();
    assert_eq!(report.directional_derivative_intervals().len(), 5);
    assert_eq!(report.problem_size().auxiliary_variables(), 5);
    assert_eq!(report.problem_size().quadratic_objective_terms(), 4);
    assert_eq!(report.problem_size().linear_objective_terms(), 1);
    assert_eq!(report.problem_size().affine_inequality_constraints(), 10);

    let side = |source: &str, expected_side: BoundSide| {
        report
            .directional_derivative_intervals()
            .iter()
            .find(|assessment| {
                assessment.source_id().as_str() == source && assessment.side() == expected_side
            })
            .unwrap()
    };
    for source in ["soft-lower", "soft-lower-duplicate"] {
        let assessment = side(source, BoundSide::Lower);
        assert!(
            (assessment.violation() - 1.0).abs() <= 1.0e-6,
            "{source} recovered derivative {}, violation {}",
            assessment.recovered_directional_derivative(),
            assessment.violation()
        );
        assert!((assessment.loss().unwrap() - 1.0).abs() <= 1.0e-5);
        assert_eq!(assessment.quadratic_penalty().unwrap().weight(), 2.0);
        assert!(assessment.linear_violation_penalty().is_none());
    }
    let upper = side("soft-upper", BoundSide::Upper);
    assert!((upper.violation() - 2.0).abs() <= 1.0e-6);
    assert!((upper.loss().unwrap() - 6.0).abs() <= 1.0e-5);
    assert_eq!(upper.linear_violation_penalty().unwrap().weight(), 3.0);

    let interval_lower = side("soft-interval", BoundSide::Lower);
    assert!((interval_lower.violation() - 1.0).abs() <= 1.0e-6);
    assert!((interval_lower.loss().unwrap() - 2.0).abs() <= 1.0e-5);
    let interval_upper = side("soft-interval", BoundSide::Upper);
    assert!(interval_upper.violation() <= 1.0e-4);
    assert!(interval_upper.slack() >= 1.9);
    assert!(interval_upper.loss().unwrap().abs() <= 1.0e-7);
    assert!((report.total_objective().unwrap() - 10.0).abs() <= 1.0e-4);
}

#[test]
fn exact_derivative_conflicts_are_preflighted_before_the_backend() {
    let location = point(0.0, 0.0, 0.0);
    let mut exact_conflict = builder();
    exact_conflict
        .add(FieldValueObservation::try_new(SourceId::new("value"), location, 0.0).unwrap())
        .unwrap();
    exact_conflict
        .add(
            TangentDirectionObservation::try_new(
                SourceId::new("tangent"),
                location,
                vector(2.0, 0.0, 0.0),
            )
            .unwrap(),
        )
        .unwrap();
    exact_conflict
        .add(
            DirectionalDerivativeInterval::try_lower(
                SourceId::new("lower"),
                location,
                vector(20.0, 0.0, 0.0),
                1.0,
            )
            .unwrap(),
        )
        .unwrap();
    let failure = exact_conflict.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
    assert_eq!(
        failure
            .report()
            .direct_input_conflict()
            .unwrap()
            .first_source()
            .as_str(),
        "lower"
    );

    let mut local_interval_conflict = builder();
    local_interval_conflict
        .add(
            DirectionalDerivativeInterval::try_lower(
                SourceId::new("lower"),
                location,
                vector(1.0, 2.0, 0.0),
                2.0,
            )
            .unwrap(),
        )
        .unwrap();
    local_interval_conflict
        .add(
            DirectionalDerivativeInterval::try_upper(
                SourceId::new("upper"),
                location,
                vector(5.0, 10.0, 0.0),
                1.0,
            )
            .unwrap(),
        )
        .unwrap();
    let failure = local_interval_conflict.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
}

fn exact_zero_derivative_fit(use_interval: bool) -> georbf::fit::FitSuccess {
    let mut builder = builder();
    for (source, location, value) in [
        ("origin", point(0.0, 0.0, 0.0), 0.0),
        ("east", point(1.0, 0.0, 0.0), 1.0),
        ("north", point(0.0, 1.0, 0.0), 2.0),
        ("up", point(0.0, 0.0, 1.0), 3.0),
    ] {
        builder
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }
    let location = point(0.2, -0.1, 0.4);
    let direction = vector(2.0, -1.0, 0.0);
    if use_interval {
        builder
            .add(
                DirectionalDerivativeInterval::try_interval(
                    SourceId::new("zero-derivative"),
                    location,
                    direction,
                    0.0,
                    0.0,
                )
                .unwrap(),
            )
            .unwrap();
    } else {
        builder
            .add(
                TangentDirectionObservation::try_new(
                    SourceId::new("zero-derivative"),
                    location,
                    direction,
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

#[test]
fn zero_width_hard_interval_and_tangent_have_equivalent_field_observables() {
    let tangent = exact_zero_derivative_fit(false);
    let interval = exact_zero_derivative_fit(true);
    for location in [
        point(0.1, 0.2, 0.3),
        point(-0.4, 0.7, 0.2),
        point(1.5, -0.2, 0.8),
    ] {
        let tangent_sample = tangent.model().evaluate(location).unwrap();
        let interval_sample = interval.model().evaluate(location).unwrap();
        assert!((tangent_sample.value() - interval_sample.value()).abs() <= 5.0e-6);
        for (left, right) in tangent_sample
            .gradient()
            .components()
            .into_iter()
            .zip(interval_sample.gradient().components())
        {
            assert!((left - right).abs() <= 5.0e-6);
        }
    }
    assert_eq!(tangent.report().hard_relations().len(), 5);
    assert!(
        tangent
            .report()
            .directional_derivative_intervals()
            .is_empty()
    );
    assert_eq!(interval.report().hard_relations().len(), 4);
    assert_eq!(
        interval.report().directional_derivative_intervals().len(),
        2
    );
}

fn duplicate_interval_fit(reverse: bool) -> georbf::fit::FitSuccess {
    let mut builder = builder();
    let mut observations = [
        ("origin", point(0.0, 0.0, 0.0), 0.0),
        ("east", point(1.0, 0.0, 0.0), 1.0),
        ("north", point(0.0, 1.0, 0.0), 2.0),
        ("up", point(0.0, 0.0, 1.0), 3.0),
    ];
    if reverse {
        observations.reverse();
    }
    for (source, location, value) in observations {
        builder
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }
    let mut relations = [
        DirectionalDerivativeInterval::try_upper(
            SourceId::new("duplicate-a"),
            point(0.2, 0.3, 0.4),
            vector(1.0, 2.0, 2.0),
            5.0,
        )
        .unwrap(),
        DirectionalDerivativeInterval::try_upper(
            SourceId::new("duplicate-b"),
            point(0.2, 0.3, 0.4),
            vector(10.0, 20.0, 20.0),
            5.0,
        )
        .unwrap(),
    ];
    if reverse {
        relations.reverse();
    }
    for relation in relations {
        builder.add(relation).unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

#[test]
fn scaled_duplicate_directions_share_hard_rows_and_reports_are_stably_ordered() {
    let forward = duplicate_interval_fit(false);
    let reverse = duplicate_interval_fit(true);
    for success in [&forward, &reverse] {
        assert_eq!(
            success
                .report()
                .problem_size()
                .affine_inequality_constraints(),
            1
        );
        assert_eq!(
            success
                .report()
                .directional_derivative_intervals()
                .iter()
                .map(|assessment| assessment.source_id().as_str())
                .collect::<Vec<_>>(),
            ["duplicate-a", "duplicate-b"]
        );
    }
    let query = point(0.6, -0.2, 0.1);
    let left = forward.model().evaluate(query).unwrap();
    let right = reverse.model().evaluate(query).unwrap();
    assert!((left.value() - right.value()).abs() <= 1.0e-10);
    assert_eq!(left.gradient(), right.gradient());
}

fn covariant_derivative_fit(length_scale: f64, field_scale: f64) -> georbf::fit::FitSuccess {
    let transformed = length_scale != 1.0 || field_scale != 1.0;
    let frame = InputCoordinateFrame::try_new(
        if transformed {
            ["north-prime", "east-prime", "up-prime"]
        } else {
            ["east", "north", "up"]
        },
        if transformed {
            Handedness::Left
        } else {
            Handedness::Right
        },
        LengthUnitLabel::new(if transformed { "scaled-m" } else { "m" }),
    )
    .unwrap();
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("field"));
    builder
        .set_global_anisotropy_metric(
            GlobalAnisotropyMetric::try_from_matrix(if transformed {
                [[0.25, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 1.0]]
            } else {
                [[4.0, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 1.0]]
            })
            .unwrap(),
        )
        .unwrap();
    builder
        .set_field_energy_normalization(
            FieldEnergyNormalization::try_new(length_scale.powi(3) / field_scale.powi(2)).unwrap(),
        )
        .unwrap();
    let transform = |[x, y, z]: [f64; 3]| {
        if transformed {
            [
                length_scale * y + 10.0,
                length_scale * x - 3.0,
                length_scale * z + 4.0,
            ]
        } else {
            [x, y, z]
        }
    };
    for (source, support, value) in [
        ("origin", [0.0, 0.0, 0.0], 0.0),
        ("east", [1.0, 0.0, 0.0], 1.0),
        ("north", [0.0, 1.0, 0.0], 2.0),
        ("up", [0.0, 0.0, 1.0], 3.0),
    ] {
        let [x, y, z] = transform(support);
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(source),
                    point(x, y, z),
                    field_scale * value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    let [x, y, z] = transform([0.2, -0.3, 0.4]);
    let derivative_scale = field_scale / length_scale;
    builder
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(x, y, z),
            if transformed {
                vector(
                    derivative_scale * 2.0,
                    derivative_scale,
                    derivative_scale * 3.0,
                )
            } else {
                vector(1.0, 2.0, 3.0)
            },
        ))
        .unwrap();
    let direction = if transformed {
        vector(0.0, 1.0, 0.0)
    } else {
        vector(1.0, 0.0, 0.0)
    };
    builder
        .add(
            DirectionalDerivativeInterval::try_lower_with_quadratic_penalty(
                SourceId::new("quadratic"),
                point(x, y, z),
                direction,
                derivative_scale * 2.0,
                QuadraticPenalty::try_new(2.0 / derivative_scale.powi(2)).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    builder
        .add(
            DirectionalDerivativeInterval::try_upper_with_linear_violation_penalty(
                SourceId::new("linear"),
                point(x, y, z),
                direction,
                0.0,
                LinearViolationPenalty::try_new(3.0 / derivative_scale).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    builder.build().unwrap().fit().unwrap()
}

#[test]
fn derivative_results_covary_under_rotation_reflection_scale_and_anisotropy() {
    let original = covariant_derivative_fit(1.0, 1.0);
    let transformed = covariant_derivative_fit(2.5, 4.0);
    let original_sample = original.model().evaluate(point(0.35, -0.2, 0.1)).unwrap();
    let transformed_sample = transformed
        .model()
        .evaluate(point(9.5, -2.125, 4.25))
        .unwrap();
    assert!((transformed_sample.value() - 4.0 * original_sample.value()).abs() <= 1.0e-6);
    let expected_gradient = [
        4.0 / 2.5 * original_sample.gradient().components()[1],
        4.0 / 2.5 * original_sample.gradient().components()[0],
        4.0 / 2.5 * original_sample.gradient().components()[2],
    ];
    for (actual, expected) in transformed_sample
        .gradient()
        .components()
        .into_iter()
        .zip(expected_gradient)
    {
        assert!((actual - expected).abs() <= 1.0e-6);
    }

    let derivative_scale = 4.0 / 2.5;
    for (left, right) in original
        .report()
        .directional_derivative_intervals()
        .iter()
        .zip(transformed.report().directional_derivative_intervals())
    {
        assert_eq!(left.source_id(), right.source_id());
        assert!((right.bound() - derivative_scale * left.bound()).abs() <= 1.0e-10);
        assert!(
            (right.recovered_directional_derivative()
                - derivative_scale * left.recovered_directional_derivative())
            .abs()
                <= 1.0e-6
        );
        assert!((right.slack() - derivative_scale * left.slack()).abs() <= 1.0e-6);
        assert!((right.violation() - derivative_scale * left.violation()).abs() <= 1.0e-6);
        assert!((right.tolerance() - derivative_scale * left.tolerance()).abs() <= 1.0e-9);
        assert!((right.loss().unwrap() - left.loss().unwrap()).abs() <= 1.0e-8);
    }
    assert!(
        (transformed.report().total_objective().unwrap()
            - original.report().total_objective().unwrap())
        .abs()
            <= 1.0e-8
    );
}
