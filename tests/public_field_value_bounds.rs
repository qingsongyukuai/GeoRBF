use georbf::diagnostics::{
    BackendRegularizationApplication, ProblemDiagnosis, RepresentationEvidenceStage,
};
use georbf::fit::{BoundActiveState, BoundSide};
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{FieldValueObservation, QuadraticPenalty, QuadraticPenaltyError};
use georbf::problem::BuildError;
use georbf::relation::{
    FieldValueBound, FieldValueBoundError, FieldValueViolationPenalty, LinearViolationPenalty,
    LinearViolationPenaltyError, SharedLevelSetBuilder,
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

fn affine_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    1.0 + x + 2.0 * y - 0.5 * z
}

#[test]
fn checked_bound_constructors_reject_non_finite_and_empty_intervals() {
    let location = point(1.0, 2.0, 3.0);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            FieldValueBound::try_lower(SourceId::new("lower"), location, value),
            Err(FieldValueBoundError::NonFiniteBound)
        );
        assert_eq!(
            FieldValueBound::try_upper(SourceId::new("upper"), location, value),
            Err(FieldValueBoundError::NonFiniteBound)
        );
        assert_eq!(
            FieldValueBound::try_interval(SourceId::new("interval-lower"), location, value, 1.0),
            Err(FieldValueBoundError::NonFiniteBound)
        );
        assert_eq!(
            FieldValueBound::try_interval(SourceId::new("interval-upper"), location, -1.0, value),
            Err(FieldValueBoundError::NonFiniteBound)
        );
    }
    assert_eq!(
        FieldValueBound::try_interval(SourceId::new("empty"), location, 2.0, 1.0),
        Err(FieldValueBoundError::EmptyInterval {
            lower: 2.0,
            upper: 1.0,
        })
    );

    let lower = FieldValueBound::try_lower(SourceId::new("signed-zero"), location, -0.0).unwrap();
    assert_eq!(lower.lower_bound().unwrap().to_bits(), 0.0_f64.to_bits());
    let interval =
        FieldValueBound::try_interval(SourceId::new("closed"), location, -0.0, 0.0).unwrap();
    assert_eq!(interval.lower_bound(), Some(0.0));
    assert_eq!(interval.upper_bound(), Some(0.0));
}

#[test]
fn close_support_positive_mode_is_retained_without_geometric_deduplication() {
    let mut problem = builder();
    let close_supports = [point(0.25, 0.5, 0.75), point(0.255, 0.5, 0.75)];
    for (index, location) in [
        point(0.0, 0.0, 0.0),
        point(1.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0),
        point(0.0, 0.0, 1.0),
        close_supports[0],
        close_supports[1],
    ]
    .into_iter()
    .enumerate()
    {
        problem
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("close-value-{index}")),
                    location,
                    affine_value(location),
                )
                .unwrap(),
            )
            .unwrap();
    }
    problem
        .add(
            FieldValueBound::try_lower(
                SourceId::new("close-support-loose-bound"),
                point(0.0, 0.0, 0.0),
                -100.0,
            )
            .unwrap(),
        )
        .unwrap();

    let success = problem
        .build()
        .unwrap()
        .fit()
        .expect("a reliably positive close-support quotient mode must be retained");
    let analysis = success.report().cubic_analysis().unwrap();
    assert_eq!(analysis.fitting_functional_count(), 6);
    assert_eq!(analysis.quotient_construction().quotient_dimension(), 2);
    let factorization = analysis.quotient_factorization();
    assert_eq!(factorization.retained_modes(), 2);
    assert_eq!(factorization.truncated_modes(), 0);
    assert_eq!(factorization.unregularized_llt_count(), 1);
    assert_eq!(factorization.full_spectrum_analysis_count(), 0);
    assert!(
        factorization
            .pivot_intervals()
            .iter()
            .all(|interval| interval.lower_bound() > 0.0)
    );
    let minimum_pivot = factorization
        .pivot_intervals()
        .iter()
        .map(|interval| interval.lower_bound())
        .fold(f64::INFINITY, f64::min);
    let maximum_pivot = factorization
        .pivot_intervals()
        .iter()
        .map(|interval| interval.upper_bound())
        .fold(0.0_f64, f64::max);
    let pivot_ratio = minimum_pivot / maximum_pivot;
    assert!(
        pivot_ratio < 1.0e-3,
        "the close-support mode should be small but certified positive: ratio={pivot_ratio:e}"
    );
    assert!(!factorization.kernel_ridge_applied());
    assert!(!factorization.gram_jitter_applied());
    assert!(!factorization.mode_truncation_applied());
    let representation = success.report().representation_evidence();
    assert!(!representation.problem_regularization_applied());
    assert!(representation.unregularized_canonical_recovery_verified());
    assert!(
        representation
            .backend_factorization_regularization()
            .iter()
            .any(|attempt| attempt.enabled()),
        "Clarabel's factorization-only stabilization must be disclosed"
    );
    assert!(
        representation
            .backend_factorization_regularization()
            .iter()
            .all(|attempt| {
                attempt.static_application() == BackendRegularizationApplication::Applied
                    && attempt.dynamic_application()
                        == BackendRegularizationApplication::EnabledApplicationNotReported
            })
    );

    for location in close_supports {
        let sample = success.model().evaluate(location).unwrap();
        assert!((sample.value() - affine_value(location)).abs() <= 1.0e-8);
    }
}

#[test]
fn violation_penalties_are_checked_and_soft_interval_sides_stay_independent() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            LinearViolationPenalty::try_new(value),
            Err(LinearViolationPenaltyError::NotFinite)
        );
    }
    for value in [-1.0, -0.0, 0.0] {
        assert_eq!(
            LinearViolationPenalty::try_new(value),
            Err(LinearViolationPenaltyError::NotPositive)
        );
    }
    assert_eq!(
        QuadraticPenalty::try_new(0.0),
        Err(QuadraticPenaltyError::NotPositive)
    );

    let bound = FieldValueBound::try_interval_with_quadratic_penalties(
        SourceId::new("soft-interval"),
        point(0.0, 0.0, 0.0),
        -1.0,
        QuadraticPenalty::try_new(2.0).unwrap(),
        1.0,
        QuadraticPenalty::try_new(3.0).unwrap(),
    )
    .unwrap();
    assert!(bound.is_soft());
    assert_eq!(bound.lower_bound(), Some(-1.0));
    assert_eq!(bound.upper_bound(), Some(1.0));

    let mixed = FieldValueBound::try_interval_with_violation_penalties(
        SourceId::new("mixed-soft-interval"),
        point(0.0, 0.0, 0.0),
        -1.0,
        FieldValueViolationPenalty::Quadratic(QuadraticPenalty::try_new(2.0).unwrap()),
        1.0,
        FieldValueViolationPenalty::Linear(LinearViolationPenalty::try_new(3.0).unwrap()),
    )
    .unwrap();
    assert!(mixed.is_soft());
}

#[test]
fn hard_bounds_enter_the_builder_as_atomic_top_level_relations() {
    let location = point(0.0, 0.0, 0.0);
    let mut builder = builder();
    builder
        .add(FieldValueObservation::try_new(SourceId::new("value"), location, 0.0).unwrap())
        .unwrap();
    builder
        .add(FieldValueBound::try_lower(SourceId::new("bound"), location, -1.0).unwrap())
        .unwrap();
    assert!(
        builder
            .add(FieldValueBound::try_upper(SourceId::new("bound"), location, 1.0).unwrap())
            .is_err()
    );

    let snapshot = builder.build().unwrap();
    assert_eq!(snapshot.observation_count(), 1);
    assert_eq!(snapshot.field_value_bound_count(), 1);
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
    builder
        .add(
            FieldValueBound::try_lower(SourceId::new("hard-lower"), point(2.0, 0.0, 0.0), 2.1)
                .unwrap(),
        )
        .unwrap();
    builder
        .add(
            FieldValueBound::try_upper(SourceId::new("hard-upper"), point(0.0, 2.0, 0.0), 5.0)
                .unwrap(),
        )
        .unwrap();
    builder
        .add(
            FieldValueBound::try_interval(
                SourceId::new("hard-interval"),
                point(0.0, 0.0, 2.0),
                5.0,
                5.9,
            )
            .unwrap(),
        )
        .unwrap();

    let success = builder.build().unwrap().fit().unwrap();
    let backend = success.report().backend_fingerprint().unwrap();
    assert_eq!(backend.enabled_features().collect::<Vec<_>>(), ["serde"]);
    let acceptance = success.report().canonical_acceptance().unwrap();
    assert!(acceptance.accepted());
    assert!(acceptance.backend_standard_form_verified());
    assert!(acceptance.backend_standard_form_residual().unwrap() <= 1.0e-8);
    assert!(acceptance.hard_affine_inequality_violation_max().unwrap() <= 1.0e-8);
    assert!(acceptance.scaling_round_trip_error().unwrap() <= 1.0e-11);
    assert!(acceptance.reduction_round_trip_error().unwrap() <= 1.0e-11);
    assert!(
        acceptance
            .backend_internal_scaling_round_trip_error()
            .unwrap()
            <= 1.0e-11
    );
    let residual = success
        .report()
        .attempts()
        .last()
        .unwrap()
        .convex_residual()
        .unwrap();
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
    let analysis = success.report().cubic_analysis().unwrap();
    assert!(
        analysis
            .reduced_condition_estimate()
            .is_some_and(|value| value >= 1.0)
    );
    for (location, expected) in [
        (point(0.0, 0.0, 0.0), 0.0),
        (point(1.0, 0.0, 0.0), 1.0),
        (point(0.0, 1.0, 0.0), 2.0),
        (point(0.0, 0.0, 1.0), 3.0),
    ] {
        let sample = success.model().evaluate(location).unwrap();
        assert!((sample.value() - expected).abs() <= 1.0e-7);
    }
    let assessments = success.report().field_value_bounds();
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
    assert_eq!(lower.semantic_role().as_str(), "field-value-bound/lower");
    assert!((lower.bound() - 2.1).abs() <= 1.0e-12);
    assert!(
        (lower.recovered_value() - 2.1).abs() <= 5.0e-4,
        "lower recovered {}, slack {}, tolerance {}",
        lower.recovered_value(),
        lower.slack(),
        lower.tolerance()
    );
    assert_eq!(lower.violation(), 0.0);
    assert_eq!(
        lower.active_state(),
        if lower.slack() <= lower.tolerance() {
            BoundActiveState::Active
        } else {
            BoundActiveState::Inactive
        }
    );
    assert!(lower.loss().is_none());

    let upper = side("hard-upper", BoundSide::Upper);
    assert!((upper.recovered_value() - 4.0).abs() <= 5.0e-4);
    assert!((upper.slack() - 1.0).abs() <= 5.0e-4);
    assert_eq!(upper.active_state(), BoundActiveState::Inactive);

    let interval_lower = side("hard-interval", BoundSide::Lower);
    assert!((interval_lower.slack() - 0.9).abs() <= 5.0e-4);
    assert_eq!(interval_lower.active_state(), BoundActiveState::Inactive);
    let interval_upper = side("hard-interval", BoundSide::Upper);
    assert!(interval_upper.violation() <= 1.0e-4);
}

#[test]
fn soft_bounds_require_normalization_and_recover_independent_physical_losses() {
    let mut missing_normalization = builder();
    missing_normalization
        .add(
            FieldValueBound::try_lower_with_quadratic_penalty(
                SourceId::new("soft"),
                point(0.0, 0.0, 0.0),
                1.0,
                QuadraticPenalty::try_new(2.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    let failure = missing_normalization.build().unwrap_err();
    assert!(
        failure
            .errors()
            .contains(&BuildError::MissingFieldEnergyNormalization)
    );

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
    for bound in [
        FieldValueBound::try_lower_with_quadratic_penalty(
            SourceId::new("soft-lower"),
            point(0.0, 0.0, 0.0),
            1.0,
            QuadraticPenalty::try_new(2.0).unwrap(),
        )
        .unwrap(),
        FieldValueBound::try_lower_with_quadratic_penalty(
            SourceId::new("soft-lower-duplicate"),
            point(0.0, 0.0, 0.0),
            1.0,
            QuadraticPenalty::try_new(2.0).unwrap(),
        )
        .unwrap(),
        FieldValueBound::try_upper_with_linear_violation_penalty(
            SourceId::new("soft-upper"),
            point(0.0, 0.0, 0.0),
            -2.0,
            LinearViolationPenalty::try_new(3.0).unwrap(),
        )
        .unwrap(),
        FieldValueBound::try_interval_with_quadratic_penalties(
            SourceId::new("soft-interval"),
            point(0.0, 0.0, 0.0),
            1.0,
            QuadraticPenalty::try_new(4.0).unwrap(),
            2.0,
            QuadraticPenalty::try_new(5.0).unwrap(),
        )
        .unwrap(),
    ] {
        builder.add(bound).unwrap();
    }

    let success = builder.build().unwrap().fit().unwrap();
    let report = success.report();
    assert_eq!(report.field_value_bounds().len(), 5);
    assert_eq!(report.problem_size().auxiliary_variables(), 5);
    assert_eq!(report.problem_size().quadratic_objective_terms(), 4);
    assert_eq!(report.problem_size().linear_objective_terms(), 1);
    assert_eq!(report.problem_size().affine_inequality_constraints(), 10);
    let ledger = report.all_source_recovery().unwrap();
    assert!(ledger.verified());
    assert_eq!(ledger.canonical_hard_relation_count(), 4);
    assert_eq!(ledger.canonical_soft_relation_count(), 5);
    assert_eq!(ledger.recovery_edge_count(), 9);
    assert_eq!(ledger.solver_relation_row_count(), 9);
    assert_eq!(ledger.recovered_sources(), ledger.participating_sources());

    let side = |source: &str, expected_side: BoundSide| {
        report
            .field_value_bounds()
            .iter()
            .find(|assessment| {
                assessment.source_id().as_str() == source && assessment.side() == expected_side
            })
            .unwrap()
    };
    for source in ["soft-lower", "soft-lower-duplicate"] {
        let assessment = side(source, BoundSide::Lower);
        assert!((assessment.violation() - 1.0).abs() <= 1.0e-6);
        assert!((assessment.loss().unwrap() - 1.0).abs() <= 1.0e-5);
        assert_eq!(assessment.quadratic_penalty().unwrap().weight(), 2.0);
        assert!(assessment.linear_violation_penalty().is_none());
        assert_eq!(assessment.active_state(), BoundActiveState::Active);
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
fn locally_provable_hard_bound_conflicts_stop_before_the_backend() {
    let location = point(0.0, 0.0, 0.0);
    let mut equality_conflict = builder();
    equality_conflict
        .add(FieldValueObservation::try_new(SourceId::new("equality"), location, 0.0).unwrap())
        .unwrap();
    equality_conflict
        .add(FieldValueBound::try_lower(SourceId::new("lower"), location, 1.0).unwrap())
        .unwrap();
    let failure = equality_conflict.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
    let evidence = failure.report().direct_input_conflict().unwrap();
    assert_eq!(evidence.first_source().as_str(), "equality");
    assert_eq!(evidence.second_source().as_str(), "lower");
    assert_eq!(
        failure
            .report()
            .conflict_witness()
            .unwrap()
            .relations()
            .iter()
            .map(|relation| (
                relation.source().source_id().as_str(),
                relation.source().semantic_role().as_str(),
            ))
            .collect::<Vec<_>>(),
        [
            ("equality", "field-value-observation/value"),
            ("lower", "field-value-bound/lower"),
        ]
    );

    let mut empty_feasible_set = builder();
    empty_feasible_set
        .add(FieldValueBound::try_lower(SourceId::new("lower"), location, 2.0).unwrap())
        .unwrap();
    empty_feasible_set
        .add(FieldValueBound::try_upper(SourceId::new("upper"), location, 1.0).unwrap())
        .unwrap();
    let failure = empty_feasible_set.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().attempts().is_empty());
    assert_eq!(failure.report().direct_input_conflicts().len(), 1);
}

#[test]
fn duplicate_hard_bounds_share_one_backend_relation_and_keep_both_sources() {
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
    for source in ["duplicate-a", "duplicate-b"] {
        builder
            .add(
                FieldValueBound::try_upper(SourceId::new(source), point(2.0, 0.0, 0.0), 3.0)
                    .unwrap(),
            )
            .unwrap();
    }

    let success = builder.build().unwrap().fit().unwrap();
    let ledger = success.report().all_source_recovery().unwrap();
    assert!(ledger.verified());
    assert_eq!(ledger.canonical_hard_relation_count(), 5);
    assert_eq!(ledger.recovery_edge_count(), 6);
    assert_eq!(ledger.recovered_sources(), ledger.participating_sources());
    assert_eq!(success.report().field_value_bounds().len(), 2);
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
            .field_value_bounds()
            .iter()
            .map(|assessment| assessment.source_id().as_str())
            .collect::<Vec<_>>(),
        ["duplicate-a", "duplicate-b"]
    );
}

#[test]
fn general_bound_infeasibility_requires_a_validated_farkas_certificate() {
    let mut level = SharedLevelSetBuilder::new(GroupId::new("one-level"));
    for (source, location) in [
        ("level/origin", point(0.0, 0.0, 0.0)),
        ("level/east", point(1.0, 0.0, 0.0)),
        ("level/north", point(0.0, 1.0, 0.0)),
        ("level/up", point(0.0, 0.0, 1.0)),
    ] {
        level.add_member(SourceId::new(source), location).unwrap();
    }
    level
        .add_member(SourceId::new("level/origin-alias"), point(0.0, 0.0, 0.0))
        .unwrap();
    let mut builder = builder();
    builder.add(level.build().unwrap()).unwrap();
    builder
        .add(FieldValueBound::try_lower(SourceId::new("lower"), point(0.0, 0.0, 0.0), 2.0).unwrap())
        .unwrap();
    builder
        .add(FieldValueBound::try_upper(SourceId::new("upper"), point(1.0, 0.0, 0.0), 1.0).unwrap())
        .unwrap();

    let failure = builder.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::InfeasibleProblem);
    let representation = failure.report().representation_evidence();
    assert_eq!(
        representation.failure_stage(),
        Some(RepresentationEvidenceStage::Backend)
    );
    assert!(representation.quotient_construction().is_some());
    assert!(representation.quotient_factorization().is_some());
    assert_eq!(representation.source_count(), 7);
    assert_eq!(representation.recovery_source_count(), 0);
    assert!(!representation.problem_regularization_applied());
    assert!(
        representation
            .backend_factorization_regularization()
            .iter()
            .any(|attempt| attempt.enabled())
    );
    assert!(!failure.report().attempts().is_empty());
    assert!(failure.report().canonical_acceptance().is_none());
    assert!(failure.report().hard_relations().is_empty());
    assert!(failure.report().field_value_bounds().is_empty());
    let certificate = failure
        .report()
        .infeasibility_certificate()
        .expect("infeasible diagnoses retain independently validated evidence");
    assert!(certificate.backend_invoked());
    assert!(certificate.finite());
    assert!((certificate.normalized_ray_norm() - 1.0).abs() <= 1.0e-12);
    assert!(certificate.stationarity_residual() <= certificate.residual_limit());
    assert!(certificate.dual_cone_violation() <= certificate.residual_limit());
    assert!(certificate.separation_margin() >= certificate.separation_limit());
    assert!(certificate.provenance_verified());
    assert!(certificate.recovery_round_trip_error() <= 1.0e-11);
    assert!(
        certificate
            .sources()
            .iter()
            .any(|source| source.source_id().as_str() == "lower")
    );
    assert!(
        certificate
            .sources()
            .iter()
            .any(|source| source.source_id().as_str() == "upper")
    );
    assert!(
        failure
            .report()
            .attempts()
            .iter()
            .any(|attempt| attempt.certificate_present())
    );

    let witness = failure
        .report()
        .conflict_witness()
        .expect("validated infeasibility is localized to original canonical sources");
    assert!(witness.backend_invoked());
    assert!(witness.provenance_verified());
    assert!(witness.canonical_residual() <= witness.residual_limit());
    assert!(witness.separation_margin() >= witness.separation_limit());
    let witness_sources = witness
        .sources()
        .iter()
        .map(|source| source.source_id().as_str())
        .collect::<Vec<_>>();
    assert!(witness_sources.contains(&"lower"));
    assert!(witness_sources.contains(&"upper"));
    assert!(witness_sources.contains(&"level/origin"));
    assert!(witness_sources.contains(&"level/origin-alias"));
    assert!(witness_sources.windows(2).all(|pair| pair[0] <= pair[1]));
    assert!(
        witness
            .relations()
            .iter()
            .all(|relation| relation.multiplier().is_finite() && relation.multiplier() != 0.0)
    );
}

fn covariant_bound_fit(
    length_scale: f64,
    field_scale: f64,
    reverse_inputs: bool,
) -> georbf::fit::FitSuccess {
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
    let mut hard = [
        ("origin", [0.0, 0.0, 0.0], 0.0),
        ("east", [1.0, 0.0, 0.0], 1.0),
        ("north", [0.0, 1.0, 0.0], 2.0),
        ("up", [0.0, 0.0, 1.0], 3.0),
    ];
    if reverse_inputs {
        hard.reverse();
    }
    for (source, support, value) in hard {
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
    let quadratic = FieldValueBound::try_lower_with_quadratic_penalty(
        SourceId::new("quadratic"),
        point(x, y, z),
        field_scale * 1.8,
        QuadraticPenalty::try_new(2.0 / field_scale.powi(2)).unwrap(),
    )
    .unwrap();
    let linear = FieldValueBound::try_upper_with_linear_violation_penalty(
        SourceId::new("linear"),
        point(x, y, z),
        field_scale * -0.2,
        LinearViolationPenalty::try_new(3.0 / field_scale).unwrap(),
    )
    .unwrap();
    if reverse_inputs {
        builder.add(linear).unwrap();
        builder.add(quadratic).unwrap();
    } else {
        builder.add(quadratic).unwrap();
        builder.add(linear).unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

#[test]
fn bound_results_are_covariant_under_frame_and_field_unit_changes() {
    let original = covariant_bound_fit(1.0, 1.0, false);
    let transformed = covariant_bound_fit(2.5, 4.0, false);
    let original_sample = original.model().evaluate(point(0.35, -0.2, 0.1)).unwrap();
    let transformed_sample = transformed
        .model()
        .evaluate(point(9.5, -2.125, 4.25))
        .unwrap();
    assert!((transformed_sample.value() - 4.0 * original_sample.value()).abs() <= 1.0e-7);
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
        assert!((actual - expected).abs() <= 1.0e-7);
    }
    for (original_bound, transformed_bound) in original
        .report()
        .field_value_bounds()
        .iter()
        .zip(transformed.report().field_value_bounds())
    {
        assert_eq!(transformed_bound.source_id(), original_bound.source_id());
        assert!((transformed_bound.bound() - 4.0 * original_bound.bound()).abs() <= 1.0e-7);
        assert!(
            (transformed_bound.recovered_value() - 4.0 * original_bound.recovered_value()).abs()
                <= 1.0e-7
        );
        assert!((transformed_bound.violation() - 4.0 * original_bound.violation()).abs() <= 1.0e-7);
        assert!(
            (transformed_bound.loss().unwrap() - original_bound.loss().unwrap()).abs() <= 1.0e-7
        );
    }
    assert!(
        (transformed.report().total_objective().unwrap()
            - original.report().total_objective().unwrap())
        .abs()
            <= 1.0e-7
    );
}

#[test]
fn input_permutation_preserves_bound_results_and_stable_report_order() {
    let baseline = covariant_bound_fit(1.0, 1.0, false);
    let permuted = covariant_bound_fit(1.0, 1.0, true);
    assert_eq!(
        baseline
            .report()
            .field_value_bounds()
            .iter()
            .map(|bound| (bound.source_id().clone(), bound.side()))
            .collect::<Vec<_>>(),
        permuted
            .report()
            .field_value_bounds()
            .iter()
            .map(|bound| (bound.source_id().clone(), bound.side()))
            .collect::<Vec<_>>()
    );
    for query in [
        point(0.2, -0.3, 0.4),
        point(-0.1, 0.25, 0.5),
        point(0.5, 0.5, 0.5),
    ] {
        let expected = baseline.model().evaluate(query).unwrap();
        let actual = permuted.model().evaluate(query).unwrap();
        assert!((actual.value() - expected.value()).abs() <= 1.0e-10);
    }
    assert!(
        (permuted.report().total_objective().unwrap()
            - baseline.report().total_objective().unwrap())
        .abs()
            <= 1.0e-10
    );
}
