use georbf::diagnostics::{SolveAttemptKind, SolveAttemptTermination};
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::KernelKind;
use georbf::observation::{FieldValueObservation, GradientObservation};
use georbf::problem::{FitConfiguration, ThreadBudget};
use georbf::{Point3, ProblemBuilder, SourceId, Vector3};
use std::num::NonZeroUsize;

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
    builder.set_fit_configuration(
        FitConfiguration::default()
            .with_thread_budget(ThreadBudget::Exact(NonZeroUsize::new(1).unwrap())),
    );

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
    let aliased_location = point(-1.0, -1.0, -1.0);
    builder
        .add(
            FieldValueObservation::try_new(
                SourceId::new("value-alias"),
                aliased_location,
                affine_value(aliased_location),
            )
            .expect("the duplicate exact fact is finite"),
        )
        .expect("a distinct source may report the same exact fact");

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
    assert_eq!(
        report.requested_thread_budget(),
        ThreadBudget::Exact(NonZeroUsize::new(1).unwrap())
    );
    let backend = report
        .backend_fingerprint()
        .expect("a successful fit records its selected backend");
    let legacy_features: [&str; 2] = backend.features();
    assert_eq!(legacy_features, ["linalg", "std"]);
    assert_eq!(
        backend.enabled_features().collect::<Vec<_>>(),
        ["linalg", "std"]
    );
    assert_eq!(backend.requested_threads(), 1);
    assert_eq!(backend.actual_threads(), 1);
    assert!(report.attempts().iter().all(|attempt| {
        attempt.backend_fingerprint().requested_threads() == 1
            && attempt.backend_fingerprint().actual_threads() == 1
    }));
    let problem_size = report.problem_size();
    assert_eq!(problem_size.input_observations(), 9);
    assert_eq!(problem_size.scalar_hard_relations(), 15);
    assert_eq!(problem_size.canonical_hard_equalities(), Some(14));
    assert_eq!(problem_size.center_coefficients(), Some(14));
    assert_eq!(problem_size.semantic_latents(), 0);
    assert_eq!(problem_size.auxiliary_variables(), 0);
    assert_eq!(problem_size.cone_blocks(), 0);
    assert_eq!(problem_size.primal_variables(), Some(18));
    assert_eq!(problem_size.equality_constraints(), Some(18));
    assert_eq!(problem_size.kkt_dimension(), Some(36));
    let accepted_attempt = report
        .attempts()
        .last()
        .expect("a successful fit retains its accepted backend attempt");
    assert_eq!(
        accepted_attempt.kind(),
        SolveAttemptKind::BunchKaufmanRefinement
    );
    assert_eq!(
        accepted_attempt.termination(),
        SolveAttemptTermination::CandidateProduced
    );
    assert_eq!(accepted_attempt.settings().kind(), accepted_attempt.kind());
    assert_eq!(
        accepted_attempt.scaling().method(),
        "block-aware Ruiz max-norm diagonal congruence"
    );
    assert_eq!(accepted_attempt.scaling().rounds(), 8);
    assert!(accepted_attempt.residual().is_some());
    assert_eq!(accepted_attempt.failure_reason(), None);
    assert_close(report.field_energy().expect("fit succeeded"), 0.0);
    assert_eq!(report.hard_relations().len(), 15);
    assert!(report.hard_relations().iter().all(|relation| {
        relation.residual().abs() <= relation.tolerance()
            && !relation.source_id().as_str().is_empty()
            && !relation.semantic_role().as_str().is_empty()
    }));
    assert_eq!(
        report
            .hard_relations()
            .iter()
            .filter(|relation| relation.source_id().as_str() == "gradient-0")
            .map(|relation| relation.semantic_role().as_str())
            .collect::<Vec<_>>(),
        [
            "gradient-observation/component/0",
            "gradient-observation/component/1",
            "gradient-observation/component/2",
        ]
    );
    assert_eq!(
        report
            .hard_relations()
            .iter()
            .filter(|relation| {
                matches!(relation.source_id().as_str(), "value-0" | "value-alias")
            })
            .count(),
        2
    );

    let cubic_analysis = report
        .cubic_analysis()
        .expect("the Cubic path retains representation rank and condition evidence");
    assert_eq!(cubic_analysis.fitting_functional_count(), 14);
    assert_eq!(cubic_analysis.polynomial_dimension(), 4);
    assert_eq!(cubic_analysis.polynomial_rank(), 4);
    let quotient = cubic_analysis.quotient_construction();
    assert_eq!(quotient.quotient_dimension(), 10);
    assert_eq!(quotient.householder_reflector_count(), 4);
    assert_eq!(quotient.congruence_pass_count(), 2);
    assert!(quotient.householder_orthogonality_error() <= 1.0e-11);
    assert!(cubic_analysis.null_space_defect() <= 1.0e-11);
    assert!(quotient.canonical_response_round_trip_error() <= 1.0e-11);
    let factorization = cubic_analysis.quotient_factorization();
    assert_eq!(factorization.quotient_dimension(), 10);
    assert_eq!(factorization.retained_modes(), 10);
    assert_eq!(factorization.truncated_modes(), 0);
    assert_eq!(factorization.unregularized_llt_count(), 1);
    assert_eq!(factorization.full_spectrum_analysis_count(), 0);
    assert!(factorization.normalized_backward_error() <= 1.0e-11);
    assert_eq!(factorization.pivot_intervals().len(), 10);
    assert!(
        factorization
            .pivot_intervals()
            .iter()
            .all(|interval| interval.lower_bound() > 0.0)
    );
    assert!(factorization.field_energy_identity_error().unwrap() <= 1.0e-11);
    assert!(factorization.side_condition_error().unwrap() <= 1.0e-11);
    assert!(factorization.recovery_round_trip_error().unwrap() <= 1.0e-11);
    assert!(factorization.canonical_response_round_trip_error().unwrap() <= 1.0e-11);
    assert!(!factorization.kernel_ridge_applied());
    assert!(!factorization.gram_jitter_applied());
    assert!(!factorization.mode_truncation_applied());
    assert!(cubic_analysis.reduced_smallest_singular_value() > 0.0);
    let backend_rank = report
        .backend_rank()
        .expect("the admitted KKT path retains rank evidence");
    assert!(backend_rank.is_full_rank());
    assert!(backend_rank.condition_estimate().is_some());
    let inertia = report
        .inertia()
        .expect("the admitted KKT path retains inertia evidence");
    assert_eq!(inertia.expected(), inertia.observed());
    let acceptance = report
        .canonical_acceptance()
        .expect("success retains physical Recover-and-Verify evidence");
    assert!(acceptance.accepted());
    assert!(acceptance.recovery_finite());
    assert!(acceptance.provenance_verified());
    assert!(acceptance.side_condition().is_some());

    fn assert_send_sync_clone<T: Send + Sync + Clone>() {}
    assert_send_sync_clone::<georbf::SolvedModel>();
    let model_debug = format!("{:?}", success.model());
    assert_eq!(model_debug, "SolvedModel { .. }");
    assert!(!model_debug.contains("coefficient"));
    assert!(!model_debug.contains("polynomial"));
}
