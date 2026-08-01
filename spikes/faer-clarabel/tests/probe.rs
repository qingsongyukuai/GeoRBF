use georbf_backend_probe::{
    CLARABEL_VERSION, CapabilityState, CapacityOutcome, FAER_VERSION, SolveAttemptTermination,
    probe_clarabel_primal_infeasible, probe_clarabel_qp, probe_clarabel_socp,
    probe_clarabel_unbounded, probe_faer,
};

const SOLVE_TOLERANCE: f64 = 1.0e-8;

#[test]
fn faer_exposes_the_required_dense_kkt_and_rank_evidence() {
    let evidence = probe_faer().expect("the faer probe should run");

    assert_eq!(evidence.version, FAER_VERSION);
    assert_eq!(evidence.kkt.solution.len(), 3);
    assert!((evidence.kkt.solution[0] - 0.5).abs() <= SOLVE_TOLERANCE);
    assert!((evidence.kkt.solution[1] + 0.5).abs() <= SOLVE_TOLERANCE);
    assert!((evidence.kkt.solution[2] - 1.0).abs() <= SOLVE_TOLERANCE);
    assert!(evidence.kkt.normalized_backward_error <= 1.0e-11);
    assert_eq!(evidence.kkt.inertia.positive, 2);
    assert_eq!(evidence.kkt.inertia.negative, 1);
    assert_eq!(evidence.kkt.inertia.zero, 0);
    assert!(!evidence.kkt.adapter_scaling_applied);

    assert!(evidence.factorizations.cholesky_spd_succeeded);
    assert!(evidence.factorizations.cholesky_indefinite_rejected);
    assert_eq!(evidence.factorizations.col_pivoted_qr_rank, 1);
    assert_eq!(evidence.factorizations.svd_rank, 1);
    assert_eq!(evidence.factorizations.singular_values.len(), 2);
    assert!(evidence.factorizations.singular_values[0] > 0.0);
    assert!(evidence.factorizations.singular_values[1] <= 1.0e-12);
    assert_eq!(
        evidence.factorizations.capacity.representable_oversize,
        CapacityOutcome::RejectedBeforeAllocation
    );
    assert_eq!(
        evidence.factorizations.capacity.arithmetic_overflow,
        CapacityOutcome::RejectedBeforeAllocation
    );
    assert_eq!(
        evidence.factorizations.capacity.state,
        CapabilityState::Ambiguous
    );
    assert!(evidence.factorizations.capacity.first_oversize_square_bytes > 8 * 1024 * 1024 * 1024);

    assert_eq!(evidence.threads.requested, 1);
    assert_eq!(evidence.threads.actual, 1);
    assert!(!evidence.threads.process_global_state_modified);
    assert!(evidence.failure_reason.is_none());
}

#[test]
fn clarabel_qp_exposes_candidates_duals_settings_and_residuals() {
    let evidence = probe_clarabel_qp().expect("the Clarabel QP probe should run");

    assert_eq!(evidence.version, CLARABEL_VERSION);
    assert_eq!(evidence.problem_class, "convex_qp");
    assert_eq!(evidence.termination, SolveAttemptTermination::Solved);
    assert!((evidence.primal[0] - 3.0 / 7.0).abs() <= SOLVE_TOLERANCE);
    assert!((evidence.primal[1] - 3.0 / 14.0).abs() <= SOLVE_TOLERANCE);
    assert_eq!(evidence.dual.len(), 5);
    assert_eq!(evidence.slack.len(), 5);
    assert!(evidence.primal_residual <= 1.0e-8);
    assert!(evidence.dual_residual <= 1.0e-8);
    assert!(evidence.absolute_gap <= 1.0e-8);
    assert!(evidence.settings.equilibration);
    assert!(evidence.settings.iterative_refinement);
    assert!(evidence.settings.static_regularization);
    assert!(evidence.settings.dynamic_regularization);
    assert_eq!(evidence.settings.max_threads, 1);
    assert_eq!(evidence.threads.actual, 1);
    assert!(!evidence.threads.process_global_state_modified);
    assert!(evidence.failure_reason.is_none());
    assert_valid_scaling(&evidence);
}

#[test]
fn clarabel_socp_exposes_a_soc_candidate_and_evidence() {
    let evidence = probe_clarabel_socp().expect("the Clarabel SOCP probe should run");

    assert_eq!(evidence.problem_class, "convex_socp");
    assert_eq!(evidence.termination, SolveAttemptTermination::Solved);
    assert!((evidence.primal[0] - 1.0).abs() <= SOLVE_TOLERANCE);
    assert!((evidence.primal[1] - 1.0).abs() <= SOLVE_TOLERANCE);
    assert!(evidence.primal_residual <= 1.0e-8);
    assert!(evidence.dual_residual <= 1.0e-8);
    assert!(evidence.absolute_gap <= 1.0e-8);
    assert_eq!(evidence.cones, vec!["second_order(3)"]);
    assert_valid_scaling(&evidence);
}

#[test]
fn clarabel_returns_independently_checkable_infeasibility_evidence() {
    let evidence = probe_clarabel_primal_infeasible().expect("the infeasible probe should run");

    assert_eq!(
        evidence.attempt.termination,
        SolveAttemptTermination::PrimalInfeasible
    );
    assert_eq!(evidence.attempt.version, CLARABEL_VERSION);
    assert_eq!(evidence.attempt.linear_solver, "qdldl");
    assert_eq!(evidence.attempt.threads.requested, 1);
    assert_eq!(evidence.attempt.threads.actual, 1);
    assert!(evidence.attempt.failure_reason.is_none());
    assert_valid_scaling(&evidence.attempt);
    assert!(evidence.attempt.primal_infeasibility_residual <= 1.0e-8);
    assert!(evidence.certificate_residual <= 1.0e-8);
    assert!(evidence.cone_violation <= 1.0e-12);
    assert!(evidence.separation_margin >= 1.0e-7);
    assert!(evidence.certificate.iter().all(|value| value.is_finite()));
}

#[test]
fn clarabel_returns_independently_checkable_unboundedness_evidence() {
    let evidence = probe_clarabel_unbounded().expect("the unbounded probe should run");

    assert_eq!(
        evidence.attempt.termination,
        SolveAttemptTermination::DualInfeasible
    );
    assert_eq!(evidence.attempt.version, CLARABEL_VERSION);
    assert_eq!(evidence.attempt.linear_solver, "qdldl");
    assert_eq!(evidence.attempt.threads.requested, 1);
    assert_eq!(evidence.attempt.threads.actual, 1);
    assert!(evidence.attempt.failure_reason.is_none());
    assert_valid_scaling(&evidence.attempt);
    assert!(evidence.attempt.dual_infeasibility_residual <= 1.0e-8);
    assert!(evidence.certificate_residual <= 1.0e-8);
    assert!(evidence.cone_violation <= 1.0e-12);
    assert!(evidence.descent_margin >= 1.0e-7);
    assert!(evidence.certificate.iter().all(|value| value.is_finite()));
}

fn assert_valid_scaling(evidence: &georbf_backend_probe::ClarabelEvidence) {
    assert_eq!(evidence.scaling.variable.len(), evidence.primal.len());
    assert_eq!(
        evidence.scaling.inverse_variable.len(),
        evidence.primal.len()
    );
    assert_eq!(evidence.scaling.constraint.len(), evidence.slack.len());
    assert_eq!(
        evidence.scaling.inverse_constraint.len(),
        evidence.slack.len()
    );
    assert!(evidence.scaling.objective.is_finite());
    assert!(evidence.scaling.objective > 0.0);
    for (forward, inverse) in evidence
        .scaling
        .variable
        .iter()
        .zip(&evidence.scaling.inverse_variable)
        .chain(
            evidence
                .scaling
                .constraint
                .iter()
                .zip(&evidence.scaling.inverse_constraint),
        )
    {
        assert!(forward.is_finite());
        assert!(inverse.is_finite());
        assert!(*forward > 0.0);
        assert!(*inverse > 0.0);
        assert!((forward * inverse - 1.0).abs() <= 1.0e-12);
    }
}
