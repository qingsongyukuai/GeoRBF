use georbf_cubic_cpd_recovery_spike::{FailureKind, run_counterexamples};

#[test]
fn rank_deficient_polynomial_is_rejected_before_solving() {
    let evidence = run_counterexamples().expect("the counterexamples should run");

    assert_eq!(
        evidence.rank_deficient_polynomial.kind,
        FailureKind::PolynomialRankDeficient
    );
    assert!(!evidence.rank_deficient_polynomial.solver_invoked);
    assert!(
        !evidence
            .rank_deficient_polynomial
            .hidden_regularization_applied
    );
}

#[test]
fn nonpositive_reduced_pairing_is_rejected_before_solving() {
    let evidence = run_counterexamples().expect("the counterexamples should run");

    assert_eq!(
        evidence.nonpositive_reduced_pairing.kind,
        FailureKind::ReducedPairingNotPositive
    );
    assert!(!evidence.nonpositive_reduced_pairing.solver_invoked);
    assert!(
        !evidence
            .nonpositive_reduced_pairing
            .hidden_regularization_applied
    );
}

#[test]
fn broken_recovery_is_rejected_after_backend_contract_passes() {
    let evidence = run_counterexamples().expect("the counterexamples should run");
    let positive = georbf_cubic_cpd_recovery_spike::run_manufactured_experiment()
        .expect("policy is available");

    assert_eq!(
        evidence.broken_recovery.kind,
        FailureKind::RecoveryVerification
    );
    assert!(evidence.broken_recovery.solver_invoked);
    assert!(evidence.broken_recovery.backend_contract_passed);
    assert!(!evidence.broken_recovery.hidden_regularization_applied);
    assert!(evidence.broken_recovery.detected_violation > positive.acceptance.canonical);
}
