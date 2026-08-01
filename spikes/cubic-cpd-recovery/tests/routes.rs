use georbf_cubic_cpd_recovery_spike::run_manufactured_experiment;

#[test]
fn equality_augmented_kkt_recovers_the_physical_canonical_problem() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");
    let equality = &evidence.equality;

    assert_eq!(equality.inertia, equality.expected_inertia);
    assert!(equality.normalized_backward_error <= 1.0e-11);
    assert!(equality.scaling_round_trip_error <= 1.0e-11);
    assert!(equality.recovered.side_condition_violation <= 1.0e-10);
    assert!(equality.recovered.hard_violation <= 1.0e-8);
    assert!(equality.manufactured_truth_error <= 1.0e-8);
}

#[test]
fn reduced_qp_recovers_the_same_physical_canonical_problem() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");
    let qp = &evidence.qp;

    assert!(qp.scaled.primal <= 1.0e-8);
    assert!(qp.scaled.dual <= 1.0e-8);
    assert!(qp.scaled.stationarity <= 1.0e-8);
    assert!(qp.scaled.complementarity <= 1.0e-8);
    assert!(qp.scaled.relative_gap <= 1.0e-8);
    assert!(qp.recovery_round_trip_error <= 1.0e-11);
    assert!(qp.physical_slack_equation_violation <= 1.0e-8);
    assert!(qp.recovered.side_condition_violation <= 1.0e-10);
    assert!(qp.recovered.hard_violation <= 1.0e-8);
    assert!(qp.recovered.slacks.iter().all(|slack| slack.is_finite()));
    assert!(qp.manufactured_truth_error <= 1.0e-8);
}

#[test]
fn reduced_socp_recovers_the_same_physical_canonical_problem() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");
    let socp = &evidence.socp;

    assert!(socp.scaled.primal <= 1.0e-8);
    assert!(socp.scaled.dual <= 1.0e-8);
    assert!(socp.scaled.stationarity <= 1.0e-8);
    assert!(socp.scaled.complementarity <= 1.0e-8);
    assert!(socp.scaled.relative_gap <= 1.0e-8);
    assert!(socp.recovery_round_trip_error <= 1.0e-11);
    assert!(socp.physical_slack_equation_violation <= 1.0e-8);
    assert!(socp.recovered.side_condition_violation <= 1.0e-10);
    assert!(socp.recovered.hard_violation <= 1.0e-8);
    assert_eq!(socp.recovered.slacks.len(), 4);
    assert!(socp.manufactured_truth_error <= 1.0e-8);
}

#[test]
fn all_forms_agree_on_canonical_observables_and_objective() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");

    assert!(evidence.cross_route_observable_error <= 1.0e-8);
    assert!(evidence.equality.recovered.field_energy > 0.0);
    assert_eq!(evidence.equality.recovered.residuals.len(), 11);
    assert_eq!(evidence.equality.recovered.slacks.len(), 4);
    assert_eq!(evidence.qp.recovered.slacks.len(), 4);
    assert_eq!(evidence.socp.recovered.slacks.len(), 4);
}
