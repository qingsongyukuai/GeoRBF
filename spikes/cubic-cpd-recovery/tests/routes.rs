use georbf_cubic_cpd_recovery_spike::{NUMERICAL_POLICY_VERSION, run_manufactured_experiment};

#[test]
fn equality_augmented_kkt_recovers_the_physical_canonical_problem() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");
    let equality = &evidence.equality;
    let limits = evidence.acceptance;

    assert_eq!(evidence.numerical_policy_version, NUMERICAL_POLICY_VERSION);
    assert_eq!(equality.inertia, equality.expected_inertia);
    assert!(equality.normalized_backward_error <= limits.backward_error);
    assert!(equality.scaling_round_trip_error <= limits.round_trip);
    assert!(equality.recovered.side_condition_violation <= limits.side_condition);
    assert!(equality.recovered.hard_violation <= limits.canonical);
    assert!(equality.manufactured_truth_error <= limits.canonical);
}

#[test]
fn reduced_qp_recovers_the_same_physical_canonical_problem() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");
    let qp = &evidence.qp;
    let limits = evidence.acceptance;

    assert!(qp.scaled.primal <= limits.backend_residual);
    assert!(qp.scaled.dual <= limits.backend_residual);
    assert!(qp.scaled.stationarity <= limits.backend_residual);
    assert!(qp.scaled.complementarity <= limits.backend_residual);
    assert!(qp.scaled.relative_gap <= limits.backend_residual);
    assert!(qp.recovery_round_trip_error <= limits.round_trip);
    assert!(qp.physical_slack_equation_violation <= limits.canonical);
    assert!(qp.recovered.side_condition_violation <= limits.side_condition);
    assert!(qp.recovered.hard_violation <= limits.canonical);
    assert!(qp.cumulative_scaling_bounds.minimum >= limits.cumulative_scaling_minimum);
    assert!(qp.cumulative_scaling_bounds.maximum <= limits.cumulative_scaling_maximum);
    assert!(qp.recovered.slacks.iter().all(|slack| slack.is_finite()));
    assert!(qp.manufactured_truth_error <= limits.canonical);
}

#[test]
fn reduced_socp_recovers_the_same_physical_canonical_problem() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");
    let socp = &evidence.socp;
    let limits = evidence.acceptance;

    assert!(socp.scaled.primal <= limits.backend_residual);
    assert!(socp.scaled.dual <= limits.backend_residual);
    assert!(socp.scaled.stationarity <= limits.backend_residual);
    assert!(socp.scaled.complementarity <= limits.backend_residual);
    assert!(socp.scaled.relative_gap <= limits.backend_residual);
    assert!(socp.recovery_round_trip_error <= limits.round_trip);
    assert!(socp.physical_slack_equation_violation <= limits.canonical);
    assert!(socp.recovered.side_condition_violation <= limits.side_condition);
    assert!(socp.recovered.hard_violation <= limits.canonical);
    assert!(socp.cumulative_scaling_bounds.minimum >= limits.cumulative_scaling_minimum);
    assert!(socp.cumulative_scaling_bounds.maximum <= limits.cumulative_scaling_maximum);
    assert_eq!(socp.recovered.slacks.len(), 4);
    assert!(socp.manufactured_truth_error <= limits.canonical);
}

#[test]
fn all_forms_agree_on_canonical_observables_and_objective() {
    let evidence = run_manufactured_experiment().expect("the manufactured experiment should run");

    assert!(evidence.cross_route_observable_error <= evidence.acceptance.cross_route);
    assert!(evidence.equality.recovered.field_energy > 0.0);
    assert_eq!(evidence.equality.recovered.residuals.len(), 11);
    assert_eq!(evidence.equality.recovered.slacks.len(), 4);
    assert_eq!(evidence.qp.recovered.slacks.len(), 4);
    assert_eq!(evidence.socp.recovered.slacks.len(), 4);
}
