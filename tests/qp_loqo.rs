use georbf::{
    loqo_step_divisor, solve_loqo_qp, solve_loqo_qp_with_options, validate_loqo_qp, DenseMatrix,
    DenseVector, Error, LoqoKktStage, LoqoOptions, LoqoSolveErrorKind, LoqoStepError,
    LoqoStopReason,
};

fn matrix(rows: usize, columns: usize, values: &[f64]) -> DenseMatrix {
    DenseMatrix::from_row_major(rows, columns, values.to_vec()).expect("valid test matrix")
}

fn vector(values: &[f64]) -> DenseVector {
    DenseVector::from_values(values.to_vec())
}

fn assert_close(left: f64, right: f64, tolerance: f64) {
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    assert!(
        (left - right).abs() <= tolerance * scale,
        "{left:.17e} != {right:.17e}"
    );
}

#[test]
fn inactive_single_box_keeps_the_unconstrained_minimum() {
    let solution = solve_loqo_qp(
        &matrix(1, 1, &[1.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[-1.0]),
        &vector(&[2.0]),
    )
    .unwrap();

    assert_eq!(
        solution.weights().values()[0].to_bits(),
        0x3b4d_1b63_2e9f_17e0
    );
    assert_eq!(solution.objective().to_bits(), 0x36aa_79bb_3512_707b);
    assert_close(solution.weights().values()[0], 0.0, 1.0e-7);
    assert_close(solution.objective(), 0.0, 1.0e-10);
    assert!(solution.residual().accepted());
    assert_eq!(solution.stop_reason(), LoqoStopReason::SignificantFigures);
    assert_eq!(solution.trace().len(), 10);
    assert_eq!(solution.trace()[0].iteration(), 1);
    assert_eq!(solution.trace()[0].primal_objective(), 0.0625);
    assert_eq!(solution.trace()[0].dual_objective(), -199.3125);
    assert_eq!(
        solution.trace()[0].predictor_mu(),
        Some(f64::from_bits(0x3fcd_0662_af4f_8731))
    );
    assert_eq!(solution.trace()[9].iteration(), 10);
    assert_close(
        solution.trace()[9].significant_figures(),
        6.995_317_254_619_095,
        2.0e-14,
    );
    assert_eq!(solution.kkt_solves()[0].dimension(), 2);
    assert_eq!(solution.kkt_solves()[0].stage(), LoqoKktStage::Initial);
    assert_eq!(solution.kkt_solves().len(), 19);
    assert_eq!(solution.kkt_solves()[1].stage(), LoqoKktStage::Predictor(1));
    assert_eq!(solution.kkt_solves()[2].stage(), LoqoKktStage::Corrector(1));
}

#[test]
fn active_single_and_two_sided_boxes_match_the_analytic_minima() {
    let active = solve_loqo_qp(
        &matrix(1, 1, &[1.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[1.0]),
        &vector(&[1.0]),
    )
    .unwrap();
    assert_eq!(
        active.weights().values()[0].to_bits(),
        0x3ff0_0000_1b67_1ccd
    );
    assert_eq!(active.objective().to_bits(), 0x3ff0_0000_36ce_39c9);
    assert_eq!(active.trace().len(), 10);
    assert_eq!(
        active.trace()[1].predictor_mu(),
        Some(f64::from_bits(0x3f8a_7528_1e37_32a9))
    );
    assert_close(active.weights().values()[0], 1.0, 2.0e-7);
    assert_close(active.objective(), 1.0, 3.0e-7);
    assert!(active.residual().minimum_lower_slack() >= -1.0e-7);

    let two_variables = solve_loqo_qp(
        &matrix(2, 2, &[1.0, 0.0, 0.0, 2.0]),
        &matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]),
        &vector(&[0.5, -2.0]),
        &vector(&[1.0, 1.0]),
    )
    .unwrap();
    assert_eq!(
        two_variables
            .weights()
            .values()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        vec![0x3fe0_0001_5cee_f1cc, 0xbff0_0000_0ab0_9aeb]
    );
    assert_eq!(two_variables.objective().to_bits(), 0x4002_0000_6c9c_f607);
    assert_eq!(two_variables.trace().len(), 10);
    assert_close(two_variables.weights().values()[0], 0.5, 1.0e-6);
    assert_close(two_variables.weights().values()[1], -1.0, 1.0e-6);
    assert!(two_variables.residual().accepted());

    let upper_active = solve_loqo_qp(
        &matrix(1, 1, &[1.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[-2.0]),
        &vector(&[1.0]),
    )
    .unwrap();
    assert_close(upper_active.weights().values()[0], -1.0, 2.0e-7);
    assert_close(upper_active.objective(), 1.0, 3.0e-7);
    assert!(upper_active.residual().minimum_upper_slack() >= -1.0e-7);
}

#[test]
fn positivity_divisor_preserves_the_frozen_point_ninety_five_fraction() {
    assert_eq!(loqo_step_divisor(&[190.0], &[100.0]).unwrap(), 2.0);
    assert_eq!(
        loqo_step_divisor(&[95.0, -190.0], &[100.0, 200.0]).unwrap(),
        1.0
    );
    assert_eq!(loqo_step_divisor(&[0.0], &[100.0]).unwrap(), 1.0);
    assert_eq!(
        loqo_step_divisor(&[f64::NAN], &[100.0]),
        Err(LoqoStepError::NonFiniteRatio)
    );
}

#[test]
fn validation_records_llt_without_using_it_as_a_solve_gate() {
    let interpolation = matrix(2, 2, &[1.0, 0.0, 0.0, -0.25]);
    let constraints = matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]);
    let lower = vector(&[0.25, -2.0]);
    let range = vector(&[0.5, 1.0]);
    let validation = validate_loqo_qp(&interpolation, &constraints, &lower, &range);
    assert!(!validation.surfe_matrix_system_valid());
    assert!(validation.safe_shape_valid());

    let outcome = solve_loqo_qp(&interpolation, &constraints, &lower, &range);
    match outcome {
        Ok(solution) => assert!(solution.attempted()),
        Err(error) => assert!(error.attempted()),
    }
}

#[test]
fn invalid_shapes_and_nonfinite_inputs_fail_with_typed_evidence() {
    let interpolation = matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]);
    let lower = vector(&[0.0, 0.0]);
    let range = vector(&[1.0, 1.0]);

    let wrong_rows = solve_loqo_qp(
        &interpolation,
        &matrix(1, 2, &[1.0, 0.0]),
        &vector(&[0.0]),
        &vector(&[1.0]),
    )
    .unwrap_err();
    assert_eq!(wrong_rows.kind(), LoqoSolveErrorKind::ConstraintRowMismatch);
    assert!(!wrong_rows.attempted());

    let wrong_range = solve_loqo_qp(
        &interpolation,
        &matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]),
        &lower,
        &vector(&[1.0]),
    )
    .unwrap_err();
    assert_eq!(wrong_range.kind(), LoqoSolveErrorKind::RangeValueMismatch);
    assert!(!wrong_range.attempted());

    let nonfinite = solve_loqo_qp(
        &matrix(1, 1, &[f64::NAN]),
        &matrix(1, 1, &[1.0]),
        &vector(&[0.0]),
        &vector(&[1.0]),
    )
    .unwrap_err();
    assert_eq!(nonfinite.kind(), LoqoSolveErrorKind::NonFiniteInput);
    assert!(nonfinite.attempted());
    assert_eq!(nonfinite.trace().len(), 1);
    assert!(nonfinite.kkt_failure().is_some());
    assert_eq!(nonfinite.surfe_error(), Error::LoqoSolverFailure);

    let valid = validate_loqo_qp(
        &interpolation,
        &matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]),
        &lower,
        &range,
    );
    assert!(valid.safe_shape_valid());
}

#[test]
fn zero_hessian_is_attempted_and_preserves_the_frozen_finite_candidate() {
    let solution = solve_loqo_qp(
        &matrix(1, 1, &[0.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[-1.0]),
        &vector(&[2.0]),
    )
    .unwrap();
    assert!(!solution.validation().surfe_matrix_system_valid());
    assert_eq!(
        solution.weights().values()[0].to_bits(),
        0xbcc7_d498_0000_0000
    );
    assert_eq!(solution.objective().to_bits(), (-0.0_f64).to_bits());
    assert_eq!(solution.trace().len(), 10);
    assert!(solution.residual().accepted());
}

#[test]
fn ill_conditioned_hessian_is_attempted_without_a_condition_number_gate() {
    let solution = solve_loqo_qp(
        &matrix(1, 1, &[1.0e-14]),
        &matrix(1, 1, &[1.0]),
        &vector(&[0.25]),
        &vector(&[0.5]),
    )
    .unwrap();
    assert!(solution.validation().surfe_matrix_system_valid());
    assert_eq!(
        solution.weights().values()[0].to_bits(),
        0x3fe0_1792_49ba_14bb
    );
    assert_eq!(solution.objective().to_bits(), 0x3ce6_c725_68a9_5775);
    assert_eq!(solution.trace().len(), 9);
    assert!(solution.residual().accepted());
    assert!(solution.residual().complementarity() <= solution.residual().residual_limit());
}

#[test]
fn frozen_tight_indefinite_and_impossible_stops_are_typed_failures() {
    let tight = solve_loqo_qp(
        &matrix(1, 1, &[1.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[0.5]),
        &vector(&[0.0]),
    )
    .unwrap_err();
    assert_eq!(tight.kind(), LoqoSolveErrorKind::DualObjectiveAbovePrimal);
    assert_eq!(tight.trace().len(), 1);

    let indefinite = solve_loqo_qp(
        &matrix(1, 1, &[-0.25]),
        &matrix(1, 1, &[1.0]),
        &vector(&[0.25]),
        &vector(&[0.5]),
    )
    .unwrap_err();
    assert!(!indefinite.validation().surfe_matrix_system_valid());
    assert_eq!(
        indefinite.kind(),
        LoqoSolveErrorKind::DualObjectiveAbovePrimal
    );
    assert_eq!(indefinite.trace().len(), 6);
    assert_close(
        indefinite.trace()[5].primal_objective(),
        -0.131_859_279_518_660_28,
        2.0e-14,
    );
    assert_close(
        indefinite.trace()[5].dual_objective(),
        -0.129_415_997_969_703_38,
        2.0e-14,
    );

    let impossible = solve_loqo_qp(
        &matrix(1, 1, &[1.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[0.0]),
        &vector(&[-1.0]),
    )
    .unwrap_err();
    assert_eq!(
        impossible.kind(),
        LoqoSolveErrorKind::DualObjectiveAbovePrimal
    );
    assert_eq!(impossible.trace().len(), 1);
    assert!(impossible.residual().unwrap().minimum_upper_slack() < 0.0);
}

#[test]
fn explicit_iteration_limit_is_failure_and_keeps_the_candidate() {
    let error = solve_loqo_qp_with_options(
        &matrix(1, 1, &[1.0]),
        &matrix(1, 1, &[1.0]),
        &vector(&[-1.0]),
        &vector(&[2.0]),
        LoqoOptions { max_iterations: 0 },
    )
    .unwrap_err();
    assert_eq!(error.kind(), LoqoSolveErrorKind::IterationLimit);
    assert!(error.attempted());
    assert!(error.candidate_weights().is_some());
    assert_eq!(error.trace().len(), 1);
}
