use georbf::{
    predictor_corrector_step_length, solve_predictor_corrector_qp,
    solve_predictor_corrector_qp_with_options, validate_predictor_corrector_qp, ConstraintSystem,
    DenseMatrix, DenseVector, Error, QpOptions, QpSolveErrorKind, QpStopReason,
};

fn matrix(rows: usize, columns: usize, values: &[f64]) -> DenseMatrix {
    DenseMatrix::from_row_major(rows, columns, values.to_vec()).expect("valid test matrix")
}

fn vector(values: &[f64]) -> DenseVector {
    DenseVector::from_values(values.to_vec())
}

fn system(rows: usize, columns: usize, values: &[f64], rhs: &[f64]) -> ConstraintSystem {
    ConstraintSystem::new(matrix(rows, columns, values), vector(rhs))
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn assert_close(left: f64, right: f64, tolerance: f64) {
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    assert!(
        (left - right).abs() <= tolerance * scale,
        "{left:.17e} != {right:.17e}"
    );
}

#[test]
fn inactive_and_active_bounds_match_frozen_weights_objective_and_mu_trace() {
    let interpolation = matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]);
    let equality = system(1, 2, &[1.0, 1.0], &[1.0]);

    let inactive = solve_predictor_corrector_qp(
        &interpolation,
        &equality,
        &system(1, 2, &[1.0, 0.0], &[0.2]),
    )
    .unwrap();
    assert_eq!(
        bits(inactive.weights().values()),
        vec![0x3fe0_0000_0000_000c, 0x3fdf_ffff_ffff_ffe8]
    );
    assert_eq!(
        inactive.stop_reason(),
        QpStopReason::ComplementarityTolerance
    );
    assert_eq!(inactive.trace().len(), 12);
    assert_close(inactive.trace()[0].mu(), 1_001_414.679_609_126_1, 2.0e-13);
    assert_close(inactive.trace()[10].mu(), 4.664_932_464_530_536e-6, 2.0e-13);
    assert_close(
        inactive.trace()[11].mu(),
        1.566_468_512_521_041e-15,
        2.0e-13,
    );
    assert_eq!(inactive.objective().to_bits(), 0x3fe0_0000_0000_0000);
    assert_eq!(inactive.residual().equality_linf(), 0.0);
    assert_close(inactive.residual().minimum_inequality_slack(), 0.3, 2.0e-14);
    assert!(inactive.residual().accepted());

    let active = solve_predictor_corrector_qp(
        &interpolation,
        &equality,
        &system(1, 2, &[1.0, 0.0], &[0.8]),
    )
    .unwrap();
    assert_eq!(
        bits(active.weights().values()),
        vec![0x3fe9_9999_9999_99a5, 0x3fc9_9999_9999_996b]
    );
    assert_eq!(active.trace().len(), 12);
    assert_close(active.trace()[0].mu(), 1_001_414.272_315_62, 2.0e-13);
    assert!(active.residual().minimum_inequality_slack() >= -2.0e-14);
    assert_close(active.objective(), 0.68, 2.0e-14);
}

#[test]
fn no_equality_and_mixed_qps_preserve_empty_blocks_and_multiple_inequalities() {
    let no_equality = solve_predictor_corrector_qp(
        &matrix(1, 1, &[1.0]),
        &system(0, 1, &[], &[]),
        &system(1, 1, &[1.0], &[1.0]),
    )
    .unwrap();
    assert_eq!(
        bits(no_equality.weights().values()),
        vec![0x3ff0_0000_0000_3bd2]
    );
    assert_eq!(no_equality.dual_equality().len(), 0);
    assert_eq!(no_equality.trace().len(), 11);
    assert_close(
        no_equality.trace()[10].mu(),
        6.800_999_406_677_297e-12,
        2.0e-13,
    );

    let mixed = solve_predictor_corrector_qp(
        &matrix(3, 3, &[2.0, 0.5, 0.0, 0.5, 1.5, 0.25, 0.0, 0.25, 1.0]),
        &system(1, 3, &[1.0, 1.0, 1.0], &[1.0]),
        &system(2, 3, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0], &[0.1, 0.2]),
    )
    .unwrap();
    assert_eq!(
        bits(mixed.weights().values()),
        vec![
            0x3fce_a5db_e911_9e9d,
            0x3fcc_d856_9d6f_bf12,
            0x3fe1_2073_5e5f_a895,
        ]
    );
    assert_eq!(mixed.trace().len(), 14);
    assert_close(mixed.trace()[13].mu(), 3.567_704_156_896_468_3e-10, 2.0e-13);
    assert!(mixed.residual().accepted());
}

#[test]
fn step_length_keeps_frozen_ratio_sign_threshold_and_unit_cap() {
    assert_eq!(
        predictor_corrector_step_length(&[2.0, 4.0], &[4.0, -1.0], &[3.0, -2.0], &[6.0, -1.0])
            .unwrap(),
        0.5
    );
    assert_eq!(
        predictor_corrector_step_length(&[1.0], &[1.0e14], &[1.0], &[0.0]).unwrap(),
        1.0
    );
    assert_eq!(
        predictor_corrector_step_length(&[-2.0], &[4.0], &[-3.0], &[6.0]).unwrap(),
        0.5
    );
}

#[test]
fn validation_records_frozen_llt_result_without_gating_attempted_solve() {
    let interpolation = matrix(2, 2, &[1.0, 0.0, 0.0, -0.25]);
    let equality = system(1, 2, &[1.0, 1.0], &[1.0]);
    let inequality = system(1, 2, &[1.0, 0.0], &[0.2]);
    let validation = validate_predictor_corrector_qp(&interpolation, &equality, &inequality);
    assert!(!validation.surfe_matrix_system_valid());
    assert!(validation.safe_shape_valid());

    let solution = solve_predictor_corrector_qp(&interpolation, &equality, &inequality).unwrap();
    assert!(solution.attempted());
    assert_eq!(
        bits(solution.weights().values()),
        vec![0x3fc9_9999_9999_9a99, 0x3fe9_9999_9999_995a]
    );
    assert_eq!(solution.trace().len(), 12);
}

#[test]
fn ill_conditioned_system_is_attempted_and_matches_frozen_prediction() {
    let solution = solve_predictor_corrector_qp(
        &matrix(2, 2, &[1.0e-14, 0.0, 0.0, 1.0]),
        &system(1, 2, &[1.0, 1.0], &[1.0]),
        &system(1, 2, &[1.0, 0.0], &[0.25]),
    )
    .unwrap();
    assert_eq!(
        bits(solution.weights().values()),
        vec![0x3ff0_0000_0106_9008, 0xbe30_6900_826c_0000]
    );
    assert!(solution.residual().accepted());
    assert_eq!(solution.trace().len(), 11);
}

#[test]
fn degenerate_zero_hessian_preserves_frozen_finite_nonunique_candidate() {
    let solution = solve_predictor_corrector_qp(
        &matrix(2, 2, &[0.0, 0.0, 0.0, 0.0]),
        &system(1, 2, &[1.0, 0.0], &[1.0]),
        &system(1, 2, &[0.0, 1.0], &[0.0]),
    )
    .unwrap();
    assert!(!solution.validation().surfe_matrix_system_valid());
    assert_eq!(bits(solution.weights().values()), bits(&[1.0, 1000.0]));
    assert_eq!(solution.objective(), 0.0);
    assert_eq!(solution.trace().len(), 2);
    assert_eq!(solution.trace()[0].mu(), 1_000_000.0);
    assert_eq!(solution.trace()[1].mu(), 0.0);
    assert!(solution.residual().accepted());
}

#[test]
fn frozen_infeasible_success_is_observed_then_safely_rejected() {
    let error = solve_predictor_corrector_qp(
        &matrix(1, 1, &[1.0]),
        &system(1, 1, &[1.0], &[0.0]),
        &system(1, 1, &[1.0], &[1.0]),
    )
    .unwrap_err();
    assert_eq!(error.kind(), QpSolveErrorKind::InfeasibleSolution);
    assert!(error.attempted());
    assert_eq!(error.surfe_error(), Error::PredictorCorrectorSolverFailure);
    assert_eq!(
        error.stop_reason(),
        Some(QpStopReason::ComplementarityTolerance)
    );
    assert_eq!(error.trace().len(), 2);
    assert_close(error.trace()[0].mu(), 1_005_418.627_775_935_5, 2.0e-13);
    assert_close(error.trace()[1].mu(), -1_006.422_036_571_875_7, 2.0e-13);
    assert_eq!(error.candidate_weights().unwrap().values(), &[0.0]);
    assert_eq!(error.residual().unwrap().minimum_inequality_slack(), -1.0);
}

#[test]
fn sixth_iteration_mu_increase_branch_is_preserved_then_feasibility_checked() {
    let error = solve_predictor_corrector_qp(
        &matrix(
            2,
            2,
            &[
                -1.295_358_919_297_941_1,
                0.593_463_939_355_649_6,
                0.593_463_939_355_649_6,
                1.336_850_734_178_098_9,
            ],
        ),
        &system(
            1,
            2,
            &[1.246_943_181_585_774_3, 0.343_801_132_568_106_3],
            &[-0.779_386_113_706_520_8],
        ),
        &system(
            2,
            2,
            &[
                1.845_105_482_275_249_3,
                -0.706_794_151_743_586_1,
                -0.786_740_860_313_468_9,
                0.053_180_649_095_110_78,
            ],
            &[-1.059_091_143_660_698_5, 0.100_628_356_048_078_96],
        ),
    )
    .unwrap_err();
    assert_eq!(error.kind(), QpSolveErrorKind::InfeasibleSolution);
    assert_eq!(
        error.stop_reason(),
        Some(QpStopReason::ComplementarityIncreased)
    );
    assert_eq!(error.trace().len(), 7);
    assert_close(error.trace()[5].mu(), 17.105_335_153_876_93, 2.0e-13);
    assert_close(error.trace()[6].mu(), 17.386_342_564_804_785, 2.0e-13);
    assert_eq!(
        bits(error.candidate_weights().unwrap().values()),
        vec![0xbfc1_3154_0b55_8616, 0xbffc_7a15_f2da_55ca]
    );
}

#[test]
fn invalid_shapes_missing_inequalities_and_nonfinite_inputs_fail_with_evidence() {
    let interpolation = matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]);
    let equality = system(1, 3, &[1.0, 1.0, 1.0], &[1.0]);
    let inequality = system(1, 2, &[1.0, 0.0], &[0.2]);
    let shape_error =
        solve_predictor_corrector_qp(&interpolation, &equality, &inequality).unwrap_err();
    assert_eq!(shape_error.kind(), QpSolveErrorKind::EqualityColumnMismatch);
    assert!(!shape_error.attempted());

    let missing = solve_predictor_corrector_qp(
        &interpolation,
        &system(0, 2, &[], &[]),
        &system(0, 2, &[], &[]),
    )
    .unwrap_err();
    assert_eq!(missing.kind(), QpSolveErrorKind::MissingInequalities);
    assert!(!missing.attempted());

    let nonfinite = solve_predictor_corrector_qp(
        &matrix(1, 1, &[f64::NAN]),
        &system(0, 1, &[], &[]),
        &system(1, 1, &[1.0], &[0.0]),
    )
    .unwrap_err();
    assert_eq!(nonfinite.kind(), QpSolveErrorKind::NonFiniteInput);
    assert!(nonfinite.attempted());
    assert!(nonfinite.kkt_failure().is_some());
}

#[test]
fn explicit_iteration_limit_is_a_failure_not_a_relaxed_success() {
    let error = solve_predictor_corrector_qp_with_options(
        &matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]),
        &system(1, 2, &[1.0, 1.0], &[1.0]),
        &system(1, 2, &[1.0, 0.0], &[0.2]),
        QpOptions { max_iterations: 0 },
    )
    .unwrap_err();
    assert_eq!(error.kind(), QpSolveErrorKind::IterationLimit);
    assert!(error.attempted());
    assert_eq!(error.trace().len(), 1);
    assert!(error.candidate_weights().is_some());
}
