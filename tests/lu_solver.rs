use georbf::{solve_partial_pivot_lu, validate_lu_system, Error, LuSolveErrorKind};

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
fn well_conditioned_system_matches_frozen_weights_pivots_and_packed_lu() {
    let matrix = [4.0, 1.0, 1.0, 3.0];
    let rhs = [1.0, 2.0];
    let validation = validate_lu_system(2, 2, &matrix, &rhs);
    assert!(validation.surfe_matrix_system_valid());
    assert!(validation.safe_preflight_valid());

    let solution = solve_partial_pivot_lu(2, 2, &matrix, &rhs).unwrap();
    assert!(solution.attempted());
    assert_eq!(
        bits(solution.weights().values()),
        vec![0x3fb7_45d1_745d_1746, 0x3fe4_5d17_45d1_745d]
    );
    assert_eq!(solution.factorization().row_transpositions(), &[0, 1]);
    assert_eq!(solution.factorization().permutation(), &[0, 1]);
    assert_eq!(solution.factorization().exact_zero_pivot(), None);
    assert_eq!(
        bits(solution.factorization().packed_lu().data()),
        bits(&[4.0, 1.0, 0.25, 2.75])
    );
    assert_eq!(solution.residual().l2(), 0.0);
    assert_eq!(solution.residual().relative_l2(), 0.0);
    assert_eq!(solution.residual().backward_error(), 0.0);
    assert!(solution.residual().accepted());
}

#[test]
fn hilbert_system_is_attempted_without_condition_number_rejection() {
    let matrix = [
        1.0,
        1.0 / 2.0,
        1.0 / 3.0,
        1.0 / 2.0,
        1.0 / 3.0,
        1.0 / 4.0,
        1.0 / 3.0,
        1.0 / 4.0,
        1.0 / 5.0,
    ];
    // Exact binary64 values emitted by frozen Eigen for `hilbert * Ones()`.
    let rhs = [
        f64::from_bits(0x3ffd_5555_5555_5555),
        f64::from_bits(0x3ff1_5555_5555_5555),
        f64::from_bits(0x3fe9_1111_1111_1111),
    ];
    let solution = solve_partial_pivot_lu(3, 3, &matrix, &rhs).unwrap();

    assert_eq!(solution.factorization().row_transpositions(), &[0, 2, 2]);
    assert_eq!(solution.factorization().permutation(), &[0, 2, 1]);
    let frozen = [
        f64::from_bits(0x3ff0_0000_0000_0008),
        f64::from_bits(0x3fef_ffff_ffff_ffa3),
        f64::from_bits(0x3ff0_0000_0000_002d),
    ];
    assert_eq!(bits(solution.weights().values()), bits(&frozen));
    assert_eq!(
        bits(solution.factorization().packed_lu().data()),
        vec![
            0x3ff0_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3fd5_5555_5555_5555,
            0x3fd5_5555_5555_5555,
            0x3fb5_5555_5555_5556,
            0x3fb6_c16c_16c1_6c18,
            0x3fe0_0000_0000_0000,
            0x3fef_ffff_ffff_fffd,
            0xbf76_c16c_16c1_6c00,
        ]
    );
    assert!(solution.residual().l2() <= 2.0_f64.powi(-52));
    assert!(solution.residual().accepted());
}

#[test]
fn equal_absolute_pivot_uses_the_first_row_like_frozen_eigen() {
    let matrix = [0.0, 2.0, 1.0, 3.0, 1.0, 4.0, -3.0, 5.0, 2.0];
    let expected = [1.0, -2.0, 0.5];
    let rhs = [
        matrix[0] * expected[0] + matrix[1] * expected[1] + matrix[2] * expected[2],
        matrix[3] * expected[0] + matrix[4] * expected[1] + matrix[5] * expected[2],
        matrix[6] * expected[0] + matrix[7] * expected[1] + matrix[8] * expected[2],
    ];
    let solution = solve_partial_pivot_lu(3, 3, &matrix, &rhs).unwrap();

    assert_eq!(solution.factorization().row_transpositions(), &[1, 2, 2]);
    assert_eq!(solution.factorization().permutation(), &[2, 0, 1]);
    assert_eq!(bits(solution.weights().values()), bits(&expected));
    assert_eq!(solution.residual().l2(), 0.0);
}

#[test]
fn consistent_singular_system_preserves_frozen_finite_solution_semantics() {
    let matrix = [1.0, 2.0, 2.0, 4.0];
    let solution = solve_partial_pivot_lu(2, 2, &matrix, &[3.0, 6.0]).unwrap();

    assert_eq!(solution.factorization().exact_zero_pivot(), Some(1));
    assert_eq!(solution.factorization().row_transpositions(), &[1, 1]);
    assert_eq!(solution.factorization().permutation(), &[1, 0]);
    assert_eq!(bits(solution.weights().values()), bits(&[3.0, 0.0]));
    assert_eq!(solution.residual().l2(), 0.0);
    assert!(solution.residual().accepted());
}

#[test]
fn inconsistent_singular_system_fails_after_attempt_with_linear_solver_category() {
    let matrix = [1.0, 2.0, 2.0, 4.0];
    let error = solve_partial_pivot_lu(2, 2, &matrix, &[3.0, 7.0]).unwrap_err();

    assert_eq!(error.kind(), LuSolveErrorKind::SingularSystem);
    assert!(error.attempted());
    assert_eq!(error.surfe_error(), Error::LinearSolverFailure);
    let factorization = error.factorization().unwrap();
    assert_eq!(factorization.exact_zero_pivot(), Some(1));
    assert_eq!(factorization.row_transpositions(), &[1, 1]);
}

#[test]
fn non_finite_matrix_and_rhs_keep_validation_and_attempt_evidence_distinct() {
    let mut matrix = [4.0, 1.0, 1.0, 3.0];
    matrix[0] = f64::NAN;
    let matrix_validation = validate_lu_system(2, 2, &matrix, &[1.0, 2.0]);
    assert!(!matrix_validation.surfe_matrix_system_valid());
    assert!(!matrix_validation.safe_preflight_valid());
    let matrix_error = solve_partial_pivot_lu(2, 2, &matrix, &[1.0, 2.0]).unwrap_err();
    assert_eq!(matrix_error.kind(), LuSolveErrorKind::NonFiniteMatrix);
    assert!(matrix_error.attempted());
    assert!(matrix_error.factorization().is_some());
    let zero_rhs_error = solve_partial_pivot_lu(2, 2, &matrix, &[0.0, 0.0]).unwrap_err();
    assert_eq!(zero_rhs_error.kind(), LuSolveErrorKind::NonFiniteMatrix);
    assert!(zero_rhs_error.attempted());

    let finite_matrix = [4.0, 1.0, 1.0, 3.0];
    let rhs = [f64::INFINITY, 2.0];
    let rhs_validation = validate_lu_system(2, 2, &finite_matrix, &rhs);
    assert!(rhs_validation.surfe_matrix_system_valid());
    assert!(!rhs_validation.safe_preflight_valid());
    let rhs_error = solve_partial_pivot_lu(2, 2, &finite_matrix, &rhs).unwrap_err();
    assert_eq!(rhs_error.kind(), LuSolveErrorKind::NonFiniteRightHandSide);
    assert!(rhs_error.attempted());
}

#[test]
fn invalid_shapes_fail_safely_before_attempting_factorization() {
    let mismatch =
        solve_partial_pivot_lu(2, 2, &[4.0, 1.0, 1.0, 3.0], &[1.0, 2.0, 3.0]).unwrap_err();
    assert_eq!(mismatch.kind(), LuSolveErrorKind::DimensionMismatch);
    assert!(!mismatch.attempted());

    let non_square =
        solve_partial_pivot_lu(2, 3, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0], &[1.0, 2.0]).unwrap_err();
    assert_eq!(non_square.kind(), LuSolveErrorKind::NonSquareMatrix);
    assert!(!non_square.attempted());

    let storage = solve_partial_pivot_lu(2, 2, &[1.0, 0.0, 0.0], &[1.0, 2.0]).unwrap_err();
    assert_eq!(storage.kind(), LuSolveErrorKind::InvalidStorageLength);
    assert!(!storage.attempted());

    let empty = solve_partial_pivot_lu(0, 0, &[], &[]).unwrap_err();
    assert_eq!(empty.kind(), LuSolveErrorKind::EmptySystem);
    assert!(!empty.attempted());
}

#[test]
fn deterministic_diagonally_dominant_systems_have_small_backward_error() {
    let mut state = 0x5eed_1234_89ab_cdef_u64;
    for size in 1..=12 {
        let mut matrix = vec![0.0; size * size];
        let mut expected = vec![0.0; size];
        for value in &mut expected {
            *value = random_unit(&mut state);
        }
        for row in 0..size {
            let mut off_diagonal_sum = 0.0;
            for column in 0..size {
                if row != column {
                    let value = random_unit(&mut state);
                    matrix[row * size + column] = value;
                    off_diagonal_sum += value.abs();
                }
            }
            matrix[row * size + row] = off_diagonal_sum + 1.0 + row as f64 / size as f64;
        }
        let rhs = (0..size)
            .map(|row| {
                (0..size)
                    .map(|column| matrix[row * size + column] * expected[column])
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let solution = solve_partial_pivot_lu(size, size, &matrix, &rhs).unwrap();
        for (actual, expected) in solution.weights().values().iter().zip(expected) {
            assert_close(*actual, expected, 2.0e-13);
        }
        assert!(solution.residual().accepted());
        assert!(solution.residual().backward_error() <= solution.residual().acceptance_limit());
    }
}

fn random_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    let mantissa = *state >> 11;
    (mantissa as f64 / ((1_u64 << 53) as f64)) * 2.0 - 1.0
}
