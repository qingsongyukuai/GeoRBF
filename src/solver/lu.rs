//! Partial-pivot LU matching Surfe's Eigen solve path for dynamic vectors.
//!
//! Sources:
//! - `surfe_lib/matrix_solver.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`;
//! - `Eigen/src/LU/PartialPivLU.h@36b95962756c1fce8e29b1f8bc45967f30773c00`;
//! - `Eigen/src/Core/products/TriangularSolverVector.h@36b95962756c1fce8e29b1f8bc45967f30773c00`.

use crate::{DenseMatrix, DenseVector};

use super::{LuSolveError, LuSolveErrorKind};

/// Matrix/RHS facts kept separate from whether the solve was attempted.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LuValidation {
    storage_length_matches: bool,
    non_empty: bool,
    square: bool,
    right_hand_side_length_matches: bool,
    matrix_finite: bool,
    right_hand_side_finite: bool,
}

impl LuValidation {
    pub const fn storage_length_matches(self) -> bool {
        self.storage_length_matches
    }

    pub const fn is_non_empty(self) -> bool {
        self.non_empty
    }

    pub const fn is_square(self) -> bool {
        self.square
    }

    pub const fn right_hand_side_length_matches(self) -> bool {
        self.right_hand_side_length_matches
    }

    pub const fn matrix_is_finite(self) -> bool {
        self.matrix_finite
    }

    pub const fn right_hand_side_is_finite(self) -> bool {
        self.right_hand_side_finite
    }

    /// Exact observable predicate of frozen `validate_matrix_systems()` once
    /// valid dense storage has reached its constructor.
    pub const fn surfe_matrix_system_valid(self) -> bool {
        self.storage_length_matches && self.matrix_finite
    }

    /// Safe Rust preflight, retained as evidence rather than a condition-number
    /// gate. The solve still attempts finite ill-conditioned systems.
    pub const fn safe_preflight_valid(self) -> bool {
        self.storage_length_matches
            && self.non_empty
            && self.square
            && self.right_hand_side_length_matches
            && self.matrix_finite
            && self.right_hand_side_finite
    }
}

/// Pivots and packed `L/U` values from the actual attempted factorization.
#[derive(Clone, Debug, PartialEq)]
pub struct LuFactorizationEvidence {
    row_transpositions: Vec<usize>,
    permutation: Vec<usize>,
    pivot_values: Vec<f64>,
    exact_zero_pivot: Option<usize>,
    packed_lu: DenseMatrix,
}

impl LuFactorizationEvidence {
    pub fn row_transpositions(&self) -> &[usize] {
        &self.row_transpositions
    }

    /// Eigen `PermutationMatrix::indices()` representation of the row swaps.
    pub fn permutation(&self) -> &[usize] {
        &self.permutation
    }

    pub fn pivot_values(&self) -> &[f64] {
        &self.pivot_values
    }

    pub const fn exact_zero_pivot(&self) -> Option<usize> {
        self.exact_zero_pivot
    }

    pub const fn packed_lu(&self) -> &DenseMatrix {
        &self.packed_lu
    }
}

/// Scale-aware post-solve evidence; no condition-number estimate is used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuResidualEvidence {
    l2: f64,
    relative_l2: f64,
    linf: f64,
    backward_error: f64,
    acceptance_limit: f64,
    accepted: bool,
}

impl LuResidualEvidence {
    pub const fn l2(self) -> f64 {
        self.l2
    }

    pub const fn relative_l2(self) -> f64 {
        self.relative_l2
    }

    pub const fn linf(self) -> f64 {
        self.linf
    }

    pub const fn backward_error(self) -> f64 {
        self.backward_error
    }

    pub const fn acceptance_limit(self) -> f64 {
        self.acceptance_limit
    }

    pub const fn accepted(self) -> bool {
        self.accepted
    }
}

/// Successful finite LU weights with pivot and residual evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct LuSolution {
    validation: LuValidation,
    weights: DenseVector,
    factorization: LuFactorizationEvidence,
    residual: LuResidualEvidence,
}

impl LuSolution {
    pub const fn attempted(&self) -> bool {
        true
    }

    pub const fn validation(&self) -> LuValidation {
        self.validation
    }

    pub const fn weights(&self) -> &DenseVector {
        &self.weights
    }

    pub const fn factorization(&self) -> &LuFactorizationEvidence {
        &self.factorization
    }

    pub const fn residual(&self) -> LuResidualEvidence {
        self.residual
    }
}

/// Inspect an explicit row-major system without solving it.
pub fn validate_lu_system(
    rows: usize,
    columns: usize,
    matrix_row_major: &[f64],
    right_hand_side: &[f64],
) -> LuValidation {
    LuValidation {
        storage_length_matches: rows
            .checked_mul(columns)
            .is_some_and(|length| length == matrix_row_major.len()),
        non_empty: rows != 0 && columns != 0,
        square: rows == columns,
        right_hand_side_length_matches: right_hand_side.len() == rows,
        matrix_finite: matrix_row_major.iter().all(|value| value.is_finite()),
        right_hand_side_finite: right_hand_side.iter().all(|value| value.is_finite()),
    }
}

/// Solve an explicit row-major system with frozen partial-pivot and triangular
/// update order. Ill-conditioned matrices are always attempted.
pub fn solve_partial_pivot_lu(
    rows: usize,
    columns: usize,
    matrix_row_major: &[f64],
    right_hand_side: &[f64],
) -> Result<LuSolution, LuSolveError> {
    let validation = validate_lu_system(rows, columns, matrix_row_major, right_hand_side);
    if !validation.storage_length_matches {
        return Err(early_error(
            LuSolveErrorKind::InvalidStorageLength,
            validation,
        ));
    }
    if !validation.non_empty {
        return Err(early_error(LuSolveErrorKind::EmptySystem, validation));
    }
    if !validation.square {
        return Err(early_error(LuSolveErrorKind::NonSquareMatrix, validation));
    }
    if !validation.right_hand_side_length_matches {
        return Err(early_error(LuSolveErrorKind::DimensionMismatch, validation));
    }

    let mut matrix = DenseMatrix::zeros(rows, columns);
    let mut row = 0;
    while row < rows {
        let mut column = 0;
        while column < columns {
            matrix.set(row, column, matrix_row_major[row * columns + column]);
            column += 1;
        }
        row += 1;
    }
    solve_dense_partial_pivot_lu_with_validation(&matrix, right_hand_side, validation)
}

/// Solve T17 dense assembly values without copying through a public raw-matrix
/// constructor.
pub fn solve_dense_partial_pivot_lu(
    matrix: &DenseMatrix,
    right_hand_side: &DenseVector,
) -> Result<LuSolution, LuSolveError> {
    let validation = validate_lu_system(
        matrix.rows(),
        matrix.cols(),
        matrix.data(),
        right_hand_side.values(),
    );
    if !validation.non_empty {
        return Err(early_error(LuSolveErrorKind::EmptySystem, validation));
    }
    if !validation.square {
        return Err(early_error(LuSolveErrorKind::NonSquareMatrix, validation));
    }
    if !validation.right_hand_side_length_matches {
        return Err(early_error(LuSolveErrorKind::DimensionMismatch, validation));
    }
    solve_dense_partial_pivot_lu_with_validation(matrix, right_hand_side.values(), validation)
}

fn solve_dense_partial_pivot_lu_with_validation(
    matrix: &DenseMatrix,
    right_hand_side: &[f64],
    validation: LuValidation,
) -> Result<LuSolution, LuSolveError> {
    let factorization = factorize(matrix);
    let mut weights = right_hand_side.to_vec();

    let mut index = 0;
    while index < factorization.row_transpositions.len() {
        weights.swap(index, factorization.row_transpositions[index]);
        index += 1;
    }

    // Eigen's dynamic column-major UnitLower vector solve updates the tail of
    // the RHS after each nonzero entry.
    let size = weights.len();
    let mut column = 0;
    while column < size {
        if weights[column] != 0.0 {
            let multiplier = weights[column];
            let mut row = column + 1;
            while row < size {
                weights[row] -= multiplier * packed(&factorization, row, column);
                row += 1;
            }
        }
        column += 1;
    }

    // Frozen Eigen uses the analogous column update for the Upper solve and
    // intentionally skips division when the current RHS entry is exactly 0.
    let mut remaining = size;
    while remaining != 0 {
        let pivot = remaining - 1;
        if weights[pivot] != 0.0 {
            weights[pivot] /= packed(&factorization, pivot, pivot);
            let solved = weights[pivot];
            let mut row = 0;
            while row < pivot {
                weights[row] -= solved * packed(&factorization, row, pivot);
                row += 1;
            }
        }
        remaining -= 1;
    }

    if !validation.matrix_finite {
        return Err(attempted_error(
            LuSolveErrorKind::NonFiniteMatrix,
            validation,
            factorization,
            None,
            None,
        ));
    }
    if !validation.right_hand_side_finite {
        return Err(attempted_error(
            LuSolveErrorKind::NonFiniteRightHandSide,
            validation,
            factorization,
            None,
            None,
        ));
    }
    if weights.iter().any(|value| !value.is_finite()) {
        let kind = if factorization.exact_zero_pivot.is_some() {
            LuSolveErrorKind::SingularSystem
        } else {
            LuSolveErrorKind::NonFiniteSolution
        };
        return Err(attempted_error(kind, validation, factorization, None, None));
    }

    let residual = residual_evidence(matrix, right_hand_side, &weights);
    let candidate_weights = DenseVector::from_values(weights);
    if !residual.l2.is_finite()
        || !residual.relative_l2.is_finite()
        || !residual.linf.is_finite()
        || !residual.backward_error.is_finite()
    {
        return Err(attempted_error(
            LuSolveErrorKind::NonFiniteResidual,
            validation,
            factorization,
            Some(candidate_weights),
            Some(residual),
        ));
    }
    if !residual.accepted {
        return Err(attempted_error(
            LuSolveErrorKind::ResidualTooLarge,
            validation,
            factorization,
            Some(candidate_weights),
            Some(residual),
        ));
    }

    Ok(LuSolution {
        validation,
        weights: candidate_weights,
        factorization,
        residual,
    })
}

fn factorize(matrix: &DenseMatrix) -> LuFactorizationEvidence {
    let size = matrix.rows();
    let mut packed_lu = matrix.clone();
    let mut row_transpositions = Vec::with_capacity(size);
    let mut pivot_values = Vec::with_capacity(size);
    let mut exact_zero_pivot = None;

    let mut pivot_column = 0;
    while pivot_column < size {
        let mut pivot_row = pivot_column;
        let mut biggest = packed_lu
            .get(pivot_column, pivot_column)
            .expect("validated square matrix")
            .abs();
        let mut candidate = pivot_column + 1;
        while candidate < size {
            let score = packed_lu
                .get(candidate, pivot_column)
                .expect("validated square matrix")
                .abs();
            if score > biggest {
                biggest = score;
                pivot_row = candidate;
            }
            candidate += 1;
        }
        row_transpositions.push(pivot_row);

        if biggest != 0.0 {
            if pivot_row != pivot_column {
                swap_rows(&mut packed_lu, pivot_column, pivot_row);
            }
            let pivot = packed_lu
                .get(pivot_column, pivot_column)
                .expect("validated square matrix");
            let mut row = pivot_column + 1;
            while row < size {
                let value = packed_lu
                    .get(row, pivot_column)
                    .expect("validated square matrix")
                    / pivot;
                packed_lu.set(row, pivot_column, value);
                row += 1;
            }
        } else if exact_zero_pivot.is_none() {
            exact_zero_pivot = Some(pivot_column);
        }

        let mut row = pivot_column + 1;
        while row < size {
            let lower = packed_lu
                .get(row, pivot_column)
                .expect("validated square matrix");
            let mut column = pivot_column + 1;
            while column < size {
                let value = packed_lu.get(row, column).expect("validated square matrix")
                    - lower
                        * packed_lu
                            .get(pivot_column, column)
                            .expect("validated square matrix");
                packed_lu.set(row, column, value);
                column += 1;
            }
            row += 1;
        }
        pivot_values.push(
            packed_lu
                .get(pivot_column, pivot_column)
                .expect("validated square matrix"),
        );
        pivot_column += 1;
    }

    let permutation = eigen_permutation_indices(&row_transpositions);
    LuFactorizationEvidence {
        row_transpositions,
        permutation,
        pivot_values,
        exact_zero_pivot,
        packed_lu,
    }
}

fn swap_rows(matrix: &mut DenseMatrix, first: usize, second: usize) {
    let mut column = 0;
    while column < matrix.cols() {
        let left = matrix.get(first, column).expect("row in bounds");
        let right = matrix.get(second, column).expect("row in bounds");
        matrix.set(first, column, right);
        matrix.set(second, column, left);
        column += 1;
    }
}

fn eigen_permutation_indices(transpositions: &[usize]) -> Vec<usize> {
    let mut row_order = (0..transpositions.len()).collect::<Vec<_>>();
    let mut index = 0;
    while index < transpositions.len() {
        row_order.swap(index, transpositions[index]);
        index += 1;
    }
    let mut inverse = vec![0; row_order.len()];
    let mut position = 0;
    while position < row_order.len() {
        inverse[row_order[position]] = position;
        position += 1;
    }
    inverse
}

fn packed(evidence: &LuFactorizationEvidence, row: usize, column: usize) -> f64 {
    evidence
        .packed_lu
        .get(row, column)
        .expect("factorization indices are in bounds")
}

fn residual_evidence(
    matrix: &DenseMatrix,
    right_hand_side: &[f64],
    weights: &[f64],
) -> LuResidualEvidence {
    let mut squared_l2 = 0.0;
    let mut linf = 0.0_f64;
    let mut matrix_linf = 0.0_f64;
    let mut row = 0;
    while row < matrix.rows() {
        let mut prediction = 0.0;
        let mut row_sum = 0.0;
        let mut column = 0;
        while column < matrix.cols() {
            let coefficient = matrix.get(row, column).expect("matrix index in bounds");
            prediction += coefficient * weights[column];
            row_sum += coefficient.abs();
            column += 1;
        }
        let residual = prediction - right_hand_side[row];
        squared_l2 += residual * residual;
        linf = linf.max(residual.abs());
        matrix_linf = matrix_linf.max(row_sum);
        row += 1;
    }

    let mut rhs_squared_l2 = 0.0;
    let mut rhs_linf = 0.0_f64;
    for value in right_hand_side {
        rhs_squared_l2 += value * value;
        rhs_linf = rhs_linf.max(value.abs());
    }
    let weights_linf = weights
        .iter()
        .fold(0.0_f64, |maximum, value| maximum.max(value.abs()));
    let l2 = squared_l2.sqrt();
    let rhs_l2 = rhs_squared_l2.sqrt();
    let relative_l2 = if rhs_l2 == 0.0 {
        if l2 == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        l2 / rhs_l2
    };
    let denominator = matrix_linf * weights_linf + rhs_linf;
    let backward_error = if denominator == 0.0 {
        if linf == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        linf / denominator
    };
    let acceptance_limit = 64.0 * f64::EPSILON * matrix.rows().max(1) as f64;
    let accepted = backward_error.is_finite() && backward_error <= acceptance_limit;
    LuResidualEvidence {
        l2,
        relative_l2,
        linf,
        backward_error,
        acceptance_limit,
        accepted,
    }
}

fn early_error(kind: LuSolveErrorKind, validation: LuValidation) -> LuSolveError {
    LuSolveError {
        kind,
        attempted: false,
        validation,
        factorization: None,
        candidate_weights: None,
        residual: None,
    }
}

fn attempted_error(
    kind: LuSolveErrorKind,
    validation: LuValidation,
    factorization: LuFactorizationEvidence,
    candidate_weights: Option<DenseVector>,
    residual: Option<LuResidualEvidence>,
) -> LuSolveError {
    LuSolveError {
        kind,
        attempted: true,
        validation,
        factorization: Some(Box::new(factorization)),
        candidate_weights: candidate_weights.map(Box::new),
        residual,
    }
}
