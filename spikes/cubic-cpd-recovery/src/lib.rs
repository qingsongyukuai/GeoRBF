use std::cmp::Ordering;
use std::error::Error;
use std::fmt::{Display, Formatter};

use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSettings, DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SolverStatus,
    ZeroConeT,
};
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::householder::{
    apply_block_householder_sequence_on_the_left_in_place_scratch,
    apply_block_householder_sequence_on_the_left_in_place_with_conj,
    apply_block_householder_sequence_transpose_on_the_left_in_place_scratch,
    apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj,
};
use faer::linalg::solvers::Qr;
use faer::linalg::solvers::Solve;
use faer::prelude::*;
use faer::{Conj, Side, get_global_parallelism};

const POLYNOMIAL_DIMENSION: usize = 4;

pub const FAER_VERSION: &str = "0.24.4";
pub const CLARABEL_VERSION: &str = "0.11.1";

#[derive(Debug, Clone)]
pub struct ExperimentEvidence {
    pub cpd: CpdEvidence,
    pub equality: EqualityEvidence,
    pub qp: ConvexRouteEvidence,
    pub socp: ConvexRouteEvidence,
    pub cross_route_observable_error: f64,
}

#[derive(Debug, Clone)]
pub struct CpdEvidence {
    pub functional_count: usize,
    pub polynomial_dimension: usize,
    pub polynomial_rank: usize,
    pub singular_values: Vec<f64>,
    pub null_space_defect: f64,
    pub reduced_symmetry_defect: f64,
    pub symmetry_defect_limit: f64,
    pub reduced_smallest_eigenvalue: f64,
    pub affine_reproduction_error: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InertiaEvidence {
    pub positive: usize,
    pub negative: usize,
    pub zero: usize,
}

#[derive(Debug, Clone)]
pub struct CanonicalObservables {
    pub field_coefficients: Vec<f64>,
    pub polynomial_coefficients: [f64; POLYNOMIAL_DIMENSION],
    pub semantic_latents: Vec<f64>,
    pub functional_values: Vec<f64>,
    pub residuals: Vec<f64>,
    pub slacks: Vec<f64>,
    pub side_condition_violation: f64,
    pub hard_violation: f64,
    pub field_energy: f64,
    pub objective: f64,
}

#[derive(Debug, Clone)]
pub struct EqualityEvidence {
    pub inertia: InertiaEvidence,
    pub expected_inertia: InertiaEvidence,
    pub normalized_backward_error: f64,
    pub scaling_round_trip_error: f64,
    pub manufactured_truth_error: f64,
    pub recovered: CanonicalObservables,
}

#[derive(Debug, Clone)]
pub struct ScaledConvexResiduals {
    pub primal: f64,
    pub dual: f64,
    pub stationarity: f64,
    pub complementarity: f64,
    pub relative_gap: f64,
}

#[derive(Debug, Clone)]
pub struct ConvexRouteEvidence {
    pub problem_class: &'static str,
    pub scaled: ScaledConvexResiduals,
    pub recovery_round_trip_error: f64,
    pub manufactured_truth_error: f64,
    pub recovered: CanonicalObservables,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FailureKind {
    PolynomialRankDeficient,
    ReducedPairingNotPositive,
    RecoveryVerification,
}

#[derive(Debug, Clone)]
pub struct FailureEvidence {
    pub kind: FailureKind,
    pub solver_invoked: bool,
    pub backend_contract_passed: bool,
    pub hidden_regularization_applied: bool,
    pub detected_violation: f64,
}

#[derive(Debug, Clone)]
pub struct CounterexampleEvidence {
    pub rank_deficient_polynomial: FailureEvidence,
    pub nonpositive_reduced_pairing: FailureEvidence,
    pub broken_recovery: FailureEvidence,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExperimentError {
    detail: String,
}

impl ExperimentError {
    fn new(detail: impl Into<String>) -> Self {
        Self {
            detail: detail.into(),
        }
    }
}

impl Display for ExperimentError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.detail)
    }
}

impl Error for ExperimentError {}

#[derive(Debug, Clone, Copy)]
struct Functional {
    point: [f64; 3],
    value_weight: f64,
    derivative_weight: [f64; 3],
}

impl Functional {
    fn polynomial_pairing(self) -> [f64; POLYNOMIAL_DIMENSION] {
        [
            self.value_weight,
            self.value_weight * self.point[0] + self.derivative_weight[0],
            self.value_weight * self.point[1] + self.derivative_weight[1],
            self.value_weight * self.point[2] + self.derivative_weight[2],
        ]
    }
}

#[derive(Debug, Clone)]
struct DenseMatrix {
    rows: usize,
    columns: usize,
    values: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
struct CoordinateNormalization {
    center: [f64; 3],
    length: f64,
}

impl CoordinateNormalization {
    fn from_functionals(functionals: &[Functional]) -> Self {
        let mut minimum = [f64::INFINITY; 3];
        let mut maximum = [f64::NEG_INFINITY; 3];
        for functional in functionals {
            for axis in 0..3 {
                minimum[axis] = minimum[axis].min(functional.point[axis]);
                maximum[axis] = maximum[axis].max(functional.point[axis]);
            }
        }
        let center = [
            0.5 * (minimum[0] + maximum[0]),
            0.5 * (minimum[1] + maximum[1]),
            0.5 * (minimum[2] + maximum[2]),
        ];
        let support_radius = functionals
            .iter()
            .map(|functional| {
                let displacement = subtract(functional.point, center);
                dot(displacement, displacement).sqrt()
            })
            .fold(0.0_f64, f64::max);
        let length = if support_radius > 0.0 {
            support_radius
        } else {
            1.0
        };
        Self { center, length }
    }

    fn polynomial_pairing(self, functional: Functional) -> [f64; POLYNOMIAL_DIMENSION] {
        let physical = functional.polynomial_pairing();
        [
            physical[0],
            (physical[1] - self.center[0] * physical[0]) / self.length,
            (physical[2] - self.center[1] * physical[0]) / self.length,
            (physical[3] - self.center[2] * physical[0]) / self.length,
        ]
    }

    fn to_standard_polynomial(
        self,
        physical: [f64; POLYNOMIAL_DIMENSION],
    ) -> [f64; POLYNOMIAL_DIMENSION] {
        [
            physical[0]
                + physical[1] * self.center[0]
                + physical[2] * self.center[1]
                + physical[3] * self.center[2],
            physical[1] * self.length,
            physical[2] * self.length,
            physical[3] * self.length,
        ]
    }

    fn to_physical_polynomial(
        self,
        standard: [f64; POLYNOMIAL_DIMENSION],
    ) -> [f64; POLYNOMIAL_DIMENSION] {
        let linear = [
            standard[1] / self.length,
            standard[2] / self.length,
            standard[3] / self.length,
        ];
        [
            standard[0]
                - linear[0] * self.center[0]
                - linear[1] * self.center[1]
                - linear[2] * self.center[2],
            linear[0],
            linear[1],
            linear[2],
        ]
    }
}

#[derive(Debug, Clone)]
struct CanonicalProblem {
    kernel: DenseMatrix,
    polynomial: DenseMatrix,
    normalization: CoordinateNormalization,
    truth_coefficients: Vec<f64>,
    truth_polynomial: [f64; POLYNOMIAL_DIMENSION],
    truth_latent: f64,
    truth_values: Vec<f64>,
}

impl DenseMatrix {
    fn from_fn(rows: usize, columns: usize, mut value: impl FnMut(usize, usize) -> f64) -> Self {
        let mut values = Vec::with_capacity(rows * columns);
        for row in 0..rows {
            for column in 0..columns {
                values.push(value(row, column));
            }
        }
        Self {
            rows,
            columns,
            values,
        }
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.columns + column]
    }

    fn set(&mut self, row: usize, column: usize, value: f64) {
        self.values[row * self.columns + column] = value;
    }

    fn multiply_vector(&self, vector: &[f64]) -> Vec<f64> {
        assert_eq!(self.columns, vector.len());
        (0..self.rows)
            .map(|row| {
                (0..self.columns)
                    .map(|column| self.get(row, column) * vector[column])
                    .sum()
            })
            .collect()
    }

    fn to_faer(&self) -> Mat<f64> {
        Mat::from_fn(self.rows, self.columns, |row, column| self.get(row, column))
    }

    fn to_csc(&self) -> CscMatrix<f64> {
        let rows = self
            .values
            .chunks(self.columns)
            .map(|row| row.to_vec())
            .collect::<Vec<_>>();
        CscMatrix::from(&rows)
    }
}

#[derive(Debug, Clone)]
struct HouseholderNullSpace {
    factor: Qr<f64>,
    ambient_dimension: usize,
    polynomial_rank: usize,
}

impl HouseholderNullSpace {
    fn new(polynomial: &DenseMatrix, polynomial_rank: usize) -> Self {
        Self {
            factor: polynomial.to_faer().qr(),
            ambient_dimension: polynomial.rows,
            polynomial_rank,
        }
    }

    fn reduced_dimension(&self) -> usize {
        self.ambient_dimension - self.polynomial_rank
    }

    fn expand(&self, reduced: &[f64]) -> Vec<f64> {
        assert_eq!(reduced.len(), self.reduced_dimension());
        let mut embedded = Mat::<f64>::zeros(self.ambient_dimension, 1);
        for (index, value) in reduced.iter().enumerate() {
            embedded[(self.polynomial_rank + index, 0)] = *value;
        }
        let mut memory = MemBuffer::new(
            apply_block_householder_sequence_on_the_left_in_place_scratch::<f64>(
                self.ambient_dimension,
                self.factor.Q_coeff().nrows(),
                1,
            ),
        );
        apply_block_householder_sequence_on_the_left_in_place_with_conj(
            self.factor.Q_basis(),
            self.factor.Q_coeff(),
            Conj::No,
            embedded.as_mut(),
            get_global_parallelism(),
            MemStack::new(&mut memory),
        );
        (0..self.ambient_dimension)
            .map(|row| embedded[(row, 0)])
            .collect()
    }

    fn project(&self, ambient: &[f64]) -> Vec<f64> {
        assert_eq!(ambient.len(), self.ambient_dimension);
        let mut transformed = Mat::<f64>::from_fn(self.ambient_dimension, 1, |row, _| ambient[row]);
        let mut memory = MemBuffer::new(
            apply_block_householder_sequence_transpose_on_the_left_in_place_scratch::<f64>(
                self.ambient_dimension,
                self.factor.Q_coeff().nrows(),
                1,
            ),
        );
        apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj(
            self.factor.Q_basis(),
            self.factor.Q_coeff(),
            Conj::No,
            transformed.as_mut(),
            get_global_parallelism(),
            MemStack::new(&mut memory),
        );
        (self.polynomial_rank..self.ambient_dimension)
            .map(|row| transformed[(row, 0)])
            .collect()
    }
}

pub fn run_manufactured_experiment() -> Result<ExperimentEvidence, ExperimentError> {
    let functionals = manufactured_functionals();
    let normalization = CoordinateNormalization::from_functionals(&functionals);
    let polynomial = assemble_polynomial_pairing(&functionals, normalization);
    let singular_values = polynomial
        .to_faer()
        .singular_values()
        .map_err(|error| ExperimentError::new(format!("polynomial SVD failed: {error:?}")))?;
    let largest = singular_values.first().copied().unwrap_or(0.0);
    let threshold =
        4096.0 * f64::EPSILON * functionals.len().max(POLYNOMIAL_DIMENSION) as f64 * largest;
    let polynomial_rank = singular_values
        .iter()
        .filter(|singular_value| **singular_value >= threshold)
        .count();
    if polynomial_rank != POLYNOMIAL_DIMENSION {
        return Err(ExperimentError::new(format!(
            "expected rank(P)=4 before solving, observed {polynomial_rank}"
        )));
    }
    let kernel = assemble_cubic_pairing(&functionals);
    let null_space = HouseholderNullSpace::new(&polynomial, polynomial_rank);
    let (mut reduced, null_space_defect) =
        materialize_reduced_pairing(&kernel, &polynomial, &null_space);
    let reduced_symmetry_defect = normalized_symmetry_defect(&reduced);
    let symmetry_defect_limit = 256.0 * f64::EPSILON * reduced.rows.max(reduced.columns) as f64;
    if null_space_defect > 1.0e-10 {
        return Err(ExperimentError::new(format!(
            "Householder null-space side-condition defect {null_space_defect:e}"
        )));
    }
    if reduced_symmetry_defect > symmetry_defect_limit {
        return Err(ExperimentError::new(format!(
            "reduced symmetry defect {reduced_symmetry_defect:e} exceeds {symmetry_defect_limit:e}"
        )));
    }
    symmetrize(&mut reduced);
    let reduced_eigenvalues = reduced
        .to_faer()
        .self_adjoint_eigenvalues(Side::Lower)
        .map_err(|error| ExperimentError::new(format!("reduced eigensolve failed: {error:?}")))?;
    let reduced_smallest_eigenvalue = reduced_eigenvalues[0];
    reduced
        .to_faer()
        .llt(Side::Lower)
        .map_err(|_| ExperimentError::new("reduced Cubic pairing is not strictly positive"))?;
    let affine_reproduction_error = affine_reproduction_error(&kernel, &polynomial);
    let problem = manufacture_canonical_problem(kernel, polynomial, normalization);
    let equality = solve_equality(&problem)?;
    let qp = solve_qp(&problem, &null_space, &reduced)?;
    let socp = solve_socp(&problem, &null_space, &reduced)?;
    let cross_route_observable_error =
        canonical_observable_difference(&equality.recovered, &qp.recovered).max(
            canonical_observable_difference(&equality.recovered, &socp.recovered),
        );

    Ok(ExperimentEvidence {
        cpd: CpdEvidence {
            functional_count: functionals.len(),
            polynomial_dimension: POLYNOMIAL_DIMENSION,
            polynomial_rank,
            singular_values,
            null_space_defect,
            reduced_symmetry_defect,
            symmetry_defect_limit,
            reduced_smallest_eigenvalue,
            affine_reproduction_error,
        },
        equality,
        qp,
        socp,
        cross_route_observable_error,
    })
}

pub fn run_counterexamples() -> Result<CounterexampleEvidence, ExperimentError> {
    let mut flattened_functionals = manufactured_functionals();
    for functional in &mut flattened_functionals {
        functional.point[2] = 0.0;
        functional.derivative_weight[2] = 0.0;
    }
    let flattened_normalization = CoordinateNormalization::from_functionals(&flattened_functionals);
    let flattened_polynomial =
        assemble_polynomial_pairing(&flattened_functionals, flattened_normalization);
    let flattened_rank = numerical_rank(&flattened_polynomial)?;
    let rank_deficient_polynomial = FailureEvidence {
        kind: FailureKind::PolynomialRankDeficient,
        solver_invoked: false,
        backend_contract_passed: false,
        hidden_regularization_applied: false,
        detected_violation: (POLYNOMIAL_DIMENSION - flattened_rank) as f64,
    };

    let functionals = manufactured_functionals();
    let normalization = CoordinateNormalization::from_functionals(&functionals);
    let polynomial = assemble_polynomial_pairing(&functionals, normalization);
    let polynomial_rank = numerical_rank(&polynomial)?;
    let kernel = assemble_cubic_pairing(&functionals);
    let null_space = HouseholderNullSpace::new(&polynomial, polynomial_rank);
    let mut first_reduced_unit = vec![0.0; null_space.reduced_dimension()];
    first_reduced_unit[0] = 1.0;
    let first_null_vector = null_space.expand(&first_reduced_unit);
    let first_energy = dot_product(
        &first_null_vector,
        &kernel.multiply_vector(&first_null_vector),
    );
    let damaged_kernel = DenseMatrix::from_fn(kernel.rows, kernel.columns, |row, column| {
        kernel.get(row, column)
            - 2.0 * first_energy * first_null_vector[row] * first_null_vector[column]
    });
    let (mut damaged_reduced, _) =
        materialize_reduced_pairing(&damaged_kernel, &polynomial, &null_space);
    symmetrize(&mut damaged_reduced);
    let smallest_damaged_eigenvalue = damaged_reduced
        .to_faer()
        .self_adjoint_eigenvalues(Side::Lower)
        .map_err(|error| {
            ExperimentError::new(format!("damaged reduced eigensolve failed: {error:?}"))
        })?[0];
    let nonpositive_reduced_pairing = FailureEvidence {
        kind: FailureKind::ReducedPairingNotPositive,
        solver_invoked: false,
        backend_contract_passed: false,
        hidden_regularization_applied: false,
        detected_violation: (-smallest_damaged_eigenvalue).max(0.0),
    };

    let (mut reduced, _) = materialize_reduced_pairing(&kernel, &polynomial, &null_space);
    symmetrize(&mut reduced);
    let problem = manufacture_canonical_problem(kernel, polynomial, normalization);
    let qp = solve_qp(&problem, &null_space, &reduced)?;
    let backend_contract_passed = qp.scaled.primal <= 1.0e-8
        && qp.scaled.dual <= 1.0e-8
        && qp.scaled.stationarity <= 1.0e-8
        && qp.scaled.complementarity <= 1.0e-8
        && qp.scaled.relative_gap <= 1.0e-8;
    let mut corrupted_coefficients = qp.recovered.field_coefficients.clone();
    corrupted_coefficients[0] += 0.01;
    let standard_polynomial =
        normalization.to_standard_polynomial(qp.recovered.polynomial_coefficients);
    let mut corrupted_primal = corrupted_coefficients.clone();
    corrupted_primal.extend_from_slice(&standard_polynomial);
    corrupted_primal.push(qp.recovered.semantic_latents[0]);
    let corrupted = recover_canonical(&problem, &corrupted_primal, qp.recovered.slacks.clone());
    let projected = null_space.project(&corrupted_coefficients);
    let round_trip = null_space.expand(&projected);
    let broken_map_round_trip = corrupted_coefficients
        .iter()
        .zip(round_trip)
        .map(|(actual, recovered)| (actual - recovered).abs() / actual.abs().max(1.0))
        .fold(0.0_f64, f64::max);
    let broken_recovery = FailureEvidence {
        kind: FailureKind::RecoveryVerification,
        solver_invoked: true,
        backend_contract_passed,
        hidden_regularization_applied: false,
        detected_violation: corrupted
            .hard_violation
            .max(corrupted.side_condition_violation)
            .max(broken_map_round_trip),
    };

    Ok(CounterexampleEvidence {
        rank_deficient_polynomial,
        nonpositive_reduced_pairing,
        broken_recovery,
    })
}

fn manufactured_functionals() -> Vec<Functional> {
    vec![
        Functional {
            point: [-1.0, -1.0, -1.0],
            value_weight: 1.0,
            derivative_weight: [0.0; 3],
        },
        Functional {
            point: [1.0, -1.0, -1.0],
            value_weight: 1.0,
            derivative_weight: [0.0; 3],
        },
        Functional {
            point: [-1.0, 1.0, -1.0],
            value_weight: 1.0,
            derivative_weight: [0.0; 3],
        },
        Functional {
            point: [-1.0, -1.0, 1.0],
            value_weight: 1.0,
            derivative_weight: [0.0; 3],
        },
        Functional {
            point: [1.0, 1.0, 1.0],
            value_weight: 1.0,
            derivative_weight: [0.0; 3],
        },
        Functional {
            point: [0.25, -0.5, 0.75],
            value_weight: 0.0,
            derivative_weight: [1.0, 0.0, 0.0],
        },
        Functional {
            point: [-0.75, 0.25, 0.5],
            value_weight: 0.0,
            derivative_weight: [0.0, 1.0, 0.0],
        },
        Functional {
            point: [0.5, 0.75, -0.25],
            value_weight: 0.0,
            derivative_weight: [0.0, 0.0, 1.0],
        },
        Functional {
            point: [0.0, 0.0, 0.0],
            value_weight: 1.0,
            derivative_weight: [0.5, -0.25, 0.125],
        },
        Functional {
            point: [-0.5, 0.625, -0.75],
            value_weight: 0.0,
            derivative_weight: [1.0, 1.0, 1.0],
        },
    ]
}

fn assemble_polynomial_pairing(
    functionals: &[Functional],
    normalization: CoordinateNormalization,
) -> DenseMatrix {
    DenseMatrix::from_fn(functionals.len(), POLYNOMIAL_DIMENSION, |row, column| {
        normalization.polynomial_pairing(functionals[row])[column]
    })
}

fn numerical_rank(matrix: &DenseMatrix) -> Result<usize, ExperimentError> {
    let singular_values = matrix
        .to_faer()
        .singular_values()
        .map_err(|error| ExperimentError::new(format!("rank SVD failed: {error:?}")))?;
    let largest = singular_values.first().copied().unwrap_or(0.0);
    let threshold = 4096.0 * f64::EPSILON * matrix.rows.max(matrix.columns) as f64 * largest;
    Ok(singular_values
        .iter()
        .filter(|singular_value| **singular_value >= threshold)
        .count())
}

fn assemble_cubic_pairing(functionals: &[Functional]) -> DenseMatrix {
    let mut pairing = DenseMatrix::from_fn(functionals.len(), functionals.len(), |_, _| 0.0);
    for row in 0..functionals.len() {
        for column in row..functionals.len() {
            let value = cubic_functional_pairing(functionals[row], functionals[column]);
            pairing.set(row, column, value);
            pairing.set(column, row, value);
        }
    }
    pairing
}

fn cubic_functional_pairing(left: Functional, right: Functional) -> f64 {
    let displacement = subtract(left.point, right.point);
    let radius = dot(displacement, displacement).sqrt();
    if radius == 0.0 {
        return 0.0;
    }

    let gradient_x = scale(displacement, 3.0 * radius);
    let gradient_y = scale(gradient_x, -1.0);
    let mixed_derivative = DenseMatrix::from_fn(3, 3, |row, column| {
        -3.0 * (if row == column { radius } else { 0.0 }
            + displacement[row] * displacement[column] / radius)
    });

    left.value_weight * right.value_weight * radius.powi(3)
        + left.value_weight * dot(right.derivative_weight, gradient_y)
        + right.value_weight * dot(left.derivative_weight, gradient_x)
        + dot(
            left.derivative_weight,
            multiply_three_by_three(&mixed_derivative, right.derivative_weight),
        )
}

fn materialize_reduced_pairing(
    kernel: &DenseMatrix,
    polynomial: &DenseMatrix,
    null_space: &HouseholderNullSpace,
) -> (DenseMatrix, f64) {
    let reduced_dimension = null_space.reduced_dimension();
    let mut null_space_defect = 0.0_f64;
    let reduced = DenseMatrix::from_fn(reduced_dimension, reduced_dimension, |row, column| {
        let mut unit = vec![0.0; reduced_dimension];
        unit[column] = 1.0;
        let null_column = null_space.expand(&unit);
        for polynomial_column in 0..POLYNOMIAL_DIMENSION {
            let side_value = (0..polynomial.rows)
                .map(|index| polynomial.get(index, polynomial_column) * null_column[index])
                .sum::<f64>();
            null_space_defect = null_space_defect.max(side_value.abs());
        }
        let projected = null_space.project(&kernel.multiply_vector(&null_column));
        projected[row]
    });
    (reduced, null_space_defect)
}

fn normalized_symmetry_defect(matrix: &DenseMatrix) -> f64 {
    let mut defect = 0.0_f64;
    let mut scale = 1.0_f64;
    for row in 0..matrix.rows {
        for column in 0..matrix.columns {
            defect = defect.max((matrix.get(row, column) - matrix.get(column, row)).abs());
            scale = scale.max(matrix.get(row, column).abs());
        }
    }
    defect / scale
}

fn symmetrize(matrix: &mut DenseMatrix) {
    for row in 0..matrix.rows {
        for column in (row + 1)..matrix.columns {
            let average = 0.5 * (matrix.get(row, column) + matrix.get(column, row));
            matrix.set(row, column, average);
            matrix.set(column, row, average);
        }
    }
}

fn affine_reproduction_error(kernel: &DenseMatrix, polynomial: &DenseMatrix) -> f64 {
    let augmented_dimension = kernel.rows + POLYNOMIAL_DIMENSION;
    let augmented =
        Mat::<f64>::from_fn(
            augmented_dimension,
            augmented_dimension,
            |row, column| match (row < kernel.rows, column < kernel.columns) {
                (true, true) => kernel.get(row, column),
                (true, false) => polynomial.get(row, column - kernel.columns),
                (false, true) => polynomial.get(column, row - kernel.rows),
                (false, false) => 0.0,
            },
        );
    let factor = augmented.lblt(Side::Lower);
    let mut maximum_error = 0.0_f64;
    for basis in 0..POLYNOMIAL_DIMENSION {
        let rhs = Mat::<f64>::from_fn(augmented_dimension, 1, |row, _| {
            if row < polynomial.rows {
                polynomial.get(row, basis)
            } else {
                0.0
            }
        });
        let solution = factor.solve(&rhs);
        for row in 0..kernel.rows {
            maximum_error = maximum_error.max(solution[(row, 0)].abs());
        }
        for coefficient in 0..POLYNOMIAL_DIMENSION {
            let expected = if coefficient == basis { 1.0 } else { 0.0 };
            maximum_error =
                maximum_error.max((solution[(kernel.rows + coefficient, 0)] - expected).abs());
        }
    }
    maximum_error
}

fn manufacture_canonical_problem(
    kernel: DenseMatrix,
    polynomial: DenseMatrix,
    normalization: CoordinateNormalization,
) -> CanonicalProblem {
    let truth_coefficients = vec![
        0.195, -0.105, -0.17, -0.10, 0.10, -0.07, 0.12, -0.05, 0.08, -0.04,
    ];
    let kernel_values = kernel.multiply_vector(&truth_coefficients);
    let physical_polynomial = [0.6, 0.5 * (kernel_values[0] - kernel_values[1]), 0.4, 0.15];
    let standard_polynomial = normalization.to_standard_polynomial(physical_polynomial);
    let polynomial_values = polynomial.multiply_vector(&standard_polynomial);
    let truth_values = kernel_values
        .iter()
        .zip(polynomial_values)
        .map(|(kernel_value, polynomial_value)| kernel_value + polynomial_value)
        .collect::<Vec<_>>();
    let truth_latent = 0.5 * (truth_values[0] + truth_values[1]);

    CanonicalProblem {
        kernel,
        polynomial,
        normalization,
        truth_coefficients,
        truth_polynomial: physical_polynomial,
        truth_latent,
        truth_values,
    }
}

fn solve_equality(problem: &CanonicalProblem) -> Result<EqualityEvidence, ExperimentError> {
    let (hessian, constraints, constraint_rhs) = equality_primal_form(problem);
    let primal_dimension = hessian.rows;
    let constraint_dimension = constraints.rows;
    let kkt_dimension = primal_dimension + constraint_dimension;
    let kkt = DenseMatrix::from_fn(kkt_dimension, kkt_dimension, |row, column| {
        match (row < primal_dimension, column < primal_dimension) {
            (true, true) => hessian.get(row, column),
            (true, false) => constraints.get(column - primal_dimension, row),
            (false, true) => constraints.get(row - primal_dimension, column),
            (false, false) => 0.0,
        }
    });
    let rhs = (0..kkt_dimension)
        .map(|row| {
            if row < primal_dimension {
                0.0
            } else {
                constraint_rhs[row - primal_dimension]
            }
        })
        .collect::<Vec<_>>();
    let scaled = scale_symmetric_kkt(&kkt, &rhs)?;
    let factor = scaled.matrix.to_faer().lblt(Side::Lower);
    let scaled_rhs = Mat::<f64>::from_fn(kkt_dimension, 1, |row, _| scaled.rhs[row]);
    let scaled_solution_matrix = factor.solve(&scaled_rhs);
    let scaled_solution = (0..kkt_dimension)
        .map(|row| scaled_solution_matrix[(row, 0)])
        .collect::<Vec<_>>();
    let solution = scaled_solution
        .iter()
        .zip(&scaled.factors)
        .map(|(value, factor)| value * factor)
        .collect::<Vec<_>>();
    let scaled_backward_error =
        normalized_backward_error(&scaled.matrix, &scaled_solution, &scaled.rhs);
    let physical_backward_error = normalized_backward_error(&kkt, &solution, &rhs);
    let normalized_backward_error = scaled_backward_error.max(physical_backward_error);
    let inertia = inertia_from_lblt(&factor, kkt_dimension);
    let expected_inertia = InertiaEvidence {
        positive: primal_dimension,
        negative: constraint_dimension,
        zero: 0,
    };

    let primal = &solution[..primal_dimension];
    let recovered = recover_canonical(problem, primal, Vec::new());
    let scaling_round_trip_error = scaling_round_trip_error(
        primal,
        &scaled.factors[..primal_dimension],
        problem.normalization,
        recovered.polynomial_coefficients,
    );
    let manufactured_truth_error = manufactured_truth_error(problem, &recovered);

    Ok(EqualityEvidence {
        inertia,
        expected_inertia,
        normalized_backward_error,
        scaling_round_trip_error,
        manufactured_truth_error,
        recovered,
    })
}

fn equality_primal_form(problem: &CanonicalProblem) -> (DenseMatrix, DenseMatrix, Vec<f64>) {
    let coefficient_count = problem.kernel.rows;
    let primal_dimension = coefficient_count + POLYNOMIAL_DIMENSION + 1;
    let constraint_dimension = POLYNOMIAL_DIMENSION + 3 + (coefficient_count - 2);
    let hessian = DenseMatrix::from_fn(primal_dimension, primal_dimension, |row, column| {
        if row < coefficient_count && column < coefficient_count {
            problem.kernel.get(row, column)
        } else {
            0.0
        }
    });
    let latent_column = coefficient_count + POLYNOMIAL_DIMENSION;
    let constraints =
        DenseMatrix::from_fn(constraint_dimension, primal_dimension, |row, column| {
            if row < POLYNOMIAL_DIMENSION {
                if column < coefficient_count {
                    problem.polynomial.get(column, row)
                } else {
                    0.0
                }
            } else if row == POLYNOMIAL_DIMENSION || row == POLYNOMIAL_DIMENSION + 1 {
                let functional = row - POLYNOMIAL_DIMENSION;
                if column < coefficient_count {
                    problem.kernel.get(functional, column)
                } else if column < latent_column {
                    problem
                        .polynomial
                        .get(functional, column - coefficient_count)
                } else {
                    -1.0
                }
            } else if row == POLYNOMIAL_DIMENSION + 2 {
                if column == latent_column { 1.0 } else { 0.0 }
            } else {
                let functional = row - (POLYNOMIAL_DIMENSION + 1);
                if column < coefficient_count {
                    problem.kernel.get(functional, column)
                } else if column < latent_column {
                    problem
                        .polynomial
                        .get(functional, column - coefficient_count)
                } else {
                    0.0
                }
            }
        });
    let constraint_rhs = (0..constraint_dimension)
        .map(|row| match row.cmp(&(POLYNOMIAL_DIMENSION + 2)) {
            Ordering::Less => 0.0,
            Ordering::Equal => problem.truth_latent,
            Ordering::Greater => {
                let functional = row - (POLYNOMIAL_DIMENSION + 1);
                problem.truth_values[functional]
            }
        })
        .collect();
    (hessian, constraints, constraint_rhs)
}

fn solve_qp(
    problem: &CanonicalProblem,
    null_space: &HouseholderNullSpace,
    reduced_hessian: &DenseMatrix,
) -> Result<ConvexRouteEvidence, ExperimentError> {
    let (equality_rows, equality_rhs) = reduced_equalities(problem, null_space);
    let reduced_dimension = null_space.reduced_dimension();
    let variable_count = reduced_dimension + POLYNOMIAL_DIMENSION + 1;
    let hessian = DenseMatrix::from_fn(variable_count, variable_count, |row, column| {
        if row < reduced_dimension && column < reduced_dimension {
            reduced_hessian.get(row, column)
        } else {
            0.0
        }
    });
    let linear_objective = vec![0.0; variable_count];
    let affine_bound_row = reduced_functional_row(problem, null_space, 4);
    let constraint_matrix =
        DenseMatrix::from_fn(equality_rows.rows + 1, variable_count, |row, column| {
            if row < equality_rows.rows {
                equality_rows.get(row, column)
            } else {
                affine_bound_row[column]
            }
        });
    let mut constraint_rhs = equality_rhs;
    constraint_rhs.push(problem.truth_values[4] + 0.75);
    let cone_blocks = vec![
        ConeBlock::Zero(equality_rows.rows),
        ConeBlock::Nonnegative(1),
    ];
    let scaled = scale_conic_form(
        &hessian,
        &linear_objective,
        &constraint_matrix,
        &constraint_rhs,
        &cone_blocks,
    )?;
    let p = scaled.hessian.to_csc().to_triu();
    let a = scaled.constraints.to_csc();
    let cones = [ZeroConeT(equality_rows.rows), NonnegativeConeT(1)];
    let settings = clarabel_settings();
    let mut solver = DefaultSolver::new(
        &p,
        &scaled.linear_objective,
        &a,
        &scaled.constraint_rhs,
        &cones,
        settings,
    )
    .map_err(|error| ExperimentError::new(format!("Clarabel QP setup failed: {error}")))?;
    solver.solve();
    if !matches!(
        solver.solution.status,
        SolverStatus::Solved | SolverStatus::AlmostSolved
    ) {
        return Err(ExperimentError::new(format!(
            "Clarabel QP returned {}",
            solver.solution.status
        )));
    }
    let scaled_primal = solver.solution.x.clone();
    let scaled_dual = solver.solution.z.clone();
    let scaled_slack = solver.solution.s.clone();
    let residuals = scaled_convex_residuals(
        &scaled.hessian,
        &scaled.linear_objective,
        &scaled.constraints,
        &scaled.constraint_rhs,
        &cone_blocks,
        &scaled_primal,
        &scaled_dual,
        &scaled_slack,
    );
    let physical_reduced = scaled_primal
        .iter()
        .zip(&scaled.variable_factors)
        .map(|(value, factor)| value * factor)
        .collect::<Vec<_>>();
    let physical_slacks = scaled_slack
        .iter()
        .zip(&scaled.constraint_factors)
        .map(|(value, factor)| value / factor)
        .collect::<Vec<_>>();
    let coefficients = null_space.expand(&physical_reduced[..reduced_dimension]);
    let mut full_primal = coefficients.clone();
    full_primal.extend_from_slice(&physical_reduced[reduced_dimension..]);
    let mut recovered = recover_canonical(
        problem,
        &full_primal,
        vec![physical_slacks[equality_rows.rows]],
    );
    let affine_violation =
        (recovered.functional_values[4] - (problem.truth_values[4] + 0.75)).max(0.0);
    recovered.hard_violation = recovered.hard_violation.max(affine_violation);
    let recovery_round_trip_error = conic_recovery_round_trip_error(
        null_space,
        &physical_reduced,
        &scaled_primal,
        &physical_slacks,
        &scaled_slack,
        &scaled.variable_factors,
        &scaled.constraint_factors,
        problem.normalization,
        recovered.polynomial_coefficients,
    );
    let manufactured_truth_error = manufactured_truth_error(problem, &recovered);

    Ok(ConvexRouteEvidence {
        problem_class: "reduced_qp",
        scaled: residuals,
        recovery_round_trip_error,
        manufactured_truth_error,
        recovered,
    })
}

fn solve_socp(
    problem: &CanonicalProblem,
    null_space: &HouseholderNullSpace,
    reduced_hessian: &DenseMatrix,
) -> Result<ConvexRouteEvidence, ExperimentError> {
    let (equality_rows, equality_rhs) = reduced_equalities(problem, null_space);
    let reduced_dimension = null_space.reduced_dimension();
    let variable_count = reduced_dimension + POLYNOMIAL_DIMENSION + 1;
    let hessian = DenseMatrix::from_fn(variable_count, variable_count, |row, column| {
        if row < reduced_dimension && column < reduced_dimension {
            reduced_hessian.get(row, column)
        } else {
            0.0
        }
    });
    let linear_objective = vec![0.0; variable_count];
    let cone_rows = [
        reduced_functional_row(problem, null_space, 4),
        reduced_functional_row(problem, null_space, 5),
        reduced_functional_row(problem, null_space, 6),
    ];
    let constraint_matrix =
        DenseMatrix::from_fn(equality_rows.rows + 3, variable_count, |row, column| {
            if row < equality_rows.rows {
                equality_rows.get(row, column)
            } else {
                cone_rows[row - equality_rows.rows][column]
            }
        });
    let mut constraint_rhs = equality_rhs;
    constraint_rhs.extend([
        problem.truth_values[4] + 1.0,
        problem.truth_values[5] + 0.2,
        problem.truth_values[6] - 0.1,
    ]);
    let cone_blocks = vec![
        ConeBlock::Zero(equality_rows.rows),
        ConeBlock::SecondOrder(3),
    ];
    let scaled = scale_conic_form(
        &hessian,
        &linear_objective,
        &constraint_matrix,
        &constraint_rhs,
        &cone_blocks,
    )?;
    let p = scaled.hessian.to_csc().to_triu();
    let a = scaled.constraints.to_csc();
    let cones = [ZeroConeT(equality_rows.rows), SecondOrderConeT(3)];
    let settings = clarabel_settings();
    let mut solver = DefaultSolver::new(
        &p,
        &scaled.linear_objective,
        &a,
        &scaled.constraint_rhs,
        &cones,
        settings,
    )
    .map_err(|error| ExperimentError::new(format!("Clarabel SOCP setup failed: {error}")))?;
    solver.solve();
    if !matches!(
        solver.solution.status,
        SolverStatus::Solved | SolverStatus::AlmostSolved
    ) {
        return Err(ExperimentError::new(format!(
            "Clarabel SOCP returned {}",
            solver.solution.status
        )));
    }
    let scaled_primal = solver.solution.x.clone();
    let scaled_dual = solver.solution.z.clone();
    let scaled_slack = solver.solution.s.clone();
    let residuals = scaled_convex_residuals(
        &scaled.hessian,
        &scaled.linear_objective,
        &scaled.constraints,
        &scaled.constraint_rhs,
        &cone_blocks,
        &scaled_primal,
        &scaled_dual,
        &scaled_slack,
    );
    let physical_reduced = scaled_primal
        .iter()
        .zip(&scaled.variable_factors)
        .map(|(value, factor)| value * factor)
        .collect::<Vec<_>>();
    let physical_slacks = scaled_slack
        .iter()
        .zip(&scaled.constraint_factors)
        .map(|(value, factor)| value / factor)
        .collect::<Vec<_>>();
    let coefficients = null_space.expand(&physical_reduced[..reduced_dimension]);
    let mut full_primal = coefficients.clone();
    full_primal.extend_from_slice(&physical_reduced[reduced_dimension..]);
    let canonical_slacks = physical_slacks[equality_rows.rows..].to_vec();
    let mut recovered = recover_canonical(problem, &full_primal, canonical_slacks.clone());
    let cone_violation =
        (canonical_slacks[1].hypot(canonical_slacks[2]) - canonical_slacks[0]).max(0.0);
    recovered.hard_violation = recovered.hard_violation.max(cone_violation);
    let recovery_round_trip_error = conic_recovery_round_trip_error(
        null_space,
        &physical_reduced,
        &scaled_primal,
        &physical_slacks,
        &scaled_slack,
        &scaled.variable_factors,
        &scaled.constraint_factors,
        problem.normalization,
        recovered.polynomial_coefficients,
    );
    let manufactured_truth_error = manufactured_truth_error(problem, &recovered);

    Ok(ConvexRouteEvidence {
        problem_class: "reduced_socp",
        scaled: residuals,
        recovery_round_trip_error,
        manufactured_truth_error,
        recovered,
    })
}

fn reduced_equalities(
    problem: &CanonicalProblem,
    null_space: &HouseholderNullSpace,
) -> (DenseMatrix, Vec<f64>) {
    let (_, full_constraints, full_rhs) = equality_primal_form(problem);
    let reduced_dimension = null_space.reduced_dimension();
    let reduced_variable_count = reduced_dimension + POLYNOMIAL_DIMENSION + 1;
    let rows = full_constraints.rows - POLYNOMIAL_DIMENSION;
    let reduced = DenseMatrix::from_fn(rows, reduced_variable_count, |row, column| {
        let full_row = row + POLYNOMIAL_DIMENSION;
        if column < reduced_dimension {
            let ambient = (0..problem.kernel.rows)
                .map(|index| full_constraints.get(full_row, index))
                .collect::<Vec<_>>();
            null_space.project(&ambient)[column]
        } else {
            let full_column = problem.kernel.rows + (column - reduced_dimension);
            full_constraints.get(full_row, full_column)
        }
    });
    (reduced, full_rhs[POLYNOMIAL_DIMENSION..].to_vec())
}

fn reduced_functional_row(
    problem: &CanonicalProblem,
    null_space: &HouseholderNullSpace,
    functional: usize,
) -> Vec<f64> {
    let reduced_dimension = null_space.reduced_dimension();
    let mut row = null_space.project(
        &(0..problem.kernel.columns)
            .map(|column| problem.kernel.get(functional, column))
            .collect::<Vec<_>>(),
    );
    row.extend((0..POLYNOMIAL_DIMENSION).map(|column| problem.polynomial.get(functional, column)));
    row.push(0.0);
    assert_eq!(row.len(), reduced_dimension + POLYNOMIAL_DIMENSION + 1);
    row
}

#[derive(Debug, Clone, Copy)]
enum ConeBlock {
    Zero(usize),
    Nonnegative(usize),
    SecondOrder(usize),
}

impl ConeBlock {
    fn dimension(self) -> usize {
        match self {
            Self::Zero(dimension) | Self::Nonnegative(dimension) | Self::SecondOrder(dimension) => {
                dimension
            }
        }
    }

    fn uses_common_scaling(self) -> bool {
        matches!(self, Self::SecondOrder(_))
    }
}

#[derive(Debug, Clone)]
struct ScaledConicForm {
    hessian: DenseMatrix,
    linear_objective: Vec<f64>,
    constraints: DenseMatrix,
    constraint_rhs: Vec<f64>,
    variable_factors: Vec<f64>,
    constraint_factors: Vec<f64>,
}

fn scale_conic_form(
    hessian: &DenseMatrix,
    linear_objective: &[f64],
    constraints: &DenseMatrix,
    constraint_rhs: &[f64],
    cone_blocks: &[ConeBlock],
) -> Result<ScaledConicForm, ExperimentError> {
    let mut scaled_hessian = hessian.clone();
    let mut scaled_objective = linear_objective.to_vec();
    let mut scaled_constraints = constraints.clone();
    let mut scaled_rhs = constraint_rhs.to_vec();
    let mut variable_factors = vec![1.0; hessian.columns];
    let mut constraint_factors = vec![1.0; constraints.rows];
    for _ in 0..8 {
        let variable_round = (0..hessian.columns)
            .map(|column| {
                let hessian_norm = (0..hessian.rows)
                    .map(|row| scaled_hessian.get(row, column).abs())
                    .fold(0.0_f64, f64::max);
                let constraint_norm = (0..constraints.rows)
                    .map(|row| scaled_constraints.get(row, column).abs())
                    .fold(0.0_f64, f64::max);
                bounded_power_of_two_factor(
                    hessian_norm.max(constraint_norm),
                    variable_factors[column],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        for row in 0..scaled_hessian.rows {
            for column in 0..scaled_hessian.columns {
                scaled_hessian.set(
                    row,
                    column,
                    scaled_hessian.get(row, column) * variable_round[row] * variable_round[column],
                );
            }
        }
        for column in 0..scaled_objective.len() {
            scaled_objective[column] *= variable_round[column];
            variable_factors[column] *= variable_round[column];
        }
        for row in 0..scaled_constraints.rows {
            for (column, factor) in variable_round.iter().copied().enumerate() {
                scaled_constraints.set(row, column, scaled_constraints.get(row, column) * factor);
            }
        }

        let mut row_round = vec![1.0; scaled_constraints.rows];
        let mut start = 0;
        for block in cone_blocks {
            let end = start + block.dimension();
            if block.uses_common_scaling() {
                let mut norm = 0.0_f64;
                for row in start..end {
                    for column in 0..scaled_constraints.columns {
                        norm = norm.max(scaled_constraints.get(row, column).abs());
                    }
                }
                let cumulative = constraint_factors[start..end]
                    .iter()
                    .copied()
                    .fold(1.0_f64, f64::max);
                let factor = bounded_power_of_two_factor(norm, cumulative)?;
                row_round[start..end].fill(factor);
            } else {
                for row in start..end {
                    let norm = (0..scaled_constraints.columns)
                        .map(|column| scaled_constraints.get(row, column).abs())
                        .fold(0.0_f64, f64::max);
                    row_round[row] = bounded_power_of_two_factor(norm, constraint_factors[row])?;
                }
            }
            start = end;
        }
        for row in 0..scaled_constraints.rows {
            constraint_factors[row] *= row_round[row];
            scaled_rhs[row] *= row_round[row];
            for column in 0..scaled_constraints.columns {
                scaled_constraints.set(
                    row,
                    column,
                    scaled_constraints.get(row, column) * row_round[row],
                );
            }
        }
    }
    Ok(ScaledConicForm {
        hessian: scaled_hessian,
        linear_objective: scaled_objective,
        constraints: scaled_constraints,
        constraint_rhs: scaled_rhs,
        variable_factors,
        constraint_factors,
    })
}

fn clarabel_settings() -> DefaultSettings<f64> {
    DefaultSettings {
        verbose: false,
        max_threads: 1,
        direct_solve_method: "qdldl".into(),
        tol_gap_abs: 1.0e-10,
        tol_gap_rel: 1.0e-10,
        tol_feas: 1.0e-10,
        ..DefaultSettings::default()
    }
}

#[allow(clippy::too_many_arguments)]
fn scaled_convex_residuals(
    hessian: &DenseMatrix,
    linear_objective: &[f64],
    constraints: &DenseMatrix,
    constraint_rhs: &[f64],
    cone_blocks: &[ConeBlock],
    primal: &[f64],
    dual: &[f64],
    slack: &[f64],
) -> ScaledConvexResiduals {
    let affine = constraints.multiply_vector(primal);
    let equation_residual = affine
        .iter()
        .zip(slack)
        .zip(constraint_rhs)
        .map(|((affine, slack), rhs)| (affine + slack - rhs).abs())
        .fold(0.0_f64, f64::max);
    let primal_scale = affine
        .iter()
        .chain(slack)
        .chain(constraint_rhs)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let primal_cone_violation = cone_violation(slack, cone_blocks);
    let primal_residual =
        (equation_residual / primal_scale).max(primal_cone_violation / primal_scale);

    let hessian_product = hessian.multiply_vector(primal);
    let mut stationarity_vector = hessian_product
        .iter()
        .zip(linear_objective)
        .map(|(quadratic, linear)| quadratic + linear)
        .collect::<Vec<_>>();
    for (column, stationarity) in stationarity_vector.iter_mut().enumerate() {
        *stationarity += (0..constraints.rows)
            .map(|row| constraints.get(row, column) * dual[row])
            .sum::<f64>();
    }
    let stationarity_scale = hessian_product
        .iter()
        .chain(linear_objective)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let stationarity = stationarity_vector
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max)
        / stationarity_scale;
    let dual_residual = cone_violation(dual, cone_blocks)
        / dual.iter().map(|value| value.abs()).fold(1.0_f64, f64::max);
    let primal_objective =
        0.5 * dot_product(primal, &hessian_product) + dot_product(linear_objective, primal);
    let dual_objective =
        -0.5 * dot_product(primal, &hessian_product) - dot_product(constraint_rhs, dual);
    let complementarity =
        dot_product(slack, dual).abs() / (1.0 + primal_objective.abs().max(dual_objective.abs()));
    let relative_gap = (primal_objective - dual_objective).abs()
        / (1.0 + primal_objective.abs().max(dual_objective.abs()));
    ScaledConvexResiduals {
        primal: primal_residual,
        dual: dual_residual,
        stationarity,
        complementarity,
        relative_gap,
    }
}

fn cone_violation(vector: &[f64], cone_blocks: &[ConeBlock]) -> f64 {
    let mut violation = 0.0_f64;
    let mut start = 0;
    for block in cone_blocks {
        let end = start + block.dimension();
        match block {
            ConeBlock::Zero(_) => {}
            ConeBlock::Nonnegative(_) => {
                violation = vector[start..end]
                    .iter()
                    .map(|value| (-value).max(0.0))
                    .fold(violation, f64::max);
            }
            ConeBlock::SecondOrder(_) => {
                let tail_norm = vector[(start + 1)..end]
                    .iter()
                    .map(|value| value * value)
                    .sum::<f64>()
                    .sqrt();
                violation = violation.max((tail_norm - vector[start]).max(0.0));
            }
        }
        start = end;
    }
    violation
}

#[allow(clippy::too_many_arguments)]
fn conic_recovery_round_trip_error(
    null_space: &HouseholderNullSpace,
    physical_reduced: &[f64],
    scaled_primal: &[f64],
    physical_slacks: &[f64],
    scaled_slacks: &[f64],
    variable_factors: &[f64],
    constraint_factors: &[f64],
    normalization: CoordinateNormalization,
    physical_polynomial: [f64; POLYNOMIAL_DIMENSION],
) -> f64 {
    let variable_error = physical_reduced
        .iter()
        .zip(scaled_primal)
        .zip(variable_factors)
        .map(|((physical, scaled), factor)| {
            (physical - factor * scaled).abs() / physical.abs().max(1.0)
        })
        .fold(0.0_f64, f64::max);
    let slack_error = physical_slacks
        .iter()
        .zip(scaled_slacks)
        .zip(constraint_factors)
        .map(|((physical, scaled), factor)| {
            (physical - scaled / factor).abs() / physical.abs().max(1.0)
        })
        .fold(0.0_f64, f64::max);
    let reduced_dimension = null_space.reduced_dimension();
    let coefficients = null_space.expand(&physical_reduced[..reduced_dimension]);
    let projected = null_space.project(&coefficients);
    let reduction_error = projected
        .iter()
        .zip(&physical_reduced[..reduced_dimension])
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(0.0_f64, f64::max);
    let polynomial_round_trip = normalization
        .to_physical_polynomial(normalization.to_standard_polynomial(physical_polynomial));
    polynomial_round_trip
        .iter()
        .zip(physical_polynomial)
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(
            variable_error.max(slack_error).max(reduction_error),
            f64::max,
        )
}

#[derive(Debug, Clone)]
struct SymmetricScaling {
    matrix: DenseMatrix,
    rhs: Vec<f64>,
    factors: Vec<f64>,
}

fn scale_symmetric_kkt(
    matrix: &DenseMatrix,
    rhs: &[f64],
) -> Result<SymmetricScaling, ExperimentError> {
    let mut scaled_matrix = matrix.clone();
    let mut scaled_rhs = rhs.to_vec();
    let mut factors = vec![1.0; matrix.rows];
    for _ in 0..8 {
        let round_factors = (0..matrix.rows)
            .map(|row| {
                let norm = (0..matrix.columns)
                    .map(|column| scaled_matrix.get(row, column).abs())
                    .fold(0.0_f64, f64::max);
                bounded_power_of_two_factor(norm, factors[row])
            })
            .collect::<Result<Vec<_>, _>>()?;
        for row in 0..matrix.rows {
            factors[row] *= round_factors[row];
            scaled_rhs[row] *= round_factors[row];
        }
        for row in 0..matrix.rows {
            for column in 0..matrix.columns {
                scaled_matrix.set(
                    row,
                    column,
                    scaled_matrix.get(row, column) * round_factors[row] * round_factors[column],
                );
            }
        }
    }
    Ok(SymmetricScaling {
        matrix: scaled_matrix,
        rhs: scaled_rhs,
        factors,
    })
}

fn bounded_power_of_two_factor(norm: f64, cumulative: f64) -> Result<f64, ExperimentError> {
    if !norm.is_finite() || norm == 0.0 {
        return Err(ExperimentError::new(format!(
            "Ruiz scaling encountered invalid norm {norm}"
        )));
    }
    let proposed_exponent = (-0.5 * norm.log2()).round() as i32;
    let round_exponent = proposed_exponent.clamp(-8, 8);
    let cumulative_exponent = cumulative.log2().round() as i32;
    let bounded_exponent =
        round_exponent.clamp(-32 - cumulative_exponent, 32 - cumulative_exponent);
    Ok(2.0_f64.powi(bounded_exponent))
}

fn inertia_from_lblt(
    factor: &faer::linalg::solvers::Lblt<f64>,
    dimension: usize,
) -> InertiaEvidence {
    let diagonal = (0..dimension)
        .map(|index| factor.B_diag()[index])
        .collect::<Vec<_>>();
    let subdiagonal = (0..dimension)
        .map(|index| factor.B_subdiag()[index])
        .collect::<Vec<_>>();
    let scale = diagonal
        .iter()
        .chain(&subdiagonal)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let tolerance = 4096.0 * f64::EPSILON * dimension as f64 * scale;
    let mut eigenvalues = Vec::with_capacity(dimension);
    let mut index = 0;
    while index < dimension {
        if index + 1 < dimension && subdiagonal[index] != 0.0 {
            let radius = (diagonal[index] - diagonal[index + 1]).hypot(2.0 * subdiagonal[index]);
            eigenvalues.push(0.5 * (diagonal[index] + diagonal[index + 1] + radius));
            eigenvalues.push(0.5 * (diagonal[index] + diagonal[index + 1] - radius));
            index += 2;
        } else {
            eigenvalues.push(diagonal[index]);
            index += 1;
        }
    }
    InertiaEvidence {
        positive: eigenvalues
            .iter()
            .filter(|value| **value > tolerance)
            .count(),
        negative: eigenvalues
            .iter()
            .filter(|value| **value < -tolerance)
            .count(),
        zero: eigenvalues
            .iter()
            .filter(|value| value.abs() <= tolerance)
            .count(),
    }
}

fn normalized_backward_error(matrix: &DenseMatrix, solution: &[f64], rhs: &[f64]) -> f64 {
    let residual = matrix
        .multiply_vector(solution)
        .iter()
        .zip(rhs)
        .map(|(actual, expected)| (actual - expected).abs())
        .fold(0.0_f64, f64::max);
    let matrix_norm = (0..matrix.rows)
        .map(|row| {
            (0..matrix.columns)
                .map(|column| matrix.get(row, column).abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let solution_norm = solution
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let rhs_norm = rhs.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    residual / (matrix_norm * solution_norm + rhs_norm).max(f64::MIN_POSITIVE)
}

fn recover_canonical(
    problem: &CanonicalProblem,
    primal: &[f64],
    slacks: Vec<f64>,
) -> CanonicalObservables {
    let coefficient_count = problem.kernel.rows;
    let field_coefficients = primal[..coefficient_count].to_vec();
    let standard_polynomial = [
        primal[coefficient_count],
        primal[coefficient_count + 1],
        primal[coefficient_count + 2],
        primal[coefficient_count + 3],
    ];
    let polynomial_coefficients = problem
        .normalization
        .to_physical_polynomial(standard_polynomial);
    let latent = primal[coefficient_count + POLYNOMIAL_DIMENSION];
    let functional_values = problem
        .kernel
        .multiply_vector(&field_coefficients)
        .iter()
        .zip(problem.polynomial.multiply_vector(&standard_polynomial))
        .map(|(kernel_value, polynomial_value)| kernel_value + polynomial_value)
        .collect::<Vec<_>>();
    let mut residuals = vec![
        functional_values[0] - latent,
        functional_values[1] - latent,
        latent - problem.truth_latent,
    ];
    residuals.extend(
        functional_values[2..]
            .iter()
            .zip(&problem.truth_values[2..])
            .map(|(actual, expected)| actual - expected),
    );
    let side_condition_violation = (0..POLYNOMIAL_DIMENSION)
        .map(|column| {
            (0..coefficient_count)
                .map(|row| problem.polynomial.get(row, column) * field_coefficients[row])
                .sum::<f64>()
                .abs()
        })
        .fold(0.0_f64, f64::max);
    let hard_violation = residuals
        .iter()
        .map(|value| value.abs())
        .fold(side_condition_violation, f64::max);
    let field_energy = dot_product(
        &field_coefficients,
        &problem.kernel.multiply_vector(&field_coefficients),
    );
    CanonicalObservables {
        field_coefficients,
        polynomial_coefficients,
        semantic_latents: vec![latent],
        functional_values,
        residuals,
        slacks,
        side_condition_violation,
        hard_violation,
        field_energy,
        objective: 0.5 * field_energy,
    }
}

fn scaling_round_trip_error(
    physical_primal: &[f64],
    factors: &[f64],
    normalization: CoordinateNormalization,
    physical_polynomial: [f64; POLYNOMIAL_DIMENSION],
) -> f64 {
    let variable_error = physical_primal
        .iter()
        .zip(factors)
        .map(|(value, factor)| {
            let round_trip = factor * (value / factor);
            (round_trip - value).abs() / value.abs().max(1.0)
        })
        .fold(0.0_f64, f64::max);
    let polynomial_round_trip = normalization
        .to_physical_polynomial(normalization.to_standard_polynomial(physical_polynomial));
    polynomial_round_trip
        .iter()
        .zip(physical_polynomial)
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(variable_error, f64::max)
}

fn manufactured_truth_error(problem: &CanonicalProblem, recovered: &CanonicalObservables) -> f64 {
    let observable_error = recovered
        .field_coefficients
        .iter()
        .zip(&problem.truth_coefficients)
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .chain(
            recovered
                .polynomial_coefficients
                .iter()
                .zip(problem.truth_polynomial)
                .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0)),
        )
        .chain(std::iter::once(
            (recovered.semantic_latents[0] - problem.truth_latent).abs()
                / problem.truth_latent.abs().max(1.0),
        ))
        .chain(
            recovered
                .functional_values
                .iter()
                .zip(&problem.truth_values)
                .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0)),
        )
        .fold(0.0_f64, f64::max);
    let truth_energy = dot_product(
        &problem.truth_coefficients,
        &problem.kernel.multiply_vector(&problem.truth_coefficients),
    );
    observable_error
        .max((recovered.field_energy - truth_energy).abs() / truth_energy.abs().max(1.0))
        .max((recovered.objective - 0.5 * truth_energy).abs() / (0.5 * truth_energy).abs().max(1.0))
}

fn canonical_observable_difference(
    left: &CanonicalObservables,
    right: &CanonicalObservables,
) -> f64 {
    left.field_coefficients
        .iter()
        .zip(&right.field_coefficients)
        .chain(
            left.polynomial_coefficients
                .iter()
                .zip(&right.polynomial_coefficients),
        )
        .chain(left.semantic_latents.iter().zip(&right.semantic_latents))
        .chain(left.functional_values.iter().zip(&right.functional_values))
        .chain(left.residuals.iter().zip(&right.residuals))
        .map(|(left, right)| (left - right).abs() / left.abs().max(right.abs()).max(1.0))
        .chain([
            (left.field_energy - right.field_energy).abs()
                / left
                    .field_energy
                    .abs()
                    .max(right.field_energy.abs())
                    .max(1.0),
            (left.objective - right.objective).abs()
                / left.objective.abs().max(right.objective.abs()).max(1.0),
        ])
        .fold(0.0_f64, f64::max)
}

fn dot_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum()
}

fn multiply_three_by_three(matrix: &DenseMatrix, vector: [f64; 3]) -> [f64; 3] {
    let product = matrix.multiply_vector(&vector);
    [product[0], product[1], product[2]]
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}
