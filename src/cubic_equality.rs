use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::householder::{
    apply_block_householder_sequence_on_the_left_in_place_scratch,
    apply_block_householder_sequence_on_the_left_in_place_with_conj,
};
use faer::linalg::solvers::Qr;
use faer::prelude::*;
use faer::{Conj, Side};

use crate::capacity::{CapacityExceededEvidence, plan_equality_capacity};
use crate::cubic::{CubicKernel, GlobalAnisotropyMetric};
use crate::faer_backend;
use crate::functional::{
    CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, UsageProvenance,
};
use crate::kkt::{EqualityKktSystem, KktFailure, KktSolveEvidence, solve_equality_kkt};
use crate::math::dot3;
use crate::numerical::{EQUALITY_KKT_POLICY_V1, SpectralRankDecision};

const POLYNOMIAL_DIMENSION: usize = 4;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct DenseMatrix {
    rows: usize,
    columns: usize,
    values: Vec<f64>,
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

    pub(crate) fn shape(&self) -> (usize, usize) {
        (self.rows, self.columns)
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.columns + column]
    }

    fn set(&mut self, row: usize, column: usize, value: f64) {
        self.values[row * self.columns + column] = value;
    }

    fn values(&self) -> &[f64] {
        &self.values
    }

    fn multiply_vector(&self, vector: &[f64]) -> Vec<f64> {
        debug_assert_eq!(self.columns, vector.len());
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
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PolynomialCoordinates {
    center: [f64; 3],
    length: f64,
    degenerate_extent: bool,
}

impl PolynomialCoordinates {
    fn from_functionals(
        functionals: &[CanonicalFunctional],
    ) -> Result<Self, RepresentationFailure> {
        let mut minimum = [f64::INFINITY; 3];
        let mut maximum = [f64::NEG_INFINITY; 3];
        for functional in functionals {
            for term in functional.terms() {
                for axis in 0..3 {
                    minimum[axis] = minimum[axis].min(term.support()[axis]);
                    maximum[axis] = maximum[axis].max(term.support()[axis]);
                }
            }
        }
        let center = std::array::from_fn(|axis| 0.5 * minimum[axis] + 0.5 * maximum[axis]);
        if center.iter().any(|value| !value.is_finite()) {
            return Err(RepresentationFailure::InvalidSolveCoordinateTransform {
                reason: "non-finite bounding-box center",
                solver_invoked: false,
            });
        }
        let length = functionals
            .iter()
            .flat_map(CanonicalFunctional::terms)
            .map(|term| stable_norm(subtract(term.support(), center)))
            .fold(0.0_f64, f64::max);
        if !length.is_finite() {
            return Err(RepresentationFailure::InvalidSolveCoordinateTransform {
                reason: "non-finite characteristic length",
                solver_invoked: false,
            });
        }
        let degenerate_extent = length == 0.0;
        let coordinates = Self {
            center,
            length: if degenerate_extent { 1.0 } else { length },
            degenerate_extent,
        };
        let kernel_scale = coordinates.length.powi(3);
        if !kernel_scale.is_finite() || kernel_scale <= 0.0 {
            return Err(RepresentationFailure::InvalidSolveCoordinateTransform {
                reason: "non-invertible Cubic field recovery scale",
                solver_invoked: false,
            });
        }
        Ok(coordinates)
    }

    pub(crate) fn center(self) -> [f64; 3] {
        self.center
    }

    pub(crate) fn length(self) -> f64 {
        self.length
    }

    pub(crate) fn degenerate_extent(self) -> bool {
        self.degenerate_extent
    }

    pub(crate) fn to_standard_functional(
        self,
        physical: &CanonicalFunctional,
    ) -> CanonicalFunctional {
        let terms = physical
            .terms()
            .iter()
            .map(|term| {
                FunctionalTerm::new(
                    std::array::from_fn(|axis| {
                        (term.support()[axis] - self.center[axis]) / self.length
                    }),
                    term.value_coefficient(),
                    term.gradient_coefficient().map(|value| value / self.length),
                )
            })
            .collect();
        CanonicalFunctional::new(physical.dimension(), terms)
            .expect("an invertible similarity transform preserves a nonzero finite functional")
    }

    pub(crate) fn to_physical_field_coefficients(self, standard: &[f64]) -> Vec<f64> {
        let kernel_scale = self.length.powi(3);
        standard
            .iter()
            .map(|coefficient| coefficient / kernel_scale)
            .collect()
    }

    fn to_standard_field_coefficients(self, physical: &[f64]) -> Vec<f64> {
        let kernel_scale = self.length.powi(3);
        physical
            .iter()
            .map(|coefficient| coefficient * kernel_scale)
            .collect()
    }

    fn is_valid_recovery_map(self) -> bool {
        self.center.iter().all(|value| value.is_finite())
            && self.length.is_finite()
            && self.length > 0.0
    }

    fn pairing(self, functional: &CanonicalFunctional) -> [f64; POLYNOMIAL_DIMENSION] {
        functional
            .terms()
            .iter()
            .fold([0.0; 4], |mut pairing, term| {
                let value = term.value_coefficient();
                let support = term.support();
                let gradient = term.gradient_coefficient();
                pairing[0] += value;
                for axis in 0..3 {
                    pairing[axis + 1] += (value * (support[axis] - self.center[axis])
                        + gradient[axis])
                        / self.length;
                }
                pairing
            })
    }

    pub(crate) fn to_standard(self, physical: [f64; 4]) -> [f64; 4] {
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

    pub(crate) fn to_physical(self, standard: [f64; 4]) -> [f64; 4] {
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

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CpdEvidence {
    pub(crate) fitting_functional_count: usize,
    pub(crate) polynomial_dimension: usize,
    pub(crate) polynomial_rank: usize,
    pub(crate) singular_values: Vec<f64>,
    pub(crate) polynomial_rrqr_ratio: f64,
    pub(crate) polynomial_svd_ratio: f64,
    pub(crate) polynomial_rank_reject_ratio: f64,
    pub(crate) polynomial_rank_accept_ratio: f64,
    pub(crate) null_space_defect: f64,
    pub(crate) reduced_symmetry_defect: f64,
    pub(crate) symmetry_defect_limit: f64,
    pub(crate) reduced_smallest_eigenvalue: f64,
    pub(crate) affine_reproduction_error: f64,
    pub(crate) solve_coordinate_center: [f64; 3],
    pub(crate) solve_coordinate_length: f64,
    pub(crate) degenerate_extent: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PolynomialRankEvidence {
    pub(crate) rrqr_ratio: f64,
    pub(crate) svd_ratio: f64,
    pub(crate) reject_ratio: f64,
    pub(crate) accept_ratio: f64,
    pub(crate) backend_invoked: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CubicFieldEnergyNormalization {
    factor: f64,
}

impl CubicFieldEnergyNormalization {
    fn all_hard() -> Self {
        Self { factor: 1.0 }
    }

    pub(crate) fn factor(self) -> f64 {
        self.factor
    }

    pub(crate) fn for_rescaled_units(
        self,
        length_scale: f64,
        field_scale: f64,
    ) -> Result<Self, FieldEnergyNormalizationError> {
        if !length_scale.is_finite() || length_scale <= 0.0 {
            return Err(FieldEnergyNormalizationError::InvalidLengthScale);
        }
        if !field_scale.is_finite() || field_scale <= 0.0 {
            return Err(FieldEnergyNormalizationError::InvalidFieldScale);
        }
        let factor = self.factor * length_scale.powi(3) / field_scale.powi(2);
        if !factor.is_finite() || factor <= 0.0 {
            return Err(FieldEnergyNormalizationError::NonFiniteRescaledFactor);
        }
        Ok(Self { factor })
    }

    pub(crate) fn dimensionless_energy(self, native_cubic_energy: f64) -> f64 {
        self.factor * native_cubic_energy
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FieldEnergyNormalizationError {
    InvalidLengthScale,
    InvalidFieldScale,
    NonFiniteRescaledFactor,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RepresentationFailure {
    EmptyRepresenterSpan,
    Capacity(Box<CapacityExceededEvidence>),
    InvalidSolveCoordinateTransform {
        reason: &'static str,
        solver_invoked: bool,
    },
    PolynomialRankDeficient {
        rank: usize,
        solver_invoked: bool,
        hidden_regularization_applied: bool,
    },
    PolynomialRankGrayZone {
        evidence: PolynomialRankEvidence,
    },
    NumericalDecisionGrayZone {
        stage: &'static str,
        solver_invoked: bool,
        hidden_regularization_applied: bool,
    },
    NullSpaceWorkspaceAllocation,
    NullSpaceDefect {
        observed: f64,
        limit: f64,
    },
    ReducedSymmetryContract {
        observed: f64,
        limit: f64,
    },
    ReducedPairingNotPositive {
        smallest_eigenvalue: f64,
        solver_invoked: bool,
        hidden_regularization_applied: bool,
    },
    AffineReproductionBackend(Box<KktFailure>),
    AffineReproductionContract {
        observed: f64,
        limit: f64,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct CubicRepresentation {
    fitting_uses: Vec<FunctionalUse>,
    metric: GlobalAnisotropyMetric,
    coordinates: PolynomialCoordinates,
    kernel: DenseMatrix,
    polynomial: DenseMatrix,
    evidence: CpdEvidence,
    field_energy_normalization: CubicFieldEnergyNormalization,
}

impl CubicRepresentation {
    pub(crate) fn new(
        fitting_uses: Vec<FunctionalUse>,
        metric: GlobalAnisotropyMetric,
    ) -> Result<Self, RepresentationFailure> {
        if fitting_uses.is_empty() {
            return Err(RepresentationFailure::EmptyRepresenterSpan);
        }
        let augmented_dimension = fitting_uses
            .len()
            .checked_add(POLYNOMIAL_DIMENSION)
            .unwrap_or(usize::MAX);
        plan_equality_capacity(augmented_dimension, augmented_dimension)
            .map_err(|failure| RepresentationFailure::Capacity(Box::new(failure)))?;
        let functionals = fitting_uses
            .iter()
            .map(|usage| usage.functional().clone())
            .collect::<Vec<_>>();
        let coordinates = PolynomialCoordinates::from_functionals(&functionals)?;
        let standard_functionals = functionals
            .iter()
            .map(|functional| coordinates.to_standard_functional(functional))
            .collect::<Vec<_>>();
        let polynomial = DenseMatrix::from_fn(
            standard_functionals.len(),
            POLYNOMIAL_DIMENSION,
            |row, column| {
                standard_functionals[row].evaluate_affine(
                    if column == 0 { 1.0 } else { 0.0 },
                    std::array::from_fn(|axis| if axis + 1 == column { 1.0 } else { 0.0 }),
                )
            },
        );
        let (singular_values, polynomial_rank, polynomial_rank_evidence) =
            verify_polynomial_rank(&polynomial)?;
        let kernel = assemble_kernel_pairing(&standard_functionals, &metric);
        let null_space = HouseholderNullSpace::new(&polynomial, polynomial_rank);
        let (mut reduced, null_space_defect) =
            materialize_reduced_pairing(&kernel, &polynomial, &null_space)?;
        if null_space_defect > EQUALITY_KKT_POLICY_V1.null_space_defect_limit {
            return Err(RepresentationFailure::NullSpaceDefect {
                observed: null_space_defect,
                limit: EQUALITY_KKT_POLICY_V1.null_space_defect_limit,
            });
        }
        let reduced_symmetry_defect = normalized_symmetry_defect(&reduced);
        let symmetry_defect_limit = EQUALITY_KKT_POLICY_V1.reduced_symmetry_multiplier
            * f64::EPSILON
            * reduced.rows.max(reduced.columns) as f64;
        if reduced_symmetry_defect > symmetry_defect_limit {
            return Err(RepresentationFailure::ReducedSymmetryContract {
                observed: reduced_symmetry_defect,
                limit: symmetry_defect_limit,
            });
        }
        symmetrize(&mut reduced);
        let reduced_smallest_eigenvalue = verify_reduced_pairing(&reduced)?;
        let affine_reproduction_error = affine_reproduction_error(&kernel, &polynomial)?;
        if affine_reproduction_error > EQUALITY_KKT_POLICY_V1.affine_reproduction_limit {
            return Err(RepresentationFailure::AffineReproductionContract {
                observed: affine_reproduction_error,
                limit: EQUALITY_KKT_POLICY_V1.affine_reproduction_limit,
            });
        }

        Ok(Self {
            fitting_uses,
            metric,
            coordinates,
            kernel,
            polynomial,
            evidence: CpdEvidence {
                fitting_functional_count: functionals.len(),
                polynomial_dimension: POLYNOMIAL_DIMENSION,
                polynomial_rank,
                singular_values,
                polynomial_rrqr_ratio: polynomial_rank_evidence.rrqr_ratio,
                polynomial_svd_ratio: polynomial_rank_evidence.svd_ratio,
                polynomial_rank_reject_ratio: polynomial_rank_evidence.reject_ratio,
                polynomial_rank_accept_ratio: polynomial_rank_evidence.accept_ratio,
                null_space_defect,
                reduced_symmetry_defect,
                symmetry_defect_limit,
                reduced_smallest_eigenvalue,
                affine_reproduction_error,
                solve_coordinate_center: coordinates.center(),
                solve_coordinate_length: coordinates.length(),
                degenerate_extent: coordinates.degenerate_extent(),
            },
            field_energy_normalization: CubicFieldEnergyNormalization::all_hard(),
        })
    }

    pub(crate) fn evidence(&self) -> &CpdEvidence {
        &self.evidence
    }

    pub(crate) fn kernel_pairing(&self) -> &DenseMatrix {
        &self.kernel
    }

    pub(crate) fn polynomial_pairing(&self) -> &DenseMatrix {
        &self.polynomial
    }

    pub(crate) fn polynomial_coordinates(&self) -> PolynomialCoordinates {
        self.coordinates
    }

    pub(crate) fn solve_coordinate_transform(&self) -> PolynomialCoordinates {
        self.coordinates
    }

    pub(crate) fn field_energy_normalization(&self) -> CubicFieldEnergyNormalization {
        self.field_energy_normalization
    }
}

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

    fn expand(&self, reduced: &[f64]) -> Result<Vec<f64>, RepresentationFailure> {
        debug_assert_eq!(reduced.len(), self.reduced_dimension());
        let mut embedded = Mat::<f64>::zeros(self.ambient_dimension, 1);
        for (index, value) in reduced.iter().enumerate() {
            embedded[(self.polynomial_rank + index, 0)] = *value;
        }
        let requirement = apply_block_householder_sequence_on_the_left_in_place_scratch::<f64>(
            self.ambient_dimension,
            self.factor.Q_coeff().nrows(),
            1,
        );
        let mut memory = MemBuffer::try_new(requirement)
            .map_err(|_| RepresentationFailure::NullSpaceWorkspaceAllocation)?;
        apply_block_householder_sequence_on_the_left_in_place_with_conj(
            self.factor.Q_basis(),
            self.factor.Q_coeff(),
            Conj::No,
            embedded.as_mut(),
            faer_backend::parallelism(),
            MemStack::new(&mut memory),
        );
        Ok((0..self.ambient_dimension)
            .map(|row| embedded[(row, 0)])
            .collect())
    }
}

fn verify_polynomial_rank(
    polynomial: &DenseMatrix,
) -> Result<(Vec<f64>, usize, PolynomialRankEvidence), RepresentationFailure> {
    let exact_nonzero_columns = (0..polynomial.columns)
        .filter(|column| (0..polynomial.rows).any(|row| polynomial.get(row, *column) != 0.0))
        .count();
    if exact_nonzero_columns != POLYNOMIAL_DIMENSION {
        return Err(RepresentationFailure::PolynomialRankDeficient {
            rank: exact_nonzero_columns,
            solver_invoked: false,
            hidden_regularization_applied: false,
        });
    }
    let qr = polynomial.to_faer().col_piv_qr();
    let rrqr_largest = (0..POLYNOMIAL_DIMENSION)
        .map(|index| qr.R()[(index, index)].abs())
        .fold(0.0_f64, f64::max);
    let rrqr_smallest = (0..POLYNOMIAL_DIMENSION)
        .map(|index| qr.R()[(index, index)].abs())
        .fold(f64::INFINITY, f64::min);
    let rrqr_ratio = if rrqr_largest > 0.0 {
        rrqr_smallest / rrqr_largest
    } else {
        0.0
    };
    let singular_values = polynomial.to_faer().singular_values().map_err(|_| {
        RepresentationFailure::NumericalDecisionGrayZone {
            stage: "polynomial-rank",
            solver_invoked: false,
            hidden_regularization_applied: false,
        }
    })?;
    let largest = singular_values.first().copied().unwrap_or(0.0);
    let smallest = singular_values.last().copied().unwrap_or(0.0);
    let dimension = polynomial.rows.max(polynomial.columns);
    let (reject, _) = EQUALITY_KKT_POLICY_V1.spectral_thresholds(dimension, largest);
    let (reject_ratio, accept_ratio) = EQUALITY_KKT_POLICY_V1.spectral_ratio_thresholds(dimension);
    let svd_ratio = if largest > 0.0 {
        smallest / largest
    } else {
        0.0
    };
    let rank = singular_values
        .iter()
        .filter(|singular_value| **singular_value > reject)
        .count();
    match EQUALITY_KKT_POLICY_V1.classify_spectral_ratio(dimension, svd_ratio) {
        SpectralRankDecision::Reject => {
            return Err(RepresentationFailure::PolynomialRankDeficient {
                rank,
                solver_invoked: false,
                hidden_regularization_applied: false,
            });
        }
        SpectralRankDecision::GrayZone => {
            return Err(RepresentationFailure::PolynomialRankGrayZone {
                evidence: PolynomialRankEvidence {
                    rrqr_ratio,
                    svd_ratio,
                    reject_ratio,
                    accept_ratio,
                    backend_invoked: false,
                },
            });
        }
        SpectralRankDecision::Accept => {}
    }
    Ok((
        singular_values,
        rank,
        PolynomialRankEvidence {
            rrqr_ratio,
            svd_ratio,
            reject_ratio,
            accept_ratio,
            backend_invoked: false,
        },
    ))
}

fn assemble_kernel_pairing(
    functionals: &[CanonicalFunctional],
    metric: &GlobalAnisotropyMetric,
) -> DenseMatrix {
    let mut kernel = DenseMatrix::from_fn(functionals.len(), functionals.len(), |_, _| 0.0);
    for row in 0..functionals.len() {
        for column in row..functionals.len() {
            let value = CubicKernel::pairing(&functionals[row], &functionals[column], metric);
            kernel.set(row, column, value);
            kernel.set(column, row, value);
        }
    }
    kernel
}

fn materialize_reduced_pairing(
    kernel: &DenseMatrix,
    polynomial: &DenseMatrix,
    null_space: &HouseholderNullSpace,
) -> Result<(DenseMatrix, f64), RepresentationFailure> {
    let reduced_dimension = null_space.reduced_dimension();
    let mut columns = Vec::with_capacity(reduced_dimension);
    let mut null_space_defect = 0.0_f64;
    for column in 0..reduced_dimension {
        let mut unit = vec![0.0; reduced_dimension];
        unit[column] = 1.0;
        let expanded = null_space.expand(&unit)?;
        for polynomial_column in 0..POLYNOMIAL_DIMENSION {
            let side_value = (0..polynomial.rows)
                .map(|row| polynomial.get(row, polynomial_column) * expanded[row])
                .sum::<f64>();
            null_space_defect = null_space_defect.max(side_value.abs());
        }
        columns.push(expanded);
    }
    let kernel_columns = columns
        .iter()
        .map(|column| kernel.multiply_vector(column))
        .collect::<Vec<_>>();
    let reduced = DenseMatrix::from_fn(reduced_dimension, reduced_dimension, |row, column| {
        dot_product(&columns[row], &kernel_columns[column])
    });
    Ok((reduced, null_space_defect))
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

fn verify_reduced_pairing(matrix: &DenseMatrix) -> Result<f64, RepresentationFailure> {
    if matrix.rows == 0 {
        return Ok(f64::INFINITY);
    }
    let eigenvalues = matrix
        .to_faer()
        .self_adjoint_eigenvalues(Side::Lower)
        .map_err(|_| RepresentationFailure::NumericalDecisionGrayZone {
            stage: "reduced-positivity",
            solver_invoked: false,
            hidden_regularization_applied: false,
        })?;
    let smallest = eigenvalues[0];
    let largest = *eigenvalues
        .last()
        .expect("the nonempty reduced pairing has an eigenvalue");
    let (reject, accept) = EQUALITY_KKT_POLICY_V1.spectral_thresholds(matrix.rows, largest);
    if smallest <= reject {
        return Err(RepresentationFailure::ReducedPairingNotPositive {
            smallest_eigenvalue: smallest,
            solver_invoked: false,
            hidden_regularization_applied: false,
        });
    }
    if smallest < accept {
        return Err(RepresentationFailure::NumericalDecisionGrayZone {
            stage: "reduced-positivity",
            solver_invoked: false,
            hidden_regularization_applied: false,
        });
    }
    matrix.to_faer().llt(Side::Lower).map_err(|_| {
        RepresentationFailure::NumericalDecisionGrayZone {
            stage: "reduced-cholesky",
            solver_invoked: false,
            hidden_regularization_applied: false,
        }
    })?;
    Ok(smallest)
}

fn affine_reproduction_error(
    kernel: &DenseMatrix,
    polynomial: &DenseMatrix,
) -> Result<f64, RepresentationFailure> {
    let equality_jacobian =
        DenseMatrix::from_fn(POLYNOMIAL_DIMENSION, kernel.rows, |row, column| {
            polynomial.get(column, row)
        });
    let equality_rhs = [0.0; POLYNOMIAL_DIMENSION];
    let mut maximum_error = 0.0_f64;
    for basis in 0..POLYNOMIAL_DIMENSION {
        let stationarity_rhs = (0..kernel.rows)
            .map(|row| polynomial.get(row, basis))
            .collect::<Vec<_>>();
        let evidence = solve_equality_kkt(&EqualityKktSystem {
            primal_variables: kernel.rows,
            equality_constraints: POLYNOMIAL_DIMENSION,
            hessian: kernel.values(),
            equality_jacobian: equality_jacobian.values(),
            stationarity_rhs: &stationarity_rhs,
            equality_rhs: &equality_rhs,
        })
        .map_err(|failure| RepresentationFailure::AffineReproductionBackend(Box::new(failure)))?;
        maximum_error =
            evidence
                .candidate
                .iter()
                .map(|coefficient| coefficient.abs())
                .chain(evidence.equality_multipliers.iter().enumerate().map(
                    |(coefficient, actual)| {
                        let expected = if coefficient == basis { 1.0 } else { 0.0 };
                        (actual - expected).abs()
                    },
                ))
                .fold(maximum_error, f64::max);
    }
    Ok(maximum_error)
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct HardEquality {
    functional: FunctionalUse,
    target: f64,
}

impl HardEquality {
    pub(crate) fn new(functional: FunctionalUse, target: f64) -> Self {
        Self { functional, target }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct EqualityAssemblyEvidence {
    pub(crate) primal_variables: usize,
    pub(crate) field_coefficients: usize,
    pub(crate) polynomial_coefficients: usize,
    pub(crate) side_conditions: usize,
    pub(crate) hard_equalities: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FieldSample {
    pub(crate) value: f64,
    pub(crate) gradient: [f64; 3],
}

#[derive(Debug, Clone)]
pub(crate) struct RecoveredCubicField {
    representers: Vec<CanonicalFunctional>,
    metric: GlobalAnisotropyMetric,
    coefficients: Vec<f64>,
    physical_polynomial: [f64; 4],
}

impl RecoveredCubicField {
    pub(crate) fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    pub(crate) fn physical_polynomial(&self) -> [f64; 4] {
        self.physical_polynomial
    }

    pub(crate) fn sample(&self, point: [f64; 3]) -> FieldSample {
        let mut value = self.physical_polynomial[0]
            + self.physical_polynomial[1] * point[0]
            + self.physical_polynomial[2] * point[1]
            + self.physical_polynomial[3] * point[2];
        let mut gradient = [
            self.physical_polynomial[1],
            self.physical_polynomial[2],
            self.physical_polynomial[3],
        ];
        for (coefficient, functional) in self.coefficients.iter().zip(&self.representers) {
            for term in functional.terms() {
                let jet = CubicKernel::jet(point, term.support(), &self.metric);
                let derivative = term.gradient_coefficient();
                let representer_value =
                    term.value_coefficient() * jet.value() + dot3(derivative, jet.gradient_y());
                value += coefficient * representer_value;
                let mixed_derivative = jet.mixed_xy();
                for axis in 0..3 {
                    let representer_derivative = term.value_coefficient() * jet.gradient_x()[axis]
                        + dot3(mixed_derivative[axis], derivative);
                    gradient[axis] += coefficient * representer_derivative;
                }
            }
        }
        FieldSample { value, gradient }
    }

    fn evaluate_functional(&self, functional: &CanonicalFunctional) -> f64 {
        functional
            .terms()
            .iter()
            .map(|term| {
                let sample = self.sample(term.support());
                term.value_coefficient() * sample.value
                    + dot3(term.gradient_coefficient(), sample.gradient)
            })
            .sum()
    }

    fn native_cubic_energy(&self) -> f64 {
        self.coefficients
            .iter()
            .enumerate()
            .flat_map(|(row, left)| {
                self.coefficients
                    .iter()
                    .enumerate()
                    .map(move |(column, right)| {
                        left * right
                            * CubicKernel::pairing(
                                &self.representers[row],
                                &self.representers[column],
                                &self.metric,
                            )
                    })
            })
            .sum()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveredHardEquality {
    pub(crate) usage: FunctionalUse,
    pub(crate) target: f64,
    pub(crate) value: f64,
    pub(crate) residual: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct FunctionalViolationEnvelope {
    pub(crate) field_value: f64,
    pub(crate) field_value_per_length: f64,
}

impl FunctionalViolationEnvelope {
    fn from_dimensioned_residuals(
        residuals: impl IntoIterator<Item = (FunctionalDimension, f64)>,
    ) -> Self {
        residuals
            .into_iter()
            .fold(Self::default(), |mut envelope, (dimension, residual)| {
                match dimension {
                    FunctionalDimension::FieldValue => {
                        envelope.field_value = envelope.field_value.max(residual.abs());
                    }
                    FunctionalDimension::FieldValuePerLength => {
                        envelope.field_value_per_length =
                            envelope.field_value_per_length.max(residual.abs());
                    }
                }
                envelope
            })
    }

    fn is_within_policy(self) -> bool {
        self.field_value <= EQUALITY_KKT_POLICY_V1.field_value_recovery_limit
            && self.field_value_per_length <= EQUALITY_KKT_POLICY_V1.field_derivative_recovery_limit
    }
}

#[derive(Debug, Clone)]
pub(crate) struct CubicEqualitySolution {
    pub(crate) representation: CpdEvidence,
    pub(crate) assembly: EqualityAssemblyEvidence,
    pub(crate) backend: KktSolveEvidence,
    pub(crate) field: RecoveredCubicField,
    pub(crate) hard_equalities: Vec<RecoveredHardEquality>,
    pub(crate) side_condition_violation: f64,
    pub(crate) hard_equality_violations: FunctionalViolationEnvelope,
    pub(crate) polynomial_round_trip_error: f64,
    pub(crate) field_coefficient_round_trip_error: f64,
    pub(crate) field_energy_round_trip_error: f64,
    pub(crate) recovery_finite: bool,
    pub(crate) provenance_verified: bool,
    pub(crate) semantic_latent_count: usize,
    pub(crate) field_energy: f64,
    pub(crate) total_objective: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RecoveryVerificationFailureReason {
    InvalidRecoveryMap,
    ProvenanceMismatch,
    NonFiniteRecoveredQuantity,
    SideConditionViolation,
    HardEqualityViolation,
    PolynomialRoundTripViolation,
    FieldCoefficientRoundTripViolation,
    FieldEnergyRoundTripViolation,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveryVerificationFailureEvidence {
    pub(crate) reasons: Vec<RecoveryVerificationFailureReason>,
    pub(crate) side_condition_violation: Option<f64>,
    pub(crate) hard_equality_violations: Option<FunctionalViolationEnvelope>,
    pub(crate) polynomial_round_trip_error: Option<f64>,
    pub(crate) field_coefficient_round_trip_error: Option<f64>,
    pub(crate) field_energy_round_trip_error: Option<f64>,
    pub(crate) no_model_produced: bool,
}

impl RecoveryVerificationFailureEvidence {
    fn early(reason: RecoveryVerificationFailureReason) -> Self {
        Self {
            reasons: vec![reason],
            side_condition_violation: None,
            hard_equality_violations: None,
            polynomial_round_trip_error: None,
            field_coefficient_round_trip_error: None,
            field_energy_round_trip_error: None,
            no_model_produced: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct RecoveryProvenanceMap {
    usages: Vec<UsageProvenance>,
}

impl RecoveryProvenanceMap {
    fn from_representation(representation: &CubicRepresentation) -> Self {
        Self {
            usages: representation
                .fitting_uses
                .iter()
                .map(|usage| usage.provenance().clone())
                .collect(),
        }
    }

    fn verifies(&self, representation: &CubicRepresentation) -> bool {
        self.usages.len() == representation.fitting_uses.len()
            && self
                .usages
                .iter()
                .zip(&representation.fitting_uses)
                .all(|(expected, actual)| expected == actual.provenance())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CubicEqualityFailure {
    EmptyEqualitySet,
    NonFiniteTarget { equality: usize },
    Representation(Box<RepresentationFailure>),
    Backend(Box<KktFailure>),
    RecoveryVerification(Box<RecoveryVerificationFailureEvidence>),
}

pub(crate) struct CubicEqualityCore;

impl CubicEqualityCore {
    pub(crate) fn solve(
        equalities: Vec<HardEquality>,
        metric: GlobalAnisotropyMetric,
    ) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
        if equalities.is_empty() {
            return Err(CubicEqualityFailure::EmptyEqualitySet);
        }
        for (index, equality) in equalities.iter().enumerate() {
            if !equality.target.is_finite() {
                return Err(CubicEqualityFailure::NonFiniteTarget { equality: index });
            }
        }
        let targets = equalities
            .iter()
            .map(|equality| equality.target)
            .collect::<Vec<_>>();
        let representation = CubicRepresentation::new(
            equalities
                .into_iter()
                .map(|equality| equality.functional)
                .collect(),
            metric,
        )
        .map_err(|failure| CubicEqualityFailure::Representation(Box::new(failure)))?;
        let recovery_provenance = RecoveryProvenanceMap::from_representation(&representation);
        let (assembly, backend) = solve_standard_form(&representation, &targets)?;
        recover_and_verify(
            representation,
            targets,
            assembly,
            backend,
            recovery_provenance,
        )
    }
}

fn solve_standard_form(
    representation: &CubicRepresentation,
    targets: &[f64],
) -> Result<(EqualityAssemblyEvidence, KktSolveEvidence), CubicEqualityFailure> {
    let coefficient_count = representation.kernel.rows;
    let primal_variables = coefficient_count + POLYNOMIAL_DIMENSION;
    let equality_constraints = POLYNOMIAL_DIMENSION + targets.len();
    let hessian = DenseMatrix::from_fn(primal_variables, primal_variables, |row, column| {
        if row < coefficient_count && column < coefficient_count {
            representation.kernel.get(row, column)
        } else {
            0.0
        }
    });
    let equality_jacobian =
        DenseMatrix::from_fn(equality_constraints, primal_variables, |row, column| {
            if row < POLYNOMIAL_DIMENSION {
                if column < coefficient_count {
                    representation.polynomial.get(column, row)
                } else {
                    0.0
                }
            } else {
                let functional = row - POLYNOMIAL_DIMENSION;
                if column < coefficient_count {
                    representation.kernel.get(functional, column)
                } else {
                    representation
                        .polynomial
                        .get(functional, column - coefficient_count)
                }
            }
        });
    let stationarity_rhs = vec![0.0; primal_variables];
    let equality_rhs = std::iter::repeat_n(0.0, POLYNOMIAL_DIMENSION)
        .chain(targets.iter().copied())
        .collect::<Vec<_>>();
    let backend = solve_equality_kkt(&EqualityKktSystem {
        primal_variables,
        equality_constraints,
        hessian: hessian.values(),
        equality_jacobian: equality_jacobian.values(),
        stationarity_rhs: &stationarity_rhs,
        equality_rhs: &equality_rhs,
    })
    .map_err(|failure| CubicEqualityFailure::Backend(Box::new(failure)))?;
    Ok((
        EqualityAssemblyEvidence {
            primal_variables,
            field_coefficients: coefficient_count,
            polynomial_coefficients: POLYNOMIAL_DIMENSION,
            side_conditions: POLYNOMIAL_DIMENSION,
            hard_equalities: targets.len(),
        },
        backend,
    ))
}

fn recover_and_verify(
    representation: CubicRepresentation,
    targets: Vec<f64>,
    assembly: EqualityAssemblyEvidence,
    backend: KktSolveEvidence,
    recovery_provenance: RecoveryProvenanceMap,
) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
    if !representation.coordinates.is_valid_recovery_map() {
        return Err(CubicEqualityFailure::RecoveryVerification(Box::new(
            RecoveryVerificationFailureEvidence::early(
                RecoveryVerificationFailureReason::InvalidRecoveryMap,
            ),
        )));
    }
    let provenance_verified = recovery_provenance.verifies(&representation);
    if !provenance_verified {
        return Err(CubicEqualityFailure::RecoveryVerification(Box::new(
            RecoveryVerificationFailureEvidence::early(
                RecoveryVerificationFailureReason::ProvenanceMismatch,
            ),
        )));
    }

    let coefficient_count = assembly.field_coefficients;
    let standard_coefficients = &backend.candidate[..coefficient_count];
    let coefficients = representation
        .coordinates
        .to_physical_field_coefficients(standard_coefficients);
    let standard_polynomial: [f64; 4] = backend.candidate
        [coefficient_count..coefficient_count + POLYNOMIAL_DIMENSION]
        .try_into()
        .expect("the augmented primal retains exactly four polynomial coefficients");
    let physical_polynomial = representation.coordinates.to_physical(standard_polynomial);
    let field = RecoveredCubicField {
        representers: representation
            .fitting_uses
            .iter()
            .map(|usage| usage.functional().clone())
            .collect(),
        metric: representation.metric.clone(),
        coefficients,
        physical_polynomial,
    };
    let hard_equalities = representation
        .fitting_uses
        .iter()
        .cloned()
        .zip(&targets)
        .map(|(usage, target)| {
            let value = field.evaluate_functional(usage.functional());
            RecoveredHardEquality {
                usage,
                target: *target,
                value,
                residual: value - *target,
            }
        })
        .collect::<Vec<_>>();
    let side_condition_violation = (0..POLYNOMIAL_DIMENSION)
        .map(|column| {
            representation
                .fitting_uses
                .iter()
                .zip(&field.coefficients)
                .map(|(usage, coefficient)| {
                    representation.coordinates.pairing(usage.functional())[column] * coefficient
                })
                .sum::<f64>()
                .abs()
        })
        .fold(0.0_f64, f64::max);
    let hard_equality_violations = FunctionalViolationEnvelope::from_dimensioned_residuals(
        hard_equalities
            .iter()
            .map(|equality| (equality.usage.functional().dimension(), equality.residual)),
    );
    let recovered_standard_polynomial = representation
        .coordinates
        .to_standard(field.physical_polynomial);
    let polynomial_round_trip_error =
        relative_slice_error(&recovered_standard_polynomial, &standard_polynomial);
    let recovered_standard_coefficients = representation
        .coordinates
        .to_standard_field_coefficients(&field.coefficients);
    let field_coefficient_round_trip_error =
        relative_slice_error(&recovered_standard_coefficients, standard_coefficients);

    let native_cubic_energy = field.native_cubic_energy();
    let field_energy = representation
        .field_energy_normalization
        .dimensionless_energy(native_cubic_energy);
    let standard_energy = dot_product(
        standard_coefficients,
        &representation.kernel.multiply_vector(standard_coefficients),
    );
    let recovered_energy = standard_energy / representation.coordinates.length().powi(3);
    let field_energy_round_trip_error =
        (field_energy - recovered_energy).abs() / recovered_energy.abs().max(1.0);
    let total_objective = 0.5 * field_energy;
    let recovery_finite = field.coefficients.iter().all(|value| value.is_finite())
        && field
            .physical_polynomial
            .iter()
            .all(|value| value.is_finite())
        && hard_equalities
            .iter()
            .all(|equality| equality.value.is_finite() && equality.residual.is_finite())
        && side_condition_violation.is_finite()
        && polynomial_round_trip_error.is_finite()
        && field_coefficient_round_trip_error.is_finite()
        && field_energy.is_finite()
        && field_energy_round_trip_error.is_finite()
        && total_objective.is_finite();

    let mut reasons = Vec::new();
    if !recovery_finite {
        reasons.push(RecoveryVerificationFailureReason::NonFiniteRecoveredQuantity);
    }
    if side_condition_violation > EQUALITY_KKT_POLICY_V1.side_condition_limit {
        reasons.push(RecoveryVerificationFailureReason::SideConditionViolation);
    }
    if !hard_equality_violations.is_within_policy() {
        reasons.push(RecoveryVerificationFailureReason::HardEqualityViolation);
    }
    if polynomial_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::PolynomialRoundTripViolation);
    }
    if field_coefficient_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::FieldCoefficientRoundTripViolation);
    }
    if field_energy_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::FieldEnergyRoundTripViolation);
    }
    if !reasons.is_empty() {
        return Err(CubicEqualityFailure::RecoveryVerification(Box::new(
            RecoveryVerificationFailureEvidence {
                reasons,
                side_condition_violation: Some(side_condition_violation),
                hard_equality_violations: Some(hard_equality_violations),
                polynomial_round_trip_error: Some(polynomial_round_trip_error),
                field_coefficient_round_trip_error: Some(field_coefficient_round_trip_error),
                field_energy_round_trip_error: Some(field_energy_round_trip_error),
                no_model_produced: true,
            },
        )));
    }

    Ok(CubicEqualitySolution {
        representation: representation.evidence.clone(),
        assembly,
        backend,
        field,
        hard_equalities,
        side_condition_violation,
        hard_equality_violations,
        polynomial_round_trip_error,
        field_coefficient_round_trip_error,
        field_energy_round_trip_error,
        recovery_finite,
        provenance_verified,
        semantic_latent_count: 0,
        field_energy,
        total_objective,
    })
}

fn relative_slice_error(actual: &[f64], expected: &[f64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(0.0_f64, f64::max)
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn stable_norm(vector: [f64; 3]) -> f64 {
    let scale = vector
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        0.0
    } else if !scale.is_finite() {
        f64::INFINITY
    } else {
        scale
            * vector
                .map(|value| (value / scale).powi(2))
                .into_iter()
                .sum::<f64>()
                .sqrt()
    }
}

fn dot_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic::GlobalAnisotropyMetric;
    use crate::functional::{
        CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, GroupId,
        RelationId, SemanticRolePath, SourceId, UsageProvenance,
    };

    const TRUTH_COEFFICIENTS: [f64; 10] = [
        0.195, -0.105, -0.17, -0.10, 0.10, -0.07, 0.12, -0.05, 0.08, -0.04,
    ];
    // Independently evaluated with CPython Decimal at 120 digits from the
    // issue #14 manufactured declaration. These are inputs to the product
    // solve, never values generated by the product implementation.
    const MANUFACTURED_TARGETS: [f64; 10] = [
        0.204_930_731_175_453_18,
        0.204_930_731_175_453_18,
        -1.670_953_367_752_023_9,
        -2.060_930_563_063_515_7,
        3.156_001_299_339_634_2,
        1.033_400_763_358_782_7,
        0.613_739_534_748_635_2,
        0.764_288_970_221_123_3,
        0.840_004_985_206_390_4,
        0.891_736_553_208_148_7,
    ];
    const TRUTH_POLYNOMIAL: [f64; 4] = [0.6, 1.117_364_783_083_071_7, 0.4, 0.15];

    fn functional(
        support: [f64; 3],
        value_coefficient: f64,
        gradient_coefficient: [f64; 3],
    ) -> CanonicalFunctional {
        CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![FunctionalTerm::new(
                support,
                value_coefficient,
                gradient_coefficient,
            )],
        )
        .expect("manufactured functional is nonzero and finite")
    }

    fn manufactured_functionals() -> Vec<CanonicalFunctional> {
        vec![
            functional([-1.0, -1.0, -1.0], 1.0, [0.0; 3]),
            functional([1.0, -1.0, -1.0], 1.0, [0.0; 3]),
            functional([-1.0, 1.0, -1.0], 1.0, [0.0; 3]),
            functional([-1.0, -1.0, 1.0], 1.0, [0.0; 3]),
            functional([1.0, 1.0, 1.0], 1.0, [0.0; 3]),
            functional([0.25, -0.5, 0.75], 0.0, [1.0, 0.0, 0.0]),
            functional([-0.75, 0.25, 0.5], 0.0, [0.0, 1.0, 0.0]),
            functional([0.5, 0.75, -0.25], 0.0, [0.0, 0.0, 1.0]),
            functional([0.0, 0.0, 0.0], 1.0, [0.5, -0.25, 0.125]),
            functional([-0.5, 0.625, -0.75], 0.0, [1.0, 1.0, 1.0]),
        ]
    }

    fn usages() -> Vec<FunctionalUse> {
        manufactured_functionals()
            .into_iter()
            .enumerate()
            .map(|(index, functional)| {
                FunctionalUse::new(
                    functional,
                    UsageProvenance::new(
                        SourceId::new(format!("issue-17-manufactured-{index}")),
                        Some(GroupId::new("issue-17-manufactured")),
                        RelationId::new(format!("equality-{index}")),
                        SemanticRolePath::new(format!("hard-equality/{index}")),
                    ),
                )
            })
            .collect()
    }

    fn assert_close(actual: f64, expected: f64, limit: f64) {
        assert!(
            (actual - expected).abs() <= limit * expected.abs().max(1.0),
            "actual={actual:e}, expected={expected:e}, relative limit={limit:e}"
        );
    }

    #[test]
    fn cubic_representation_retains_full_pi1_and_passes_cpd_preflight() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representer span is Cubic-admissible");
        let evidence = representation.evidence();

        assert_eq!(evidence.fitting_functional_count, 10);
        assert_eq!(evidence.polynomial_dimension, 4);
        assert_eq!(evidence.polynomial_rank, 4);
        assert!(evidence.polynomial_rrqr_ratio > evidence.polynomial_rank_accept_ratio);
        assert!(evidence.polynomial_svd_ratio > evidence.polynomial_rank_accept_ratio);
        assert_eq!(representation.kernel_pairing().shape(), (10, 10));
        assert_eq!(representation.polynomial_pairing().shape(), (10, 4));
        assert!(evidence.null_space_defect <= 1.0e-12);
        assert!(evidence.reduced_smallest_eigenvalue > 0.0);
        assert!(evidence.reduced_symmetry_defect <= evidence.symmetry_defect_limit);
        assert!(evidence.affine_reproduction_error <= 1.0e-11);
        let normalization = representation.field_energy_normalization();
        assert_eq!(normalization.factor(), 1.0);
        let transformed = normalization
            .for_rescaled_units(2.0, 4.0)
            .expect("finite positive unit scales preserve normalization validity");
        assert_close(transformed.factor(), 0.5, 1.0e-15);
        assert_close(
            transformed.dimensionless_energy(6.0),
            normalization.dimensionless_energy(3.0),
            1.0e-15,
        );

        let standard = representation
            .polynomial_coordinates()
            .to_standard(TRUTH_POLYNOMIAL);
        let round_trip = representation
            .polynomial_coordinates()
            .to_physical(standard);
        for (actual, expected) in round_trip.into_iter().zip(TRUTH_POLYNOMIAL) {
            assert_close(actual, expected, 1.0e-15);
        }
    }

    #[test]
    fn solve_coordinates_are_deterministic_and_recover_physical_field_coefficients() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representer span is Cubic-admissible");
        let transform = representation.solve_coordinate_transform();

        assert_eq!(transform.center(), [0.0; 3]);
        assert_close(transform.length(), 3.0_f64.sqrt(), 1.0e-15);
        assert!(!transform.degenerate_extent());

        let physical = manufactured_functionals()[5].clone();
        let standard = transform.to_standard_functional(&physical);
        let physical_term = physical.terms()[0];
        let standard_term = standard.terms()[0];
        for (actual, expected) in standard_term
            .support()
            .into_iter()
            .zip(physical_term.support().map(|value| value / 3.0_f64.sqrt()))
        {
            assert_close(actual, expected, 1.0e-15);
        }
        for (actual, expected) in standard_term.gradient_coefficient().into_iter().zip(
            physical_term
                .gradient_coefficient()
                .map(|value| value / 3.0_f64.sqrt()),
        ) {
            assert_close(actual, expected, 1.0e-15);
        }

        let recovered = transform.to_physical_field_coefficients(&[3.0_f64.sqrt().powi(3)]);
        assert_close(recovered[0], 1.0, 1.0e-15);
    }

    #[test]
    fn cubic_equality_recovers_the_manufactured_field_through_the_full_augmented_kkt() {
        let fitting_uses = usages();
        let expected_uses = fitting_uses.clone();
        let equalities = fitting_uses
            .into_iter()
            .zip(MANUFACTURED_TARGETS)
            .map(|(functional, target)| HardEquality::new(functional, target))
            .collect();
        let solution = CubicEqualityCore::solve(equalities, GlobalAnisotropyMetric::identity())
            .expect("the valid manufactured field should recover");

        assert_eq!(solution.assembly.primal_variables, 14);
        assert_eq!(solution.assembly.field_coefficients, 10);
        assert_eq!(solution.assembly.polynomial_coefficients, 4);
        assert_eq!(solution.assembly.side_conditions, 4);
        assert_eq!(solution.assembly.hard_equalities, 10);
        assert_eq!(solution.backend.capacity.kkt_dimension, 28);
        assert!(solution.backend.normalized_backward_error <= 1.0e-11);
        for ((recovered, expected_usage), expected_target) in solution
            .hard_equalities
            .iter()
            .zip(&expected_uses)
            .zip(MANUFACTURED_TARGETS)
        {
            assert_eq!(&recovered.usage, expected_usage);
            assert_eq!(
                recovered.usage.functional().dimension(),
                FunctionalDimension::FieldValue
            );
            assert_eq!(recovered.target, expected_target);
            assert!(recovered.residual.abs() <= 1.0e-8);
        }
        for (actual, expected) in solution
            .field
            .coefficients()
            .iter()
            .copied()
            .zip(TRUTH_COEFFICIENTS)
        {
            assert_close(actual, expected, 1.0e-8);
        }
        for (actual, expected) in solution
            .field
            .physical_polynomial()
            .into_iter()
            .zip(TRUTH_POLYNOMIAL)
        {
            assert_close(actual, expected, 1.0e-8);
        }
        assert!(solution.side_condition_violation <= 1.0e-10);
        assert!(solution.hard_equality_violations.field_value <= 1.0e-8);
        assert_eq!(
            solution.hard_equality_violations.field_value_per_length,
            0.0
        );
        assert!(solution.polynomial_round_trip_error <= 1.0e-11);
        assert!(solution.field_coefficient_round_trip_error <= 1.0e-11);
        assert!(solution.field_energy_round_trip_error <= 1.0e-11);
        assert!(solution.recovery_finite);
        assert!(solution.provenance_verified);
        assert_eq!(solution.semantic_latent_count, 0);
        assert_close(solution.field_energy, 0.818_826_203_475_800_5, 1.0e-8);
        assert_close(solution.total_objective, 0.409_413_101_737_900_25, 1.0e-8);

        let sample = solution.field.sample([0.2, -0.3, 0.4]);
        assert_close(sample.value, 0.334_298_988_103_683_4, 1.0e-8);
        for (actual, expected) in sample.gradient.into_iter().zip([
            1.044_736_659_872_982_4,
            0.933_851_542_772_418_8,
            -0.142_406_691_059_823_7,
        ]) {
            assert_close(actual, expected, 1.0e-8);
        }
    }

    #[test]
    fn rank_deficient_pi1_is_rejected_without_hidden_repair() {
        let flattened = usages()
            .into_iter()
            .map(|usage| {
                let term = usage.functional().terms()[0];
                let mut support = term.support();
                let mut gradient = term.gradient_coefficient();
                support[2] = 0.0;
                gradient[2] = 0.0;
                if term.value_coefficient() == 0.0 && gradient == [0.0; 3] {
                    gradient[0] = 1.0;
                }
                FunctionalUse::new(
                    functional(support, term.value_coefficient(), gradient),
                    usage.provenance().clone(),
                )
            })
            .collect();

        let failure = match CubicRepresentation::new(flattened, GlobalAnisotropyMetric::identity())
        {
            Ok(_) => panic!("a rank-three Pi1 pairing must not reach a solver"),
            Err(failure) => failure,
        };
        assert_eq!(
            failure,
            RepresentationFailure::PolynomialRankDeficient {
                rank: 3,
                solver_invoked: false,
                hidden_regularization_applied: false,
            }
        );
    }

    #[test]
    fn polynomial_rank_band_returns_stable_gray_zone_evidence() {
        let gray_ratio = 1.0e-13;
        let polynomial = DenseMatrix::from_fn(4, 4, |row, column| {
            if row == column {
                if row == 3 { gray_ratio } else { 1.0 }
            } else {
                0.0
            }
        });

        let failure = verify_polynomial_rank(&polynomial)
            .expect_err("a ratio between the reject and accept bands is not guessed");
        match failure {
            RepresentationFailure::PolynomialRankGrayZone { evidence } => {
                assert!(evidence.rrqr_ratio > evidence.reject_ratio);
                assert!(evidence.rrqr_ratio < evidence.accept_ratio);
                assert!(evidence.svd_ratio > evidence.reject_ratio);
                assert!(evidence.svd_ratio < evidence.accept_ratio);
                assert!(!evidence.backend_invoked);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn unrepresentable_similarity_scale_fails_before_assembly_or_backend_invocation() {
        let extreme = [
            [1.0e308, 1.0e308, 1.0e308],
            [-1.0e308, 1.0e308, -1.0e308],
            [1.0e308, -1.0e308, -1.0e308],
            [-1.0e308, -1.0e308, 1.0e308],
        ]
        .into_iter()
        .enumerate()
        .map(|(index, support)| {
            FunctionalUse::new(
                functional(support, 1.0, [0.0; 3]),
                UsageProvenance::new(
                    SourceId::new(format!("extreme-{index}")),
                    None,
                    RelationId::new(format!("extreme-relation-{index}")),
                    SemanticRolePath::new(format!("hard-equality/{index}")),
                ),
            )
        })
        .collect();

        let failure = CubicRepresentation::new(extreme, GlobalAnisotropyMetric::identity())
            .expect_err("an infinite derived radius has no reversible f64 recovery map");
        assert_eq!(
            failure,
            RepresentationFailure::InvalidSolveCoordinateTransform {
                reason: "non-invertible Cubic field recovery scale",
                solver_invoked: false,
            }
        );
    }

    #[test]
    fn recovery_violation_envelopes_keep_value_and_derivative_units_separate() {
        let envelope = FunctionalViolationEnvelope::from_dimensioned_residuals([
            (FunctionalDimension::FieldValue, 2.0e-9),
            (FunctionalDimension::FieldValuePerLength, -7.0e-9),
        ]);

        assert_eq!(envelope.field_value, 2.0e-9);
        assert_eq!(envelope.field_value_per_length, 7.0e-9);
        assert!(envelope.is_within_policy());
    }

    #[test]
    fn damaged_coordinate_map_is_a_recovery_failure_not_a_backend_contract_failure() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representation is valid");
        let targets = MANUFACTURED_TARGETS.to_vec();
        let recovery_provenance = RecoveryProvenanceMap::from_representation(&representation);
        let (assembly, backend) = solve_standard_form(&representation, &targets)
            .expect("the undamaged standard form should produce a backend candidate");
        let mut damaged = representation;
        damaged.coordinates.length = 0.0;

        let failure = recover_and_verify(damaged, targets, assembly, backend, recovery_provenance)
            .expect_err("an invalid inverse coordinate map must fail during recovery");

        match failure {
            CubicEqualityFailure::RecoveryVerification(evidence) => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::InvalidRecoveryMap]
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_provenance_map_fails_recovery_without_returning_a_partial_model() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representation is valid");
        let targets = MANUFACTURED_TARGETS.to_vec();
        let mut recovery_provenance = RecoveryProvenanceMap::from_representation(&representation);
        recovery_provenance.usages[0] = UsageProvenance::new(
            SourceId::new("corrupted-source"),
            None,
            RelationId::new("corrupted-relation"),
            SemanticRolePath::new("corrupted-role"),
        );
        let (assembly, backend) = solve_standard_form(&representation, &targets)
            .expect("the standard form should still satisfy its backend contract");

        let failure = recover_and_verify(
            representation,
            targets,
            assembly,
            backend,
            recovery_provenance,
        )
        .expect_err("canonical provenance must round-trip before a model exists");

        match failure {
            CubicEqualityFailure::RecoveryVerification(evidence) => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::ProvenanceMismatch]
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }
}
