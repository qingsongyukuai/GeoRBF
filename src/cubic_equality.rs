use faer::Conj;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::householder::{
    apply_block_householder_sequence_on_the_left_in_place_scratch,
    apply_block_householder_sequence_on_the_left_in_place_with_conj,
};
use faer::prelude::*;
#[cfg(test)]
use std::cell::RefCell;

use crate::capacity::{CapacityExceededEvidence, plan_equality_capacity};
use crate::cubic::{CubicKernel, GlobalAnisotropyMetric};
use crate::faer_backend;
use crate::functional::{
    CanonicalFunctional, DerivedBlockId, DerivedColumnId, DerivedRowId, FunctionalDimension,
    FunctionalTerm, FunctionalUse, GroupId, ResidualId, SourceId, UsageProvenance,
};
use crate::geometry::FieldUnitLabel;
use crate::kernel::FieldEnergyNormalization;
use crate::kkt::{EqualityKktSystem, KktFailure, KktSolveEvidence, solve_equality_kkt};
use crate::math::dot3;
use crate::numerical::{
    EQUALITY_KKT_POLICY_V1, SpectralAnalysisFailure, SpectralRankDecision, analyze_spectral_rank,
};

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
pub(crate) struct CubicSolveCoordinateTransform {
    center: [f64; 3],
    length: f64,
    degenerate_extent: bool,
}

impl CubicSolveCoordinateTransform {
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
                reason: SolveCoordinateTransformFailureReason::BoundingBoxCenterNotFinite,
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
                reason: SolveCoordinateTransformFailureReason::CharacteristicLengthNotFinite,
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
                reason: SolveCoordinateTransformFailureReason::FieldRecoveryScaleNotInvertible,
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
    ) -> Result<CanonicalFunctional, RepresentationFailure> {
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
        CanonicalFunctional::new(physical.dimension(), terms).map_err(|_| {
            RepresentationFailure::InvalidSolveCoordinateTransform {
                reason: SolveCoordinateTransformFailureReason::StandardFunctionalNotFinite,
                solver_invoked: false,
            }
        })
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

    fn to_standard_tolerance(
        self,
        _dimension: FunctionalDimension,
        physical_tolerance: f64,
    ) -> f64 {
        // Functionals transform contragrediently, so their scalar output and
        // tolerance retain their physical units in standard coordinates.
        physical_tolerance
    }

    fn to_physical_tolerance(
        self,
        _dimension: FunctionalDimension,
        standard_tolerance: f64,
    ) -> f64 {
        standard_tolerance
    }

    fn to_physical_side_condition(self, standard: [f64; 4]) -> [f64; 4] {
        let field_scale = self.length.powi(3);
        [
            standard[0] / field_scale,
            (self.center[0] * standard[0] + self.length * standard[1]) / field_scale,
            (self.center[1] * standard[0] + self.length * standard[2]) / field_scale,
            (self.center[2] * standard[0] + self.length * standard[3]) / field_scale,
        ]
    }

    fn to_standard_side_condition(self, physical: [f64; 4]) -> [f64; 4] {
        let field_scale = self.length.powi(3);
        let constant = field_scale * physical[0];
        [
            constant,
            (field_scale * physical[1] - self.center[0] * constant) / self.length,
            (field_scale * physical[2] - self.center[1] * constant) / self.length,
            (field_scale * physical[3] - self.center[2] * constant) / self.length,
        ]
    }

    fn to_physical_side_condition_tolerances(self, standard: [f64; 4]) -> [f64; 4] {
        let field_scale = self.length.powi(3);
        [
            standard[0] / field_scale,
            (self.center[0].abs() * standard[0] + self.length * standard[1]) / field_scale,
            (self.center[1].abs() * standard[0] + self.length * standard[2]) / field_scale,
            (self.center[2].abs() * standard[0] + self.length * standard[3]) / field_scale,
        ]
    }

    fn is_valid_recovery_map(self) -> bool {
        self.center.iter().all(|value| value.is_finite())
            && self.length.is_finite()
            && self.length > 0.0
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
    pub(crate) reduced_smallest_singular_value: f64,
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
pub(crate) struct AnalysisExecutionEvidence {
    pub(crate) solver_invoked: bool,
    pub(crate) hidden_regularization_applied: bool,
}

impl AnalysisExecutionEvidence {
    fn pre_backend() -> Self {
        Self {
            solver_invoked: false,
            hidden_regularization_applied: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct VerifiedCanonicalMode {
    pub(crate) residual: f64,
    pub(crate) execution: AnalysisExecutionEvidence,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RepresentationFailure {
    EmptyRepresenterSpan,
    Capacity(Box<CapacityExceededEvidence>),
    InvalidSolveCoordinateTransform {
        reason: SolveCoordinateTransformFailureReason,
        solver_invoked: bool,
    },
    PolynomialRankDeficient {
        rank: Option<usize>,
        mode: VerifiedCanonicalMode,
    },
    PolynomialRankGrayZone {
        evidence: PolynomialRankEvidence,
    },
    ReducedPairingGrayZone(AnalysisExecutionEvidence),
    AlgebraicAnalysisFailure {
        stage: AlgebraicAnalysisStage,
        solver_invoked: bool,
    },
    AlgebraicAnalysisWorkspaceAllocation {
        stage: AlgebraicAnalysisStage,
        bytes: u64,
        alignment: usize,
        solver_invoked: bool,
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
        classification: ReducedPairingFailureClassification,
        rank: usize,
        negative_pivots: usize,
        execution: AnalysisExecutionEvidence,
    },
    AffineReproductionBackend(Box<KktFailure>),
    AffineReproductionContract {
        observed: f64,
        limit: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlgebraicAnalysisStage {
    PolynomialRank,
    ReducedCholesky,
    ReducedInertia,
    ReducedSpectrum,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReducedPairingFailureClassification {
    RankDeficient,
    NegativeCurvature,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SolveCoordinateTransformFailureReason {
    BoundingBoxCenterNotFinite,
    CharacteristicLengthNotFinite,
    FieldRecoveryScaleNotInvertible,
    StandardFunctionalNotFinite,
}

#[derive(Debug, Clone)]
pub(crate) struct CubicRepresentation {
    fitting_uses: Vec<FunctionalUse>,
    metric: GlobalAnisotropyMetric,
    coordinates: CubicSolveCoordinateTransform,
    kernel: DenseMatrix,
    polynomial: DenseMatrix,
    evidence: CpdEvidence,
    field_energy_normalization: FieldEnergyNormalization,
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
        let (coordinates, standard_functionals, polynomial) =
            assemble_polynomial_pairing(&functionals)?;
        let (singular_values, polynomial_rank, polynomial_rank_evidence) =
            verify_polynomial_rank(&polynomial, &functionals, &coordinates)?;
        let kernel = assemble_kernel_pairing(&standard_functionals, &metric);
        let null_space = HouseholderNullSpace::new(&polynomial, polynomial_rank)?;
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
        let reduced_smallest_singular_value = verify_reduced_pairing(&reduced)?;
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
                reduced_smallest_singular_value,
                affine_reproduction_error,
                solve_coordinate_center: coordinates.center(),
                solve_coordinate_length: coordinates.length(),
                degenerate_extent: coordinates.degenerate_extent(),
            },
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
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

    pub(crate) fn solve_coordinate_transform(&self) -> CubicSolveCoordinateTransform {
        self.coordinates
    }

    pub(crate) fn field_energy_normalization(&self) -> FieldEnergyNormalization {
        self.field_energy_normalization
    }

    fn set_field_energy_normalization(&mut self, normalization: FieldEnergyNormalization) {
        self.field_energy_normalization = normalization;
    }
}

pub(crate) fn preflight_polynomial_analysis_failure(
    fitting_uses: &[FunctionalUse],
) -> Option<RepresentationFailure> {
    if fitting_uses.is_empty() {
        return None;
    }
    let functionals = fitting_uses
        .iter()
        .map(|usage| usage.functional().clone())
        .collect::<Vec<_>>();
    let (coordinates, _, polynomial) = match assemble_polynomial_pairing(&functionals) {
        Ok(pairing) => pairing,
        Err(failure) => return Some(failure),
    };
    match verify_polynomial_rank(&polynomial, &functionals, &coordinates) {
        Ok(_) => None,
        Err(failure) => Some(failure),
    }
}

pub(crate) fn canonical_fitting_uses(
    equalities: &[CanonicalHardEquality],
    soft_equalities: &[CanonicalSoftEquality],
) -> Vec<FunctionalUse> {
    let mut fitting_uses = Vec::<FunctionalUse>::new();
    for usage in equalities
        .iter()
        .filter(|equality| {
            equality.participation() == CanonicalEqualityParticipation::SolverConstraint
        })
        .filter_map(CanonicalHardEquality::field)
        .chain(soft_equalities.iter().map(CanonicalSoftEquality::field))
    {
        if !fitting_uses
            .iter()
            .any(|existing| existing.functional() == usage.functional())
        {
            fitting_uses.push(usage.clone());
        }
    }
    fitting_uses
}

fn assemble_polynomial_pairing(
    functionals: &[CanonicalFunctional],
) -> Result<
    (
        CubicSolveCoordinateTransform,
        Vec<CanonicalFunctional>,
        DenseMatrix,
    ),
    RepresentationFailure,
> {
    let coordinates = CubicSolveCoordinateTransform::from_functionals(functionals)?;
    let standard_functionals = functionals
        .iter()
        .map(|functional| coordinates.to_standard_functional(functional))
        .collect::<Result<Vec<_>, _>>()?;
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
    Ok((coordinates, standard_functionals, polynomial))
}

struct HouseholderNullSpace {
    basis: Mat<f64>,
    coefficients: Mat<f64>,
    ambient_dimension: usize,
    polynomial_rank: usize,
}

impl HouseholderNullSpace {
    fn new(
        polynomial: &DenseMatrix,
        polynomial_rank: usize,
    ) -> Result<Self, RepresentationFailure> {
        let matrix = polynomial.to_faer();
        let factors = faer_backend::householder_qr(matrix.as_ref())
            .map_err(|_| RepresentationFailure::NullSpaceWorkspaceAllocation)?;
        Ok(Self {
            basis: factors.basis,
            coefficients: factors.coefficients,
            ambient_dimension: polynomial.rows,
            polynomial_rank,
        })
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
            self.coefficients.nrows(),
            1,
        );
        let mut memory = MemBuffer::try_new(requirement)
            .map_err(|_| RepresentationFailure::NullSpaceWorkspaceAllocation)?;
        apply_block_householder_sequence_on_the_left_in_place_with_conj(
            self.basis.as_ref(),
            self.coefficients.as_ref(),
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
    functionals: &[CanonicalFunctional],
    coordinates: &CubicSolveCoordinateTransform,
) -> Result<(Vec<f64>, usize, PolynomialRankEvidence), RepresentationFailure> {
    if let Some(column) = (0..polynomial.columns)
        .find(|column| (0..polynomial.rows).all(|row| polynomial.get(row, *column) == 0.0))
    {
        let nonzero_columns = (0..polynomial.columns)
            .filter(|candidate| {
                (0..polynomial.rows).any(|row| polynomial.get(row, *candidate) != 0.0)
            })
            .count();
        let structurally_proven_rank = (polynomial.rows == 1 || nonzero_columns == 1).then_some(1);
        let mut standard_mode = [0.0; POLYNOMIAL_DIMENSION];
        standard_mode[column] = 1.0;
        return polynomial_rank_failure(
            structurally_proven_rank,
            functionals,
            coordinates,
            standard_mode,
        );
    }
    let matrix = polynomial.to_faer();
    let analysis = analyze_spectral_rank(matrix.as_ref()).map_err(|failure| match failure {
        SpectralAnalysisFailure::WorkspaceAllocation(failure) => {
            RepresentationFailure::AlgebraicAnalysisWorkspaceAllocation {
                stage: AlgebraicAnalysisStage::PolynomialRank,
                bytes: failure.bytes,
                alignment: failure.alignment,
                solver_invoked: false,
            }
        }
        SpectralAnalysisFailure::NumericalError => {
            RepresentationFailure::AlgebraicAnalysisFailure {
                stage: AlgebraicAnalysisStage::PolynomialRank,
                solver_invoked: false,
            }
        }
    })?;
    if polynomial.rows < POLYNOMIAL_DIMENSION {
        let standard_mode =
            smallest_right_mode(matrix.as_ref(), AlgebraicAnalysisStage::PolynomialRank)?;
        return polynomial_rank_failure(
            Some(analysis.rank),
            functionals,
            coordinates,
            standard_mode
                .try_into()
                .expect("the Cubic polynomial pairing has exactly four columns"),
        );
    }
    match analysis.decision {
        SpectralRankDecision::Reject => {
            let standard_mode =
                smallest_right_mode(matrix.as_ref(), AlgebraicAnalysisStage::PolynomialRank)?;
            return polynomial_rank_failure(
                Some(analysis.rank),
                functionals,
                coordinates,
                standard_mode
                    .try_into()
                    .expect("the Cubic polynomial pairing has exactly four columns"),
            );
        }
        SpectralRankDecision::GrayZone => {
            return Err(RepresentationFailure::PolynomialRankGrayZone {
                evidence: PolynomialRankEvidence {
                    rrqr_ratio: analysis.rrqr_ratio,
                    svd_ratio: analysis.svd_ratio,
                    reject_ratio: analysis.reject_ratio,
                    accept_ratio: analysis.accept_ratio,
                    backend_invoked: false,
                },
            });
        }
        SpectralRankDecision::Accept => {}
    }
    Ok((
        analysis.singular_values,
        analysis.rank,
        PolynomialRankEvidence {
            rrqr_ratio: analysis.rrqr_ratio,
            svd_ratio: analysis.svd_ratio,
            reject_ratio: analysis.reject_ratio,
            accept_ratio: analysis.accept_ratio,
            backend_invoked: false,
        },
    ))
}

fn polynomial_rank_failure(
    rank: Option<usize>,
    functionals: &[CanonicalFunctional],
    coordinates: &CubicSolveCoordinateTransform,
    standard_mode: [f64; POLYNOMIAL_DIMENSION],
) -> Result<(Vec<f64>, usize, PolynomialRankEvidence), RepresentationFailure> {
    Err(RepresentationFailure::PolynomialRankDeficient {
        rank,
        mode: VerifiedCanonicalMode {
            residual: canonical_polynomial_mode_residual(functionals, coordinates, standard_mode)?,
            execution: AnalysisExecutionEvidence::pre_backend(),
        },
    })
}

fn smallest_right_mode(
    matrix: faer::MatRef<'_, f64>,
    stage: AlgebraicAnalysisStage,
) -> Result<Vec<f64>, RepresentationFailure> {
    faer_backend::smallest_right_singular_vector(matrix).map_err(|failure| match failure {
        faer_backend::DecompositionFailure::WorkspaceAllocation(failure) => {
            RepresentationFailure::AlgebraicAnalysisWorkspaceAllocation {
                stage,
                bytes: failure.bytes,
                alignment: failure.alignment,
                solver_invoked: false,
            }
        }
        faer_backend::DecompositionFailure::NumericalError => {
            RepresentationFailure::AlgebraicAnalysisFailure {
                stage,
                solver_invoked: false,
            }
        }
    })
}

fn canonical_polynomial_mode_residual(
    functionals: &[CanonicalFunctional],
    coordinates: &CubicSolveCoordinateTransform,
    mut standard_mode: [f64; POLYNOMIAL_DIMENSION],
) -> Result<f64, RepresentationFailure> {
    let scale = standard_mode
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return Err(RepresentationFailure::AlgebraicAnalysisFailure {
            stage: AlgebraicAnalysisStage::PolynomialRank,
            solver_invoked: false,
        });
    }
    for coefficient in &mut standard_mode {
        *coefficient /= scale;
    }
    let physical_mode = coordinates.to_physical(standard_mode);
    let recovered_standard_mode = coordinates.to_standard(physical_mode);
    let round_trip_error = relative_slice_error(&standard_mode, &recovered_standard_mode);
    let functional_residual = functionals
        .iter()
        .map(|functional| {
            let standard = coordinates
                .to_standard_functional(functional)
                .expect("an analyzed canonical functional has a finite standard form");
            let response = functional
                .evaluate_affine(
                    physical_mode[0],
                    [physical_mode[1], physical_mode[2], physical_mode[3]],
                )
                .abs();
            let response_scale = standard
                .terms()
                .iter()
                .map(|term| {
                    let support = term.support();
                    (term.value_coefficient()
                        * (standard_mode[0]
                            + standard_mode[1] * support[0]
                            + standard_mode[2] * support[1]
                            + standard_mode[3] * support[2]))
                        .abs()
                        + term
                            .gradient_coefficient()
                            .into_iter()
                            .zip(standard_mode[1..].iter().copied())
                            .map(|(left, right)| (left * right).abs())
                            .sum::<f64>()
                })
                .sum::<f64>()
                .max(1.0);
            response / response_scale
        })
        .fold(0.0_f64, f64::max);
    let residual = functional_residual.max(round_trip_error);
    if !residual.is_finite() {
        return Err(RepresentationFailure::AlgebraicAnalysisFailure {
            stage: AlgebraicAnalysisStage::PolynomialRank,
            solver_invoked: false,
        });
    }
    if residual > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        return Err(RepresentationFailure::AlgebraicAnalysisFailure {
            stage: AlgebraicAnalysisStage::PolynomialRank,
            solver_invoked: false,
        });
    }
    Ok(residual)
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
    let faer_matrix = matrix.to_faer();
    let cholesky = faer_backend::cholesky_minimum_diagonal(faer_matrix.as_ref());
    let failed_cholesky_inertia = match cholesky {
        Ok(_) => None,
        Err(faer_backend::CholeskyFailure::WorkspaceAllocation(failure)) => {
            return Err(
                RepresentationFailure::AlgebraicAnalysisWorkspaceAllocation {
                    stage: AlgebraicAnalysisStage::ReducedCholesky,
                    bytes: failure.bytes,
                    alignment: failure.alignment,
                    solver_invoked: false,
                },
            );
        }
        Err(faer_backend::CholeskyFailure::NonPositivePivot) => Some(
            faer_backend::bunch_kaufman_inertia(faer_matrix.as_ref()).map_err(|failure| {
                match failure {
                    faer_backend::DecompositionFailure::WorkspaceAllocation(failure) => {
                        RepresentationFailure::AlgebraicAnalysisWorkspaceAllocation {
                            stage: AlgebraicAnalysisStage::ReducedInertia,
                            bytes: failure.bytes,
                            alignment: failure.alignment,
                            solver_invoked: false,
                        }
                    }
                    faer_backend::DecompositionFailure::NumericalError => {
                        RepresentationFailure::AlgebraicAnalysisFailure {
                            stage: AlgebraicAnalysisStage::ReducedInertia,
                            solver_invoked: false,
                        }
                    }
                }
            })?,
        ),
    };
    let spectrum =
        analyze_spectral_rank(faer_matrix.as_ref()).map_err(|failure| match failure {
            SpectralAnalysisFailure::WorkspaceAllocation(failure) => {
                RepresentationFailure::AlgebraicAnalysisWorkspaceAllocation {
                    stage: AlgebraicAnalysisStage::ReducedSpectrum,
                    bytes: failure.bytes,
                    alignment: failure.alignment,
                    solver_invoked: false,
                }
            }
            SpectralAnalysisFailure::NumericalError => {
                RepresentationFailure::AlgebraicAnalysisFailure {
                    stage: AlgebraicAnalysisStage::ReducedSpectrum,
                    solver_invoked: false,
                }
            }
        })?;
    match spectrum.decision {
        SpectralRankDecision::Reject => Err(RepresentationFailure::ReducedPairingNotPositive {
            classification: ReducedPairingFailureClassification::RankDeficient,
            rank: spectrum.rank,
            negative_pivots: failed_cholesky_inertia
                .map(|inertia| inertia.negative)
                .unwrap_or(0),
            execution: AnalysisExecutionEvidence::pre_backend(),
        }),
        SpectralRankDecision::GrayZone => Err(RepresentationFailure::ReducedPairingGrayZone(
            AnalysisExecutionEvidence::pre_backend(),
        )),
        SpectralRankDecision::Accept => {
            if let Some(inertia) = failed_cholesky_inertia {
                if inertia.negative > 0 {
                    return Err(RepresentationFailure::ReducedPairingNotPositive {
                        classification: ReducedPairingFailureClassification::NegativeCurvature,
                        rank: spectrum.rank,
                        negative_pivots: inertia.negative,
                        execution: AnalysisExecutionEvidence::pre_backend(),
                    });
                }
                return Err(RepresentationFailure::AlgebraicAnalysisFailure {
                    stage: AlgebraicAnalysisStage::ReducedInertia,
                    solver_invoked: false,
                });
            }
            Ok(*spectrum
                .singular_values
                .last()
                .expect("a nonempty reduced pairing has a singular value"))
        }
    }
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

    pub(crate) fn usage(&self) -> &FunctionalUse {
        &self.functional
    }

    pub(crate) fn target(&self) -> f64 {
        self.target
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CanonicalEqualityParticipation {
    SolverConstraint,
    VerificationOnly,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SemanticLatentCoefficient {
    pub(crate) latent: usize,
    pub(crate) coefficient: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SemanticLatentDefinition {
    pub(crate) group_id: GroupId,
    pub(crate) field_unit: FieldUnitLabel,
    pub(crate) member_source_ids: Vec<SourceId>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardEquality {
    field: Option<FunctionalUse>,
    latent_coefficients: Vec<SemanticLatentCoefficient>,
    provenance: UsageProvenance,
    dimension: FunctionalDimension,
    target: f64,
    participation: CanonicalEqualityParticipation,
}

impl CanonicalHardEquality {
    pub(crate) fn new(
        field: Option<FunctionalUse>,
        latent_coefficients: Vec<SemanticLatentCoefficient>,
        provenance: UsageProvenance,
        dimension: FunctionalDimension,
        target: f64,
        participation: CanonicalEqualityParticipation,
    ) -> Self {
        Self {
            field,
            latent_coefficients,
            provenance,
            dimension,
            target,
            participation,
        }
    }

    fn from_field_only(equality: HardEquality) -> Self {
        let provenance = equality.functional.provenance().clone();
        let dimension = equality.functional.functional().dimension();
        Self::new(
            Some(equality.functional),
            Vec::new(),
            provenance,
            dimension,
            equality.target,
            CanonicalEqualityParticipation::SolverConstraint,
        )
    }

    pub(crate) fn field(&self) -> Option<&FunctionalUse> {
        self.field.as_ref()
    }

    pub(crate) fn latent_coefficients(&self) -> &[SemanticLatentCoefficient] {
        &self.latent_coefficients
    }

    pub(crate) fn provenance(&self) -> &UsageProvenance {
        &self.provenance
    }

    pub(crate) fn dimension(&self) -> FunctionalDimension {
        self.dimension
    }

    pub(crate) fn target(&self) -> f64 {
        self.target
    }

    pub(crate) fn participation(&self) -> CanonicalEqualityParticipation {
        self.participation
    }

    pub(crate) fn promote_to_solver_constraint(&mut self) {
        self.participation = CanonicalEqualityParticipation::SolverConstraint;
    }

    fn constant_shift_response(&self) -> f64 {
        self.field
            .as_ref()
            .map(|usage| {
                usage
                    .functional()
                    .terms()
                    .iter()
                    .map(|term| term.value_coefficient())
                    .sum::<f64>()
            })
            .unwrap_or(0.0)
            + self
                .latent_coefficients
                .iter()
                .map(|term| term.coefficient)
                .sum::<f64>()
    }

    fn evaluate(&self, field: &RecoveredCubicField, latents: &[f64]) -> f64 {
        self.field
            .as_ref()
            .map(|usage| field.evaluate_functional(usage.functional()))
            .unwrap_or(0.0)
            + self
                .latent_coefficients
                .iter()
                .map(|term| term.coefficient * latents[term.latent])
                .sum::<f64>()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CanonicalSoftLoss {
    QuadraticPenalty {
        weight: f64,
    },
    StandardDeviation {
        standard_deviation: f64,
    },
    Covariance {
        dimension: usize,
        covariance: Vec<f64>,
        whitening: Vec<f64>,
        inverse_whitening: Vec<f64>,
        precision: Vec<f64>,
    },
}

impl CanonicalSoftLoss {
    pub(crate) fn covariance(dimension: usize, covariance: Vec<f64>) -> Self {
        let (whitening, inverse_whitening, precision) =
            covariance_transforms(dimension, &covariance)
                .unwrap_or_else(|| (Vec::new(), Vec::new(), Vec::new()));
        Self::Covariance {
            dimension,
            covariance,
            whitening,
            inverse_whitening,
            precision,
        }
    }

    pub(crate) fn precision_matrix(&self, dimension: usize) -> Vec<f64> {
        match self {
            Self::QuadraticPenalty { weight } => diagonal_matrix(dimension, *weight),
            Self::StandardDeviation { standard_deviation } => {
                diagonal_matrix(dimension, 1.0 / standard_deviation.powi(2))
            }
            Self::Covariance { precision, .. } => precision.clone(),
        }
    }

    fn whitening_matrix(&self, dimension: usize) -> Vec<f64> {
        match self {
            Self::QuadraticPenalty { weight } => diagonal_matrix(dimension, weight.sqrt()),
            Self::StandardDeviation { standard_deviation } => {
                diagonal_matrix(dimension, 1.0 / standard_deviation)
            }
            Self::Covariance { whitening, .. } => whitening.clone(),
        }
    }

    fn inverse_whitening_matrix(&self, dimension: usize) -> Vec<f64> {
        match self {
            Self::QuadraticPenalty { weight } => diagonal_matrix(dimension, 1.0 / weight.sqrt()),
            Self::StandardDeviation { standard_deviation } => {
                diagonal_matrix(dimension, *standard_deviation)
            }
            Self::Covariance {
                inverse_whitening, ..
            } => inverse_whitening.clone(),
        }
    }

    fn is_valid(&self, dimension: usize) -> bool {
        match self {
            Self::QuadraticPenalty { weight } => weight.is_finite() && *weight > 0.0,
            Self::StandardDeviation { standard_deviation } => {
                standard_deviation.is_finite() && *standard_deviation > 0.0
            }
            Self::Covariance {
                dimension: covariance_dimension,
                covariance,
                whitening,
                inverse_whitening,
                precision,
            } => {
                *covariance_dimension == dimension
                    && dimension > 0
                    && [covariance, whitening, inverse_whitening, precision]
                        .into_iter()
                        .all(|matrix| {
                            matrix.len() == dimension * dimension
                                && matrix.iter().all(|entry| entry.is_finite())
                        })
            }
        }
    }

    fn residual_reference_scale(&self) -> f64 {
        match self {
            Self::QuadraticPenalty { weight } => 1.0 / weight.sqrt(),
            Self::StandardDeviation { standard_deviation } => *standard_deviation,
            Self::Covariance {
                dimension,
                covariance,
                ..
            } => (0..*dimension)
                .map(|index| covariance[index * dimension + index].sqrt())
                .fold(0.0_f64, f64::max),
        }
    }

    pub(crate) fn covariance_entries(&self) -> Option<(usize, &[f64])> {
        match self {
            Self::Covariance {
                dimension,
                covariance,
                ..
            } => Some((*dimension, covariance)),
            _ => None,
        }
    }
}

fn diagonal_matrix(dimension: usize, diagonal: f64) -> Vec<f64> {
    (0..dimension * dimension)
        .map(|index| {
            if index / dimension == index % dimension {
                diagonal
            } else {
                0.0
            }
        })
        .collect()
}

fn covariance_transforms(
    dimension: usize,
    covariance: &[f64],
) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    if dimension == 0 || covariance.len() != dimension * dimension {
        return None;
    }
    let scale = covariance
        .iter()
        .map(|entry| entry.abs())
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let root_scale = scale.sqrt();
    let mut lower = vec![0.0; dimension * dimension];
    for row in 0..dimension {
        for column in 0..=row {
            let product = (0..column)
                .map(|index| lower[row * dimension + index] * lower[column * dimension + index])
                .sum::<f64>();
            let remainder = covariance[row * dimension + column] / scale - product;
            let value = if row == column {
                if remainder <= 0.0 {
                    return None;
                }
                remainder.sqrt()
            } else {
                remainder / lower[column * dimension + column]
            };
            if !value.is_finite() {
                return None;
            }
            lower[row * dimension + column] = value;
        }
    }
    let inverse_whitening = lower
        .iter()
        .map(|entry| root_scale * entry)
        .collect::<Vec<_>>();
    let mut whitening = vec![0.0; dimension * dimension];
    for column in 0..dimension {
        for row in 0..dimension {
            let rhs = if row == column { 1.0 } else { 0.0 };
            let prior = (0..row)
                .map(|index| lower[row * dimension + index] * whitening[index * dimension + column])
                .sum::<f64>();
            whitening[row * dimension + column] =
                (rhs / root_scale - prior) / lower[row * dimension + row];
        }
    }
    let precision = (0..dimension * dimension)
        .map(|index| {
            let row = index / dimension;
            let column = index % dimension;
            (0..dimension)
                .map(|inner| {
                    whitening[inner * dimension + row] * whitening[inner * dimension + column]
                })
                .sum::<f64>()
        })
        .collect::<Vec<_>>();
    if whitening
        .iter()
        .chain(&inverse_whitening)
        .chain(&precision)
        .all(|entry| entry.is_finite())
    {
        Some((whitening, inverse_whitening, precision))
    } else {
        None
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftEquality {
    field: FunctionalUse,
    provenance: UsageProvenance,
    dimension: FunctionalDimension,
    target: f64,
}

impl CanonicalSoftEquality {
    pub(crate) fn new(field: FunctionalUse, target: f64) -> Self {
        let provenance = field.provenance().clone();
        let dimension = field.functional().dimension();
        Self {
            field,
            provenance,
            dimension,
            target,
        }
    }

    pub(crate) fn field(&self) -> &FunctionalUse {
        &self.field
    }

    pub(crate) fn provenance(&self) -> &UsageProvenance {
        &self.provenance
    }

    pub(crate) fn dimension(&self) -> FunctionalDimension {
        self.dimension
    }

    pub(crate) fn target(&self) -> f64 {
        self.target
    }

    fn evaluate(&self, field: &RecoveredCubicField) -> f64 {
        field.evaluate_functional(self.field.functional())
    }

    fn constant_shift_response(&self) -> f64 {
        self.field
            .functional()
            .terms()
            .iter()
            .map(|term| term.value_coefficient())
            .sum()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftObjective {
    residuals: Vec<ResidualId>,
    loss: CanonicalSoftLoss,
    covariance_group: Option<GroupId>,
}

impl CanonicalSoftObjective {
    pub(crate) fn new(residual: ResidualId, loss: CanonicalSoftLoss) -> Self {
        Self::new_block(vec![residual], loss, None)
    }

    pub(crate) fn new_block(
        residuals: Vec<ResidualId>,
        loss: CanonicalSoftLoss,
        covariance_group: Option<GroupId>,
    ) -> Self {
        Self {
            residuals,
            loss,
            covariance_group,
        }
    }

    pub(crate) fn residuals(&self) -> &[ResidualId] {
        &self.residuals
    }

    pub(crate) fn loss(&self) -> &CanonicalSoftLoss {
        &self.loss
    }

    pub(crate) fn covariance_group(&self) -> Option<&GroupId> {
        self.covariance_group.as_ref()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CubicCanonicalProblem {
    pub(crate) equalities: Vec<CanonicalHardEquality>,
    pub(crate) soft_equalities: Vec<CanonicalSoftEquality>,
    pub(crate) soft_objectives: Vec<CanonicalSoftObjective>,
    pub(crate) semantic_latents: Vec<SemanticLatentDefinition>,
    pub(crate) field_energy_normalization: FieldEnergyNormalization,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct EqualityAssemblyEvidence {
    pub(crate) primal_variables: usize,
    pub(crate) field_coefficients: usize,
    pub(crate) polynomial_coefficients: usize,
    pub(crate) semantic_latents: usize,
    pub(crate) side_conditions: usize,
    pub(crate) hard_equalities: usize,
    pub(crate) canonical_hard_equalities: usize,
    hard_equality_rows: Vec<AssembledHardEqualityRow>,
    soft_objective_blocks: Vec<AssembledSoftObjectiveBlock>,
}

#[derive(Debug, Clone, PartialEq)]
struct AssembledHardEqualityRow {
    kkt_equality_row: usize,
    canonical_index: usize,
    solver_index: usize,
    provenance: UsageProvenance,
    derived_block: DerivedBlockId,
    residual: ResidualId,
    derived_row: DerivedRowId,
    derived_column: Option<DerivedColumnId>,
    standard_jacobian_row: Vec<f64>,
    rhs: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct AssembledSoftObjectiveBlock {
    objective_index: usize,
    canonical_indices: Vec<usize>,
    fitting_functional_indices: Vec<usize>,
    provenances: Vec<UsageProvenance>,
    residuals: Vec<ResidualId>,
    targets: Vec<f64>,
    canonical_precision: Vec<f64>,
    standard_precision: Vec<f64>,
    whitening: Vec<f64>,
    inverse_whitening: Vec<f64>,
    covariance_group: Option<GroupId>,
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

    pub(crate) fn evaluate_functional(&self, functional: &CanonicalFunctional) -> f64 {
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
    pub(crate) provenance: UsageProvenance,
    pub(crate) dimension: FunctionalDimension,
    pub(crate) target: f64,
    pub(crate) value: f64,
    pub(crate) residual: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveredSoftEquality {
    pub(crate) provenance: UsageProvenance,
    pub(crate) dimension: FunctionalDimension,
    pub(crate) target: f64,
    pub(crate) value: f64,
    pub(crate) residual: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveredSoftObjective {
    pub(crate) canonical_indices: Vec<usize>,
    pub(crate) loss: CanonicalSoftLoss,
    pub(crate) covariance_group: Option<GroupId>,
    pub(crate) whitened_residual: Vec<f64>,
    pub(crate) objective_contribution: f64,
    pub(crate) whitening_round_trip_error: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveredSemanticLatent {
    pub(crate) group_id: GroupId,
    pub(crate) field_unit: FieldUnitLabel,
    pub(crate) member_source_ids: Vec<SourceId>,
    pub(crate) value: f64,
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
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CanonicalRelationToleranceEvidence {
    pub(crate) dimension: FunctionalDimension,
    pub(crate) characteristic_scale: f64,
    pub(crate) relation_reference_scale: f64,
    pub(crate) physical_tolerance: f64,
    pub(crate) standard_tolerance: f64,
    pub(crate) scaled_kkt_tolerance: Option<f64>,
    pub(crate) recovered_physical_tolerance: f64,
    pub(crate) round_trip_error: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct CanonicalRelationTolerancePlan {
    dimension: FunctionalDimension,
    characteristic_scale: f64,
    relation_reference_scale: f64,
    physical_tolerance: f64,
    standard_tolerance: f64,
    kkt_row: Option<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PhysicalSideConditionEvidence {
    pub(crate) components: [f64; 4],
    pub(crate) physical_tolerances: [f64; 4],
    pub(crate) standard_components: [f64; 4],
    pub(crate) recovered_standard_components: [f64; 4],
    pub(crate) round_trip_error: f64,
}

impl PhysicalSideConditionEvidence {
    fn is_within_policy(self) -> bool {
        self.components
            .into_iter()
            .zip(self.physical_tolerances)
            .all(|(component, tolerance)| component.abs() <= tolerance)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct CubicEqualitySolution {
    pub(crate) representation: CpdEvidence,
    pub(crate) assembly: EqualityAssemblyEvidence,
    pub(crate) backend: KktSolveEvidence,
    pub(crate) field: RecoveredCubicField,
    pub(crate) semantic_latents: Vec<RecoveredSemanticLatent>,
    pub(crate) hard_equalities: Vec<RecoveredHardEquality>,
    pub(crate) soft_equalities: Vec<RecoveredSoftEquality>,
    pub(crate) soft_objectives: Vec<RecoveredSoftObjective>,
    pub(crate) side_condition: PhysicalSideConditionEvidence,
    pub(crate) hard_equality_violations: FunctionalViolationEnvelope,
    pub(crate) relation_tolerances: Vec<CanonicalRelationToleranceEvidence>,
    pub(crate) tolerance_round_trip_error: f64,
    pub(crate) polynomial_round_trip_error: f64,
    pub(crate) field_coefficient_round_trip_error: f64,
    pub(crate) field_energy_round_trip_error: f64,
    pub(crate) whitening_round_trip_error: f64,
    pub(crate) objective_round_trip_error: f64,
    pub(crate) objective_verified: bool,
    pub(crate) recovery_finite: bool,
    pub(crate) provenance_verified: bool,
    pub(crate) semantic_latent_count: usize,
    pub(crate) field_energy: f64,
    pub(crate) total_objective: f64,
}

/// A physical Recover-and-Verify rejection reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RecoveryVerificationFailureReason {
    /// A stored coordinate or scaling recovery map was invalid.
    InvalidRecoveryMap,
    /// Canonical usage provenance no longer matched assembly evidence.
    ProvenanceMismatch,
    /// At least one recovered physical quantity was non-finite.
    NonFiniteRecoveredQuantity,
    /// The recovered Cubic Π₁ side condition exceeded its tolerance.
    SideConditionViolation,
    /// The side-condition recovery map exceeded its round-trip limit.
    SideConditionRoundTripViolation,
    /// At least one recovered hard equality exceeded physical tolerance.
    HardEqualityViolation,
    /// Polynomial recovery exceeded its round-trip limit.
    PolynomialRoundTripViolation,
    /// Field-coefficient recovery exceeded its round-trip limit.
    FieldCoefficientRoundTripViolation,
    /// FieldEnergy recovery exceeded its round-trip limit.
    FieldEnergyRoundTripViolation,
    /// A whitening map failed to recover its original physical residual.
    WhiteningRoundTripViolation,
    /// The independently recovered physical objective disagreed with the
    /// standard-form objective.
    ObjectiveRoundTripViolation,
    /// Relation-tolerance recovery exceeded its round-trip limit.
    ToleranceRoundTripViolation,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveryVerificationFailureEvidence {
    pub(crate) reasons: Vec<RecoveryVerificationFailureReason>,
    pub(crate) side_condition: Option<PhysicalSideConditionEvidence>,
    pub(crate) hard_equalities: Option<Vec<RecoveredHardEquality>>,
    pub(crate) soft_equalities: Option<Vec<RecoveredSoftEquality>>,
    pub(crate) soft_objectives: Option<Vec<RecoveredSoftObjective>>,
    pub(crate) relation_tolerances: Option<Vec<CanonicalRelationToleranceEvidence>>,
    pub(crate) hard_equality_violations: Option<FunctionalViolationEnvelope>,
    pub(crate) polynomial_round_trip_error: Option<f64>,
    pub(crate) field_coefficient_round_trip_error: Option<f64>,
    pub(crate) field_energy_round_trip_error: Option<f64>,
    pub(crate) whitening_round_trip_error: Option<f64>,
    pub(crate) objective_round_trip_error: Option<f64>,
    pub(crate) tolerance_round_trip_error: Option<f64>,
    pub(crate) recovery_finite: Option<bool>,
    pub(crate) provenance_verified: Option<bool>,
    pub(crate) no_model_produced: bool,
}

impl RecoveryVerificationFailureEvidence {
    fn early(reason: RecoveryVerificationFailureReason) -> Self {
        Self {
            reasons: vec![reason],
            side_condition: None,
            hard_equalities: None,
            soft_equalities: None,
            soft_objectives: None,
            relation_tolerances: None,
            hard_equality_violations: None,
            polynomial_round_trip_error: None,
            field_coefficient_round_trip_error: None,
            field_energy_round_trip_error: None,
            whitening_round_trip_error: None,
            objective_round_trip_error: None,
            tolerance_round_trip_error: None,
            recovery_finite: None,
            provenance_verified: None,
            no_model_produced: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CubicEqualityFailure {
    EmptyEqualitySet,
    NonFiniteTarget {
        equality: usize,
    },
    Representation(Box<RepresentationFailure>),
    Backend {
        failure: Box<KktFailure>,
        representation: Box<CpdEvidence>,
    },
    RecoveryVerification {
        evidence: Box<RecoveryVerificationFailureEvidence>,
        representation: Box<CpdEvidence>,
        backend: Box<KktSolveEvidence>,
    },
}

pub(crate) struct CubicEqualityCore;

#[cfg(test)]
thread_local! {
    static INJECTED_KKT_FAILURE: RefCell<Option<KktFailure>> = const { RefCell::new(None) };
}

#[cfg(test)]
pub(crate) fn inject_kkt_failure_once(failure: KktFailure) {
    INJECTED_KKT_FAILURE.with(|slot| {
        assert!(
            slot.replace(Some(failure)).is_none(),
            "only one KKT failure may be injected per fit"
        );
    });
}

#[cfg(test)]
fn take_injected_kkt_failure() -> Option<KktFailure> {
    INJECTED_KKT_FAILURE.with(|slot| slot.borrow_mut().take())
}

impl CubicEqualityCore {
    pub(crate) fn solve(
        equalities: Vec<HardEquality>,
        metric: GlobalAnisotropyMetric,
    ) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
        Self::solve_canonical(
            CubicCanonicalProblem {
                equalities: equalities
                    .into_iter()
                    .map(CanonicalHardEquality::from_field_only)
                    .collect(),
                soft_equalities: Vec::new(),
                soft_objectives: Vec::new(),
                semantic_latents: Vec::new(),
                field_energy_normalization: FieldEnergyNormalization::all_hard(),
            },
            metric,
        )
    }

    pub(crate) fn solve_canonical(
        problem: CubicCanonicalProblem,
        metric: GlobalAnisotropyMetric,
    ) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
        if problem.equalities.is_empty() && problem.soft_equalities.is_empty() {
            return Err(CubicEqualityFailure::EmptyEqualitySet);
        }
        for (index, equality) in problem.equalities.iter().enumerate() {
            if !equality.target.is_finite()
                || equality.latent_coefficients.iter().any(|term| {
                    !term.coefficient.is_finite() || term.latent >= problem.semantic_latents.len()
                })
            {
                return Err(CubicEqualityFailure::NonFiniteTarget { equality: index });
            }
        }
        let objective_residuals = problem
            .soft_objectives
            .iter()
            .flat_map(|objective| objective.residuals.iter())
            .collect::<Vec<_>>();
        if problem
            .soft_equalities
            .iter()
            .any(|relation| !relation.target.is_finite())
            || objective_residuals.len() != problem.soft_equalities.len()
            || problem
                .soft_equalities
                .iter()
                .zip(objective_residuals)
                .any(|(relation, residual)| relation.provenance.residual() != residual)
            || problem
                .soft_objectives
                .iter()
                .any(|objective| !objective.loss.is_valid(objective.residuals.len()))
        {
            return Err(CubicEqualityFailure::NonFiniteTarget {
                equality: problem.equalities.len(),
            });
        }
        let fitting_uses = canonical_fitting_uses(&problem.equalities, &problem.soft_equalities);
        let mut representation = CubicRepresentation::new(fitting_uses, metric)
            .map_err(|failure| CubicEqualityFailure::Representation(Box::new(failure)))?;
        representation.set_field_energy_normalization(problem.field_energy_normalization);
        #[cfg(test)]
        if let Some(failure) = take_injected_kkt_failure() {
            return Err(CubicEqualityFailure::Backend {
                failure: Box::new(failure),
                representation: Box::new(representation.evidence.clone()),
            });
        }
        let (assembly, backend) = solve_standard_form(&representation, &problem)?;
        recover_and_verify(representation, problem, assembly, backend)
    }
}

fn solve_standard_form(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
) -> Result<(EqualityAssemblyEvidence, KktSolveEvidence), CubicEqualityFailure> {
    let solver_equalities = problem
        .equalities
        .iter()
        .enumerate()
        .filter(|(_, equality)| {
            equality.participation == CanonicalEqualityParticipation::SolverConstraint
        })
        .collect::<Vec<_>>();
    let coefficient_count = representation.kernel.rows;
    let latent_offset = coefficient_count + POLYNOMIAL_DIMENSION;
    let primal_variables = latent_offset + problem.semantic_latents.len();
    let equality_constraints = POLYNOMIAL_DIMENSION + solver_equalities.len();
    let mut canonical_offset = 0;
    let soft_objective_blocks = problem
        .soft_objectives
        .iter()
        .enumerate()
        .map(|(objective_index, objective)| {
            let dimension = objective.residuals.len();
            let canonical_indices =
                (canonical_offset..canonical_offset + dimension).collect::<Vec<_>>();
            canonical_offset += dimension;
            let relations = canonical_indices
                .iter()
                .map(|index| &problem.soft_equalities[*index])
                .collect::<Vec<_>>();
            let canonical_precision = objective.loss.precision_matrix(dimension);
            AssembledSoftObjectiveBlock {
                objective_index,
                canonical_indices,
                fitting_functional_indices: relations
                    .iter()
                    .map(|relation| {
                        representation
                            .fitting_uses
                            .iter()
                            .position(|usage| usage.functional() == relation.field.functional())
                            .expect("every soft field functional enters the representer span")
                    })
                    .collect(),
                provenances: relations
                    .iter()
                    .map(|relation| relation.provenance.clone())
                    .collect(),
                residuals: objective.residuals.clone(),
                targets: relations.iter().map(|relation| relation.target).collect(),
                standard_precision: canonical_precision.clone(),
                canonical_precision,
                whitening: objective.loss.whitening_matrix(dimension),
                inverse_whitening: objective.loss.inverse_whitening_matrix(dimension),
                covariance_group: objective.covariance_group.clone(),
            }
        })
        .collect::<Vec<_>>();
    let field_energy_scale = representation.field_energy_normalization.factor()
        / representation.coordinates.length().powi(3);
    let hessian = DenseMatrix::from_fn(primal_variables, primal_variables, |row, column| {
        let field_energy = if row < coefficient_count && column < coefficient_count {
            field_energy_scale * representation.kernel.get(row, column)
        } else {
            0.0
        };
        field_energy
            + soft_objective_blocks
                .iter()
                .map(|objective| {
                    let dimension = objective.canonical_indices.len();
                    (0..dimension)
                        .flat_map(|left| {
                            (0..dimension).map(move |right| {
                                objective.standard_precision[left * dimension + right]
                                    * standard_field_affine_coefficient(
                                        representation,
                                        row,
                                        objective.fitting_functional_indices[left],
                                    )
                                    * standard_field_affine_coefficient(
                                        representation,
                                        column,
                                        objective.fitting_functional_indices[right],
                                    )
                            })
                        })
                        .sum::<f64>()
                })
                .sum::<f64>()
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
                let (_, equality) = solver_equalities[row - POLYNOMIAL_DIMENSION];
                standard_affine_coefficient(
                    representation,
                    column,
                    equality.field.as_ref(),
                    &equality.latent_coefficients,
                )
            }
        });
    let stationarity_rhs = (0..primal_variables)
        .map(|column| {
            soft_objective_blocks
                .iter()
                .map(|objective| {
                    let dimension = objective.canonical_indices.len();
                    (0..dimension)
                        .flat_map(|left| {
                            (0..dimension).map(move |right| {
                                standard_field_affine_coefficient(
                                    representation,
                                    column,
                                    objective.fitting_functional_indices[left],
                                ) * objective.standard_precision[left * dimension + right]
                                    * objective.targets[right]
                            })
                        })
                        .sum::<f64>()
                })
                .sum()
        })
        .collect::<Vec<_>>();
    let equality_rhs = std::iter::repeat_n(0.0, POLYNOMIAL_DIMENSION)
        .chain(
            solver_equalities
                .iter()
                .map(|(_, equality)| equality.target),
        )
        .collect::<Vec<_>>();
    let hard_equality_rows = solver_equalities
        .iter()
        .enumerate()
        .map(|(solver_index, (canonical_index, equality))| {
            let kkt_equality_row = POLYNOMIAL_DIMENSION + solver_index;
            AssembledHardEqualityRow {
                kkt_equality_row,
                canonical_index: *canonical_index,
                solver_index,
                provenance: equality.provenance.clone(),
                derived_block: DerivedBlockId::from_residual(equality.provenance.residual()),
                residual: equality.provenance.residual().clone(),
                derived_row: DerivedRowId::from_residual(equality.provenance.residual()),
                derived_column: equality
                    .field
                    .as_ref()
                    .map(|_| DerivedColumnId::from_residual(equality.provenance.residual())),
                standard_jacobian_row: (0..primal_variables)
                    .map(|column| equality_jacobian.get(kkt_equality_row, column))
                    .collect(),
                rhs: equality_rhs[kkt_equality_row],
            }
        })
        .collect();
    let backend = solve_equality_kkt(&EqualityKktSystem {
        primal_variables,
        equality_constraints,
        hessian: hessian.values(),
        equality_jacobian: equality_jacobian.values(),
        stationarity_rhs: &stationarity_rhs,
        equality_rhs: &equality_rhs,
    })
    .map_err(|failure| CubicEqualityFailure::Backend {
        failure: Box::new(failure),
        representation: Box::new(representation.evidence.clone()),
    })?;
    Ok((
        EqualityAssemblyEvidence {
            primal_variables,
            field_coefficients: coefficient_count,
            polynomial_coefficients: POLYNOMIAL_DIMENSION,
            semantic_latents: problem.semantic_latents.len(),
            side_conditions: POLYNOMIAL_DIMENSION,
            hard_equalities: solver_equalities.len(),
            canonical_hard_equalities: problem.equalities.len(),
            hard_equality_rows,
            soft_objective_blocks,
        },
        backend,
    ))
}

fn standard_affine_row(
    representation: &CubicRepresentation,
    primal_variables: usize,
    field: Option<&FunctionalUse>,
    latent_coefficients: &[SemanticLatentCoefficient],
) -> Vec<f64> {
    (0..primal_variables)
        .map(|column| {
            standard_affine_coefficient(representation, column, field, latent_coefficients)
        })
        .collect()
}

fn standard_affine_coefficient(
    representation: &CubicRepresentation,
    column: usize,
    field: Option<&FunctionalUse>,
    latent_coefficients: &[SemanticLatentCoefficient],
) -> f64 {
    let coefficient_count = representation.kernel.rows;
    let latent_offset = coefficient_count + POLYNOMIAL_DIMENSION;
    if column < latent_offset {
        field
            .map(|field| {
                let functional = representation
                    .fitting_uses
                    .iter()
                    .position(|use_| use_.functional() == field.functional())
                    .expect("every objective or hard field functional enters the representer span");
                standard_field_affine_coefficient(representation, column, functional)
            })
            .unwrap_or(0.0)
    } else {
        latent_coefficients
            .iter()
            .find(|term| term.latent == column - latent_offset)
            .map(|term| term.coefficient)
            .unwrap_or(0.0)
    }
}

fn standard_field_affine_coefficient(
    representation: &CubicRepresentation,
    column: usize,
    fitting_functional_index: usize,
) -> f64 {
    let coefficient_count = representation.kernel.rows;
    if column < coefficient_count {
        representation.kernel.get(fitting_functional_index, column)
    } else if column < coefficient_count + POLYNOMIAL_DIMENSION {
        representation
            .polynomial
            .get(fitting_functional_index, column - coefficient_count)
    } else {
        0.0
    }
}

fn standard_field_affine_value(
    representation: &CubicRepresentation,
    fitting_functional_index: usize,
    candidate: &[f64],
) -> f64 {
    candidate
        .iter()
        .enumerate()
        .map(|(column, value)| {
            standard_field_affine_coefficient(representation, column, fitting_functional_index)
                * value
        })
        .sum()
}

fn verifies_assembled_provenance_and_rows(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    assembly: &EqualityAssemblyEvidence,
    backend: &KktSolveEvidence,
) -> bool {
    let solver_equalities = problem
        .equalities
        .iter()
        .enumerate()
        .filter(|(_, equality)| {
            equality.participation == CanonicalEqualityParticipation::SolverConstraint
        })
        .collect::<Vec<_>>();
    if assembly.hard_equality_rows.len() != solver_equalities.len()
        || assembly.soft_objective_blocks.len() != problem.soft_objectives.len()
        || backend.equality_multipliers.len() != assembly.side_conditions + assembly.hard_equalities
    {
        return false;
    }
    let hard_rows_verified = assembly
        .hard_equality_rows
        .iter()
        .zip(solver_equalities)
        .enumerate()
        .all(|(solver_index, (row, (canonical_index, equality)))| {
            let expected = expected_hard_equality_row(representation, problem, equality);
            row.canonical_index == canonical_index
                && row.solver_index == solver_index
                && row.kkt_equality_row == POLYNOMIAL_DIMENSION + solver_index
                && row.provenance == equality.provenance
                && row.derived_block
                    == DerivedBlockId::from_residual(equality.provenance.residual())
                && row.residual == *equality.provenance.residual()
                && row.derived_row == DerivedRowId::from_residual(equality.provenance.residual())
                && row.derived_column
                    == equality
                        .field
                        .as_ref()
                        .map(|_| DerivedColumnId::from_residual(equality.provenance.residual()))
                && row.rhs == equality.target
                && row.standard_jacobian_row == expected
        });
    let soft_rows_verified = assembly
        .soft_objective_blocks
        .iter()
        .zip(&problem.soft_objectives)
        .enumerate()
        .all(|(objective_index, (block, objective))| {
            let dimension = objective.residuals.len();
            block.objective_index == objective_index
                && block.canonical_indices.len() == dimension
                && block.fitting_functional_indices.len() == dimension
                && block.provenances.len() == dimension
                && block.residuals == objective.residuals
                && block.targets.len() == dimension
                && block.covariance_group == objective.covariance_group
                && block.canonical_precision == objective.loss.precision_matrix(dimension)
                && block.standard_precision == block.canonical_precision
                && block.whitening == objective.loss.whitening_matrix(dimension)
                && block.inverse_whitening == objective.loss.inverse_whitening_matrix(dimension)
                && block
                    .canonical_indices
                    .iter()
                    .enumerate()
                    .all(|(component, canonical_index)| {
                        problem
                            .soft_equalities
                            .get(*canonical_index)
                            .is_some_and(|relation| {
                                representation
                                    .fitting_uses
                                    .get(block.fitting_functional_indices[component])
                                    .is_some_and(|usage| {
                                        usage.functional() == relation.field.functional()
                                    })
                                    && block.provenances[component] == relation.provenance
                                    && block.residuals[component] == *relation.provenance.residual()
                                    && block.targets[component] == relation.target
                            })
                    })
        });
    hard_rows_verified && soft_rows_verified
}

fn expected_hard_equality_row(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    equality: &CanonicalHardEquality,
) -> Vec<f64> {
    standard_affine_row(
        representation,
        representation.kernel.rows + POLYNOMIAL_DIMENSION + problem.semantic_latents.len(),
        equality.field.as_ref(),
        &equality.latent_coefficients,
    )
}

fn recover_and_verify(
    representation: CubicRepresentation,
    problem: CubicCanonicalProblem,
    assembly: EqualityAssemblyEvidence,
    backend: KktSolveEvidence,
) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
    if !representation.coordinates.is_valid_recovery_map() {
        return Err(CubicEqualityFailure::RecoveryVerification {
            evidence: Box::new(RecoveryVerificationFailureEvidence::early(
                RecoveryVerificationFailureReason::InvalidRecoveryMap,
            )),
            representation: Box::new(representation.evidence.clone()),
            backend: Box::new(backend),
        });
    }
    let provenance_verified =
        verifies_assembled_provenance_and_rows(&representation, &problem, &assembly, &backend);
    if !provenance_verified {
        return Err(CubicEqualityFailure::RecoveryVerification {
            evidence: Box::new(RecoveryVerificationFailureEvidence::early(
                RecoveryVerificationFailureReason::ProvenanceMismatch,
            )),
            representation: Box::new(representation.evidence.clone()),
            backend: Box::new(backend),
        });
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
    let latent_offset = coefficient_count + POLYNOMIAL_DIMENSION;
    let latent_values =
        backend.candidate[latent_offset..latent_offset + problem.semantic_latents.len()].to_vec();
    let semantic_latents = problem
        .semantic_latents
        .iter()
        .zip(&latent_values)
        .map(|(definition, value)| RecoveredSemanticLatent {
            group_id: definition.group_id.clone(),
            field_unit: definition.field_unit.clone(),
            member_source_ids: definition.member_source_ids.clone(),
            value: *value,
        })
        .collect::<Vec<_>>();
    let hard_equalities = problem
        .equalities
        .iter()
        .map(|equality| {
            let value = equality.evaluate(&field, &latent_values);
            RecoveredHardEquality {
                provenance: equality.provenance.clone(),
                dimension: equality.dimension,
                target: equality.target,
                value,
                residual: value - equality.target,
            }
        })
        .collect::<Vec<_>>();
    let soft_equalities = problem
        .soft_equalities
        .iter()
        .map(|relation| {
            let value = relation.evaluate(&field);
            let residual = value - relation.target;
            RecoveredSoftEquality {
                provenance: relation.provenance.clone(),
                dimension: relation.dimension,
                target: relation.target,
                value,
                residual,
            }
        })
        .collect::<Vec<_>>();
    let soft_objectives = problem
        .soft_objectives
        .iter()
        .zip(&assembly.soft_objective_blocks)
        .map(|(objective, block)| {
            let residual = block
                .canonical_indices
                .iter()
                .map(|index| soft_equalities[*index].residual)
                .collect::<Vec<_>>();
            let dimension = residual.len();
            let whitened_residual = matrix_vector_product(&block.whitening, dimension, &residual);
            let recovered_residual =
                matrix_vector_product(&block.inverse_whitening, dimension, &whitened_residual);
            let whitening_round_trip_error = relative_slice_error(&recovered_residual, &residual);
            RecoveredSoftObjective {
                canonical_indices: block.canonical_indices.clone(),
                loss: objective.loss.clone(),
                covariance_group: objective.covariance_group.clone(),
                objective_contribution: 0.5
                    * whitened_residual
                        .iter()
                        .map(|component| component * component)
                        .sum::<f64>(),
                whitened_residual,
                whitening_round_trip_error,
            }
        })
        .collect::<Vec<_>>();
    let standard_side_components = std::array::from_fn(|column| {
        (0..coefficient_count)
            .map(|row| representation.polynomial.get(row, column) * standard_coefficients[row])
            .sum::<f64>()
    });
    let mapped_physical_side = representation
        .coordinates
        .to_physical_side_condition(standard_side_components);
    let physical_side_components = std::array::from_fn(|column| {
        representation
            .fitting_uses
            .iter()
            .zip(&field.coefficients)
            .map(|(usage, coefficient)| {
                let pairing = usage.functional().evaluate_affine(
                    if column == 0 { 1.0 } else { 0.0 },
                    std::array::from_fn(|axis| if axis + 1 == column { 1.0 } else { 0.0 }),
                );
                pairing * coefficient
            })
            .sum::<f64>()
    });
    let recovered_standard_side = representation
        .coordinates
        .to_standard_side_condition(physical_side_components);
    let side_condition_round_trip_error =
        relative_slice_error(&mapped_physical_side, &physical_side_components).max(
            relative_slice_error(&recovered_standard_side, &standard_side_components),
        );
    let side_condition = PhysicalSideConditionEvidence {
        components: physical_side_components,
        physical_tolerances: representation
            .coordinates
            .to_physical_side_condition_tolerances(
                [EQUALITY_KKT_POLICY_V1.side_condition_limit; POLYNOMIAL_DIMENSION],
            ),
        standard_components: standard_side_components,
        recovered_standard_components: recovered_standard_side,
        round_trip_error: side_condition_round_trip_error,
    };
    let hard_equality_violations = FunctionalViolationEnvelope::from_dimensioned_residuals(
        hard_equalities
            .iter()
            .map(|equality| (equality.dimension, equality.residual)),
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
    let field_energy = representation.field_energy_normalization.factor() * native_cubic_energy;
    let standard_energy = dot_product(
        standard_coefficients,
        &representation.kernel.multiply_vector(standard_coefficients),
    );
    let recovered_energy = representation.field_energy_normalization.factor() * standard_energy
        / representation.coordinates.length().powi(3);
    let field_energy_round_trip_error =
        (field_energy - recovered_energy).abs() / recovered_energy.abs().max(1.0);
    let soft_loss = soft_objectives
        .iter()
        .map(|objective| objective.objective_contribution)
        .sum::<f64>();
    let total_objective = 0.5 * field_energy + soft_loss;
    let standard_soft_loss = assembly
        .soft_objective_blocks
        .iter()
        .map(|objective| {
            let dimension = objective.canonical_indices.len();
            let residual = objective
                .fitting_functional_indices
                .iter()
                .zip(&objective.targets)
                .map(|(functional, target)| {
                    standard_field_affine_value(&representation, *functional, &backend.candidate)
                        - target
                })
                .collect::<Vec<_>>();
            let weighted =
                matrix_vector_product(&objective.standard_precision, dimension, &residual);
            0.5 * residual
                .iter()
                .zip(weighted)
                .map(|(left, right)| left * right)
                .sum::<f64>()
        })
        .sum::<f64>();
    let standard_total_objective = 0.5 * recovered_energy + standard_soft_loss;
    let whitening_round_trip_error = soft_objectives
        .iter()
        .map(|objective| objective.whitening_round_trip_error)
        .fold(0.0_f64, f64::max);
    let objective_round_trip_error = (total_objective - standard_total_objective).abs()
        / standard_total_objective.abs().max(1.0);
    let objective_verified =
        objective_round_trip_error <= EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit;
    let relation_tolerances =
        canonical_relation_tolerances(&representation, &problem, field_energy, &assembly, &backend);
    let tolerance_round_trip_error = relation_tolerances
        .iter()
        .map(|evidence| evidence.round_trip_error)
        .fold(0.0_f64, f64::max);
    let recovery_finite = field.coefficients.iter().all(|value| value.is_finite())
        && field
            .physical_polynomial
            .iter()
            .all(|value| value.is_finite())
        && hard_equalities
            .iter()
            .all(|equality| equality.value.is_finite() && equality.residual.is_finite())
        && soft_equalities
            .iter()
            .all(|equality| equality.value.is_finite() && equality.residual.is_finite())
        && soft_objectives.iter().all(|objective| {
            objective.objective_contribution.is_finite()
                && objective.whitening_round_trip_error.is_finite()
                && objective
                    .whitened_residual
                    .iter()
                    .all(|component| component.is_finite())
        })
        && semantic_latents
            .iter()
            .all(|latent| latent.value.is_finite())
        && side_condition
            .components
            .into_iter()
            .chain(side_condition.physical_tolerances)
            .chain(side_condition.standard_components)
            .chain(side_condition.recovered_standard_components)
            .chain([side_condition.round_trip_error])
            .all(f64::is_finite)
        && polynomial_round_trip_error.is_finite()
        && field_coefficient_round_trip_error.is_finite()
        && field_energy.is_finite()
        && field_energy_round_trip_error.is_finite()
        && objective_round_trip_error.is_finite()
        && relation_tolerances.iter().all(|evidence| {
            evidence.characteristic_scale.is_finite()
                && evidence.relation_reference_scale.is_finite()
                && evidence.physical_tolerance.is_finite()
                && evidence.standard_tolerance.is_finite()
                && evidence.scaled_kkt_tolerance.is_none_or(f64::is_finite)
                && evidence.recovered_physical_tolerance.is_finite()
                && evidence.round_trip_error.is_finite()
        })
        && total_objective.is_finite();

    let mut reasons = Vec::new();
    if !recovery_finite {
        reasons.push(RecoveryVerificationFailureReason::NonFiniteRecoveredQuantity);
    }
    if !side_condition.is_within_policy() {
        reasons.push(RecoveryVerificationFailureReason::SideConditionViolation);
    }
    if side_condition.round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::SideConditionRoundTripViolation);
    }
    if hard_equalities
        .iter()
        .zip(&relation_tolerances)
        .any(|(equality, tolerance)| equality.residual.abs() > tolerance.physical_tolerance)
    {
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
    if whitening_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::WhiteningRoundTripViolation);
    }
    if !objective_verified {
        reasons.push(RecoveryVerificationFailureReason::ObjectiveRoundTripViolation);
    }
    if tolerance_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::ToleranceRoundTripViolation);
    }
    if !reasons.is_empty() {
        return Err(CubicEqualityFailure::RecoveryVerification {
            evidence: Box::new(RecoveryVerificationFailureEvidence {
                reasons,
                side_condition: Some(side_condition),
                hard_equalities: Some(hard_equalities),
                soft_equalities: Some(soft_equalities),
                soft_objectives: Some(soft_objectives),
                relation_tolerances: Some(relation_tolerances),
                hard_equality_violations: Some(hard_equality_violations),
                polynomial_round_trip_error: Some(polynomial_round_trip_error),
                field_coefficient_round_trip_error: Some(field_coefficient_round_trip_error),
                field_energy_round_trip_error: Some(field_energy_round_trip_error),
                whitening_round_trip_error: Some(whitening_round_trip_error),
                objective_round_trip_error: Some(objective_round_trip_error),
                tolerance_round_trip_error: Some(tolerance_round_trip_error),
                recovery_finite: Some(recovery_finite),
                provenance_verified: Some(provenance_verified),
                no_model_produced: true,
            }),
            representation: Box::new(representation.evidence.clone()),
            backend: Box::new(backend),
        });
    }

    Ok(CubicEqualitySolution {
        representation: representation.evidence.clone(),
        assembly,
        backend,
        field,
        semantic_latents,
        hard_equalities,
        soft_equalities,
        soft_objectives,
        side_condition,
        hard_equality_violations,
        relation_tolerances,
        tolerance_round_trip_error,
        polynomial_round_trip_error,
        field_coefficient_round_trip_error,
        field_energy_round_trip_error,
        whitening_round_trip_error,
        objective_round_trip_error,
        objective_verified,
        recovery_finite,
        provenance_verified,
        semantic_latent_count: problem.semantic_latents.len(),
        field_energy,
        total_objective,
    })
}

fn canonical_relation_tolerances(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    field_energy: f64,
    assembly: &EqualityAssemblyEvidence,
    backend: &KktSolveEvidence,
) -> Vec<CanonicalRelationToleranceEvidence> {
    let field_value_gauge_offset = canonical_gauge_offset(problem, FunctionalDimension::FieldValue);
    let derivative_gauge_offset =
        canonical_gauge_offset(problem, FunctionalDimension::FieldValuePerLength);
    let value_implied_field_scale = problem
        .equalities
        .iter()
        .filter(|equality| equality.dimension == FunctionalDimension::FieldValue)
        .map(|equality| {
            (equality.target - equality.constant_shift_response() * field_value_gauge_offset).abs()
        })
        .chain(
            problem
                .soft_equalities
                .iter()
                .filter(|equality| equality.dimension == FunctionalDimension::FieldValue)
                .map(|equality| {
                    (equality.target
                        - equality.constant_shift_response() * field_value_gauge_offset)
                        .abs()
                }),
        )
        .fold(0.0_f64, f64::max);
    let derivative_implied_field_scale = representation.coordinates.length()
        * problem
            .equalities
            .iter()
            .filter(|equality| equality.dimension == FunctionalDimension::FieldValuePerLength)
            .map(|equality| equality.target.abs())
            .chain(
                problem
                    .soft_equalities
                    .iter()
                    .filter(|equality| {
                        equality.dimension == FunctionalDimension::FieldValuePerLength
                    })
                    .map(|equality| equality.target.abs()),
            )
            .fold(0.0_f64, f64::max);
    let mut soft_offset = 0;
    let soft_loss_implied_field_scale = problem
        .soft_objectives
        .iter()
        .map(|objective| {
            let relation = &problem.soft_equalities[soft_offset];
            soft_offset += objective.residuals.len();
            match relation.dimension {
                FunctionalDimension::FieldValue => objective.loss.residual_reference_scale(),
                FunctionalDimension::FieldValuePerLength => {
                    representation.coordinates.length() * objective.loss.residual_reference_scale()
                }
            }
        })
        .fold(0.0_f64, f64::max);
    let native_field_energy = field_energy / problem.field_energy_normalization.factor();
    let field_scale = (native_field_energy.abs() * representation.coordinates.length().powi(3))
        .sqrt()
        .max(value_implied_field_scale)
        .max(derivative_implied_field_scale)
        .max(soft_loss_implied_field_scale);
    let mut standard_by_kkt_row = vec![0.0; backend.scaling.cumulative_exponents.len()];
    let mut tolerance_plans = problem
        .equalities
        .iter()
        .enumerate()
        .map(|(index, equality)| {
            let dimension = equality.dimension;
            let (characteristic_scale, gauge_offset) = tolerance_scales_for_dimension(
                dimension,
                field_scale,
                representation.coordinates.length(),
                field_value_gauge_offset,
                derivative_gauge_offset,
            );
            let relation_reference_scale =
                (equality.target - equality.constant_shift_response() * gauge_offset).abs();
            let physical_tolerance = EQUALITY_KKT_POLICY_V1
                .canonical_characteristic_tolerance_multiplier
                * characteristic_scale
                + EQUALITY_KKT_POLICY_V1.canonical_relation_reference_tolerance_multiplier
                    * relation_reference_scale;
            let standard_tolerance = representation
                .coordinates
                .to_standard_tolerance(dimension, physical_tolerance);
            let kkt_row = assembly
                .hard_equality_rows
                .iter()
                .find(|row| row.canonical_index == index)
                .map(|row| assembly.primal_variables + row.kkt_equality_row);
            if let Some(kkt_row) = kkt_row {
                standard_by_kkt_row[kkt_row] = standard_tolerance;
            }
            CanonicalRelationTolerancePlan {
                dimension,
                characteristic_scale,
                relation_reference_scale,
                physical_tolerance,
                standard_tolerance,
                kkt_row,
            }
        })
        .collect::<Vec<_>>();
    let scaled = backend
        .scaling
        .scale_residual_or_tolerance(&standard_by_kkt_row);
    let recovered = backend.scaling.recover_residual_or_tolerance(&scaled);
    tolerance_plans
        .drain(..)
        .map(|plan| {
            let scaled_kkt_tolerance = plan.kkt_row.map(|kkt_row| scaled[kkt_row]);
            let recovered_standard_tolerance = plan
                .kkt_row
                .map(|kkt_row| recovered[kkt_row])
                .unwrap_or(plan.standard_tolerance);
            let recovered_physical_tolerance = representation
                .coordinates
                .to_physical_tolerance(plan.dimension, recovered_standard_tolerance);
            CanonicalRelationToleranceEvidence {
                dimension: plan.dimension,
                characteristic_scale: plan.characteristic_scale,
                relation_reference_scale: plan.relation_reference_scale,
                physical_tolerance: plan.physical_tolerance,
                standard_tolerance: plan.standard_tolerance,
                scaled_kkt_tolerance,
                recovered_physical_tolerance,
                round_trip_error: (recovered_physical_tolerance - plan.physical_tolerance).abs()
                    / plan.physical_tolerance.abs().max(1.0),
            }
        })
        .collect()
}

fn tolerance_scales_for_dimension(
    dimension: FunctionalDimension,
    field_scale: f64,
    length: f64,
    field_value_gauge_offset: f64,
    derivative_gauge_offset: f64,
) -> (f64, f64) {
    match dimension {
        FunctionalDimension::FieldValue => (field_scale, field_value_gauge_offset),
        FunctionalDimension::FieldValuePerLength => (field_scale / length, derivative_gauge_offset),
    }
}

fn canonical_gauge_offset(problem: &CubicCanonicalProblem, dimension: FunctionalDimension) -> f64 {
    let response_targets = problem
        .equalities
        .iter()
        .filter(|equality| equality.dimension == dimension)
        .map(|equality| (equality.constant_shift_response(), equality.target))
        .chain(
            problem
                .soft_equalities
                .iter()
                .filter(|equality| equality.dimension == dimension)
                .map(|equality| (equality.constant_shift_response(), equality.target)),
        )
        .collect::<Vec<_>>();
    let response_scale = response_targets
        .iter()
        .map(|(response, _)| response.abs())
        .fold(0.0_f64, f64::max);
    let target_scale = response_targets
        .iter()
        .map(|(_, target)| target.abs())
        .fold(0.0_f64, f64::max);
    if response_scale == 0.0 || target_scale == 0.0 {
        return 0.0;
    }
    let numerator = response_targets
        .iter()
        .map(|(response, target)| (response / response_scale) * (target / target_scale))
        .sum::<f64>();
    let denominator = response_targets
        .iter()
        .map(|(response, _)| (response / response_scale).powi(2))
        .sum::<f64>();
    target_scale / response_scale * numerator / denominator
}

fn relative_slice_error(actual: &[f64], expected: &[f64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(0.0_f64, f64::max)
}

fn matrix_vector_product(matrix: &[f64], dimension: usize, vector: &[f64]) -> Vec<f64> {
    debug_assert_eq!(matrix.len(), dimension * dimension);
    debug_assert_eq!(vector.len(), dimension);
    (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix[row * dimension + column] * vector[column])
                .sum()
        })
        .collect()
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
                        ResidualId::new(format!("residual-{index}")),
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

    fn canonical_problem(uses: Vec<FunctionalUse>, targets: Vec<f64>) -> CubicCanonicalProblem {
        CubicCanonicalProblem {
            equalities: uses
                .into_iter()
                .zip(targets)
                .map(|(usage, target)| {
                    CanonicalHardEquality::from_field_only(HardEquality::new(usage, target))
                })
                .collect(),
            soft_equalities: Vec::new(),
            soft_objectives: Vec::new(),
            semantic_latents: Vec::new(),
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
        }
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
        assert!(evidence.reduced_smallest_singular_value > 0.0);
        assert!(evidence.reduced_symmetry_defect <= evidence.symmetry_defect_limit);
        assert!(evidence.affine_reproduction_error <= 1.0e-11);
        let normalization = representation.field_energy_normalization();
        assert_eq!(normalization.factor(), 1.0);
        let transformed = FieldEnergyNormalization::try_new(
            normalization.factor() * 2.0_f64.powi(3) / 4.0_f64.powi(2),
        )
        .expect("finite positive unit scales preserve normalization validity");
        assert_close(transformed.factor(), 0.5, 1.0e-15);
        assert_close(
            transformed.factor() * 6.0,
            normalization.factor() * 3.0,
            1.0e-15,
        );

        let standard = representation
            .solve_coordinate_transform()
            .to_standard(TRUTH_POLYNOMIAL);
        let round_trip = representation
            .solve_coordinate_transform()
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
        let standard = transform
            .to_standard_functional(&physical)
            .expect("the manufactured functional has a finite standard form");
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
            assert_eq!(&recovered.provenance, expected_usage.provenance());
            assert_eq!(recovered.dimension, FunctionalDimension::FieldValue);
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
        assert!(solution.side_condition.is_within_policy());
        assert!(solution.side_condition.round_trip_error <= 1.0e-11);
        assert!(solution.hard_equality_violations.field_value <= 1.0e-8);
        assert_eq!(
            solution.hard_equality_violations.field_value_per_length,
            0.0
        );
        assert!(solution.polynomial_round_trip_error <= 1.0e-11);
        assert!(solution.field_coefficient_round_trip_error <= 1.0e-11);
        assert!(solution.field_energy_round_trip_error <= 1.0e-11);
        assert_eq!(
            solution.relation_tolerances.len(),
            MANUFACTURED_TARGETS.len()
        );
        assert!(solution.relation_tolerances.iter().all(|tolerance| {
            assert_close(
                tolerance.physical_tolerance,
                1.0e-10 * tolerance.characteristic_scale
                    + 1.0e-8 * tolerance.relation_reference_scale,
                1.0e-15,
            );
            tolerance.standard_tolerance == tolerance.physical_tolerance
                && tolerance.recovered_physical_tolerance == tolerance.physical_tolerance
                && tolerance.round_trip_error <= 1.0e-11
        }));
        assert!(solution.tolerance_round_trip_error <= 1.0e-11);
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
    fn semantic_latent_is_solved_in_the_kkt_and_verified_from_the_candidate() {
        let group_id = GroupId::new("manufactured-shared-level");
        let mut uses = usages();
        let member = uses.remove(0);
        let member_provenance = member.provenance().clone();
        let mut equalities = vec![CanonicalHardEquality::new(
            Some(member),
            vec![SemanticLatentCoefficient {
                latent: 0,
                coefficient: -1.0,
            }],
            member_provenance,
            FunctionalDimension::FieldValue,
            0.0,
            CanonicalEqualityParticipation::SolverConstraint,
        )];
        equalities.extend(
            uses.into_iter()
                .zip(MANUFACTURED_TARGETS[1..].iter().copied())
                .map(|(usage, target)| {
                    CanonicalHardEquality::from_field_only(HardEquality::new(usage, target))
                }),
        );
        let gauge_provenance = UsageProvenance::new(
            SourceId::new("manufactured-latent-gauge"),
            Some(group_id.clone()),
            RelationId::new("manufactured-latent-gauge/relation"),
            ResidualId::new("manufactured-latent-gauge/residual"),
            SemanticRolePath::new("additive-field-gauge/level-set"),
        );
        equalities.push(CanonicalHardEquality::new(
            None,
            vec![SemanticLatentCoefficient {
                latent: 0,
                coefficient: 1.0,
            }],
            gauge_provenance.clone(),
            FunctionalDimension::FieldValue,
            MANUFACTURED_TARGETS[0],
            CanonicalEqualityParticipation::SolverConstraint,
        ));
        let problem = CubicCanonicalProblem {
            equalities: equalities.clone(),
            soft_equalities: Vec::new(),
            soft_objectives: Vec::new(),
            semantic_latents: vec![SemanticLatentDefinition {
                group_id: group_id.clone(),
                field_unit: FieldUnitLabel::new("manufactured-unit"),
                member_source_ids: vec![SourceId::new("issue-17-manufactured-0")],
            }],
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
        };

        let solution =
            CubicEqualityCore::solve_canonical(problem.clone(), GlobalAnisotropyMetric::identity())
                .expect("the explicit semantic-latent KKT should recover");
        assert_eq!(solution.assembly.primal_variables, 15);
        assert_eq!(solution.assembly.semantic_latents, 1);
        assert_eq!(solution.assembly.hard_equalities, 11);
        assert_eq!(solution.assembly.canonical_hard_equalities, 11);
        assert_eq!(solution.backend.capacity.kkt_dimension, 30);
        assert_eq!(solution.semantic_latent_count, 1);
        assert_eq!(solution.semantic_latents[0].group_id, group_id);
        assert_eq!(
            solution.semantic_latents[0].field_unit.as_str(),
            "manufactured-unit"
        );
        assert_close(
            solution.semantic_latents[0].value,
            MANUFACTURED_TARGETS[0],
            1.0e-8,
        );
        assert_eq!(
            solution.hard_equalities.last().unwrap().provenance,
            gauge_provenance
        );

        let mut inconsistent = problem;
        let conflicting_provenance = UsageProvenance::new(
            SourceId::new("manufactured-conflicting-gauge"),
            Some(GroupId::new("manufactured-shared-level")),
            RelationId::new("manufactured-conflicting-gauge/relation"),
            ResidualId::new("manufactured-conflicting-gauge/residual"),
            SemanticRolePath::new("additive-field-gauge/level-set"),
        );
        inconsistent.equalities.push(CanonicalHardEquality::new(
            None,
            vec![SemanticLatentCoefficient {
                latent: 0,
                coefficient: 1.0,
            }],
            conflicting_provenance,
            FunctionalDimension::FieldValue,
            MANUFACTURED_TARGETS[0] + 1.0,
            CanonicalEqualityParticipation::VerificationOnly,
        ));
        let failure =
            CubicEqualityCore::solve_canonical(inconsistent, GlobalAnisotropyMetric::identity())
                .expect_err("an incompatible verification-only gauge must reject the candidate");
        match failure {
            CubicEqualityFailure::RecoveryVerification { evidence, .. } => {
                assert!(
                    evidence
                        .reasons
                        .contains(&RecoveryVerificationFailureReason::HardEqualityViolation)
                );
                assert_eq!(evidence.hard_equalities.as_ref().unwrap().len(), 12);
                assert_eq!(evidence.relation_tolerances.as_ref().unwrap().len(), 12);
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
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
                rank: None,
                mode: VerifiedCanonicalMode {
                    residual: 0.0,
                    execution: AnalysisExecutionEvidence::pre_backend(),
                },
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

        let coordinates = CubicSolveCoordinateTransform {
            center: [0.0; 3],
            length: 1.0,
            degenerate_extent: false,
        };
        let failure = verify_polynomial_rank(&polynomial, &[], &coordinates)
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
    fn capacity_preflight_polynomial_analysis_retains_gray_zone_failure() {
        let functionals = [
            functional([0.0; 3], 1.0, [0.0; 3]),
            functional([0.0; 3], 0.0, [1.0, 0.0, 0.0]),
            functional([0.0; 3], 0.0, [0.0, 1.0, 0.0]),
            functional([0.0; 3], 0.0, [0.0, 0.0, 1.0e-13]),
        ];
        let fitting_uses = functionals
            .into_iter()
            .enumerate()
            .map(|(index, functional)| {
                FunctionalUse::new(
                    functional,
                    UsageProvenance::new(
                        SourceId::new(format!("gray-{index}")),
                        None,
                        RelationId::new(format!("gray-relation-{index}")),
                        ResidualId::new(format!("gray-residual-{index}")),
                        SemanticRolePath::new(format!("gray/{index}")),
                    ),
                )
            })
            .collect::<Vec<_>>();

        assert!(matches!(
            preflight_polynomial_analysis_failure(&fitting_uses),
            Some(RepresentationFailure::PolynomialRankGrayZone { .. })
        ));
    }

    #[test]
    fn exact_zero_column_does_not_invent_the_rank_of_the_remaining_pairing() {
        let functionals = [0.0, 1.0, 2.0]
            .into_iter()
            .map(|coordinate| functional([coordinate, coordinate, 0.0], 1.0, [0.0; 3]))
            .collect::<Vec<_>>();
        let polynomial = DenseMatrix::from_fn(3, POLYNOMIAL_DIMENSION, |row, column| {
            let coordinate = row as f64;
            match column {
                0 => 1.0,
                1 | 2 => coordinate,
                3 => 0.0,
                _ => unreachable!(),
            }
        });
        let coordinates = CubicSolveCoordinateTransform {
            center: [0.0; 3],
            length: 1.0,
            degenerate_extent: false,
        };

        let failure = verify_polynomial_rank(&polynomial, &functionals, &coordinates)
            .expect_err("the exact zero column proves a missing affine mode");
        assert!(matches!(
            failure,
            RepresentationFailure::PolynomialRankDeficient { rank: None, .. }
        ));
    }

    #[test]
    fn reduced_spd_fallback_distinguishes_negative_rank_and_gray_evidence() {
        let negative = DenseMatrix::from_fn(2, 2, |row, column| {
            if row == column {
                if row == 0 { 1.0 } else { -1.0 }
            } else {
                0.0
            }
        });
        match verify_reduced_pairing(&negative)
            .expect_err("negative curvature must fail after the Cholesky-first path")
        {
            RepresentationFailure::ReducedPairingNotPositive {
                classification,
                rank,
                negative_pivots,
                ..
            } => {
                assert_eq!(
                    classification,
                    ReducedPairingFailureClassification::NegativeCurvature
                );
                assert_eq!(rank, 2);
                assert_eq!(negative_pivots, 1);
            }
            other => panic!("unexpected negative-curvature failure: {other:?}"),
        }

        let two_by_two_pivot =
            DenseMatrix::from_fn(2, 2, |row, column| if row != column { 1.0 } else { 0.0 });
        match verify_reduced_pairing(&two_by_two_pivot)
            .expect_err("an indefinite two-by-two pivot must retain its inertia")
        {
            RepresentationFailure::ReducedPairingNotPositive {
                classification,
                negative_pivots,
                ..
            } => {
                assert_eq!(
                    classification,
                    ReducedPairingFailureClassification::NegativeCurvature
                );
                assert_eq!(negative_pivots, 1);
            }
            other => panic!("unexpected two-by-two inertia failure: {other:?}"),
        }

        let rank_deficient =
            DenseMatrix::from_fn(2, 2, |row, column| (row == column && row == 0) as u8 as f64);
        match verify_reduced_pairing(&rank_deficient)
            .expect_err("a zero reduced mode must be rank deficient")
        {
            RepresentationFailure::ReducedPairingNotPositive {
                classification,
                rank,
                ..
            } => {
                assert_eq!(
                    classification,
                    ReducedPairingFailureClassification::RankDeficient
                );
                assert_eq!(rank, 1);
            }
            other => panic!("unexpected rank failure: {other:?}"),
        }

        let gray = DenseMatrix::from_fn(2, 2, |row, column| {
            if row == column {
                if row == 0 { 1.0 } else { 1.0e-13 }
            } else {
                0.0
            }
        });
        assert_eq!(
            verify_reduced_pairing(&gray)
                .expect_err("a reduced mode inside the spectral band remains undecided"),
            RepresentationFailure::ReducedPairingGrayZone(AnalysisExecutionEvidence::pre_backend())
        );
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
                    ResidualId::new(format!("extreme-residual-{index}")),
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
                reason: SolveCoordinateTransformFailureReason::FieldRecoveryScaleNotInvertible,
                solver_invoked: false,
            }
        );
    }

    #[test]
    fn normalized_gradient_overflow_is_a_typed_preassembly_rejection() {
        let tiny = 1.0e-100;
        let uses = [
            [0.0, 0.0, 0.0],
            [tiny, 0.0, 0.0],
            [0.0, tiny, 0.0],
            [0.0, 0.0, tiny],
        ]
        .into_iter()
        .enumerate()
        .map(|(index, support)| {
            FunctionalUse::new(
                functional(support, 1.0, [1.0e250, 0.0, 0.0]),
                UsageProvenance::new(
                    SourceId::new(format!("overflow-{index}")),
                    None,
                    RelationId::new(format!("overflow-relation-{index}")),
                    ResidualId::new(format!("overflow-residual-{index}")),
                    SemanticRolePath::new(format!("hard-equality/{index}")),
                ),
            )
        })
        .collect();

        let failure = CubicRepresentation::new(uses, GlobalAnisotropyMetric::identity())
            .expect_err("a derived nonfinite functional must fail closed without panicking");
        assert_eq!(
            failure,
            RepresentationFailure::InvalidSolveCoordinateTransform {
                reason: SolveCoordinateTransformFailureReason::StandardFunctionalNotFinite,
                solver_invoked: false,
            }
        );
    }

    #[test]
    fn canonical_relation_tolerances_are_invariant_to_additive_field_gauge() {
        let original_uses = usages();
        let original = CubicEqualityCore::solve(
            original_uses
                .clone()
                .into_iter()
                .zip(MANUFACTURED_TARGETS)
                .map(|(usage, target)| HardEquality::new(usage, target))
                .collect(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect("the original gauge should solve");
        let gauge_shift = 137.0;
        let shifted = CubicEqualityCore::solve(
            original_uses
                .into_iter()
                .zip(MANUFACTURED_TARGETS)
                .map(|(usage, target)| {
                    let constant_response = usage
                        .functional()
                        .terms()
                        .iter()
                        .map(|term| term.value_coefficient())
                        .sum::<f64>();
                    HardEquality::new(usage, target + gauge_shift * constant_response)
                })
                .collect(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect("an additive gauge shift should preserve the physical solve");

        for (actual, expected) in shifted
            .relation_tolerances
            .iter()
            .zip(&original.relation_tolerances)
        {
            assert_close(
                actual.characteristic_scale,
                expected.characteristic_scale,
                1.0e-12,
            );
            assert_close(
                actual.relation_reference_scale,
                expected.relation_reference_scale,
                1.0e-12,
            );
            assert_close(
                actual.physical_tolerance,
                expected.physical_tolerance,
                1.0e-12,
            );
        }
    }

    #[test]
    fn physical_side_condition_recovery_handles_a_translated_solve_frame() {
        let translation = [10.0, -20.0, 5.0];
        let shifted_uses = usages()
            .into_iter()
            .map(|usage| {
                let terms = usage
                    .functional()
                    .terms()
                    .iter()
                    .map(|term| {
                        FunctionalTerm::new(
                            std::array::from_fn(|axis| term.support()[axis] + translation[axis]),
                            term.value_coefficient(),
                            term.gradient_coefficient(),
                        )
                    })
                    .collect();
                FunctionalUse::new(
                    CanonicalFunctional::new(usage.functional().dimension(), terms)
                        .expect("translation preserves a canonical functional"),
                    usage.provenance().clone(),
                )
            })
            .zip(MANUFACTURED_TARGETS)
            .map(|(usage, target)| HardEquality::new(usage, target))
            .collect();
        let solution = CubicEqualityCore::solve(shifted_uses, GlobalAnisotropyMetric::identity())
            .expect("translation should preserve the manufactured equality solve");

        assert_eq!(solution.representation.solve_coordinate_center, translation);
        assert!(solution.side_condition.is_within_policy());
        assert!(solution.side_condition.round_trip_error <= 1.0e-11);
        assert_ne!(
            solution.side_condition.physical_tolerances[1],
            EQUALITY_KKT_POLICY_V1.side_condition_limit
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
    }

    #[test]
    fn damaged_coordinate_map_is_a_recovery_failure_not_a_backend_contract_failure() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representation is valid");
        let problem = canonical_problem(usages(), MANUFACTURED_TARGETS.to_vec());
        let (assembly, backend) = solve_standard_form(&representation, &problem)
            .expect("the undamaged standard form should produce a backend candidate");
        let mut damaged = representation;
        damaged.coordinates.length = 0.0;

        let failure = recover_and_verify(damaged, problem, assembly, backend)
            .expect_err("an invalid inverse coordinate map must fail during recovery");

        match failure {
            CubicEqualityFailure::RecoveryVerification {
                evidence, backend, ..
            } => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::InvalidRecoveryMap]
                );
                assert!(evidence.no_model_produced);
                assert!(!backend.attempts.is_empty());
                assert!(backend.attempts.iter().any(|attempt| {
                    attempt.termination == crate::kkt::SolveAttemptTermination::CandidateProduced
                }));
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_provenance_map_fails_recovery_without_returning_a_partial_model() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representation is valid");
        let problem = canonical_problem(usages(), MANUFACTURED_TARGETS.to_vec());
        let (mut assembly, backend) = solve_standard_form(&representation, &problem)
            .expect("the standard form should still satisfy its backend contract");
        assembly.hard_equality_rows[0].provenance = UsageProvenance::new(
            SourceId::new("corrupted-source"),
            None,
            RelationId::new("corrupted-relation"),
            ResidualId::new("corrupted-residual"),
            SemanticRolePath::new("corrupted-role"),
        );
        let failure = recover_and_verify(representation, problem, assembly, backend)
            .expect_err("canonical provenance must round-trip before a model exists");

        match failure {
            CubicEqualityFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::ProvenanceMismatch]
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_kkt_row_target_association_fails_provenance_recovery() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representation is valid");
        let problem = canonical_problem(usages(), MANUFACTURED_TARGETS.to_vec());
        let (mut assembly, backend) = solve_standard_form(&representation, &problem)
            .expect("the standard form should produce a verified backend candidate");
        assembly.hard_equality_rows.swap(0, 1);

        let failure = recover_and_verify(representation, problem, assembly, backend)
            .expect_err("a KKT row/relation reassociation must not produce a model");
        match failure {
            CubicEqualityFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::ProvenanceMismatch]
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_soft_objective_provenance_is_rejected_even_for_an_exact_fit() {
        let hard_uses = usages();
        let soft_use = hard_uses.last().unwrap().clone();
        let soft_target = *MANUFACTURED_TARGETS.last().unwrap();
        let mut problem = canonical_problem(hard_uses, MANUFACTURED_TARGETS.to_vec());
        problem.field_energy_normalization = FieldEnergyNormalization::try_new(3.0).unwrap();
        let soft = CanonicalSoftEquality::new(soft_use, soft_target);
        problem.soft_objectives.push(CanonicalSoftObjective::new(
            soft.provenance().residual().clone(),
            CanonicalSoftLoss::QuadraticPenalty { weight: 2.0 },
        ));
        problem.soft_equalities.push(soft);
        let mut representation = CubicRepresentation::new(
            canonical_fitting_uses(&problem.equalities, &problem.soft_equalities),
            GlobalAnisotropyMetric::identity(),
        )
        .expect("the mixed hard/soft representation is valid");
        representation.set_field_energy_normalization(problem.field_energy_normalization);
        let (mut assembly, backend) = solve_standard_form(&representation, &problem)
            .expect("the undamaged objective should produce a backend candidate");
        let residual = standard_field_affine_value(
            &representation,
            assembly.soft_objective_blocks[0].fitting_functional_indices[0],
            &backend.candidate,
        ) - soft_target;
        assert!(residual.abs() <= 1.0e-11, "residual={residual:e}");
        assembly.soft_objective_blocks[0].standard_precision[0] *= 2.0;

        let failure = recover_and_verify(representation, problem, assembly, backend)
            .expect_err("a damaged objective recovery map must not produce a model");
        match failure {
            CubicEqualityFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::ProvenanceMismatch]
                );
                assert!(evidence.objective_round_trip_error.is_none());
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_whitening_recovery_map_is_rejected_without_a_model() {
        let hard_uses = usages();
        let soft_uses = hard_uses[8..].to_vec();
        let mut problem = canonical_problem(hard_uses, MANUFACTURED_TARGETS.to_vec());
        problem.field_energy_normalization = FieldEnergyNormalization::try_new(3.0).unwrap();
        let soft_equalities = soft_uses
            .into_iter()
            .zip([MANUFACTURED_TARGETS[8] + 1.0, MANUFACTURED_TARGETS[9] - 2.0])
            .map(|(usage, target)| CanonicalSoftEquality::new(usage, target))
            .collect::<Vec<_>>();
        problem
            .soft_objectives
            .push(CanonicalSoftObjective::new_block(
                soft_equalities
                    .iter()
                    .map(|equality| equality.provenance().residual().clone())
                    .collect(),
                CanonicalSoftLoss::covariance(2, vec![1.0, 0.25, 0.25, 2.0]),
                Some(GroupId::new("damaged-whitening")),
            ));
        problem.soft_equalities = soft_equalities;
        let mut representation = CubicRepresentation::new(
            canonical_fitting_uses(&problem.equalities, &problem.soft_equalities),
            GlobalAnisotropyMetric::identity(),
        )
        .expect("the mixed hard/soft representation is valid");
        representation.set_field_energy_normalization(problem.field_energy_normalization);
        let (mut assembly, backend) = solve_standard_form(&representation, &problem)
            .expect("the undamaged whitening map should produce a candidate");
        if let CanonicalSoftLoss::Covariance {
            inverse_whitening, ..
        } = &mut problem.soft_objectives[0].loss
        {
            inverse_whitening[0] *= 2.0;
        } else {
            panic!("the fixture uses covariance whitening");
        }
        assembly.soft_objective_blocks[0].inverse_whitening[0] *= 2.0;

        let failure = recover_and_verify(representation, problem, assembly, backend)
            .expect_err("a non-invertible whitening recovery map must not produce a model");
        match failure {
            CubicEqualityFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::WhiteningRoundTripViolation]
                );
                assert!(
                    evidence.whitening_round_trip_error.unwrap()
                        > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }
}
