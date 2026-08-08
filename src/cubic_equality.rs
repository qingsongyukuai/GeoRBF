use faer::Accum;
use faer::Conj;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::householder::{
    apply_block_householder_sequence_on_the_left_in_place_scratch,
    apply_block_householder_sequence_on_the_left_in_place_with_conj,
    apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj,
};
use faer::linalg::{matmul::matmul, triangular_solve::solve_lower_triangular_in_place};
use faer::prelude::*;
#[cfg(test)]
use std::cell::RefCell;

use crate::capacity::{CapacityExceededEvidence, plan_equality_capacity};
use crate::cubic::{CubicKernel, GlobalAnisotropyMetric};
use crate::cubic_solver_form::{
    AllSourceRecoveryLedger, CanonicalCubicFieldForm, CanonicalCubicSolverForm,
    CanonicalHardConflictWitness, CanonicalHardRecoveryGraph, CubicFieldCoordinateLayout,
};
use crate::faer_backend;
use crate::functional::{
    CanonicalFunctional, DerivedBlockId, DerivedColumnId, DerivedRowId, FunctionalDimension,
    FunctionalRepresenterSpan, FunctionalTerm, FunctionalUse, GroupId, ResidualId, SourceId,
    UsageProvenance,
};
use crate::geometry::FieldUnitLabel;
use crate::kernel::FieldEnergyNormalization;
use crate::kkt::{EqualityKktSystem, KktFailure, KktSolveEvidence, solve_equality_kkt};
use crate::math::dot3;
use crate::numerical::{
    EXECUTED_NUMERICAL_POLICY, SpectralAnalysisFailure, SpectralRankDecision, analyze_spectral_rank,
};
use crate::precision_rescue::{
    CertifiedDoubleDouble, DOUBLE_DOUBLE_PRECISION_BITS, DoubleDouble, PrecisionRescueConclusion,
    classify_symmetric_schur, cubic_jet_dd, cubic_pairing_dd_certified,
};

pub(crate) const POLYNOMIAL_DIMENSION: usize = 4;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct DenseMatrix {
    rows: usize,
    columns: usize,
    values: Vec<f64>,
}

impl DenseMatrix {
    pub(crate) fn from_fn(
        rows: usize,
        columns: usize,
        mut value: impl FnMut(usize, usize) -> f64,
    ) -> Self {
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

    pub(crate) fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.columns + column]
    }

    pub(crate) fn set(&mut self, row: usize, column: usize, value: f64) {
        self.values[row * self.columns + column] = value;
    }

    pub(crate) fn values(&self) -> &[f64] {
        &self.values
    }

    pub(crate) fn multiply_vector(&self, vector: &[f64]) -> Vec<f64> {
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

    pub(crate) fn to_standard_field_coefficients(self, physical: &[f64]) -> Vec<f64> {
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

    pub(crate) fn to_physical_side_condition(self, standard: [f64; 4]) -> [f64; 4] {
        let field_scale = self.length.powi(3);
        [
            standard[0] / field_scale,
            (self.center[0] * standard[0] + self.length * standard[1]) / field_scale,
            (self.center[1] * standard[0] + self.length * standard[2]) / field_scale,
            (self.center[2] * standard[0] + self.length * standard[3]) / field_scale,
        ]
    }

    pub(crate) fn to_standard_side_condition(self, physical: [f64; 4]) -> [f64; 4] {
        let field_scale = self.length.powi(3);
        let constant = field_scale * physical[0];
        [
            constant,
            (field_scale * physical[1] - self.center[0] * constant) / self.length,
            (field_scale * physical[2] - self.center[1] * constant) / self.length,
            (field_scale * physical[3] - self.center[2] * constant) / self.length,
        ]
    }

    pub(crate) fn to_physical_side_condition_tolerances(self, standard: [f64; 4]) -> [f64; 4] {
        let field_scale = self.length.powi(3);
        [
            standard[0] / field_scale,
            (self.center[0].abs() * standard[0] + self.length * standard[1]) / field_scale,
            (self.center[1].abs() * standard[0] + self.length * standard[2]) / field_scale,
            (self.center[2].abs() * standard[0] + self.length * standard[3]) / field_scale,
        ]
    }

    pub(crate) fn is_valid_recovery_map(self) -> bool {
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
    pub(crate) quotient_construction: QuotientConstructionEvidence,
    pub(crate) quotient_factorization: QuotientFactorizationEvidence,
    pub(crate) singular_values: Vec<f64>,
    pub(crate) polynomial_rrqr_ratio: f64,
    pub(crate) polynomial_svd_ratio: f64,
    pub(crate) polynomial_rank_reject_ratio: f64,
    pub(crate) polynomial_rank_accept_ratio: f64,
    pub(crate) polynomial_precision_rescue: Option<PrecisionRescueEvidence>,
    pub(crate) reduced_symmetry_defect: f64,
    pub(crate) symmetry_defect_limit: f64,
    pub(crate) reduced_largest_singular_value: f64,
    pub(crate) reduced_smallest_singular_value: f64,
    pub(crate) affine_reproduction_error: f64,
    pub(crate) solve_coordinate_center: [f64; 3],
    pub(crate) solve_coordinate_length: f64,
    pub(crate) degenerate_extent: bool,
    pub(crate) problem_regularization_applied: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RepresentationBuildStage {
    SourceParticipation,
    RepresenterAssembly,
    PolynomialPairing,
    HouseholderQuotient,
    QuotientFactorization,
    ResponseAssembly,
    Backend,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RepresentationBuildEvidence {
    pub(crate) failure_stage: RepresentationBuildStage,
    pub(crate) last_completed_stage: RepresentationBuildStage,
    pub(crate) representer_count: Option<usize>,
    pub(crate) polynomial_dimension: Option<usize>,
    pub(crate) polynomial_rank: Option<usize>,
    pub(crate) quotient_construction: Option<QuotientConstructionEvidence>,
    pub(crate) quotient_factorization: Option<QuotientFactorizationEvidence>,
    pub(crate) retained_modes: Option<usize>,
    pub(crate) truncated_modes: Option<usize>,
    pub(crate) problem_regularization_applied: bool,
}

impl RepresentationBuildEvidence {
    fn new() -> Self {
        Self {
            failure_stage: RepresentationBuildStage::RepresenterAssembly,
            last_completed_stage: RepresentationBuildStage::SourceParticipation,
            representer_count: None,
            polynomial_dimension: None,
            polynomial_rank: None,
            quotient_construction: None,
            quotient_factorization: None,
            retained_modes: None,
            truncated_modes: None,
            problem_regularization_applied: false,
        }
    }

    pub(crate) fn completed(evidence: &CpdEvidence) -> Self {
        Self {
            failure_stage: RepresentationBuildStage::Backend,
            last_completed_stage: RepresentationBuildStage::ResponseAssembly,
            representer_count: Some(evidence.fitting_functional_count),
            polynomial_dimension: Some(evidence.polynomial_dimension),
            polynomial_rank: Some(evidence.polynomial_rank),
            quotient_construction: Some(evidence.quotient_construction),
            quotient_factorization: Some(evidence.quotient_factorization.clone()),
            retained_modes: Some(evidence.quotient_factorization.retained_modes),
            truncated_modes: Some(evidence.quotient_factorization.truncated_modes),
            problem_regularization_applied: evidence.problem_regularization_applied,
        }
    }

    pub(crate) fn response_assembly_failed(evidence: &CpdEvidence) -> Self {
        let mut build = Self::completed(evidence);
        build.failure_stage = RepresentationBuildStage::ResponseAssembly;
        build.last_completed_stage = RepresentationBuildStage::QuotientFactorization;
        build
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PolynomialRankEvidence {
    pub(crate) rrqr_ratio: f64,
    pub(crate) svd_ratio: f64,
    pub(crate) reject_ratio: f64,
    pub(crate) accept_ratio: f64,
    pub(crate) backend_invoked: bool,
    pub(crate) precision_rescue: Option<PrecisionRescueEvidence>,
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
    pub(crate) precision_rescue: Option<PrecisionRescueEvidence>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RepresentationFailure {
    AuditedBuild {
        evidence: Box<RepresentationBuildEvidence>,
        failure: Box<RepresentationFailure>,
    },
    AuditedFactorization {
        evidence: Box<QuotientFactorizationEvidence>,
        failure: Box<RepresentationFailure>,
    },
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
    PolynomialNegativeCurvature {
        evidence: PrecisionRescueEvidence,
    },
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
    HouseholderOrthogonalityContract {
        observed: f64,
        limit: f64,
    },
    CanonicalResponseRoundTripContract {
        observed: f64,
        limit: f64,
    },
    QuotientPivotRequiresPrecisionRescue {
        quotient_dimension: usize,
        pivot_index: usize,
        interval: Option<OutwardRoundedInterval>,
        execution: AnalysisExecutionEvidence,
    },
    QuotientPrecisionRescueGrayZone {
        quotient_dimension: usize,
        evidence: PrecisionRescueEvidence,
        execution: AnalysisExecutionEvidence,
    },
    QuotientRankDeficient {
        quotient_dimension: usize,
        evidence: PrecisionRescueEvidence,
        mode: VerifiedCanonicalMode,
    },
    QuotientNegativeCurvature {
        quotient_dimension: usize,
        evidence: PrecisionRescueEvidence,
        execution: AnalysisExecutionEvidence,
    },
    QuotientFactorizationNotPositive {
        quotient_dimension: usize,
        pivot_index: usize,
        interval: OutwardRoundedInterval,
        execution: AnalysisExecutionEvidence,
    },
    QuotientLltContract {
        observed: f64,
        limit: f64,
    },
    QuotientFieldEnergyIdentityContract {
        observed: f64,
        limit: f64,
    },
    QuotientSideConditionContract {
        observed: f64,
        limit: f64,
    },
    QuotientRecoveryRoundTripContract {
        observed: f64,
        limit: f64,
    },
    QuotientResponseRoundTripContract {
        observed: f64,
        limit: f64,
    },
    ReducedSymmetryContract {
        observed: f64,
        limit: f64,
    },
    AffineReproductionBackend(Box<KktFailure>),
    AffineReproductionContract {
        observed: f64,
        limit: f64,
    },
}

impl RepresentationFailure {
    pub(crate) fn audited(self, evidence: RepresentationBuildEvidence) -> RepresentationFailure {
        match self {
            Self::AuditedBuild { .. } => self,
            _ => Self::AuditedBuild {
                evidence: Box::new(evidence),
                failure: Box::new(self),
            },
        }
    }

    pub(crate) fn root_cause(&self) -> &RepresentationFailure {
        match self {
            Self::AuditedBuild { failure, .. } | Self::AuditedFactorization { failure, .. } => {
                failure.root_cause()
            }
            _ => self,
        }
    }

    pub(crate) fn build_evidence(&self) -> Option<&RepresentationBuildEvidence> {
        match self {
            Self::AuditedBuild { evidence, .. } => Some(evidence),
            _ => None,
        }
    }

    fn with_factorization_evidence(
        self,
        evidence: QuotientFactorizationEvidence,
    ) -> RepresentationFailure {
        Self::AuditedFactorization {
            evidence: Box::new(evidence),
            failure: Box::new(self),
        }
    }

    fn factorization_evidence(&self) -> Option<&QuotientFactorizationEvidence> {
        match self {
            Self::AuditedBuild { failure, .. } => failure.factorization_evidence(),
            Self::AuditedFactorization { evidence, .. } => Some(evidence),
            _ => None,
        }
    }

    fn completed_quotient_modes(&self) -> Option<usize> {
        match self.root_cause() {
            Self::QuotientPivotRequiresPrecisionRescue { pivot_index, .. }
            | Self::QuotientFactorizationNotPositive { pivot_index, .. } => Some(*pivot_index),
            Self::QuotientPrecisionRescueGrayZone { evidence, .. }
            | Self::QuotientRankDeficient { evidence, .. }
            | Self::QuotientNegativeCurvature { evidence, .. } => Some(evidence.first_mode),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlgebraicAnalysisStage {
    PolynomialRank,
    ReducedCholesky,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SolveCoordinateTransformFailureReason {
    BoundingBoxCenterNotFinite,
    CharacteristicLengthNotFinite,
    FieldRecoveryScaleNotInvertible,
    StandardFunctionalNotFinite,
}

/// Cubic Quotient Representation Module.
///
/// Callers construct the representation with [`Self::build`], obtain a
/// canonical functional's solver-coordinate response with [`Self::response`],
/// and recover solver field coordinates with [`Self::recover`]. Numerical
/// matrices, Householder storage, and backend realization stay behind this
/// crate-private interface.
#[derive(Debug, Clone)]
pub(crate) struct CubicRepresentation {
    fitting_uses: Vec<FunctionalUse>,
    metric: GlobalAnisotropyMetric,
    coordinates: CubicSolveCoordinateTransform,
    kernel: DenseMatrix,
    polynomial: DenseMatrix,
    null_space: HouseholderNullSpace,
    energy_basis: EnergyOrthonormalQuotientBasis,
    evidence: CpdEvidence,
    field_energy_normalization: FieldEnergyNormalization,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CubicFunctionalResponse {
    pub(crate) standard_field: Vec<f64>,
    pub(crate) quotient_field: Vec<f64>,
    pub(crate) polynomial: [f64; POLYNOMIAL_DIMENSION],
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum CubicSolverFieldCoordinates<'a> {
    Standard(&'a [f64]),
    Quotient(&'a [f64]),
}

#[derive(Debug, Clone)]
pub(crate) struct CubicRepresentationRecovery {
    pub(crate) field: RecoveredCubicField,
    pub(crate) standard_coefficients: Vec<f64>,
    pub(crate) recovered_solver_coordinates: Vec<f64>,
    pub(crate) solver_round_trip_error: f64,
    pub(crate) side_condition: PhysicalSideConditionEvidence,
    pub(crate) polynomial_round_trip_error: f64,
    pub(crate) field_coefficient_round_trip_error: f64,
    pub(crate) field_energy: f64,
    pub(crate) recovered_energy: f64,
    pub(crate) field_energy_round_trip_error: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CubicRepresentationRecoveryFailure {
    InvalidRecoveryMap,
    Representation(RepresentationFailure),
}

impl CubicRepresentation {
    pub(crate) fn audit_response_assembly_failure(
        &self,
        failure: RepresentationFailure,
    ) -> RepresentationFailure {
        failure.audited(RepresentationBuildEvidence::response_assembly_failed(
            &self.evidence,
        ))
    }

    pub(crate) fn build(
        fitting_uses: Vec<FunctionalUse>,
        metric: GlobalAnisotropyMetric,
        field_energy_normalization: FieldEnergyNormalization,
    ) -> Result<(Self, CanonicalCubicFieldForm), RepresentationFailure> {
        let representation =
            Self::new_with_audit(fitting_uses, metric, field_energy_normalization)?;
        let field_form = representation.solver_field_form();
        Ok((representation, field_form))
    }

    fn new(
        fitting_uses: Vec<FunctionalUse>,
        metric: GlobalAnisotropyMetric,
    ) -> Result<Self, RepresentationFailure> {
        Self::new_with_normalization(fitting_uses, metric, FieldEnergyNormalization::all_hard())
    }

    fn new_with_normalization(
        fitting_uses: Vec<FunctionalUse>,
        metric: GlobalAnisotropyMetric,
        field_energy_normalization: FieldEnergyNormalization,
    ) -> Result<Self, RepresentationFailure> {
        Self::new_with_audit(fitting_uses, metric, field_energy_normalization)
            .map_err(|failure| failure.root_cause().clone())
    }

    fn new_with_audit(
        fitting_uses: Vec<FunctionalUse>,
        metric: GlobalAnisotropyMetric,
        field_energy_normalization: FieldEnergyNormalization,
    ) -> Result<Self, RepresentationFailure> {
        let mut build_evidence = RepresentationBuildEvidence::new();
        if fitting_uses.is_empty() {
            return Err(RepresentationFailure::EmptyRepresenterSpan.audited(build_evidence));
        }
        let functionals = fitting_uses
            .iter()
            .map(|usage| usage.functional().clone())
            .collect::<Vec<_>>();
        let (coordinates, standard_functionals, polynomial) =
            assemble_polynomial_pairing(&functionals)
                .map_err(|failure| failure.audited(build_evidence.clone()))?;
        build_evidence.representer_count = Some(functionals.len());
        build_evidence.failure_stage = RepresentationBuildStage::PolynomialPairing;
        build_evidence.last_completed_stage = RepresentationBuildStage::RepresenterAssembly;
        build_evidence.polynomial_dimension = Some(polynomial.columns);
        let (singular_values, polynomial_rank, polynomial_rank_evidence) =
            verify_polynomial_rank(&polynomial, &functionals, &coordinates)
                .map_err(|failure| failure.audited(build_evidence.clone()))?;
        build_evidence.failure_stage = RepresentationBuildStage::HouseholderQuotient;
        build_evidence.last_completed_stage = RepresentationBuildStage::PolynomialPairing;
        build_evidence.polynomial_rank = Some(polynomial_rank);
        let kernel = assemble_kernel_pairing(&standard_functionals, &metric);
        let null_space = HouseholderNullSpace::new(&polynomial, polynomial_rank)
            .map_err(|failure| failure.audited(build_evidence.clone()))?;
        let (mut reduced, trailing_polynomial, canonical_responses, quotient_construction) =
            implicit_quotient_congruence(&kernel, &polynomial, &null_space)
                .map_err(|failure| failure.audited(build_evidence.clone()))?;
        build_evidence.quotient_construction = Some(quotient_construction);
        let null_space_defect = quotient_construction.null_space_defect;
        if null_space_defect > EXECUTED_NUMERICAL_POLICY.null_space_defect_limit {
            return Err(RepresentationFailure::NullSpaceDefect {
                observed: null_space_defect,
                limit: EXECUTED_NUMERICAL_POLICY.null_space_defect_limit,
            }
            .audited(build_evidence));
        }
        if quotient_construction.householder_orthogonality_error
            > EXECUTED_NUMERICAL_POLICY.quotient_householder_orthogonality_limit
        {
            return Err(RepresentationFailure::HouseholderOrthogonalityContract {
                observed: quotient_construction.householder_orthogonality_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_householder_orthogonality_limit,
            }
            .audited(build_evidence));
        }
        if quotient_construction.canonical_response_round_trip_error
            > EXECUTED_NUMERICAL_POLICY.quotient_canonical_response_round_trip_limit
        {
            return Err(RepresentationFailure::CanonicalResponseRoundTripContract {
                observed: quotient_construction.canonical_response_round_trip_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_canonical_response_round_trip_limit,
            }
            .audited(build_evidence));
        }
        let reduced_symmetry_defect = normalized_symmetry_defect(&reduced);
        let symmetry_defect_limit = EXECUTED_NUMERICAL_POLICY.reduced_symmetry_multiplier
            * f64::EPSILON
            * reduced.rows.max(reduced.columns) as f64;
        if reduced_symmetry_defect > symmetry_defect_limit {
            return Err(RepresentationFailure::ReducedSymmetryContract {
                observed: reduced_symmetry_defect,
                limit: symmetry_defect_limit,
            }
            .audited(build_evidence));
        }
        build_evidence.failure_stage = RepresentationBuildStage::QuotientFactorization;
        build_evidence.last_completed_stage = RepresentationBuildStage::HouseholderQuotient;
        symmetrize(&mut reduced);
        let field_energy_scale = field_energy_normalization.factor() / coordinates.length().powi(3);
        let quotient_energy_gram =
            DenseMatrix::from_fn(reduced.rows, reduced.columns, |row, column| {
                field_energy_scale * reduced.get(row, column)
            });
        let rescue_source = CanonicalPrecisionRescueSource {
            functionals: &standard_functionals,
            metric: &metric,
            null_space: &null_space,
            field_energy_scale,
        };
        let energy_basis = EnergyOrthonormalQuotientBasis::factor(
            &quotient_energy_gram,
            &trailing_polynomial,
            &canonical_responses,
            Some(&rescue_source),
        )
        .map_err(|failure| {
            if let Some(evidence) = failure.factorization_evidence() {
                build_evidence.quotient_factorization = Some(evidence.clone());
                build_evidence.retained_modes = Some(evidence.retained_modes);
                build_evidence.truncated_modes = Some(evidence.truncated_modes);
            } else if let Some(completed_modes) = failure.completed_quotient_modes() {
                build_evidence.retained_modes = Some(completed_modes);
                build_evidence.truncated_modes = Some(0);
            }
            failure.audited(build_evidence.clone())
        })?;
        let quotient_factorization = energy_basis.evidence.clone();
        build_evidence.failure_stage = RepresentationBuildStage::ResponseAssembly;
        build_evidence.last_completed_stage = RepresentationBuildStage::QuotientFactorization;
        build_evidence.quotient_factorization = Some(quotient_factorization.clone());
        build_evidence.retained_modes = Some(quotient_factorization.retained_modes);
        build_evidence.truncated_modes = Some(quotient_factorization.truncated_modes);
        let (reduced_largest_singular_value, reduced_smallest_singular_value) =
            quotient_mode_risk_estimates(&reduced, &energy_basis, field_energy_scale);
        let affine_reproduction_error = affine_reproduction_error(&kernel, &polynomial)
            .map_err(|failure| failure.audited(build_evidence.clone()))?;
        if affine_reproduction_error > EXECUTED_NUMERICAL_POLICY.affine_reproduction_limit {
            return Err(RepresentationFailure::AffineReproductionContract {
                observed: affine_reproduction_error,
                limit: EXECUTED_NUMERICAL_POLICY.affine_reproduction_limit,
            }
            .audited(build_evidence));
        }
        Ok(Self {
            fitting_uses,
            metric,
            coordinates,
            kernel,
            polynomial,
            null_space,
            energy_basis,
            evidence: CpdEvidence {
                fitting_functional_count: functionals.len(),
                polynomial_dimension: POLYNOMIAL_DIMENSION,
                polynomial_rank,
                quotient_construction,
                quotient_factorization,
                singular_values,
                polynomial_rrqr_ratio: polynomial_rank_evidence.rrqr_ratio,
                polynomial_svd_ratio: polynomial_rank_evidence.svd_ratio,
                polynomial_rank_reject_ratio: polynomial_rank_evidence.reject_ratio,
                polynomial_rank_accept_ratio: polynomial_rank_evidence.accept_ratio,
                polynomial_precision_rescue: polynomial_rank_evidence.precision_rescue,
                reduced_symmetry_defect,
                symmetry_defect_limit,
                reduced_largest_singular_value,
                reduced_smallest_singular_value,
                affine_reproduction_error,
                solve_coordinate_center: coordinates.center(),
                solve_coordinate_length: coordinates.length(),
                degenerate_extent: coordinates.degenerate_extent(),
                problem_regularization_applied: false,
            },
            field_energy_normalization,
        })
    }

    fn standard_functional_row(
        &self,
        functional: &CanonicalFunctional,
    ) -> Result<(Vec<f64>, [f64; POLYNOMIAL_DIMENSION]), RepresentationFailure> {
        let standard = self.coordinates.to_standard_functional(functional)?;
        let ambient = self
            .fitting_uses
            .iter()
            .map(|usage| {
                let basis = self
                    .coordinates
                    .to_standard_functional(usage.functional())
                    .expect("a retained basis functional has a finite standard form");
                CubicKernel::pairing(&standard, &basis, &self.metric)
            })
            .collect();
        let polynomial = std::array::from_fn(|column| {
            standard.evaluate_affine(
                if column == 0 { 1.0 } else { 0.0 },
                std::array::from_fn(|axis| if axis + 1 == column { 1.0 } else { 0.0 }),
            )
        });
        Ok((ambient, polynomial))
    }

    fn solver_field_form(&self) -> CanonicalCubicFieldForm {
        let standard_field_variables = self.fitting_uses.len();
        let quotient_field_variables = self.null_space.reduced_dimension();
        let field_energy_scale =
            self.field_energy_normalization.factor() / self.coordinates.length().powi(3);
        let standard_field_energy = (0..standard_field_variables)
            .flat_map(|row| {
                (0..standard_field_variables)
                    .map(move |column| field_energy_scale * self.kernel.get(row, column))
            })
            .collect();
        let quotient_field_energy = (0..quotient_field_variables)
            .flat_map(|row| {
                (0..quotient_field_variables).map(move |column| f64::from(row == column))
            })
            .collect();
        let standard_side_conditions = (0..POLYNOMIAL_DIMENSION)
            .flat_map(|polynomial| {
                (0..standard_field_variables)
                    .map(move |field| self.polynomial.get(field, polynomial))
            })
            .collect();
        CanonicalCubicFieldForm::new(
            standard_field_variables,
            quotient_field_variables,
            POLYNOMIAL_DIMENSION,
            standard_field_energy,
            quotient_field_energy,
            standard_side_conditions,
            self.coordinates.length(),
            self.evidence.clone(),
        )
    }

    pub(crate) fn response(
        &self,
        functional: &CanonicalFunctional,
    ) -> Result<CubicFunctionalResponse, RepresentationFailure> {
        #[cfg(test)]
        if take_injected_response_assembly_failure() {
            return Err(RepresentationFailure::QuotientResponseRoundTripContract {
                observed: f64::INFINITY,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_basis_response_round_trip_limit,
            });
        }
        let (standard_field, polynomial) = self.standard_functional_row(functional)?;
        let quotient_response = self.null_space.project(&standard_field)?;
        let quotient_field = self.energy_basis.response(&quotient_response)?;
        Ok(CubicFunctionalResponse {
            standard_field,
            quotient_field,
            polynomial,
        })
    }

    pub(crate) fn recover(
        &self,
        solver_coordinates: CubicSolverFieldCoordinates<'_>,
        standard_polynomial: [f64; POLYNOMIAL_DIMENSION],
    ) -> Result<CubicRepresentationRecovery, CubicRepresentationRecoveryFailure> {
        if !self.coordinates.is_valid_recovery_map() {
            return Err(CubicRepresentationRecoveryFailure::InvalidRecoveryMap);
        }
        let (standard_coefficients, recovered_solver_coordinates, solver_round_trip_error) =
            match solver_coordinates {
                CubicSolverFieldCoordinates::Standard(coefficients) => {
                    (coefficients.to_vec(), coefficients.to_vec(), 0.0)
                }
                CubicSolverFieldCoordinates::Quotient(coefficients) => {
                    let reduced = self
                        .energy_basis
                        .to_householder_coordinates(coefficients)
                        .map_err(CubicRepresentationRecoveryFailure::Representation)?;
                    let standard = self
                        .null_space
                        .expand(&reduced)
                        .map_err(CubicRepresentationRecoveryFailure::Representation)?;
                    let recovered_reduced = self
                        .null_space
                        .project(&standard)
                        .map_err(CubicRepresentationRecoveryFailure::Representation)?;
                    let recovered = self
                        .energy_basis
                        .to_solver_coordinates(&recovered_reduced)
                        .map_err(CubicRepresentationRecoveryFailure::Representation)?;
                    let round_trip_error = relative_slice_error(&recovered, coefficients);
                    (standard, recovered, round_trip_error)
                }
            };
        let field = RecoveredCubicField::from_standard_candidate(
            self,
            &standard_coefficients,
            standard_polynomial,
        );
        let standard_side_components = std::array::from_fn(|column| {
            (0..standard_coefficients.len())
                .map(|row| self.polynomial.get(row, column) * standard_coefficients[row])
                .sum::<f64>()
        });
        let mapped_physical_side = self
            .coordinates
            .to_physical_side_condition(standard_side_components);
        let physical_side_components = std::array::from_fn(|column| {
            self.fitting_uses
                .iter()
                .zip(field.coefficients())
                .map(|(usage, coefficient)| {
                    usage.functional().evaluate_affine(
                        if column == 0 { 1.0 } else { 0.0 },
                        std::array::from_fn(|axis| if axis + 1 == column { 1.0 } else { 0.0 }),
                    ) * coefficient
                })
                .sum::<f64>()
        });
        let recovered_standard_side = self
            .coordinates
            .to_standard_side_condition(physical_side_components);
        let side_condition = PhysicalSideConditionEvidence {
            components: physical_side_components,
            physical_tolerances: self.coordinates.to_physical_side_condition_tolerances(
                [EXECUTED_NUMERICAL_POLICY.side_condition_limit; POLYNOMIAL_DIMENSION],
            ),
            standard_components: standard_side_components,
            recovered_standard_components: recovered_standard_side,
            round_trip_error: relative_slice_error(
                &mapped_physical_side,
                &physical_side_components,
            )
            .max(relative_slice_error(
                &recovered_standard_side,
                &standard_side_components,
            )),
        };
        let recovered_standard_polynomial =
            self.coordinates.to_standard(field.physical_polynomial());
        let polynomial_round_trip_error =
            relative_slice_error(&recovered_standard_polynomial, &standard_polynomial);
        let recovered_standard_coefficients = self
            .coordinates
            .to_standard_field_coefficients(field.coefficients());
        let field_coefficient_round_trip_error =
            relative_slice_error(&recovered_standard_coefficients, &standard_coefficients);
        let field_energy = self.field_energy_normalization.factor() * field.native_cubic_energy();
        let standard_energy = dot_product(
            &standard_coefficients,
            &self.kernel.multiply_vector(&standard_coefficients),
        );
        let recovered_energy = self.field_energy_normalization.factor() * standard_energy
            / self.coordinates.length().powi(3);
        let field_energy_round_trip_error =
            (field_energy - recovered_energy).abs() / recovered_energy.abs().max(1.0);
        Ok(CubicRepresentationRecovery {
            field,
            standard_coefficients,
            recovered_solver_coordinates,
            solver_round_trip_error,
            side_condition,
            polynomial_round_trip_error,
            field_coefficient_round_trip_error,
            field_energy,
            recovered_energy,
            field_energy_round_trip_error,
        })
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
    let mut evidence = RepresentationBuildEvidence::new();
    let (coordinates, _, polynomial) = match assemble_polynomial_pairing(&functionals) {
        Ok(pairing) => pairing,
        Err(failure) => return Some(failure.audited(evidence)),
    };
    evidence.representer_count = Some(functionals.len());
    evidence.failure_stage = RepresentationBuildStage::PolynomialPairing;
    evidence.last_completed_stage = RepresentationBuildStage::RepresenterAssembly;
    evidence.polynomial_dimension = Some(polynomial.columns);
    match verify_polynomial_rank(&polynomial, &functionals, &coordinates) {
        Ok(_) => None,
        Err(failure) => Some(failure.audited(evidence)),
    }
}

pub(crate) fn canonical_fitting_uses(
    equalities: &[CanonicalHardEquality],
    soft_equalities: &[CanonicalSoftEquality],
    affine_inequalities: &[CanonicalAffineInequality],
) -> Vec<FunctionalUse> {
    let mut fitting_uses = Vec::<FunctionalUse>::new();
    for usage in equalities
        .iter()
        .filter(|equality| {
            equality.participation() == CanonicalEqualityParticipation::SolverConstraint
        })
        .filter_map(CanonicalHardEquality::field)
        .chain(soft_equalities.iter().map(CanonicalSoftEquality::field))
        .chain(
            affine_inequalities
                .iter()
                .filter_map(CanonicalAffineInequality::field),
        )
    {
        let normal_basis = (usage.representer_span()
            == FunctionalRepresenterSpan::CompleteGradientAtSupport)
            .then(|| {
                usage
                    .functional()
                    .terms()
                    .first()
                    .map(|term| {
                        (0..3)
                            .map(|axis| {
                                let functional = CanonicalFunctional::new(
                                    FunctionalDimension::FieldValuePerLength,
                                    vec![FunctionalTerm::new(
                                        term.support(),
                                        0.0,
                                        std::array::from_fn(|component| {
                                            if component == axis { 1.0 } else { 0.0 }
                                        }),
                                    )],
                                )
                                .expect("a gradient-component basis functional is valid");
                                FunctionalUse::new(functional, usage.provenance().clone())
                            })
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default()
            })
            .unwrap_or_else(|| vec![usage.clone()]);
        for basis_use in normal_basis {
            if !fitting_uses
                .iter()
                .any(|existing| existing.functional() == basis_use.functional())
            {
                fitting_uses.push(basis_use);
            }
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

#[derive(Debug, Clone)]
pub(crate) struct HouseholderNullSpace {
    basis: Mat<f64>,
    coefficients: Mat<f64>,
    ambient_dimension: usize,
    polynomial_rank: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HouseholderDirection {
    Forward,
    Transpose,
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

    pub(crate) fn reduced_dimension(&self) -> usize {
        self.ambient_dimension - self.polynomial_rank
    }

    fn reflector_count(&self) -> usize {
        self.coefficients.ncols()
    }

    fn apply_on_left(
        &self,
        matrix: MatMut<'_, f64>,
        direction: HouseholderDirection,
    ) -> Result<(), RepresentationFailure> {
        let rhs_columns = matrix.ncols();
        let requirement = apply_block_householder_sequence_on_the_left_in_place_scratch::<f64>(
            self.ambient_dimension,
            self.coefficients.nrows(),
            rhs_columns,
        );
        let mut memory = MemBuffer::try_new(requirement)
            .map_err(|_| RepresentationFailure::NullSpaceWorkspaceAllocation)?;
        match direction {
            HouseholderDirection::Forward => {
                apply_block_householder_sequence_on_the_left_in_place_with_conj(
                    self.basis.as_ref(),
                    self.coefficients.as_ref(),
                    Conj::No,
                    matrix,
                    faer_backend::parallelism(),
                    MemStack::new(&mut memory),
                );
            }
            HouseholderDirection::Transpose => {
                apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj(
                    self.basis.as_ref(),
                    self.coefficients.as_ref(),
                    Conj::No,
                    matrix,
                    faer_backend::parallelism(),
                    MemStack::new(&mut memory),
                );
            }
        }
        Ok(())
    }

    fn apply_on_right(&self, mut matrix: MatMut<'_, f64>) -> Result<(), RepresentationFailure> {
        self.apply_on_left(
            matrix.rb_mut().transpose_mut(),
            HouseholderDirection::Transpose,
        )
    }

    pub(crate) fn expand(&self, reduced: &[f64]) -> Result<Vec<f64>, RepresentationFailure> {
        debug_assert_eq!(reduced.len(), self.reduced_dimension());
        let mut embedded = Mat::<f64>::zeros(self.ambient_dimension, 1);
        for (index, value) in reduced.iter().enumerate() {
            embedded[(self.polynomial_rank + index, 0)] = *value;
        }
        self.apply_on_left(embedded.as_mut(), HouseholderDirection::Forward)?;
        Ok((0..self.ambient_dimension)
            .map(|row| embedded[(row, 0)])
            .collect())
    }

    pub(crate) fn project(&self, ambient: &[f64]) -> Result<Vec<f64>, RepresentationFailure> {
        debug_assert_eq!(ambient.len(), self.ambient_dimension);
        let mut transformed = Mat::<f64>::from_fn(self.ambient_dimension, 1, |row, _| ambient[row]);
        self.apply_on_left(transformed.as_mut(), HouseholderDirection::Transpose)?;
        Ok((self.polynomial_rank..self.ambient_dimension)
            .map(|row| transformed[(row, 0)])
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
            if functionals.len() != polynomial.rows {
                return Err(RepresentationFailure::PolynomialRankGrayZone {
                    evidence: PolynomialRankEvidence {
                        rrqr_ratio: analysis.rrqr_ratio,
                        svd_ratio: analysis.svd_ratio,
                        reject_ratio: analysis.reject_ratio,
                        accept_ratio: analysis.accept_ratio,
                        backend_invoked: false,
                        precision_rescue: None,
                    },
                });
            }
            let rescue = rescue_polynomial_rank(functionals, coordinates)?;
            let rescue_evidence = PrecisionRescueEvidence {
                first_mode: 0,
                mode_count: POLYNOMIAL_DIMENSION,
                precision_bits: DOUBLE_DOUBLE_PRECISION_BITS,
                conclusion: rescue.conclusion,
            };
            match rescue.conclusion {
                PrecisionRescueConclusion::Positive => {}
                PrecisionRescueConclusion::AlgebraicZero => {
                    let standard_mode = smallest_right_mode(
                        matrix.as_ref(),
                        AlgebraicAnalysisStage::PolynomialRank,
                    )?;
                    let mut failure = polynomial_rank_failure(
                        Some(analysis.rank),
                        functionals,
                        coordinates,
                        standard_mode
                            .try_into()
                            .expect("the Cubic polynomial pairing has exactly four columns"),
                    )
                    .expect_err("an algebraic zero is a rank deficiency");
                    if let RepresentationFailure::PolynomialRankDeficient { mode, .. } =
                        &mut failure
                    {
                        mode.precision_rescue = Some(rescue_evidence);
                    }
                    return Err(failure);
                }
                PrecisionRescueConclusion::NegativeCurvature => {
                    return Err(RepresentationFailure::PolynomialNegativeCurvature {
                        evidence: rescue_evidence,
                    });
                }
                PrecisionRescueConclusion::GrayZone
                | PrecisionRescueConclusion::CapacityExceeded => {
                    return Err(RepresentationFailure::PolynomialRankGrayZone {
                        evidence: PolynomialRankEvidence {
                            rrqr_ratio: analysis.rrqr_ratio,
                            svd_ratio: analysis.svd_ratio,
                            reject_ratio: analysis.reject_ratio,
                            accept_ratio: analysis.accept_ratio,
                            backend_invoked: false,
                            precision_rescue: Some(rescue_evidence),
                        },
                    });
                }
            }
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
            precision_rescue: (analysis.decision == SpectralRankDecision::GrayZone).then_some(
                PrecisionRescueEvidence {
                    first_mode: 0,
                    mode_count: POLYNOMIAL_DIMENSION,
                    precision_bits: DOUBLE_DOUBLE_PRECISION_BITS,
                    conclusion: PrecisionRescueConclusion::Positive,
                },
            ),
        },
    ))
}

fn rescue_polynomial_rank(
    functionals: &[CanonicalFunctional],
    coordinates: &CubicSolveCoordinateTransform,
) -> Result<crate::precision_rescue::SymmetricRescueResult, RepresentationFailure> {
    let standard_functionals = functionals
        .iter()
        .map(|functional| coordinates.to_standard_functional(functional))
        .collect::<Result<Vec<_>, _>>()?;
    let mut responses = Vec::with_capacity(standard_functionals.len() * POLYNOMIAL_DIMENSION);
    for functional in &standard_functionals {
        for column in 0..POLYNOMIAL_DIMENSION {
            responses.push(affine_response_dd(functional, column));
        }
    }
    let functional_count = standard_functionals.len();
    let gram = (0..POLYNOMIAL_DIMENSION)
        .flat_map(|row| {
            let responses = &responses;
            (0..POLYNOMIAL_DIMENSION).map(move |column| {
                let mut value = DoubleDouble::from(0.0);
                let mut absolute_scale = 0.0;
                let mut propagated_error = 0.0;
                for index in 0..functional_count {
                    let left = responses[index * POLYNOMIAL_DIMENSION + row];
                    let right = responses[index * POLYNOMIAL_DIMENSION + column];
                    let contribution = left.value * right.value;
                    value += contribution;
                    absolute_scale += contribution.to_f64().abs();
                    propagated_error += left.value.to_f64().abs() * right.error
                        + right.value.to_f64().abs() * left.error
                        + left.error * right.error;
                }
                let mut certified =
                    CertifiedDoubleDouble::new(value, absolute_scale, functional_count * 2);
                certified.error += propagated_error;
                certified
            })
        })
        .collect::<Vec<_>>();
    Ok(classify_symmetric_schur(
        &gram,
        POLYNOMIAL_DIMENSION,
        |mode| {
            standard_functionals.iter().all(|functional| {
                (0..POLYNOMIAL_DIMENSION)
                    .fold(DoubleDouble::from(0.0), |sum, column| {
                        sum + affine_response_dd(functional, column).value * mode[column]
                    })
                    .is_zero()
            })
        },
    ))
}

fn affine_response_dd(functional: &CanonicalFunctional, column: usize) -> CertifiedDoubleDouble {
    let mut value = DoubleDouble::from(0.0);
    let mut absolute_scale = 0.0;
    for term in functional.terms() {
        let contribution = if column == 0 {
            DoubleDouble::from(term.value_coefficient())
        } else {
            DoubleDouble::from(term.value_coefficient())
                * DoubleDouble::from(term.support()[column - 1])
                + DoubleDouble::from(term.gradient_coefficient()[column - 1])
        };
        value += contribution;
        absolute_scale += contribution.to_f64().abs();
    }
    CertifiedDoubleDouble::new(value, absolute_scale, functional.terms().len() * 3)
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
            precision_rescue: None,
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
    if residual > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct QuotientConstructionEvidence {
    pub(crate) quotient_dimension: usize,
    pub(crate) householder_reflector_count: usize,
    pub(crate) congruence_pass_count: usize,
    pub(crate) householder_orthogonality_error: f64,
    pub(crate) null_space_defect: f64,
    pub(crate) canonical_response_round_trip_error: f64,
}

fn implicit_quotient_congruence(
    kernel: &DenseMatrix,
    polynomial: &DenseMatrix,
    null_space: &HouseholderNullSpace,
) -> Result<
    (
        DenseMatrix,
        DenseMatrix,
        DenseMatrix,
        QuotientConstructionEvidence,
    ),
    RepresentationFailure,
> {
    let ambient_dimension = kernel.rows;
    let reduced_dimension = null_space.reduced_dimension();
    let polynomial_rank = null_space.polynomial_rank;

    let mut q1 = Mat::<f64>::from_fn(ambient_dimension, polynomial_rank, |row, column| {
        usize::from(row == column) as f64
    });
    null_space.apply_on_left(q1.as_mut(), HouseholderDirection::Forward)?;
    let householder_orthogonality_error = (0..polynomial_rank)
        .flat_map(|row| {
            let q1 = &q1;
            (0..polynomial_rank).map(move |column| {
                let product = (0..ambient_dimension)
                    .map(|index| q1[(index, row)] * q1[(index, column)])
                    .sum::<f64>();
                (product - f64::from(row == column)).abs()
            })
        })
        .fold(0.0_f64, f64::max);

    let mut transformed_polynomial = polynomial.to_faer();
    null_space.apply_on_left(
        transformed_polynomial.as_mut(),
        HouseholderDirection::Transpose,
    )?;
    let null_space_defect = (polynomial_rank..ambient_dimension)
        .flat_map(|row| {
            let transformed_polynomial = &transformed_polynomial;
            (0..POLYNOMIAL_DIMENSION).map(move |column| transformed_polynomial[(row, column)].abs())
        })
        .fold(0.0_f64, f64::max);
    let trailing_polynomial =
        DenseMatrix::from_fn(reduced_dimension, POLYNOMIAL_DIMENSION, |row, column| {
            transformed_polynomial[(polynomial_rank + row, column)]
        });

    let mut congruence_pass_count = 0;
    let mut transformed_kernel = kernel.to_faer();
    null_space.apply_on_left(transformed_kernel.as_mut(), HouseholderDirection::Transpose)?;
    congruence_pass_count += 1;

    let canonical_response_round_trip_error = if reduced_dimension == 0 {
        0.0
    } else {
        let response_column = (0..ambient_dimension)
            .max_by(|left, right| {
                let norm = |column: usize| {
                    (polynomial_rank..ambient_dimension)
                        .map(|row| transformed_kernel[(row, column)].powi(2))
                        .sum::<f64>()
                };
                norm(*left).total_cmp(&norm(*right))
            })
            .expect("a nonempty representer span has a canonical response");
        let quotient_response = (polynomial_rank..ambient_dimension)
            .map(|row| transformed_kernel[(row, response_column)])
            .collect::<Vec<_>>();
        let expanded = null_space.expand(&quotient_response)?;
        let recovered = null_space.project(&expanded)?;
        relative_slice_error(&recovered, &quotient_response)
    };
    let canonical_responses =
        DenseMatrix::from_fn(reduced_dimension, ambient_dimension, |row, column| {
            transformed_kernel[(polynomial_rank + row, column)]
        });

    null_space.apply_on_right(transformed_kernel.as_mut())?;
    congruence_pass_count += 1;
    let reduced = DenseMatrix::from_fn(reduced_dimension, reduced_dimension, |row, column| {
        transformed_kernel[(polynomial_rank + row, polynomial_rank + column)]
    });
    Ok((
        reduced,
        trailing_polynomial,
        canonical_responses,
        QuotientConstructionEvidence {
            quotient_dimension: reduced_dimension,
            householder_reflector_count: null_space.reflector_count(),
            congruence_pass_count,
            householder_orthogonality_error,
            null_space_defect,
            canonical_response_round_trip_error,
        },
    ))
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OutwardRoundedInterval {
    pub(crate) lower: f64,
    pub(crate) upper: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct QuotientFactorizationEvidence {
    pub(crate) quotient_dimension: usize,
    pub(crate) retained_modes: usize,
    pub(crate) truncated_modes: usize,
    pub(crate) unregularized_llt_count: usize,
    pub(crate) full_spectrum_analysis_count: usize,
    pub(crate) normalized_backward_error: f64,
    pub(crate) pivot_intervals: Vec<OutwardRoundedInterval>,
    pub(crate) field_energy_identity_error: Option<f64>,
    pub(crate) side_condition_error: Option<f64>,
    pub(crate) recovery_round_trip_error: Option<f64>,
    pub(crate) canonical_response_round_trip_error: Option<f64>,
    pub(crate) kernel_ridge_applied: bool,
    pub(crate) gram_jitter_applied: bool,
    pub(crate) mode_truncation_applied: bool,
    pub(crate) precision_rescue: Option<PrecisionRescueEvidence>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PrecisionRescueEvidence {
    pub(crate) first_mode: usize,
    pub(crate) mode_count: usize,
    pub(crate) precision_bits: u32,
    pub(crate) conclusion: PrecisionRescueConclusion,
}

#[derive(Debug, Clone)]
pub(crate) struct EnergyOrthonormalQuotientBasis {
    lower: DenseMatrix,
    permutation: Vec<usize>,
    evidence: QuotientFactorizationEvidence,
}

trait PrecisionRescueSource {
    fn pairing_modes(&self, left: &[f64], right: &[f64]) -> CertifiedDoubleDouble;

    fn algebraic_zero(&self, _mode: &[DoubleDouble]) -> bool {
        false
    }
}

struct DirectPrecisionRescueSource {
    gram: DenseMatrix,
    algebraic_zero: bool,
}

impl DirectPrecisionRescueSource {
    #[cfg(test)]
    fn from_gram(gram: &DenseMatrix) -> Self {
        Self {
            gram: gram.clone(),
            algebraic_zero: false,
        }
    }

    #[cfg(test)]
    fn with_algebraic_zero(gram: &DenseMatrix) -> Self {
        Self {
            gram: gram.clone(),
            algebraic_zero: true,
        }
    }
}

impl PrecisionRescueSource for DirectPrecisionRescueSource {
    fn pairing_modes(&self, left: &[f64], right: &[f64]) -> CertifiedDoubleDouble {
        let mut value = DoubleDouble::from(0.0);
        let mut absolute_scale = 0.0;
        for (row, left) in left.iter().copied().enumerate() {
            for (column, right) in right.iter().copied().enumerate() {
                let contribution = DoubleDouble::from(left)
                    * DoubleDouble::from(self.gram.get(row, column))
                    * DoubleDouble::from(right);
                value += contribution;
                absolute_scale += contribution.to_f64().abs();
            }
        }
        CertifiedDoubleDouble::new(value, absolute_scale, left.len() * right.len() * 3)
    }

    fn algebraic_zero(&self, _mode: &[DoubleDouble]) -> bool {
        self.algebraic_zero
    }
}

struct CanonicalPrecisionRescueSource<'a> {
    functionals: &'a [CanonicalFunctional],
    metric: &'a GlobalAnisotropyMetric,
    null_space: &'a HouseholderNullSpace,
    field_energy_scale: f64,
}

impl PrecisionRescueSource for CanonicalPrecisionRescueSource<'_> {
    fn pairing_modes(&self, left: &[f64], right: &[f64]) -> CertifiedDoubleDouble {
        let row_weights = self
            .null_space
            .expand(left)
            .expect("an analyzed Householder basis remains expandable");
        let column_weights = self
            .null_space
            .expand(right)
            .expect("an analyzed Householder basis remains expandable");
        let mut pairing = DoubleDouble::from(0.0);
        let mut propagated_error = 0.0;
        let mut absolute_scale = 0.0;
        let mut operations = 0;
        for (left, left_weight) in self.functionals.iter().zip(row_weights) {
            if left_weight == 0.0 {
                continue;
            }
            for (right, right_weight) in self.functionals.iter().zip(&column_weights) {
                if *right_weight != 0.0 {
                    let canonical = cubic_pairing_dd_certified(left, right, self.metric);
                    let weight =
                        DoubleDouble::from(left_weight) * DoubleDouble::from(*right_weight);
                    let contribution = weight * canonical.value;
                    pairing += contribution;
                    propagated_error += weight.to_f64().abs() * canonical.error;
                    absolute_scale += contribution.to_f64().abs();
                    operations += 3;
                }
            }
        }
        let mut certified = CertifiedDoubleDouble::new(
            DoubleDouble::from(self.field_energy_scale) * pairing,
            self.field_energy_scale.abs() * absolute_scale,
            operations,
        );
        certified.error += self.field_energy_scale.abs() * propagated_error;
        certified
    }

    fn algebraic_zero(&self, mode: &[DoubleDouble]) -> bool {
        let quotient_mode = mode.iter().map(|value| value.to_f64()).collect::<Vec<_>>();
        let ambient = self
            .null_space
            .expand(&quotient_mode)
            .expect("an analyzed Householder basis remains expandable");
        let mut combined = Vec::<([f64; 3], [DoubleDouble; 4])>::new();
        for (functional, weight) in self.functionals.iter().zip(ambient) {
            for term in functional.terms() {
                let index = combined
                    .iter()
                    .position(|(support, _)| *support == term.support())
                    .unwrap_or_else(|| {
                        combined.push((term.support(), [DoubleDouble::from(0.0); 4]));
                        combined.len() - 1
                    });
                combined[index].1[0] +=
                    DoubleDouble::from(weight) * DoubleDouble::from(term.value_coefficient());
                for axis in 0..3 {
                    combined[index].1[axis + 1] += DoubleDouble::from(weight)
                        * DoubleDouble::from(term.gradient_coefficient()[axis]);
                }
            }
        }
        combined
            .iter()
            .all(|(_, coefficients)| coefficients.iter().all(|value| value.is_zero()))
    }
}

impl EnergyOrthonormalQuotientBasis {
    fn factor(
        gram: &DenseMatrix,
        trailing_polynomial: &DenseMatrix,
        canonical_responses: &DenseMatrix,
        rescue_source: Option<&dyn PrecisionRescueSource>,
    ) -> Result<Self, RepresentationFailure> {
        debug_assert_eq!(gram.rows, gram.columns);
        debug_assert_eq!(trailing_polynomial.rows, gram.rows);
        debug_assert_eq!(canonical_responses.rows, gram.rows);
        let dimension = gram.rows;
        if dimension == 0 {
            return Ok(Self {
                lower: DenseMatrix::from_fn(0, 0, |_, _| unreachable!()),
                permutation: Vec::new(),
                evidence: QuotientFactorizationEvidence {
                    quotient_dimension: 0,
                    retained_modes: 0,
                    truncated_modes: 0,
                    unregularized_llt_count: 0,
                    full_spectrum_analysis_count: 0,
                    normalized_backward_error: 0.0,
                    pivot_intervals: Vec::new(),
                    field_energy_identity_error: Some(0.0),
                    side_condition_error: Some(0.0),
                    recovery_round_trip_error: Some(0.0),
                    canonical_response_round_trip_error: Some(0.0),
                    kernel_ridge_applied: false,
                    gram_jitter_applied: false,
                    mode_truncation_applied: false,
                    precision_rescue: None,
                },
            });
        }

        let factorization = match faer_backend::unregularized_llt(gram.to_faer().as_ref()) {
            Ok(factors) => {
                if factors.dynamic_regularization_count != 0 {
                    return Err(RepresentationFailure::QuotientLltContract {
                        observed: f64::INFINITY,
                        limit: EXECUTED_NUMERICAL_POLICY.quotient_llt_backward_error_limit,
                    });
                }
                let lower = retained_lower_triangle(factors.lower.as_ref());
                let pivot_intervals = (0..dimension)
                    .map(|index| outward_rounded_pivot_interval(gram, &lower, index))
                    .collect::<Vec<_>>();
                if let Some((pivot_index, interval)) = pivot_intervals
                    .iter()
                    .copied()
                    .enumerate()
                    .find(|(_, interval)| interval.lower <= 0.0)
                {
                    if interval.upper <= 0.0 && rescue_source.is_none() {
                        return Err(RepresentationFailure::QuotientFactorizationNotPositive {
                            quotient_dimension: dimension,
                            pivot_index,
                            interval,
                            execution: AnalysisExecutionEvidence::pre_backend(),
                        });
                    }
                    rescue_quotient_factor(gram, &lower, pivot_index, rescue_source)?
                } else {
                    RescuedFactorization {
                        lower,
                        permutation: (0..dimension).collect(),
                        pivot_intervals,
                        precision_rescue: None,
                    }
                }
            }
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
            Err(faer_backend::CholeskyFailure::NonPositivePivot { index }) => {
                let lower = recomputed_lower_for_rescue(gram, index)?;
                let interval = outward_rounded_pivot_interval(gram, &lower, index);
                if interval.upper <= 0.0 && rescue_source.is_none() {
                    return Err(RepresentationFailure::QuotientFactorizationNotPositive {
                        quotient_dimension: dimension,
                        pivot_index: index,
                        interval,
                        execution: AnalysisExecutionEvidence::pre_backend(),
                    });
                }
                rescue_quotient_factor(gram, &lower, index, rescue_source)?
            }
        };
        let RescuedFactorization {
            lower,
            permutation,
            pivot_intervals,
            precision_rescue,
        } = factorization;

        let permuted_gram = permute_symmetric(gram, &permutation);
        let normalized_backward_error = normalized_llt_backward_error(&permuted_gram, &lower);
        let factorization_evidence = QuotientFactorizationEvidence {
            quotient_dimension: dimension,
            retained_modes: dimension,
            truncated_modes: 0,
            unregularized_llt_count: 1,
            full_spectrum_analysis_count: 0,
            normalized_backward_error,
            pivot_intervals,
            field_energy_identity_error: None,
            side_condition_error: None,
            recovery_round_trip_error: None,
            canonical_response_round_trip_error: None,
            kernel_ridge_applied: false,
            gram_jitter_applied: false,
            mode_truncation_applied: false,
            precision_rescue,
        };
        if !normalized_backward_error.is_finite()
            || normalized_backward_error
                > EXECUTED_NUMERICAL_POLICY.quotient_llt_backward_error_limit
        {
            return Err(RepresentationFailure::QuotientLltContract {
                observed: normalized_backward_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_llt_backward_error_limit,
            }
            .with_factorization_evidence(factorization_evidence));
        }

        let provisional = Self {
            lower,
            permutation,
            evidence: factorization_evidence,
        };

        let solver_probe = (0..dimension)
            .map(|index| {
                let magnitude = (index % 7 + 1) as f64;
                if index % 2 == 0 {
                    magnitude
                } else {
                    -magnitude
                }
            })
            .collect::<Vec<_>>();
        let householder_probe = provisional
            .to_householder_coordinates(&solver_probe)
            .map_err(|failure| failure.with_factorization_evidence(provisional.evidence.clone()))?;
        let recovered_probe = provisional
            .to_solver_coordinates(&householder_probe)
            .map_err(|failure| failure.with_factorization_evidence(provisional.evidence.clone()))?;
        let recovery_round_trip_error = relative_slice_error(&recovered_probe, &solver_probe);
        if recovery_round_trip_error
            > EXECUTED_NUMERICAL_POLICY.quotient_basis_recovery_round_trip_limit
        {
            let mut evidence = provisional.evidence.clone();
            evidence.recovery_round_trip_error = Some(recovery_round_trip_error);
            return Err(RepresentationFailure::QuotientRecoveryRoundTripContract {
                observed: recovery_round_trip_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_basis_recovery_round_trip_limit,
            }
            .with_factorization_evidence(evidence));
        }

        let field_energy_identity_error =
            quotient_field_energy_identity_error(&permuted_gram, &provisional.lower);
        if field_energy_identity_error
            > EXECUTED_NUMERICAL_POLICY.quotient_field_energy_identity_limit
        {
            let mut evidence = provisional.evidence.clone();
            evidence.recovery_round_trip_error = Some(recovery_round_trip_error);
            evidence.field_energy_identity_error = Some(field_energy_identity_error);
            return Err(RepresentationFailure::QuotientFieldEnergyIdentityContract {
                observed: field_energy_identity_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_field_energy_identity_limit,
            }
            .with_factorization_evidence(evidence));
        }

        let side_condition_error = (0..trailing_polynomial.columns)
            .map(|column| {
                let response = (0..dimension)
                    .map(|row| trailing_polynomial.get(row, column))
                    .collect::<Vec<_>>();
                provisional
                    .response(&response)
                    .expect("a verified basis produces finite polynomial responses")
                    .into_iter()
                    .map(f64::abs)
                    .fold(0.0_f64, f64::max)
            })
            .fold(0.0_f64, f64::max);
        if side_condition_error > EXECUTED_NUMERICAL_POLICY.quotient_basis_side_condition_limit {
            let mut evidence = provisional.evidence.clone();
            evidence.recovery_round_trip_error = Some(recovery_round_trip_error);
            evidence.field_energy_identity_error = Some(field_energy_identity_error);
            evidence.side_condition_error = Some(side_condition_error);
            return Err(RepresentationFailure::QuotientSideConditionContract {
                observed: side_condition_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_basis_side_condition_limit,
            }
            .with_factorization_evidence(evidence));
        }

        let canonical_response_round_trip_error = (0..canonical_responses.columns)
            .map(|column| {
                let quotient_response = (0..dimension)
                    .map(|row| canonical_responses.get(row, column))
                    .collect::<Vec<_>>();
                let energy_response = provisional.response(&quotient_response)?;
                let recovered_response = provisional
                    .unpermute_from_factor(&provisional.multiply_lower(&energy_response));
                Ok(relative_slice_error(
                    &recovered_response,
                    &quotient_response,
                ))
            })
            .collect::<Result<Vec<_>, RepresentationFailure>>()
            .map_err(|failure| failure.with_factorization_evidence(provisional.evidence.clone()))?
            .into_iter()
            .fold(0.0_f64, f64::max);
        if canonical_response_round_trip_error
            > EXECUTED_NUMERICAL_POLICY.quotient_basis_response_round_trip_limit
        {
            let mut evidence = provisional.evidence.clone();
            evidence.recovery_round_trip_error = Some(recovery_round_trip_error);
            evidence.field_energy_identity_error = Some(field_energy_identity_error);
            evidence.side_condition_error = Some(side_condition_error);
            evidence.canonical_response_round_trip_error =
                Some(canonical_response_round_trip_error);
            return Err(RepresentationFailure::QuotientResponseRoundTripContract {
                observed: canonical_response_round_trip_error,
                limit: EXECUTED_NUMERICAL_POLICY.quotient_basis_response_round_trip_limit,
            }
            .with_factorization_evidence(evidence));
        }

        Ok(Self {
            evidence: QuotientFactorizationEvidence {
                field_energy_identity_error: Some(field_energy_identity_error),
                side_condition_error: Some(side_condition_error),
                recovery_round_trip_error: Some(recovery_round_trip_error),
                canonical_response_round_trip_error: Some(canonical_response_round_trip_error),
                ..provisional.evidence
            },
            ..provisional
        })
    }

    fn response(&self, quotient_response: &[f64]) -> Result<Vec<f64>, RepresentationFailure> {
        let response = self.solve_lower(&self.permute_to_factor(quotient_response));
        finite_basis_coordinates(response)
    }

    fn to_householder_coordinates(
        &self,
        solver_coordinates: &[f64],
    ) -> Result<Vec<f64>, RepresentationFailure> {
        let coordinates = self.solve_lower_transpose(solver_coordinates);
        finite_basis_coordinates(self.unpermute_from_factor(&coordinates))
    }

    fn to_solver_coordinates(
        &self,
        householder_coordinates: &[f64],
    ) -> Result<Vec<f64>, RepresentationFailure> {
        finite_basis_coordinates(
            self.multiply_lower_transpose(&self.permute_to_factor(householder_coordinates)),
        )
    }

    fn solve_lower(&self, right_hand_side: &[f64]) -> Vec<f64> {
        debug_assert_eq!(right_hand_side.len(), self.lower.rows);
        let mut solution = vec![0.0; right_hand_side.len()];
        for row in 0..self.lower.rows {
            let residual = right_hand_side[row]
                - (0..row)
                    .map(|column| self.lower.get(row, column) * solution[column])
                    .sum::<f64>();
            solution[row] = residual / self.lower.get(row, row);
        }
        solution
    }

    fn solve_lower_transpose(&self, right_hand_side: &[f64]) -> Vec<f64> {
        debug_assert_eq!(right_hand_side.len(), self.lower.rows);
        let mut solution = vec![0.0; right_hand_side.len()];
        for row in (0..self.lower.rows).rev() {
            let residual = right_hand_side[row]
                - ((row + 1)..self.lower.rows)
                    .map(|column| self.lower.get(column, row) * solution[column])
                    .sum::<f64>();
            solution[row] = residual / self.lower.get(row, row);
        }
        solution
    }

    fn solve_energy_gram(&self, right_hand_side: &[f64]) -> Vec<f64> {
        let lower_solution = self.solve_lower(&self.permute_to_factor(right_hand_side));
        self.unpermute_from_factor(&self.solve_lower_transpose(&lower_solution))
    }

    fn permute_to_factor(&self, vector: &[f64]) -> Vec<f64> {
        self.permutation
            .iter()
            .map(|index| vector[*index])
            .collect()
    }

    fn unpermute_from_factor(&self, vector: &[f64]) -> Vec<f64> {
        let mut original = vec![0.0; vector.len()];
        for (factor_index, original_index) in self.permutation.iter().copied().enumerate() {
            original[original_index] = vector[factor_index];
        }
        original
    }

    fn multiply_lower(&self, vector: &[f64]) -> Vec<f64> {
        debug_assert_eq!(vector.len(), self.lower.rows);
        (0..self.lower.rows)
            .map(|row| {
                (0..=row)
                    .map(|column| self.lower.get(row, column) * vector[column])
                    .sum()
            })
            .collect()
    }

    fn multiply_lower_transpose(&self, vector: &[f64]) -> Vec<f64> {
        debug_assert_eq!(vector.len(), self.lower.rows);
        (0..self.lower.rows)
            .map(|row| {
                (row..self.lower.rows)
                    .map(|column| self.lower.get(column, row) * vector[column])
                    .sum()
            })
            .collect()
    }
}

struct RescuedFactorization {
    lower: DenseMatrix,
    permutation: Vec<usize>,
    pivot_intervals: Vec<OutwardRoundedInterval>,
    precision_rescue: Option<PrecisionRescueEvidence>,
}

fn rescue_quotient_factor(
    gram: &DenseMatrix,
    prefix_lower: &DenseMatrix,
    first_mode: usize,
    source: Option<&dyn PrecisionRescueSource>,
) -> Result<RescuedFactorization, RepresentationFailure> {
    let dimension = gram.rows;
    let Some(source) = source else {
        return Err(
            RepresentationFailure::QuotientPivotRequiresPrecisionRescue {
                quotient_dimension: dimension,
                pivot_index: first_mode,
                interval: Some(outward_rounded_pivot_interval(
                    gram,
                    prefix_lower,
                    first_mode,
                )),
                execution: AnalysisExecutionEvidence::pre_backend(),
            },
        );
    };
    let (mut lower, mut permutation, first_mode, stable_intervals) =
        isolate_symmetric_stable_prefix(gram);
    let mode_count = dimension - first_mode;
    if mode_count == 0 {
        return Ok(RescuedFactorization {
            lower,
            permutation,
            pivot_intervals: stable_intervals,
            precision_rescue: None,
        });
    }

    let ambiguity_modes = (first_mode..dimension)
        .map(|factor_mode| {
            let mut mode = vec![0.0; dimension];
            mode[permutation[factor_mode]] = 1.0;
            for row in (0..first_mode).rev() {
                let tail = (first_mode..dimension)
                    .map(|column| lower.get(column, row) * mode[permutation[column]])
                    .sum::<f64>();
                mode[permutation[row]] = -tail / lower.get(row, row);
            }
            mode
        })
        .collect::<Vec<_>>();
    let mut schur =
        vec![CertifiedDoubleDouble::new(DoubleDouble::from(0.0), 0.0, 1); mode_count * mode_count];
    for row in 0..mode_count {
        for column in row..mode_count {
            let entry = source.pairing_modes(&ambiguity_modes[row], &ambiguity_modes[column]);
            schur[row * mode_count + column] = entry;
            schur[column * mode_count + row] = entry;
        }
    }
    let rescue = classify_symmetric_schur(&schur, mode_count, |local_mode| {
        let mut full_mode = vec![DoubleDouble::from(0.0); dimension];
        for (mode, coefficient) in ambiguity_modes.iter().zip(local_mode) {
            for (combined, component) in full_mode.iter_mut().zip(mode) {
                *combined += *coefficient * DoubleDouble::from(*component);
            }
        }
        source.algebraic_zero(&full_mode)
    });
    let evidence = PrecisionRescueEvidence {
        first_mode,
        mode_count,
        precision_bits: DOUBLE_DOUBLE_PRECISION_BITS,
        conclusion: rescue.conclusion,
    };
    match rescue.conclusion {
        PrecisionRescueConclusion::Positive => {}
        PrecisionRescueConclusion::AlgebraicZero => {
            return Err(RepresentationFailure::QuotientRankDeficient {
                quotient_dimension: dimension,
                evidence,
                mode: VerifiedCanonicalMode {
                    residual: 0.0,
                    execution: AnalysisExecutionEvidence::pre_backend(),
                    precision_rescue: Some(evidence),
                },
            });
        }
        PrecisionRescueConclusion::NegativeCurvature => {
            return Err(RepresentationFailure::QuotientNegativeCurvature {
                quotient_dimension: dimension,
                evidence,
                execution: AnalysisExecutionEvidence::pre_backend(),
            });
        }
        PrecisionRescueConclusion::GrayZone | PrecisionRescueConclusion::CapacityExceeded => {
            return Err(RepresentationFailure::QuotientPrecisionRescueGrayZone {
                quotient_dimension: dimension,
                evidence,
                execution: AnalysisExecutionEvidence::pre_backend(),
            });
        }
    }

    if (0..dimension).all(|index| prefix_lower.get(index, index) > 0.0) {
        let mut intervals = (0..dimension)
            .map(|index| outward_rounded_pivot_interval(gram, prefix_lower, index))
            .collect::<Vec<_>>();
        let rescued_lower = rescue
            .pivot_lower_bounds
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        for interval in &mut intervals {
            if interval.lower <= 0.0 {
                *interval = OutwardRoundedInterval {
                    lower: rescued_lower,
                    upper: next_up(rescued_lower),
                };
            }
        }
        return Ok(RescuedFactorization {
            lower: prefix_lower.clone(),
            permutation: (0..dimension).collect(),
            pivot_intervals: intervals,
            precision_rescue: Some(evidence),
        });
    }

    let old_permutation = permutation.clone();
    let old_lower = lower.clone();
    for (local_factor, local_original) in rescue.permutation.iter().copied().enumerate() {
        permutation[first_mode + local_factor] = old_permutation[first_mode + local_original];
    }
    for local_factor in 0..mode_count {
        let local_original = rescue.permutation[local_factor];
        for column in 0..first_mode {
            lower.set(
                first_mode + local_factor,
                column,
                old_lower.get(first_mode + local_original, column),
            );
        }
    }
    for row in 0..mode_count {
        for column in 0..=row {
            lower.set(
                first_mode + row,
                first_mode + column,
                rescue.lower[row * mode_count + column].to_f64(),
            );
        }
    }
    let mut intervals = stable_intervals;
    intervals.extend(
        rescue
            .pivot_lower_bounds
            .iter()
            .map(|lower| OutwardRoundedInterval {
                lower: *lower,
                upper: next_up(*lower),
            }),
    );
    Ok(RescuedFactorization {
        lower,
        permutation,
        pivot_intervals: intervals,
        precision_rescue: Some(evidence),
    })
}

fn isolate_symmetric_stable_prefix(
    gram: &DenseMatrix,
) -> (DenseMatrix, Vec<usize>, usize, Vec<OutwardRoundedInterval>) {
    let dimension = gram.rows;
    let mut permutation = (0..dimension).collect::<Vec<_>>();
    let mut lower = DenseMatrix::from_fn(dimension, dimension, |_, _| 0.0);
    let mut intervals = Vec::new();
    for pivot in 0..dimension {
        let selected = (pivot..dimension)
            .max_by(|left, right| {
                let diagonal = |factor_index: usize| {
                    let original = permutation[factor_index];
                    gram.get(original, original)
                        - (0..pivot)
                            .map(|column| lower.get(factor_index, column).powi(2))
                            .sum::<f64>()
                };
                diagonal(*left).total_cmp(&diagonal(*right))
            })
            .expect("the remaining symmetric block is nonempty");
        if selected != pivot {
            permutation.swap(pivot, selected);
            for column in 0..pivot {
                let value = lower.get(pivot, column);
                lower.set(pivot, column, lower.get(selected, column));
                lower.set(selected, column, value);
            }
        }
        let permuted_gram = permute_symmetric(gram, &permutation);
        let interval = outward_rounded_pivot_interval(&permuted_gram, &lower, pivot);
        if interval.lower <= 0.0 {
            return (lower, permutation, pivot, intervals);
        }
        intervals.push(interval);
        let pivot_value = permuted_gram.get(pivot, pivot)
            - (0..pivot)
                .map(|column| lower.get(pivot, column).powi(2))
                .sum::<f64>();
        lower.set(pivot, pivot, pivot_value.sqrt());
        for row in (pivot + 1)..dimension {
            let residual = permuted_gram.get(row, pivot)
                - (0..pivot)
                    .map(|column| lower.get(row, column) * lower.get(pivot, column))
                    .sum::<f64>();
            lower.set(row, pivot, residual / lower.get(pivot, pivot));
        }
    }
    (lower, permutation, dimension, intervals)
}

fn recomputed_lower_for_rescue(
    gram: &DenseMatrix,
    first_mode: usize,
) -> Result<DenseMatrix, RepresentationFailure> {
    let mut lower = DenseMatrix::from_fn(gram.rows, gram.columns, |_, _| 0.0);
    for column in 0..first_mode {
        for row in column..gram.rows {
            let residual = gram.get(row, column)
                - (0..column)
                    .map(|index| lower.get(row, index) * lower.get(column, index))
                    .sum::<f64>();
            let coordinate = if row == column {
                if residual <= 0.0 || !residual.is_finite() {
                    return Err(RepresentationFailure::QuotientLltContract {
                        observed: f64::INFINITY,
                        limit: EXECUTED_NUMERICAL_POLICY.quotient_llt_backward_error_limit,
                    });
                }
                residual.sqrt()
            } else {
                residual / lower.get(column, column)
            };
            if !coordinate.is_finite() {
                return Err(RepresentationFailure::QuotientLltContract {
                    observed: f64::INFINITY,
                    limit: EXECUTED_NUMERICAL_POLICY.quotient_llt_backward_error_limit,
                });
            }
            lower.set(row, column, coordinate);
        }
    }
    Ok(lower)
}

fn permute_symmetric(matrix: &DenseMatrix, permutation: &[usize]) -> DenseMatrix {
    DenseMatrix::from_fn(matrix.rows, matrix.columns, |row, column| {
        matrix.get(permutation[row], permutation[column])
    })
}

fn retained_lower_triangle(matrix: faer::MatRef<'_, f64>) -> DenseMatrix {
    DenseMatrix::from_fn(matrix.nrows(), matrix.ncols(), |row, column| {
        if column <= row {
            matrix[(row, column)]
        } else {
            0.0
        }
    })
}

fn finite_basis_coordinates(coordinates: Vec<f64>) -> Result<Vec<f64>, RepresentationFailure> {
    if coordinates.iter().all(|value| value.is_finite()) {
        Ok(coordinates)
    } else {
        Err(RepresentationFailure::AlgebraicAnalysisFailure {
            stage: AlgebraicAnalysisStage::ReducedCholesky,
            solver_invoked: false,
        })
    }
}

fn outward_rounded_pivot_interval(
    gram: &DenseMatrix,
    lower: &DenseMatrix,
    pivot: usize,
) -> OutwardRoundedInterval {
    let mut sum_lower = 0.0;
    let mut sum_upper = 0.0;
    for column in 0..pivot {
        let value = lower.get(pivot, column);
        let rounded_square = value * value;
        let square_lower = if rounded_square == 0.0 {
            0.0
        } else {
            next_down(rounded_square)
        };
        let square_upper = if value == 0.0 {
            0.0
        } else {
            next_up(rounded_square)
        };
        sum_lower = if sum_lower == 0.0 && square_lower == 0.0 {
            0.0
        } else {
            next_down(sum_lower + square_lower)
        };
        sum_upper = if sum_upper == 0.0 && square_upper == 0.0 {
            0.0
        } else {
            next_up(sum_upper + square_upper)
        };
    }
    if sum_lower == 0.0 && sum_upper == 0.0 {
        let diagonal = gram.get(pivot, pivot);
        OutwardRoundedInterval {
            lower: diagonal,
            upper: diagonal,
        }
    } else {
        OutwardRoundedInterval {
            lower: next_down(gram.get(pivot, pivot) - sum_upper),
            upper: next_up(gram.get(pivot, pivot) - sum_lower),
        }
    }
}

fn normalized_llt_backward_error(gram: &DenseMatrix, lower: &DenseMatrix) -> f64 {
    let lower = lower.to_faer();
    let absolute_lower = Mat::<f64>::from_fn(lower.nrows(), lower.ncols(), |row, column| {
        lower[(row, column)].abs()
    });
    let mut product = Mat::<f64>::zeros(lower.nrows(), lower.ncols());
    let mut absolute_product = Mat::<f64>::zeros(lower.nrows(), lower.ncols());
    matmul(
        product.as_mut(),
        Accum::Replace,
        lower.as_ref(),
        lower.as_ref().transpose(),
        1.0,
        faer_backend::parallelism(),
    );
    matmul(
        absolute_product.as_mut(),
        Accum::Replace,
        absolute_lower.as_ref(),
        absolute_lower.as_ref().transpose(),
        1.0,
        faer_backend::parallelism(),
    );
    let residual_norm = (0..gram.rows)
        .map(|row| {
            (0..gram.columns)
                .map(|column| (gram.get(row, column) - product[(row, column)]).abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let gram_norm = (0..gram.rows)
        .map(|row| {
            (0..gram.columns)
                .map(|column| gram.get(row, column).abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let factor_norm = (0..gram.rows)
        .map(|row| {
            (0..gram.columns)
                .map(|column| absolute_product[(row, column)].abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    residual_norm / (gram_norm + factor_norm)
}

fn quotient_field_energy_identity_error(gram: &DenseMatrix, lower: &DenseMatrix) -> f64 {
    let lower = lower.to_faer();
    let mut transformed = gram.to_faer();
    solve_lower_triangular_in_place(
        lower.as_ref(),
        transformed.as_mut(),
        faer_backend::parallelism(),
    );
    solve_lower_triangular_in_place(
        lower.as_ref(),
        transformed.as_mut().transpose_mut(),
        faer_backend::parallelism(),
    );
    (0..gram.rows)
        .flat_map(|row| {
            let transformed = &transformed;
            (0..gram.columns)
                .map(move |column| (transformed[(row, column)] - f64::from(row == column)).abs())
        })
        .fold(0.0_f64, f64::max)
}

fn quotient_mode_risk_estimates(
    reduced: &DenseMatrix,
    basis: &EnergyOrthonormalQuotientBasis,
    field_energy_scale: f64,
) -> (f64, f64) {
    if reduced.rows == 0 {
        return (f64::INFINITY, f64::INFINITY);
    }
    let seed = (0..reduced.rows)
        .map(|index| {
            let magnitude = (index % 7 + 1) as f64;
            if index % 2 == 0 {
                magnitude
            } else {
                -magnitude
            }
        })
        .collect::<Vec<_>>();
    let largest = fixed_iteration_rayleigh_estimate(reduced, seed.clone(), |vector| {
        reduced.multiply_vector(vector)
    });
    let smallest = fixed_iteration_rayleigh_estimate(reduced, seed, |vector| {
        basis
            .solve_energy_gram(vector)
            .into_iter()
            .map(|value| field_energy_scale * value)
            .collect()
    });
    let pivot_fallback = basis
        .evidence
        .pivot_intervals
        .iter()
        .map(|interval| interval.lower / field_energy_scale)
        .fold(f64::INFINITY, f64::min);
    (
        largest.filter(|value| *value > 0.0).unwrap_or_else(|| {
            (0..reduced.rows)
                .map(|row| {
                    (0..reduced.columns)
                        .map(|column| reduced.get(row, column).abs())
                        .sum::<f64>()
                })
                .fold(0.0_f64, f64::max)
        }),
        smallest
            .filter(|value| *value > 0.0)
            .unwrap_or(pivot_fallback),
    )
}

fn fixed_iteration_rayleigh_estimate(
    matrix: &DenseMatrix,
    mut vector: Vec<f64>,
    mut apply_iteration_operator: impl FnMut(&[f64]) -> Vec<f64>,
) -> Option<f64> {
    const ITERATIONS: usize = 12;
    for _ in 0..ITERATIONS {
        let mut next = apply_iteration_operator(&vector);
        let scale = next.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
        if !scale.is_finite() || scale == 0.0 {
            return None;
        }
        next.iter_mut().for_each(|value| *value /= scale);
        vector = next;
    }
    let product = matrix.multiply_vector(&vector);
    let denominator = dot_product(&vector, &vector);
    let estimate = dot_product(&vector, &product) / denominator;
    estimate.is_finite().then_some(estimate.abs())
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits - 1 } else { bits + 1 })
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
    source_recoveries: Vec<CanonicalHardSourceRecovery>,
    dimension: FunctionalDimension,
    target: f64,
    participation: CanonicalEqualityParticipation,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardSourceRecovery {
    pub(crate) provenance: UsageProvenance,
    pub(crate) coefficient: f64,
    pub(crate) target: f64,
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
            source_recoveries: vec![CanonicalHardSourceRecovery {
                provenance: provenance.clone(),
                coefficient: 1.0,
                target,
            }],
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

    pub(crate) fn source_recoveries(&self) -> &[CanonicalHardSourceRecovery] {
        &self.source_recoveries
    }

    pub(crate) fn add_source_recovery(
        &mut self,
        provenance: UsageProvenance,
        coefficient: f64,
        target: f64,
    ) {
        self.source_recoveries.push(CanonicalHardSourceRecovery {
            provenance,
            coefficient,
            target,
        });
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

    pub(crate) fn constant_shift_response(&self) -> f64 {
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

    pub(crate) fn evaluate(&self, field: &RecoveredCubicField, latents: &[f64]) -> f64 {
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
pub(crate) struct CanonicalAffineInequality {
    field: Option<FunctionalUse>,
    latent_coefficients: Vec<SemanticLatentCoefficient>,
    provenance: UsageProvenance,
    source_provenances: Vec<UsageProvenance>,
    dimension: FunctionalDimension,
    sense: CanonicalInequalitySense,
    bound: f64,
    violation_channel: Option<CanonicalViolationChannel>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum CanonicalInequalitySense {
    Lower,
    Upper,
}

impl CanonicalInequalitySense {
    pub(crate) fn upper_form_multiplier(self) -> f64 {
        match self {
            Self::Lower => -1.0,
            Self::Upper => 1.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CanonicalViolationLoss {
    QuadraticPenalty { weight: f64 },
    LinearViolationPenalty { weight: f64 },
}

impl CanonicalViolationLoss {
    pub(crate) fn weight(self) -> f64 {
        match self {
            Self::QuadraticPenalty { weight } | Self::LinearViolationPenalty { weight } => weight,
        }
    }

    pub(crate) fn is_valid(self) -> bool {
        self.weight().is_finite() && self.weight() > 0.0
    }

    pub(crate) fn objective_contribution(self, violation: f64) -> f64 {
        match self {
            Self::QuadraticPenalty { weight } => 0.5 * weight * violation.powi(2),
            Self::LinearViolationPenalty { weight } => weight * violation,
        }
    }

    pub(crate) fn residual_reference_scale(self) -> f64 {
        match self {
            Self::QuadraticPenalty { weight } => 1.0 / weight.sqrt(),
            Self::LinearViolationPenalty { weight } => 1.0 / weight,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalViolationChannel {
    residual: ResidualId,
    loss: CanonicalViolationLoss,
}

impl CanonicalViolationChannel {
    pub(crate) fn new(residual: ResidualId, loss: CanonicalViolationLoss) -> Self {
        Self { residual, loss }
    }

    pub(crate) fn residual(&self) -> &ResidualId {
        &self.residual
    }

    pub(crate) fn loss(&self) -> CanonicalViolationLoss {
        self.loss
    }
}

impl CanonicalAffineInequality {
    pub(crate) fn upper_bound(
        field: Option<FunctionalUse>,
        latent_coefficients: Vec<SemanticLatentCoefficient>,
        provenance: UsageProvenance,
        dimension: FunctionalDimension,
        upper_bound: f64,
    ) -> Self {
        Self::new(
            field,
            latent_coefficients,
            provenance,
            dimension,
            CanonicalInequalitySense::Upper,
            upper_bound,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        field: Option<FunctionalUse>,
        latent_coefficients: Vec<SemanticLatentCoefficient>,
        provenance: UsageProvenance,
        dimension: FunctionalDimension,
        sense: CanonicalInequalitySense,
        bound: f64,
        violation_channel: Option<CanonicalViolationChannel>,
    ) -> Self {
        Self {
            field,
            latent_coefficients,
            source_provenances: vec![provenance.clone()],
            provenance,
            dimension,
            sense,
            bound,
            violation_channel,
        }
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

    pub(crate) fn source_provenances(&self) -> &[UsageProvenance] {
        &self.source_provenances
    }

    pub(crate) fn add_source_provenance(&mut self, provenance: UsageProvenance) {
        debug_assert!(!self.source_provenances.contains(&provenance));
        self.source_provenances.push(provenance);
        self.source_provenances
            .sort_by(|left, right| left.source().cmp(right.source()));
    }

    pub(crate) fn dimension(&self) -> FunctionalDimension {
        self.dimension
    }

    pub(crate) fn bound(&self) -> f64 {
        self.bound
    }

    pub(crate) fn sense(&self) -> CanonicalInequalitySense {
        self.sense
    }

    pub(crate) fn violation_channel(&self) -> Option<&CanonicalViolationChannel> {
        self.violation_channel.as_ref()
    }

    pub(crate) fn upper_form_bound(&self) -> f64 {
        self.sense.upper_form_multiplier() * self.bound
    }

    pub(crate) fn constant_shift_response(&self) -> f64 {
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

    pub(crate) fn evaluate(&self, field: &RecoveredCubicField, latents: &[f64]) -> f64 {
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

    pub(crate) fn physical_margin(&self, value: f64) -> f64 {
        match self.sense {
            CanonicalInequalitySense::Lower => value - self.bound,
            CanonicalInequalitySense::Upper => self.bound - value,
        }
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

    pub(crate) fn whitening_matrix(&self, dimension: usize) -> Vec<f64> {
        match self {
            Self::QuadraticPenalty { weight } => diagonal_matrix(dimension, weight.sqrt()),
            Self::StandardDeviation { standard_deviation } => {
                diagonal_matrix(dimension, 1.0 / standard_deviation)
            }
            Self::Covariance { whitening, .. } => whitening.clone(),
        }
    }

    pub(crate) fn inverse_whitening_matrix(&self, dimension: usize) -> Vec<f64> {
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

    pub(crate) fn is_valid(&self, dimension: usize) -> bool {
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

    pub(crate) fn residual_reference_scale(&self) -> f64 {
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

    pub(crate) fn evaluate(&self, field: &RecoveredCubicField) -> f64 {
        field.evaluate_functional(self.field.functional())
    }

    pub(crate) fn constant_shift_response(&self) -> f64 {
        self.field
            .functional()
            .terms()
            .iter()
            .map(|term| term.value_coefficient())
            .sum()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CanonicalSoftResidualMemberKind {
    FieldValue,
    Gradient,
    Tangent,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CanonicalHardResidualBlockKind {
    NormalProjection,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CanonicalHardResidualBlock {
    canonical_indices: Vec<usize>,
    kind: CanonicalHardResidualBlockKind,
}

impl CanonicalHardResidualBlock {
    pub(crate) fn normal_projection(canonical_indices: Vec<usize>) -> Self {
        Self {
            canonical_indices,
            kind: CanonicalHardResidualBlockKind::NormalProjection,
        }
    }

    pub(crate) fn canonical_indices(&self) -> &[usize] {
        &self.canonical_indices
    }

    pub(crate) fn is_valid(&self, equality_count: usize) -> bool {
        match self.kind {
            CanonicalHardResidualBlockKind::NormalProjection => {
                matches!(self.canonical_indices.len(), 2 | 3)
                    && self
                        .canonical_indices
                        .iter()
                        .all(|index| *index < equality_count)
            }
        }
    }
}

impl CanonicalSoftResidualMemberKind {
    pub(crate) fn component_count(self) -> usize {
        match self {
            Self::FieldValue | Self::Tangent => 1,
            Self::Gradient => 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CanonicalSoftResidualBlockKind {
    Independent(CanonicalSoftResidualMemberKind),
    NormalProjection,
    CovarianceGroup {
        members: Vec<CanonicalSoftResidualMemberKind>,
    },
}

impl CanonicalSoftResidualBlockKind {
    pub(crate) fn is_valid(
        &self,
        residual_count: usize,
        covariance_group: Option<&GroupId>,
    ) -> bool {
        match self {
            Self::Independent(member) => {
                covariance_group.is_none() && member.component_count() == residual_count
            }
            Self::NormalProjection => covariance_group.is_none() && matches!(residual_count, 2 | 3),
            Self::CovarianceGroup { members } => {
                covariance_group.is_some()
                    && !members.is_empty()
                    && members.iter().try_fold(0_usize, |count, member| {
                        count.checked_add(member.component_count())
                    }) == Some(residual_count)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftObjective {
    residuals: Vec<ResidualId>,
    loss: CanonicalSoftLoss,
    covariance_group: Option<GroupId>,
    block_kind: CanonicalSoftResidualBlockKind,
}

impl CanonicalSoftObjective {
    pub(crate) fn new(residual: ResidualId, loss: CanonicalSoftLoss) -> Self {
        Self::new_block(
            vec![residual],
            loss,
            None,
            CanonicalSoftResidualBlockKind::Independent(
                CanonicalSoftResidualMemberKind::FieldValue,
            ),
        )
    }

    pub(crate) fn new_block(
        residuals: Vec<ResidualId>,
        loss: CanonicalSoftLoss,
        covariance_group: Option<GroupId>,
        block_kind: CanonicalSoftResidualBlockKind,
    ) -> Self {
        Self {
            residuals,
            loss,
            covariance_group,
            block_kind,
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

    pub(crate) fn block_kind(&self) -> &CanonicalSoftResidualBlockKind {
        &self.block_kind
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CubicCanonicalProblem {
    pub(crate) equalities: Vec<CanonicalHardEquality>,
    pub(crate) hard_residual_blocks: Vec<CanonicalHardResidualBlock>,
    pub(crate) affine_inequalities: Vec<CanonicalAffineInequality>,
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
    standard_rows: Vec<Vec<f64>>,
    provenances: Vec<UsageProvenance>,
    residuals: Vec<ResidualId>,
    targets: Vec<f64>,
    canonical_precision: Vec<f64>,
    standard_precision: Vec<f64>,
    whitening: Vec<f64>,
    inverse_whitening: Vec<f64>,
    covariance_group: Option<GroupId>,
    block_kind: CanonicalSoftResidualBlockKind,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FieldSample {
    pub(crate) value: f64,
    pub(crate) gradient: [f64; 3],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum QuerySampleFailure {
    NonFiniteResult,
    NumericalIndeterminate,
}

#[derive(Debug, Clone, Copy)]
struct DoubleDoubleQuerySample {
    components: [DoubleDouble; 4],
    error_bounds: [f64; 4],
}

#[derive(Debug, Clone, Copy, Default)]
struct CompensatedQuerySum {
    sum: f64,
    correction: f64,
    absolute_scale: f64,
    operations: usize,
    has_positive: bool,
    has_negative: bool,
}

impl CompensatedQuerySum {
    fn add(&mut self, value: f64, operations: usize) {
        self.add_with_scale(value, value.abs(), operations);
    }

    fn add_with_scale(&mut self, value: f64, absolute_scale: f64, operations: usize) {
        self.operations = self.operations.saturating_add(operations.max(1));
        self.has_positive |= value.is_sign_positive() && value != 0.0;
        self.has_negative |= value.is_sign_negative() && value != 0.0;
        self.absolute_scale += absolute_scale;
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }

    fn error_bound(self) -> f64 {
        let relative_bound = 8.0 * self.operations as f64 * f64::EPSILON;
        if relative_bound >= 1.0 {
            f64::INFINITY
        } else {
            relative_bound / (1.0 - relative_bound) * self.absolute_scale
        }
    }

    fn has_cancellation(self) -> bool {
        self.has_positive && self.has_negative
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RecoveredCubicField {
    representers: Vec<CanonicalFunctional>,
    metric: GlobalAnisotropyMetric,
    coefficients: Vec<f64>,
    physical_polynomial: [f64; 4],
    query_coordinates: CubicSolveCoordinateTransform,
    query_standard_polynomial: [f64; 4],
    query_component_field_scales: Option<[f64; 4]>,
}

impl RecoveredCubicField {
    pub(crate) fn from_standard_candidate(
        representation: &CubicRepresentation,
        standard_coefficients: &[f64],
        standard_polynomial: [f64; POLYNOMIAL_DIMENSION],
    ) -> Self {
        Self {
            representers: representation
                .fitting_uses
                .iter()
                .map(|usage| usage.functional().clone())
                .collect(),
            metric: representation.metric.clone(),
            coefficients: representation
                .coordinates
                .to_physical_field_coefficients(standard_coefficients),
            physical_polynomial: representation.coordinates.to_physical(standard_polynomial),
            query_coordinates: representation.coordinates,
            query_standard_polynomial: standard_polynomial,
            query_component_field_scales: None,
        }
    }

    pub(crate) fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    pub(crate) fn physical_polynomial(&self) -> [f64; 4] {
        self.physical_polynomial
    }

    pub(crate) fn finalize_verified_query_representation(
        &mut self,
        field_scale: f64,
        characteristic_length: f64,
        basis_round_trip_error: f64,
    ) -> Result<f64, QuerySampleFailure> {
        debug_assert!(field_scale.is_finite() && field_scale >= 0.0);
        debug_assert!(characteristic_length.is_finite() && characteristic_length > 0.0);
        debug_assert!(
            basis_round_trip_error.is_finite()
                && basis_round_trip_error <= EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit
        );
        self.query_component_field_scales = Some([
            field_scale,
            field_scale / characteristic_length,
            field_scale / characteristic_length,
            field_scale / characteristic_length,
        ]);
        let mut round_trip_error = 0.0_f64;
        for functional in &self.representers {
            let recovered_response = self.evaluate_functional(functional);
            let mut query_response = CompensatedQuerySum::default();
            for term in functional.terms() {
                let sample = self.reliable_sample(term.support())?;
                query_response.add(
                    term.value_coefficient() * sample.value
                        + dot3(term.gradient_coefficient(), sample.gradient),
                    16,
                );
            }
            round_trip_error = round_trip_error.max(
                (query_response.value() - recovered_response).abs()
                    / recovered_response.abs().max(1.0),
            );
        }
        Ok(round_trip_error)
    }

    pub(crate) fn reliable_sample(
        &self,
        point: [f64; 3],
    ) -> Result<FieldSample, QuerySampleFailure> {
        let component_field_scales = self
            .query_component_field_scales
            .expect("only a verified recovered field can reach SolvedModel");
        let primary = self.compensated_sample(point);
        let primary_values = [
            primary[0].value(),
            primary[1].value(),
            primary[2].value(),
            primary[3].value(),
        ];
        if primary_values
            .iter()
            .zip(primary)
            .zip(component_field_scales)
            .all(|((value, sum), field_scale)| {
                value.is_finite()
                    && sum.error_bound() <= query_reliability_envelope(field_scale, *value)
            })
        {
            return Ok(FieldSample {
                value: primary_values[0],
                gradient: [primary_values[1], primary_values[2], primary_values[3]],
            });
        }

        let upgraded = self.double_double_sample(point);
        let upgraded_values = upgraded.components.map(DoubleDouble::to_f64);
        if upgraded_values
            .iter()
            .zip(upgraded.error_bounds)
            .zip(component_field_scales)
            .all(|((value, error), field_scale)| {
                value.is_finite() && error <= query_reliability_envelope(field_scale, *value)
            })
        {
            return Ok(FieldSample {
                value: upgraded_values[0],
                gradient: [upgraded_values[1], upgraded_values[2], upgraded_values[3]],
            });
        }

        if primary
            .into_iter()
            .zip(upgraded_values)
            .any(|(sum, value)| !value.is_finite() && !sum.has_cancellation())
        {
            Err(QuerySampleFailure::NonFiniteResult)
        } else {
            Err(QuerySampleFailure::NumericalIndeterminate)
        }
    }

    fn compensated_sample(&self, point: [f64; 3]) -> [CompensatedQuerySum; 4] {
        // Keep the f64 and double-double expansions arithmetically independent:
        // sharing a numeric abstraction here would make the rescue path repeat
        // primary-rounding defects instead of recomputing the same field.
        let mut sums = [CompensatedQuerySum::default(); 4];
        sums[0].add(self.query_standard_polynomial[0], 1);
        for axis in 0..3 {
            let standard_coordinate = (point[axis] - self.query_coordinates.center()[axis])
                / self.query_coordinates.length();
            sums[0].add(
                self.query_standard_polynomial[axis + 1] * standard_coordinate,
                4,
            );
            sums[axis + 1].add(
                self.query_standard_polynomial[axis + 1] / self.query_coordinates.length(),
                2,
            );
        }
        for (coefficient, functional) in self.coefficients.iter().zip(&self.representers) {
            for term in functional.terms() {
                let jet = CubicKernel::jet(point, term.support(), &self.metric);
                let derivative = term.gradient_coefficient();
                let value_term = term.value_coefficient() * jet.value();
                let derivative_terms: [f64; 3] =
                    std::array::from_fn(|axis| derivative[axis] * jet.gradient_y()[axis]);
                sums[0].add_with_scale(
                    coefficient * (value_term + derivative_terms.into_iter().sum::<f64>()),
                    coefficient.abs()
                        * (value_term.abs()
                            + derivative_terms.into_iter().map(f64::abs).sum::<f64>()),
                    128,
                );
                let mixed_derivative = jet.mixed_xy();
                for axis in 0..3 {
                    let value_derivative_term = term.value_coefficient() * jet.gradient_x()[axis];
                    let mixed_terms: [f64; 3] = std::array::from_fn(|column| {
                        mixed_derivative[axis][column] * derivative[column]
                    });
                    sums[axis + 1].add_with_scale(
                        coefficient
                            * (value_derivative_term + mixed_terms.into_iter().sum::<f64>()),
                        coefficient.abs()
                            * (value_derivative_term.abs()
                                + mixed_terms.into_iter().map(f64::abs).sum::<f64>()),
                        128,
                    );
                }
            }
        }
        sums
    }

    fn double_double_sample(&self, point: [f64; 3]) -> DoubleDoubleQuerySample {
        let mut sums = [DoubleDouble::from(0.0); 4];
        let mut absolute_scales = [0.0_f64; 4];
        let mut operations = [0_usize; 4];
        let mut add = |component: usize, contribution: DoubleDouble, operation_count: usize| {
            sums[component] += contribution;
            absolute_scales[component] += contribution.to_f64().abs();
            operations[component] += operation_count;
        };
        add(0, DoubleDouble::from(self.query_standard_polynomial[0]), 1);
        for (axis, coordinate) in point.into_iter().enumerate() {
            let standard_coordinate = (DoubleDouble::from(coordinate)
                - DoubleDouble::from(self.query_coordinates.center()[axis]))
                / DoubleDouble::from(self.query_coordinates.length());
            add(
                0,
                DoubleDouble::from(self.query_standard_polynomial[axis + 1]) * standard_coordinate,
                4,
            );
            add(
                axis + 1,
                DoubleDouble::from(self.query_standard_polynomial[axis + 1])
                    / DoubleDouble::from(self.query_coordinates.length()),
                2,
            );
        }
        for (coefficient, representer) in self.coefficients.iter().zip(&self.representers) {
            let coefficient = DoubleDouble::from(*coefficient);
            for term in representer.terms() {
                let jet = cubic_jet_dd(point, term.support(), self.metric.matrix());
                let derivative = term.gradient_coefficient().map(DoubleDouble::from);
                add(
                    0,
                    coefficient
                        * (DoubleDouble::from(term.value_coefficient()) * jet.value()
                            + dot3_double_double(derivative, jet.gradient_y())),
                    128,
                );
                let mixed = jet.mixed_xy();
                for (axis, mixed_row) in mixed.into_iter().enumerate() {
                    add(
                        axis + 1,
                        coefficient
                            * (DoubleDouble::from(term.value_coefficient())
                                * jet.gradient_x()[axis]
                                + dot3_double_double(mixed_row, derivative)),
                        128,
                    );
                }
            }
        }
        let errors = std::array::from_fn(|component| {
            8.0 * operations[component] as f64 * f64::EPSILON.powi(2) * absolute_scales[component]
        });
        DoubleDoubleQuerySample {
            components: sums,
            error_bounds: errors,
        }
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

    pub(crate) fn native_cubic_energy(&self) -> f64 {
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

fn query_reliability_envelope(field_scale: f64, sample_reference_scale: f64) -> f64 {
    EXECUTED_NUMERICAL_POLICY.query_field_scale_tolerance_multiplier * field_scale
        + EXECUTED_NUMERICAL_POLICY.query_sample_reference_tolerance_multiplier
            * sample_reference_scale.abs()
}

fn dot3_double_double(left: [DoubleDouble; 3], right: [DoubleDouble; 3]) -> DoubleDouble {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
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
    pub(crate) block_kind: CanonicalSoftResidualBlockKind,
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
    pub(crate) fn from_dimensioned_residuals(
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
    pub(crate) fn is_within_policy(self) -> bool {
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
    pub(crate) hard_recovery: CanonicalHardRecoveryGraph,
    pub(crate) all_source_recovery: AllSourceRecoveryLedger,
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
    pub(crate) query_response_round_trip_error: f64,
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
    /// One or more participating SourceIds lacked a complete recovery edge.
    SourceCoverageMismatch,
    /// At least one recovered physical quantity was non-finite.
    NonFiniteRecoveredQuantity,
    /// The recovered Cubic Π₁ side condition exceeded its tolerance.
    SideConditionViolation,
    /// The side-condition recovery map exceeded its round-trip limit.
    SideConditionRoundTripViolation,
    /// At least one recovered hard equality exceeded physical tolerance.
    HardEqualityViolation,
    /// At least one recovered affine inequality exceeded physical tolerance.
    AffineInequalityViolation,
    /// A recovered inequality slack disagreed with the backend equation.
    BackendSlackMismatch,
    /// The reduced/full field map exceeded its round-trip limit.
    ReductionRoundTripViolation,
    /// The QP scaling map exceeded its round-trip limit.
    ScalingRoundTripViolation,
    /// Polynomial recovery exceeded its round-trip limit.
    PolynomialRoundTripViolation,
    /// Field-coefficient recovery exceeded its round-trip limit.
    FieldCoefficientRoundTripViolation,
    /// Recovered-basis responses disagreed with the finalized query path.
    QueryRepresentationRoundTripViolation,
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
    AffineInequalityRequiresConvexQp,
    NonFiniteTarget {
        equality: usize,
    },
    Representation(Box<RepresentationFailure>),
    DirectInputConflict {
        evidence: CanonicalHardConflictWitness,
        representation: Box<CpdEvidence>,
    },
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
    static INJECTED_RESPONSE_ASSEMBLY_FAILURE: RefCell<bool> = const { RefCell::new(false) };
}

#[cfg(test)]
pub(crate) fn inject_response_assembly_failure_once() {
    INJECTED_RESPONSE_ASSEMBLY_FAILURE.with(|slot| {
        assert!(
            !slot.replace(true),
            "only one response failure may be injected per fit"
        );
    });
}

#[cfg(test)]
fn take_injected_response_assembly_failure() -> bool {
    INJECTED_RESPONSE_ASSEMBLY_FAILURE.with(|slot| slot.replace(false))
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
                hard_residual_blocks: Vec::new(),
                affine_inequalities: Vec::new(),
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
        if !problem.affine_inequalities.is_empty() {
            return Err(CubicEqualityFailure::AffineInequalityRequiresConvexQp);
        }
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
        if problem
            .hard_residual_blocks
            .iter()
            .any(|block| !block.is_valid(problem.equalities.len()))
        {
            return Err(CubicEqualityFailure::NonFiniteTarget {
                equality: problem.equalities.len(),
            });
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
            || problem.soft_objectives.iter().any(|objective| {
                !objective.loss.is_valid(objective.residuals.len())
                    || !objective.block_kind.is_valid(
                        objective.residuals.len(),
                        objective.covariance_group.as_ref(),
                    )
            })
        {
            return Err(CubicEqualityFailure::NonFiniteTarget {
                equality: problem.equalities.len(),
            });
        }
        let fitting_uses = canonical_fitting_uses(
            &problem.equalities,
            &problem.soft_equalities,
            &problem.affine_inequalities,
        );
        let primal_variables = fitting_uses
            .len()
            .checked_add(POLYNOMIAL_DIMENSION)
            .and_then(|count| count.checked_add(problem.semantic_latents.len()))
            .unwrap_or(usize::MAX);
        let (representation, field_form) =
            CubicRepresentation::build(fitting_uses, metric, problem.field_energy_normalization)
                .map_err(|failure| CubicEqualityFailure::Representation(Box::new(failure)))?;
        let solver_form = CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
            .map_err(|failure| {
                CubicEqualityFailure::Representation(Box::new(
                    representation.audit_response_assembly_failure(failure),
                ))
            })?;
        if solver_form.verifies_hard_conflict_witness() {
            if let Some(evidence) = solver_form.hard_recovery.conflict_witness.clone() {
                return Err(CubicEqualityFailure::DirectInputConflict {
                    evidence,
                    representation: Box::new(solver_form.representation_evidence.clone()),
                });
            }
        }
        let equality_constraints = solver_form
            .solver_hard_rows()
            .count()
            .checked_add(POLYNOMIAL_DIMENSION)
            .unwrap_or(usize::MAX);
        plan_equality_capacity_after_representation(
            primal_variables,
            equality_constraints,
            &solver_form.representation_evidence,
        )?;
        #[cfg(test)]
        if let Some(failure) = take_injected_kkt_failure() {
            return Err(CubicEqualityFailure::Backend {
                failure: Box::new(failure),
                representation: Box::new(solver_form.representation_evidence.clone()),
            });
        }
        let (assembly, backend) = solve_standard_form(&solver_form)?;
        recover_and_verify(representation, solver_form, problem, assembly, backend)
    }
}

fn plan_equality_capacity_after_representation(
    primal_variables: usize,
    equality_constraints: usize,
    representation: &CpdEvidence,
) -> Result<(), CubicEqualityFailure> {
    plan_equality_capacity(primal_variables, equality_constraints)
        .map(|_| ())
        .map_err(|failure| {
            let evidence = RepresentationBuildEvidence::completed(representation);
            CubicEqualityFailure::Representation(Box::new(
                RepresentationFailure::Capacity(Box::new(failure)).audited(evidence),
            ))
        })
}

fn solve_standard_form(
    form: &CanonicalCubicSolverForm,
) -> Result<(EqualityAssemblyEvidence, KktSolveEvidence), CubicEqualityFailure> {
    let coordinate_layout = CubicFieldCoordinateLayout::Standard;
    let variable_layout = form.variable_layout(coordinate_layout, 0);
    let solver_equalities = form.solver_hard_rows().collect::<Vec<_>>();
    let coefficient_count = variable_layout.field;
    let primal_variables = variable_layout
        .variables()
        .expect("the checked Equality variable layout is finite");
    let equality_constraints = POLYNOMIAL_DIMENSION + solver_equalities.len();
    let solver_equality_rows = solver_equalities
        .iter()
        .map(|equality| {
            equality
                .row
                .coefficients(coordinate_layout, variable_layout)
        })
        .collect::<Vec<_>>();
    let soft_objective_blocks = form
        .soft_objectives
        .iter()
        .map(|objective| {
            let relations = objective
                .canonical_indices
                .iter()
                .map(|index| &form.soft_rows[*index].row)
                .collect::<Vec<_>>();
            AssembledSoftObjectiveBlock {
                objective_index: objective.objective_index,
                canonical_indices: objective.canonical_indices.clone(),
                standard_rows: relations
                    .iter()
                    .map(|relation| relation.coefficients(coordinate_layout, variable_layout))
                    .collect(),
                provenances: relations
                    .iter()
                    .map(|relation| relation.provenance.clone())
                    .collect(),
                residuals: objective.residuals.clone(),
                targets: relations.iter().map(|relation| relation.target).collect(),
                standard_precision: objective.precision.clone(),
                canonical_precision: objective.precision.clone(),
                whitening: objective.whitening.clone(),
                inverse_whitening: objective.inverse_whitening.clone(),
                covariance_group: objective.covariance_group.clone(),
                block_kind: objective.block_kind.clone(),
            }
        })
        .collect::<Vec<_>>();
    let hessian = DenseMatrix::from_fn(primal_variables, primal_variables, |row, column| {
        let field_energy = if row < coefficient_count && column < coefficient_count {
            form.field_energy(coordinate_layout)[row * coefficient_count + column]
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
                                    * objective.standard_rows[left][row]
                                    * objective.standard_rows[right][column]
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
                    form.standard_side_conditions[row * coefficient_count + column]
                } else {
                    0.0
                }
            } else {
                solver_equality_rows[row - POLYNOMIAL_DIMENSION][column]
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
                                objective.standard_rows[left][column]
                                    * objective.standard_precision[left * dimension + right]
                                    * objective.targets[right]
                            })
                        })
                        .sum::<f64>()
                })
                .sum()
        })
        .collect::<Vec<_>>();
    let equality_rhs = std::iter::repeat_n(0.0, POLYNOMIAL_DIMENSION)
        .chain(solver_equalities.iter().map(|equality| equality.row.target))
        .collect::<Vec<_>>();
    let hard_equality_rows = solver_equalities
        .iter()
        .enumerate()
        .map(|(solver_index, equality)| {
            let kkt_equality_row = POLYNOMIAL_DIMENSION + solver_index;
            AssembledHardEqualityRow {
                kkt_equality_row,
                canonical_index: equality.row.canonical_index,
                solver_index,
                provenance: equality.row.provenance.clone(),
                derived_block: equality.row.derived_block.clone(),
                residual: equality.row.residual.clone(),
                derived_row: equality.row.derived_row.clone(),
                derived_column: equality.row.derived_column.clone(),
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
        representation: Box::new(form.representation_evidence.clone()),
    })?;
    Ok((
        EqualityAssemblyEvidence {
            primal_variables,
            field_coefficients: coefficient_count,
            polynomial_coefficients: POLYNOMIAL_DIMENSION,
            semantic_latents: form.semantic_latents,
            side_conditions: POLYNOMIAL_DIMENSION,
            hard_equalities: solver_equalities.len(),
            canonical_hard_equalities: form.hard_rows.len(),
            hard_equality_rows,
            soft_objective_blocks,
        },
        backend,
    ))
}

fn verifies_assembled_provenance_and_rows(
    solver_form: &CanonicalCubicSolverForm,
    assembly: &EqualityAssemblyEvidence,
    backend: &KktSolveEvidence,
) -> bool {
    let variable_layout = solver_form.variable_layout(CubicFieldCoordinateLayout::Standard, 0);
    let solver_equalities = solver_form.solver_hard_rows().collect::<Vec<_>>();
    if !solver_form.verifies_hard_recovery()
        || !solver_form.verifies_soft_recovery()
        || assembly.hard_equality_rows.len() != solver_equalities.len()
        || assembly.soft_objective_blocks.len() != solver_form.soft_objectives.len()
        || backend.equality_multipliers.len() != assembly.side_conditions + assembly.hard_equalities
    {
        return false;
    }
    let hard_rows_verified = assembly
        .hard_equality_rows
        .iter()
        .zip(solver_equalities)
        .enumerate()
        .all(|(solver_index, (row, equality))| {
            let expected = equality
                .row
                .coefficients(CubicFieldCoordinateLayout::Standard, variable_layout);
            row.canonical_index == equality.row.canonical_index
                && row.solver_index == solver_index
                && row.kkt_equality_row == POLYNOMIAL_DIMENSION + solver_index
                && row.provenance == equality.row.provenance
                && row.derived_block == equality.row.derived_block
                && row.residual == equality.row.residual
                && row.derived_row == equality.row.derived_row
                && row.derived_column == equality.row.derived_column
                && row.rhs == equality.row.target
                && row.standard_jacobian_row == expected
        });
    let soft_rows_verified = assembly
        .soft_objective_blocks
        .iter()
        .zip(&solver_form.soft_objectives)
        .enumerate()
        .all(|(objective_index, (block, objective))| {
            let dimension = objective.residuals.len();
            block.objective_index == objective_index
                && objective.objective_index == objective_index
                && block.canonical_indices.len() == dimension
                && block.standard_rows.len() == dimension
                && block.provenances.len() == dimension
                && block.residuals == objective.residuals
                && block.targets.len() == dimension
                && block.covariance_group == objective.covariance_group
                && block.block_kind == objective.block_kind
                && objective.precision == objective.loss.precision_matrix(dimension)
                && objective.whitening == objective.loss.whitening_matrix(dimension)
                && objective.inverse_whitening == objective.loss.inverse_whitening_matrix(dimension)
                && block.canonical_precision == objective.precision
                && block.standard_precision == block.canonical_precision
                && block.whitening == objective.whitening
                && block.inverse_whitening == objective.inverse_whitening
                && block
                    .canonical_indices
                    .iter()
                    .enumerate()
                    .all(|(component, canonical_index)| {
                        solver_form
                            .soft_rows
                            .get(*canonical_index)
                            .is_some_and(|relation| {
                                block.standard_rows[component]
                                    == relation.row.coefficients(
                                        CubicFieldCoordinateLayout::Standard,
                                        variable_layout,
                                    )
                                    && block.provenances[component] == relation.row.provenance
                                    && block.residuals[component] == relation.row.residual
                                    && block.targets[component] == relation.row.target
                            })
                    })
        });
    hard_rows_verified && soft_rows_verified
}

fn recover_and_verify(
    representation: CubicRepresentation,
    solver_form: CanonicalCubicSolverForm,
    problem: CubicCanonicalProblem,
    assembly: EqualityAssemblyEvidence,
    backend: KktSolveEvidence,
) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
    let representation_evidence = solver_form.representation_evidence.clone();
    let characteristic_length = solver_form.characteristic_length;
    let provenance_verified =
        verifies_assembled_provenance_and_rows(&solver_form, &assembly, &backend);
    if !provenance_verified {
        return Err(CubicEqualityFailure::RecoveryVerification {
            evidence: Box::new(RecoveryVerificationFailureEvidence::early(
                RecoveryVerificationFailureReason::ProvenanceMismatch,
            )),
            representation: Box::new(representation_evidence.clone()),
            backend: Box::new(backend),
        });
    }

    let coefficient_count = assembly.field_coefficients;
    let standard_polynomial: [f64; 4] = backend.candidate
        [coefficient_count..coefficient_count + POLYNOMIAL_DIMENSION]
        .try_into()
        .expect("the augmented primal retains exactly four polynomial coefficients");
    let recovered_representation = match representation.recover(
        CubicSolverFieldCoordinates::Standard(&backend.candidate[..coefficient_count]),
        standard_polynomial,
    ) {
        Ok(recovered) => recovered,
        Err(CubicRepresentationRecoveryFailure::InvalidRecoveryMap) => {
            return Err(CubicEqualityFailure::RecoveryVerification {
                evidence: Box::new(RecoveryVerificationFailureEvidence::early(
                    RecoveryVerificationFailureReason::InvalidRecoveryMap,
                )),
                representation: Box::new(representation_evidence),
                backend: Box::new(backend),
            });
        }
        Err(CubicRepresentationRecoveryFailure::Representation(failure)) => {
            return Err(CubicEqualityFailure::Representation(Box::new(failure)));
        }
    };
    let mut field = recovered_representation.field;
    let side_condition = recovered_representation.side_condition;
    let polynomial_round_trip_error = recovered_representation.polynomial_round_trip_error;
    let field_coefficient_round_trip_error =
        recovered_representation.field_coefficient_round_trip_error;
    let field_energy = recovered_representation.field_energy;
    let recovered_energy = recovered_representation.recovered_energy;
    let field_energy_round_trip_error = recovered_representation.field_energy_round_trip_error;
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
            let whitened_residual =
                dense_matrix_vector_product(&block.whitening, dimension, dimension, &residual);
            let recovered_residual = dense_matrix_vector_product(
                &block.inverse_whitening,
                dimension,
                dimension,
                &whitened_residual,
            );
            let whitening_round_trip_error = relative_slice_error(&recovered_residual, &residual);
            RecoveredSoftObjective {
                canonical_indices: block.canonical_indices.clone(),
                loss: objective.loss.clone(),
                covariance_group: objective.covariance_group.clone(),
                block_kind: objective.block_kind.clone(),
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
    let hard_equality_violations = FunctionalViolationEnvelope::from_dimensioned_residuals(
        hard_equalities
            .iter()
            .map(|equality| (equality.dimension, equality.residual)),
    );
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
                .standard_rows
                .iter()
                .zip(&objective.targets)
                .map(|(row, target)| dot_product(row, &backend.candidate) - target)
                .collect::<Vec<_>>();
            let weighted = dense_matrix_vector_product(
                &objective.standard_precision,
                dimension,
                dimension,
                &residual,
            );
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
        objective_round_trip_error <= EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit;
    let relation_tolerances = canonical_relation_tolerances(
        characteristic_length,
        &problem,
        field_energy,
        &assembly,
        &backend,
    );
    let all_source_recovery = solver_form.verify_all_source_recovery(
        &hard_equalities,
        &relation_tolerances,
        &[],
        &soft_equalities,
        &soft_objectives,
    );
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
    if !all_source_recovery.verified {
        reasons.push(RecoveryVerificationFailureReason::SourceCoverageMismatch);
    }
    if !side_condition.is_within_policy() {
        reasons.push(RecoveryVerificationFailureReason::SideConditionViolation);
    }
    if side_condition.round_trip_error > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::SideConditionRoundTripViolation);
    }
    if !hard_residuals_within_tolerance(&problem, &hard_equalities, &relation_tolerances) {
        reasons.push(RecoveryVerificationFailureReason::HardEqualityViolation);
    }
    if polynomial_round_trip_error > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::PolynomialRoundTripViolation);
    }
    if field_coefficient_round_trip_error > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::FieldCoefficientRoundTripViolation);
    }
    if field_energy_round_trip_error > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::FieldEnergyRoundTripViolation);
    }
    if whitening_round_trip_error > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::WhiteningRoundTripViolation);
    }
    if !objective_verified {
        reasons.push(RecoveryVerificationFailureReason::ObjectiveRoundTripViolation);
    }
    if tolerance_round_trip_error > EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit {
        reasons.push(RecoveryVerificationFailureReason::ToleranceRoundTripViolation);
    }
    let query_field_scale =
        canonical_characteristic_field_scale(&problem, characteristic_length, field_energy);
    let query_response_round_trip_error = if reasons.is_empty() {
        match field.finalize_verified_query_representation(
            query_field_scale,
            characteristic_length,
            polynomial_round_trip_error.max(field_coefficient_round_trip_error),
        ) {
            Ok(error) if error <= EXECUTED_NUMERICAL_POLICY.recovery_round_trip_limit => {
                Some(error)
            }
            Ok(_) | Err(_) => {
                reasons
                    .push(RecoveryVerificationFailureReason::QueryRepresentationRoundTripViolation);
                None
            }
        }
    } else {
        None
    };
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
            representation: Box::new(representation_evidence.clone()),
            backend: Box::new(backend),
        });
    }

    Ok(CubicEqualitySolution {
        representation: representation_evidence,
        assembly,
        hard_recovery: solver_form.hard_recovery.clone(),
        all_source_recovery,
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
        query_response_round_trip_error: query_response_round_trip_error
            .expect("accepted recovery verified the query response round trip"),
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

fn hard_residuals_within_tolerance(
    problem: &CubicCanonicalProblem,
    recovered: &[RecoveredHardEquality],
    tolerances: &[CanonicalRelationToleranceEvidence],
) -> bool {
    let mut block_membership = vec![false; recovered.len()];
    for block in &problem.hard_residual_blocks {
        for index in block.canonical_indices() {
            block_membership[*index] = true;
        }
        let scale = block
            .canonical_indices()
            .iter()
            .map(|index| recovered[*index].residual.abs())
            .fold(0.0_f64, f64::max);
        let residual_norm = if scale == 0.0 {
            0.0
        } else {
            scale
                * block
                    .canonical_indices()
                    .iter()
                    .map(|index| (recovered[*index].residual / scale).powi(2))
                    .sum::<f64>()
                    .sqrt()
        };
        let tolerance = block
            .canonical_indices()
            .iter()
            .map(|index| tolerances[*index].physical_tolerance)
            .fold(0.0_f64, f64::max);
        if residual_norm > tolerance {
            return false;
        }
    }
    recovered.iter().zip(tolerances).zip(block_membership).all(
        |((relation, tolerance), in_block)| {
            in_block || relation.residual.abs() <= tolerance.physical_tolerance
        },
    )
}

fn canonical_relation_tolerances(
    characteristic_length: f64,
    problem: &CubicCanonicalProblem,
    field_energy: f64,
    assembly: &EqualityAssemblyEvidence,
    backend: &KktSolveEvidence,
) -> Vec<CanonicalRelationToleranceEvidence> {
    let field_value_gauge_offset = canonical_gauge_offset(problem, FunctionalDimension::FieldValue);
    let derivative_gauge_offset =
        canonical_gauge_offset(problem, FunctionalDimension::FieldValuePerLength);
    let field_scale =
        canonical_characteristic_field_scale(problem, characteristic_length, field_energy);
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
                characteristic_length,
                field_value_gauge_offset,
                derivative_gauge_offset,
            );
            let relation_reference_scale =
                (equality.target - equality.constant_shift_response() * gauge_offset).abs();
            let physical_tolerance = EXECUTED_NUMERICAL_POLICY
                .canonical_characteristic_tolerance_multiplier
                * characteristic_scale
                + EXECUTED_NUMERICAL_POLICY.canonical_relation_reference_tolerance_multiplier
                    * relation_reference_scale;
            let standard_tolerance = physical_tolerance;
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
            let recovered_physical_tolerance = recovered_standard_tolerance;
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

pub(crate) fn canonical_gauge_offset(
    problem: &CubicCanonicalProblem,
    dimension: FunctionalDimension,
) -> f64 {
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
        .chain(
            problem
                .affine_inequalities
                .iter()
                .filter(|inequality| inequality.dimension == dimension)
                .map(|inequality| (inequality.constant_shift_response(), inequality.bound)),
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

pub(crate) fn canonical_characteristic_field_scale(
    problem: &CubicCanonicalProblem,
    length: f64,
    field_energy: f64,
) -> f64 {
    let native_energy = field_energy / problem.field_energy_normalization.factor();
    let energy_scale = (native_energy.abs() * length.powi(3)).sqrt();
    let relation_scale = problem
        .equalities
        .iter()
        .map(|relation| {
            let gauge = canonical_gauge_offset(problem, relation.dimension());
            let scale = (relation.target() - relation.constant_shift_response() * gauge).abs();
            match relation.dimension() {
                FunctionalDimension::FieldValue => scale,
                FunctionalDimension::FieldValuePerLength => length * scale,
            }
        })
        .chain(problem.soft_equalities.iter().map(|relation| {
            let gauge = canonical_gauge_offset(problem, relation.dimension());
            let scale = (relation.target() - relation.constant_shift_response() * gauge).abs();
            match relation.dimension() {
                FunctionalDimension::FieldValue => scale,
                FunctionalDimension::FieldValuePerLength => length * scale,
            }
        }))
        .chain(problem.affine_inequalities.iter().map(|relation| {
            let gauge = canonical_gauge_offset(problem, relation.dimension());
            let scale = (relation.bound() - relation.constant_shift_response() * gauge).abs();
            match relation.dimension() {
                FunctionalDimension::FieldValue => scale,
                FunctionalDimension::FieldValuePerLength => length * scale,
            }
        }))
        .fold(0.0_f64, f64::max);
    let soft_loss_scale = problem
        .soft_objectives
        .iter()
        .map(|objective| {
            let reference_scale = objective.loss().residual_reference_scale();
            let is_derivative = objective.residuals().first().is_some_and(|residual| {
                problem.soft_equalities.iter().any(|relation| {
                    relation.provenance().residual() == residual
                        && relation.dimension() == FunctionalDimension::FieldValuePerLength
                })
            });
            if is_derivative {
                length * reference_scale
            } else {
                reference_scale
            }
        })
        .chain(problem.affine_inequalities.iter().filter_map(|relation| {
            relation.violation_channel().map(|channel| {
                let reference_scale = channel.loss().residual_reference_scale();
                match relation.dimension() {
                    FunctionalDimension::FieldValue => reference_scale,
                    FunctionalDimension::FieldValuePerLength => length * reference_scale,
                }
            })
        }))
        .fold(0.0_f64, f64::max);
    energy_scale.max(relation_scale).max(soft_loss_scale)
}

pub(crate) fn relative_slice_error(actual: &[f64], expected: &[f64]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(0.0_f64, f64::max)
}

pub(crate) fn dense_matrix_vector_product(
    matrix: &[f64],
    rows: usize,
    columns: usize,
    vector: &[f64],
) -> Vec<f64> {
    debug_assert_eq!(matrix.len(), rows * columns);
    debug_assert_eq!(vector.len(), columns);
    (0..rows)
        .map(|row| {
            (0..columns)
                .map(|column| matrix[row * columns + column] * vector[column])
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

pub(crate) fn dot_product(left: &[f64], right: &[f64]) -> f64 {
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
            hard_residual_blocks: Vec::new(),
            affine_inequalities: Vec::new(),
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
        let evidence = &representation.evidence;

        assert_eq!(evidence.fitting_functional_count, 10);
        assert_eq!(evidence.polynomial_dimension, 4);
        assert_eq!(evidence.polynomial_rank, 4);
        assert!(evidence.polynomial_rrqr_ratio > evidence.polynomial_rank_accept_ratio);
        assert!(evidence.polynomial_svd_ratio > evidence.polynomial_rank_accept_ratio);
        assert_eq!(representation.kernel.shape(), (10, 10));
        assert_eq!(representation.polynomial.shape(), (10, 4));
        assert!(evidence.quotient_construction.null_space_defect <= 1.0e-12);
        assert!(evidence.reduced_smallest_singular_value > 0.0);
        assert!(evidence.reduced_symmetry_defect <= evidence.symmetry_defect_limit);
        assert!(evidence.affine_reproduction_error <= 1.0e-11);
        let normalization = representation.field_energy_normalization;
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

        let standard = representation.coordinates.to_standard(TRUTH_POLYNOMIAL);
        let round_trip = representation.coordinates.to_physical(standard);
        for (actual, expected) in round_trip.into_iter().zip(TRUTH_POLYNOMIAL) {
            assert_close(actual, expected, 1.0e-15);
        }
    }

    #[test]
    fn implicit_quotient_retains_nearby_same_direction_functionals() {
        let mut distinct_uses = usages()[..4].to_vec();
        let nearby_x = f64::from_bits(0.25_f64.to_bits() + 1);
        for (index, support) in [[0.25, 0.5, 0.75], [nearby_x, 0.5, 0.75]]
            .into_iter()
            .enumerate()
        {
            distinct_uses.push(FunctionalUse::new(
                functional(support, 0.0, [1.0, 0.0, 0.0]),
                UsageProvenance::new(
                    SourceId::new(format!("nearby-direction-{index}")),
                    None,
                    RelationId::new(format!("nearby-direction-relation-{index}")),
                    ResidualId::new(format!("nearby-direction-residual-{index}")),
                    SemanticRolePath::new("nearby-same-direction"),
                ),
            ));
        }
        let problem = canonical_problem(distinct_uses, vec![0.0; 6]);

        let fitting_uses = canonical_fitting_uses(
            &problem.equalities,
            &problem.soft_equalities,
            &problem.affine_inequalities,
        );
        assert_eq!(fitting_uses.len(), 6);
        let functionals = fitting_uses
            .iter()
            .map(|usage| usage.functional().clone())
            .collect::<Vec<_>>();
        let (_, _, polynomial) = assemble_polynomial_pairing(&functionals)
            .expect("nearby finite supports preserve the complete polynomial pairing");
        let complement = HouseholderNullSpace::new(&polynomial, POLYNOMIAL_DIMENSION)
            .expect("the full-rank pairing has an implicit complement");

        assert_eq!(polynomial.shape(), (6, POLYNOMIAL_DIMENSION));
        assert_eq!(complement.reduced_dimension(), 2);
    }

    #[test]
    fn energy_orthonormal_quotient_is_verified_without_regularization_or_spectrum_gate() {
        let problem = canonical_problem(usages(), MANUFACTURED_TARGETS.to_vec());
        let (representation, field_form) = CubicRepresentation::build(
            canonical_fitting_uses(
                &problem.equalities,
                &problem.soft_equalities,
                &problem.affine_inequalities,
            ),
            GlobalAnisotropyMetric::identity(),
            FieldEnergyNormalization::all_hard(),
        )
        .expect("the manufactured quotient has a verified energy basis");
        let factorization = &representation.evidence.quotient_factorization;

        assert_eq!(factorization.quotient_dimension, 6);
        assert_eq!(factorization.retained_modes, 6);
        assert_eq!(factorization.truncated_modes, 0);
        assert_eq!(factorization.unregularized_llt_count, 1);
        assert_eq!(factorization.full_spectrum_analysis_count, 0);
        assert!(factorization.normalized_backward_error <= 1.0e-11);
        assert_eq!(factorization.pivot_intervals.len(), 6);
        assert!(
            factorization
                .pivot_intervals
                .iter()
                .all(|interval| interval.lower > 0.0)
        );
        assert!(
            factorization
                .field_energy_identity_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(
            factorization
                .side_condition_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(
            factorization
                .recovery_round_trip_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(
            factorization
                .canonical_response_round_trip_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(!factorization.kernel_ridge_applied);
        assert!(!factorization.gram_jitter_applied);
        assert!(!factorization.mode_truncation_applied);

        let solver_form = CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
            .expect("the canonical solver form consumes the energy basis");
        let quotient_dimension = representation.null_space.reduced_dimension();
        let quotient_field_energy = solver_form.field_energy(CubicFieldCoordinateLayout::Quotient);
        assert_eq!(quotient_field_energy.len(), quotient_dimension.pow(2));
        for row in 0..quotient_dimension {
            for column in 0..quotient_dimension {
                assert_eq!(
                    quotient_field_energy[row * quotient_dimension + column],
                    f64::from(row == column)
                );
            }
        }
    }

    #[test]
    fn cancellation_ambiguous_pivot_is_rescued_and_reattached() {
        let next_after_one = f64::from_bits(1.0_f64.to_bits() + 1);
        let gram = DenseMatrix::from_fn(2, 2, |row, column| match (row, column) {
            (0, 0) | (0, 1) | (1, 0) => 1.0,
            (1, 1) => next_after_one,
            _ => unreachable!(),
        });

        let basis = EnergyOrthonormalQuotientBasis::factor(
            &gram,
            &DenseMatrix::from_fn(2, POLYNOMIAL_DIMENSION, |_, _| 0.0),
            &gram,
            Some(&DirectPrecisionRescueSource::from_gram(&gram)),
        )
        .expect("double-double proves and reattaches the cancellation-scale mode");
        let rescue = basis
            .evidence
            .precision_rescue
            .expect("the actual upgraded range is recorded");
        assert_eq!(rescue.first_mode, 1);
        assert_eq!(rescue.mode_count, 1);
        assert_eq!(rescue.precision_bits, 106);
        assert_eq!(rescue.conclusion, PrecisionRescueConclusion::Positive);
        assert_eq!(basis.evidence.retained_modes, 2);
        assert!(basis.evidence.pivot_intervals[1].lower > 0.0);
        assert!(
            basis
                .evidence
                .field_energy_identity_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(
            basis
                .evidence
                .side_condition_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(
            basis
                .evidence
                .recovery_round_trip_error
                .is_some_and(|error| error <= 1.0e-11)
        );
        assert!(
            basis
                .evidence
                .canonical_response_round_trip_error
                .is_some_and(|error| error <= 1.0e-11)
        );
    }

    #[test]
    fn quotient_rescue_keeps_zero_negative_and_gray_conclusions_distinct() {
        let trailing = DenseMatrix::from_fn(1, POLYNOMIAL_DIMENSION, |_, _| 0.0);
        let zero = DenseMatrix::from_fn(1, 1, |_, _| 0.0);

        let algebraic = DirectPrecisionRescueSource::with_algebraic_zero(&zero);
        assert!(matches!(
            EnergyOrthonormalQuotientBasis::factor(&zero, &trailing, &zero, Some(&algebraic)),
            Err(RepresentationFailure::QuotientRankDeficient {
                evidence: PrecisionRescueEvidence {
                    conclusion: PrecisionRescueConclusion::AlgebraicZero,
                    ..
                },
                ..
            })
        ));

        let unresolved = DirectPrecisionRescueSource::from_gram(&zero);
        assert!(matches!(
            EnergyOrthonormalQuotientBasis::factor(&zero, &trailing, &zero, Some(&unresolved)),
            Err(RepresentationFailure::QuotientPrecisionRescueGrayZone {
                evidence: PrecisionRescueEvidence {
                    conclusion: PrecisionRescueConclusion::GrayZone,
                    ..
                },
                ..
            })
        ));

        let negative = DenseMatrix::from_fn(1, 1, |_, _| -1.0);
        let corrupted = DirectPrecisionRescueSource::from_gram(&negative);
        assert!(matches!(
            EnergyOrthonormalQuotientBasis::factor(
                &negative,
                &trailing,
                &negative,
                Some(&corrupted),
            ),
            Err(RepresentationFailure::QuotientNegativeCurvature { .. })
        ));
    }

    #[test]
    fn quotient_rescue_rejects_all_65_ambiguous_modes_without_truncation() {
        let dimension = 65;
        let zero = DenseMatrix::from_fn(dimension, dimension, |_, _| 0.0);
        let source = DirectPrecisionRescueSource::from_gram(&zero);
        let failure = EnergyOrthonormalQuotientBasis::factor(
            &zero,
            &DenseMatrix::from_fn(dimension, POLYNOMIAL_DIMENSION, |_, _| 0.0),
            &zero,
            Some(&source),
        )
        .expect_err("65 unresolved modes exceed the bounded rescue policy");
        match failure {
            RepresentationFailure::QuotientPrecisionRescueGrayZone {
                quotient_dimension,
                evidence,
                execution,
            } => {
                assert_eq!(quotient_dimension, dimension);
                assert_eq!(evidence.first_mode, 0);
                assert_eq!(evidence.mode_count, dimension);
                assert_eq!(evidence.precision_bits, 106);
                assert_eq!(
                    evidence.conclusion,
                    PrecisionRescueConclusion::CapacityExceeded,
                );
                assert_eq!(execution, AnalysisExecutionEvidence::pre_backend());
            }
            other => panic!("unexpected bounded-rescue result: {other:?}"),
        }
    }

    #[test]
    fn solve_coordinates_are_deterministic_and_recover_physical_field_coefficients() {
        let representation = CubicRepresentation::new(usages(), GlobalAnisotropyMetric::identity())
            .expect("the manufactured representer span is Cubic-admissible");
        let transform = representation.coordinates;

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
        let capacity_failure = plan_equality_capacity_after_representation(
            usize::MAX,
            usize::MAX,
            &solution.representation,
        )
        .expect_err("overflowing capacity must retain completed representation evidence");
        let CubicEqualityFailure::Representation(failure) = capacity_failure else {
            panic!("capacity rejection must remain a representation failure");
        };
        let build = failure
            .build_evidence()
            .expect("completed representation evidence must survive capacity rejection");
        assert_eq!(build.failure_stage, RepresentationBuildStage::Backend);
        assert_eq!(
            build.last_completed_stage,
            RepresentationBuildStage::ResponseAssembly
        );
        assert_eq!(
            build
                .quotient_construction
                .as_ref()
                .map(|evidence| evidence.quotient_dimension),
            Some(6)
        );
        assert!(build.quotient_factorization.is_some());
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
            hard_residual_blocks: Vec::new(),
            affine_inequalities: Vec::new(),
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
            CubicEqualityFailure::DirectInputConflict { evidence, .. } => {
                assert_eq!(evidence.canonical_residual, 0.0);
                assert!(evidence.separation_margin > 0.0);
                assert_eq!(
                    evidence
                        .relations
                        .iter()
                        .map(|relation| relation.provenance.source().as_str())
                        .collect::<Vec<_>>(),
                    [
                        "manufactured-conflicting-gauge",
                        "manufactured-latent-gauge"
                    ]
                );
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
                    precision_rescue: None,
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
    fn preflight_polynomial_analysis_rescues_a_small_positive_pi1_mode() {
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

        assert_eq!(preflight_polynomial_analysis_failure(&fitting_uses), None);

        let functionals = fitting_uses
            .iter()
            .map(|usage| usage.functional().clone())
            .collect::<Vec<_>>();
        let (coordinates, _, polynomial) =
            assemble_polynomial_pairing(&functionals).expect("the canonical Pi1 pairing assembles");
        let (_, rank, evidence) = verify_polynomial_rank(&polynomial, &functionals, &coordinates)
            .expect("double-double proves the small fourth Pi1 mode positive");
        assert_eq!(rank, POLYNOMIAL_DIMENSION);
        let rescue = evidence
            .precision_rescue
            .expect("the actual polynomial upgrade is recorded");
        assert_eq!(rescue.first_mode, 0);
        assert_eq!(rescue.mode_count, POLYNOMIAL_DIMENSION);
        assert_eq!(rescue.precision_bits, 106);
        assert_eq!(rescue.conclusion, PrecisionRescueConclusion::Positive);
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
    fn verified_llt_retains_small_positive_modes_and_distinguishes_nonpositive_pivots() {
        let trailing_polynomial = DenseMatrix::from_fn(2, POLYNOMIAL_DIMENSION, |_, _| 0.0);
        let small_positive = DenseMatrix::from_fn(2, 2, |row, column| {
            if row == column {
                if row == 0 { 1.0 } else { 1.0e-30 }
            } else {
                0.0
            }
        });
        let basis = EnergyOrthonormalQuotientBasis::factor(
            &small_positive,
            &trailing_polynomial,
            &small_positive,
            None,
        )
        .expect("a reliably positive mode is retained regardless of relative scale");
        assert_eq!(basis.evidence.retained_modes, 2);
        assert_eq!(basis.evidence.truncated_modes, 0);
        assert_eq!(basis.evidence.full_spectrum_analysis_count, 0);

        for (gram, expected_pivot) in [
            (
                DenseMatrix::from_fn(2, 2, |row, column| {
                    if row == column {
                        if row == 0 { 1.0 } else { -1.0 }
                    } else {
                        0.0
                    }
                }),
                -1.0,
            ),
            (
                DenseMatrix::from_fn(2, 2, |row, column| f64::from(row == column && row == 0)),
                0.0,
            ),
        ] {
            match EnergyOrthonormalQuotientBasis::factor(&gram, &trailing_polynomial, &gram, None) {
                Err(RepresentationFailure::QuotientFactorizationNotPositive {
                    quotient_dimension,
                    pivot_index,
                    interval,
                    ..
                }) => {
                    assert_eq!(quotient_dimension, 2);
                    assert_eq!(pivot_index, 1);
                    assert_eq!(interval.lower, expected_pivot);
                    assert_eq!(interval.upper, expected_pivot);
                }
                other => panic!("unexpected non-positive factorization result: {other:?}"),
            }
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
            EXECUTED_NUMERICAL_POLICY.side_condition_limit
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
        let field_form = representation.solver_field_form();
        let solver_form = CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
            .expect("the canonical solver form is valid");
        let (assembly, backend) = solve_standard_form(&solver_form)
            .expect("the undamaged standard form should produce a backend candidate");
        let mut damaged = representation;
        damaged.coordinates.length = 0.0;

        let failure = recover_and_verify(damaged, solver_form, problem, assembly, backend)
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
        let field_form = representation.solver_field_form();
        let solver_form = CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
            .expect("the canonical solver form is valid");
        let (mut assembly, backend) = solve_standard_form(&solver_form)
            .expect("the standard form should still satisfy its backend contract");
        assembly.hard_equality_rows[0].provenance = UsageProvenance::new(
            SourceId::new("corrupted-source"),
            None,
            RelationId::new("corrupted-relation"),
            ResidualId::new("corrupted-residual"),
            SemanticRolePath::new("corrupted-role"),
        );
        let failure = recover_and_verify(representation, solver_form, problem, assembly, backend)
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
        let field_form = representation.solver_field_form();
        let solver_form = CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
            .expect("the canonical solver form is valid");
        let (mut assembly, backend) = solve_standard_form(&solver_form)
            .expect("the standard form should produce a verified backend candidate");
        assembly.hard_equality_rows.swap(0, 1);

        let failure = recover_and_verify(representation, solver_form, problem, assembly, backend)
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
        let (representation, field_form) = CubicRepresentation::build(
            canonical_fitting_uses(
                &problem.equalities,
                &problem.soft_equalities,
                &problem.affine_inequalities,
            ),
            GlobalAnisotropyMetric::identity(),
            problem.field_energy_normalization,
        )
        .expect("the mixed hard/soft representation is valid");
        let solver_form = CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
            .expect("the canonical solver form is valid");
        let (mut assembly, backend) = solve_standard_form(&solver_form)
            .expect("the undamaged objective should produce a backend candidate");
        let residual = dot_product(
            &assembly.soft_objective_blocks[0].standard_rows[0],
            &backend.candidate,
        ) - soft_target;
        assert!(residual.abs() <= 1.0e-11, "residual={residual:e}");
        assembly.soft_objective_blocks[0].standard_precision[0] *= 2.0;

        let failure = recover_and_verify(representation, solver_form, problem, assembly, backend)
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
    fn damaged_whitening_objective_association_is_rejected_without_a_model() {
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
                CanonicalSoftResidualBlockKind::CovarianceGroup {
                    members: vec![
                        CanonicalSoftResidualMemberKind::FieldValue,
                        CanonicalSoftResidualMemberKind::FieldValue,
                    ],
                },
            ));
        problem.soft_equalities = soft_equalities;
        let (representation, field_form) = CubicRepresentation::build(
            canonical_fitting_uses(
                &problem.equalities,
                &problem.soft_equalities,
                &problem.affine_inequalities,
            ),
            GlobalAnisotropyMetric::identity(),
            problem.field_energy_normalization,
        )
        .expect("the mixed hard/soft representation is valid");
        let mut solver_form =
            CanonicalCubicSolverForm::assemble(&representation, field_form, &problem)
                .expect("the canonical solver form is valid");
        let (mut assembly, backend) = solve_standard_form(&solver_form)
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
        solver_form.soft_objectives[0].inverse_whitening[0] *= 2.0;
        if let CanonicalSoftLoss::Covariance {
            inverse_whitening, ..
        } = &mut solver_form.soft_objectives[0].loss
        {
            inverse_whitening[0] *= 2.0;
        } else {
            panic!("the fixture uses covariance whitening");
        }

        let failure = recover_and_verify(representation, solver_form, problem, assembly, backend)
            .expect_err("a non-invertible whitening recovery map must not produce a model");
        match failure {
            CubicEqualityFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![RecoveryVerificationFailureReason::ProvenanceMismatch]
                );
                assert_eq!(evidence.whitening_round_trip_error, None);
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }
}
