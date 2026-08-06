//! Strict fit outcomes and physical-unit fit reports.

use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fmt;

use crate::capacity::{
    CapacityExceededEvidence, CapacityExceededReason, EqualityCapacityShape, SourceStorageShape,
    plan_equality_capacity, plan_equality_capacity_for,
};
use crate::clarabel_backend::{
    ClarabelAttemptEvidence, ClarabelAttemptProfile, ClarabelTermination,
};
use crate::cubic_equality::{
    AlgebraicAnalysisStage as InternalCubicAnalysisStage, CanonicalAffineInequality,
    CanonicalEqualityParticipation, CanonicalHardEquality, CanonicalHardResidualBlock,
    CanonicalInequalitySense, CanonicalRelationToleranceEvidence, CanonicalSoftEquality,
    CanonicalSoftLoss, CanonicalSoftObjective, CanonicalSoftResidualBlockKind,
    CanonicalSoftResidualMemberKind, CanonicalViolationChannel, CanonicalViolationLoss,
    CpdEvidence, CubicCanonicalProblem, CubicEqualityFailure, CubicEqualitySolution,
    PhysicalSideConditionEvidence, RecoveryVerificationFailureEvidence, RepresentationFailure,
    SemanticLatentCoefficient, SemanticLatentDefinition,
    SolveCoordinateTransformFailureReason as InternalSolveCoordinateFailure,
    canonical_fitting_uses, preflight_polynomial_analysis_failure,
};
use crate::cubic_execution::{
    CubicExecutionCore, CubicExecutionFailure, CubicExecutionSolution, QpAttemptFailureReason,
    QpAttemptRecord, RecoveredAffineInequality, ValidatedInfeasibilityEvidence,
    ValidatedRecessionEvidence,
};
use crate::diagnostics::{
    AnalysisContractQuantity, AnalysisFailureEvidence, AnalysisFailureStage,
    AttemptFailureCategory, AttemptFailureEvidence, BackendAttemptSettings, BackendFingerprint,
    BackendFingerprintParts, BackendInputField, CanonicalAcceptanceEvidence,
    CanonicalAcceptanceEvidenceParts, CanonicalEvidenceSource, CapacityEvidence,
    CapacityFailureKind, ConvexResidualEvidence as PublicConvexResidualEvidence,
    ConvexResidualEvidenceParts, CubicAnalysisEvidence, CubicAnalysisEvidenceParts,
    CubicLltPivotInterval, CubicLltPivotIntervalParts, CubicQuotientConstructionEvidence,
    CubicQuotientConstructionEvidenceParts, CubicQuotientFactorizationEvidence,
    CubicQuotientFactorizationEvidenceParts, DirectInputConflictEvidence, InertiaCounts,
    InertiaEvidence, InfeasibilityCertificateEvidence, InfeasibilityCertificateEvidenceParts,
    InterpretableRankDeficiencyEvidence, InterpretableRankDeficiencyEvidenceParts,
    LinearResidualEvidence, ProblemDiagnosis, RankDecision, RankDeficiencyConcept, RankEvidence,
    RankEvidenceDomain, RankEvidenceParts, RecessionRayEvidence, RecessionRayEvidenceParts,
    RecoveryVerificationEvidence, RecoveryVerificationEvidenceParts, RelationGraphConflictEvidence,
    ResidualDimension, ScalingFailureReason, ScalingSummary, SharedLevelSetConflictSourceEvidence,
    SharedLevelSetRelationConflictEvidence, SideConditionEvidence, SolveAttemptKind,
    SolveAttemptRecord, SolveAttemptRecordParts, SolveAttemptTermination,
    SolveCoordinateFailureReason, UnidentifiedAdditiveGaugeEvidence,
    UninformativeSharedLevelSetEvidence, UnresolvedAxialNormalEvidence,
};
use crate::functional::{
    CanonicalFunctional, FunctionalDimension, FunctionalRepresenterSpan, FunctionalTerm,
    FunctionalUse, GroupId, RelationId, ResidualId, SemanticRolePath, SourceId, UsageProvenance,
};
use crate::geometry::{FieldUnitLabel, Point3, Vector3};
use crate::kernel::{FieldEnergyNormalization, KernelConfig};
use crate::kkt::{
    AlgebraicAnalysisPhase as InternalBackendAnalysisPhase,
    AlgebraicRankEvidence as InternalRankEvidence,
    BackendAttemptSettings as InternalAttemptSettings,
    BackendContractViolationReason as InternalBackendContractReason,
    BackendFingerprint as InternalBackendFingerprint, Inertia as InternalInertia,
    KktAttemptFailureReason, KktAttemptKind, KktAttemptRecord, KktFailure,
    KktInputField as InternalBackendInputField, NumericalFailureReason,
    RankClassification as InternalRankDecision, RuizScalingFailure as InternalScalingFailure,
    SolveAttemptTermination as InternalTermination, WorkspacePhase as InternalWorkspacePhase,
};
use crate::model::SolvedModel;
use crate::numerical::NumericalPolicyId;
use crate::observation::{
    CovarianceGroupMember, CovarianceMatrix, FieldValueConfiguration, GradientConfiguration,
    MinimumNormalSlope, MinimumNormalSlopeConfiguration, MinimumNormalSlopeEnforcement,
    NormalDirectionConfiguration, NormalDirectionEnforcement, ObservationInput, QuadraticPenalty,
    StandardDeviation, TangentConfiguration,
};
use crate::problem::{ProblemSnapshot, ThreadBudget};
use crate::relation::{
    AdditiveFieldGaugeReference, AffineBoundConfiguration, AffineBoundSide,
    FieldSeparationInterval, LinearViolationPenalty, MinimumFieldOffset, MinimumFieldSeparation,
    PointToLevelSetRelation, PointToLevelSetSide, PolaritySelection, SharedLevelSetRelationInput,
    SharedLevelSetRelationKind, SharedLevelSetRelationOrientation, StratigraphicFieldDirection,
};

/// Successful fit output: one accepted model and its complete report.
#[derive(Debug)]
pub struct FitSuccess {
    model: SolvedModel,
    report: FitReport,
}

impl FitSuccess {
    /// Returns the immutable accepted model.
    pub fn model(&self) -> &SolvedModel {
        &self.model
    }

    /// Returns the report for the successful fit.
    pub fn report(&self) -> &FitReport {
        &self.report
    }

    /// Splits the outcome into owning model and report values.
    pub fn into_parts(self) -> (SolvedModel, FitReport) {
        (self.model, self.report)
    }
}

/// Failed fit output. It never contains a model or candidate field.
#[derive(Debug)]
pub struct FitFailure {
    diagnosis: ProblemDiagnosis,
    report: Box<FitReport>,
}

impl FitFailure {
    /// Returns GeoRBF's structured semantic diagnosis.
    pub fn diagnosis(&self) -> ProblemDiagnosis {
        self.diagnosis
    }

    /// Returns the report accumulated for the failed fit.
    pub fn report(&self) -> &FitReport {
        &self.report
    }
}

impl fmt::Display for FitFailure {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "GeoRBF fit failed: {:?}", self.diagnosis)
    }
}

impl Error for FitFailure {}

/// Public problem-size evidence recorded before solver presolve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProblemSize {
    input_observations: usize,
    scalar_hard_relations: usize,
    scalar_soft_relations: usize,
    canonical_hard_equalities: Option<usize>,
    canonical_soft_equalities: Option<usize>,
    center_coefficients: Option<usize>,
    semantic_latents: usize,
    auxiliary_variables: usize,
    quadratic_objective_terms: usize,
    linear_objective_terms: usize,
    affine_inequality_constraints: usize,
    cone_blocks: usize,
    primal_variables: Option<usize>,
    equality_constraints: Option<usize>,
    kkt_dimension: Option<usize>,
}

struct CubicProblemSizeParts {
    input_observations: usize,
    scalar_hard_relations: usize,
    scalar_soft_relations: usize,
    canonical_hard_equalities: usize,
    canonical_soft_equalities: usize,
    quadratic_objective_terms: usize,
    linear_objective_terms: usize,
    affine_inequality_constraints: usize,
    center_coefficients: usize,
    semantic_latents: usize,
    solver_hard_equalities: usize,
    auxiliary_variables: usize,
}

impl ProblemSize {
    fn cubic_equality(parts: CubicProblemSizeParts) -> Option<Self> {
        let primal_variables = parts
            .center_coefficients
            .checked_add(4)?
            .checked_add(parts.semantic_latents)?;
        let equality_constraints = parts.solver_hard_equalities.checked_add(4)?;
        let kkt_dimension = primal_variables.checked_add(equality_constraints)?;
        Some(Self {
            input_observations: parts.input_observations,
            scalar_hard_relations: parts.scalar_hard_relations,
            scalar_soft_relations: parts.scalar_soft_relations,
            canonical_hard_equalities: Some(parts.canonical_hard_equalities),
            canonical_soft_equalities: Some(parts.canonical_soft_equalities),
            center_coefficients: Some(parts.center_coefficients),
            semantic_latents: parts.semantic_latents,
            auxiliary_variables: parts.auxiliary_variables,
            quadratic_objective_terms: parts.quadratic_objective_terms,
            linear_objective_terms: parts.linear_objective_terms,
            affine_inequality_constraints: parts.affine_inequality_constraints,
            cone_blocks: 0,
            primal_variables: Some(primal_variables),
            equality_constraints: Some(equality_constraints),
            kkt_dimension: Some(kkt_dimension),
        })
    }

    fn cubic_qp(parts: CubicProblemSizeParts) -> Option<Self> {
        let primal_variables = parts
            .center_coefficients
            .checked_add(parts.semantic_latents)?
            .checked_add(parts.auxiliary_variables)?;
        Some(Self {
            input_observations: parts.input_observations,
            scalar_hard_relations: parts.scalar_hard_relations,
            scalar_soft_relations: parts.scalar_soft_relations,
            canonical_hard_equalities: Some(parts.canonical_hard_equalities),
            canonical_soft_equalities: Some(parts.canonical_soft_equalities),
            center_coefficients: Some(parts.center_coefficients),
            semantic_latents: parts.semantic_latents,
            auxiliary_variables: parts.auxiliary_variables,
            quadratic_objective_terms: parts.quadratic_objective_terms,
            linear_objective_terms: parts.linear_objective_terms,
            affine_inequality_constraints: parts.affine_inequality_constraints,
            cone_blocks: 0,
            primal_variables: Some(primal_variables),
            equality_constraints: Some(parts.solver_hard_equalities),
            kkt_dimension: None,
        })
    }

    /// Returns the independent top-level input count for audit context.
    pub fn input_observations(self) -> usize {
        self.input_observations
    }

    /// Returns scalar hard components across every caller-owned source.
    pub fn scalar_hard_relations(self) -> usize {
        self.scalar_hard_relations
    }

    /// Returns independently retained scalar soft residual channels.
    pub fn scalar_soft_relations(self) -> usize {
        self.scalar_soft_relations
    }

    /// Returns the exact hard equalities retained after duplicate merging, or
    /// `None` when fit stopped before canonical lowering.
    pub fn canonical_hard_equalities(self) -> Option<usize> {
        self.canonical_hard_equalities
    }

    /// Returns canonical soft equality residuals, or `None` before lowering.
    pub fn canonical_soft_equalities(self) -> Option<usize> {
        self.canonical_soft_equalities
    }

    /// Returns the representer-center coefficient count, or `None` before
    /// representation planning.
    pub fn center_coefficients(self) -> Option<usize> {
        self.center_coefficients
    }

    /// Returns the semantic latent count.
    pub fn semantic_latents(self) -> usize {
        self.semantic_latents
    }

    /// Returns the auxiliary-variable count.
    pub fn auxiliary_variables(self) -> usize {
        self.auxiliary_variables
    }

    /// Returns quadratic soft-loss terms in the physical objective.
    pub fn quadratic_objective_terms(self) -> usize {
        self.quadratic_objective_terms
    }

    /// Returns linear soft-violation terms in the physical objective.
    pub fn linear_objective_terms(self) -> usize {
        self.linear_objective_terms
    }

    /// Returns scalar affine-inequality rows before backend presolve.
    pub fn affine_inequality_constraints(self) -> usize {
        self.affine_inequality_constraints
    }

    /// Returns the conic block count.
    pub fn cone_blocks(self) -> usize {
        self.cone_blocks
    }

    /// Returns the backend-standard-form primal dimension, or `None` when no
    /// exact standard form was planned.
    pub fn primal_variables(self) -> Option<usize> {
        self.primal_variables
    }

    /// Returns the backend-standard-form equality dimension, or `None` when no
    /// exact standard form was planned.
    pub fn equality_constraints(self) -> Option<usize> {
        self.equality_constraints
    }

    /// Returns the actual symmetric augmented KKT dimension, or `None` when no
    /// exact KKT was planned.
    pub fn kkt_dimension(self) -> Option<usize> {
        self.kkt_dimension
    }
}

/// Physical recovery assessment for one scalar hard-relation component.
#[derive(Debug, Clone, PartialEq)]
pub struct HardRelationAssessment {
    source_id: SourceId,
    group_id: Option<GroupId>,
    semantic_role: SemanticRolePath,
    dimension: ResidualDimension,
    target: f64,
    recovered_value: f64,
    residual: f64,
    tolerance: f64,
    characteristic_scale: f64,
    relation_reference_scale: f64,
    standard_tolerance: f64,
    scaled_kkt_tolerance: Option<f64>,
    recovered_physical_tolerance: f64,
    tolerance_round_trip_error: f64,
}

/// Which side of a Field Value Bound one assessment describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum BoundSide {
    /// A lower relation of the form `bound <= field`.
    Lower,
    /// An upper relation of the form `field <= bound`.
    Upper,
}

/// Whether a recovered bound side lies on its accepted activity envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BoundActiveState {
    /// The satisfaction slack is no larger than the physical relation tolerance.
    Active,
    /// The relation has physical satisfaction slack above its tolerance.
    Inactive,
}

/// Physical recovery assessment for one independently identified bound side.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldValueBoundAssessment {
    source_id: SourceId,
    semantic_role: SemanticRolePath,
    side: BoundSide,
    bound: f64,
    recovered_value: f64,
    slack: f64,
    violation: f64,
    tolerance: f64,
    active_state: BoundActiveState,
    quadratic_penalty: Option<QuadraticPenalty>,
    linear_violation_penalty: Option<LinearViolationPenalty>,
    loss: Option<f64>,
}

/// Physical recovery assessment for one Directional Derivative Interval side.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectionalDerivativeIntervalAssessment {
    source_id: SourceId,
    semantic_role: SemanticRolePath,
    direction: Vector3,
    side: BoundSide,
    bound: f64,
    recovered_directional_derivative: f64,
    slack: f64,
    violation: f64,
    tolerance: f64,
    active_state: BoundActiveState,
    quadratic_penalty: Option<QuadraticPenalty>,
    linear_violation_penalty: Option<LinearViolationPenalty>,
    loss: Option<f64>,
}

/// Physical recovery assessment for one Field Separation Interval side.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldSeparationIntervalAssessment {
    relation: FieldSeparationInterval,
    semantic_role: SemanticRolePath,
    side: BoundSide,
    bound: f64,
    recovered_reference_value: f64,
    recovered_target_value: f64,
    recovered_field_separation: f64,
    slack: f64,
    violation: f64,
    tolerance: f64,
    active_state: BoundActiveState,
    quadratic_penalty: Option<QuadraticPenalty>,
    linear_violation_penalty: Option<LinearViolationPenalty>,
    loss: Option<f64>,
}

/// Physical recovery assessment for one Point to Level Set Relation.
#[derive(Debug, Clone, PartialEq)]
pub struct PointToLevelSetRelationAssessment {
    relation: PointToLevelSetRelation,
    semantic_role: SemanticRolePath,
    recovered_point_value: f64,
    recovered_level_value: f64,
    recovered_field_offset: f64,
    slack: f64,
    violation: f64,
    tolerance: f64,
    active_state: BoundActiveState,
    quadratic_penalty: Option<QuadraticPenalty>,
    linear_violation_penalty: Option<LinearViolationPenalty>,
    loss: Option<f64>,
}

impl FieldSeparationIntervalAssessment {
    /// Returns the caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        self.relation.source_id()
    }

    /// Returns the stable lower- or upper-side semantic role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the ordered reference shared level set.
    pub fn reference_group_id(&self) -> &GroupId {
        self.relation.reference_group_id()
    }

    /// Returns the ordered target shared level set.
    pub fn target_group_id(&self) -> &GroupId {
        self.relation.target_group_id()
    }

    /// Returns which side of the signed interval this assessment describes.
    pub fn side(&self) -> BoundSide {
        self.side
    }

    /// Returns this side's finite bound in field-value units.
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Returns the recovered reference-group field value.
    pub fn recovered_reference_value(&self) -> f64 {
        self.recovered_reference_value
    }

    /// Returns the recovered target-group field value.
    pub fn recovered_target_value(&self) -> f64 {
        self.recovered_target_value
    }

    /// Returns recovered `target - reference` in field-value units.
    pub fn recovered_field_separation(&self) -> f64 {
        self.recovered_field_separation
    }

    /// Returns nonnegative physical satisfaction slack.
    pub fn slack(&self) -> f64 {
        self.slack
    }

    /// Returns nonnegative physical violation.
    pub fn violation(&self) -> f64 {
        self.violation
    }

    /// Returns the physical acceptance tolerance for this side.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns the stable physical active state.
    pub fn active_state(&self) -> BoundActiveState {
        self.active_state
    }

    /// Returns the configured quadratic violation penalty when present.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the configured linear violation penalty when present.
    pub fn linear_violation_penalty(&self) -> Option<LinearViolationPenalty> {
        self.linear_violation_penalty
    }

    /// Returns this soft side's objective contribution; hard sides return `None`.
    pub fn loss(&self) -> Option<f64> {
        self.loss
    }
}

impl PointToLevelSetRelationAssessment {
    /// Returns the caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        self.relation.source_id()
    }

    /// Returns the stable semantic role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the finite sampled location.
    pub fn location(&self) -> Point3 {
        self.relation.location()
    }

    /// Returns the referenced shared level set.
    pub fn group_id(&self) -> &GroupId {
        self.relation.group_id()
    }

    /// Returns the explicitly declared increasing or decreasing side.
    pub fn side(&self) -> PointToLevelSetSide {
        self.relation.side()
    }

    /// Returns the strictly positive field-value offset.
    pub fn minimum_offset(&self) -> MinimumFieldOffset {
        self.relation.minimum_offset()
    }

    /// Returns the recovered field value at the sampled point.
    pub fn recovered_point_value(&self) -> f64 {
        self.recovered_point_value
    }

    /// Returns the recovered shared-level field value.
    pub fn recovered_level_value(&self) -> f64 {
        self.recovered_level_value
    }

    /// Returns the oriented point-to-level difference in field-value units.
    pub fn recovered_field_offset(&self) -> f64 {
        self.recovered_field_offset
    }

    /// Returns nonnegative physical satisfaction slack.
    pub fn slack(&self) -> f64 {
        self.slack
    }

    /// Returns nonnegative physical violation.
    pub fn violation(&self) -> f64 {
        self.violation
    }

    /// Returns the physical acceptance tolerance.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns the stable physical active state.
    pub fn active_state(&self) -> BoundActiveState {
        self.active_state
    }

    /// Returns the configured quadratic violation penalty when present.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the configured linear violation penalty when present.
    pub fn linear_violation_penalty(&self) -> Option<LinearViolationPenalty> {
        self.linear_violation_penalty
    }

    /// Returns this soft relation's objective contribution; hard returns `None`.
    pub fn loss(&self) -> Option<f64> {
        self.loss
    }
}

#[derive(Debug, Clone, PartialEq)]
enum RecoveredSharedLevelSetRelationValues {
    StratigraphicAge { younger: f64, older: f64 },
    FieldLevelOrder { lower: f64, upper: f64 },
}

/// Physical recovery assessment for one relation between shared level sets.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedLevelSetRelationAssessment {
    relation: SharedLevelSetRelationInput,
    orientation: SharedLevelSetRelationOrientation,
    semantic_role: SemanticRolePath,
    recovered_values: RecoveredSharedLevelSetRelationValues,
    recovered_field_separation: f64,
    slack: f64,
    violation: f64,
    tolerance: f64,
    active_state: BoundActiveState,
    quadratic_penalty: Option<QuadraticPenalty>,
    linear_violation_penalty: Option<LinearViolationPenalty>,
    loss: Option<f64>,
}

impl SharedLevelSetRelationAssessment {
    /// Returns the caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        self.relation.source_id()
    }

    /// Returns the stable semantic role of this relation family.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the caller-declared relation kind.
    pub fn kind(&self) -> SharedLevelSetRelationKind {
        self.relation.kind()
    }

    /// Returns the geologically younger group for a stratigraphic age relation.
    pub fn younger_group_id(&self) -> Option<&GroupId> {
        self.relation.younger_group_id()
    }

    /// Returns the geologically older group for a stratigraphic age relation.
    pub fn older_group_id(&self) -> Option<&GroupId> {
        self.relation.older_group_id()
    }

    /// Returns the declared lower group for a Field Level Order.
    pub fn lower_group_id(&self) -> Option<&GroupId> {
        self.relation.lower_group_id()
    }

    /// Returns the declared upper group for a Field Level Order.
    pub fn upper_group_id(&self) -> Option<&GroupId> {
        self.relation.upper_group_id()
    }

    /// Returns the group on the lower-field-value side after semantic orientation.
    pub fn lower_field_group_id(&self) -> &GroupId {
        &self.orientation.lower_field_group_id
    }

    /// Returns the group on the upper-field-value side after semantic orientation.
    pub fn upper_field_group_id(&self) -> &GroupId {
        &self.orientation.upper_field_group_id
    }

    /// Returns the explicit age-to-field mapping for age relations.
    pub fn field_direction(&self) -> Option<StratigraphicFieldDirection> {
        self.orientation.field_direction
    }

    /// Returns the strict minimum for age relations; non-strict order returns `None`.
    pub fn minimum_separation(&self) -> Option<MinimumFieldSeparation> {
        self.orientation.minimum_separation
    }

    /// Returns the recovered younger-group value for a stratigraphic age relation.
    pub fn recovered_younger_value(&self) -> Option<f64> {
        match self.recovered_values {
            RecoveredSharedLevelSetRelationValues::StratigraphicAge { younger, .. } => {
                Some(younger)
            }
            RecoveredSharedLevelSetRelationValues::FieldLevelOrder { .. } => None,
        }
    }

    /// Returns the recovered older-group value for a stratigraphic age relation.
    pub fn recovered_older_value(&self) -> Option<f64> {
        match self.recovered_values {
            RecoveredSharedLevelSetRelationValues::StratigraphicAge { older, .. } => Some(older),
            RecoveredSharedLevelSetRelationValues::FieldLevelOrder { .. } => None,
        }
    }

    /// Returns the recovered lower-group value for a Field Level Order.
    pub fn recovered_lower_value(&self) -> Option<f64> {
        match self.recovered_values {
            RecoveredSharedLevelSetRelationValues::FieldLevelOrder { lower, .. } => Some(lower),
            RecoveredSharedLevelSetRelationValues::StratigraphicAge { .. } => None,
        }
    }

    /// Returns the recovered upper-group value for a Field Level Order.
    pub fn recovered_upper_value(&self) -> Option<f64> {
        match self.recovered_values {
            RecoveredSharedLevelSetRelationValues::FieldLevelOrder { upper, .. } => Some(upper),
            RecoveredSharedLevelSetRelationValues::StratigraphicAge { .. } => None,
        }
    }

    /// Returns upper-field value minus lower-field value in physical field units.
    pub fn recovered_field_separation(&self) -> f64 {
        self.recovered_field_separation
    }

    /// Returns the nonnegative satisfaction slack in field-value units.
    pub fn slack(&self) -> f64 {
        self.slack
    }

    /// Returns the nonnegative violation in field-value units.
    pub fn violation(&self) -> f64 {
        self.violation
    }

    /// Returns the physical acceptance tolerance.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns whether the relation is active at the accepted tolerance.
    pub fn active_state(&self) -> BoundActiveState {
        self.active_state
    }

    /// Returns the configured quadratic violation penalty when present.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the configured linear violation penalty when present.
    pub fn linear_violation_penalty(&self) -> Option<LinearViolationPenalty> {
        self.linear_violation_penalty
    }

    /// Returns this soft relation's objective contribution; hard relations return `None`.
    pub fn loss(&self) -> Option<f64> {
        self.loss
    }
}

impl DirectionalDerivativeIntervalAssessment {
    /// Returns the caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable lower- or upper-side semantic role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the oriented unit direction used by this derivative relation.
    pub fn direction(&self) -> Vector3 {
        self.direction
    }

    /// Returns which side of the derivative interval this assessment describes.
    pub fn side(&self) -> BoundSide {
        self.side
    }

    /// Returns the finite bound in field-value-per-length units.
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Returns the independently recovered directional derivative.
    pub fn recovered_directional_derivative(&self) -> f64 {
        self.recovered_directional_derivative
    }

    /// Returns nonnegative physical satisfaction slack.
    pub fn slack(&self) -> f64 {
        self.slack
    }

    /// Returns nonnegative physical violation.
    pub fn violation(&self) -> f64 {
        self.violation
    }

    /// Returns the physical acceptance tolerance for this side.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns the stable active/inactive state derived in physical units.
    pub fn active_state(&self) -> BoundActiveState {
        self.active_state
    }

    /// Returns the configured quadratic violation penalty when present.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the configured linear violation penalty when present.
    pub fn linear_violation_penalty(&self) -> Option<LinearViolationPenalty> {
        self.linear_violation_penalty
    }

    /// Returns this soft side's objective contribution; hard sides return `None`.
    pub fn loss(&self) -> Option<f64> {
        self.loss
    }
}

impl FieldValueBoundAssessment {
    /// Returns the caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable lower- or upper-side semantic role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns which side of the bound this assessment describes.
    pub fn side(&self) -> BoundSide {
        self.side
    }

    /// Returns the finite bound in field-value units.
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Returns the independently recovered field value at the bound support.
    pub fn recovered_value(&self) -> f64 {
        self.recovered_value
    }

    /// Returns nonnegative physical satisfaction slack in field-value units.
    pub fn slack(&self) -> f64 {
        self.slack
    }

    /// Returns nonnegative physical violation in field-value units.
    pub fn violation(&self) -> f64 {
        self.violation
    }

    /// Returns the physical acceptance tolerance for this side.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns the stable active/inactive state derived in physical units.
    pub fn active_state(&self) -> BoundActiveState {
        self.active_state
    }

    /// Returns the configured quadratic violation penalty when present.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the configured linear violation penalty when present.
    pub fn linear_violation_penalty(&self) -> Option<LinearViolationPenalty> {
        self.linear_violation_penalty
    }

    /// Returns this soft side's objective contribution; hard sides return `None`.
    pub fn loss(&self) -> Option<f64> {
        self.loss
    }
}

/// Recovered physical residual and objective contribution for one soft field
/// value.
#[derive(Debug, Clone, PartialEq)]
pub struct SoftFieldValueAssessment {
    source_id: SourceId,
    semantic_role: SemanticRolePath,
    target: f64,
    recovered_value: f64,
    residual: f64,
    quadratic_penalty: Option<QuadraticPenalty>,
    standard_deviation: Option<StandardDeviation>,
    loss: f64,
}

impl SoftFieldValueAssessment {
    /// Returns the caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable field-value residual role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the observed field value in the caller's units.
    pub fn target(&self) -> f64 {
        self.target
    }

    /// Returns the independently recovered field value in caller units.
    pub fn recovered_value(&self) -> f64 {
        self.recovered_value
    }

    /// Returns recovered value minus target in field-value units.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Returns the non-statistical penalty when that path configured this
    /// relation.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the statistical standard deviation when that path configured
    /// this relation.
    pub fn standard_deviation(&self) -> Option<StandardDeviation> {
        self.standard_deviation
    }

    /// Returns this independent residual's dimensionless objective
    /// contribution.
    pub fn loss(&self) -> f64 {
        self.loss
    }
}

/// Recovered complete-gradient vector residual and its one block-level loss.
#[derive(Debug, Clone, PartialEq)]
pub struct SoftGradientAssessment {
    source_id: SourceId,
    semantic_role: SemanticRolePath,
    target: Vector3,
    recovered_gradient: Vector3,
    residual: Vector3,
    quadratic_penalty: Option<QuadraticPenalty>,
    standard_deviation: Option<StandardDeviation>,
    covariance: Option<CovarianceMatrix>,
    whitened_residual: Box<[f64]>,
    whitening_round_trip_error: f64,
    loss: f64,
}

impl SoftGradientAssessment {
    /// Returns the caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable ordered-vector residual role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the observed complete gradient in caller units.
    pub fn target(&self) -> Vector3 {
        self.target
    }

    /// Returns the independently recovered complete gradient.
    pub fn recovered_gradient(&self) -> Vector3 {
        self.recovered_gradient
    }

    /// Returns recovered gradient minus target in field-per-length units.
    pub fn residual(&self) -> Vector3 {
        self.residual
    }

    /// Returns the Euclidean vector penalty when configured.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the isotropic statistical standard deviation when configured.
    pub fn standard_deviation(&self) -> Option<StandardDeviation> {
        self.standard_deviation
    }

    /// Returns the statistical covariance when configured.
    pub fn covariance(&self) -> Option<&CovarianceMatrix> {
        self.covariance.as_ref()
    }

    /// Returns the derived whitened residual in canonical component order.
    pub fn whitened_residual(&self) -> &[f64] {
        &self.whitened_residual
    }

    /// Returns the whitening forward/inverse recovery error.
    pub fn whitening_round_trip_error(&self) -> f64 {
        self.whitening_round_trip_error
    }

    /// Returns this vector block's dimensionless objective contribution.
    pub fn loss(&self) -> f64 {
        self.loss
    }
}

/// Recovered physical observables for one Directed Normal semantic relation.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectedNormalAssessment {
    source_id: SourceId,
    direction_semantic_role: SemanticRolePath,
    slope_semantic_role: SemanticRolePath,
    direction: Vector3,
    recovered_gradient: Vector3,
    projection_residual: Vector3,
    projection_residual_norm: f64,
    recovered_slope: f64,
    minimum_slope: MinimumNormalSlope,
    slope_slack: f64,
    slope_violation: f64,
    direction_tolerance: Option<f64>,
    slope_tolerance: f64,
    slope_active_state: BoundActiveState,
    direction_enforcement: NormalDirectionEnforcement,
    minimum_slope_enforcement: MinimumNormalSlopeEnforcement,
    direction_loss: Option<f64>,
    slope_loss: Option<f64>,
    input_axis: Option<Vector3>,
    polarity_resolution_source_id: Option<SourceId>,
    polarity_selection: Option<PolaritySelection>,
}

impl DirectedNormalAssessment {
    /// Returns the caller identity of the Directed or resolved Axial input.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable rotation-invariant direction-channel role.
    pub fn direction_semantic_role(&self) -> &SemanticRolePath {
        &self.direction_semantic_role
    }

    /// Returns the stable independent minimum-slope channel role.
    pub fn slope_semantic_role(&self) -> &SemanticRolePath {
        &self.slope_semantic_role
    }

    /// Returns the resolved oriented physical unit normal.
    pub fn direction(&self) -> Vector3 {
        self.direction
    }

    /// Returns the complete recovered gradient in physical input units.
    pub fn recovered_gradient(&self) -> Vector3 {
        self.recovered_gradient
    }

    /// Returns `(I - n n^T) grad(f)` in physical derivative units.
    pub fn projection_residual(&self) -> Vector3 {
        self.projection_residual
    }

    /// Returns the Euclidean norm of the rotation-invariant projection residual.
    pub fn projection_residual_norm(&self) -> f64 {
        self.projection_residual_norm
    }

    /// Returns `n^T grad(f)` in physical derivative units.
    pub fn recovered_slope(&self) -> f64 {
        self.recovered_slope
    }

    /// Returns the caller's finite positive minimum normal slope.
    pub fn minimum_slope(&self) -> MinimumNormalSlope {
        self.minimum_slope
    }

    /// Returns nonnegative physical satisfaction slack for the slope lower bound.
    pub fn slope_slack(&self) -> f64 {
        self.slope_slack
    }

    /// Returns nonnegative physical violation of the minimum slope.
    pub fn slope_violation(&self) -> f64 {
        self.slope_violation
    }

    /// Returns the hard projection tolerance, or `None` for a soft direction channel.
    pub fn direction_tolerance(&self) -> Option<f64> {
        self.direction_tolerance
    }

    /// Returns the physical acceptance tolerance for the slope relation.
    pub fn slope_tolerance(&self) -> f64 {
        self.slope_tolerance
    }

    /// Returns the stable active state of the minimum-slope lower bound.
    pub fn slope_active_state(&self) -> BoundActiveState {
        self.slope_active_state
    }

    /// Returns the direction channel's Euclidean quadratic penalty when present.
    pub fn direction_quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.direction_enforcement.quadratic_penalty()
    }

    /// Returns the direction channel's isotropic statistical scale when present.
    pub fn direction_standard_deviation(&self) -> Option<StandardDeviation> {
        self.direction_enforcement.standard_deviation()
    }

    /// Returns the slope channel's quadratic violation penalty when present.
    pub fn slope_quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.minimum_slope_enforcement.quadratic_penalty()
    }

    /// Returns the slope channel's linear violation penalty when present.
    pub fn slope_linear_violation_penalty(&self) -> Option<LinearViolationPenalty> {
        self.minimum_slope_enforcement.linear_violation_penalty()
    }

    /// Returns the direction channel's independent objective contribution.
    pub fn direction_loss(&self) -> Option<f64> {
        self.direction_loss
    }

    /// Returns the slope channel's independent objective contribution.
    pub fn slope_loss(&self) -> Option<f64> {
        self.slope_loss
    }

    /// Returns the retained normalized Axial input orientation, when applicable.
    pub fn input_axis(&self) -> Option<Vector3> {
        self.input_axis
    }

    /// Returns independent Polarity Resolution provenance for a resolved Axial input.
    pub fn polarity_resolution_source_id(&self) -> Option<&SourceId> {
        self.polarity_resolution_source_id.as_ref()
    }

    /// Returns the explicit Axial polarity selection, when applicable.
    pub fn polarity_selection(&self) -> Option<PolaritySelection> {
        self.polarity_selection
    }
}

/// Recovered scalar directional-derivative residual for one soft Tangent.
#[derive(Debug, Clone, PartialEq)]
pub struct SoftTangentAssessment {
    source_id: SourceId,
    semantic_role: SemanticRolePath,
    recovered_directional_derivative: f64,
    residual: f64,
    quadratic_penalty: Option<QuadraticPenalty>,
    standard_deviation: Option<StandardDeviation>,
    loss: f64,
}

impl SoftTangentAssessment {
    /// Returns the caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable directional-derivative residual role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the Tangent target, which is exactly zero.
    pub fn target(&self) -> f64 {
        0.0
    }

    /// Returns the independently recovered directional derivative.
    pub fn recovered_directional_derivative(&self) -> f64 {
        self.recovered_directional_derivative
    }

    /// Returns recovered directional derivative minus zero target.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Returns the non-statistical penalty when configured.
    pub fn quadratic_penalty(&self) -> Option<QuadraticPenalty> {
        self.quadratic_penalty
    }

    /// Returns the statistical standard deviation when configured.
    pub fn standard_deviation(&self) -> Option<StandardDeviation> {
        self.standard_deviation
    }

    /// Returns this scalar residual's dimensionless objective contribution.
    pub fn loss(&self) -> f64 {
        self.loss
    }
}

/// One original-unit member residual recovered from a covariance group.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceGroupMemberAssessment {
    source_id: SourceId,
    semantic_role: SemanticRolePath,
    dimension: ResidualDimension,
    target_components: Box<[f64]>,
    recovered_components: Box<[f64]>,
    residual_components: Box<[f64]>,
}

impl CovarianceGroupMemberAssessment {
    /// Returns the stable caller-owned member identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the stable scalar or ordered-vector member role.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the common physical residual dimension.
    pub fn dimension(&self) -> ResidualDimension {
        self.dimension
    }

    /// Returns target components in the group's explicit component order.
    pub fn target_components(&self) -> &[f64] {
        &self.target_components
    }

    /// Returns independently recovered components in original physical units.
    pub fn recovered_components(&self) -> &[f64] {
        &self.recovered_components
    }

    /// Returns recovered-minus-target components in original physical units.
    pub fn residual_components(&self) -> &[f64] {
        &self.residual_components
    }
}

/// Recovered observables and the only identifiable objective contribution for
/// one named covariance group.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceGroupAssessment {
    group_id: GroupId,
    covariance: CovarianceMatrix,
    members: Vec<CovarianceGroupMemberAssessment>,
    whitened_residual: Box<[f64]>,
    whitening_round_trip_error: f64,
    objective_contribution: f64,
}

impl CovarianceGroupAssessment {
    /// Returns the stable caller-owned group identity.
    pub fn group_id(&self) -> &GroupId {
        &self.group_id
    }

    /// Returns the checked covariance in explicit flattened member order.
    pub fn covariance(&self) -> &CovarianceMatrix {
        &self.covariance
    }

    /// Returns original member blocks without invented member-level losses.
    pub fn members(&self) -> &[CovarianceGroupMemberAssessment] {
        &self.members
    }

    /// Returns the derived whitened residual vector.
    pub fn whitened_residual(&self) -> &[f64] {
        &self.whitened_residual
    }

    /// Returns the whitening forward/inverse recovery error.
    pub fn whitening_round_trip_error(&self) -> f64 {
        self.whitening_round_trip_error
    }

    /// Returns the unique group-level dimensionless objective contribution.
    pub fn objective_contribution(&self) -> f64 {
        self.objective_contribution
    }
}

impl HardRelationAssessment {
    /// Returns the caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the referenced semantic group, when this relation uses one.
    pub fn group_id(&self) -> Option<&GroupId> {
        self.group_id.as_ref()
    }

    /// Returns the stable semantic component path.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns this residual's physical dimension.
    pub fn dimension(&self) -> ResidualDimension {
        self.dimension
    }

    /// Returns the hard target in physical units.
    pub fn target(&self) -> f64 {
        self.target
    }

    /// Returns the independently recovered model value in physical units.
    pub fn recovered_value(&self) -> f64 {
        self.recovered_value
    }

    /// Returns recovered value minus target, in physical units.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Returns the versioned physical acceptance tolerance.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Returns the gauge-invariant characteristic scale used by policy.
    pub fn characteristic_scale(&self) -> f64 {
        self.characteristic_scale
    }

    /// Returns the relation-local physical reference scale.
    pub fn relation_reference_scale(&self) -> f64 {
        self.relation_reference_scale
    }

    /// Returns the tolerance in normalized representation coordinates.
    pub fn standard_tolerance(&self) -> f64 {
        self.standard_tolerance
    }

    /// Returns the tolerance after backend KKT scaling, when this relation
    /// participated in the solver system.
    pub fn scaled_kkt_tolerance(&self) -> Option<f64> {
        self.scaled_kkt_tolerance
    }

    /// Returns the independently recovered physical tolerance.
    pub fn recovered_physical_tolerance(&self) -> f64 {
        self.recovered_physical_tolerance
    }

    /// Returns the tolerance forward/inverse round-trip error.
    pub fn tolerance_round_trip_error(&self) -> f64 {
        self.tolerance_round_trip_error
    }
}

/// One recovered shared-level semantic latent in physical field units.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedLevelValue {
    group_id: GroupId,
    value: f64,
    field_unit: FieldUnitLabel,
    member_source_ids: Box<[SourceId]>,
}

impl SharedLevelValue {
    /// Returns the stable identity of the recovered semantic latent.
    pub fn group_id(&self) -> &GroupId {
        &self.group_id
    }

    /// Returns the recovered shared field value.
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Returns the problem's caller-declared field-unit label.
    pub fn field_unit(&self) -> &FieldUnitLabel {
        &self.field_unit
    }

    /// Returns complete member provenance in stable SourceId order.
    pub fn member_source_ids(&self) -> &[SourceId] {
        &self.member_source_ids
    }
}

/// Typed audit report shared by successful and failed fits.
#[derive(Debug, Clone)]
pub struct FitReport {
    problem_size: ProblemSize,
    resolved_kernel: KernelConfig,
    field_energy_normalization: FieldEnergyNormalization,
    numerical_policy: NumericalPolicyId,
    requested_thread_budget: ThreadBudget,
    hard_relations: Vec<HardRelationAssessment>,
    field_value_bounds: Vec<FieldValueBoundAssessment>,
    directional_derivative_intervals: Vec<DirectionalDerivativeIntervalAssessment>,
    field_separation_intervals: Vec<FieldSeparationIntervalAssessment>,
    point_to_level_set_relations: Vec<PointToLevelSetRelationAssessment>,
    shared_level_set_relations: Vec<SharedLevelSetRelationAssessment>,
    soft_field_values: Vec<SoftFieldValueAssessment>,
    soft_gradients: Vec<SoftGradientAssessment>,
    directed_normals: Vec<DirectedNormalAssessment>,
    soft_tangents: Vec<SoftTangentAssessment>,
    covariance_groups: Vec<CovarianceGroupAssessment>,
    shared_level_values: Vec<SharedLevelValue>,
    field_energy: Option<f64>,
    total_objective: Option<f64>,
    backend_fingerprint: Option<BackendFingerprint>,
    attempts: Vec<SolveAttemptRecord>,
    recovery_verification: Option<RecoveryVerificationEvidence>,
    direct_input_conflicts: Vec<DirectInputConflictEvidence>,
    relation_graph_conflicts: Vec<RelationGraphConflictEvidence>,
    shared_level_set_relation_conflicts: Vec<SharedLevelSetRelationConflictEvidence>,
    execution_failure: Option<AttemptFailureEvidence>,
    cubic_analysis: Option<CubicAnalysisEvidence>,
    backend_rank: Option<RankEvidence>,
    interpretable_rank_deficiency: Option<InterpretableRankDeficiencyEvidence>,
    inertia: Option<InertiaEvidence>,
    canonical_acceptance: Option<CanonicalAcceptanceEvidence>,
    capacity: Option<CapacityEvidence>,
    analysis_failure: Option<AnalysisFailureEvidence>,
    infeasibility_certificate: Option<InfeasibilityCertificateEvidence>,
    recession_ray: Option<RecessionRayEvidence>,
    unidentified_additive_gauge: Option<UnidentifiedAdditiveGaugeEvidence>,
    uninformative_shared_level_sets: Vec<UninformativeSharedLevelSetEvidence>,
    unresolved_axial_normals: Vec<UnresolvedAxialNormalEvidence>,
}

impl FitReport {
    /// Returns the canonical problem size used for this fit.
    pub fn problem_size(&self) -> ProblemSize {
        self.problem_size
    }

    /// Returns the resolved kernel configuration.
    pub fn resolved_kernel(&self) -> &KernelConfig {
        &self.resolved_kernel
    }

    /// Returns the resolved physical FieldEnergy normalization.
    pub fn field_energy_normalization(&self) -> FieldEnergyNormalization {
        self.field_energy_normalization
    }

    /// Returns the versioned numerical policy identity.
    pub fn numerical_policy(&self) -> NumericalPolicyId {
        self.numerical_policy
    }

    /// Returns the caller's resource request stored in the snapshot.
    pub fn requested_thread_budget(&self) -> ThreadBudget {
        self.requested_thread_budget
    }

    /// Returns physical hard-relation recovery assessments in stable order.
    pub fn hard_relations(&self) -> &[HardRelationAssessment] {
        &self.hard_relations
    }

    /// Returns lower/upper Field Value Bound assessments in SourceId/role order.
    pub fn field_value_bounds(&self) -> &[FieldValueBoundAssessment] {
        &self.field_value_bounds
    }

    /// Returns Directional Derivative Interval sides in SourceId/role order.
    pub fn directional_derivative_intervals(&self) -> &[DirectionalDerivativeIntervalAssessment] {
        &self.directional_derivative_intervals
    }

    /// Returns Field Separation Interval sides in SourceId/role order.
    pub fn field_separation_intervals(&self) -> &[FieldSeparationIntervalAssessment] {
        &self.field_separation_intervals
    }

    /// Returns Point to Level Set Relation assessments in SourceId order.
    pub fn point_to_level_set_relations(&self) -> &[PointToLevelSetRelationAssessment] {
        &self.point_to_level_set_relations
    }

    /// Returns shared-level-set relation assessments in stable SourceId order.
    pub fn shared_level_set_relations(&self) -> &[SharedLevelSetRelationAssessment] {
        &self.shared_level_set_relations
    }

    /// Returns soft Field Value assessments in stable SourceId order.
    pub fn soft_field_values(&self) -> &[SoftFieldValueAssessment] {
        &self.soft_field_values
    }

    /// Returns complete soft Gradient assessments in stable SourceId order.
    pub fn soft_gradients(&self) -> &[SoftGradientAssessment] {
        &self.soft_gradients
    }

    /// Returns Directed and explicitly resolved Axial Normal assessments.
    pub fn directed_normals(&self) -> &[DirectedNormalAssessment] {
        &self.directed_normals
    }

    /// Returns soft Tangent assessments in stable SourceId order.
    pub fn soft_tangents(&self) -> &[SoftTangentAssessment] {
        &self.soft_tangents
    }

    /// Returns named covariance-group assessments in stable GroupId order.
    pub fn covariance_groups(&self) -> &[CovarianceGroupAssessment] {
        &self.covariance_groups
    }

    /// Returns every recovered shared-level semantic latent in GroupId order.
    pub fn shared_level_values(&self) -> &[SharedLevelValue] {
        &self.shared_level_values
    }

    /// Returns accepted FieldEnergy, or `None` when no model was accepted.
    pub fn field_energy(&self) -> Option<f64> {
        self.field_energy
    }

    /// Returns the accepted total objective, or `None` on failure.
    pub fn total_objective(&self) -> Option<f64> {
        self.total_objective
    }

    /// Returns the selected backend fingerprint when a candidate was accepted.
    pub fn backend_fingerprint(&self) -> Option<&BackendFingerprint> {
        self.backend_fingerprint.as_ref()
    }

    /// Returns bounded attempt records in execution order.
    pub fn attempts(&self) -> &[SolveAttemptRecord] {
        &self.attempts
    }

    /// Returns physical rejection evidence when Recover and Verify failed.
    pub fn recovery_verification(&self) -> Option<&RecoveryVerificationEvidence> {
        self.recovery_verification.as_ref()
    }

    /// Returns stable source evidence for a direct hard-input conflict.
    pub fn direct_input_conflict(&self) -> Option<&DirectInputConflictEvidence> {
        self.direct_input_conflicts.first()
    }

    /// Returns every direct hard-input conflict in stable semantic/source order.
    pub fn direct_input_conflicts(&self) -> &[DirectInputConflictEvidence] {
        &self.direct_input_conflicts
    }

    /// Returns complete path provenance for a hard relation-graph conflict.
    pub fn relation_graph_conflict(&self) -> Option<&RelationGraphConflictEvidence> {
        self.relation_graph_conflicts.first()
    }

    /// Returns every hard relation-graph conflict in stable semantic/source order.
    pub fn relation_graph_conflicts(&self) -> &[RelationGraphConflictEvidence] {
        &self.relation_graph_conflicts
    }

    /// Returns the first impossible hard shared-level-set relation proof.
    pub fn shared_level_set_relation_conflict(
        &self,
    ) -> Option<&SharedLevelSetRelationConflictEvidence> {
        self.shared_level_set_relation_conflicts.first()
    }

    /// Returns every impossible hard shared-level-set relation proof in stable order.
    pub fn shared_level_set_relation_conflicts(&self) -> &[SharedLevelSetRelationConflictEvidence] {
        &self.shared_level_set_relation_conflicts
    }

    /// Returns the terminal backend-contract or numerical failure evidence.
    pub fn execution_failure(&self) -> Option<AttemptFailureEvidence> {
        self.execution_failure
    }

    /// Returns complete Cubic representation analysis when construction succeeded.
    pub fn cubic_analysis(&self) -> Option<&CubicAnalysisEvidence> {
        self.cubic_analysis.as_ref()
    }

    /// Returns representation or backend rank evidence at the terminal boundary.
    pub fn rank_evidence(&self) -> Option<&RankEvidence> {
        self.backend_rank.as_ref()
    }

    /// Returns backend KKT rank evidence when that analysis was reached.
    pub fn backend_rank(&self) -> Option<&RankEvidence> {
        self.backend_rank
            .as_ref()
            .filter(|evidence| evidence.domain() == RankEvidenceDomain::BackendKkt)
    }

    /// Returns canonical mode evidence when a rank loss has semantic meaning.
    pub fn interpretable_rank_deficiency(&self) -> Option<&InterpretableRankDeficiencyEvidence> {
        self.interpretable_rank_deficiency.as_ref()
    }

    /// Returns expected and observed Equality KKT inertia when available.
    pub fn inertia(&self) -> Option<InertiaEvidence> {
        self.inertia
    }

    /// Returns physical Recover-and-Verify acceptance evidence when reached.
    pub fn canonical_acceptance(&self) -> Option<&CanonicalAcceptanceEvidence> {
        self.canonical_acceptance.as_ref()
    }

    /// Returns checked capacity evidence when planning rejected the fit.
    pub fn capacity(&self) -> Option<CapacityEvidence> {
        self.capacity
    }

    /// Returns structured numerical-analysis failure evidence when available.
    pub fn analysis_failure(&self) -> Option<&AnalysisFailureEvidence> {
        self.analysis_failure.as_ref()
    }

    /// Returns independently validated Farkas-ray evidence for infeasibility.
    pub fn infeasibility_certificate(&self) -> Option<InfeasibilityCertificateEvidence> {
        self.infeasibility_certificate.clone()
    }

    /// Returns independently validated recession-ray evidence for unboundedness.
    pub fn recession_ray(&self) -> Option<RecessionRayEvidence> {
        self.recession_ray.clone()
    }

    /// Returns structural evidence for a missing additive-field representative.
    pub fn unidentified_additive_gauge(&self) -> Option<&UnidentifiedAdditiveGaugeEvidence> {
        self.unidentified_additive_gauge.as_ref()
    }

    /// Returns structural evidence for a disconnected one-member shared level set.
    pub fn uninformative_shared_level_set(&self) -> Option<&UninformativeSharedLevelSetEvidence> {
        self.uninformative_shared_level_sets.first()
    }

    /// Returns every disconnected one-member shared level set in stable group order.
    pub fn uninformative_shared_level_sets(&self) -> &[UninformativeSharedLevelSetEvidence] {
        &self.uninformative_shared_level_sets
    }

    /// Returns every unresolved Axial Normal in stable SourceId order.
    pub fn unresolved_axial_normals(&self) -> &[UnresolvedAxialNormalEvidence] {
        &self.unresolved_axial_normals
    }
}

pub(crate) fn fit_snapshot(snapshot: &ProblemSnapshot) -> Result<FitSuccess, FitFailure> {
    let scalar_relation_counts = scalar_relation_counts(snapshot).unwrap_or(ScalarRelationCounts {
        hard: usize::MAX,
        soft: usize::MAX,
    });
    let scalar_relation_count = scalar_relation_counts.total().unwrap_or(usize::MAX);
    let quadratic_objective_terms = quadratic_objective_term_count(snapshot).unwrap_or(usize::MAX);
    let linear_objective_terms = linear_objective_term_count(snapshot).unwrap_or(usize::MAX);
    let input_observations = snapshot
        .inner
        .observations
        .len()
        .checked_add(snapshot.inner.directed_normals.len())
        .and_then(|count| count.checked_add(snapshot.inner.axial_normals.len()))
        .and_then(|count| count.checked_add(snapshot.inner.polarity_resolutions.len()))
        .and_then(|count| count.checked_add(snapshot.inner.field_value_bounds.len()))
        .and_then(|count| count.checked_add(snapshot.inner.directional_derivative_intervals.len()))
        .and_then(|count| count.checked_add(snapshot.inner.field_separation_intervals.len()))
        .and_then(|count| count.checked_add(snapshot.inner.point_to_level_set_relations.len()))
        .and_then(|count| count.checked_add(snapshot.inner.shared_level_set_relations.len()))
        .unwrap_or(usize::MAX);
    let source_identifier_bytes = source_identifier_bytes(snapshot).unwrap_or(usize::MAX);
    let conservative_problem_size = conservative_problem_size(
        input_observations,
        scalar_relation_counts,
        snapshot.inner.shared_level_sets.len(),
        quadratic_objective_terms,
        linear_objective_terms,
        snapshot
            .inner
            .field_value_bounds
            .iter()
            .map(|bound| {
                usize::from(bound.lower().is_some()) + usize::from(bound.upper().is_some())
            })
            .chain(
                snapshot
                    .inner
                    .directional_derivative_intervals
                    .iter()
                    .map(|interval| {
                        usize::from(interval.lower().is_some())
                            + usize::from(interval.upper().is_some())
                    }),
            )
            .chain(std::iter::repeat_n(
                2_usize,
                snapshot.inner.field_separation_intervals.len(),
            ))
            .chain(std::iter::repeat_n(
                1_usize,
                snapshot.inner.point_to_level_set_relations.len(),
            ))
            .chain(std::iter::repeat_n(
                1_usize,
                snapshot.inner.shared_level_set_relations.len(),
            ))
            .chain(std::iter::repeat_n(
                1_usize,
                snapshot.inner.directed_normals.len() + snapshot.inner.axial_normals.len(),
            ))
            .sum(),
    );
    let mut preflight_report = empty_report(snapshot, conservative_problem_size);
    let resolved_axial_sources = snapshot
        .inner
        .polarity_resolutions
        .iter()
        .map(|resolution| resolution.axial_normal_source_id())
        .collect::<std::collections::BTreeSet<_>>();
    preflight_report.unresolved_axial_normals = snapshot
        .inner
        .axial_normals
        .iter()
        .filter(|normal| !resolved_axial_sources.contains(normal.source_id()))
        .map(|normal| {
            UnresolvedAxialNormalEvidence::new(normal.source_id().clone(), normal.input_axis())
        })
        .collect();
    let source_lifecycle_capacity =
        plan_source_lifecycle_capacity(scalar_relation_count, source_identifier_bytes);
    let source_preflight_capacity =
        plan_source_preflight_capacity(scalar_relation_count, source_identifier_bytes);
    preflight_report.uninformative_shared_level_sets = snapshot
        .inner
        .shared_level_sets
        .iter()
        .filter(|group| {
            group.members().len() == 1 && !snapshot_references_group(snapshot, group.group_id())
        })
        .map(|group| {
            UninformativeSharedLevelSetEvidence::new(
                group.group_id().clone(),
                group.members()[0].source_id().clone(),
                false,
            )
        })
        .collect();
    preflight_report
        .uninformative_shared_level_sets
        .sort_by(|left, right| {
            left.group_id()
                .cmp(right.group_id())
                .then_with(|| left.member_source_id().cmp(right.member_source_id()))
        });
    let has_absolute_reference = snapshot
        .inner
        .observations
        .iter()
        .any(|observation| matches!(observation, ObservationInput::FieldValue(_)))
        || snapshot.inner.covariance_groups.iter().any(|group| {
            group
                .members()
                .iter()
                .any(|member| matches!(member, CovarianceGroupMember::FieldValue(_)))
        })
        || !snapshot.inner.additive_field_gauges.is_empty()
        || !snapshot.inner.field_value_bounds.is_empty();
    if !has_absolute_reference {
        let mut source_ids = snapshot
            .inner
            .observations
            .iter()
            .map(|observation| observation.source_id().clone())
            .chain(
                snapshot
                    .inner
                    .covariance_groups
                    .iter()
                    .flat_map(|group| group.members())
                    .map(|member| member.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .shared_level_sets
                    .iter()
                    .flat_map(|group| group.members())
                    .map(|member| member.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .directional_derivative_intervals
                    .iter()
                    .map(|interval| interval.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .field_separation_intervals
                    .iter()
                    .map(|interval| interval.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .point_to_level_set_relations
                    .iter()
                    .map(|relation| relation.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .shared_level_set_relations
                    .iter()
                    .map(|relation| relation.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .directed_normals
                    .iter()
                    .map(|normal| normal.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .axial_normals
                    .iter()
                    .map(|normal| normal.source_id().clone()),
            )
            .chain(
                snapshot
                    .inner
                    .polarity_resolutions
                    .iter()
                    .map(|resolution| resolution.source_id().clone()),
            )
            .collect::<Vec<_>>();
        source_ids.sort();
        let group_ids = snapshot
            .inner
            .shared_level_sets
            .iter()
            .map(|group| group.group_id().clone())
            .collect();
        preflight_report.unidentified_additive_gauge = Some(
            UnidentifiedAdditiveGaugeEvidence::new(source_ids, group_ids, false),
        );
    }
    preflight_report.shared_level_set_relation_conflicts =
        preflight_shared_level_set_relation_conflicts(snapshot);
    if !preflight_report
        .shared_level_set_relation_conflicts
        .is_empty()
    {
        return Err(FitFailure {
            diagnosis: primary_preflight_diagnosis(&preflight_report)
                .expect("a shared-level-set conflict always supplies a diagnosis"),
            report: Box::new(preflight_report),
        });
    }
    if source_lifecycle_capacity.is_err() && source_preflight_capacity.is_err() {
        let evidence =
            source_lifecycle_capacity.expect_err("the source lifecycle failure was checked above");
        preflight_report.capacity = Some(public_capacity(&evidence));
        return Err(FitFailure {
            diagnosis: primary_preflight_diagnosis(&preflight_report)
                .expect("source lifecycle capacity always supplies a diagnosis"),
            report: Box::new(preflight_report),
        });
    }
    let lowering = lower_snapshot(snapshot);
    if let Some(conflict) = canonical_affine_value_cycle_conflict(snapshot, &lowering) {
        preflight_report
            .shared_level_set_relation_conflicts
            .push(conflict);
        preflight_report
            .shared_level_set_relation_conflicts
            .sort_by(|left, right| {
                left.source_ids()
                    .cmp(right.source_ids())
                    .then_with(|| left.group_ids().cmp(right.group_ids()))
            });
        return Err(FitFailure {
            diagnosis: primary_preflight_diagnosis(&preflight_report)
                .expect("a canonical affine cycle conflict always supplies a diagnosis"),
            report: Box::new(preflight_report),
        });
    }
    let fitting_uses = canonical_fitting_uses(
        &lowering.canonical_equalities,
        &lowering.canonical_soft_equalities,
        &lowering.canonical_affine_inequalities,
    );
    let problem_size_parts = CubicProblemSizeParts {
        input_observations,
        scalar_hard_relations: scalar_relation_counts.hard,
        scalar_soft_relations: scalar_relation_counts.soft,
        canonical_hard_equalities: lowering.canonical_equalities.len(),
        canonical_soft_equalities: lowering.canonical_soft_equalities.len(),
        quadratic_objective_terms,
        linear_objective_terms,
        affine_inequality_constraints: lowering.canonical_affine_inequalities.len()
            + lowering
                .canonical_affine_inequalities
                .iter()
                .filter(|relation| relation.violation_channel().is_some())
                .count(),
        center_coefficients: fitting_uses.len(),
        semantic_latents: lowering.semantic_latents.len(),
        solver_hard_equalities: lowering.solver_equality_count(),
        auxiliary_variables: lowering
            .canonical_affine_inequalities
            .iter()
            .filter(|relation| relation.violation_channel().is_some())
            .count(),
    };
    let exact_problem_size = if lowering.canonical_affine_inequalities.is_empty() {
        ProblemSize::cubic_equality(problem_size_parts)
    } else {
        ProblemSize::cubic_qp(problem_size_parts)
    };
    if let Some(problem_size) = exact_problem_size {
        preflight_report.problem_size = problem_size;
    }
    preflight_report.direct_input_conflicts = lowering.direct_input_conflicts.clone();
    preflight_report.relation_graph_conflicts = lowering.relation_graph_conflicts.clone();
    let capacity_failure = match source_lifecycle_capacity {
        Ok(()) if lowering.canonical_affine_inequalities.is_empty() => plan_snapshot_capacity(
            &lowering,
            fitting_uses.len(),
            scalar_relation_count,
            source_identifier_bytes,
        )
        .err(),
        Ok(()) => None,
        Err(evidence) => Some(evidence),
    };
    if let Some(evidence) = capacity_failure {
        preflight_report.capacity = Some(public_capacity(&evidence));
        if let Some(failure) = preflight_polynomial_analysis_failure(&fitting_uses) {
            retain_representation_failure(
                &mut preflight_report,
                &failure,
                &lowering.source_relations,
            );
        }
    }
    if let Some(diagnosis) = primary_preflight_diagnosis(&preflight_report) {
        return Err(FitFailure {
            diagnosis,
            report: Box::new(preflight_report),
        });
    }
    let problem_size = exact_problem_size
        .expect("the successful conservative capacity plan proves exact dimensions representable");
    preflight_report.problem_size = problem_size;
    let base_report = || preflight_report.clone();
    fit_snapshot_after_preflight(snapshot, lowering, problem_size, base_report)
}

fn primary_preflight_diagnosis(report: &FitReport) -> Option<ProblemDiagnosis> {
    if !report.unresolved_axial_normals.is_empty() {
        return Some(ProblemDiagnosis::UnresolvedSemantics);
    }
    if !report.uninformative_shared_level_sets.is_empty() {
        return Some(ProblemDiagnosis::UninformativeSharedLevelSet);
    }
    if !report.direct_input_conflicts.is_empty()
        || !report.relation_graph_conflicts.is_empty()
        || !report.shared_level_set_relation_conflicts.is_empty()
    {
        return Some(ProblemDiagnosis::DirectInputConflict);
    }
    if report.unidentified_additive_gauge.is_some() {
        return Some(ProblemDiagnosis::UnidentifiedAdditiveGauge);
    }
    if report.interpretable_rank_deficiency.is_some() {
        return Some(ProblemDiagnosis::UnidentifiedFieldMode);
    }
    report
        .capacity
        .is_some()
        .then_some(ProblemDiagnosis::CapacityExceeded)
}

fn empty_report(snapshot: &ProblemSnapshot, problem_size: ProblemSize) -> FitReport {
    FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: snapshot.inner.fit_configuration.numerical_policy(),
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations: Vec::new(),
        field_value_bounds: Vec::new(),
        directional_derivative_intervals: Vec::new(),
        field_separation_intervals: Vec::new(),
        point_to_level_set_relations: Vec::new(),
        shared_level_set_relations: Vec::new(),
        soft_field_values: Vec::new(),
        soft_gradients: Vec::new(),
        directed_normals: Vec::new(),
        soft_tangents: Vec::new(),
        covariance_groups: Vec::new(),
        shared_level_values: Vec::new(),
        field_energy: None,
        total_objective: None,
        backend_fingerprint: None,
        attempts: Vec::new(),
        recovery_verification: None,
        direct_input_conflicts: Vec::new(),
        relation_graph_conflicts: Vec::new(),
        shared_level_set_relation_conflicts: Vec::new(),
        execution_failure: None,
        cubic_analysis: None,
        backend_rank: None,
        interpretable_rank_deficiency: None,
        inertia: None,
        canonical_acceptance: None,
        capacity: None,
        analysis_failure: None,
        infeasibility_certificate: None,
        recession_ray: None,
        unidentified_additive_gauge: None,
        uninformative_shared_level_sets: Vec::new(),
        unresolved_axial_normals: Vec::new(),
    }
}

#[derive(Debug, Clone, Copy)]
struct ScalarRelationCounts {
    hard: usize,
    soft: usize,
}

impl ScalarRelationCounts {
    fn total(self) -> Option<usize> {
        self.hard.checked_add(self.soft)
    }
}

fn scalar_relation_counts(snapshot: &ProblemSnapshot) -> Option<ScalarRelationCounts> {
    let (observation_hard, observation_soft) = snapshot.inner.observations.iter().try_fold(
        (0_usize, 0_usize),
        |(hard, soft), observation| match observation {
            ObservationInput::FieldValue(value) if value.configuration().is_soft() => {
                Some((hard, soft.checked_add(1)?))
            }
            ObservationInput::Gradient(gradient) if gradient.configuration().is_soft() => {
                Some((hard, soft.checked_add(3)?))
            }
            ObservationInput::TangentDirection(tangent) if tangent.configuration().is_soft() => {
                Some((hard, soft.checked_add(1)?))
            }
            ObservationInput::FieldValue(_) | ObservationInput::TangentDirection(_) => {
                Some((hard.checked_add(1)?, soft))
            }
            ObservationInput::Gradient(_) => Some((hard.checked_add(3)?, soft)),
        },
    )?;
    let group_hard = snapshot
        .inner
        .shared_level_sets
        .iter()
        .try_fold(0_usize, |count, group| {
            count.checked_add(group.members().len())
        })?;
    let covariance_group_soft = snapshot
        .inner
        .covariance_groups
        .iter()
        .try_fold(0_usize, |count, group| {
            count.checked_add(group.scalar_residual_count())
        })?;
    let (normal_hard, normal_soft) = snapshot
        .inner
        .directed_normals
        .iter()
        .map(|normal| {
            (
                normal.direction(),
                normal.direction_enforcement(),
                normal.minimum_slope_enforcement(),
            )
        })
        .chain(snapshot.inner.axial_normals.iter().map(|normal| {
            (
                normal.axis(),
                normal.direction_enforcement(),
                normal.minimum_slope_enforcement(),
            )
        }))
        .try_fold(
            (0_usize, 0_usize),
            |(hard, soft), (normal, direction, slope)| {
                let projection_components = normal_projection_component_count(normal);
                let (hard, soft) = if direction.is_soft() {
                    (hard, soft.checked_add(projection_components)?)
                } else {
                    (hard.checked_add(projection_components)?, soft)
                };
                if slope.is_soft() {
                    Some((hard, soft.checked_add(1)?))
                } else {
                    Some((hard.checked_add(1)?, soft))
                }
            },
        )?;
    let (bound_hard, bound_soft) =
        affine_bound_sides(snapshot).try_fold((0_usize, 0_usize), |(hard, soft), side| {
            if side.configuration.is_soft() {
                Some((hard, soft.checked_add(1)?))
            } else {
                Some((hard.checked_add(1)?, soft))
            }
        })?;
    let (level_relation_hard, level_relation_soft) = snapshot
        .inner
        .shared_level_set_relations
        .iter()
        .try_fold((0_usize, 0_usize), |(hard, soft), relation| {
            if relation.is_soft() {
                Some((hard, soft.checked_add(1)?))
            } else {
                Some((hard.checked_add(1)?, soft))
            }
        })?;
    let (point_relation_hard, point_relation_soft) = snapshot
        .inner
        .point_to_level_set_relations
        .iter()
        .try_fold((0_usize, 0_usize), |(hard, soft), relation| {
            if relation.is_soft() {
                Some((hard, soft.checked_add(1)?))
            } else {
                Some((hard.checked_add(1)?, soft))
            }
        })?;
    Some(ScalarRelationCounts {
        hard: observation_hard
            .checked_add(group_hard)?
            .checked_add(snapshot.inner.additive_field_gauges.len())?
            .checked_add(bound_hard)?
            .checked_add(point_relation_hard)?
            .checked_add(level_relation_hard)?
            .checked_add(normal_hard)?,
        soft: observation_soft
            .checked_add(covariance_group_soft)?
            .checked_add(bound_soft)?
            .checked_add(point_relation_soft)?
            .checked_add(level_relation_soft)?
            .checked_add(normal_soft)?,
    })
}

fn quadratic_objective_term_count(snapshot: &ProblemSnapshot) -> Option<usize> {
    let independent = snapshot
        .inner
        .observations
        .iter()
        .filter(|observation| match observation {
            ObservationInput::FieldValue(value) => value.configuration().is_soft(),
            ObservationInput::Gradient(gradient) => gradient.configuration().is_soft(),
            ObservationInput::TangentDirection(tangent) => tangent.configuration().is_soft(),
        })
        .count();
    let bound_terms = affine_bound_sides(snapshot)
        .filter(|side| {
            matches!(
                side.configuration,
                AffineBoundConfiguration::QuadraticPenalty(_)
            )
        })
        .count();
    let level_relation_terms = snapshot
        .inner
        .shared_level_set_relations
        .iter()
        .filter(|relation| {
            matches!(
                relation.configuration(),
                AffineBoundConfiguration::QuadraticPenalty(_)
            )
        })
        .count();
    let point_relation_terms = snapshot
        .inner
        .point_to_level_set_relations
        .iter()
        .filter(|relation| {
            matches!(
                relation.configuration(),
                AffineBoundConfiguration::QuadraticPenalty(_)
            )
        })
        .count();
    let normal_terms = snapshot
        .inner
        .directed_normals
        .iter()
        .map(|normal| {
            usize::from(normal.direction_enforcement().is_soft())
                + usize::from(
                    normal
                        .minimum_slope_enforcement()
                        .quadratic_penalty()
                        .is_some(),
                )
        })
        .chain(snapshot.inner.axial_normals.iter().map(|normal| {
            usize::from(normal.direction_enforcement().is_soft())
                + usize::from(
                    normal
                        .minimum_slope_enforcement()
                        .quadratic_penalty()
                        .is_some(),
                )
        }))
        .try_fold(0_usize, usize::checked_add)?;
    independent
        .checked_add(snapshot.inner.covariance_groups.len())?
        .checked_add(bound_terms)?
        .checked_add(point_relation_terms)?
        .checked_add(level_relation_terms)?
        .checked_add(normal_terms)
}

fn linear_objective_term_count(snapshot: &ProblemSnapshot) -> Option<usize> {
    let bound_terms = affine_bound_sides(snapshot)
        .filter(|side| {
            matches!(
                side.configuration,
                AffineBoundConfiguration::LinearViolationPenalty(_)
            )
        })
        .count();
    let normal_slope_terms = snapshot
        .inner
        .directed_normals
        .iter()
        .filter(|normal| {
            normal
                .minimum_slope_enforcement()
                .linear_violation_penalty()
                .is_some()
        })
        .count()
        .checked_add(
            snapshot
                .inner
                .axial_normals
                .iter()
                .filter(|normal| {
                    normal
                        .minimum_slope_enforcement()
                        .linear_violation_penalty()
                        .is_some()
                })
                .count(),
        )?;
    bound_terms
        .checked_add(
            snapshot
                .inner
                .point_to_level_set_relations
                .iter()
                .filter(|relation| {
                    matches!(
                        relation.configuration(),
                        AffineBoundConfiguration::LinearViolationPenalty(_)
                    )
                })
                .count(),
        )?
        .checked_add(
            snapshot
                .inner
                .shared_level_set_relations
                .iter()
                .filter(|relation| {
                    matches!(
                        relation.configuration(),
                        AffineBoundConfiguration::LinearViolationPenalty(_)
                    )
                })
                .count(),
        )?
        .checked_add(normal_slope_terms)
}

fn affine_bound_sides(snapshot: &ProblemSnapshot) -> impl Iterator<Item = &AffineBoundSide> {
    snapshot
        .inner
        .field_value_bounds
        .iter()
        .flat_map(|bound| bound.lower().into_iter().chain(bound.upper()))
        .chain(
            snapshot
                .inner
                .directional_derivative_intervals
                .iter()
                .flat_map(|interval| interval.lower().into_iter().chain(interval.upper())),
        )
        .chain(
            snapshot
                .inner
                .field_separation_intervals
                .iter()
                .flat_map(|interval| [interval.lower(), interval.upper()]),
        )
}

fn source_identifier_bytes(snapshot: &ProblemSnapshot) -> Option<usize> {
    let observation_bytes =
        snapshot
            .inner
            .observations
            .iter()
            .try_fold(0_usize, |bytes, observation| {
                let multiplicity = match observation {
                    ObservationInput::FieldValue(_) | ObservationInput::TangentDirection(_) => 1,
                    ObservationInput::Gradient(_) => 3,
                };
                observation
                    .source_id()
                    .as_str()
                    .len()
                    .checked_mul(multiplicity)
                    .and_then(|source_bytes| bytes.checked_add(source_bytes))
            })?;
    let shared_level_bytes =
        snapshot
            .inner
            .shared_level_sets
            .iter()
            .try_fold(0_usize, |bytes, group| {
                group.members().iter().try_fold(bytes, |bytes, member| {
                    bytes
                        .checked_add(group.group_id().as_str().len())
                        .and_then(|bytes| bytes.checked_add(member.source_id().as_str().len()))
                })
            })?;
    let covariance_group_bytes =
        snapshot
            .inner
            .covariance_groups
            .iter()
            .try_fold(0_usize, |bytes, group| {
                group.members().iter().try_fold(bytes, |bytes, member| {
                    let multiplicity = member.scalar_residual_count();
                    group
                        .group_id()
                        .as_str()
                        .len()
                        .checked_add(member.source_id().as_str().len())?
                        .checked_mul(multiplicity)
                        .and_then(|member_bytes| bytes.checked_add(member_bytes))
                })
            })?;
    let gauge_bytes =
        snapshot
            .inner
            .additive_field_gauges
            .iter()
            .try_fold(0_usize, |bytes, gauge| {
                let group_bytes = match gauge.reference() {
                    AdditiveFieldGaugeReference::Point(_) => 0,
                    AdditiveFieldGaugeReference::LevelSet(group_id) => group_id.as_str().len(),
                };
                bytes
                    .checked_add(gauge.source_id().as_str().len())
                    .and_then(|bytes| bytes.checked_add(group_bytes))
            })?;
    let bound_bytes =
        snapshot
            .inner
            .field_value_bounds
            .iter()
            .try_fold(0_usize, |bytes, bound| {
                let multiplicity = usize::from(bound.lower().is_some())
                    .checked_add(usize::from(bound.upper().is_some()))?;
                bound
                    .source_id()
                    .as_str()
                    .len()
                    .checked_mul(multiplicity)
                    .and_then(|source_bytes| bytes.checked_add(source_bytes))
            })?;
    let derivative_interval_bytes = snapshot
        .inner
        .directional_derivative_intervals
        .iter()
        .try_fold(0_usize, |bytes, interval| {
            let multiplicity = usize::from(interval.lower().is_some())
                .checked_add(usize::from(interval.upper().is_some()))?;
            interval
                .source_id()
                .as_str()
                .len()
                .checked_mul(multiplicity)
                .and_then(|source_bytes| bytes.checked_add(source_bytes))
        })?;
    let shared_level_set_relation_bytes = snapshot
        .inner
        .shared_level_set_relations
        .iter()
        .try_fold(0_usize, |bytes, relation| {
            relation.declared_group_ids().into_iter().try_fold(
                bytes.checked_add(relation.source_id().as_str().len())?,
                |bytes, group_id| bytes.checked_add(group_id.as_str().len()),
            )
        })?;
    let field_separation_bytes =
        snapshot
            .inner
            .field_separation_intervals
            .iter()
            .try_fold(0_usize, |bytes, relation| {
                bytes
                    .checked_add(relation.source_id().as_str().len().checked_mul(2)?)?
                    .checked_add(
                        relation
                            .reference_group_id()
                            .as_str()
                            .len()
                            .checked_mul(2)?,
                    )?
                    .checked_add(relation.target_group_id().as_str().len().checked_mul(2)?)
            })?;
    let point_to_level_set_bytes = snapshot
        .inner
        .point_to_level_set_relations
        .iter()
        .try_fold(0_usize, |bytes, relation| {
            bytes
                .checked_add(relation.source_id().as_str().len())?
                .checked_add(relation.group_id().as_str().len())
        })?;
    let normal_bytes = snapshot
        .inner
        .directed_normals
        .iter()
        .map(|normal| normal.source_id())
        .chain(
            snapshot
                .inner
                .axial_normals
                .iter()
                .map(|normal| normal.source_id()),
        )
        .try_fold(0_usize, |bytes, source_id| {
            source_id
                .as_str()
                .len()
                .checked_mul(4)
                .and_then(|source_bytes| bytes.checked_add(source_bytes))
        })?;
    let polarity_resolution_bytes =
        snapshot
            .inner
            .polarity_resolutions
            .iter()
            .try_fold(0_usize, |bytes, resolution| {
                bytes
                    .checked_add(resolution.source_id().as_str().len())?
                    .checked_add(resolution.axial_normal_source_id().as_str().len())
            })?;
    observation_bytes
        .checked_add(shared_level_bytes)?
        .checked_add(covariance_group_bytes)?
        .checked_add(gauge_bytes)?
        .checked_add(bound_bytes)?
        .checked_add(derivative_interval_bytes)?
        .checked_add(field_separation_bytes)?
        .checked_add(point_to_level_set_bytes)?
        .checked_add(shared_level_set_relation_bytes)?
        .checked_add(normal_bytes)?
        .checked_add(polarity_resolution_bytes)
}

fn snapshot_references_group(snapshot: &ProblemSnapshot, group_id: &GroupId) -> bool {
    snapshot.inner.additive_field_gauges.iter().any(|gauge| {
        matches!(
            gauge.reference(),
            AdditiveFieldGaugeReference::LevelSet(referenced) if referenced == group_id
        )
    }) || snapshot
        .inner
        .shared_level_set_relations
        .iter()
        .any(|relation| relation.declared_group_ids().contains(&group_id))
        || snapshot
            .inner
            .field_separation_intervals
            .iter()
            .any(|relation| {
                relation.reference_group_id() == group_id || relation.target_group_id() == group_id
            })
        || snapshot
            .inner
            .point_to_level_set_relations
            .iter()
            .any(|relation| relation.group_id() == group_id)
}

fn conservative_problem_size(
    input_observations: usize,
    scalar_relations: ScalarRelationCounts,
    semantic_latents: usize,
    quadratic_objective_terms: usize,
    linear_objective_terms: usize,
    affine_inequality_constraints: usize,
) -> ProblemSize {
    ProblemSize {
        input_observations,
        scalar_hard_relations: scalar_relations.hard,
        scalar_soft_relations: scalar_relations.soft,
        canonical_hard_equalities: None,
        canonical_soft_equalities: None,
        center_coefficients: None,
        semantic_latents,
        auxiliary_variables: 0,
        quadratic_objective_terms,
        linear_objective_terms,
        affine_inequality_constraints,
        cone_blocks: 0,
        primal_variables: None,
        equality_constraints: None,
        kkt_dimension: None,
    }
}

fn plan_snapshot_capacity(
    lowering: &EqualityLowering,
    fitting_functional_count: usize,
    scalar_relations: usize,
    source_identifier_bytes: usize,
) -> Result<(), CapacityExceededEvidence> {
    let Some(primal_variables) = fitting_functional_count
        .checked_add(4)
        .and_then(|value| value.checked_add(lowering.semantic_latents.len()))
    else {
        return Err(plan_equality_capacity(usize::MAX, usize::MAX)
            .expect_err("maximal dimensions must overflow the capacity plan"));
    };
    let Some(equality_constraints) = lowering.solver_equality_count().checked_add(4) else {
        return Err(plan_equality_capacity(usize::MAX, usize::MAX)
            .expect_err("maximal dimensions must overflow the capacity plan"));
    };
    let Some(canonical_relations) = lowering
        .canonical_equalities
        .len()
        .checked_add(lowering.canonical_soft_equalities.len())
        .and_then(|value| value.checked_add(4))
    else {
        return Err(plan_equality_capacity(usize::MAX, usize::MAX)
            .expect_err("maximal dimensions must overflow the capacity plan"));
    };
    plan_equality_capacity_for(EqualityCapacityShape {
        primal_variables,
        equality_constraints,
        canonical_relations,
        source_lowering: SourceStorageShape {
            relations: scalar_relations,
            identifier_bytes: source_identifier_bytes,
        },
        source_report: SourceStorageShape {
            relations: scalar_relations,
            identifier_bytes: source_identifier_bytes,
        },
    })
    .map(|_| ())
}

fn plan_source_lifecycle_capacity(
    scalar_relations: usize,
    source_identifier_bytes: usize,
) -> Result<(), CapacityExceededEvidence> {
    plan_equality_capacity_for(EqualityCapacityShape {
        primal_variables: 0,
        equality_constraints: 0,
        canonical_relations: 0,
        source_lowering: SourceStorageShape {
            relations: scalar_relations,
            identifier_bytes: source_identifier_bytes,
        },
        source_report: SourceStorageShape {
            relations: scalar_relations,
            identifier_bytes: source_identifier_bytes,
        },
    })
    .map(|_| ())
}

fn plan_source_preflight_capacity(
    scalar_relations: usize,
    source_identifier_bytes: usize,
) -> Result<(), CapacityExceededEvidence> {
    plan_equality_capacity_for(EqualityCapacityShape {
        primal_variables: 0,
        equality_constraints: 0,
        canonical_relations: 0,
        source_lowering: SourceStorageShape {
            relations: scalar_relations,
            identifier_bytes: source_identifier_bytes,
        },
        source_report: SourceStorageShape::default(),
    })
    .map(|_| ())
}

fn fit_snapshot_after_preflight(
    snapshot: &ProblemSnapshot,
    lowering: EqualityLowering,
    problem_size: ProblemSize,
    base_report: impl Fn() -> FitReport,
) -> Result<FitSuccess, FitFailure> {
    if !lowering.canonical_affine_inequalities.is_empty() {
        return fit_snapshot_after_preflight_qp(snapshot, lowering, problem_size, base_report());
    }
    let solution = match CubicExecutionCore::solve_equality_production(
        CubicCanonicalProblem {
            equalities: lowering.canonical_equalities.clone(),
            hard_residual_blocks: lowering.canonical_hard_residual_blocks.clone(),
            affine_inequalities: Vec::new(),
            soft_equalities: lowering.canonical_soft_equalities.clone(),
            soft_objectives: lowering.canonical_soft_objectives.clone(),
            semantic_latents: lowering.semantic_latents.clone(),
            field_energy_normalization: snapshot.inner.field_energy_normalization,
        },
        snapshot.inner.global_anisotropy_metric.as_cubic_metric(),
    ) {
        Ok(solution) => solution,
        Err(failure) => {
            return Err(FitFailure {
                diagnosis: diagnose(&failure),
                report: Box::new(failure_report(
                    base_report(),
                    &failure,
                    &lowering.source_relations,
                )),
            });
        }
    };
    let shared_level_values = solution
        .semantic_latents
        .iter()
        .map(|latent| SharedLevelValue {
            group_id: latent.group_id.clone(),
            value: latent.value,
            field_unit: latent.field_unit.clone(),
            member_source_ids: latent.member_source_ids.clone().into(),
        })
        .collect::<Vec<_>>();
    let report = success_report(
        snapshot,
        problem_size,
        &solution,
        &lowering.source_relations,
        &shared_level_values,
    );
    let model = SolvedModel::new(
        snapshot.clone(),
        solution.field,
        shared_level_values
            .iter()
            .map(|level| (level.group_id.clone(), level.value))
            .collect(),
    );
    Ok(FitSuccess { model, report })
}

fn fit_snapshot_after_preflight_qp(
    snapshot: &ProblemSnapshot,
    lowering: EqualityLowering,
    problem_size: ProblemSize,
    base_report: FitReport,
) -> Result<FitSuccess, FitFailure> {
    let solution = match CubicExecutionCore::solve(
        CubicCanonicalProblem {
            equalities: lowering.canonical_equalities.clone(),
            hard_residual_blocks: lowering.canonical_hard_residual_blocks.clone(),
            affine_inequalities: lowering.canonical_affine_inequalities.clone(),
            soft_equalities: lowering.canonical_soft_equalities.clone(),
            soft_objectives: lowering.canonical_soft_objectives.clone(),
            semantic_latents: lowering.semantic_latents.clone(),
            field_energy_normalization: snapshot.inner.field_energy_normalization,
        },
        snapshot.inner.global_anisotropy_metric.as_cubic_metric(),
    ) {
        Ok(solution) => solution,
        Err(failure) => {
            let diagnosis = diagnose_qp(&failure);
            return Err(FitFailure {
                diagnosis,
                report: Box::new(qp_failure_report(
                    base_report,
                    &failure,
                    &lowering.source_relations,
                )),
            });
        }
    };
    let shared_level_values = solution
        .semantic_latents
        .iter()
        .map(|latent| SharedLevelValue {
            group_id: latent.group_id.clone(),
            value: latent.value,
            field_unit: latent.field_unit.clone(),
            member_source_ids: latent.member_source_ids.clone().into(),
        })
        .collect::<Vec<_>>();
    let report = success_report_qp(
        snapshot,
        problem_size,
        &solution,
        &lowering.source_relations,
        &lowering.source_bound_relations,
        &lowering.resolved_normals,
        &shared_level_values,
    );
    let model = SolvedModel::new(
        snapshot.clone(),
        solution.field,
        shared_level_values
            .iter()
            .map(|level| (level.group_id.clone(), level.value))
            .collect(),
    );
    Ok(FitSuccess { model, report })
}

#[derive(Debug, Clone)]
struct ScalarObservation {
    source_id: SourceId,
    group_id: Option<GroupId>,
    support: [f64; 3],
    functional: ScalarFunctionalDescriptor,
    semantic_role: SemanticRolePath,
    target: f64,
}

#[derive(Debug, Clone)]
struct ObservationResidualBlock {
    components: Vec<ScalarObservation>,
    configuration: ResidualBlockConfiguration,
    kind: CanonicalSoftResidualBlockKind,
}

#[derive(Debug, Clone)]
struct ResolvedNormalInput {
    source_id: SourceId,
    location: Point3,
    direction: Vector3,
    direction_enforcement: NormalDirectionEnforcement,
    minimum_slope: MinimumNormalSlope,
    minimum_slope_enforcement: MinimumNormalSlopeEnforcement,
    input_axis: Option<Vector3>,
    polarity_resolution_source_id: Option<SourceId>,
    polarity_selection: Option<PolaritySelection>,
}

fn resolved_normal_inputs(snapshot: &ProblemSnapshot) -> Vec<ResolvedNormalInput> {
    let mut normals = snapshot
        .inner
        .directed_normals
        .iter()
        .map(|normal| ResolvedNormalInput {
            source_id: normal.source_id().clone(),
            location: normal.location(),
            direction: normal.direction(),
            direction_enforcement: normal.direction_enforcement(),
            minimum_slope: normal.minimum_slope(),
            minimum_slope_enforcement: normal.minimum_slope_enforcement(),
            input_axis: None,
            polarity_resolution_source_id: None,
            polarity_selection: None,
        })
        .collect::<Vec<_>>();
    normals.extend(snapshot.inner.axial_normals.iter().filter_map(|normal| {
        let resolution = snapshot
            .inner
            .polarity_resolutions
            .iter()
            .find(|resolution| resolution.axial_normal_source_id() == normal.source_id())?;
        let direction = match resolution.selection() {
            PolaritySelection::AlongInputAxis => normal.input_axis(),
            PolaritySelection::AgainstInputAxis => {
                let [x, y, z] = normal.input_axis().components();
                Vector3::try_new(-x, -y, -z)
                    .expect("reversing a finite unit normal keeps it finite")
            }
        };
        Some(ResolvedNormalInput {
            source_id: normal.source_id().clone(),
            location: normal.location(),
            direction,
            direction_enforcement: normal.direction_enforcement(),
            minimum_slope: normal.minimum_slope(),
            minimum_slope_enforcement: normal.minimum_slope_enforcement(),
            input_axis: Some(normal.input_axis()),
            polarity_resolution_source_id: Some(resolution.source_id().clone()),
            polarity_selection: Some(resolution.selection()),
        })
    }));
    normals.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    normals
}

#[derive(Debug, Clone)]
enum ResidualBlockConfiguration {
    Hard,
    QuadraticPenalty(QuadraticPenalty),
    StandardDeviation(StandardDeviation),
    Covariance(CovarianceMatrix),
}

#[derive(Debug, Clone, Copy)]
struct ScalarFunctionalDescriptor {
    dimension: FunctionalDimension,
    value_coefficient: f64,
    gradient_coefficient: [f64; 3],
}

impl ScalarFunctionalDescriptor {
    fn field_value() -> Self {
        Self {
            dimension: FunctionalDimension::FieldValue,
            value_coefficient: 1.0,
            gradient_coefficient: [0.0; 3],
        }
    }

    fn gradient_component(axis: usize) -> Self {
        Self {
            dimension: FunctionalDimension::FieldValuePerLength,
            value_coefficient: 0.0,
            gradient_coefficient: std::array::from_fn(
                |component| {
                    if component == axis { 1.0 } else { 0.0 }
                },
            ),
        }
    }

    fn directional_derivative(direction: [f64; 3]) -> Self {
        Self {
            dimension: FunctionalDimension::FieldValuePerLength,
            value_coefficient: 0.0,
            gradient_coefficient: direction,
        }
    }
}

fn observation_residual_blocks(observations: &[ObservationInput]) -> Vec<ObservationResidualBlock> {
    observations
        .iter()
        .map(|observation| match observation {
            ObservationInput::FieldValue(observation) => ObservationResidualBlock {
                components: vec![ScalarObservation {
                    source_id: observation.source_id().clone(),
                    group_id: None,
                    support: observation.location().components(),
                    functional: ScalarFunctionalDescriptor::field_value(),
                    semantic_role: SemanticRolePath::new("field-value-observation/value"),
                    target: observation.value(),
                }],
                configuration: match observation.configuration() {
                    FieldValueConfiguration::Hard => ResidualBlockConfiguration::Hard,
                    FieldValueConfiguration::QuadraticPenalty(penalty) => {
                        ResidualBlockConfiguration::QuadraticPenalty(penalty)
                    }
                    FieldValueConfiguration::StandardDeviation(standard_deviation) => {
                        ResidualBlockConfiguration::StandardDeviation(standard_deviation)
                    }
                },
                kind: CanonicalSoftResidualBlockKind::Independent(
                    CanonicalSoftResidualMemberKind::FieldValue,
                ),
            },
            ObservationInput::Gradient(observation) => ObservationResidualBlock {
                components: observation
                    .gradient()
                    .components()
                    .into_iter()
                    .enumerate()
                    .map(|(axis, target)| ScalarObservation {
                        source_id: observation.source_id().clone(),
                        group_id: None,
                        support: observation.location().components(),
                        functional: ScalarFunctionalDescriptor::gradient_component(axis),
                        semantic_role: SemanticRolePath::new(format!(
                            "gradient-observation/component/{axis}"
                        )),
                        target,
                    })
                    .collect(),
                configuration: match observation.configuration() {
                    GradientConfiguration::Hard => ResidualBlockConfiguration::Hard,
                    GradientConfiguration::QuadraticPenalty(penalty) => {
                        ResidualBlockConfiguration::QuadraticPenalty(*penalty)
                    }
                    GradientConfiguration::StandardDeviation(standard_deviation) => {
                        ResidualBlockConfiguration::StandardDeviation(*standard_deviation)
                    }
                    GradientConfiguration::Covariance(covariance) => {
                        ResidualBlockConfiguration::Covariance(covariance.clone())
                    }
                },
                kind: CanonicalSoftResidualBlockKind::Independent(
                    CanonicalSoftResidualMemberKind::Gradient,
                ),
            },
            ObservationInput::TangentDirection(observation) => ObservationResidualBlock {
                components: vec![ScalarObservation {
                    source_id: observation.source_id().clone(),
                    group_id: None,
                    support: observation.location().components(),
                    functional: ScalarFunctionalDescriptor::directional_derivative(
                        observation.direction().components(),
                    ),
                    semantic_role: SemanticRolePath::new(
                        "tangent-direction-observation/directional-derivative",
                    ),
                    target: 0.0,
                }],
                configuration: match observation.configuration() {
                    TangentConfiguration::Hard => ResidualBlockConfiguration::Hard,
                    TangentConfiguration::QuadraticPenalty(penalty) => {
                        ResidualBlockConfiguration::QuadraticPenalty(penalty)
                    }
                    TangentConfiguration::StandardDeviation(standard_deviation) => {
                        ResidualBlockConfiguration::StandardDeviation(standard_deviation)
                    }
                },
                kind: CanonicalSoftResidualBlockKind::Independent(
                    CanonicalSoftResidualMemberKind::Tangent,
                ),
            },
        })
        .collect()
}

fn normal_projection_components(normal: &ResolvedNormalInput) -> Vec<(usize, ScalarObservation)> {
    let direction = normal.direction.components();
    (0..3)
        .filter_map(|axis| {
            let coefficients = std::array::from_fn(|component| {
                let identity = if component == axis { 1.0 } else { 0.0 };
                crate::math::canonical_zero(identity - direction[axis] * direction[component])
            });
            coefficients
                .iter()
                .any(|coefficient| *coefficient != 0.0)
                .then(|| {
                    (
                        axis,
                        ScalarObservation {
                            source_id: normal.source_id.clone(),
                            group_id: None,
                            support: normal.location.components(),
                            functional: ScalarFunctionalDescriptor::directional_derivative(
                                coefficients,
                            ),
                            semantic_role: SemanticRolePath::new(format!(
                                "directed-normal/direction-projection/component/{axis}"
                            )),
                            target: 0.0,
                        },
                    )
                })
        })
        .collect()
}

fn normal_projection_component_count(direction: Vector3) -> usize {
    let direction = direction.components();
    (0..3)
        .filter(|axis| {
            (0..3).any(|component| {
                let identity = if component == *axis { 1.0 } else { 0.0 };
                crate::math::canonical_zero(identity - direction[*axis] * direction[component])
                    != 0.0
            })
        })
        .count()
}

fn covariance_group_components(
    group_id: &GroupId,
    members: &[CovarianceGroupMember],
) -> Vec<ScalarObservation> {
    members
        .iter()
        .flat_map(|member| match member {
            CovarianceGroupMember::FieldValue(observation) => vec![ScalarObservation {
                source_id: observation.source_id().clone(),
                group_id: Some(group_id.clone()),
                support: observation.location().components(),
                functional: ScalarFunctionalDescriptor::field_value(),
                semantic_role: SemanticRolePath::new("field-value-observation/value"),
                target: observation.value(),
            }],
            CovarianceGroupMember::Gradient(observation) => observation
                .gradient()
                .components()
                .into_iter()
                .enumerate()
                .map(|(axis, target)| ScalarObservation {
                    source_id: observation.source_id().clone(),
                    group_id: Some(group_id.clone()),
                    support: observation.location().components(),
                    functional: ScalarFunctionalDescriptor::gradient_component(axis),
                    semantic_role: SemanticRolePath::new(format!(
                        "gradient-observation/component/{axis}"
                    )),
                    target,
                })
                .collect(),
            CovarianceGroupMember::Tangent(observation) => vec![ScalarObservation {
                source_id: observation.source_id().clone(),
                group_id: Some(group_id.clone()),
                support: observation.location().components(),
                functional: ScalarFunctionalDescriptor::directional_derivative(
                    observation.direction().components(),
                ),
                semantic_role: SemanticRolePath::new(
                    "tangent-direction-observation/directional-derivative",
                ),
                target: 0.0,
            }],
        })
        .collect()
}

fn covariance_group_member_kinds(
    members: &[CovarianceGroupMember],
) -> Vec<CanonicalSoftResidualMemberKind> {
    members
        .iter()
        .map(|member| match member {
            CovarianceGroupMember::FieldValue(_) => CanonicalSoftResidualMemberKind::FieldValue,
            CovarianceGroupMember::Gradient(_) => CanonicalSoftResidualMemberKind::Gradient,
            CovarianceGroupMember::Tangent(_) => CanonicalSoftResidualMemberKind::Tangent,
        })
        .collect()
}

#[derive(Debug, Clone)]
struct SourceHardRelation {
    equality: CanonicalHardEquality,
    canonical_index: usize,
    kind: SourceHardRelationKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceHardRelationKind {
    Scalar,
    NormalProjection { normal_index: usize },
}

#[derive(Debug, Clone)]
struct SourceBoundRelation {
    inequality: CanonicalAffineInequality,
    canonical_index: usize,
    kind: SourceBoundKind,
}

#[derive(Debug, Clone)]
enum SourceBoundKind {
    FieldValue,
    DirectionalDerivative {
        direction: Vector3,
    },
    FieldSeparationInterval {
        relation: FieldSeparationInterval,
    },
    PointToLevelSetRelation {
        relation: PointToLevelSetRelation,
    },
    SharedLevelSetRelation {
        relation: SharedLevelSetRelationInput,
        orientation: SharedLevelSetRelationOrientation,
    },
    NormalSlope {
        normal_index: usize,
    },
}

#[derive(Debug, Clone)]
struct EqualityLowering {
    source_relations: Vec<SourceHardRelation>,
    source_bound_relations: Vec<SourceBoundRelation>,
    canonical_equalities: Vec<CanonicalHardEquality>,
    canonical_hard_residual_blocks: Vec<CanonicalHardResidualBlock>,
    canonical_affine_inequalities: Vec<CanonicalAffineInequality>,
    canonical_soft_equalities: Vec<CanonicalSoftEquality>,
    canonical_soft_objectives: Vec<CanonicalSoftObjective>,
    canonical_index_by_key: BTreeMap<CanonicalEqualityKey, (usize, f64)>,
    canonical_bound_index_by_key: BTreeMap<CanonicalBoundKey, usize>,
    direct_input_conflicts: Vec<DirectInputConflictEvidence>,
    relation_graph_conflicts: Vec<RelationGraphConflictEvidence>,
    semantic_latents: Vec<SemanticLatentDefinition>,
    resolved_normals: Vec<ResolvedNormalInput>,
}

impl EqualityLowering {
    fn new(resolved_normals: Vec<ResolvedNormalInput>) -> Self {
        Self {
            source_relations: Vec::new(),
            source_bound_relations: Vec::new(),
            canonical_equalities: Vec::new(),
            canonical_hard_residual_blocks: Vec::new(),
            canonical_affine_inequalities: Vec::new(),
            canonical_soft_equalities: Vec::new(),
            canonical_soft_objectives: Vec::new(),
            canonical_index_by_key: BTreeMap::new(),
            canonical_bound_index_by_key: BTreeMap::new(),
            direct_input_conflicts: Vec::new(),
            relation_graph_conflicts: Vec::new(),
            semantic_latents: Vec::new(),
            resolved_normals,
        }
    }

    fn push_soft_block(
        &mut self,
        equalities: Vec<CanonicalSoftEquality>,
        loss: CanonicalSoftLoss,
        covariance_group: Option<GroupId>,
        block_kind: CanonicalSoftResidualBlockKind,
    ) {
        self.canonical_soft_objectives
            .push(CanonicalSoftObjective::new_block(
                equalities
                    .iter()
                    .map(|equality| equality.provenance().residual().clone())
                    .collect(),
                loss,
                covariance_group,
                block_kind,
            ));
        self.canonical_soft_equalities.extend(equalities);
    }

    fn push_source(&mut self, equality: CanonicalHardEquality) {
        self.push_source_with_kind(equality, SourceHardRelationKind::Scalar);
    }

    fn push_source_with_kind(
        &mut self,
        equality: CanonicalHardEquality,
        kind: SourceHardRelationKind,
    ) -> usize {
        let (key, normalized_target) = normalized_equality_key(&equality);
        let canonical_index =
            if let Some((index, first_target)) = self.canonical_index_by_key.get(&key).copied() {
                if first_target != normalized_target {
                    let first = &self.canonical_equalities[index];
                    self.direct_input_conflicts
                        .push(DirectInputConflictEvidence::new(
                            first.provenance().source().clone(),
                            equality.provenance().source().clone(),
                            equality.provenance().semantic_role().clone(),
                            first.target(),
                            equality.target(),
                        ));
                    self.direct_input_conflicts.sort_by(|left, right| {
                        left.semantic_role()
                            .cmp(right.semantic_role())
                            .then_with(|| left.first_source().cmp(right.first_source()))
                            .then_with(|| left.second_source().cmp(right.second_source()))
                    });
                }
                if equality.participation() == CanonicalEqualityParticipation::SolverConstraint {
                    self.canonical_equalities[index].promote_to_solver_constraint();
                }
                index
            } else {
                let index = self.canonical_equalities.len();
                self.canonical_equalities.push(equality.clone());
                self.canonical_index_by_key
                    .insert(key, (index, normalized_target));
                index
            };
        self.source_relations.push(SourceHardRelation {
            equality,
            canonical_index,
            kind,
        });
        canonical_index
    }

    fn push_bound(&mut self, inequality: CanonicalAffineInequality, kind: SourceBoundKind) {
        let canonical_index = if inequality.violation_channel().is_some() {
            let index = self.canonical_affine_inequalities.len();
            self.canonical_affine_inequalities.push(inequality.clone());
            index
        } else {
            let key = CanonicalBoundKey {
                dimension: inequality.dimension(),
                field_terms: inequality
                    .field()
                    .into_iter()
                    .flat_map(|field| field.functional().terms())
                    .map(|term| {
                        (
                            canonical_support_bits(term.support()),
                            canonical_scalar_bits(term.value_coefficient()),
                            term.gradient_coefficient().map(canonical_scalar_bits),
                        )
                    })
                    .collect(),
                latent_coefficients: inequality
                    .latent_coefficients()
                    .iter()
                    .map(|term| (term.latent, canonical_scalar_bits(term.coefficient)))
                    .collect::<Vec<_>>(),
                side: inequality.sense(),
                bound: canonical_scalar_bits(inequality.bound()),
            };
            if let Some(index) = self.canonical_bound_index_by_key.get(&key) {
                self.canonical_affine_inequalities[*index]
                    .add_source_provenance(inequality.provenance().clone());
                *index
            } else {
                let index = self.canonical_affine_inequalities.len();
                self.canonical_affine_inequalities.push(inequality.clone());
                self.canonical_bound_index_by_key.insert(key, index);
                index
            }
        };
        self.source_bound_relations.push(SourceBoundRelation {
            inequality,
            canonical_index,
            kind,
        });
    }

    fn record_graph_conflict(
        &mut self,
        outcome: &CanonicalValueEdgeOutcome,
        semantic_role: &SemanticRolePath,
    ) {
        let CanonicalValueEdgeOutcome::Conflict {
            proof_source_ids,
            proof_group_ids,
            first_absolute,
            second_absolute,
        } = outcome
        else {
            return;
        };
        self.relation_graph_conflicts
            .push(RelationGraphConflictEvidence::new(
                proof_source_ids.clone(),
                proof_group_ids.clone(),
                semantic_role.clone(),
                first_absolute.source_id.clone(),
                first_absolute.target,
                second_absolute.source_id.clone(),
                second_absolute.target,
            ));
        self.relation_graph_conflicts.sort_by(|left, right| {
            left.semantic_role()
                .cmp(right.semantic_role())
                .then_with(|| left.source_ids().cmp(right.source_ids()))
        });
    }

    fn record_direct_derivative_bound_conflicts(&mut self, facts: &[HardDerivativeBoundFact]) {
        let mut by_functional =
            BTreeMap::<CanonicalEqualityKey, Vec<&HardDerivativeBoundFact>>::new();
        for fact in facts {
            by_functional
                .entry(fact.functional_key.clone())
                .or_default()
                .push(fact);
        }
        for (functional_key, functional_facts) in by_functional {
            for (index, left) in functional_facts.iter().enumerate() {
                for right in &functional_facts[index + 1..] {
                    let incompatible = match (left.sense, right.sense) {
                        (CanonicalInequalitySense::Lower, CanonicalInequalitySense::Upper) => {
                            left.bound > right.bound
                        }
                        (CanonicalInequalitySense::Upper, CanonicalInequalitySense::Lower) => {
                            right.bound > left.bound
                        }
                        _ => false,
                    };
                    if incompatible {
                        let (first, second) = if left.source_id <= right.source_id {
                            (left, right)
                        } else {
                            (right, left)
                        };
                        self.direct_input_conflicts
                            .push(DirectInputConflictEvidence::new(
                                first.source_id.clone(),
                                second.source_id.clone(),
                                second.semantic_role.clone(),
                                first.bound,
                                second.bound,
                            ));
                    }
                }
            }
            for source_relation in &self.source_relations {
                let (equality_key, equality_target) =
                    normalized_equality_key(&source_relation.equality);
                if equality_key != functional_key {
                    continue;
                }
                for fact in &functional_facts {
                    let incompatible = match fact.sense {
                        CanonicalInequalitySense::Lower => equality_target < fact.bound,
                        CanonicalInequalitySense::Upper => equality_target > fact.bound,
                    };
                    if incompatible {
                        let equality_source = source_relation.equality.provenance().source();
                        let (first_source, first_target, second_source, second_target) =
                            if equality_source <= &fact.source_id {
                                (
                                    equality_source.clone(),
                                    equality_target,
                                    fact.source_id.clone(),
                                    fact.bound,
                                )
                            } else {
                                (
                                    fact.source_id.clone(),
                                    fact.bound,
                                    equality_source.clone(),
                                    equality_target,
                                )
                            };
                        self.direct_input_conflicts
                            .push(DirectInputConflictEvidence::new(
                                first_source,
                                second_source,
                                fact.semantic_role.clone(),
                                first_target,
                                second_target,
                            ));
                    }
                }
            }
        }
        self.direct_input_conflicts.sort_by(|left, right| {
            left.semantic_role()
                .cmp(right.semantic_role())
                .then_with(|| left.first_source().cmp(right.first_source()))
                .then_with(|| left.second_source().cmp(right.second_source()))
        });
    }

    fn solver_equality_count(&self) -> usize {
        self.canonical_equalities
            .iter()
            .filter(|equality| {
                equality.participation() == CanonicalEqualityParticipation::SolverConstraint
            })
            .count()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CanonicalEqualityKey {
    dimension: FunctionalDimension,
    field_terms: Vec<([u64; 3], u64, [u64; 3])>,
    latent_coefficients: Vec<(usize, u64)>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CanonicalBoundKey {
    dimension: FunctionalDimension,
    field_terms: Vec<([u64; 3], u64, [u64; 3])>,
    latent_coefficients: Vec<(usize, u64)>,
    side: CanonicalInequalitySense,
    bound: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum CanonicalValueNode {
    Group(GroupId),
    Support([u64; 3]),
}

#[derive(Debug, Clone, PartialEq)]
struct ComponentAbsoluteTarget {
    node: usize,
    source_id: SourceId,
    target: f64,
}

#[derive(Debug, Clone)]
struct AbsoluteGroupProof {
    target: f64,
    source_ids: Vec<SourceId>,
}

#[derive(Debug, Default)]
/// Equality components retain original absolute anchors, never derived offsets.
struct CanonicalValueConstraintForest {
    node_index: BTreeMap<CanonicalValueNode, usize>,
    nodes: Vec<CanonicalValueNode>,
    parent: Vec<usize>,
    absolute_target: Vec<Option<ComponentAbsoluteTarget>>,
    adjacency: Vec<Vec<(usize, SourceId)>>,
}

impl CanonicalValueConstraintForest {
    fn absolute_target_for_support(
        &mut self,
        support: [u64; 3],
    ) -> Option<ComponentAbsoluteTarget> {
        let index = *self.node_index.get(&CanonicalValueNode::Support(support))?;
        let root = self.root(index);
        self.absolute_target[root].clone()
    }

    fn absolute_proof_for_group(&mut self, group_id: &GroupId) -> Option<AbsoluteGroupProof> {
        let node = *self
            .node_index
            .get(&CanonicalValueNode::Group(group_id.clone()))?;
        let root = self.root(node);
        let absolute = self.absolute_target[root].clone()?;
        let mut source_ids = self.proof_source_ids(node, absolute.node);
        source_ids.push(absolute.source_id);
        Some(AbsoluteGroupProof {
            target: absolute.target,
            source_ids,
        })
    }

    fn add_member_equality(
        &mut self,
        group_id: &GroupId,
        support: [f64; 3],
        source_id: &SourceId,
    ) -> CanonicalValueEdgeOutcome {
        self.add_equality_edge(
            CanonicalValueNode::Group(group_id.clone()),
            CanonicalValueNode::Support(canonical_support_bits(support)),
            source_id,
        )
    }

    fn add_absolute_support(
        &mut self,
        support: [f64; 3],
        value: f64,
        source_id: &SourceId,
    ) -> CanonicalValueEdgeOutcome {
        self.add_absolute_target(
            CanonicalValueNode::Support(canonical_support_bits(support)),
            value,
            source_id,
        )
    }

    fn add_absolute_group(
        &mut self,
        group_id: &GroupId,
        value: f64,
        source_id: &SourceId,
    ) -> CanonicalValueEdgeOutcome {
        self.add_absolute_target(
            CanonicalValueNode::Group(group_id.clone()),
            value,
            source_id,
        )
    }

    fn add_equality_edge(
        &mut self,
        left: CanonicalValueNode,
        right: CanonicalValueNode,
        source_id: &SourceId,
    ) -> CanonicalValueEdgeOutcome {
        let left = self.intern(left);
        let right = self.intern(right);
        let left_root = self.root(left);
        let right_root = self.root(right);
        if left_root == right_root {
            return CanonicalValueEdgeOutcome::Redundant;
        }
        let left_absolute = self.absolute_target[left_root].clone();
        let right_absolute = self.absolute_target[right_root].clone();
        let both_sides_anchored = left_absolute.is_some() && right_absolute.is_some();
        if let (Some(left_absolute), Some(right_absolute)) = (&left_absolute, &right_absolute) {
            if left_absolute.target != right_absolute.target {
                return self.equality_conflict(
                    left,
                    right,
                    source_id,
                    left_absolute,
                    right_absolute,
                );
            }
        }
        self.parent[right_root] = left_root;
        self.absolute_target[left_root] = left_absolute.or(right_absolute);
        self.absolute_target[right_root] = None;
        self.adjacency[left].push((right, source_id.clone()));
        self.adjacency[right].push((left, source_id.clone()));
        if both_sides_anchored {
            CanonicalValueEdgeOutcome::Redundant
        } else {
            CanonicalValueEdgeOutcome::Independent
        }
    }

    fn add_absolute_target(
        &mut self,
        node: CanonicalValueNode,
        target: f64,
        source_id: &SourceId,
    ) -> CanonicalValueEdgeOutcome {
        let node = self.intern(node);
        let root = self.root(node);
        if let Some(existing) = &self.absolute_target[root] {
            return if existing.target == target {
                CanonicalValueEdgeOutcome::Redundant
            } else {
                let first_absolute = existing.clone();
                let second_absolute = ComponentAbsoluteTarget {
                    node,
                    source_id: source_id.clone(),
                    target,
                };
                CanonicalValueEdgeOutcome::Conflict {
                    proof_source_ids: self
                        .absolute_conflict_source_ids(&first_absolute, &second_absolute),
                    proof_group_ids: self.proof_group_ids(first_absolute.node, node),
                    first_absolute,
                    second_absolute,
                }
            };
        }
        self.absolute_target[root] = Some(ComponentAbsoluteTarget {
            node,
            source_id: source_id.clone(),
            target,
        });
        CanonicalValueEdgeOutcome::Independent
    }

    fn intern(&mut self, node: CanonicalValueNode) -> usize {
        if let Some(index) = self.node_index.get(&node) {
            return *index;
        }
        let index = self.parent.len();
        self.nodes.push(node.clone());
        self.parent.push(index);
        self.absolute_target.push(None);
        self.adjacency.push(Vec::new());
        self.node_index.insert(node, index);
        index
    }

    fn root(&mut self, index: usize) -> usize {
        let parent = self.parent[index];
        if parent == index {
            return index;
        }
        let root = self.root(parent);
        self.parent[index] = root;
        root
    }

    fn proof_source_ids(&self, start: usize, end: usize) -> Vec<SourceId> {
        let path = self.proof_path(start, end);
        path.windows(2)
            .map(|nodes| {
                self.adjacency[nodes[0]]
                    .iter()
                    .find(|(neighbor, _)| *neighbor == nodes[1])
                    .expect("the recovered proof path retains every edge source")
                    .1
                    .clone()
            })
            .collect()
    }

    fn absolute_conflict_source_ids(
        &self,
        first: &ComponentAbsoluteTarget,
        second: &ComponentAbsoluteTarget,
    ) -> Vec<SourceId> {
        let mut sources = vec![first.source_id.clone()];
        sources.extend(self.proof_source_ids(first.node, second.node));
        sources.push(second.source_id.clone());
        sources
    }

    fn equality_conflict(
        &self,
        left: usize,
        right: usize,
        closing_source: &SourceId,
        left_absolute: &ComponentAbsoluteTarget,
        right_absolute: &ComponentAbsoluteTarget,
    ) -> CanonicalValueEdgeOutcome {
        let mut proof_source_ids = vec![left_absolute.source_id.clone()];
        proof_source_ids.extend(self.proof_source_ids(left_absolute.node, left));
        proof_source_ids.push(closing_source.clone());
        proof_source_ids.extend(self.proof_source_ids(right, right_absolute.node));
        proof_source_ids.push(right_absolute.source_id.clone());

        let mut proof_nodes = self.proof_path(left_absolute.node, left);
        proof_nodes.extend(self.proof_path(right, right_absolute.node));
        CanonicalValueEdgeOutcome::Conflict {
            proof_source_ids,
            proof_group_ids: self.group_ids(proof_nodes),
            first_absolute: left_absolute.clone(),
            second_absolute: right_absolute.clone(),
        }
    }

    fn proof_group_ids(&self, start: usize, end: usize) -> Vec<GroupId> {
        self.group_ids(self.proof_path(start, end))
    }

    fn group_ids(&self, nodes: Vec<usize>) -> Vec<GroupId> {
        let mut groups = nodes
            .into_iter()
            .filter_map(|index| match &self.nodes[index] {
                CanonicalValueNode::Group(group_id) => Some(group_id.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();
        groups.sort();
        groups.dedup();
        groups
    }

    fn proof_path(&self, start: usize, end: usize) -> Vec<usize> {
        let mut predecessor = vec![None::<usize>; self.adjacency.len()];
        let mut queue = VecDeque::from([start]);
        predecessor[start] = Some(start);
        while let Some(node) = queue.pop_front() {
            if node == end {
                break;
            }
            for (neighbor, _) in &self.adjacency[node] {
                if predecessor[*neighbor].is_none() {
                    predecessor[*neighbor] = Some(node);
                    queue.push_back(*neighbor);
                }
            }
        }
        let mut path = vec![end];
        let mut cursor = end;
        while cursor != start {
            cursor = predecessor[cursor].expect("equal union-find roots imply an adjacency path");
            path.push(cursor);
        }
        path.reverse();
        path
    }
}

#[derive(Debug, Clone, PartialEq)]
enum CanonicalValueEdgeOutcome {
    Independent,
    Redundant,
    Conflict {
        proof_source_ids: Vec<SourceId>,
        proof_group_ids: Vec<GroupId>,
        first_absolute: ComponentAbsoluteTarget,
        second_absolute: ComponentAbsoluteTarget,
    },
}

fn canonical_support_bits(support: [f64; 3]) -> [u64; 3] {
    support.map(canonical_scalar_bits)
}

fn canonical_scalar_bits(value: f64) -> u64 {
    if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    }
}

fn normalized_equality_key(equality: &CanonicalHardEquality) -> (CanonicalEqualityKey, f64) {
    let sign = canonical_normalization_sign(equality.field(), equality.latent_coefficients());
    let dimension = equality.dimension();
    let field_terms = equality
        .field()
        .into_iter()
        .flat_map(|field| field.functional().terms())
        .map(|term| {
            (
                term.support().map(f64::to_bits),
                signed_coefficient_bits(sign, term.value_coefficient()),
                term.gradient_coefficient()
                    .map(|coefficient| signed_coefficient_bits(sign, coefficient)),
            )
        })
        .collect();
    let mut latent_coefficients = equality
        .latent_coefficients()
        .iter()
        .map(|term| (term.latent, signed_coefficient_bits(sign, term.coefficient)))
        .collect::<Vec<_>>();
    latent_coefficients.sort_by_key(|(latent, _)| *latent);
    (
        CanonicalEqualityKey {
            dimension,
            field_terms,
            latent_coefficients,
        },
        sign * equality.target(),
    )
}

fn canonical_normalization_sign(
    field: Option<&FunctionalUse>,
    latent_coefficients: &[SemanticLatentCoefficient],
) -> f64 {
    let first_coefficient = field
        .into_iter()
        .flat_map(|field| field.functional().terms())
        .flat_map(|term| {
            std::iter::once(term.value_coefficient()).chain(term.gradient_coefficient())
        })
        .chain(latent_coefficients.iter().map(|term| term.coefficient))
        .find(|coefficient| *coefficient != 0.0)
        .expect("a canonical affine functional has a field or semantic-latent coefficient");
    if first_coefficient.is_sign_negative() {
        -1.0
    } else {
        1.0
    }
}

fn affine_bound_inequality(
    source_id: &SourceId,
    semantic_role: &'static str,
    functional: CanonicalFunctional,
    dimension: FunctionalDimension,
    sense: CanonicalInequalitySense,
    side: AffineBoundSide,
) -> CanonicalAffineInequality {
    affine_bound_inequality_with_representer_span(
        source_id,
        semantic_role,
        functional,
        dimension,
        sense,
        side,
        FunctionalRepresenterSpan::Native,
    )
}

fn affine_bound_inequality_with_representer_span(
    source_id: &SourceId,
    semantic_role: &'static str,
    functional: CanonicalFunctional,
    dimension: FunctionalDimension,
    sense: CanonicalInequalitySense,
    side: AffineBoundSide,
    representer_span: FunctionalRepresenterSpan,
) -> CanonicalAffineInequality {
    let provenance = relation_provenance(
        source_id.clone(),
        None,
        SemanticRolePath::new(semantic_role),
    );
    let violation_channel = violation_channel(side.configuration, &provenance);
    CanonicalAffineInequality::new(
        Some(FunctionalUse::with_representer_span(
            functional,
            provenance.clone(),
            representer_span,
        )),
        Vec::new(),
        provenance,
        dimension,
        sense,
        side.bound,
        violation_channel,
    )
}

fn violation_channel(
    configuration: AffineBoundConfiguration,
    provenance: &UsageProvenance,
) -> Option<CanonicalViolationChannel> {
    match configuration {
        AffineBoundConfiguration::Hard => None,
        AffineBoundConfiguration::QuadraticPenalty(penalty) => {
            Some(CanonicalViolationChannel::new(
                provenance.residual().clone(),
                CanonicalViolationLoss::QuadraticPenalty {
                    weight: penalty.weight(),
                },
            ))
        }
        AffineBoundConfiguration::LinearViolationPenalty(penalty) => {
            Some(CanonicalViolationChannel::new(
                provenance.residual().clone(),
                CanonicalViolationLoss::LinearViolationPenalty {
                    weight: penalty.weight(),
                },
            ))
        }
    }
}

#[derive(Debug, Clone)]
struct HardSharedLevelEdge {
    source_id: SourceId,
    from: GroupId,
    to: GroupId,
    semantic_role: SemanticRolePath,
    minimum_difference: f64,
}

fn hard_shared_level_edge(
    relation: &SharedLevelSetRelationInput,
    direction: Option<StratigraphicFieldDirection>,
) -> Option<HardSharedLevelEdge> {
    if !matches!(relation.configuration(), AffineBoundConfiguration::Hard) {
        return None;
    }
    let orientation = relation.orientation(direction);
    let minimum_difference = orientation.required_difference();
    Some(HardSharedLevelEdge {
        source_id: relation.source_id().clone(),
        from: orientation.lower_field_group_id,
        to: orientation.upper_field_group_id,
        semantic_role: relation.semantic_role(),
        minimum_difference,
    })
}

fn shared_level_set_relation_conflicts(
    relations: &[SharedLevelSetRelationInput],
    direction: Option<StratigraphicFieldDirection>,
) -> Vec<SharedLevelSetRelationConflictEvidence> {
    let mut edges = relations
        .iter()
        .filter_map(|relation| hard_shared_level_edge(relation, direction))
        .collect::<Vec<_>>();
    edges.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let mut adjacency = BTreeMap::<GroupId, Vec<usize>>::new();
    for (index, edge) in edges.iter().enumerate() {
        adjacency.entry(edge.from.clone()).or_default().push(index);
    }

    let mut seen = std::collections::BTreeSet::<Vec<SourceId>>::new();
    let mut conflicts = Vec::new();
    for (edge_index, edge) in edges.iter().enumerate() {
        let proof = if edge.from == edge.to {
            Some(vec![edge_index])
        } else if edge.minimum_difference > 0.0 {
            shared_level_path(&edges, &adjacency, &edge.to, &edge.from).map(|mut path| {
                path.push(edge_index);
                path
            })
        } else {
            None
        };
        let Some(proof) = proof else {
            continue;
        };
        let source_provenance = proof
            .iter()
            .map(|index| {
                let edge = &edges[*index];
                SharedLevelSetConflictSourceEvidence::new(
                    edge.source_id.clone(),
                    vec![edge.from.clone(), edge.to.clone()],
                    edge.semantic_role.clone(),
                )
            })
            .collect::<Vec<_>>();
        let mut source_ids = source_provenance
            .iter()
            .map(|source| source.source_id().clone())
            .collect::<Vec<_>>();
        source_ids.sort();
        source_ids.dedup();
        if !seen.insert(source_ids.clone()) {
            continue;
        }
        conflicts.push(SharedLevelSetRelationConflictEvidence::new(
            source_provenance,
        ));
    }
    conflicts.sort_by(|left, right| {
        left.source_ids()
            .cmp(right.source_ids())
            .then_with(|| left.group_ids().cmp(right.group_ids()))
    });
    conflicts
}

fn absolute_shared_level_set_relation_conflicts(
    snapshot: &ProblemSnapshot,
    value_constraints: &mut CanonicalValueConstraintForest,
) -> Vec<SharedLevelSetRelationConflictEvidence> {
    let mut edges = snapshot
        .inner
        .shared_level_set_relations
        .iter()
        .filter_map(|relation| {
            hard_shared_level_edge(relation, snapshot.inner.stratigraphic_field_direction)
        })
        .collect::<Vec<_>>();
    edges.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let provenance_by_source = shared_level_set_conflict_provenance(snapshot, &edges);
    let mut absolutes = BTreeMap::<GroupId, AbsoluteGroupProof>::new();
    for group in &snapshot.inner.shared_level_sets {
        if let Some(proof) = value_constraints.absolute_proof_for_group(group.group_id()) {
            absolutes.insert(group.group_id().clone(), proof);
        }
    }

    let mut seen = std::collections::BTreeSet::<Vec<SourceId>>::new();
    let mut conflicts = Vec::new();
    for (lower_group, lower_absolute) in &absolutes {
        let paths = maximum_shared_level_paths(&edges, lower_group);
        for (upper_group, (minimum_difference, edge_path)) in paths {
            if &upper_group == lower_group || edge_path.is_empty() {
                continue;
            }
            let Some(upper_absolute) = absolutes.get(&upper_group) else {
                continue;
            };
            if upper_absolute.target - lower_absolute.target >= minimum_difference {
                continue;
            }
            let mut source_ids = lower_absolute.source_ids.clone();
            for edge_index in edge_path {
                let edge = &edges[edge_index];
                source_ids.push(edge.source_id.clone());
            }
            source_ids.extend(upper_absolute.source_ids.iter().cloned());
            source_ids.sort();
            source_ids.dedup();
            if !seen.insert(source_ids.clone()) {
                continue;
            }
            let source_provenance = source_ids
                .iter()
                .map(|source_id| {
                    provenance_by_source
                        .get(source_id)
                        .cloned()
                        .expect("every absolute relation proof source retains original provenance")
                })
                .collect();
            conflicts.push(SharedLevelSetRelationConflictEvidence::new(
                source_provenance,
            ));
        }
    }
    conflicts.sort_by(|left, right| {
        left.source_ids()
            .cmp(right.source_ids())
            .then_with(|| left.group_ids().cmp(right.group_ids()))
    });
    conflicts
}

fn shared_level_set_conflict_provenance(
    snapshot: &ProblemSnapshot,
    edges: &[HardSharedLevelEdge],
) -> BTreeMap<SourceId, SharedLevelSetConflictSourceEvidence> {
    let mut provenance = BTreeMap::new();
    for block in observation_residual_blocks(&snapshot.inner.observations) {
        if !matches!(block.configuration, ResidualBlockConfiguration::Hard) {
            continue;
        }
        for observation in block.components {
            if observation.functional.dimension != FunctionalDimension::FieldValue {
                continue;
            }
            provenance.insert(
                observation.source_id.clone(),
                SharedLevelSetConflictSourceEvidence::new(
                    observation.source_id,
                    observation.group_id.into_iter().collect(),
                    observation.semantic_role,
                ),
            );
        }
    }
    for group in &snapshot.inner.shared_level_sets {
        for member in group.members() {
            provenance.insert(
                member.source_id().clone(),
                SharedLevelSetConflictSourceEvidence::new(
                    member.source_id().clone(),
                    vec![group.group_id().clone()],
                    SemanticRolePath::new("shared-level-set/member/value"),
                ),
            );
        }
    }
    for gauge in &snapshot.inner.additive_field_gauges {
        let (group_ids, semantic_role) = match gauge.reference() {
            AdditiveFieldGaugeReference::Point(_) => (
                Vec::new(),
                SemanticRolePath::new("additive-field-gauge/point"),
            ),
            AdditiveFieldGaugeReference::LevelSet(group_id) => (
                vec![group_id.clone()],
                SemanticRolePath::new("additive-field-gauge/level-set"),
            ),
        };
        provenance.insert(
            gauge.source_id().clone(),
            SharedLevelSetConflictSourceEvidence::new(
                gauge.source_id().clone(),
                group_ids,
                semantic_role,
            ),
        );
    }
    for edge in edges {
        provenance.insert(
            edge.source_id.clone(),
            SharedLevelSetConflictSourceEvidence::new(
                edge.source_id.clone(),
                vec![edge.from.clone(), edge.to.clone()],
                edge.semantic_role.clone(),
            ),
        );
    }
    provenance
}

fn preflight_shared_level_set_relation_conflicts(
    snapshot: &ProblemSnapshot,
) -> Vec<SharedLevelSetRelationConflictEvidence> {
    let mut conflicts = shared_level_set_relation_conflicts(
        &snapshot.inner.shared_level_set_relations,
        snapshot.inner.stratigraphic_field_direction,
    );
    let mut value_constraints = CanonicalValueConstraintForest::default();
    for block in observation_residual_blocks(&snapshot.inner.observations) {
        if !matches!(block.configuration, ResidualBlockConfiguration::Hard) {
            continue;
        }
        for observation in block.components {
            if observation.functional.dimension == FunctionalDimension::FieldValue {
                value_constraints.add_absolute_support(
                    observation.support,
                    observation.target,
                    &observation.source_id,
                );
            }
        }
    }
    for group in &snapshot.inner.shared_level_sets {
        for member in group.members() {
            value_constraints.add_member_equality(
                group.group_id(),
                member.location().components(),
                member.source_id(),
            );
        }
    }
    for gauge in &snapshot.inner.additive_field_gauges {
        match gauge.reference() {
            AdditiveFieldGaugeReference::Point(point) => {
                value_constraints.add_absolute_support(
                    point.components(),
                    gauge.value(),
                    gauge.source_id(),
                );
            }
            AdditiveFieldGaugeReference::LevelSet(group_id) => {
                value_constraints.add_absolute_group(group_id, gauge.value(), gauge.source_id());
            }
        }
    }
    conflicts.extend(absolute_shared_level_set_relation_conflicts(
        snapshot,
        &mut value_constraints,
    ));
    conflicts.sort_by(|left, right| {
        left.source_ids()
            .cmp(right.source_ids())
            .then_with(|| left.group_ids().cmp(right.group_ids()))
    });
    conflicts.dedup_by(|left, right| left.source_ids() == right.source_ids());
    conflicts
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum AffineValueNode {
    Anchor,
    Group(GroupId),
    Support([u64; 3]),
}

#[derive(Debug, Clone)]
struct HardAffineValueEdge {
    provenance: UsageProvenance,
    from: AffineValueNode,
    to: AffineValueNode,
    minimum_difference: ExactDyadic,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// An exact signed integer multiple of `2^-1074`, the finest binary64 unit.
struct ExactDyadic {
    negative: bool,
    magnitude: Vec<u64>,
}

impl ExactDyadic {
    fn zero() -> Self {
        Self {
            negative: false,
            magnitude: Vec::new(),
        }
    }

    fn from_f64(value: f64) -> Self {
        debug_assert!(value.is_finite());
        if value == 0.0 {
            return Self::zero();
        }
        let bits = value.to_bits();
        let negative = bits >> 63 != 0;
        let raw_exponent = ((bits >> 52) & 0x7ff) as usize;
        let fraction = bits & ((1_u64 << 52) - 1);
        let (significand, shift) = if raw_exponent == 0 {
            (fraction, 0)
        } else {
            ((1_u64 << 52) | fraction, raw_exponent - 1)
        };
        let limb = shift / 64;
        let bit = shift % 64;
        let mut magnitude = vec![0_u64; limb + 2];
        magnitude[limb] = significand << bit;
        if bit != 0 {
            magnitude[limb + 1] = significand >> (64 - bit);
        }
        Self::normalized(negative, magnitude)
    }

    fn negated(&self) -> Self {
        if self.magnitude.is_empty() {
            return Self::zero();
        }
        Self {
            negative: !self.negative,
            magnitude: self.magnitude.clone(),
        }
    }

    fn add(&self, other: &Self) -> Self {
        if self.negative == other.negative {
            return Self::normalized(
                self.negative,
                add_magnitudes(&self.magnitude, &other.magnitude),
            );
        }
        match compare_magnitudes(&self.magnitude, &other.magnitude) {
            std::cmp::Ordering::Greater => Self::normalized(
                self.negative,
                subtract_magnitudes(&self.magnitude, &other.magnitude),
            ),
            std::cmp::Ordering::Less => Self::normalized(
                other.negative,
                subtract_magnitudes(&other.magnitude, &self.magnitude),
            ),
            std::cmp::Ordering::Equal => Self::zero(),
        }
    }

    fn is_greater_than(&self, other: &Self) -> bool {
        if self.negative != other.negative {
            return other.negative;
        }
        match compare_magnitudes(&self.magnitude, &other.magnitude) {
            std::cmp::Ordering::Greater => !self.negative,
            std::cmp::Ordering::Less => self.negative,
            std::cmp::Ordering::Equal => false,
        }
    }

    fn is_positive(&self) -> bool {
        !self.negative && !self.magnitude.is_empty()
    }

    fn normalized(negative: bool, mut magnitude: Vec<u64>) -> Self {
        while magnitude.last() == Some(&0) {
            magnitude.pop();
        }
        Self {
            negative: negative && !magnitude.is_empty(),
            magnitude,
        }
    }
}

fn compare_magnitudes(left: &[u64], right: &[u64]) -> std::cmp::Ordering {
    left.len()
        .cmp(&right.len())
        .then_with(|| left.iter().rev().cmp(right.iter().rev()))
}

fn add_magnitudes(left: &[u64], right: &[u64]) -> Vec<u64> {
    let mut sum = Vec::with_capacity(left.len().max(right.len()) + 1);
    let mut carry = 0_u128;
    for index in 0..left.len().max(right.len()) {
        let total = u128::from(left.get(index).copied().unwrap_or(0))
            + u128::from(right.get(index).copied().unwrap_or(0))
            + carry;
        sum.push(total as u64);
        carry = total >> 64;
    }
    if carry != 0 {
        sum.push(carry as u64);
    }
    sum
}

fn subtract_magnitudes(larger: &[u64], smaller: &[u64]) -> Vec<u64> {
    debug_assert!(compare_magnitudes(larger, smaller) != std::cmp::Ordering::Less);
    let mut difference = Vec::with_capacity(larger.len());
    let mut borrow = false;
    for (index, larger_limb) in larger.iter().copied().enumerate() {
        let (without_smaller, first_borrow) =
            larger_limb.overflowing_sub(smaller.get(index).copied().unwrap_or(0));
        let (limb, second_borrow) = without_smaller.overflowing_sub(u64::from(borrow));
        difference.push(limb);
        borrow = first_borrow || second_borrow;
    }
    debug_assert!(!borrow);
    difference
}

fn canonical_affine_value_cycle_conflict(
    snapshot: &ProblemSnapshot,
    lowering: &EqualityLowering,
) -> Option<SharedLevelSetRelationConflictEvidence> {
    if snapshot.inner.field_separation_intervals.is_empty()
        && snapshot.inner.point_to_level_set_relations.is_empty()
    {
        return None;
    }
    let issue_33_sources = snapshot
        .inner
        .field_separation_intervals
        .iter()
        .map(|relation| relation.source_id().clone())
        .chain(
            snapshot
                .inner
                .point_to_level_set_relations
                .iter()
                .map(|relation| relation.source_id().clone()),
        )
        .collect::<std::collections::BTreeSet<_>>();
    let mut edges = lowering
        .canonical_equalities
        .iter()
        .filter_map(|equality| canonical_exact_value_edge(equality, lowering))
        .flat_map(|edge| {
            let reverse = HardAffineValueEdge {
                provenance: edge.provenance.clone(),
                from: edge.to.clone(),
                to: edge.from.clone(),
                minimum_difference: edge.minimum_difference.negated(),
            };
            [edge, reverse]
        })
        .chain(
            lowering
                .source_bound_relations
                .iter()
                .filter_map(|relation| canonical_affine_value_edge(&relation.inequality, lowering)),
        )
        .collect::<Vec<_>>();

    edges.sort_by(|left, right| {
        left.provenance
            .source()
            .cmp(right.provenance.source())
            .then_with(|| {
                left.provenance
                    .semantic_role()
                    .cmp(right.provenance.semantic_role())
            })
            .then_with(|| left.from.cmp(&right.from))
            .then_with(|| left.to.cmp(&right.to))
    });
    let nodes = edges
        .iter()
        .flat_map(|edge| [edge.from.clone(), edge.to.clone()])
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let proof_edges = positive_canonical_affine_cycle(&edges, &nodes, &issue_33_sources)?;
    let source_provenance = proof_edges
        .into_iter()
        .map(|edge_index| {
            let edge = &edges[edge_index];
            SharedLevelSetConflictSourceEvidence::new(
                edge.provenance.source().clone(),
                edge.provenance.groups().to_vec(),
                edge.provenance.semantic_role().clone(),
            )
        })
        .collect();
    Some(SharedLevelSetRelationConflictEvidence::new(
        source_provenance,
    ))
}

fn canonical_exact_value_edge(
    equality: &CanonicalHardEquality,
    lowering: &EqualityLowering,
) -> Option<HardAffineValueEdge> {
    if equality.dimension() != FunctionalDimension::FieldValue {
        return None;
    }
    let (from, to) = canonical_unit_difference_nodes(
        equality.field(),
        equality.latent_coefficients(),
        1.0,
        lowering,
    )?;
    Some(HardAffineValueEdge {
        provenance: equality.provenance().clone(),
        from,
        to,
        minimum_difference: ExactDyadic::from_f64(equality.target()),
    })
}

fn canonical_affine_value_edge(
    inequality: &CanonicalAffineInequality,
    lowering: &EqualityLowering,
) -> Option<HardAffineValueEdge> {
    if inequality.dimension() != FunctionalDimension::FieldValue
        || inequality.violation_channel().is_some()
    {
        return None;
    }
    let multiplier = match inequality.sense() {
        CanonicalInequalitySense::Lower => 1.0,
        CanonicalInequalitySense::Upper => -1.0,
    };
    let (from, to) = canonical_unit_difference_nodes(
        inequality.field(),
        inequality.latent_coefficients(),
        multiplier,
        lowering,
    )?;
    Some(HardAffineValueEdge {
        provenance: inequality.provenance().clone(),
        from,
        to,
        minimum_difference: ExactDyadic::from_f64(multiplier * inequality.bound()),
    })
}

fn canonical_unit_difference_nodes(
    field: Option<&FunctionalUse>,
    latent_coefficients: &[SemanticLatentCoefficient],
    multiplier: f64,
    lowering: &EqualityLowering,
) -> Option<(AffineValueNode, AffineValueNode)> {
    let mut coefficients = Vec::new();
    for term in field.into_iter().flat_map(|use_| use_.functional().terms()) {
        if term.gradient_coefficient() != [0.0; 3] {
            return None;
        }
        coefficients.push((
            AffineValueNode::Support(canonical_support_bits(term.support())),
            multiplier * term.value_coefficient(),
        ));
    }
    for term in latent_coefficients {
        coefficients.push((
            AffineValueNode::Group(lowering.semantic_latents.get(term.latent)?.group_id.clone()),
            multiplier * term.coefficient,
        ));
    }
    coefficients.sort_by(|left, right| left.0.cmp(&right.0));
    match coefficients.as_slice() {
        [(node, 1.0)] => Some((AffineValueNode::Anchor, node.clone())),
        [(node, -1.0)] => Some((node.clone(), AffineValueNode::Anchor)),
        [(left, -1.0), (right, 1.0)] => Some((left.clone(), right.clone())),
        [(left, 1.0), (right, -1.0)] => Some((right.clone(), left.clone())),
        _ => None,
    }
}

fn positive_canonical_affine_cycle(
    edges: &[HardAffineValueEdge],
    nodes: &[AffineValueNode],
    issue_33_sources: &std::collections::BTreeSet<SourceId>,
) -> Option<Vec<usize>> {
    let node_index = nodes
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, node)| (node, index))
        .collect::<BTreeMap<_, _>>();
    let mut adjacency = vec![Vec::new(); nodes.len()];
    let mut reverse_adjacency = vec![Vec::new(); nodes.len()];
    for edge in edges {
        let from = node_index[&edge.from];
        let to = node_index[&edge.to];
        adjacency[from].push(to);
        reverse_adjacency[to].push(from);
    }
    let components = strongly_connected_components(&adjacency, &reverse_adjacency);
    let selected_components = edges
        .iter()
        .filter(|edge| issue_33_sources.contains(edge.provenance.source()))
        .filter_map(|edge| {
            let from_component = components[node_index[&edge.from]];
            (from_component == components[node_index[&edge.to]]).then_some(from_component)
        })
        .collect::<std::collections::BTreeSet<_>>();
    let proof_edge_indices = edges
        .iter()
        .enumerate()
        .filter_map(|(edge_index, edge)| {
            let from_component = components[node_index[&edge.from]];
            (from_component == components[node_index[&edge.to]]
                && selected_components.contains(&from_component))
            .then_some(edge_index)
        })
        .collect::<Vec<_>>();
    if proof_edge_indices.is_empty() {
        return None;
    }

    let mut distances = vec![ExactDyadic::zero(); nodes.len()];
    let mut predecessor = vec![None::<usize>; nodes.len()];
    let mut updated = None;
    for _ in 0..nodes.len() {
        updated = None;
        for edge_index in &proof_edge_indices {
            let edge = &edges[*edge_index];
            let source = node_index[&edge.from];
            let destination = node_index[&edge.to];
            let candidate = distances[source].add(&edge.minimum_difference);
            if !candidate.is_greater_than(&distances[destination]) {
                continue;
            }
            distances[destination] = candidate;
            predecessor[destination] = Some(*edge_index);
            updated = Some(destination);
        }
    }
    let mut cycle_node = updated?;
    for _ in 0..nodes.len() {
        let edge_index = predecessor[cycle_node]?;
        cycle_node = node_index[&edges[edge_index].from];
    }
    let start = cycle_node;
    let mut cycle = Vec::new();
    loop {
        let edge_index = predecessor[cycle_node]?;
        cycle.push(edge_index);
        cycle_node = node_index[&edges[edge_index].from];
        if cycle_node == start {
            break;
        }
        if cycle.len() > nodes.len() {
            return None;
        }
    }
    let cycle_weight = cycle
        .iter()
        .fold(ExactDyadic::zero(), |weight, edge_index| {
            weight.add(&edges[*edge_index].minimum_difference)
        });
    cycle_weight.is_positive().then_some(cycle)
}

fn strongly_connected_components(adjacency: &[Vec<usize>], reverse: &[Vec<usize>]) -> Vec<usize> {
    let mut visited = vec![false; adjacency.len()];
    let mut finish_order = Vec::with_capacity(adjacency.len());
    for start in 0..adjacency.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut stack = vec![(start, 0_usize)];
        while let Some((node, next_neighbor)) = stack.last_mut() {
            if let Some(neighbor) = adjacency[*node].get(*next_neighbor).copied() {
                *next_neighbor += 1;
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push((neighbor, 0));
                }
                continue;
            }
            let (finished, _) = stack.pop().expect("the DFS stack has one current node");
            finish_order.push(finished);
        }
    }

    let mut components = vec![usize::MAX; adjacency.len()];
    let mut component = 0_usize;
    for start in finish_order.into_iter().rev() {
        if components[start] != usize::MAX {
            continue;
        }
        components[start] = component;
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            for neighbor in &reverse[node] {
                if components[*neighbor] == usize::MAX {
                    components[*neighbor] = component;
                    stack.push(*neighbor);
                }
            }
        }
        component += 1;
    }
    components
}

fn maximum_shared_level_paths(
    edges: &[HardSharedLevelEdge],
    start: &GroupId,
) -> BTreeMap<GroupId, (f64, Vec<usize>)> {
    let mut paths = BTreeMap::from([(start.clone(), (0.0, Vec::new()))]);
    for _ in 0..edges.len() {
        let mut changed = false;
        for (edge_index, edge) in edges.iter().enumerate() {
            let Some((distance, path)) = paths.get(&edge.from).cloned() else {
                continue;
            };
            let candidate = distance + edge.minimum_difference;
            let should_update = paths
                .get(&edge.to)
                .is_none_or(|(known, _)| candidate > *known);
            if should_update {
                let mut candidate_path = path;
                candidate_path.push(edge_index);
                paths.insert(edge.to.clone(), (candidate, candidate_path));
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    paths
}

fn shared_level_path(
    edges: &[HardSharedLevelEdge],
    adjacency: &BTreeMap<GroupId, Vec<usize>>,
    start: &GroupId,
    end: &GroupId,
) -> Option<Vec<usize>> {
    let mut predecessor = BTreeMap::<GroupId, (GroupId, usize)>::new();
    let mut queue = VecDeque::from([start.clone()]);
    let mut visited = std::collections::BTreeSet::from([start.clone()]);
    while let Some(node) = queue.pop_front() {
        if &node == end {
            break;
        }
        for edge_index in adjacency.get(&node).into_iter().flatten() {
            let next = edges[*edge_index].to.clone();
            if visited.insert(next.clone()) {
                predecessor.insert(next.clone(), (node.clone(), *edge_index));
                queue.push_back(next);
            }
        }
    }
    if !visited.contains(end) {
        return None;
    }
    let mut path = Vec::new();
    let mut cursor = end.clone();
    while &cursor != start {
        let (previous, edge_index) = predecessor.get(&cursor)?.clone();
        path.push(edge_index);
        cursor = previous;
    }
    path.reverse();
    Some(path)
}

fn lower_snapshot(snapshot: &ProblemSnapshot) -> EqualityLowering {
    let mut lowering = EqualityLowering::new(resolved_normal_inputs(snapshot));
    let mut value_constraints = CanonicalValueConstraintForest::default();
    let mut hard_bound_facts = Vec::new();
    let mut hard_derivative_bound_facts = Vec::new();
    for block in observation_residual_blocks(&snapshot.inner.observations) {
        let soft_loss = match &block.configuration {
            ResidualBlockConfiguration::Hard => None,
            ResidualBlockConfiguration::QuadraticPenalty(penalty) => {
                Some(CanonicalSoftLoss::QuadraticPenalty {
                    weight: penalty.weight(),
                })
            }
            ResidualBlockConfiguration::StandardDeviation(standard_deviation) => {
                Some(CanonicalSoftLoss::StandardDeviation {
                    standard_deviation: standard_deviation.value(),
                })
            }
            ResidualBlockConfiguration::Covariance(covariance) => {
                Some(CanonicalSoftLoss::covariance(
                    covariance.dimension(),
                    covariance.entries().to_vec(),
                ))
            }
        };
        if let Some(loss) = soft_loss {
            lowering.push_soft_block(
                block.components.iter().map(soft_field_equality).collect(),
                loss,
                None,
                block.kind,
            );
            continue;
        }
        for observation in block.components {
            let edge =
                (observation.functional.dimension == FunctionalDimension::FieldValue).then(|| {
                    value_constraints.add_absolute_support(
                        observation.support,
                        observation.target,
                        &observation.source_id,
                    )
                });
            if let Some(edge) = &edge {
                lowering.record_graph_conflict(edge, &observation.semantic_role);
            }
            let participation = if edge
                .as_ref()
                .is_some_and(|edge| !matches!(edge, CanonicalValueEdgeOutcome::Independent))
            {
                CanonicalEqualityParticipation::VerificationOnly
            } else {
                CanonicalEqualityParticipation::SolverConstraint
            };
            lowering.push_source(field_equality(&observation, participation));
        }
    }

    for (normal_index, normal) in lowering.resolved_normals.clone().iter().enumerate() {
        let projection_components = normal_projection_components(normal);
        match normal.direction_enforcement.configuration() {
            NormalDirectionConfiguration::Hard => {
                let dropped_axis = (0..3)
                    .max_by(|left, right| {
                        normal.direction.components()[*left]
                            .abs()
                            .total_cmp(&normal.direction.components()[*right].abs())
                    })
                    .expect("a physical normal has three components");
                let drop_one = projection_components.len() == 3;
                let mut canonical_indices = Vec::with_capacity(projection_components.len());
                for (axis, component) in &projection_components {
                    let participation = if drop_one && *axis == dropped_axis {
                        CanonicalEqualityParticipation::VerificationOnly
                    } else {
                        CanonicalEqualityParticipation::SolverConstraint
                    };
                    let equality = field_equality_with_representer_span(
                        component,
                        participation,
                        FunctionalRepresenterSpan::CompleteGradientAtSupport,
                    );
                    canonical_indices.push(lowering.push_source_with_kind(
                        equality,
                        SourceHardRelationKind::NormalProjection { normal_index },
                    ));
                }
                lowering.canonical_hard_residual_blocks.push(
                    CanonicalHardResidualBlock::normal_projection(canonical_indices),
                );
            }
            NormalDirectionConfiguration::QuadraticPenalty(penalty) => {
                lowering.push_soft_block(
                    projection_components
                        .iter()
                        .map(|(_, component)| {
                            soft_field_equality_with_representer_span(
                                component,
                                FunctionalRepresenterSpan::CompleteGradientAtSupport,
                            )
                        })
                        .collect(),
                    CanonicalSoftLoss::QuadraticPenalty {
                        weight: penalty.weight(),
                    },
                    None,
                    CanonicalSoftResidualBlockKind::NormalProjection,
                );
            }
            NormalDirectionConfiguration::StandardDeviation(standard_deviation) => {
                lowering.push_soft_block(
                    projection_components
                        .iter()
                        .map(|(_, component)| {
                            soft_field_equality_with_representer_span(
                                component,
                                FunctionalRepresenterSpan::CompleteGradientAtSupport,
                            )
                        })
                        .collect(),
                    CanonicalSoftLoss::StandardDeviation {
                        standard_deviation: standard_deviation.value(),
                    },
                    None,
                    CanonicalSoftResidualBlockKind::NormalProjection,
                );
            }
        }

        let configuration = match normal.minimum_slope_enforcement.configuration() {
            MinimumNormalSlopeConfiguration::Hard => AffineBoundConfiguration::Hard,
            MinimumNormalSlopeConfiguration::QuadraticPenalty(penalty) => {
                AffineBoundConfiguration::QuadraticPenalty(penalty)
            }
            MinimumNormalSlopeConfiguration::LinearViolationPenalty(penalty) => {
                AffineBoundConfiguration::LinearViolationPenalty(penalty)
            }
        };
        let functional = CanonicalFunctional::new(
            FunctionalDimension::FieldValuePerLength,
            vec![FunctionalTerm::new(
                normal.location.components(),
                0.0,
                normal.direction.components(),
            )],
        )
        .expect("a checked normal slope lowers to one finite derivative functional");
        let inequality = affine_bound_inequality_with_representer_span(
            &normal.source_id,
            "directed-normal/minimum-slope",
            functional,
            FunctionalDimension::FieldValuePerLength,
            CanonicalInequalitySense::Lower,
            AffineBoundSide {
                bound: normal.minimum_slope.value(),
                configuration,
            },
            FunctionalRepresenterSpan::CompleteGradientAtSupport,
        );
        if matches!(configuration, AffineBoundConfiguration::Hard) {
            let (functional_key, sign) = canonical_inequality_functional_key(&inequality);
            hard_derivative_bound_facts.push(HardDerivativeBoundFact {
                source_id: normal.source_id.clone(),
                functional_key,
                sense: if sign < 0.0 {
                    CanonicalInequalitySense::Upper
                } else {
                    CanonicalInequalitySense::Lower
                },
                bound: sign * normal.minimum_slope.value(),
                semantic_role: inequality.provenance().semantic_role().clone(),
            });
        }
        lowering.push_bound(inequality, SourceBoundKind::NormalSlope { normal_index });
    }

    for group in &snapshot.inner.covariance_groups {
        let components = covariance_group_components(group.group_id(), group.members());
        lowering.push_soft_block(
            components.iter().map(soft_field_equality).collect(),
            CanonicalSoftLoss::covariance(
                group.covariance().dimension(),
                group.covariance().entries().to_vec(),
            ),
            Some(group.group_id().clone()),
            CanonicalSoftResidualBlockKind::CovarianceGroup {
                members: covariance_group_member_kinds(group.members()),
            },
        );
    }

    for bound in &snapshot.inner.field_value_bounds {
        for (sense, side, role) in [
            (
                CanonicalInequalitySense::Lower,
                bound.lower(),
                "field-value-bound/lower",
            ),
            (
                CanonicalInequalitySense::Upper,
                bound.upper(),
                "field-value-bound/upper",
            ),
        ] {
            let Some(side) = side else {
                continue;
            };
            let functional = CanonicalFunctional::new(
                FunctionalDimension::FieldValue,
                vec![FunctionalTerm::new(
                    bound.location().components(),
                    1.0,
                    [0.0; 3],
                )],
            )
            .expect("a checked Field Value Bound lowers to one finite value functional");
            let inequality = affine_bound_inequality(
                bound.source_id(),
                role,
                functional,
                FunctionalDimension::FieldValue,
                sense,
                *side,
            );
            if matches!(side.configuration, AffineBoundConfiguration::Hard) {
                hard_bound_facts.push(HardBoundFact {
                    source_id: bound.source_id().clone(),
                    support: bound.location().components(),
                    sense,
                    bound: side.bound,
                    semantic_role: inequality.provenance().semantic_role().clone(),
                });
            }
            lowering.push_bound(inequality, SourceBoundKind::FieldValue);
        }
    }

    for interval in &snapshot.inner.directional_derivative_intervals {
        for (sense, side, role) in [
            (
                CanonicalInequalitySense::Lower,
                interval.lower(),
                "directional-derivative-interval/lower",
            ),
            (
                CanonicalInequalitySense::Upper,
                interval.upper(),
                "directional-derivative-interval/upper",
            ),
        ] {
            let Some(side) = side else {
                continue;
            };
            let functional = CanonicalFunctional::new(
                FunctionalDimension::FieldValuePerLength,
                vec![FunctionalTerm::new(
                    interval.location().components(),
                    0.0,
                    interval.direction().components(),
                )],
            )
            .expect("a checked Directional Derivative Interval lowers to one finite functional");
            let inequality = affine_bound_inequality(
                interval.source_id(),
                role,
                functional,
                FunctionalDimension::FieldValuePerLength,
                sense,
                *side,
            );
            if matches!(side.configuration, AffineBoundConfiguration::Hard) {
                let (functional_key, sign) = canonical_inequality_functional_key(&inequality);
                hard_derivative_bound_facts.push(HardDerivativeBoundFact {
                    source_id: interval.source_id().clone(),
                    functional_key,
                    sense: if sign < 0.0 {
                        match sense {
                            CanonicalInequalitySense::Lower => CanonicalInequalitySense::Upper,
                            CanonicalInequalitySense::Upper => CanonicalInequalitySense::Lower,
                        }
                    } else {
                        sense
                    },
                    bound: sign * side.bound,
                    semantic_role: inequality.provenance().semantic_role().clone(),
                });
            }
            lowering.push_bound(
                inequality,
                SourceBoundKind::DirectionalDerivative {
                    direction: interval.direction(),
                },
            );
        }
    }

    let latent_index_by_group = snapshot
        .inner
        .shared_level_sets
        .iter()
        .enumerate()
        .map(|(index, group)| (group.group_id().clone(), index))
        .collect::<BTreeMap<_, _>>();
    for group in &snapshot.inner.shared_level_sets {
        let latent = lowering.semantic_latents.len();
        debug_assert_eq!(latent_index_by_group[group.group_id()], latent);
        lowering.semantic_latents.push(SemanticLatentDefinition {
            group_id: group.group_id().clone(),
            field_unit: snapshot.inner.field_unit.clone(),
            member_source_ids: group
                .members()
                .iter()
                .map(|member| member.source_id().clone())
                .collect(),
        });
        for member in group.members() {
            let role = SemanticRolePath::new("shared-level-set/member/value");
            let edge = value_constraints.add_member_equality(
                group.group_id(),
                member.location().components(),
                member.source_id(),
            );
            lowering.record_graph_conflict(&edge, &role);
            let provenance = relation_provenance(
                member.source_id().clone(),
                Some(group.group_id().clone()),
                role,
            );
            let field = CanonicalFunctional::new(
                FunctionalDimension::FieldValue,
                vec![FunctionalTerm::new(
                    member.location().components(),
                    1.0,
                    [0.0; 3],
                )],
            )
            .expect("a shared-level member lowers to one finite value functional");
            lowering.push_source(CanonicalHardEquality::new(
                Some(FunctionalUse::new(field, provenance.clone())),
                vec![SemanticLatentCoefficient {
                    latent,
                    coefficient: -1.0,
                }],
                provenance,
                FunctionalDimension::FieldValue,
                0.0,
                if matches!(edge, CanonicalValueEdgeOutcome::Independent) {
                    CanonicalEqualityParticipation::SolverConstraint
                } else {
                    CanonicalEqualityParticipation::VerificationOnly
                },
            ));
        }
    }

    for relation in &snapshot.inner.shared_level_set_relations {
        let orientation = relation.orientation(snapshot.inner.stratigraphic_field_direction);
        let lower_latent = latent_index_by_group[&orientation.lower_field_group_id];
        let upper_latent = latent_index_by_group[&orientation.upper_field_group_id];
        let semantic_role = relation.semantic_role();
        let declared_group_ids = relation.declared_group_ids();
        let provenance = relation_provenance_for_groups(
            relation.source_id().clone(),
            declared_group_ids.into_iter().cloned().collect(),
            semantic_role,
        );
        let inequality = CanonicalAffineInequality::new(
            None,
            vec![
                SemanticLatentCoefficient {
                    latent: lower_latent,
                    coefficient: -1.0,
                },
                SemanticLatentCoefficient {
                    latent: upper_latent,
                    coefficient: 1.0,
                },
            ],
            provenance.clone(),
            FunctionalDimension::FieldValue,
            CanonicalInequalitySense::Lower,
            orientation.required_difference(),
            violation_channel(relation.configuration(), &provenance),
        );
        lowering.push_bound(
            inequality,
            SourceBoundKind::SharedLevelSetRelation {
                relation: relation.clone(),
                orientation,
            },
        );
    }

    for interval in &snapshot.inner.field_separation_intervals {
        let reference_latent = latent_index_by_group[interval.reference_group_id()];
        let target_latent = latent_index_by_group[interval.target_group_id()];
        for (sense, side, role) in [
            (
                CanonicalInequalitySense::Lower,
                interval.lower(),
                "field-separation-interval/lower",
            ),
            (
                CanonicalInequalitySense::Upper,
                interval.upper(),
                "field-separation-interval/upper",
            ),
        ] {
            let provenance = relation_provenance_for_groups(
                interval.source_id().clone(),
                vec![
                    interval.reference_group_id().clone(),
                    interval.target_group_id().clone(),
                ],
                SemanticRolePath::new(role),
            );
            let inequality = CanonicalAffineInequality::new(
                None,
                vec![
                    SemanticLatentCoefficient {
                        latent: reference_latent,
                        coefficient: -1.0,
                    },
                    SemanticLatentCoefficient {
                        latent: target_latent,
                        coefficient: 1.0,
                    },
                ],
                provenance.clone(),
                FunctionalDimension::FieldValue,
                sense,
                side.bound,
                violation_channel(side.configuration, &provenance),
            );
            lowering.push_bound(
                inequality,
                SourceBoundKind::FieldSeparationInterval {
                    relation: interval.clone(),
                },
            );
        }
    }

    for relation in &snapshot.inner.point_to_level_set_relations {
        let latent = latent_index_by_group[relation.group_id()];
        let orientation = relation.orientation();
        let sense = if orientation.is_lower_bounded() {
            CanonicalInequalitySense::Lower
        } else {
            CanonicalInequalitySense::Upper
        };
        let provenance = relation_provenance_for_groups(
            relation.source_id().clone(),
            vec![relation.group_id().clone()],
            orientation.semantic_role(),
        );
        let functional = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![FunctionalTerm::new(
                relation.location().components(),
                1.0,
                [0.0; 3],
            )],
        )
        .expect("a checked Point to Level Set Relation lowers to a finite value functional");
        let inequality = CanonicalAffineInequality::new(
            Some(FunctionalUse::new(functional, provenance.clone())),
            vec![SemanticLatentCoefficient {
                latent,
                coefficient: -1.0,
            }],
            provenance.clone(),
            FunctionalDimension::FieldValue,
            sense,
            orientation.bound(relation.minimum_offset()),
            violation_channel(relation.configuration(), &provenance),
        );
        lowering.push_bound(
            inequality,
            SourceBoundKind::PointToLevelSetRelation {
                relation: relation.clone(),
            },
        );
    }

    let has_absolute_observation = snapshot
        .inner
        .observations
        .iter()
        .any(|observation| matches!(observation, ObservationInput::FieldValue(_)))
        || snapshot.inner.covariance_groups.iter().any(|group| {
            group
                .members()
                .iter()
                .any(|member| matches!(member, CovarianceGroupMember::FieldValue(_)))
        })
        || !snapshot.inner.field_value_bounds.is_empty();
    let primary_gauge = (!has_absolute_observation)
        .then(|| snapshot.inner.additive_field_gauges.first())
        .flatten()
        .map(|gauge| gauge.source_id());
    for gauge in &snapshot.inner.additive_field_gauges {
        let participation = if primary_gauge == Some(gauge.source_id()) {
            CanonicalEqualityParticipation::SolverConstraint
        } else {
            CanonicalEqualityParticipation::VerificationOnly
        };
        match gauge.reference() {
            AdditiveFieldGaugeReference::Point(point) => {
                let role = SemanticRolePath::new("additive-field-gauge/point");
                let edge = value_constraints.add_absolute_support(
                    point.components(),
                    gauge.value(),
                    gauge.source_id(),
                );
                lowering.record_graph_conflict(&edge, &role);
                let participation = if participation
                    == CanonicalEqualityParticipation::SolverConstraint
                    && matches!(edge, CanonicalValueEdgeOutcome::Independent)
                {
                    CanonicalEqualityParticipation::SolverConstraint
                } else {
                    CanonicalEqualityParticipation::VerificationOnly
                };
                lowering.push_source(field_equality(
                    &ScalarObservation {
                        source_id: gauge.source_id().clone(),
                        group_id: None,
                        support: point.components(),
                        functional: ScalarFunctionalDescriptor::field_value(),
                        semantic_role: role,
                        target: gauge.value(),
                    },
                    participation,
                ));
            }
            AdditiveFieldGaugeReference::LevelSet(group_id) => {
                let role = SemanticRolePath::new("additive-field-gauge/level-set");
                let edge = value_constraints.add_absolute_group(
                    group_id,
                    gauge.value(),
                    gauge.source_id(),
                );
                lowering.record_graph_conflict(&edge, &role);
                let participation = if participation
                    == CanonicalEqualityParticipation::SolverConstraint
                    && matches!(edge, CanonicalValueEdgeOutcome::Independent)
                {
                    CanonicalEqualityParticipation::SolverConstraint
                } else {
                    CanonicalEqualityParticipation::VerificationOnly
                };
                let provenance =
                    relation_provenance(gauge.source_id().clone(), Some(group_id.clone()), role);
                lowering.push_source(CanonicalHardEquality::new(
                    None,
                    vec![SemanticLatentCoefficient {
                        latent: latent_index_by_group[group_id],
                        coefficient: 1.0,
                    }],
                    provenance,
                    FunctionalDimension::FieldValue,
                    gauge.value(),
                    participation,
                ));
            }
        }
    }
    record_direct_bound_conflicts(&mut lowering, &mut value_constraints, &hard_bound_facts);
    lowering.record_direct_derivative_bound_conflicts(&hard_derivative_bound_facts);
    lowering
}

#[derive(Debug, Clone)]
struct HardBoundFact {
    source_id: SourceId,
    support: [f64; 3],
    sense: CanonicalInequalitySense,
    bound: f64,
    semantic_role: SemanticRolePath,
}

#[derive(Debug, Clone)]
struct HardDerivativeBoundFact {
    source_id: SourceId,
    functional_key: CanonicalEqualityKey,
    sense: CanonicalInequalitySense,
    bound: f64,
    semantic_role: SemanticRolePath,
}

fn canonical_inequality_functional_key(
    inequality: &CanonicalAffineInequality,
) -> (CanonicalEqualityKey, f64) {
    let sign = canonical_normalization_sign(inequality.field(), inequality.latent_coefficients());
    let equality = CanonicalHardEquality::new(
        inequality.field().cloned(),
        inequality.latent_coefficients().to_vec(),
        inequality.provenance().clone(),
        inequality.dimension(),
        inequality.bound(),
        CanonicalEqualityParticipation::VerificationOnly,
    );
    (normalized_equality_key(&equality).0, sign)
}

fn record_direct_bound_conflicts(
    lowering: &mut EqualityLowering,
    value_constraints: &mut CanonicalValueConstraintForest,
    facts: &[HardBoundFact],
) {
    let mut by_support = BTreeMap::<[u64; 3], Vec<&HardBoundFact>>::new();
    for fact in facts {
        by_support
            .entry(canonical_support_bits(fact.support))
            .or_default()
            .push(fact);
    }
    for (support, support_facts) in by_support {
        for (index, left) in support_facts.iter().enumerate() {
            for right in &support_facts[index + 1..] {
                let incompatible = match (left.sense, right.sense) {
                    (CanonicalInequalitySense::Lower, CanonicalInequalitySense::Upper) => {
                        left.bound > right.bound
                    }
                    (CanonicalInequalitySense::Upper, CanonicalInequalitySense::Lower) => {
                        right.bound > left.bound
                    }
                    _ => false,
                };
                if incompatible {
                    push_direct_bound_conflict(lowering, left, right);
                }
            }
        }
        if let Some(absolute) = value_constraints.absolute_target_for_support(support) {
            for fact in support_facts {
                let incompatible = match fact.sense {
                    CanonicalInequalitySense::Lower => absolute.target < fact.bound,
                    CanonicalInequalitySense::Upper => absolute.target > fact.bound,
                };
                if incompatible {
                    let (first_source, first_target, second_source, second_target) =
                        if absolute.source_id <= fact.source_id {
                            (
                                absolute.source_id.clone(),
                                absolute.target,
                                fact.source_id.clone(),
                                fact.bound,
                            )
                        } else {
                            (
                                fact.source_id.clone(),
                                fact.bound,
                                absolute.source_id.clone(),
                                absolute.target,
                            )
                        };
                    lowering
                        .direct_input_conflicts
                        .push(DirectInputConflictEvidence::new(
                            first_source,
                            second_source,
                            fact.semantic_role.clone(),
                            first_target,
                            second_target,
                        ));
                }
            }
        }
    }
    lowering.direct_input_conflicts.sort_by(|left, right| {
        left.semantic_role()
            .cmp(right.semantic_role())
            .then_with(|| left.first_source().cmp(right.first_source()))
            .then_with(|| left.second_source().cmp(right.second_source()))
    });
}

fn push_direct_bound_conflict(
    lowering: &mut EqualityLowering,
    left: &HardBoundFact,
    right: &HardBoundFact,
) {
    let (first, second) = if left.source_id <= right.source_id {
        (left, right)
    } else {
        (right, left)
    };
    lowering
        .direct_input_conflicts
        .push(DirectInputConflictEvidence::new(
            first.source_id.clone(),
            second.source_id.clone(),
            second.semantic_role.clone(),
            first.bound,
            second.bound,
        ));
}

fn field_equality(
    observation: &ScalarObservation,
    participation: CanonicalEqualityParticipation,
) -> CanonicalHardEquality {
    field_equality_with_representer_span(
        observation,
        participation,
        FunctionalRepresenterSpan::Native,
    )
}

fn field_equality_with_representer_span(
    observation: &ScalarObservation,
    participation: CanonicalEqualityParticipation,
    representer_span: FunctionalRepresenterSpan,
) -> CanonicalHardEquality {
    let provenance = relation_provenance(
        observation.source_id.clone(),
        observation.group_id.clone(),
        observation.semantic_role.clone(),
    );
    let functional = CanonicalFunctional::new(
        observation.functional.dimension,
        vec![FunctionalTerm::new(
            observation.support,
            observation.functional.value_coefficient,
            observation.functional.gradient_coefficient,
        )],
    )
    .expect("checked public observations lower to a finite nonzero functional");
    CanonicalHardEquality::new(
        Some(FunctionalUse::with_representer_span(
            functional,
            provenance.clone(),
            representer_span,
        )),
        Vec::new(),
        provenance,
        observation.functional.dimension,
        observation.target,
        participation,
    )
}

fn soft_field_equality(observation: &ScalarObservation) -> CanonicalSoftEquality {
    soft_field_equality_with_representer_span(observation, FunctionalRepresenterSpan::Native)
}

fn soft_field_equality_with_representer_span(
    observation: &ScalarObservation,
    representer_span: FunctionalRepresenterSpan,
) -> CanonicalSoftEquality {
    let provenance = relation_provenance(
        observation.source_id.clone(),
        observation.group_id.clone(),
        observation.semantic_role.clone(),
    );
    let functional = CanonicalFunctional::new(
        observation.functional.dimension,
        vec![FunctionalTerm::new(
            observation.support,
            observation.functional.value_coefficient,
            observation.functional.gradient_coefficient,
        )],
    )
    .expect("checked public observations lower to a finite nonzero functional");
    CanonicalSoftEquality::new(
        FunctionalUse::with_representer_span(functional, provenance, representer_span),
        observation.target,
    )
}

fn relation_provenance(
    source_id: SourceId,
    group_id: Option<GroupId>,
    semantic_role: SemanticRolePath,
) -> UsageProvenance {
    let relation_id = RelationId::new(format!("{}:{}", source_id.as_str(), semantic_role.as_str()));
    let residual_id = ResidualId::new(format!("{}/residual", relation_id.as_str()));
    UsageProvenance::new(source_id, group_id, relation_id, residual_id, semantic_role)
}

fn relation_provenance_for_groups(
    source_id: SourceId,
    group_ids: Vec<GroupId>,
    semantic_role: SemanticRolePath,
) -> UsageProvenance {
    let relation_id = RelationId::new(format!("{}:{}", source_id.as_str(), semantic_role.as_str()));
    let residual_id = ResidualId::new(format!("{}/residual", relation_id.as_str()));
    UsageProvenance::new_with_groups(
        source_id,
        group_ids,
        relation_id,
        residual_id,
        semantic_role,
    )
}

fn signed_coefficient_bits(sign: f64, coefficient: f64) -> u64 {
    let value = sign * coefficient;
    if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    }
}

fn directed_normal_assessments(
    solution: &CubicExecutionSolution,
    source_relations: &[SourceHardRelation],
    source_bound_relations: &[SourceBoundRelation],
    resolved_normals: &[ResolvedNormalInput],
) -> Vec<DirectedNormalAssessment> {
    let qp = solution
        .qp
        .as_ref()
        .expect("normal slope relations select the QP execution route");
    resolved_normals
        .iter()
        .enumerate()
        .map(|(normal_index, normal)| {
            let sample = solution.field.sample(normal.location.components());
            let direction = normal.direction.components();
            let recovered_slope = direction
                .into_iter()
                .zip(sample.gradient)
                .map(|(direction, gradient)| direction * gradient)
                .sum::<f64>();
            let projection = std::array::from_fn(|axis| {
                crate::math::canonical_zero(
                    sample.gradient[axis] - recovered_slope * direction[axis],
                )
            });
            let projection_residual = Vector3::try_new(projection[0], projection[1], projection[2])
                .expect("accepted recovered projection components are finite");
            let projection_residual_norm = stable_vector_norm(projection);
            let recovered_bound = source_bound_relations
                .iter()
                .filter_map(|source_relation| match source_relation.kind {
                    SourceBoundKind::NormalSlope {
                        normal_index: candidate,
                    } if candidate == normal_index => {
                        Some(&solution.affine_inequalities[source_relation.canonical_index])
                    }
                    _ => None,
                })
                .next()
                .expect("every resolved normal lowers one slope relation");
            let direction_loss = solution.soft_objectives.iter().find_map(|objective| {
                if objective.block_kind != CanonicalSoftResidualBlockKind::NormalProjection {
                    return None;
                }
                let relation = objective
                    .canonical_indices
                    .first()
                    .map(|index| &solution.soft_equalities[*index])?;
                (relation.provenance.source() == &normal.source_id)
                    .then_some(objective.objective_contribution)
            });
            let direction_tolerance = matches!(
                normal.direction_enforcement.configuration(),
                NormalDirectionConfiguration::Hard
            )
            .then(|| {
                source_relations
                    .iter()
                    .filter(|relation| {
                        relation.equality.provenance().source() == &normal.source_id
                            && relation.kind
                                == SourceHardRelationKind::NormalProjection { normal_index }
                    })
                    .map(|relation| {
                        qp.hard_relation_tolerances[relation.canonical_index].physical_tolerance
                    })
                    .fold(0.0_f64, f64::max)
            });
            DirectedNormalAssessment {
                source_id: normal.source_id.clone(),
                direction_semantic_role: SemanticRolePath::new(
                    "directed-normal/direction-projection",
                ),
                slope_semantic_role: source_bound_relations
                    .iter()
                    .find_map(|source_relation| match source_relation.kind {
                        SourceBoundKind::NormalSlope {
                            normal_index: candidate,
                        } if candidate == normal_index => Some(
                            source_relation
                                .inequality
                                .provenance()
                                .semantic_role()
                                .clone(),
                        ),
                        _ => None,
                    })
                    .expect("every resolved normal retains slope role provenance"),
                direction: normal.direction,
                recovered_gradient: Vector3::try_new(
                    sample.gradient[0],
                    sample.gradient[1],
                    sample.gradient[2],
                )
                .expect("an accepted recovered gradient is finite"),
                projection_residual,
                projection_residual_norm,
                recovered_slope,
                minimum_slope: normal.minimum_slope,
                slope_slack: recovered_bound.slack,
                slope_violation: recovered_bound.violation,
                direction_tolerance,
                slope_tolerance: recovered_bound.tolerance,
                slope_active_state: public_bound_active_state(recovered_bound),
                direction_enforcement: normal.direction_enforcement,
                minimum_slope_enforcement: normal.minimum_slope_enforcement,
                direction_loss,
                slope_loss: recovered_bound.objective_contribution,
                input_axis: normal.input_axis,
                polarity_resolution_source_id: normal.polarity_resolution_source_id.clone(),
                polarity_selection: normal.polarity_selection,
            }
        })
        .collect()
}

fn stable_vector_norm(vector: [f64; 3]) -> f64 {
    let scale = vector
        .iter()
        .map(|component| component.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    scale
        * vector
            .into_iter()
            .map(|component| (component / scale).powi(2))
            .sum::<f64>()
            .sqrt()
}

fn success_report_qp(
    snapshot: &ProblemSnapshot,
    problem_size: ProblemSize,
    solution: &CubicExecutionSolution,
    source_relations: &[SourceHardRelation],
    source_bound_relations: &[SourceBoundRelation],
    resolved_normals: &[ResolvedNormalInput],
    shared_level_values: &[SharedLevelValue],
) -> FitReport {
    let qp = solution
        .qp
        .as_ref()
        .expect("an affine-inequality execution retains QP evidence");
    let mut hard_relations = source_relations
        .iter()
        .map(|source_relation| {
            hard_relation_assessment(
                source_relation,
                solution.hard_equalities[source_relation.canonical_index].value,
                qp.hard_relation_tolerances[source_relation.canonical_index],
            )
        })
        .collect::<Vec<_>>();
    hard_relations.sort_by(|left, right| {
        left.source_id
            .cmp(&right.source_id)
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut field_value_bounds = source_bound_relations
        .iter()
        .filter_map(|source_relation| {
            matches!(&source_relation.kind, SourceBoundKind::FieldValue).then(|| {
                field_value_bound_assessment(
                    source_relation,
                    &solution.affine_inequalities[source_relation.canonical_index],
                )
            })
        })
        .collect::<Vec<_>>();
    field_value_bounds.sort_by(|left, right| {
        left.source_id
            .cmp(&right.source_id)
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut directional_derivative_intervals = source_bound_relations
        .iter()
        .filter_map(|source_relation| {
            let SourceBoundKind::DirectionalDerivative { direction } = &source_relation.kind else {
                return None;
            };
            Some(directional_derivative_interval_assessment(
                source_relation,
                &solution.affine_inequalities[source_relation.canonical_index],
                *direction,
            ))
        })
        .collect::<Vec<_>>();
    directional_derivative_intervals.sort_by(|left, right| {
        left.source_id
            .cmp(&right.source_id)
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut field_separation_intervals = source_bound_relations
        .iter()
        .filter_map(|source_relation| {
            field_separation_interval_assessment(
                source_relation,
                &solution.affine_inequalities[source_relation.canonical_index],
                shared_level_values,
            )
        })
        .collect::<Vec<_>>();
    field_separation_intervals.sort_by(|left, right| {
        left.source_id()
            .cmp(right.source_id())
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut point_to_level_set_relations = source_bound_relations
        .iter()
        .filter_map(|source_relation| {
            point_to_level_set_relation_assessment(
                source_relation,
                &solution.affine_inequalities[source_relation.canonical_index],
                shared_level_values,
            )
        })
        .collect::<Vec<_>>();
    point_to_level_set_relations.sort_by(|left, right| {
        left.source_id()
            .cmp(right.source_id())
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut shared_level_set_relations = source_bound_relations
        .iter()
        .filter_map(|source_relation| {
            shared_level_set_relation_assessment(
                source_relation,
                &solution.affine_inequalities[source_relation.canonical_index],
                shared_level_values,
            )
        })
        .collect::<Vec<_>>();
    shared_level_set_relations.sort_by(|left, right| {
        left.source_id()
            .cmp(right.source_id())
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut soft_field_values = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| soft_field_value_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    soft_field_values.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let mut soft_gradients = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| soft_gradient_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    soft_gradients.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let mut soft_tangents = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| soft_tangent_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    soft_tangents.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let mut covariance_groups = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| covariance_group_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    covariance_groups.sort_by(|left, right| left.group_id.cmp(&right.group_id));
    let directed_normals = directed_normal_assessments(
        solution,
        source_relations,
        source_bound_relations,
        resolved_normals,
    );
    let accepted_backend = &qp.attempts[qp.accepted_attempt].backend;
    FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: snapshot.inner.fit_configuration.numerical_policy(),
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations,
        field_value_bounds,
        directional_derivative_intervals,
        field_separation_intervals,
        point_to_level_set_relations,
        shared_level_set_relations,
        soft_field_values,
        soft_gradients,
        directed_normals,
        soft_tangents,
        covariance_groups,
        shared_level_values: shared_level_values.to_vec(),
        field_energy: Some(solution.field_energy),
        total_objective: Some(solution.total_objective),
        backend_fingerprint: Some(public_qp_backend_fingerprint(accepted_backend)),
        attempts: public_qp_attempts(&qp.attempts),
        recovery_verification: None,
        direct_input_conflicts: Vec::new(),
        relation_graph_conflicts: Vec::new(),
        shared_level_set_relation_conflicts: Vec::new(),
        execution_failure: None,
        cubic_analysis: Some(public_cubic_analysis(&solution.representation)),
        backend_rank: None,
        interpretable_rank_deficiency: None,
        inertia: None,
        canonical_acceptance: Some(public_qp_success_acceptance(solution)),
        capacity: None,
        analysis_failure: None,
        infeasibility_certificate: None,
        recession_ray: None,
        unidentified_additive_gauge: None,
        uninformative_shared_level_sets: Vec::new(),
        unresolved_axial_normals: Vec::new(),
    }
}

fn field_value_bound_assessment(
    source_relation: &SourceBoundRelation,
    recovered: &RecoveredAffineInequality,
) -> FieldValueBoundAssessment {
    let (quadratic_penalty, linear_violation_penalty) =
        public_bound_penalties(&source_relation.inequality);
    FieldValueBoundAssessment {
        source_id: source_relation.inequality.provenance().source().clone(),
        semantic_role: source_relation
            .inequality
            .provenance()
            .semantic_role()
            .clone(),
        side: public_bound_side(source_relation.inequality.sense()),
        bound: source_relation.inequality.bound(),
        recovered_value: recovered.value,
        slack: recovered.slack,
        violation: recovered.violation,
        tolerance: recovered.tolerance,
        active_state: public_bound_active_state(recovered),
        quadratic_penalty,
        linear_violation_penalty,
        loss: recovered.objective_contribution,
    }
}

fn directional_derivative_interval_assessment(
    source_relation: &SourceBoundRelation,
    recovered: &RecoveredAffineInequality,
    direction: Vector3,
) -> DirectionalDerivativeIntervalAssessment {
    let (quadratic_penalty, linear_violation_penalty) =
        public_bound_penalties(&source_relation.inequality);
    DirectionalDerivativeIntervalAssessment {
        source_id: source_relation.inequality.provenance().source().clone(),
        semantic_role: source_relation
            .inequality
            .provenance()
            .semantic_role()
            .clone(),
        direction,
        side: public_bound_side(source_relation.inequality.sense()),
        bound: source_relation.inequality.bound(),
        recovered_directional_derivative: recovered.value,
        slack: recovered.slack,
        violation: recovered.violation,
        tolerance: recovered.tolerance,
        active_state: public_bound_active_state(recovered),
        quadratic_penalty,
        linear_violation_penalty,
        loss: recovered.objective_contribution,
    }
}

fn shared_level_set_relation_assessment(
    source_relation: &SourceBoundRelation,
    recovered: &RecoveredAffineInequality,
    shared_level_values: &[SharedLevelValue],
) -> Option<SharedLevelSetRelationAssessment> {
    let SourceBoundKind::SharedLevelSetRelation {
        relation,
        orientation,
    } = &source_relation.kind
    else {
        return None;
    };
    let value = |group_id: &GroupId| {
        shared_level_values
            .iter()
            .find(|level| level.group_id() == group_id)
            .map(SharedLevelValue::value)
            .expect("every checked shared-level-set relation recovers both semantic latents")
    };
    let recovered_values = if let (Some(younger), Some(older)) =
        (relation.younger_group_id(), relation.older_group_id())
    {
        RecoveredSharedLevelSetRelationValues::StratigraphicAge {
            younger: value(younger),
            older: value(older),
        }
    } else {
        RecoveredSharedLevelSetRelationValues::FieldLevelOrder {
            lower: value(
                relation
                    .lower_group_id()
                    .expect("Field Level Order owns a lower group"),
            ),
            upper: value(
                relation
                    .upper_group_id()
                    .expect("Field Level Order owns an upper group"),
            ),
        }
    };
    let (quadratic_penalty, linear_violation_penalty) =
        public_bound_penalties(&source_relation.inequality);
    Some(SharedLevelSetRelationAssessment {
        relation: relation.clone(),
        orientation: orientation.clone(),
        semantic_role: source_relation
            .inequality
            .provenance()
            .semantic_role()
            .clone(),
        recovered_values,
        recovered_field_separation: recovered.value,
        slack: recovered.slack,
        violation: recovered.violation,
        tolerance: recovered.tolerance,
        active_state: public_bound_active_state(recovered),
        quadratic_penalty,
        linear_violation_penalty,
        loss: recovered.objective_contribution,
    })
}

fn field_separation_interval_assessment(
    source_relation: &SourceBoundRelation,
    recovered: &RecoveredAffineInequality,
    shared_level_values: &[SharedLevelValue],
) -> Option<FieldSeparationIntervalAssessment> {
    let SourceBoundKind::FieldSeparationInterval { relation } = &source_relation.kind else {
        return None;
    };
    let value = |group_id: &GroupId| {
        shared_level_values
            .iter()
            .find(|level| level.group_id() == group_id)
            .map(SharedLevelValue::value)
            .expect("every checked field-separation relation recovers both semantic latents")
    };
    let (quadratic_penalty, linear_violation_penalty) =
        public_bound_penalties(&source_relation.inequality);
    Some(FieldSeparationIntervalAssessment {
        relation: relation.clone(),
        semantic_role: source_relation
            .inequality
            .provenance()
            .semantic_role()
            .clone(),
        side: public_bound_side(source_relation.inequality.sense()),
        bound: source_relation.inequality.bound(),
        recovered_reference_value: value(relation.reference_group_id()),
        recovered_target_value: value(relation.target_group_id()),
        recovered_field_separation: recovered.value,
        slack: recovered.slack,
        violation: recovered.violation,
        tolerance: recovered.tolerance,
        active_state: public_bound_active_state(recovered),
        quadratic_penalty,
        linear_violation_penalty,
        loss: recovered.objective_contribution,
    })
}

fn point_to_level_set_relation_assessment(
    source_relation: &SourceBoundRelation,
    recovered: &RecoveredAffineInequality,
    shared_level_values: &[SharedLevelValue],
) -> Option<PointToLevelSetRelationAssessment> {
    let SourceBoundKind::PointToLevelSetRelation { relation } = &source_relation.kind else {
        return None;
    };
    let recovered_level_value = shared_level_values
        .iter()
        .find(|level| level.group_id() == relation.group_id())
        .map(SharedLevelValue::value)
        .expect("every checked point-to-level relation recovers its semantic latent");
    let recovered_point_value = recovered.value + recovered_level_value;
    let recovered_field_offset = relation
        .orientation()
        .recovered_field_offset(recovered.value);
    let (quadratic_penalty, linear_violation_penalty) =
        public_bound_penalties(&source_relation.inequality);
    Some(PointToLevelSetRelationAssessment {
        relation: relation.clone(),
        semantic_role: source_relation
            .inequality
            .provenance()
            .semantic_role()
            .clone(),
        recovered_point_value,
        recovered_level_value,
        recovered_field_offset,
        slack: recovered.slack,
        violation: recovered.violation,
        tolerance: recovered.tolerance,
        active_state: public_bound_active_state(recovered),
        quadratic_penalty,
        linear_violation_penalty,
        loss: recovered.objective_contribution,
    })
}

fn public_bound_penalties(
    inequality: &CanonicalAffineInequality,
) -> (Option<QuadraticPenalty>, Option<LinearViolationPenalty>) {
    inequality
        .violation_channel()
        .map(|channel| match channel.loss() {
            CanonicalViolationLoss::QuadraticPenalty { weight } => (
                Some(QuadraticPenalty::try_new(weight).expect("canonical penalty stays checked")),
                None,
            ),
            CanonicalViolationLoss::LinearViolationPenalty { weight } => (
                None,
                Some(
                    LinearViolationPenalty::try_new(weight)
                        .expect("canonical violation penalty stays checked"),
                ),
            ),
        })
        .unwrap_or((None, None))
}

fn public_bound_side(sense: CanonicalInequalitySense) -> BoundSide {
    match sense {
        CanonicalInequalitySense::Lower => BoundSide::Lower,
        CanonicalInequalitySense::Upper => BoundSide::Upper,
    }
}

fn public_bound_active_state(recovered: &RecoveredAffineInequality) -> BoundActiveState {
    if recovered.slack <= recovered.tolerance {
        BoundActiveState::Active
    } else {
        BoundActiveState::Inactive
    }
}

fn public_qp_success_acceptance(solution: &CubicExecutionSolution) -> CanonicalAcceptanceEvidence {
    let qp = solution.qp.as_ref().expect("QP success retains evidence");
    CanonicalAcceptanceEvidence::new(CanonicalAcceptanceEvidenceParts {
        accepted: solution.canonical_acceptance_verified,
        backend_standard_form_verified: solution.backend_standard_form_verified,
        recovery_finite: true,
        provenance_verified: qp.provenance_verified,
        side_condition: Some(public_side_condition(solution.side_condition)),
        hard_residual_maxima: Some((
            solution.hard_equality_violations.field_value,
            solution.hard_equality_violations.field_value_per_length,
        )),
        polynomial_round_trip_error: Some(qp.polynomial_round_trip_error),
        field_coefficient_round_trip_error: Some(qp.field_coefficient_round_trip_error),
        field_energy_round_trip_error: Some(qp.field_energy_round_trip_error),
        whitening_round_trip_error: Some(qp.whitening_round_trip_error),
        objective_round_trip_error: Some(qp.objective_round_trip_error),
        objective_verified: qp.objective_round_trip_error
            <= crate::numerical::EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit,
        tolerance_round_trip_error: Some(
            qp.hard_relation_tolerances
                .iter()
                .chain(&qp.affine_relation_tolerances)
                .map(|tolerance| tolerance.round_trip_error)
                .fold(0.0_f64, f64::max),
        ),
        hard_affine_inequality_violation_max: Some(
            solution
                .affine_inequalities
                .iter()
                .filter(|relation| relation.violation_loss.is_none())
                .map(|relation| relation.violation)
                .fold(0.0_f64, f64::max),
        ),
        backend_standard_form_residual: Some(qp.physical_standard_form_violation),
        physical_convex_residual: Some(PublicConvexResidualEvidence::new(
            ConvexResidualEvidenceParts {
                primal: qp.physical_residuals.primal,
                dual: qp.physical_residuals.dual,
                stationarity: qp.physical_residuals.stationarity,
                complementarity: qp.physical_residuals.complementarity,
                relative_gap: qp.physical_residuals.relative_gap,
            },
        )),
        scaling_round_trip_error: Some(qp.scaling_round_trip_error),
        reduction_round_trip_error: Some(qp.reduction_round_trip_error),
        backend_internal_scaling_round_trip_error: Some(
            qp.backend_internal_scaling_round_trip_error,
        ),
    })
}

fn success_report(
    snapshot: &ProblemSnapshot,
    planned_problem_size: ProblemSize,
    solution: &CubicEqualitySolution,
    source_relations: &[SourceHardRelation],
    shared_level_values: &[SharedLevelValue],
) -> FitReport {
    let problem_size = ProblemSize {
        input_observations: planned_problem_size.input_observations,
        scalar_hard_relations: planned_problem_size.scalar_hard_relations,
        scalar_soft_relations: planned_problem_size.scalar_soft_relations,
        canonical_hard_equalities: Some(solution.assembly.canonical_hard_equalities),
        canonical_soft_equalities: Some(solution.soft_equalities.len()),
        center_coefficients: Some(solution.assembly.field_coefficients),
        semantic_latents: solution.assembly.semantic_latents,
        auxiliary_variables: 0,
        quadratic_objective_terms: solution.soft_objectives.len(),
        linear_objective_terms: planned_problem_size.linear_objective_terms,
        affine_inequality_constraints: planned_problem_size.affine_inequality_constraints,
        cone_blocks: 0,
        primal_variables: Some(solution.assembly.primal_variables),
        equality_constraints: Some(
            solution.assembly.side_conditions + solution.assembly.hard_equalities,
        ),
        kkt_dimension: Some(solution.backend.capacity.kkt_dimension),
    };
    let mut hard_relations = source_relations
        .iter()
        .map(|source_relation| {
            let tolerance = solution.relation_tolerances[source_relation.canonical_index];
            let recovered_value = solution.hard_equalities[source_relation.canonical_index].value;
            hard_relation_assessment(source_relation, recovered_value, tolerance)
        })
        .collect::<Vec<_>>();
    hard_relations.sort_by(|left, right| {
        left.source_id
            .cmp(&right.source_id)
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut soft_field_values = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| soft_field_value_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    soft_field_values.sort_by(|left, right| {
        left.source_id
            .cmp(&right.source_id)
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    let mut soft_gradients = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| soft_gradient_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    soft_gradients.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let mut soft_tangents = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| soft_tangent_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    soft_tangents.sort_by(|left, right| left.source_id.cmp(&right.source_id));
    let mut covariance_groups = solution
        .soft_objectives
        .iter()
        .filter_map(|objective| covariance_group_assessment(objective, &solution.soft_equalities))
        .collect::<Vec<_>>();
    covariance_groups.sort_by(|left, right| left.group_id.cmp(&right.group_id));
    let backend_fingerprint = public_backend_fingerprint(&solution.backend.backend);
    let attempts = public_attempts(&solution.backend.attempts);
    FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: solution.backend.numerical_policy,
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations,
        field_value_bounds: Vec::new(),
        directional_derivative_intervals: Vec::new(),
        field_separation_intervals: Vec::new(),
        point_to_level_set_relations: Vec::new(),
        shared_level_set_relations: Vec::new(),
        soft_field_values,
        soft_gradients,
        directed_normals: Vec::new(),
        soft_tangents,
        covariance_groups,
        shared_level_values: shared_level_values.to_vec(),
        field_energy: Some(solution.field_energy),
        total_objective: Some(solution.total_objective),
        backend_fingerprint: Some(backend_fingerprint),
        attempts,
        recovery_verification: None,
        direct_input_conflicts: Vec::new(),
        relation_graph_conflicts: Vec::new(),
        shared_level_set_relation_conflicts: Vec::new(),
        execution_failure: None,
        cubic_analysis: Some(public_cubic_analysis(&solution.representation)),
        backend_rank: Some(public_rank_evidence(
            RankEvidenceDomain::BackendKkt,
            solution.backend.capacity.kkt_dimension,
            &solution.backend.rank,
        )),
        interpretable_rank_deficiency: None,
        inertia: Some(public_inertia(
            solution.backend.expected_inertia,
            solution.backend.observed_inertia,
            true,
        )),
        canonical_acceptance: Some(public_success_acceptance(solution)),
        capacity: None,
        analysis_failure: None,
        infeasibility_certificate: None,
        recession_ray: None,
        unidentified_additive_gauge: None,
        uninformative_shared_level_sets: Vec::new(),
        unresolved_axial_normals: Vec::new(),
    }
}

fn hard_relation_assessment(
    source_relation: &SourceHardRelation,
    recovered_value: f64,
    tolerance: CanonicalRelationToleranceEvidence,
) -> HardRelationAssessment {
    let equality = &source_relation.equality;
    let target = equality.target();
    HardRelationAssessment {
        source_id: equality.provenance().source().clone(),
        group_id: equality.provenance().group().cloned(),
        semantic_role: equality.provenance().semantic_role().clone(),
        dimension: match equality.dimension() {
            FunctionalDimension::FieldValue => ResidualDimension::FieldValue,
            FunctionalDimension::FieldValuePerLength => ResidualDimension::FieldValuePerLength,
        },
        target,
        recovered_value,
        residual: recovered_value - target,
        tolerance: tolerance.physical_tolerance,
        characteristic_scale: tolerance.characteristic_scale,
        relation_reference_scale: tolerance.relation_reference_scale,
        standard_tolerance: tolerance.standard_tolerance,
        scaled_kkt_tolerance: tolerance.scaled_kkt_tolerance,
        recovered_physical_tolerance: tolerance.recovered_physical_tolerance,
        tolerance_round_trip_error: tolerance.round_trip_error,
    }
}

fn soft_field_value_assessment(
    objective: &crate::cubic_equality::RecoveredSoftObjective,
    relations: &[crate::cubic_equality::RecoveredSoftEquality],
) -> Option<SoftFieldValueAssessment> {
    if objective.block_kind
        != CanonicalSoftResidualBlockKind::Independent(CanonicalSoftResidualMemberKind::FieldValue)
        || objective.canonical_indices.len() != 1
    {
        return None;
    }
    let relation = &relations[objective.canonical_indices[0]];
    let (quadratic_penalty, standard_deviation) = match &objective.loss {
        CanonicalSoftLoss::QuadraticPenalty { weight } => (
            Some(
                QuadraticPenalty::try_new(*weight)
                    .expect("canonical quadratic penalties retain checked public values"),
            ),
            None,
        ),
        CanonicalSoftLoss::StandardDeviation { standard_deviation } => (
            None,
            Some(
                StandardDeviation::try_new(*standard_deviation)
                    .expect("canonical standard deviations retain checked public values"),
            ),
        ),
        CanonicalSoftLoss::Covariance { .. } => return None,
    };
    Some(SoftFieldValueAssessment {
        source_id: relation.provenance.source().clone(),
        semantic_role: relation.provenance.semantic_role().clone(),
        target: relation.target,
        recovered_value: relation.value,
        residual: relation.residual,
        quadratic_penalty,
        standard_deviation,
        loss: objective.objective_contribution,
    })
}

fn soft_gradient_assessment(
    objective: &crate::cubic_equality::RecoveredSoftObjective,
    relations: &[crate::cubic_equality::RecoveredSoftEquality],
) -> Option<SoftGradientAssessment> {
    if objective.block_kind
        != CanonicalSoftResidualBlockKind::Independent(CanonicalSoftResidualMemberKind::Gradient)
        || objective.canonical_indices.len() != 3
    {
        return None;
    }
    let components = objective
        .canonical_indices
        .iter()
        .map(|index| &relations[*index])
        .collect::<Vec<_>>();
    if components
        .iter()
        .any(|relation| relation.provenance.source() != components[0].provenance.source())
    {
        return None;
    }
    let vector = |values: [f64; 3]| {
        Vector3::try_new(values[0], values[1], values[2])
            .expect("recovered finite gradient components form a public vector")
    };
    let target = vector(std::array::from_fn(|axis| components[axis].target));
    let recovered_gradient = vector(std::array::from_fn(|axis| components[axis].value));
    let residual = vector(std::array::from_fn(|axis| components[axis].residual));
    let (quadratic_penalty, standard_deviation, covariance) = match &objective.loss {
        CanonicalSoftLoss::QuadraticPenalty { weight } => (
            Some(QuadraticPenalty::try_new(*weight).expect("canonical penalty stays checked")),
            None,
            None,
        ),
        CanonicalSoftLoss::StandardDeviation { standard_deviation } => (
            None,
            Some(
                StandardDeviation::try_new(*standard_deviation)
                    .expect("canonical uncertainty stays checked"),
            ),
            None,
        ),
        CanonicalSoftLoss::Covariance {
            dimension,
            covariance,
            ..
        } => {
            let rows = covariance
                .chunks_exact(*dimension)
                .map(<[f64]>::to_vec)
                .collect();
            (
                None,
                None,
                Some(
                    CovarianceMatrix::try_from_rows(rows)
                        .expect("canonical covariance stays checked"),
                ),
            )
        }
    };
    Some(SoftGradientAssessment {
        source_id: components[0].provenance.source().clone(),
        semantic_role: SemanticRolePath::new("gradient-observation/vector"),
        target,
        recovered_gradient,
        residual,
        quadratic_penalty,
        standard_deviation,
        covariance,
        whitened_residual: objective.whitened_residual.clone().into(),
        whitening_round_trip_error: objective.whitening_round_trip_error,
        loss: objective.objective_contribution,
    })
}

fn soft_tangent_assessment(
    objective: &crate::cubic_equality::RecoveredSoftObjective,
    relations: &[crate::cubic_equality::RecoveredSoftEquality],
) -> Option<SoftTangentAssessment> {
    if objective.block_kind
        != CanonicalSoftResidualBlockKind::Independent(CanonicalSoftResidualMemberKind::Tangent)
        || objective.canonical_indices.len() != 1
    {
        return None;
    }
    let relation = &relations[objective.canonical_indices[0]];
    let (quadratic_penalty, standard_deviation) = match &objective.loss {
        CanonicalSoftLoss::QuadraticPenalty { weight } => (
            Some(QuadraticPenalty::try_new(*weight).expect("canonical penalty stays checked")),
            None,
        ),
        CanonicalSoftLoss::StandardDeviation { standard_deviation } => (
            None,
            Some(
                StandardDeviation::try_new(*standard_deviation)
                    .expect("canonical uncertainty stays checked"),
            ),
        ),
        CanonicalSoftLoss::Covariance { .. } => return None,
    };
    Some(SoftTangentAssessment {
        source_id: relation.provenance.source().clone(),
        semantic_role: relation.provenance.semantic_role().clone(),
        recovered_directional_derivative: relation.value,
        residual: relation.residual,
        quadratic_penalty,
        standard_deviation,
        loss: objective.objective_contribution,
    })
}

fn covariance_group_assessment(
    objective: &crate::cubic_equality::RecoveredSoftObjective,
    relations: &[crate::cubic_equality::RecoveredSoftEquality],
) -> Option<CovarianceGroupAssessment> {
    let group_id = objective.covariance_group.clone()?;
    let CanonicalSoftResidualBlockKind::CovarianceGroup {
        members: member_kinds,
    } = &objective.block_kind
    else {
        return None;
    };
    let (dimension, covariance_entries) = objective.loss.covariance_entries()?;
    let covariance = CovarianceMatrix::try_from_rows(
        covariance_entries
            .chunks_exact(dimension)
            .map(<[f64]>::to_vec)
            .collect(),
    )
    .expect("canonical group covariance stays checked");
    let components = objective
        .canonical_indices
        .iter()
        .map(|index| &relations[*index])
        .collect::<Vec<_>>();
    let mut members = Vec::with_capacity(member_kinds.len());
    let mut start = 0;
    for member_kind in member_kinds {
        let end = start + member_kind.component_count();
        let member_components = components.get(start..end)?;
        let source_id = member_components[0].provenance.source();
        if member_components
            .iter()
            .any(|relation| relation.provenance.source() != source_id)
        {
            return None;
        }
        let semantic_role = match member_kind {
            CanonicalSoftResidualMemberKind::Gradient => {
                SemanticRolePath::new("gradient-observation/vector")
            }
            CanonicalSoftResidualMemberKind::FieldValue
            | CanonicalSoftResidualMemberKind::Tangent => {
                member_components[0].provenance.semantic_role().clone()
            }
        };
        let residual_dimension = match member_components[0].dimension {
            FunctionalDimension::FieldValue => ResidualDimension::FieldValue,
            FunctionalDimension::FieldValuePerLength => ResidualDimension::FieldValuePerLength,
        };
        members.push(CovarianceGroupMemberAssessment {
            source_id: source_id.clone(),
            semantic_role,
            dimension: residual_dimension,
            target_components: member_components
                .iter()
                .map(|relation| relation.target)
                .collect(),
            recovered_components: member_components
                .iter()
                .map(|relation| relation.value)
                .collect(),
            residual_components: member_components
                .iter()
                .map(|relation| relation.residual)
                .collect(),
        });
        start = end;
    }
    if start != components.len() {
        return None;
    }
    Some(CovarianceGroupAssessment {
        group_id,
        covariance,
        members,
        whitened_residual: objective.whitened_residual.clone().into(),
        whitening_round_trip_error: objective.whitening_round_trip_error,
        objective_contribution: objective.objective_contribution,
    })
}

fn failure_report(
    mut report: FitReport,
    failure: &CubicEqualityFailure,
    source_relations: &[SourceHardRelation],
) -> FitReport {
    match failure {
        CubicEqualityFailure::Backend {
            failure,
            representation,
        } => {
            report.cubic_analysis = Some(public_cubic_analysis(representation));
            retain_kkt_failure(&mut report, failure);
        }
        CubicEqualityFailure::Representation(failure) => {
            retain_representation_failure(&mut report, failure, source_relations);
        }
        CubicEqualityFailure::RecoveryVerification {
            evidence,
            representation,
            backend,
        } => {
            report.cubic_analysis = Some(public_cubic_analysis(representation));
            report.backend_fingerprint = Some(public_backend_fingerprint(&backend.backend));
            report.attempts = public_attempts(&backend.attempts);
            report.backend_rank = Some(public_rank_evidence(
                RankEvidenceDomain::BackendKkt,
                backend.capacity.kkt_dimension,
                &backend.rank,
            ));
            report.inertia = Some(public_inertia(
                backend.expected_inertia,
                backend.observed_inertia,
                true,
            ));
            report.canonical_acceptance = public_failure_acceptance(evidence);
            report.hard_relations = public_failed_hard_relations(evidence, source_relations);
            let recovered_relations = evidence.soft_equalities.as_deref().unwrap_or_default();
            report.soft_field_values = evidence
                .soft_objectives
                .as_deref()
                .unwrap_or_default()
                .iter()
                .filter_map(|objective| soft_field_value_assessment(objective, recovered_relations))
                .collect();
            report.soft_field_values.sort_by(|left, right| {
                left.source_id
                    .cmp(&right.source_id)
                    .then_with(|| left.semantic_role.cmp(&right.semantic_role))
            });
            report.soft_gradients = evidence
                .soft_objectives
                .as_deref()
                .unwrap_or_default()
                .iter()
                .filter_map(|objective| soft_gradient_assessment(objective, recovered_relations))
                .collect();
            report
                .soft_gradients
                .sort_by(|left, right| left.source_id.cmp(&right.source_id));
            report.soft_tangents = evidence
                .soft_objectives
                .as_deref()
                .unwrap_or_default()
                .iter()
                .filter_map(|objective| soft_tangent_assessment(objective, recovered_relations))
                .collect();
            report
                .soft_tangents
                .sort_by(|left, right| left.source_id.cmp(&right.source_id));
            report.covariance_groups = evidence
                .soft_objectives
                .as_deref()
                .unwrap_or_default()
                .iter()
                .filter_map(|objective| covariance_group_assessment(objective, recovered_relations))
                .collect();
            report
                .covariance_groups
                .sort_by(|left, right| left.group_id.cmp(&right.group_id));
            report.recovery_verification = Some(public_recovery_evidence(evidence));
        }
        CubicEqualityFailure::EmptyEqualitySet
        | CubicEqualityFailure::AffineInequalityRequiresConvexQp
        | CubicEqualityFailure::NonFiniteTarget { .. } => {}
    }
    report
}

fn retain_representation_failure(
    report: &mut FitReport,
    failure: &RepresentationFailure,
    source_relations: &[SourceHardRelation],
) {
    match failure {
        RepresentationFailure::EmptyRepresenterSpan => {
            report.analysis_failure = Some(AnalysisFailureEvidence::EmptyRepresenterSpan);
        }
        RepresentationFailure::Capacity(evidence) => {
            report.capacity = Some(public_capacity(evidence));
        }
        RepresentationFailure::InvalidSolveCoordinateTransform {
            reason,
            solver_invoked,
        } => {
            report.analysis_failure =
                Some(AnalysisFailureEvidence::InvalidSolveCoordinateTransform {
                    reason: public_solve_coordinate_failure(*reason),
                    backend_invoked: *solver_invoked,
                });
        }
        RepresentationFailure::PolynomialRankDeficient { rank, mode } => {
            report.backend_rank = Some(RankEvidence::new(RankEvidenceParts {
                domain: RankEvidenceDomain::CubicPolynomialPairing,
                dimension: 4,
                rank: *rank,
                exact_zero_index: None,
                rrqr_ratio: None,
                singular_values: Vec::new(),
                svd_ratio: None,
                reject_ratio: None,
                accept_ratio: None,
                decision: RankDecision::RankDeficient,
                backend_invoked: mode.execution.solver_invoked,
            }));
            report.interpretable_rank_deficiency = Some(public_interpretable_rank_deficiency(
                RankDeficiencyConcept::CubicPi1FieldMode,
                RankEvidenceDomain::CubicPolynomialPairing,
                source_relations,
                mode.residual,
                mode.execution.solver_invoked,
                mode.execution.hidden_regularization_applied,
            ));
        }
        RepresentationFailure::PolynomialRankGrayZone { evidence } => {
            report.backend_rank = Some(RankEvidence::new(RankEvidenceParts {
                domain: RankEvidenceDomain::CubicPolynomialPairing,
                dimension: 4,
                rank: None,
                exact_zero_index: None,
                rrqr_ratio: Some(evidence.rrqr_ratio),
                singular_values: Vec::new(),
                svd_ratio: Some(evidence.svd_ratio),
                reject_ratio: Some(evidence.reject_ratio),
                accept_ratio: Some(evidence.accept_ratio),
                decision: RankDecision::NumericalDecisionGrayZone,
                backend_invoked: evidence.backend_invoked,
            }));
        }
        RepresentationFailure::QuotientPivotRequiresPrecisionRescue {
            quotient_dimension,
            pivot_index,
            interval,
            execution,
        } => {
            report.analysis_failure = Some(
                AnalysisFailureEvidence::QuotientPivotRequiresPrecisionRescue {
                    quotient_dimension: *quotient_dimension,
                    pivot_index: *pivot_index,
                    interval: interval.map(|interval| {
                        CubicLltPivotInterval::new(CubicLltPivotIntervalParts {
                            lower: interval.lower,
                            upper: interval.upper,
                        })
                    }),
                    backend_invoked: execution.solver_invoked,
                },
            );
        }
        RepresentationFailure::QuotientFactorizationNotPositive {
            quotient_dimension,
            pivot_index,
            interval,
            execution,
        } => {
            report.analysis_failure =
                Some(AnalysisFailureEvidence::QuotientFactorizationNotPositive {
                    quotient_dimension: *quotient_dimension,
                    pivot_index: *pivot_index,
                    interval: CubicLltPivotInterval::new(CubicLltPivotIntervalParts {
                        lower: interval.lower,
                        upper: interval.upper,
                    }),
                    backend_invoked: execution.solver_invoked,
                });
        }
        RepresentationFailure::AlgebraicAnalysisFailure {
            stage,
            solver_invoked,
        } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::NumericalAnalysis {
                stage: public_cubic_analysis_stage(*stage),
                backend_invoked: *solver_invoked,
            });
        }
        RepresentationFailure::AlgebraicAnalysisWorkspaceAllocation {
            stage,
            bytes,
            alignment,
            solver_invoked,
        } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::WorkspaceAllocation {
                stage: public_cubic_analysis_stage(*stage),
                bytes: *bytes,
                alignment: *alignment,
                backend_invoked: *solver_invoked,
            });
        }
        RepresentationFailure::NullSpaceWorkspaceAllocation => {
            report.analysis_failure = Some(AnalysisFailureEvidence::NullSpaceWorkspaceAllocation);
        }
        RepresentationFailure::NullSpaceDefect { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::NullSpaceDefect,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::HouseholderOrthogonalityContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::HouseholderOrthogonalityError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::CanonicalResponseRoundTripContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::CanonicalResponseRoundTripError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::QuotientLltContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::QuotientLltBackwardError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::QuotientFieldEnergyIdentityContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::QuotientFieldEnergyIdentityError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::QuotientSideConditionContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::QuotientSideConditionError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::QuotientRecoveryRoundTripContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::QuotientRecoveryRoundTripError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::QuotientResponseRoundTripContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::QuotientBasisResponseRoundTripError,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::ReducedSymmetryContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::ReducedSymmetryDefect,
                observed: *observed,
                limit: *limit,
            });
        }
        RepresentationFailure::AffineReproductionBackend(failure) => {
            retain_kkt_failure(report, failure);
        }
        RepresentationFailure::AffineReproductionContract { observed, limit } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ContractThresholdExceeded {
                quantity: AnalysisContractQuantity::AffineReproductionError,
                observed: *observed,
                limit: *limit,
            });
        }
    }
}

fn public_interpretable_rank_deficiency(
    concept: RankDeficiencyConcept,
    domain: RankEvidenceDomain,
    source_relations: &[SourceHardRelation],
    canonical_mode_residual: f64,
    backend_invoked: bool,
    hidden_regularization_applied: bool,
) -> InterpretableRankDeficiencyEvidence {
    let mut provenance = source_relations
        .iter()
        .filter(|relation| relation.equality.field().is_some())
        .map(|relation| {
            (
                relation.equality.provenance().semantic_role().clone(),
                relation.equality.provenance().source().clone(),
            )
        })
        .collect::<Vec<_>>();
    provenance.sort();
    provenance.dedup();
    InterpretableRankDeficiencyEvidence::new(InterpretableRankDeficiencyEvidenceParts {
        concept,
        domain,
        source_ids: provenance
            .iter()
            .map(|(_, source)| source.clone())
            .collect(),
        semantic_roles: provenance.into_iter().map(|(role, _)| role).collect(),
        canonical_mode_residual,
        canonical_mode_verified: true,
        backend_invoked,
        hidden_regularization_applied,
    })
}

fn retain_kkt_failure(report: &mut FitReport, failure: &KktFailure) {
    report.attempts = public_attempts(kkt_failure_attempts(failure));
    report.execution_failure = public_kkt_failure(failure);
    match failure {
        KktFailure::Capacity(evidence) => {
            report.capacity = Some(public_capacity(evidence));
        }
        KktFailure::InvalidLength {
            field,
            expected,
            actual,
        } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::InvalidBackendInputLength {
                field: public_backend_input_field(*field),
                expected: *expected,
                actual: *actual,
            });
        }
        KktFailure::NonFiniteInput { field, index } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::NonFiniteBackendInput {
                field: public_backend_input_field(*field),
                index: *index,
            });
        }
        KktFailure::WorkspaceAllocation {
            phase,
            bytes,
            alignment,
        } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::WorkspaceAllocation {
                stage: public_backend_workspace_stage(*phase),
                bytes: *bytes,
                alignment: *alignment,
                backend_invoked: matches!(phase, InternalWorkspacePhase::SvdRescue),
            });
        }
        KktFailure::Scaling(reason) => {
            report.analysis_failure = Some(AnalysisFailureEvidence::ScalingFailure {
                reason: public_scaling_failure(*reason),
            });
        }
        KktFailure::RankDeficient { evidence }
        | KktFailure::NumericalDecisionGrayZone { evidence } => {
            report.backend_rank = Some(public_rank_evidence(
                RankEvidenceDomain::BackendKkt,
                report
                    .problem_size
                    .kkt_dimension
                    .expect("backend evidence implies an exact KKT dimension"),
                evidence,
            ));
        }
        KktFailure::UnexpectedInertia {
            expected,
            observed,
            backend_invoked,
        } => {
            report.inertia = Some(public_inertia(*expected, *observed, *backend_invoked));
        }
        KktFailure::AlgebraicAnalysisFailure { phase } => {
            report.analysis_failure = Some(AnalysisFailureEvidence::NumericalAnalysis {
                stage: public_backend_analysis_stage(*phase),
                backend_invoked: false,
            });
        }
        KktFailure::BackendContractViolation {
            analysis: Some(analysis),
            ..
        }
        | KktFailure::NumericalFailure {
            analysis: Some(analysis),
            ..
        } => {
            report.backend_rank = Some(public_rank_evidence(
                RankEvidenceDomain::BackendKkt,
                report
                    .problem_size
                    .kkt_dimension
                    .expect("backend evidence implies an exact KKT dimension"),
                &analysis.rank,
            ));
            report.inertia = Some(public_inertia(
                analysis.expected_inertia,
                analysis.observed_inertia,
                true,
            ));
        }
        KktFailure::BackendContractViolation { analysis: None, .. }
        | KktFailure::NumericalFailure { analysis: None, .. } => {}
    }
}

fn public_solve_coordinate_failure(
    reason: InternalSolveCoordinateFailure,
) -> SolveCoordinateFailureReason {
    match reason {
        InternalSolveCoordinateFailure::BoundingBoxCenterNotFinite => {
            SolveCoordinateFailureReason::BoundingBoxCenterNotFinite
        }
        InternalSolveCoordinateFailure::CharacteristicLengthNotFinite => {
            SolveCoordinateFailureReason::CharacteristicLengthNotFinite
        }
        InternalSolveCoordinateFailure::FieldRecoveryScaleNotInvertible => {
            SolveCoordinateFailureReason::FieldRecoveryScaleNotInvertible
        }
        InternalSolveCoordinateFailure::StandardFunctionalNotFinite => {
            SolveCoordinateFailureReason::StandardFunctionalNotFinite
        }
    }
}

fn public_cubic_analysis_stage(stage: InternalCubicAnalysisStage) -> AnalysisFailureStage {
    match stage {
        InternalCubicAnalysisStage::PolynomialRank => AnalysisFailureStage::CubicPolynomialRank,
        InternalCubicAnalysisStage::ReducedCholesky => AnalysisFailureStage::CubicReducedCholesky,
    }
}

fn public_backend_analysis_stage(stage: InternalBackendAnalysisPhase) -> AnalysisFailureStage {
    match stage {
        InternalBackendAnalysisPhase::RankConfirmation => {
            AnalysisFailureStage::BackendRankConfirmation
        }
        InternalBackendAnalysisPhase::Inertia => AnalysisFailureStage::BackendInertia,
    }
}

fn public_backend_workspace_stage(stage: InternalWorkspacePhase) -> AnalysisFailureStage {
    match stage {
        InternalWorkspacePhase::RankAnalysis => AnalysisFailureStage::BackendRankWorkspace,
        InternalWorkspacePhase::InertiaAnalysis => AnalysisFailureStage::BackendInertiaWorkspace,
        InternalWorkspacePhase::Factor => AnalysisFailureStage::BackendFactorWorkspace,
        InternalWorkspacePhase::Solve => AnalysisFailureStage::BackendSolveWorkspace,
        InternalWorkspacePhase::SvdRescue => AnalysisFailureStage::BackendSvdRescueWorkspace,
    }
}

fn public_backend_input_field(field: InternalBackendInputField) -> BackendInputField {
    match field {
        InternalBackendInputField::Hessian => BackendInputField::Hessian,
        InternalBackendInputField::EqualityJacobian => BackendInputField::EqualityJacobian,
        InternalBackendInputField::StationarityRightHandSide => {
            BackendInputField::StationarityRightHandSide
        }
        InternalBackendInputField::EqualityRightHandSide => {
            BackendInputField::EqualityRightHandSide
        }
    }
}

fn public_scaling_failure(reason: InternalScalingFailure) -> ScalingFailureReason {
    match reason {
        InternalScalingFailure::InvalidShape => ScalingFailureReason::InvalidShape,
        InternalScalingFailure::ZeroNorm { index } => ScalingFailureReason::ZeroNorm { index },
        InternalScalingFailure::NonFiniteNorm { index } => {
            ScalingFailureReason::NonFiniteNorm { index }
        }
    }
}

fn public_capacity(evidence: &CapacityExceededEvidence) -> CapacityEvidence {
    let (kind, planned_peak_bytes) = match &evidence.reason {
        CapacityExceededReason::ArithmeticOverflow { .. } => {
            (CapacityFailureKind::ArithmeticOverflow, None)
        }
        CapacityExceededReason::LimitExceeded {
            planned_peak_bytes, ..
        } => (
            CapacityFailureKind::LimitExceeded,
            Some(*planned_peak_bytes),
        ),
        CapacityExceededReason::ConvexQpLimitExceeded {
            planned_peak_bytes, ..
        } => (
            CapacityFailureKind::LimitExceeded,
            Some(*planned_peak_bytes),
        ),
    };
    CapacityEvidence::new(
        kind,
        evidence.limit_bytes,
        planned_peak_bytes,
        evidence.large_allocation_attempted,
        evidence.backend_invocation_attempted,
    )
}

fn public_failed_hard_relations(
    evidence: &RecoveryVerificationFailureEvidence,
    source_relations: &[SourceHardRelation],
) -> Vec<HardRelationAssessment> {
    let (Some(recovered), Some(tolerances)) =
        (&evidence.hard_equalities, &evidence.relation_tolerances)
    else {
        return Vec::new();
    };
    source_relations
        .iter()
        .map(|source_relation| {
            hard_relation_assessment(
                source_relation,
                recovered[source_relation.canonical_index].value,
                tolerances[source_relation.canonical_index],
            )
        })
        .collect()
}

fn public_failure_acceptance(
    evidence: &RecoveryVerificationFailureEvidence,
) -> Option<CanonicalAcceptanceEvidence> {
    Some(CanonicalAcceptanceEvidence::new(
        CanonicalAcceptanceEvidenceParts {
            accepted: false,
            backend_standard_form_verified: true,
            recovery_finite: evidence.recovery_finite?,
            provenance_verified: evidence.provenance_verified?,
            side_condition: evidence.side_condition.map(public_side_condition),
            hard_residual_maxima: evidence
                .hard_equality_violations
                .map(|envelope| (envelope.field_value, envelope.field_value_per_length)),
            polynomial_round_trip_error: evidence.polynomial_round_trip_error,
            field_coefficient_round_trip_error: evidence.field_coefficient_round_trip_error,
            field_energy_round_trip_error: evidence.field_energy_round_trip_error,
            whitening_round_trip_error: evidence.whitening_round_trip_error,
            objective_round_trip_error: evidence.objective_round_trip_error,
            objective_verified: evidence.objective_round_trip_error.is_some()
                && !evidence
                    .reasons
                    .contains(&crate::cubic_equality::RecoveryVerificationFailureReason::ObjectiveRoundTripViolation),
            tolerance_round_trip_error: evidence.tolerance_round_trip_error,
            hard_affine_inequality_violation_max: None,
            backend_standard_form_residual: None,
            physical_convex_residual: None,
            scaling_round_trip_error: None,
            reduction_round_trip_error: None,
            backend_internal_scaling_round_trip_error: None,
        },
    ))
}

fn kkt_failure_attempts(failure: &KktFailure) -> &[KktAttemptRecord] {
    match failure {
        KktFailure::BackendContractViolation { attempts, .. }
        | KktFailure::NumericalFailure { attempts, .. } => attempts,
        _ => &[],
    }
}

fn public_kkt_failure(failure: &KktFailure) -> Option<AttemptFailureEvidence> {
    match failure {
        KktFailure::BackendContractViolation { reason, .. } => {
            Some(public_backend_contract_failure(*reason))
        }
        KktFailure::NumericalFailure {
            reason: NumericalFailureReason::BackendDecompositionFailure,
            ..
        } => Some(AttemptFailureEvidence::new(
            AttemptFailureCategory::BackendDecompositionFailure,
            None,
            None,
        )),
        _ => None,
    }
}

fn public_attempts(attempts: &[KktAttemptRecord]) -> Vec<SolveAttemptRecord> {
    attempts
        .iter()
        .map(|attempt| {
            SolveAttemptRecord::new(SolveAttemptRecordParts {
                sequence: attempt.sequence,
                kind: match attempt.kind {
                    KktAttemptKind::BunchKaufmanRefinement => {
                        SolveAttemptKind::BunchKaufmanRefinement
                    }
                    KktAttemptKind::SvdRescue => SolveAttemptKind::FullSvdRescue,
                },
                termination: match attempt.termination {
                    InternalTermination::CandidateProduced => {
                        SolveAttemptTermination::CandidateProduced
                    }
                    InternalTermination::NumericalError => SolveAttemptTermination::NumericalError,
                },
                settings: public_attempt_settings(attempt.settings),
                scaling: ScalingSummary::new(
                    attempt.scaling.method,
                    attempt.scaling.rounds,
                    attempt.scaling.saturated_outside_target,
                ),
                scaling_round_trip_error: None,
                refinement_steps: attempt.refinement_steps,
                residual: attempt.residual.map(|residual| {
                    LinearResidualEvidence::new([
                        residual.infinity_norm,
                        residual.matrix_infinity_norm,
                        residual.solution_infinity_norm,
                        residual.rhs_infinity_norm,
                        residual.normalized_backward_error,
                    ])
                }),
                convex_residual: None,
                certificate_present: attempt.certificate_present,
                failure_reason: attempt.failure_reason.map(public_attempt_failure),
                backend_fingerprint: public_backend_fingerprint(&attempt.backend),
            })
        })
        .collect()
}

fn public_qp_attempts(attempts: &[QpAttemptRecord]) -> Vec<SolveAttemptRecord> {
    attempts
        .iter()
        .map(|attempt| {
            let kind = match attempt.backend.profile {
                ClarabelAttemptProfile::Standard => SolveAttemptKind::ClarabelStandard,
                ClarabelAttemptProfile::Robust => SolveAttemptKind::ClarabelRobust,
            };
            SolveAttemptRecord::new(SolveAttemptRecordParts {
                sequence: attempt.backend.sequence,
                kind,
                termination: match attempt.backend.termination {
                    ClarabelTermination::Solved => SolveAttemptTermination::CandidateProduced,
                    ClarabelTermination::AlmostSolved => {
                        SolveAttemptTermination::ReducedAccuracyCandidateProduced
                    }
                    ClarabelTermination::PrimalInfeasible
                    | ClarabelTermination::AlmostPrimalInfeasible => {
                        SolveAttemptTermination::PrimalInfeasibilityCandidate
                    }
                    ClarabelTermination::DualInfeasible
                    | ClarabelTermination::AlmostDualInfeasible => {
                        SolveAttemptTermination::DualInfeasibilityCandidate
                    }
                    ClarabelTermination::IterationLimit | ClarabelTermination::TimeLimit => {
                        SolveAttemptTermination::LimitReached
                    }
                    ClarabelTermination::InsufficientProgress => {
                        SolveAttemptTermination::InsufficientProgress
                    }
                    ClarabelTermination::CallbackTerminated => {
                        SolveAttemptTermination::CallbackTermination
                    }
                    ClarabelTermination::NumericalError | ClarabelTermination::Unsolved => {
                        SolveAttemptTermination::NumericalError
                    }
                },
                settings: BackendAttemptSettings::clarabel(
                    kind,
                    attempt.backend.settings.all_settings.clone(),
                ),
                scaling: ScalingSummary::new(
                    "georbf-v1-block-aware-ruiz-power-of-two",
                    attempt.georbf_scaling.rounds.len(),
                    attempt.georbf_scaling.saturated_outside_target,
                ),
                scaling_round_trip_error: Some(attempt.georbf_scaling_round_trip_error),
                refinement_steps: 0,
                residual: None,
                convex_residual: attempt.residuals.map(|residual| {
                    PublicConvexResidualEvidence::new(ConvexResidualEvidenceParts {
                        primal: residual.primal,
                        dual: residual.dual,
                        stationarity: residual.stationarity,
                        complementarity: residual.complementarity,
                        relative_gap: residual.relative_gap,
                    })
                }),
                certificate_present: matches!(
                    attempt.backend.termination,
                    ClarabelTermination::PrimalInfeasible
                        | ClarabelTermination::AlmostPrimalInfeasible
                        | ClarabelTermination::DualInfeasible
                        | ClarabelTermination::AlmostDualInfeasible
                ),
                failure_reason: attempt.failure_reason.map(public_qp_attempt_failure),
                backend_fingerprint: public_qp_backend_fingerprint(&attempt.backend),
            })
        })
        .collect()
}

fn public_qp_attempt_failure(reason: QpAttemptFailureReason) -> AttemptFailureEvidence {
    let category = match reason {
        QpAttemptFailureReason::NonFiniteCandidate => AttemptFailureCategory::NonFiniteCandidate,
        QpAttemptFailureReason::BackendResidualExceeded => {
            AttemptFailureCategory::ConvexResidualExceeded
        }
        QpAttemptFailureReason::ThreadContractViolation => {
            AttemptFailureCategory::ThreadContractViolation
        }
        QpAttemptFailureReason::BackendFingerprintMismatch => {
            AttemptFailureCategory::BackendFingerprintMismatch
        }
        QpAttemptFailureReason::UnverifiedTermination => {
            AttemptFailureCategory::UnverifiedTermination
        }
        QpAttemptFailureReason::InvalidInfeasibilityCertificate => {
            AttemptFailureCategory::InvalidInfeasibilityCertificate
        }
        QpAttemptFailureReason::InvalidRecessionCertificate => {
            AttemptFailureCategory::InvalidRecessionCertificate
        }
    };
    AttemptFailureEvidence::new(category, None, None)
}

fn public_qp_backend_fingerprint(backend: &ClarabelAttemptEvidence) -> BackendFingerprint {
    BackendFingerprint::new(BackendFingerprintParts {
        schema_version: 1,
        crate_name: backend.backend.crate_name,
        crate_version: backend.backend.crate_version,
        features: backend.backend.features.to_vec(),
        algorithm: "Clarabel-QP/qdldl",
        target_arch: std::env::consts::ARCH,
        target_os: std::env::consts::OS,
        requested_threads: backend.requested_threads,
        actual_threads: backend.actual_threads,
    })
}

fn public_backend_fingerprint(backend: &InternalBackendFingerprint) -> BackendFingerprint {
    BackendFingerprint::new(BackendFingerprintParts {
        schema_version: backend.schema_version,
        crate_name: backend.crate_name,
        crate_version: backend.crate_version,
        features: backend.features.to_vec(),
        algorithm: backend.algorithm,
        target_arch: backend.target_arch,
        target_os: backend.target_os,
        requested_threads: backend.requested_threads,
        actual_threads: backend.actual_threads,
    })
}

fn public_attempt_settings(settings: InternalAttemptSettings) -> BackendAttemptSettings {
    match settings {
        InternalAttemptSettings::Lblt {
            pivoting,
            block_size,
            parallelism_threshold,
            factor_workspace_source,
            maximum_refinement_steps,
        } => BackendAttemptSettings::lblt(
            pivoting,
            block_size,
            parallelism_threshold,
            factor_workspace_source,
            maximum_refinement_steps,
        ),
        InternalAttemptSettings::FullSvd {
            settings_id,
            left_vectors,
            right_vectors,
        } => BackendAttemptSettings::full_svd(settings_id, left_vectors, right_vectors),
    }
}

fn public_attempt_failure(reason: KktAttemptFailureReason) -> AttemptFailureEvidence {
    match reason {
        KktAttemptFailureReason::BackendContract(reason) => public_backend_contract_failure(reason),
        KktAttemptFailureReason::Numerical(NumericalFailureReason::BackendDecompositionFailure) => {
            AttemptFailureEvidence::new(
                AttemptFailureCategory::BackendDecompositionFailure,
                None,
                None,
            )
        }
    }
}

fn public_backend_contract_failure(
    reason: InternalBackendContractReason,
) -> AttemptFailureEvidence {
    match reason {
        InternalBackendContractReason::NonFiniteCandidate => {
            AttemptFailureEvidence::new(AttemptFailureCategory::NonFiniteCandidate, None, None)
        }
        InternalBackendContractReason::BackwardErrorExceeded { observed, limit } => {
            AttemptFailureEvidence::new(
                AttemptFailureCategory::BackwardErrorExceeded,
                Some(observed),
                Some(limit),
            )
        }
        InternalBackendContractReason::ScalingRoundTripExceeded { observed, limit } => {
            AttemptFailureEvidence::new(
                AttemptFailureCategory::ScalingRoundTripExceeded,
                Some(observed),
                Some(limit),
            )
        }
    }
}

fn public_cubic_analysis(evidence: &CpdEvidence) -> CubicAnalysisEvidence {
    CubicAnalysisEvidence::new(CubicAnalysisEvidenceParts {
        fitting_functional_count: evidence.fitting_functional_count,
        polynomial_dimension: evidence.polynomial_dimension,
        polynomial_rank: evidence.polynomial_rank,
        quotient_construction: CubicQuotientConstructionEvidence::new(
            CubicQuotientConstructionEvidenceParts {
                quotient_dimension: evidence.quotient_construction.quotient_dimension,
                householder_reflector_count: evidence
                    .quotient_construction
                    .householder_reflector_count,
                congruence_pass_count: evidence.quotient_construction.congruence_pass_count,
                householder_orthogonality_error: evidence
                    .quotient_construction
                    .householder_orthogonality_error,
                canonical_response_round_trip_error: evidence
                    .quotient_construction
                    .canonical_response_round_trip_error,
            },
        ),
        quotient_factorization: CubicQuotientFactorizationEvidence::new(
            CubicQuotientFactorizationEvidenceParts {
                quotient_dimension: evidence.quotient_factorization.quotient_dimension,
                retained_modes: evidence.quotient_factorization.retained_modes,
                truncated_modes: evidence.quotient_factorization.truncated_modes,
                unregularized_llt_count: evidence.quotient_factorization.unregularized_llt_count,
                full_spectrum_analysis_count: evidence
                    .quotient_factorization
                    .full_spectrum_analysis_count,
                normalized_backward_error: evidence
                    .quotient_factorization
                    .normalized_backward_error,
                pivot_intervals: evidence
                    .quotient_factorization
                    .pivot_intervals
                    .iter()
                    .map(|interval| {
                        CubicLltPivotInterval::new(CubicLltPivotIntervalParts {
                            lower: interval.lower,
                            upper: interval.upper,
                        })
                    })
                    .collect(),
                field_energy_identity_error: evidence
                    .quotient_factorization
                    .field_energy_identity_error,
                side_condition_error: evidence.quotient_factorization.side_condition_error,
                recovery_round_trip_error: evidence
                    .quotient_factorization
                    .recovery_round_trip_error,
                canonical_response_round_trip_error: evidence
                    .quotient_factorization
                    .canonical_response_round_trip_error,
                kernel_ridge_applied: evidence.quotient_factorization.kernel_ridge_applied,
                gram_jitter_applied: evidence.quotient_factorization.gram_jitter_applied,
                mode_truncation_applied: evidence.quotient_factorization.mode_truncation_applied,
            },
        ),
        polynomial_singular_values: evidence.singular_values.clone(),
        polynomial_rrqr_ratio: evidence.polynomial_rrqr_ratio,
        polynomial_svd_ratio: evidence.polynomial_svd_ratio,
        polynomial_rank_reject_ratio: evidence.polynomial_rank_reject_ratio,
        polynomial_rank_accept_ratio: evidence.polynomial_rank_accept_ratio,
        null_space_defect: evidence.quotient_construction.null_space_defect,
        reduced_symmetry_defect: evidence.reduced_symmetry_defect,
        reduced_symmetry_defect_limit: evidence.symmetry_defect_limit,
        reduced_largest_singular_value: evidence.reduced_largest_singular_value,
        reduced_smallest_singular_value: evidence.reduced_smallest_singular_value,
        affine_reproduction_error: evidence.affine_reproduction_error,
        solve_coordinate_length: evidence.solve_coordinate_length,
        degenerate_extent: evidence.degenerate_extent,
    })
}

fn public_rank_evidence(
    domain: RankEvidenceDomain,
    dimension: usize,
    evidence: &InternalRankEvidence,
) -> RankEvidence {
    let decision = match evidence.classification {
        InternalRankDecision::FullRank => RankDecision::FullRank,
        InternalRankDecision::RankDeficient => RankDecision::RankDeficient,
        InternalRankDecision::NumericalDecisionGrayZone => RankDecision::NumericalDecisionGrayZone,
    };
    RankEvidence::new(RankEvidenceParts {
        domain,
        dimension,
        rank: (decision == RankDecision::FullRank).then_some(dimension),
        exact_zero_index: evidence.exact_zero_index,
        rrqr_ratio: Some(evidence.rrqr_ratio),
        singular_values: evidence.singular_values.clone(),
        svd_ratio: Some(evidence.svd_ratio),
        reject_ratio: Some(evidence.reject_ratio),
        accept_ratio: Some(evidence.accept_ratio),
        decision,
        backend_invoked: evidence.backend_invoked,
    })
}

fn public_inertia(
    expected: InternalInertia,
    observed: InternalInertia,
    backend_invoked: bool,
) -> InertiaEvidence {
    InertiaEvidence::new(
        InertiaCounts::new(expected.positive, expected.negative, expected.zero),
        InertiaCounts::new(observed.positive, observed.negative, observed.zero),
        backend_invoked,
    )
}

fn public_side_condition(evidence: PhysicalSideConditionEvidence) -> SideConditionEvidence {
    SideConditionEvidence::new(
        evidence.components,
        evidence.physical_tolerances,
        evidence.round_trip_error,
    )
}

fn public_success_acceptance(solution: &CubicEqualitySolution) -> CanonicalAcceptanceEvidence {
    CanonicalAcceptanceEvidence::new(CanonicalAcceptanceEvidenceParts {
        accepted: true,
        backend_standard_form_verified: true,
        recovery_finite: solution.recovery_finite,
        provenance_verified: solution.provenance_verified,
        side_condition: Some(public_side_condition(solution.side_condition)),
        hard_residual_maxima: Some((
            solution.hard_equality_violations.field_value,
            solution.hard_equality_violations.field_value_per_length,
        )),
        polynomial_round_trip_error: Some(solution.polynomial_round_trip_error),
        field_coefficient_round_trip_error: Some(solution.field_coefficient_round_trip_error),
        field_energy_round_trip_error: Some(solution.field_energy_round_trip_error),
        whitening_round_trip_error: Some(solution.whitening_round_trip_error),
        objective_round_trip_error: Some(solution.objective_round_trip_error),
        objective_verified: solution.objective_verified,
        tolerance_round_trip_error: Some(solution.tolerance_round_trip_error),
        hard_affine_inequality_violation_max: None,
        backend_standard_form_residual: None,
        physical_convex_residual: None,
        scaling_round_trip_error: None,
        reduction_round_trip_error: None,
        backend_internal_scaling_round_trip_error: None,
    })
}

fn public_recovery_evidence(
    evidence: &RecoveryVerificationFailureEvidence,
) -> RecoveryVerificationEvidence {
    RecoveryVerificationEvidence::new(RecoveryVerificationEvidenceParts {
        reasons: evidence.reasons.clone(),
        side_condition: evidence.side_condition.map(public_side_condition),
        hard_residual_maxima: evidence
            .hard_equality_violations
            .map(|envelope| (envelope.field_value, envelope.field_value_per_length)),
        polynomial_round_trip_error: evidence.polynomial_round_trip_error,
        field_coefficient_round_trip_error: evidence.field_coefficient_round_trip_error,
        field_energy_round_trip_error: evidence.field_energy_round_trip_error,
        whitening_round_trip_error: evidence.whitening_round_trip_error,
        objective_round_trip_error: evidence.objective_round_trip_error,
        tolerance_round_trip_error: evidence.tolerance_round_trip_error,
        backend_standard_form_residual: None,
        reduction_round_trip_error: None,
        scaling_round_trip_error: None,
        sources: evidence
            .hard_equalities
            .as_deref()
            .into_iter()
            .flatten()
            .map(|relation| &relation.provenance)
            .chain(
                evidence
                    .soft_equalities
                    .as_deref()
                    .into_iter()
                    .flatten()
                    .map(|relation| &relation.provenance),
            )
            .map(public_canonical_evidence_source)
            .collect(),
        no_model_produced: evidence.no_model_produced,
    })
}

fn diagnose(failure: &CubicEqualityFailure) -> ProblemDiagnosis {
    match failure {
        CubicEqualityFailure::EmptyEqualitySet
        | CubicEqualityFailure::AffineInequalityRequiresConvexQp
        | CubicEqualityFailure::NonFiniteTarget { .. } => ProblemDiagnosis::InvalidProblem,
        CubicEqualityFailure::Representation(failure) => diagnose_representation(failure),
        CubicEqualityFailure::Backend { failure, .. } => diagnose_kkt(failure),
        CubicEqualityFailure::RecoveryVerification { .. } => {
            ProblemDiagnosis::RecoveryVerificationFailure
        }
    }
}

fn diagnose_qp(failure: &CubicExecutionFailure) -> ProblemDiagnosis {
    match failure {
        CubicExecutionFailure::Equality(failure) => diagnose(failure),
        CubicExecutionFailure::Capacity(_) => ProblemDiagnosis::CapacityExceeded,
        CubicExecutionFailure::Representation(failure) => diagnose_representation(failure),
        CubicExecutionFailure::BackendContract { .. } => ProblemDiagnosis::BackendContractViolation,
        CubicExecutionFailure::RecoveryVerification { .. } => {
            ProblemDiagnosis::RecoveryVerificationFailure
        }
        CubicExecutionFailure::ValidatedInfeasible { .. } => ProblemDiagnosis::InfeasibleProblem,
        CubicExecutionFailure::ValidatedUnbounded { .. } => ProblemDiagnosis::UnboundedProblem,
        CubicExecutionFailure::InconsistentAttempts { .. } => {
            ProblemDiagnosis::NumericalConsistencyFailure
        }
        CubicExecutionFailure::Assembly(_)
        | CubicExecutionFailure::BackendAdapter(_)
        | CubicExecutionFailure::AttemptsExhausted { .. } => ProblemDiagnosis::NumericalFailure,
    }
}

fn qp_failure_report(
    mut report: FitReport,
    failure: &CubicExecutionFailure,
    source_relations: &[SourceHardRelation],
) -> FitReport {
    match failure {
        CubicExecutionFailure::Equality(failure) => {
            report = failure_report(report, failure, source_relations);
        }
        CubicExecutionFailure::Capacity(evidence) => {
            report.capacity = Some(public_capacity(evidence));
        }
        CubicExecutionFailure::Representation(failure) => {
            retain_representation_failure(&mut report, failure, source_relations);
        }
        CubicExecutionFailure::Assembly(_) => {}
        CubicExecutionFailure::BackendAdapter(_) => {
            report.execution_failure = Some(AttemptFailureEvidence::new(
                AttemptFailureCategory::BackendDecompositionFailure,
                None,
                None,
            ));
        }
        CubicExecutionFailure::BackendContract {
            attempts,
            observed,
            limit,
        } => {
            report.attempts = public_qp_attempts(attempts);
            report.execution_failure = Some(
                attempts
                    .last()
                    .and_then(|attempt| attempt.failure_reason)
                    .map(|reason| match reason {
                        QpAttemptFailureReason::BackendResidualExceeded => {
                            AttemptFailureEvidence::new(
                                AttemptFailureCategory::ConvexResidualExceeded,
                                Some(*observed),
                                Some(*limit),
                            )
                        }
                        other => public_qp_attempt_failure(other),
                    })
                    .unwrap_or_else(|| {
                        AttemptFailureEvidence::new(
                            AttemptFailureCategory::ConvexResidualExceeded,
                            Some(*observed),
                            Some(*limit),
                        )
                    }),
            );
            if let Some(last) = attempts.last() {
                report.backend_fingerprint = Some(public_qp_backend_fingerprint(&last.backend));
            }
        }
        CubicExecutionFailure::AttemptsExhausted { attempts } => {
            report.attempts = public_qp_attempts(attempts);
            report.execution_failure = attempts
                .last()
                .and_then(|attempt| attempt.failure_reason)
                .map(public_qp_attempt_failure);
            if let Some(last) = attempts.last() {
                report.backend_fingerprint = Some(public_qp_backend_fingerprint(&last.backend));
            }
        }
        CubicExecutionFailure::RecoveryVerification { evidence, attempts } => {
            report.recovery_verification = Some(public_qp_recovery_evidence(evidence));
            report.attempts = public_qp_attempts(attempts);
            if let Some(last) = attempts.last() {
                report.backend_fingerprint = Some(public_qp_backend_fingerprint(&last.backend));
            }
        }
        CubicExecutionFailure::ValidatedInfeasible { evidence, attempts } => {
            report.infeasibility_certificate = Some(public_infeasibility_certificate(evidence));
            report.attempts = public_qp_attempts(attempts);
            if let Some(last) = attempts.last() {
                report.backend_fingerprint = Some(public_qp_backend_fingerprint(&last.backend));
            }
        }
        CubicExecutionFailure::ValidatedUnbounded { evidence, attempts } => {
            report.recession_ray = Some(public_recession_ray(evidence));
            report.attempts = public_qp_attempts(attempts);
            if let Some(last) = attempts.last() {
                report.backend_fingerprint = Some(public_qp_backend_fingerprint(&last.backend));
            }
        }
        CubicExecutionFailure::InconsistentAttempts { attempts } => {
            report.attempts = public_qp_attempts(attempts);
            report.execution_failure = Some(AttemptFailureEvidence::new(
                AttemptFailureCategory::InconsistentValidatedConclusions,
                None,
                None,
            ));
            if let Some(last) = attempts.last() {
                report.backend_fingerprint = Some(public_qp_backend_fingerprint(&last.backend));
            }
        }
    }
    report
}

fn public_qp_recovery_evidence(
    evidence: &crate::cubic_execution::QpRecoveryFailureEvidence,
) -> RecoveryVerificationEvidence {
    use crate::cubic_equality::RecoveryVerificationFailureReason as PublicReason;
    use crate::cubic_execution::QpRecoveryFailureReason as InternalReason;

    let reasons = evidence
        .reasons
        .iter()
        .map(|reason| match reason {
            InternalReason::InvalidRecoveryMap => PublicReason::InvalidRecoveryMap,
            InternalReason::ProvenanceMismatch => PublicReason::ProvenanceMismatch,
            InternalReason::NonFiniteRecoveredQuantity => PublicReason::NonFiniteRecoveredQuantity,
            InternalReason::SideConditionViolation => PublicReason::SideConditionViolation,
            InternalReason::SideConditionRoundTripViolation => {
                PublicReason::SideConditionRoundTripViolation
            }
            InternalReason::HardEqualityViolation => PublicReason::HardEqualityViolation,
            InternalReason::AffineInequalityViolation => PublicReason::AffineInequalityViolation,
            InternalReason::BackendSlackMismatch => PublicReason::BackendSlackMismatch,
            InternalReason::ReductionRoundTripViolation => {
                PublicReason::ReductionRoundTripViolation
            }
            InternalReason::ScalingRoundTripViolation => PublicReason::ScalingRoundTripViolation,
            InternalReason::PolynomialRoundTripViolation => {
                PublicReason::PolynomialRoundTripViolation
            }
            InternalReason::FieldCoefficientRoundTripViolation => {
                PublicReason::FieldCoefficientRoundTripViolation
            }
            InternalReason::FieldEnergyRoundTripViolation => {
                PublicReason::FieldEnergyRoundTripViolation
            }
            InternalReason::WhiteningRoundTripViolation => {
                PublicReason::WhiteningRoundTripViolation
            }
            InternalReason::ObjectiveRoundTripViolation => {
                PublicReason::ObjectiveRoundTripViolation
            }
        })
        .collect();
    RecoveryVerificationEvidence::new(RecoveryVerificationEvidenceParts {
        reasons,
        side_condition: None,
        hard_residual_maxima: None,
        polynomial_round_trip_error: None,
        field_coefficient_round_trip_error: None,
        field_energy_round_trip_error: None,
        whitening_round_trip_error: None,
        objective_round_trip_error: None,
        tolerance_round_trip_error: None,
        backend_standard_form_residual: Some(evidence.backend_standard_form_violation),
        reduction_round_trip_error: Some(evidence.reduction_round_trip_error),
        scaling_round_trip_error: Some(evidence.scaling_round_trip_error),
        sources: evidence
            .sources
            .iter()
            .map(public_canonical_evidence_source)
            .collect(),
        no_model_produced: evidence.no_model_produced,
    })
}

fn public_infeasibility_certificate(
    evidence: &ValidatedInfeasibilityEvidence,
) -> InfeasibilityCertificateEvidence {
    InfeasibilityCertificateEvidence::new(InfeasibilityCertificateEvidenceParts {
        finite: evidence.finite,
        normalized_ray_norm: evidence.normalized_ray_norm,
        stationarity_residual: evidence.stationarity_residual,
        dual_cone_violation: evidence.dual_cone_violation,
        separation_margin: evidence.separation_margin,
        residual_limit: evidence.residual_limit,
        separation_limit: evidence.separation_limit,
        recovery_round_trip_error: evidence.recovery_round_trip_error,
        provenance_verified: evidence.provenance_verified,
        sources: evidence
            .sources
            .iter()
            .map(public_canonical_evidence_source)
            .collect(),
        backend_invoked: true,
    })
}

fn public_recession_ray(evidence: &ValidatedRecessionEvidence) -> RecessionRayEvidence {
    RecessionRayEvidence::new(RecessionRayEvidenceParts {
        finite: evidence.finite,
        normalized_ray_norm: evidence.normalized_ray_norm,
        hessian_null_residual: evidence.hessian_null_residual,
        constraint_ray_violation: evidence.constraint_ray_violation,
        descent_margin: evidence.descent_margin,
        residual_limit: evidence.residual_limit,
        separation_limit: evidence.separation_limit,
        recovery_round_trip_error: evidence.recovery_round_trip_error,
        provenance_verified: evidence.provenance_verified,
        sources: evidence
            .sources
            .iter()
            .map(public_canonical_evidence_source)
            .collect(),
        backend_invoked: true,
    })
}

fn public_canonical_evidence_source(provenance: &UsageProvenance) -> CanonicalEvidenceSource {
    CanonicalEvidenceSource::new(
        provenance.source().clone(),
        provenance.groups().to_vec(),
        provenance.semantic_role().clone(),
    )
}

fn diagnose_representation(failure: &RepresentationFailure) -> ProblemDiagnosis {
    match failure {
        RepresentationFailure::Capacity(_) => ProblemDiagnosis::CapacityExceeded,
        RepresentationFailure::PolynomialRankDeficient { .. } => {
            ProblemDiagnosis::UnidentifiedFieldMode
        }
        RepresentationFailure::PolynomialRankGrayZone { .. } => {
            ProblemDiagnosis::NumericalDecisionGrayZone
        }
        RepresentationFailure::AffineReproductionBackend(failure) => diagnose_kkt(failure),
        _ => ProblemDiagnosis::NumericalFailure,
    }
}

fn diagnose_kkt(failure: &KktFailure) -> ProblemDiagnosis {
    match failure {
        KktFailure::Capacity(_) => ProblemDiagnosis::CapacityExceeded,
        KktFailure::RankDeficient { .. } => ProblemDiagnosis::NumericalFailure,
        KktFailure::NumericalDecisionGrayZone { .. } => ProblemDiagnosis::NumericalDecisionGrayZone,
        KktFailure::BackendContractViolation { .. } => ProblemDiagnosis::BackendContractViolation,
        _ => ProblemDiagnosis::NumericalFailure,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic_equality::{RecoveryVerificationFailureReason, inject_kkt_failure_once};
    use crate::cubic_execution::{QpFaultInjection, inject_qp_fault_once};
    use crate::geometry::{Handedness, InputCoordinateFrame, LengthUnitLabel, Vector3};
    use crate::kkt::{EqualityKktSystem, solve_equality_kkt};
    use crate::observation::{
        DirectedNormalObservation, FieldValueObservation, MinimumNormalSlope,
    };
    use crate::relation::{
        DirectionalDerivativeInterval, FieldValueBound, MinimumFieldSeparation,
        SharedLevelSetBuilder, StratigraphicFieldDirection, YoungerThan,
    };
    use crate::{Point3, ProblemBuilder};

    fn injectable_snapshot() -> ProblemSnapshot {
        let mut builder = ProblemBuilder::new(
            InputCoordinateFrame::try_new(
                ["x", "y", "z"],
                Handedness::Right,
                LengthUnitLabel::new("m"),
            )
            .unwrap(),
            FieldUnitLabel::new("field"),
        );
        for (source, support, value) in [
            ("origin", [0.0, 0.0, 0.0], 1.0),
            ("east", [1.0, 0.0, 0.0], 2.0),
            ("north", [0.0, 1.0, 0.0], 3.0),
            ("up", [0.0, 0.0, 1.0], 4.0),
        ] {
            builder
                .add(
                    FieldValueObservation::try_new(
                        SourceId::new(source),
                        Point3::try_new(support[0], support[1], support[2]).unwrap(),
                        value,
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        builder.build().unwrap()
    }

    fn injectable_bounded_snapshot() -> ProblemSnapshot {
        let snapshot = injectable_snapshot();
        let mut builder = ProblemBuilder::new(
            snapshot.input_coordinate_frame().clone(),
            snapshot.field_unit().clone(),
        );
        for observation in &snapshot.inner.observations {
            match observation {
                ObservationInput::FieldValue(observation) => {
                    builder.add(observation.clone()).unwrap();
                }
                _ => unreachable!("the injectable fixture contains only field values"),
            }
        }
        builder
            .add(
                FieldValueBound::try_upper(
                    SourceId::new("bound"),
                    Point3::try_new(0.5, 0.5, 0.5).unwrap(),
                    10.0,
                )
                .unwrap(),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn injectable_normal_snapshot() -> ProblemSnapshot {
        let snapshot = injectable_snapshot();
        let mut builder = ProblemBuilder::new(
            snapshot.input_coordinate_frame().clone(),
            snapshot.field_unit().clone(),
        );
        for observation in &snapshot.inner.observations {
            match observation {
                ObservationInput::FieldValue(observation) => {
                    builder.add(observation.clone()).unwrap();
                }
                _ => unreachable!("the injectable fixture contains only field values"),
            }
        }
        builder
            .add(
                DirectedNormalObservation::try_new(
                    SourceId::new("normal"),
                    Point3::try_new(0.25, 0.25, 0.25).unwrap(),
                    Vector3::try_new(1.0, 2.0, 2.0).unwrap(),
                    MinimumNormalSlope::try_new(0.1).unwrap(),
                )
                .unwrap(),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn derivative_interval_snapshot(upper_bound: f64) -> ProblemSnapshot {
        let snapshot = injectable_snapshot();
        let mut builder = ProblemBuilder::new(
            snapshot.input_coordinate_frame().clone(),
            snapshot.field_unit().clone(),
        );
        for observation in &snapshot.inner.observations {
            match observation {
                ObservationInput::FieldValue(observation) => {
                    builder.add(observation.clone()).unwrap();
                }
                _ => unreachable!("the injectable fixture contains only field values"),
            }
        }
        builder
            .add(
                DirectionalDerivativeInterval::try_upper(
                    SourceId::new("derivative-interval"),
                    Point3::try_new(0.5, 0.5, 0.5).unwrap(),
                    Vector3::try_new(1.0, 2.0, 2.0).unwrap(),
                    upper_bound,
                )
                .unwrap(),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn injectable_derivative_interval_snapshot() -> ProblemSnapshot {
        derivative_interval_snapshot(10.0)
    }

    fn injectable_shared_level_set_relation_snapshot() -> ProblemSnapshot {
        let snapshot = injectable_snapshot();
        let mut builder = ProblemBuilder::new(
            snapshot.input_coordinate_frame().clone(),
            snapshot.field_unit().clone(),
        );
        builder
            .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
            .unwrap();
        for observation in &snapshot.inner.observations {
            match observation {
                ObservationInput::FieldValue(observation) => {
                    builder.add(observation.clone()).unwrap();
                }
                _ => unreachable!("the injectable fixture contains only field values"),
            }
        }
        for (group_id, source_id, support) in [
            ("older", "older/member", [0.0, 0.0, 0.0]),
            ("younger", "younger/member", [1.0, 0.0, 0.0]),
        ] {
            let mut group = SharedLevelSetBuilder::new(GroupId::new(group_id));
            group
                .add_member(
                    SourceId::new(source_id),
                    Point3::try_new(support[0], support[1], support[2]).unwrap(),
                )
                .unwrap();
            builder.add(group.build().unwrap()).unwrap();
        }
        builder
            .add(YoungerThan::hard(
                SourceId::new("age"),
                GroupId::new("younger"),
                GroupId::new("older"),
                MinimumFieldSeparation::try_new(0.5).unwrap(),
            ))
            .unwrap();
        builder.build().unwrap()
    }

    fn injectable_infeasible_derivative_interval_snapshot() -> ProblemSnapshot {
        injectable_infeasible_snapshot(true)
    }

    fn injectable_infeasible_bounded_snapshot() -> ProblemSnapshot {
        injectable_infeasible_snapshot(false)
    }

    fn injectable_infeasible_snapshot(include_derivative_interval: bool) -> ProblemSnapshot {
        let mut builder = ProblemBuilder::new(
            InputCoordinateFrame::try_new(
                ["x", "y", "z"],
                Handedness::Right,
                LengthUnitLabel::new("m"),
            )
            .unwrap(),
            FieldUnitLabel::new("field"),
        );
        let mut level = SharedLevelSetBuilder::new(GroupId::new("one-level"));
        for (source, support) in [
            ("member/origin", [0.0, 0.0, 0.0]),
            ("member/east", [1.0, 0.0, 0.0]),
            ("member/north", [0.0, 1.0, 0.0]),
            ("member/up", [0.0, 0.0, 1.0]),
        ] {
            level
                .add_member(
                    SourceId::new(source),
                    Point3::try_new(support[0], support[1], support[2]).unwrap(),
                )
                .unwrap();
        }
        builder.add(level.build().unwrap()).unwrap();
        builder
            .add(
                FieldValueBound::try_lower(
                    SourceId::new("lower"),
                    Point3::try_new(0.0, 0.0, 0.0).unwrap(),
                    2.0,
                )
                .unwrap(),
            )
            .unwrap();
        builder
            .add(
                FieldValueBound::try_upper(
                    SourceId::new("upper"),
                    Point3::try_new(1.0, 0.0, 0.0).unwrap(),
                    1.0,
                )
                .unwrap(),
            )
            .unwrap();
        if include_derivative_interval {
            builder
                .add(
                    DirectionalDerivativeInterval::try_interval(
                        SourceId::new("derivative-interval"),
                        Point3::try_new(0.5, 0.5, 0.5).unwrap(),
                        Vector3::try_new(1.0, 2.0, 2.0).unwrap(),
                        0.0,
                        0.0,
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        builder.build().unwrap()
    }

    #[test]
    fn snapshot_capacity_counts_source_relations_as_linear_report_storage() {
        let lowering = EqualityLowering::new(Vec::new());
        assert!(plan_snapshot_capacity(&lowering, 0, 10_000, 0).is_ok());
        assert!(plan_source_lifecycle_capacity(4_000_000, 0).is_ok());
        let lifecycle_evidence = plan_source_lifecycle_capacity(4_194_304, 0)
            .expect_err("source/report lifecycle storage is guarded before evidence cloning");
        assert!(!lifecycle_evidence.large_allocation_attempted);
        assert!(!lifecycle_evidence.backend_invocation_attempted);
        assert!(plan_source_preflight_capacity(4_194_304, 0).is_ok());
        let preflight_evidence = plan_source_preflight_capacity(8_388_608, 0)
            .expect_err("oversized canonical lowering is rejected before materialization");
        assert!(!preflight_evidence.large_allocation_attempted);
        assert!(!preflight_evidence.backend_invocation_attempted);
        let evidence = plan_snapshot_capacity(&lowering, 0, usize::MAX, 0)
            .expect_err("report storage arithmetic remains checked before allocation");
        assert!(!evidence.large_allocation_attempted);
        assert!(!evidence.backend_invocation_attempted);
        assert!(
            ProblemSize::cubic_equality(CubicProblemSizeParts {
                input_observations: 0,
                scalar_hard_relations: 0,
                scalar_soft_relations: 0,
                canonical_hard_equalities: 0,
                canonical_soft_equalities: 0,
                quadratic_objective_terms: 0,
                linear_objective_terms: 0,
                affine_inequality_constraints: 0,
                center_coefficients: usize::MAX,
                semantic_latents: 1,
                solver_hard_equalities: 0,
                auxiliary_variables: 0,
            })
            .is_none()
        );
    }

    #[test]
    fn incidence_identity_canonicalizes_signed_zero_coordinates() {
        assert_eq!(
            canonical_support_bits([-0.0, 1.0, 0.0]),
            canonical_support_bits([0.0, 1.0, -0.0])
        );
    }

    #[test]
    fn recovery_mapping_retains_quantified_side_condition_evidence() {
        let evidence = RecoveryVerificationFailureEvidence {
            reasons: vec![RecoveryVerificationFailureReason::SideConditionViolation],
            side_condition: Some(PhysicalSideConditionEvidence {
                components: [1.0, 2.0, 3.0, 4.0],
                physical_tolerances: [0.1, 0.2, 0.3, 0.4],
                standard_components: [5.0, 6.0, 7.0, 8.0],
                recovered_standard_components: [5.0, 6.0, 7.0, 8.0],
                round_trip_error: 9.0e-12,
            }),
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
            recovery_finite: Some(true),
            provenance_verified: Some(true),
            no_model_produced: true,
        };

        let public = public_recovery_evidence(&evidence);
        let side = public
            .side_condition()
            .expect("the complete physical side-condition evidence survives mapping");
        assert_eq!(side.components(), [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(side.tolerances(), [0.1, 0.2, 0.3, 0.4]);
        assert_eq!(side.round_trip_error(), 9.0e-12);
    }

    #[test]
    fn public_bound_fit_rejects_recovery_corruption_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::RecoveryMap);
        let failure = injectable_bounded_snapshot()
            .fit()
            .expect_err("a damaged QP recovery map must not publish a model");

        assert_eq!(
            failure.diagnosis(),
            ProblemDiagnosis::RecoveryVerificationFailure
        );
        assert!(failure.report().canonical_acceptance().is_none());
        let recovery = failure
            .report()
            .recovery_verification()
            .expect("the public report retains QP recovery evidence");
        assert_eq!(
            recovery.reasons(),
            &[RecoveryVerificationFailureReason::ReductionRoundTripViolation]
        );
        assert!(recovery.backend_standard_form_residual().is_some());
        assert!(recovery.reduction_round_trip_error().unwrap() > 1.0e-11);
        assert!(recovery.scaling_round_trip_error().is_some());
        assert!(
            recovery
                .sources()
                .iter()
                .any(|source| source.source_id().as_str() == "bound")
        );
        assert!(recovery.no_model_produced());
        assert!(!failure.report().attempts().is_empty());
    }

    #[test]
    fn public_normal_fit_rejects_recovery_corruption_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::RecoveryMap);
        let failure = injectable_normal_snapshot()
            .fit()
            .expect_err("a damaged Normal recovery map must not publish a model");

        assert_eq!(
            failure.diagnosis(),
            ProblemDiagnosis::RecoveryVerificationFailure
        );
        assert!(failure.report().canonical_acceptance().is_none());
        assert!(
            failure
                .report()
                .recovery_verification()
                .unwrap()
                .no_model_produced()
        );
    }

    #[test]
    fn public_normal_fit_rejects_backend_contract_corruption_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::BackendResidual);
        let failure = injectable_normal_snapshot()
            .fit()
            .expect_err("a damaged Normal backend residual must not publish a model");

        assert_eq!(
            failure.diagnosis(),
            ProblemDiagnosis::BackendContractViolation
        );
        assert!(failure.report().canonical_acceptance().is_none());
    }

    #[test]
    fn normal_qp_capacity_failure_is_structured_before_allocation() {
        inject_qp_fault_once(QpFaultInjection::Capacity);
        let failure = injectable_normal_snapshot()
            .fit()
            .expect_err("an oversized Normal QP plan must fail before allocation");

        assert_eq!(failure.diagnosis(), ProblemDiagnosis::CapacityExceeded);
        let capacity = failure.report().capacity().unwrap();
        assert!(!capacity.large_allocation_attempted());
        assert!(!capacity.backend_invocation_attempted());
        assert!(failure.report().attempts().is_empty());
    }

    #[test]
    fn public_derivative_interval_fit_rejects_recovery_corruption_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::RecoveryMap);
        let failure = injectable_derivative_interval_snapshot()
            .fit()
            .expect_err("a damaged derivative QP recovery map must not publish a model");

        assert_eq!(
            failure.diagnosis(),
            ProblemDiagnosis::RecoveryVerificationFailure
        );
        assert!(failure.report().canonical_acceptance().is_none());
        let recovery = failure
            .report()
            .recovery_verification()
            .expect("the report retains derivative-interval recovery evidence");
        assert_eq!(
            recovery.reasons(),
            &[RecoveryVerificationFailureReason::ReductionRoundTripViolation]
        );
        assert!(recovery.no_model_produced());
        assert!(!failure.report().attempts().is_empty());
    }

    #[test]
    fn derivative_interval_qp_capacity_failure_is_structured_before_allocation() {
        inject_qp_fault_once(QpFaultInjection::Capacity);
        let failure = injectable_derivative_interval_snapshot()
            .fit()
            .expect_err("an oversized derivative QP plan must not attempt allocation");

        assert_eq!(failure.diagnosis(), ProblemDiagnosis::CapacityExceeded);
        let capacity = failure
            .report()
            .capacity()
            .expect("the derivative QP capacity failure retains checked evidence");
        assert!(!capacity.large_allocation_attempted());
        assert!(!capacity.backend_invocation_attempted());
        assert!(failure.report().attempts().is_empty());
    }

    #[test]
    fn shared_level_set_relation_fit_rejects_recovery_corruption_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::ScalingMap);
        let failure = injectable_shared_level_set_relation_snapshot()
            .fit()
            .expect_err("a damaged shared-level recovery map must not publish a model");

        assert_eq!(
            failure.diagnosis(),
            ProblemDiagnosis::RecoveryVerificationFailure
        );
        assert!(failure.report().canonical_acceptance().is_none());
        let recovery = failure.report().recovery_verification().unwrap();
        assert_eq!(
            recovery.reasons(),
            &[RecoveryVerificationFailureReason::ScalingRoundTripViolation]
        );
        assert!(recovery.no_model_produced());
    }

    #[test]
    fn shared_level_set_relation_qp_capacity_failure_is_preallocation_evidence() {
        inject_qp_fault_once(QpFaultInjection::Capacity);
        let failure = injectable_shared_level_set_relation_snapshot()
            .fit()
            .expect_err("an oversized shared-level QP plan must fail before allocation");

        assert_eq!(failure.diagnosis(), ProblemDiagnosis::CapacityExceeded);
        let capacity = failure.report().capacity().unwrap();
        assert!(!capacity.large_allocation_attempted());
        assert!(!capacity.backend_invocation_attempted());
        assert!(failure.report().attempts().is_empty());
    }

    #[test]
    fn invalid_farkas_candidates_remain_a_numerical_failure_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::InfeasibilityCertificate);
        let failure = injectable_infeasible_bounded_snapshot()
            .fit()
            .expect_err("unvalidated backend rays must not prove infeasibility");

        assert_eq!(failure.diagnosis(), ProblemDiagnosis::NumericalFailure);
        assert!(failure.report().infeasibility_certificate().is_none());
        assert_eq!(failure.report().attempts().len(), 2);
        assert!(failure.report().attempts().iter().all(|attempt| {
            attempt.certificate_present()
                && attempt.failure_reason().is_some_and(|evidence| {
                    evidence.category() == AttemptFailureCategory::InvalidInfeasibilityCertificate
                })
        }));
    }

    #[test]
    fn invalid_derivative_farkas_candidates_remain_unverified_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::InfeasibilityCertificate);
        let failure = injectable_infeasible_derivative_interval_snapshot()
            .fit()
            .expect_err("unvalidated derivative-bound rays must not prove infeasibility");

        assert_eq!(failure.diagnosis(), ProblemDiagnosis::NumericalFailure);
        assert!(failure.report().infeasibility_certificate().is_none());
        assert_eq!(failure.report().attempts().len(), 2);
        assert!(failure.report().attempts().iter().all(|attempt| {
            attempt.certificate_present()
                && attempt.failure_reason().is_some_and(|evidence| {
                    evidence.category() == AttemptFailureCategory::InvalidInfeasibilityCertificate
                })
        }));
    }

    #[test]
    fn validated_recession_evidence_maps_to_unbounded_problem_without_a_model() {
        let evidence = crate::cubic_execution::ValidatedRecessionEvidence {
            finite: true,
            normalized_ray_norm: 1.0,
            hessian_null_residual: 0.0,
            constraint_ray_violation: 0.0,
            descent_margin: 1.0,
            residual_limit: 1.0e-8,
            separation_limit: 1.0e-7,
            recovery_round_trip_error: 0.0,
            provenance_verified: true,
            sources: vec![UsageProvenance::new(
                SourceId::new("unbounded-source"),
                None,
                RelationId::new("unbounded-relation"),
                ResidualId::new("unbounded-residual"),
                SemanticRolePath::new("recession/objective"),
            )],
        };
        let failure = CubicExecutionFailure::ValidatedUnbounded {
            evidence,
            attempts: Vec::new(),
        };
        assert_eq!(diagnose_qp(&failure), ProblemDiagnosis::UnboundedProblem);

        let snapshot = injectable_bounded_snapshot();
        let report = qp_failure_report(
            empty_report(
                &snapshot,
                conservative_problem_size(0, ScalarRelationCounts { hard: 0, soft: 0 }, 0, 0, 0, 0),
            ),
            &failure,
            &[],
        );
        let ray = report
            .recession_ray()
            .expect("the report exposes independently validated recession evidence");
        assert!(ray.finite());
        assert_eq!(ray.normalized_ray_norm(), 1.0);
        assert_eq!(ray.hessian_null_residual(), 0.0);
        assert_eq!(ray.constraint_ray_violation(), 0.0);
        assert_eq!(ray.descent_margin(), 1.0);
        assert!(ray.provenance_verified());
        assert_eq!(ray.recovery_round_trip_error(), 0.0);
        assert_eq!(ray.sources()[0].source_id().as_str(), "unbounded-source");
        assert!(report.canonical_acceptance().is_none());
    }

    #[test]
    fn almost_solved_remains_a_reduced_accuracy_candidate_under_identical_acceptance() {
        inject_qp_fault_once(QpFaultInjection::AlmostSolved);
        let success = injectable_bounded_snapshot()
            .fit()
            .expect("an almost-solved candidate that passes every canonical check is accepted");

        assert!(success.report().canonical_acceptance().unwrap().accepted());
        assert_eq!(success.report().attempts().len(), 1);
        assert_eq!(
            success.report().attempts()[0].termination(),
            SolveAttemptTermination::ReducedAccuracyCandidateProduced
        );
        let residual = success.report().attempts()[0].convex_residual().unwrap();
        assert!(residual.primal() <= 1.0e-8);
        assert!(residual.dual() <= 1.0e-8);
        assert!(residual.stationarity() <= 1.0e-8);
        assert!(residual.complementarity() <= 1.0e-8);
        assert!(residual.relative_gap() <= 1.0e-8);
        let physical_residual = success
            .report()
            .canonical_acceptance()
            .unwrap()
            .physical_convex_residual()
            .unwrap();
        assert!(physical_residual.primal() <= 1.0e-8);
        assert!(physical_residual.dual() <= 1.0e-8);
        assert!(physical_residual.stationarity() <= 1.0e-8);
        assert!(physical_residual.complementarity() <= 1.0e-8);
        assert!(physical_residual.relative_gap() <= 1.0e-8);
    }

    #[test]
    fn limit_terminations_exhaust_the_fixed_plan_without_a_candidate_or_diagnosis_upgrade() {
        inject_qp_fault_once(QpFaultInjection::Limit);
        let failure = injectable_bounded_snapshot()
            .fit()
            .expect_err("limit terminations must never publish a partial model");

        assert_eq!(failure.diagnosis(), ProblemDiagnosis::NumericalFailure);
        assert!(failure.report().canonical_acceptance().is_none());
        assert!(failure.report().infeasibility_certificate().is_none());
        assert!(failure.report().recession_ray().is_none());
        assert_eq!(failure.report().attempts().len(), 2);
        assert!(failure.report().attempts().iter().all(|attempt| {
            attempt.termination() == SolveAttemptTermination::LimitReached
                && attempt.failure_reason().is_some_and(|evidence| {
                    evidence.category() == AttemptFailureCategory::UnverifiedTermination
                })
        }));
    }

    #[test]
    fn contradictory_validated_attempts_map_to_a_distinct_consistency_failure() {
        let failure = CubicExecutionFailure::InconsistentAttempts {
            attempts: Vec::new(),
        };
        assert_eq!(
            diagnose_qp(&failure),
            ProblemDiagnosis::NumericalConsistencyFailure
        );
        let snapshot = injectable_bounded_snapshot();
        let report = qp_failure_report(
            empty_report(
                &snapshot,
                conservative_problem_size(0, ScalarRelationCounts { hard: 0, soft: 0 }, 0, 0, 0, 0),
            ),
            &failure,
            &[],
        );
        assert_eq!(
            report.execution_failure().unwrap().category(),
            AttemptFailureCategory::InconsistentValidatedConclusions
        );
        assert!(report.canonical_acceptance().is_none());
        assert!(report.infeasibility_certificate().is_none());
        assert!(report.recession_ray().is_none());
    }

    #[test]
    fn problem_snapshot_fit_maps_injected_backend_rejection_independently_of_termination() {
        let failure = solve_equality_kkt(&EqualityKktSystem {
            primal_variables: 1,
            equality_constraints: 0,
            hessian: &[1.0e-320],
            equality_jacobian: &[],
            stationarity_rhs: &[1.0],
            equality_rhs: &[],
        })
        .expect_err("non-finite backend candidates must exhaust the bounded plan");

        assert_eq!(
            diagnose_kkt(&failure),
            ProblemDiagnosis::BackendContractViolation
        );
        inject_kkt_failure_once(failure);
        let fit_failure = injectable_snapshot()
            .fit()
            .expect_err("the injected backend rejection must prevent model publication");
        assert_eq!(
            fit_failure.diagnosis(),
            ProblemDiagnosis::BackendContractViolation
        );
        assert!(fit_failure.report().canonical_acceptance().is_none());
        let attempts = fit_failure.report().attempts();
        assert_eq!(attempts.len(), 2);
        assert!(attempts.iter().all(|attempt| {
            attempt.termination() == SolveAttemptTermination::CandidateProduced
                && attempt.failure_reason().is_some_and(|evidence| {
                    evidence.category() == AttemptFailureCategory::NonFiniteCandidate
                })
        }));
    }

    #[test]
    fn uninterpreted_backend_rank_loss_remains_a_numerical_failure() {
        let failure = KktFailure::RankDeficient {
            evidence: InternalRankEvidence {
                exact_zero_index: Some(0),
                rrqr_ratio: 0.0,
                singular_values: Vec::new(),
                svd_ratio: 0.0,
                reject_ratio: 1.0e-14,
                accept_ratio: 1.0e-12,
                classification: InternalRankDecision::RankDeficient,
                backend_invoked: false,
            },
        };

        assert_eq!(diagnose_kkt(&failure), ProblemDiagnosis::NumericalFailure);
        inject_kkt_failure_once(failure);
        let fit_failure = injectable_snapshot()
            .fit()
            .expect_err("an uninterpreted KKT rank loss must not publish a model");
        assert_eq!(fit_failure.diagnosis(), ProblemDiagnosis::NumericalFailure);
        assert_eq!(
            fit_failure.report().rank_evidence().unwrap().domain(),
            RankEvidenceDomain::BackendKkt
        );
        assert!(
            fit_failure
                .report()
                .interpretable_rank_deficiency()
                .is_none()
        );

        let reduced_failure = RepresentationFailure::QuotientPivotRequiresPrecisionRescue {
            quotient_dimension: 2,
            pivot_index: 1,
            interval: Some(crate::cubic_equality::OutwardRoundedInterval {
                lower: -f64::from_bits(1),
                upper: f64::from_bits(1),
            }),
            execution: crate::cubic_equality::AnalysisExecutionEvidence {
                solver_invoked: false,
                hidden_regularization_applied: false,
            },
        };
        assert_eq!(
            diagnose_representation(&reduced_failure),
            ProblemDiagnosis::NumericalFailure
        );

        let snapshot = injectable_snapshot();
        let problem_size =
            conservative_problem_size(6, ScalarRelationCounts { hard: 6, soft: 0 }, 0, 0, 0, 0);
        let mut report = empty_report(&snapshot, problem_size);
        retain_representation_failure(&mut report, &reduced_failure, &[]);
        assert!(report.rank_evidence().is_none());
        match report.analysis_failure().unwrap() {
            AnalysisFailureEvidence::QuotientPivotRequiresPrecisionRescue {
                quotient_dimension,
                pivot_index,
                interval: Some(interval),
                backend_invoked,
            } => {
                assert_eq!(*quotient_dimension, 2);
                assert_eq!(*pivot_index, 1);
                assert!(interval.lower_bound() < 0.0);
                assert!(interval.upper_bound() > 0.0);
                assert!(!*backend_invoked);
            }
            other => panic!("unexpected public quotient failure evidence: {other:?}"),
        }

        let nonpositive_failure = RepresentationFailure::QuotientFactorizationNotPositive {
            quotient_dimension: 2,
            pivot_index: 1,
            interval: crate::cubic_equality::OutwardRoundedInterval {
                lower: -1.0,
                upper: -1.0,
            },
            execution: crate::cubic_equality::AnalysisExecutionEvidence {
                solver_invoked: false,
                hidden_regularization_applied: false,
            },
        };
        assert_eq!(
            diagnose_representation(&nonpositive_failure),
            ProblemDiagnosis::NumericalFailure
        );
        retain_representation_failure(&mut report, &nonpositive_failure, &[]);
        assert!(matches!(
            report.analysis_failure(),
            Some(AnalysisFailureEvidence::QuotientFactorizationNotPositive {
                quotient_dimension: 2,
                pivot_index: 1,
                ..
            })
        ));
    }
}
