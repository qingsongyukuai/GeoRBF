//! Strict fit outcomes and physical-unit fit reports.

use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fmt;

use crate::capacity::{CapacityExceededEvidence, CapacityExceededReason, plan_equality_capacity};
use crate::cubic_equality::{
    AlgebraicAnalysisStage as InternalCubicAnalysisStage, CanonicalEqualityParticipation,
    CanonicalHardEquality, CanonicalRelationToleranceEvidence, CpdEvidence, CubicCanonicalProblem,
    CubicEqualityCore, CubicEqualityFailure, CubicEqualitySolution, PhysicalSideConditionEvidence,
    RecoveryVerificationFailureEvidence, ReducedPairingFailureClassification,
    RepresentationFailure, SemanticLatentCoefficient, SemanticLatentDefinition,
    SolveCoordinateTransformFailureReason as InternalSolveCoordinateFailure,
};
use crate::diagnostics::{
    AnalysisContractQuantity, AnalysisFailureEvidence, AnalysisFailureStage,
    AttemptFailureCategory, AttemptFailureEvidence, BackendAttemptSettings, BackendFingerprint,
    BackendFingerprintParts, BackendInputField, CanonicalAcceptanceEvidence,
    CanonicalAcceptanceEvidenceParts, CapacityEvidence, CapacityFailureKind, CubicAnalysisEvidence,
    CubicAnalysisEvidenceParts, DirectInputConflictEvidence, InertiaCounts, InertiaEvidence,
    LinearResidualEvidence, ProblemDiagnosis, RankDecision, RankEvidence, RankEvidenceDomain,
    RankEvidenceParts, RecoveryVerificationEvidence, RecoveryVerificationEvidenceParts,
    RelationGraphConflictEvidence, ResidualDimension, ScalingFailureReason, ScalingSummary,
    SideConditionEvidence, SolveAttemptKind, SolveAttemptRecord, SolveAttemptRecordParts,
    SolveAttemptTermination, SolveCoordinateFailureReason, UnidentifiedAdditiveGaugeEvidence,
    UninformativeSharedLevelSetEvidence,
};
use crate::functional::{
    CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, GroupId, RelationId,
    ResidualId, SemanticRolePath, SourceId, UsageProvenance,
};
use crate::geometry::FieldUnitLabel;
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
use crate::observation::ObservationInput;
use crate::problem::{ProblemSnapshot, ThreadBudget};
use crate::relation::AdditiveFieldGaugeReference;

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
    canonical_hard_equalities: Option<usize>,
    center_coefficients: Option<usize>,
    semantic_latents: usize,
    auxiliary_variables: usize,
    cone_blocks: usize,
    primal_variables: Option<usize>,
    equality_constraints: Option<usize>,
    kkt_dimension: Option<usize>,
}

impl ProblemSize {
    fn cubic_equality(
        input_observations: usize,
        scalar_hard_relations: usize,
        canonical_hard_equalities: usize,
        center_coefficients: usize,
        semantic_latents: usize,
        solver_hard_equalities: usize,
    ) -> Option<Self> {
        let primal_variables = center_coefficients
            .checked_add(4)?
            .checked_add(semantic_latents)?;
        let equality_constraints = solver_hard_equalities.checked_add(4)?;
        let kkt_dimension = primal_variables.checked_add(equality_constraints)?;
        Some(Self {
            input_observations,
            scalar_hard_relations,
            canonical_hard_equalities: Some(canonical_hard_equalities),
            center_coefficients: Some(center_coefficients),
            semantic_latents,
            auxiliary_variables: 0,
            cone_blocks: 0,
            primal_variables: Some(primal_variables),
            equality_constraints: Some(equality_constraints),
            kkt_dimension: Some(kkt_dimension),
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

    /// Returns the exact hard equalities retained after duplicate merging, or
    /// `None` when fit stopped before canonical lowering.
    pub fn canonical_hard_equalities(self) -> Option<usize> {
        self.canonical_hard_equalities
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
    shared_level_values: Vec<SharedLevelValue>,
    field_energy: Option<f64>,
    total_objective: Option<f64>,
    backend_fingerprint: Option<BackendFingerprint>,
    attempts: Vec<SolveAttemptRecord>,
    recovery_verification: Option<RecoveryVerificationEvidence>,
    direct_input_conflict: Option<DirectInputConflictEvidence>,
    relation_graph_conflict: Option<RelationGraphConflictEvidence>,
    execution_failure: Option<AttemptFailureEvidence>,
    cubic_analysis: Option<CubicAnalysisEvidence>,
    backend_rank: Option<RankEvidence>,
    inertia: Option<InertiaEvidence>,
    canonical_acceptance: Option<CanonicalAcceptanceEvidence>,
    capacity: Option<CapacityEvidence>,
    analysis_failure: Option<AnalysisFailureEvidence>,
    unidentified_additive_gauge: Option<UnidentifiedAdditiveGaugeEvidence>,
    uninformative_shared_level_set: Option<UninformativeSharedLevelSetEvidence>,
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

    /// Returns the resolved dimensionless FieldEnergy normalization.
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
        self.direct_input_conflict.as_ref()
    }

    /// Returns complete path provenance for a hard relation-graph conflict.
    pub fn relation_graph_conflict(&self) -> Option<&RelationGraphConflictEvidence> {
        self.relation_graph_conflict.as_ref()
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

    /// Returns structural evidence for a missing additive-field representative.
    pub fn unidentified_additive_gauge(&self) -> Option<&UnidentifiedAdditiveGaugeEvidence> {
        self.unidentified_additive_gauge.as_ref()
    }

    /// Returns structural evidence for a disconnected one-member shared level set.
    pub fn uninformative_shared_level_set(&self) -> Option<&UninformativeSharedLevelSetEvidence> {
        self.uninformative_shared_level_set.as_ref()
    }
}

pub(crate) fn fit_snapshot(snapshot: &ProblemSnapshot) -> Result<FitSuccess, FitFailure> {
    let scalar_relation_count = scalar_relation_count(snapshot).unwrap_or(usize::MAX);
    let conservative_problem_size = conservative_problem_size(
        snapshot.inner.observations.len(),
        scalar_relation_count,
        snapshot.inner.shared_level_sets.len(),
    );
    if let Some(group) = snapshot.inner.shared_level_sets.iter().find(|group| {
        group.members().len() == 1 && !snapshot_references_group(snapshot, group.group_id())
    }) {
        let mut report = empty_report(snapshot, conservative_problem_size);
        report.uninformative_shared_level_set = Some(UninformativeSharedLevelSetEvidence::new(
            group.group_id().clone(),
            group.members()[0].source_id().clone(),
            false,
        ));
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::UninformativeSharedLevelSet,
            report: Box::new(report),
        });
    }
    let has_absolute_reference = snapshot
        .inner
        .observations
        .iter()
        .any(|observation| matches!(observation, ObservationInput::FieldValue(_)))
        || !snapshot.inner.additive_field_gauges.is_empty();
    if !has_absolute_reference {
        let mut source_ids = snapshot
            .inner
            .observations
            .iter()
            .map(|observation| observation.source_id().clone())
            .chain(
                snapshot
                    .inner
                    .shared_level_sets
                    .iter()
                    .flat_map(|group| group.members())
                    .map(|member| member.source_id().clone()),
            )
            .collect::<Vec<_>>();
        source_ids.sort();
        let group_ids = snapshot
            .inner
            .shared_level_sets
            .iter()
            .map(|group| group.group_id().clone())
            .collect();
        let mut report = empty_report(snapshot, conservative_problem_size);
        report.unidentified_additive_gauge = Some(UnidentifiedAdditiveGaugeEvidence::new(
            source_ids, group_ids, false,
        ));
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::UnidentifiedAdditiveGauge,
            report: Box::new(report),
        });
    }
    if let Err(evidence) = plan_snapshot_capacity(
        scalar_relation_count,
        snapshot.inner.shared_level_sets.len(),
    ) {
        let mut report = empty_report(snapshot, conservative_problem_size);
        report.capacity = Some(public_capacity(&evidence));
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::CapacityExceeded,
            report: Box::new(report),
        });
    }
    let lowering = lower_snapshot(snapshot);
    let problem_size = ProblemSize::cubic_equality(
        snapshot.inner.observations.len(),
        scalar_relation_count,
        lowering.canonical_equalities.len(),
        lowering.fitting_functional_count(),
        lowering.semantic_latents.len(),
        lowering.solver_equality_count(),
    )
    .expect("the conservative snapshot capacity plan proved exact dimensions representable");
    let base_report = || empty_report(snapshot, problem_size);
    if lowering.direct_input_conflict.is_some() || lowering.relation_graph_conflict.is_some() {
        let mut report = base_report();
        report.direct_input_conflict = lowering.direct_input_conflict.clone();
        report.relation_graph_conflict = lowering.relation_graph_conflict.clone();
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::DirectInputConflict,
            report: Box::new(report),
        });
    }
    fit_snapshot_after_preflight(snapshot, lowering, problem_size, base_report)
}

fn empty_report(snapshot: &ProblemSnapshot, problem_size: ProblemSize) -> FitReport {
    FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: snapshot.inner.fit_configuration.numerical_policy(),
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations: Vec::new(),
        shared_level_values: Vec::new(),
        field_energy: None,
        total_objective: None,
        backend_fingerprint: None,
        attempts: Vec::new(),
        recovery_verification: None,
        direct_input_conflict: None,
        relation_graph_conflict: None,
        execution_failure: None,
        cubic_analysis: None,
        backend_rank: None,
        inertia: None,
        canonical_acceptance: None,
        capacity: None,
        analysis_failure: None,
        unidentified_additive_gauge: None,
        uninformative_shared_level_set: None,
    }
}

fn scalar_relation_count(snapshot: &ProblemSnapshot) -> Option<usize> {
    snapshot
        .inner
        .observations
        .iter()
        .try_fold(0_usize, |count, observation| {
            count.checked_add(match observation {
                ObservationInput::FieldValue(_) => 1,
                ObservationInput::Gradient(_) => 3,
            })
        })?
        .checked_add(
            snapshot
                .inner
                .shared_level_sets
                .iter()
                .try_fold(0_usize, |count, group| {
                    count.checked_add(group.members().len())
                })?,
        )?
        .checked_add(snapshot.inner.additive_field_gauges.len())
}

fn snapshot_references_group(snapshot: &ProblemSnapshot, group_id: &GroupId) -> bool {
    snapshot.inner.additive_field_gauges.iter().any(|gauge| {
        matches!(
            gauge.reference(),
            AdditiveFieldGaugeReference::LevelSet(referenced) if referenced == group_id
        )
    })
}

fn conservative_problem_size(
    input_observations: usize,
    scalar_relations: usize,
    semantic_latents: usize,
) -> ProblemSize {
    ProblemSize {
        input_observations,
        scalar_hard_relations: scalar_relations,
        canonical_hard_equalities: None,
        center_coefficients: None,
        semantic_latents,
        auxiliary_variables: 0,
        cone_blocks: 0,
        primal_variables: None,
        equality_constraints: None,
        kkt_dimension: None,
    }
}

fn plan_snapshot_capacity(
    scalar_relations: usize,
    semantic_latents: usize,
) -> Result<(), CapacityExceededEvidence> {
    let Some(primal_variables) = scalar_relations
        .checked_add(4)
        .and_then(|value| value.checked_add(semantic_latents))
    else {
        return Err(plan_equality_capacity(usize::MAX, usize::MAX)
            .expect_err("maximal dimensions must overflow the capacity plan"));
    };
    let Some(equality_constraints) = scalar_relations.checked_add(4) else {
        return Err(plan_equality_capacity(usize::MAX, usize::MAX)
            .expect_err("maximal dimensions must overflow the capacity plan"));
    };
    plan_equality_capacity(primal_variables, equality_constraints).map(|_| ())
}

fn fit_snapshot_after_preflight(
    snapshot: &ProblemSnapshot,
    lowering: EqualityLowering,
    problem_size: ProblemSize,
    base_report: impl Fn() -> FitReport,
) -> Result<FitSuccess, FitFailure> {
    let solution = match CubicEqualityCore::solve_canonical(
        CubicCanonicalProblem {
            equalities: lowering.canonical_equalities.clone(),
            semantic_latents: lowering.semantic_latents.clone(),
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

#[derive(Debug, Clone)]
struct ScalarObservation {
    source_id: SourceId,
    group_id: Option<GroupId>,
    support: [f64; 3],
    component: DirectInputComponent,
    semantic_role: SemanticRolePath,
    target: f64,
}

impl ScalarObservation {
    fn dimension(&self) -> FunctionalDimension {
        match self.component {
            DirectInputComponent::FieldValue => FunctionalDimension::FieldValue,
            DirectInputComponent::Gradient(_) => FunctionalDimension::FieldValuePerLength,
        }
    }

    fn value_coefficient(&self) -> f64 {
        match self.component {
            DirectInputComponent::FieldValue => 1.0,
            DirectInputComponent::Gradient(_) => 0.0,
        }
    }

    fn gradient_coefficient(&self) -> [f64; 3] {
        match self.component {
            DirectInputComponent::FieldValue => [0.0; 3],
            DirectInputComponent::Gradient(axis) => {
                std::array::from_fn(|component| if component == axis { 1.0 } else { 0.0 })
            }
        }
    }
}

fn scalar_observations(observations: &[ObservationInput]) -> Vec<ScalarObservation> {
    observations
        .iter()
        .flat_map(|observation| match observation {
            ObservationInput::FieldValue(observation) => vec![ScalarObservation {
                source_id: observation.source_id().clone(),
                group_id: None,
                support: observation.location().components(),
                component: DirectInputComponent::FieldValue,
                semantic_role: SemanticRolePath::new("field-value-observation/value"),
                target: observation.value(),
            }],
            ObservationInput::Gradient(observation) => observation
                .gradient()
                .components()
                .into_iter()
                .enumerate()
                .map(|(axis, target)| ScalarObservation {
                    source_id: observation.source_id().clone(),
                    group_id: None,
                    support: observation.location().components(),
                    component: DirectInputComponent::Gradient(axis),
                    semantic_role: SemanticRolePath::new(format!(
                        "gradient-observation/component/{axis}"
                    )),
                    target,
                })
                .collect(),
        })
        .collect()
}

#[derive(Debug, Clone)]
struct SourceHardRelation {
    equality: CanonicalHardEquality,
    canonical_index: usize,
}

#[derive(Debug, Clone)]
struct EqualityLowering {
    source_relations: Vec<SourceHardRelation>,
    canonical_equalities: Vec<CanonicalHardEquality>,
    canonical_index_by_key: BTreeMap<CanonicalEqualityKey, (usize, f64)>,
    direct_input_conflict: Option<DirectInputConflictEvidence>,
    relation_graph_conflict: Option<RelationGraphConflictEvidence>,
    semantic_latents: Vec<SemanticLatentDefinition>,
}

impl EqualityLowering {
    fn new() -> Self {
        Self {
            source_relations: Vec::new(),
            canonical_equalities: Vec::new(),
            canonical_index_by_key: BTreeMap::new(),
            direct_input_conflict: None,
            relation_graph_conflict: None,
            semantic_latents: Vec::new(),
        }
    }

    fn push_source(&mut self, equality: CanonicalHardEquality) {
        let (key, normalized_target) = normalized_equality_key(&equality);
        let canonical_index =
            if let Some((index, first_target)) = self.canonical_index_by_key.get(&key).copied() {
                if first_target != normalized_target && self.direct_input_conflict.is_none() {
                    let first = &self.canonical_equalities[index];
                    self.direct_input_conflict = Some(DirectInputConflictEvidence::new(
                        first.provenance().source().clone(),
                        equality.provenance().source().clone(),
                        equality.provenance().semantic_role().clone(),
                        first.target(),
                        equality.target(),
                    ));
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
        if self.relation_graph_conflict.is_none() {
            self.relation_graph_conflict = Some(RelationGraphConflictEvidence::new(
                proof_source_ids.clone(),
                proof_group_ids.clone(),
                semantic_role.clone(),
                first_absolute.source_id.clone(),
                first_absolute.target,
                second_absolute.source_id.clone(),
                second_absolute.target,
            ));
        }
    }

    fn solver_equality_count(&self) -> usize {
        self.canonical_equalities
            .iter()
            .filter(|equality| {
                equality.participation() == CanonicalEqualityParticipation::SolverConstraint
            })
            .count()
    }

    fn fitting_functional_count(&self) -> usize {
        let mut functionals = Vec::<CanonicalFunctional>::new();
        for functional in self
            .canonical_equalities
            .iter()
            .filter(|equality| {
                equality.participation() == CanonicalEqualityParticipation::SolverConstraint
            })
            .filter_map(CanonicalHardEquality::field)
            .map(FunctionalUse::functional)
        {
            if !functionals.iter().any(|existing| existing == functional) {
                functionals.push(functional.clone());
            }
        }
        functionals.len()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CanonicalEqualityKey {
    dimension: u8,
    field_terms: Vec<([u64; 3], u64, [u64; 3])>,
    latent_coefficients: Vec<(usize, u64)>,
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
    support.map(|coordinate| {
        if coordinate == 0.0 {
            0.0_f64.to_bits()
        } else {
            coordinate.to_bits()
        }
    })
}

fn normalized_equality_key(equality: &CanonicalHardEquality) -> (CanonicalEqualityKey, f64) {
    let first_coefficient = equality
        .field()
        .into_iter()
        .flat_map(|field| field.functional().terms())
        .flat_map(|term| {
            std::iter::once(term.value_coefficient()).chain(term.gradient_coefficient())
        })
        .chain(
            equality
                .latent_coefficients()
                .iter()
                .map(|term| term.coefficient),
        )
        .find(|coefficient| *coefficient != 0.0)
        .expect("a canonical equality has a field or semantic-latent coefficient");
    let sign = if first_coefficient.is_sign_negative() {
        -1.0
    } else {
        1.0
    };
    let dimension = match equality.dimension() {
        FunctionalDimension::FieldValue => 0,
        FunctionalDimension::FieldValuePerLength => 1,
    };
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

fn lower_snapshot(snapshot: &ProblemSnapshot) -> EqualityLowering {
    let mut lowering = EqualityLowering::new();
    let mut value_constraints = CanonicalValueConstraintForest::default();
    for observation in scalar_observations(&snapshot.inner.observations) {
        let edge = (observation.component == DirectInputComponent::FieldValue).then(|| {
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

    let has_absolute_observation = snapshot
        .inner
        .observations
        .iter()
        .any(|observation| matches!(observation, ObservationInput::FieldValue(_)));
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
                        component: DirectInputComponent::FieldValue,
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
    lowering
}

fn field_equality(
    observation: &ScalarObservation,
    participation: CanonicalEqualityParticipation,
) -> CanonicalHardEquality {
    let provenance = relation_provenance(
        observation.source_id.clone(),
        observation.group_id.clone(),
        observation.semantic_role.clone(),
    );
    let functional = CanonicalFunctional::new(
        observation.dimension(),
        vec![FunctionalTerm::new(
            observation.support,
            observation.value_coefficient(),
            observation.gradient_coefficient(),
        )],
    )
    .expect("checked public observations lower to a finite nonzero functional");
    CanonicalHardEquality::new(
        Some(FunctionalUse::new(functional, provenance.clone())),
        Vec::new(),
        provenance,
        observation.dimension(),
        observation.target,
        participation,
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

fn signed_coefficient_bits(sign: f64, coefficient: f64) -> u64 {
    let value = sign * coefficient;
    if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    }
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
        canonical_hard_equalities: Some(solution.assembly.canonical_hard_equalities),
        center_coefficients: Some(solution.assembly.field_coefficients),
        semantic_latents: solution.assembly.semantic_latents,
        auxiliary_variables: 0,
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
    let backend_fingerprint = public_backend_fingerprint(&solution.backend.backend);
    let attempts = public_attempts(&solution.backend.attempts);
    FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: solution.backend.numerical_policy,
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations,
        shared_level_values: shared_level_values.to_vec(),
        field_energy: Some(solution.field_energy),
        total_objective: Some(solution.total_objective),
        backend_fingerprint: Some(backend_fingerprint),
        attempts,
        recovery_verification: None,
        direct_input_conflict: None,
        relation_graph_conflict: None,
        execution_failure: None,
        cubic_analysis: Some(public_cubic_analysis(&solution.representation)),
        backend_rank: Some(public_rank_evidence(
            RankEvidenceDomain::BackendKkt,
            solution.backend.capacity.kkt_dimension,
            &solution.backend.rank,
        )),
        inertia: Some(public_inertia(
            solution.backend.expected_inertia,
            solution.backend.observed_inertia,
            true,
        )),
        canonical_acceptance: Some(public_success_acceptance(solution)),
        capacity: None,
        analysis_failure: None,
        unidentified_additive_gauge: None,
        uninformative_shared_level_set: None,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DirectInputComponent {
    FieldValue,
    Gradient(usize),
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
            retain_representation_failure(&mut report, failure);
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
            report.recovery_verification = Some(public_recovery_evidence(evidence));
        }
        CubicEqualityFailure::EmptyEqualitySet | CubicEqualityFailure::NonFiniteTarget { .. } => {}
    }
    report
}

fn retain_representation_failure(report: &mut FitReport, failure: &RepresentationFailure) {
    let reduced_dimension = report
        .problem_size
        .canonical_hard_equalities
        .expect("representation evidence implies canonical lowering completed")
        .saturating_sub(4);
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
        RepresentationFailure::PolynomialRankDeficient { rank, .. } => {
            report.backend_rank = Some(RankEvidence::new(RankEvidenceParts {
                domain: RankEvidenceDomain::CubicPolynomialPairing,
                dimension: 4,
                rank: Some(*rank),
                exact_zero_index: None,
                rrqr_ratio: None,
                singular_values: Vec::new(),
                svd_ratio: None,
                reject_ratio: None,
                accept_ratio: None,
                decision: RankDecision::RankDeficient,
                backend_invoked: false,
            }));
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
        RepresentationFailure::ReducedPairingGrayZone { solver_invoked, .. } => {
            report.backend_rank = Some(RankEvidence::new(RankEvidenceParts {
                domain: RankEvidenceDomain::CubicReducedPairing,
                dimension: reduced_dimension,
                rank: None,
                exact_zero_index: None,
                rrqr_ratio: None,
                singular_values: Vec::new(),
                svd_ratio: None,
                reject_ratio: None,
                accept_ratio: None,
                decision: RankDecision::NumericalDecisionGrayZone,
                backend_invoked: *solver_invoked,
            }));
        }
        RepresentationFailure::ReducedPairingNotPositive {
            classification,
            rank,
            solver_invoked,
            ..
        } => {
            report.backend_rank = Some(RankEvidence::new(RankEvidenceParts {
                domain: RankEvidenceDomain::CubicReducedPairing,
                dimension: reduced_dimension,
                rank: Some(*rank),
                exact_zero_index: None,
                rrqr_ratio: None,
                singular_values: Vec::new(),
                svd_ratio: None,
                reject_ratio: None,
                accept_ratio: None,
                decision: match classification {
                    ReducedPairingFailureClassification::RankDeficient => {
                        RankDecision::RankDeficient
                    }
                    ReducedPairingFailureClassification::NegativeCurvature => {
                        RankDecision::FullRank
                    }
                },
                backend_invoked: *solver_invoked,
            }));
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
        InternalCubicAnalysisStage::ReducedInertia => AnalysisFailureStage::CubicReducedInertia,
        InternalCubicAnalysisStage::ReducedSpectrum => AnalysisFailureStage::CubicReducedSpectrum,
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
            recovery_finite: evidence.recovery_finite?,
            provenance_verified: evidence.provenance_verified?,
            side_condition: evidence.side_condition.map(public_side_condition),
            hard_residual_maxima: evidence
                .hard_equality_violations
                .map(|envelope| (envelope.field_value, envelope.field_value_per_length)),
            polynomial_round_trip_error: evidence.polynomial_round_trip_error,
            field_coefficient_round_trip_error: evidence.field_coefficient_round_trip_error,
            field_energy_round_trip_error: evidence.field_energy_round_trip_error,
            tolerance_round_trip_error: evidence.tolerance_round_trip_error,
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
                    InternalTermination::AcceptedCandidate => {
                        SolveAttemptTermination::AcceptedCandidate
                    }
                    InternalTermination::RejectedCandidate => {
                        SolveAttemptTermination::RejectedCandidate
                    }
                    InternalTermination::NumericalError => SolveAttemptTermination::NumericalError,
                },
                settings: public_attempt_settings(attempt.settings),
                scaling: ScalingSummary::new(
                    attempt.scaling.method,
                    attempt.scaling.rounds,
                    attempt.scaling.saturated_outside_target,
                ),
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
                certificate_present: attempt.certificate_present,
                failure_reason: attempt.failure_reason.map(public_attempt_failure),
                backend_fingerprint: public_backend_fingerprint(&attempt.backend),
            })
        })
        .collect()
}

fn public_backend_fingerprint(backend: &InternalBackendFingerprint) -> BackendFingerprint {
    BackendFingerprint::new(BackendFingerprintParts {
        schema_version: backend.schema_version,
        crate_name: backend.crate_name,
        crate_version: backend.crate_version,
        features: backend.features,
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
        polynomial_singular_values: evidence.singular_values.clone(),
        polynomial_rrqr_ratio: evidence.polynomial_rrqr_ratio,
        polynomial_svd_ratio: evidence.polynomial_svd_ratio,
        polynomial_rank_reject_ratio: evidence.polynomial_rank_reject_ratio,
        polynomial_rank_accept_ratio: evidence.polynomial_rank_accept_ratio,
        null_space_defect: evidence.null_space_defect,
        reduced_symmetry_defect: evidence.reduced_symmetry_defect,
        reduced_symmetry_defect_limit: evidence.symmetry_defect_limit,
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
        tolerance_round_trip_error: Some(solution.tolerance_round_trip_error),
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
        tolerance_round_trip_error: evidence.tolerance_round_trip_error,
        no_model_produced: evidence.no_model_produced,
    })
}

fn diagnose(failure: &CubicEqualityFailure) -> ProblemDiagnosis {
    match failure {
        CubicEqualityFailure::EmptyEqualitySet | CubicEqualityFailure::NonFiniteTarget { .. } => {
            ProblemDiagnosis::InvalidProblem
        }
        CubicEqualityFailure::Representation(failure) => diagnose_representation(failure),
        CubicEqualityFailure::Backend { failure, .. } => diagnose_kkt(failure),
        CubicEqualityFailure::RecoveryVerification { .. } => {
            ProblemDiagnosis::RecoveryVerificationFailure
        }
    }
}

fn diagnose_representation(failure: &RepresentationFailure) -> ProblemDiagnosis {
    match failure {
        RepresentationFailure::Capacity(_) => ProblemDiagnosis::CapacityExceeded,
        RepresentationFailure::PolynomialRankDeficient { .. } => {
            ProblemDiagnosis::UnidentifiedFieldMode
        }
        RepresentationFailure::PolynomialRankGrayZone { .. }
        | RepresentationFailure::ReducedPairingGrayZone { .. } => {
            ProblemDiagnosis::NumericalDecisionGrayZone
        }
        RepresentationFailure::ReducedPairingNotPositive { classification, .. } => {
            match classification {
                ReducedPairingFailureClassification::RankDeficient => {
                    ProblemDiagnosis::UnidentifiedFieldMode
                }
                ReducedPairingFailureClassification::NegativeCurvature => {
                    ProblemDiagnosis::NumericalFailure
                }
            }
        }
        RepresentationFailure::AffineReproductionBackend(failure) => diagnose_kkt(failure),
        _ => ProblemDiagnosis::NumericalFailure,
    }
}

fn diagnose_kkt(failure: &KktFailure) -> ProblemDiagnosis {
    match failure {
        KktFailure::Capacity(_) => ProblemDiagnosis::CapacityExceeded,
        KktFailure::RankDeficient { .. } => ProblemDiagnosis::UnidentifiedFieldMode,
        KktFailure::NumericalDecisionGrayZone { .. } => ProblemDiagnosis::NumericalDecisionGrayZone,
        KktFailure::BackendContractViolation { .. } => ProblemDiagnosis::BackendContractViolation,
        _ => ProblemDiagnosis::NumericalFailure,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic_equality::RecoveryVerificationFailureReason;

    #[test]
    fn snapshot_capacity_preflight_covers_all_source_relations_before_lowering() {
        assert!(plan_snapshot_capacity(2_000, 8).is_ok());
        let evidence = plan_snapshot_capacity(10_000, 0)
            .expect_err("all caller relations must count even if presolve could remove rows");
        assert!(!evidence.large_allocation_attempted);
        assert!(!evidence.backend_invocation_attempted);
        assert!(ProblemSize::cubic_equality(0, 0, 0, usize::MAX, 1, 0).is_none());
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
            relation_tolerances: None,
            hard_equality_violations: None,
            polynomial_round_trip_error: None,
            field_coefficient_round_trip_error: None,
            field_energy_round_trip_error: None,
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
}
