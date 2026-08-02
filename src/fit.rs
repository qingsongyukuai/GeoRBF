//! Strict fit outcomes and physical-unit fit reports.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use crate::capacity::{CapacityExceededEvidence, CapacityExceededReason};
use crate::cubic_equality::{
    AlgebraicAnalysisStage as InternalCubicAnalysisStage, CanonicalRelationToleranceEvidence,
    CpdEvidence, CubicEqualityCore, CubicEqualityFailure, CubicEqualitySolution, HardEquality,
    PhysicalSideConditionEvidence, RecoveryVerificationFailureEvidence,
    ReducedPairingFailureClassification, RepresentationFailure,
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
    ResidualDimension, ScalingFailureReason, ScalingSummary, SideConditionEvidence,
    SolveAttemptKind, SolveAttemptRecord, SolveAttemptRecordParts, SolveAttemptTermination,
    SolveCoordinateFailureReason, UnidentifiedAdditiveGaugeEvidence,
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
use crate::relation::{AdditiveFieldGaugeReference, SharedLevelSetInput};

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
    canonical_hard_equalities: usize,
    center_coefficients: usize,
    semantic_latents: usize,
    auxiliary_variables: usize,
    cone_blocks: usize,
    primal_variables: usize,
    equality_constraints: usize,
    kkt_dimension: usize,
}

impl ProblemSize {
    fn cubic_equality(
        input_observations: usize,
        scalar_hard_relations: usize,
        canonical_hard_equalities: usize,
        semantic_latents: usize,
    ) -> Self {
        let center_coefficients = canonical_hard_equalities;
        let primal_variables = center_coefficients + 4;
        let equality_constraints = canonical_hard_equalities + 4;
        Self {
            input_observations,
            scalar_hard_relations,
            canonical_hard_equalities,
            center_coefficients,
            semantic_latents,
            auxiliary_variables: 0,
            cone_blocks: 0,
            primal_variables,
            equality_constraints,
            kkt_dimension: primal_variables + equality_constraints,
        }
    }

    /// Returns the independent top-level input count for audit context.
    pub fn input_observations(self) -> usize {
        self.input_observations
    }

    /// Returns scalar hard components across every caller-owned source.
    pub fn scalar_hard_relations(self) -> usize {
        self.scalar_hard_relations
    }

    /// Returns the exact hard equalities retained after duplicate merging.
    pub fn canonical_hard_equalities(self) -> usize {
        self.canonical_hard_equalities
    }

    /// Returns the representer-center coefficient count.
    pub fn center_coefficients(self) -> usize {
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

    /// Returns the backend-standard-form primal dimension.
    pub fn primal_variables(self) -> usize {
        self.primal_variables
    }

    /// Returns the backend-standard-form equality dimension.
    pub fn equality_constraints(self) -> usize {
        self.equality_constraints
    }

    /// Returns the actual symmetric augmented KKT dimension.
    pub fn kkt_dimension(self) -> usize {
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
    scaled_kkt_tolerance: f64,
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

    /// Returns the tolerance after backend KKT scaling.
    pub fn scaled_kkt_tolerance(&self) -> f64 {
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
    let scalar_inputs = scalar_inputs(snapshot);
    let (level_gauges, level_gauge_conflict) = level_gauge_targets(snapshot);
    let lowering = lower_snapshot(
        &scalar_inputs,
        &snapshot.inner.shared_level_sets,
        &level_gauges,
        &snapshot.inner.field_unit,
    );
    let scalar_relation_count = scalar_observations(&snapshot.inner.observations).len()
        + snapshot
            .inner
            .shared_level_sets
            .iter()
            .map(|group| group.members().len())
            .sum::<usize>()
        + snapshot.inner.additive_field_gauges.len();
    let problem_size = ProblemSize::cubic_equality(
        snapshot.inner.source_count,
        scalar_relation_count,
        lowering.canonical_equalities.len(),
        lowering.semantic_latents.len(),
    );
    let base_report = || FitReport {
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
        execution_failure: None,
        cubic_analysis: None,
        backend_rank: None,
        inertia: None,
        canonical_acceptance: None,
        capacity: None,
        analysis_failure: None,
        unidentified_additive_gauge: None,
        uninformative_shared_level_set: None,
    };
    if let Some(conflict) = level_gauge_conflict.or_else(|| lowering.direct_input_conflict.clone())
    {
        let mut report = base_report();
        report.direct_input_conflict = Some(conflict);
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::DirectInputConflict,
            report: Box::new(report),
        });
    }
    let referenced_groups = level_gauges
        .keys()
        .collect::<std::collections::BTreeSet<_>>();
    if let Some(group) = snapshot.inner.shared_level_sets.iter().find(|group| {
        let distinct_locations = group
            .members()
            .iter()
            .map(|member| {
                member.location().components().map(|coordinate| {
                    if coordinate == 0.0 {
                        0.0_f64.to_bits()
                    } else {
                        coordinate.to_bits()
                    }
                })
            })
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        distinct_locations < 2 && !referenced_groups.contains(group.group_id())
    }) {
        let mut report = base_report();
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
        let mut report = base_report();
        report.unidentified_additive_gauge = Some(UnidentifiedAdditiveGaugeEvidence::new(
            source_ids, group_ids, false,
        ));
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::UnidentifiedAdditiveGauge,
            report: Box::new(report),
        });
    }
    let solution = match CubicEqualityCore::solve(
        lowering.canonical_equalities.clone(),
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
    let shared_level_values = recover_shared_level_values(&lowering.semantic_latents, &solution);
    let report = success_report(
        snapshot,
        problem_size,
        &solution,
        &lowering.public_source_relations,
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
    equality: HardEquality,
    canonical_index: usize,
}

#[derive(Debug, Clone)]
struct EqualityLowering {
    source_relations: Vec<SourceHardRelation>,
    public_source_relations: Vec<SourceHardRelation>,
    canonical_equalities: Vec<HardEquality>,
    canonical_index_by_key: BTreeMap<FunctionalKey, (usize, f64)>,
    direct_input_conflict: Option<DirectInputConflictEvidence>,
    semantic_latents: Vec<CanonicalSharedLevelLatent>,
}

#[derive(Debug, Clone)]
struct CanonicalSharedLevelLatent {
    group_id: GroupId,
    field_unit: FieldUnitLabel,
    members: Vec<CanonicalSharedLevelMember>,
}

#[derive(Debug, Clone)]
struct CanonicalSharedLevelMember {
    source_id: SourceId,
    support: [f64; 3],
}

fn lower_observations(observations: &[ScalarObservation]) -> EqualityLowering {
    let mut lowering = EqualityLowering {
        public_source_relations: Vec::new(),
        source_relations: Vec::new(),
        canonical_equalities: Vec::new(),
        canonical_index_by_key: BTreeMap::new(),
        direct_input_conflict: None,
        semantic_latents: Vec::new(),
    };
    for observation in observations {
        lowering.push(hard_equality(observation), true);
    }
    lowering
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct FunctionalKey {
    dimension: u8,
    terms: Vec<([u64; 3], u64, [u64; 3])>,
}

impl EqualityLowering {
    fn push(&mut self, equality: HardEquality, public: bool) {
        let (key, normalized_target) = normalized_equality_key(&equality);
        let canonical_index =
            if let Some((index, first_target)) = self.canonical_index_by_key.get(&key).copied() {
                if first_target != normalized_target && self.direct_input_conflict.is_none() {
                    let first = &self.canonical_equalities[index];
                    self.direct_input_conflict = Some(DirectInputConflictEvidence::new(
                        first.usage().provenance().source().clone(),
                        equality.usage().provenance().source().clone(),
                        equality.usage().provenance().semantic_role().clone(),
                        first.target(),
                        equality.target(),
                    ));
                }
                index
            } else {
                let index = self.canonical_equalities.len();
                self.canonical_equalities.push(equality.clone());
                self.canonical_index_by_key
                    .insert(key, (index, normalized_target));
                index
            };
        let relation = SourceHardRelation {
            equality,
            canonical_index,
        };
        self.source_relations.push(relation.clone());
        if public {
            self.public_source_relations.push(relation);
        }
    }
}

fn normalized_equality_key(equality: &HardEquality) -> (FunctionalKey, f64) {
    let functional = equality.usage().functional();
    let first_coefficient = functional
        .terms()
        .iter()
        .flat_map(|term| {
            std::iter::once(term.value_coefficient()).chain(term.gradient_coefficient())
        })
        .find(|coefficient| *coefficient != 0.0)
        .expect("canonical functionals are nonzero");
    let sign = if first_coefficient.is_sign_negative() {
        -1.0
    } else {
        1.0
    };
    let dimension = match functional.dimension() {
        FunctionalDimension::FieldValue => 0,
        FunctionalDimension::FieldValuePerLength => 1,
    };
    let terms = functional
        .terms()
        .iter()
        .map(|term| {
            (
                term.support().map(f64::to_bits),
                signed_coefficient_bits(sign, term.value_coefficient()),
                term.gradient_coefficient()
                    .map(|coefficient| signed_coefficient_bits(sign, coefficient)),
            )
        })
        .collect();
    (FunctionalKey { dimension, terms }, sign * equality.target())
}

fn signed_coefficient_bits(sign: f64, coefficient: f64) -> u64 {
    let value = sign * coefficient;
    if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    }
}

#[derive(Debug, Clone)]
struct LevelGaugeTarget {
    source_id: SourceId,
    value: f64,
}

fn level_gauge_targets(
    snapshot: &ProblemSnapshot,
) -> (
    BTreeMap<GroupId, LevelGaugeTarget>,
    Option<DirectInputConflictEvidence>,
) {
    let mut targets = BTreeMap::<GroupId, LevelGaugeTarget>::new();
    for gauge in &snapshot.inner.additive_field_gauges {
        let AdditiveFieldGaugeReference::LevelSet(group_id) = gauge.reference() else {
            continue;
        };
        if let Some(first) = targets.get(group_id).cloned() {
            if first.value != gauge.value() {
                return (
                    targets,
                    Some(DirectInputConflictEvidence::new(
                        first.source_id.clone(),
                        gauge.source_id().clone(),
                        SemanticRolePath::new("additive-field-gauge/level-set"),
                        first.value,
                        gauge.value(),
                    )),
                );
            }
        } else {
            targets.insert(
                group_id.clone(),
                LevelGaugeTarget {
                    source_id: gauge.source_id().clone(),
                    value: gauge.value(),
                },
            );
        }
    }
    (targets, None)
}

fn scalar_inputs(snapshot: &ProblemSnapshot) -> Vec<ScalarObservation> {
    let mut inputs = scalar_observations(&snapshot.inner.observations);
    inputs.extend(
        snapshot
            .inner
            .additive_field_gauges
            .iter()
            .filter_map(|gauge| match gauge.reference() {
                AdditiveFieldGaugeReference::Point(point) => Some(ScalarObservation {
                    source_id: gauge.source_id().clone(),
                    group_id: None,
                    support: point.components(),
                    component: DirectInputComponent::FieldValue,
                    semantic_role: SemanticRolePath::new("additive-field-gauge/point"),
                    target: gauge.value(),
                }),
                AdditiveFieldGaugeReference::LevelSet(_) => None,
            }),
    );
    inputs.sort_by(|left, right| {
        left.source_id
            .cmp(&right.source_id)
            .then_with(|| left.semantic_role.cmp(&right.semantic_role))
    });
    inputs
}

fn lower_snapshot(
    scalar_inputs: &[ScalarObservation],
    groups: &[SharedLevelSetInput],
    level_gauges: &BTreeMap<GroupId, LevelGaugeTarget>,
    field_unit: &FieldUnitLabel,
) -> EqualityLowering {
    let mut lowering = lower_observations(scalar_inputs);
    for group in groups {
        lowering.semantic_latents.push(CanonicalSharedLevelLatent {
            group_id: group.group_id().clone(),
            field_unit: field_unit.clone(),
            members: group
                .members()
                .iter()
                .map(|member| CanonicalSharedLevelMember {
                    source_id: member.source_id().clone(),
                    support: member.location().components(),
                })
                .collect(),
        });
        if let Some(gauge) = level_gauges.get(group.group_id()) {
            for member in group.members() {
                let equality = hard_equality(&ScalarObservation {
                    source_id: member.source_id().clone(),
                    group_id: Some(group.group_id().clone()),
                    support: member.location().components(),
                    component: DirectInputComponent::FieldValue,
                    semantic_role: SemanticRolePath::new("shared-level-set/member/value"),
                    target: gauge.value,
                });
                lowering.push(equality, false);
            }
            continue;
        }

        for (index, member) in group.members()[1..].iter().enumerate() {
            let previous_count = index + 1;
            let mut terms = group.members()[..previous_count]
                .iter()
                .map(|previous| {
                    FunctionalTerm::new(
                        previous.location().components(),
                        -1.0 / previous_count as f64,
                        [0.0; 3],
                    )
                })
                .collect::<Vec<_>>();
            terms.push(FunctionalTerm::new(
                member.location().components(),
                1.0,
                [0.0; 3],
            ));
            let Ok(functional) = CanonicalFunctional::new(FunctionalDimension::FieldValue, terms)
            else {
                continue;
            };
            let semantic_role = SemanticRolePath::new("shared-level-set/member-equality");
            let relation_id = RelationId::new(format!(
                "{}:{}:{}",
                group.group_id().as_str(),
                member.source_id().as_str(),
                semantic_role.as_str()
            ));
            let residual_id = ResidualId::new(format!("{}/residual", relation_id.as_str()));
            let equality = HardEquality::new(
                FunctionalUse::new(
                    functional,
                    UsageProvenance::new(
                        member.source_id().clone(),
                        Some(group.group_id().clone()),
                        relation_id,
                        residual_id,
                        semantic_role,
                    ),
                ),
                0.0,
            );
            lowering.push(equality, false);
        }
    }
    lowering
}

fn hard_equality(observation: &ScalarObservation) -> HardEquality {
    let semantic_role = observation.semantic_role.clone();
    let functional = CanonicalFunctional::new(
        observation.dimension(),
        vec![FunctionalTerm::new(
            observation.support,
            observation.value_coefficient(),
            observation.gradient_coefficient(),
        )],
    )
    .expect("checked public observations lower to a finite nonzero functional");
    let relation_id = RelationId::new(format!(
        "{}:{}",
        observation.source_id.as_str(),
        semantic_role.as_str()
    ));
    let residual_id = ResidualId::new(format!("{}/residual", relation_id.as_str()));
    HardEquality::new(
        FunctionalUse::new(
            functional,
            UsageProvenance::new(
                observation.source_id.clone(),
                observation.group_id.clone(),
                relation_id,
                residual_id,
                semantic_role,
            ),
        ),
        observation.target,
    )
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
        canonical_hard_equalities: solution.assembly.hard_equalities,
        center_coefficients: solution.assembly.field_coefficients,
        semantic_latents: planned_problem_size.semantic_latents,
        auxiliary_variables: 0,
        cone_blocks: 0,
        primal_variables: solution.assembly.primal_variables,
        equality_constraints: solution.assembly.side_conditions + solution.assembly.hard_equalities,
        kkt_dimension: solution.backend.capacity.kkt_dimension,
    };
    let mut hard_relations = source_relations
        .iter()
        .map(|source_relation| {
            let usage = source_relation.equality.usage();
            let tolerance = solution.relation_tolerances[source_relation.canonical_index];
            let recovered_value = solution.field.evaluate_functional(usage.functional());
            hard_relation_assessment(source_relation, recovered_value, tolerance)
        })
        .collect::<Vec<_>>();
    if !shared_level_values.is_empty() {
        let tolerance = solution
            .relation_tolerances
            .iter()
            .copied()
            .filter(|evidence| evidence.dimension == FunctionalDimension::FieldValue)
            .max_by(|left, right| left.physical_tolerance.total_cmp(&right.physical_tolerance))
            .expect("a shared level set contributes a field-value equality");
        for (group, recovered) in snapshot
            .inner
            .shared_level_sets
            .iter()
            .zip(shared_level_values)
        {
            for member in group.members() {
                let field_value = solution.field.sample(member.location().components()).value;
                hard_relations.push(hard_relation_assessment_from_parts(
                    member.source_id().clone(),
                    Some(group.group_id().clone()),
                    SemanticRolePath::new("shared-level-set/member/value"),
                    recovered.value,
                    field_value,
                    tolerance,
                ));
            }
            hard_relations.extend(snapshot.inner.additive_field_gauges.iter().filter_map(
                |gauge| match gauge.reference() {
                    AdditiveFieldGaugeReference::LevelSet(group_id)
                        if group_id == group.group_id() =>
                    {
                        Some(hard_relation_assessment_from_parts(
                            gauge.source_id().clone(),
                            Some(group.group_id().clone()),
                            SemanticRolePath::new("additive-field-gauge/level-set"),
                            gauge.value(),
                            recovered.value,
                            tolerance,
                        ))
                    }
                    _ => None,
                },
            ));
        }
    }
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
    let usage = source_relation.equality.usage();
    let target = source_relation.equality.target();
    HardRelationAssessment {
        source_id: usage.provenance().source().clone(),
        group_id: usage.provenance().group().cloned(),
        semantic_role: usage.provenance().semantic_role().clone(),
        dimension: match usage.functional().dimension() {
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

fn hard_relation_assessment_from_parts(
    source_id: SourceId,
    group_id: Option<GroupId>,
    semantic_role: SemanticRolePath,
    target: f64,
    recovered_value: f64,
    tolerance: CanonicalRelationToleranceEvidence,
) -> HardRelationAssessment {
    HardRelationAssessment {
        source_id,
        group_id,
        semantic_role,
        dimension: ResidualDimension::FieldValue,
        target,
        recovered_value,
        residual: recovered_value - target,
        tolerance: tolerance.physical_tolerance,
        characteristic_scale: tolerance.characteristic_scale,
        relation_reference_scale: target.abs(),
        standard_tolerance: tolerance.standard_tolerance,
        scaled_kkt_tolerance: tolerance.scaled_kkt_tolerance,
        recovered_physical_tolerance: tolerance.recovered_physical_tolerance,
        tolerance_round_trip_error: tolerance.round_trip_error,
    }
}

fn recover_shared_level_values(
    latents: &[CanonicalSharedLevelLatent],
    solution: &CubicEqualitySolution,
) -> Vec<SharedLevelValue> {
    latents
        .iter()
        .map(|latent| {
            let values = latent
                .members
                .iter()
                .map(|member| solution.field.sample(member.support).value)
                .collect::<Vec<_>>();
            let origin = values[0];
            let value = origin
                + values.iter().map(|value| *value - origin).sum::<f64>() / values.len() as f64;
            SharedLevelValue {
                group_id: latent.group_id.clone(),
                value,
                field_unit: latent.field_unit.clone(),
                member_source_ids: latent
                    .members
                    .iter()
                    .map(|member| member.source_id.clone())
                    .collect::<Vec<_>>()
                    .into(),
            }
        })
        .collect()
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
                report.problem_size.kkt_dimension,
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
                report.problem_size.kkt_dimension,
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
