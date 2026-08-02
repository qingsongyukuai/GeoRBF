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
    SolveCoordinateFailureReason,
};
use crate::functional::{
    CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, RelationId,
    ResidualId, SemanticRolePath, SourceId, UsageProvenance,
};
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
    ) -> Self {
        let center_coefficients = canonical_hard_equalities;
        let primal_variables = center_coefficients + 4;
        let equality_constraints = canonical_hard_equalities + 4;
        Self {
            input_observations,
            scalar_hard_relations,
            canonical_hard_equalities,
            center_coefficients,
            semantic_latents: 0,
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

/// Typed audit report shared by successful and failed fits.
#[derive(Debug, Clone)]
pub struct FitReport {
    problem_size: ProblemSize,
    resolved_kernel: KernelConfig,
    field_energy_normalization: FieldEnergyNormalization,
    numerical_policy: NumericalPolicyId,
    requested_thread_budget: ThreadBudget,
    hard_relations: Vec<HardRelationAssessment>,
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
}

pub(crate) fn fit_snapshot(snapshot: &ProblemSnapshot) -> Result<FitSuccess, FitFailure> {
    let scalar_observations = scalar_observations(&snapshot.inner.observations);
    let lowering = lower_observations(&scalar_observations);
    let problem_size = ProblemSize::cubic_equality(
        snapshot.inner.observations.len(),
        lowering.source_relations.len(),
        lowering.canonical_equalities.len(),
    );
    let base_report = || FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: snapshot.inner.fit_configuration.numerical_policy(),
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations: Vec::new(),
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
    };
    if let Some(conflict) = direct_input_conflict(&scalar_observations) {
        let mut report = base_report();
        report.direct_input_conflict = Some(conflict);
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::DirectInputConflict,
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
    let report = success_report(
        snapshot,
        problem_size,
        &solution,
        &lowering.source_relations,
    );
    let model = SolvedModel::new(snapshot.clone(), solution.field);
    Ok(FitSuccess { model, report })
}

#[derive(Debug, Clone)]
struct ScalarObservation {
    source_id: SourceId,
    support: [f64; 3],
    component: DirectInputComponent,
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
                support: observation.location().components(),
                component: DirectInputComponent::FieldValue,
                target: observation.value(),
            }],
            ObservationInput::Gradient(observation) => observation
                .gradient()
                .components()
                .into_iter()
                .enumerate()
                .map(|(axis, target)| ScalarObservation {
                    source_id: observation.source_id().clone(),
                    support: observation.location().components(),
                    component: DirectInputComponent::Gradient(axis),
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
    canonical_equalities: Vec<HardEquality>,
}

fn lower_observations(observations: &[ScalarObservation]) -> EqualityLowering {
    let mut canonical_index_by_key = BTreeMap::new();
    let mut canonical_equalities = Vec::new();
    let source_relations = observations
        .iter()
        .map(|observation| {
            let equality = hard_equality(observation);
            let next_index = canonical_equalities.len();
            let canonical_index = *canonical_index_by_key
                .entry(direct_input_key(observation.support, observation.component))
                .or_insert_with(|| {
                    canonical_equalities.push(equality.clone());
                    next_index
                });
            SourceHardRelation {
                equality,
                canonical_index,
            }
        })
        .collect();
    EqualityLowering {
        source_relations,
        canonical_equalities,
    }
}

fn hard_equality(observation: &ScalarObservation) -> HardEquality {
    let semantic_role = observation.component.semantic_role();
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
                None,
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
) -> FitReport {
    let problem_size = ProblemSize {
        input_observations: planned_problem_size.input_observations,
        scalar_hard_relations: planned_problem_size.scalar_hard_relations,
        canonical_hard_equalities: solution.assembly.hard_equalities,
        center_coefficients: solution.assembly.field_coefficients,
        semantic_latents: solution.semantic_latent_count,
        auxiliary_variables: 0,
        cone_blocks: 0,
        primal_variables: solution.assembly.primal_variables,
        equality_constraints: solution.assembly.side_conditions + solution.assembly.hard_equalities,
        kkt_dimension: solution.backend.capacity.kkt_dimension,
    };
    let hard_relations = source_relations
        .iter()
        .map(|source_relation| {
            let usage = source_relation.equality.usage();
            let tolerance = solution.relation_tolerances[source_relation.canonical_index];
            let recovered_value = solution.field.evaluate_functional(usage.functional());
            hard_relation_assessment(source_relation, recovered_value, tolerance)
        })
        .collect();
    let backend_fingerprint = public_backend_fingerprint(&solution.backend.backend);
    let attempts = public_attempts(&solution.backend.attempts);
    FitReport {
        problem_size,
        resolved_kernel: snapshot.inner.resolved_kernel.clone(),
        field_energy_normalization: snapshot.inner.field_energy_normalization,
        numerical_policy: solution.backend.numerical_policy,
        requested_thread_budget: snapshot.inner.fit_configuration.thread_budget(),
        hard_relations,
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

fn direct_input_conflict(
    observations: &[ScalarObservation],
) -> Option<DirectInputConflictEvidence> {
    let mut first_by_key = BTreeMap::new();
    for observation in observations {
        if let Some(conflict) = register_direct_input(&mut first_by_key, observation) {
            return Some(conflict);
        }
    }
    None
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DirectInputComponent {
    FieldValue,
    Gradient(usize),
}

impl DirectInputComponent {
    fn semantic_role(self) -> SemanticRolePath {
        match self {
            Self::FieldValue => SemanticRolePath::new("field-value-observation/value"),
            Self::Gradient(axis) => {
                SemanticRolePath::new(format!("gradient-observation/component/{axis}"))
            }
        }
    }
}

fn register_direct_input(
    first_by_key: &mut BTreeMap<([u64; 3], DirectInputComponent), (SourceId, f64)>,
    observation: &ScalarObservation,
) -> Option<DirectInputConflictEvidence> {
    let key = direct_input_key(observation.support, observation.component);
    match first_by_key.entry(key) {
        std::collections::btree_map::Entry::Vacant(entry) => {
            entry.insert((observation.source_id.clone(), observation.target));
            None
        }
        std::collections::btree_map::Entry::Occupied(entry) => {
            let (first_source, first_target) = entry.get();
            (*first_target != observation.target).then(|| {
                DirectInputConflictEvidence::new(
                    first_source.clone(),
                    observation.source_id.clone(),
                    observation.component.semantic_role(),
                    *first_target,
                    observation.target,
                )
            })
        }
    }
}

fn direct_input_key(
    location: [f64; 3],
    component: DirectInputComponent,
) -> ([u64; 3], DirectInputComponent) {
    let location = location.map(|coordinate| {
        if coordinate == 0.0 {
            0.0_f64.to_bits()
        } else {
            coordinate.to_bits()
        }
    });
    (location, component)
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
