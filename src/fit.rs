//! Strict fit outcomes and physical-unit fit reports.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use crate::cubic_equality::{
    CubicEqualityCore, CubicEqualityFailure, CubicEqualitySolution, HardEquality,
    RecoveryVerificationFailureEvidence, ReducedPairingFailureClassification,
    RepresentationFailure,
};
use crate::diagnostics::{
    AttemptFailureCategory, AttemptFailureEvidence, BackendAttemptSettings, BackendFingerprint,
    BackendFingerprintParts, DirectInputConflictEvidence, LinearResidualEvidence, ProblemDiagnosis,
    RecoveryVerificationEvidence, ResidualDimension, ScalingSummary, SolveAttemptKind,
    SolveAttemptRecord, SolveAttemptRecordParts, SolveAttemptTermination,
};
use crate::functional::{
    CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, RelationId,
    ResidualId, SemanticRolePath, SourceId, UsageProvenance,
};
use crate::kernel::{FieldEnergyNormalization, KernelConfig};
use crate::kkt::{
    BackendAttemptSettings as InternalAttemptSettings,
    BackendContractViolationReason as InternalBackendContractReason,
    BackendFingerprint as InternalBackendFingerprint, KktAttemptFailureReason, KktAttemptKind,
    KktAttemptRecord, KktFailure, NumericalFailureReason,
    SolveAttemptTermination as InternalTermination,
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
    center_coefficients: usize,
    semantic_latents: usize,
    auxiliary_variables: usize,
    cone_blocks: usize,
    primal_variables: usize,
    equality_constraints: usize,
    kkt_dimension: usize,
}

impl ProblemSize {
    fn cubic_equality(input_observations: usize, scalar_hard_relations: usize) -> Self {
        let center_coefficients = scalar_hard_relations;
        let primal_variables = center_coefficients + 4;
        let equality_constraints = scalar_hard_relations + 4;
        Self {
            input_observations,
            scalar_hard_relations,
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

    /// Returns the number of scalar hard equality components after lowering.
    pub fn scalar_hard_relations(self) -> usize {
        self.scalar_hard_relations
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
}

pub(crate) fn fit_snapshot(snapshot: &ProblemSnapshot) -> Result<FitSuccess, FitFailure> {
    let equalities = lower_observations(&snapshot.inner.observations);
    let problem_size =
        ProblemSize::cubic_equality(snapshot.inner.observations.len(), equalities.len());
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
    };
    if let Some(conflict) = direct_input_conflict(&snapshot.inner.observations) {
        let mut report = base_report();
        report.direct_input_conflict = Some(conflict);
        return Err(FitFailure {
            diagnosis: ProblemDiagnosis::DirectInputConflict,
            report: Box::new(report),
        });
    }
    let solution = match CubicEqualityCore::solve(
        equalities,
        snapshot.inner.global_anisotropy_metric.as_cubic_metric(),
    ) {
        Ok(solution) => solution,
        Err(failure) => {
            return Err(FitFailure {
                diagnosis: diagnose(&failure),
                report: Box::new(failure_report(base_report(), &failure)),
            });
        }
    };
    let report = success_report(snapshot, problem_size, &solution);
    let model = SolvedModel::new(snapshot.clone(), solution.field);
    Ok(FitSuccess { model, report })
}

fn lower_observations(observations: &[ObservationInput]) -> Vec<HardEquality> {
    observations
        .iter()
        .flat_map(|observation| match observation {
            ObservationInput::FieldValue(observation) => vec![hard_equality(
                observation.source_id(),
                observation.location().components(),
                FunctionalDimension::FieldValue,
                1.0,
                [0.0; 3],
                observation.value(),
                "field-value-observation/value".into(),
            )],
            ObservationInput::Gradient(observation) => {
                let target = observation.gradient().components();
                (0..3)
                    .map(|axis| {
                        hard_equality(
                            observation.source_id(),
                            observation.location().components(),
                            FunctionalDimension::FieldValuePerLength,
                            0.0,
                            std::array::from_fn(
                                |component| {
                                    if component == axis { 1.0 } else { 0.0 }
                                },
                            ),
                            target[axis],
                            format!("gradient-observation/component/{axis}"),
                        )
                    })
                    .collect()
            }
        })
        .collect()
}

fn hard_equality(
    source_id: &SourceId,
    support: [f64; 3],
    dimension: FunctionalDimension,
    value_coefficient: f64,
    gradient_coefficient: [f64; 3],
    target: f64,
    semantic_role: String,
) -> HardEquality {
    let functional = CanonicalFunctional::new(
        dimension,
        vec![FunctionalTerm::new(
            support,
            value_coefficient,
            gradient_coefficient,
        )],
    )
    .expect("checked public observations lower to a finite nonzero functional");
    let relation_id = RelationId::new(format!("{}:{semantic_role}", source_id.as_str()));
    let residual_id = ResidualId::new(format!("{}/residual", relation_id.as_str()));
    HardEquality::new(
        FunctionalUse::new(
            functional,
            UsageProvenance::new(
                source_id.clone(),
                None,
                relation_id,
                residual_id,
                SemanticRolePath::new(semantic_role),
            ),
        ),
        target,
    )
}

fn success_report(
    snapshot: &ProblemSnapshot,
    planned_problem_size: ProblemSize,
    solution: &CubicEqualitySolution,
) -> FitReport {
    let problem_size = ProblemSize {
        input_observations: planned_problem_size.input_observations,
        scalar_hard_relations: solution.assembly.hard_equalities,
        center_coefficients: solution.assembly.field_coefficients,
        semantic_latents: solution.semantic_latent_count,
        auxiliary_variables: 0,
        cone_blocks: 0,
        primal_variables: solution.assembly.primal_variables,
        equality_constraints: solution.assembly.side_conditions + solution.assembly.hard_equalities,
        kkt_dimension: solution.backend.capacity.kkt_dimension,
    };
    let hard_relations = solution
        .hard_equalities
        .iter()
        .zip(&solution.relation_tolerances)
        .map(|(relation, tolerance)| HardRelationAssessment {
            source_id: relation.usage.provenance().source().clone(),
            semantic_role: relation.usage.provenance().semantic_role().clone(),
            dimension: match relation.usage.functional().dimension() {
                FunctionalDimension::FieldValue => ResidualDimension::FieldValue,
                FunctionalDimension::FieldValuePerLength => ResidualDimension::FieldValuePerLength,
            },
            target: relation.target,
            recovered_value: relation.value,
            residual: relation.residual,
            tolerance: tolerance.physical_tolerance,
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
    }
}

fn direct_input_conflict(observations: &[ObservationInput]) -> Option<DirectInputConflictEvidence> {
    let mut first_by_key = BTreeMap::new();
    for observation in observations {
        match observation {
            ObservationInput::FieldValue(observation) => {
                if let Some(conflict) = register_direct_input(
                    &mut first_by_key,
                    observation.location().components(),
                    DirectInputComponent::FieldValue,
                    observation.source_id(),
                    observation.value(),
                ) {
                    return Some(conflict);
                }
            }
            ObservationInput::Gradient(observation) => {
                for (axis, target) in observation.gradient().components().into_iter().enumerate() {
                    if let Some(conflict) = register_direct_input(
                        &mut first_by_key,
                        observation.location().components(),
                        DirectInputComponent::Gradient(axis),
                        observation.source_id(),
                        target,
                    ) {
                        return Some(conflict);
                    }
                }
            }
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
    location: [f64; 3],
    component: DirectInputComponent,
    source: &SourceId,
    target: f64,
) -> Option<DirectInputConflictEvidence> {
    let location = location.map(|coordinate| {
        if coordinate == 0.0 {
            0.0_f64.to_bits()
        } else {
            coordinate.to_bits()
        }
    });
    match first_by_key.entry((location, component)) {
        std::collections::btree_map::Entry::Vacant(entry) => {
            entry.insert((source.clone(), target));
            None
        }
        std::collections::btree_map::Entry::Occupied(entry) => {
            let (first_source, first_target) = entry.get();
            (*first_target != target).then(|| {
                DirectInputConflictEvidence::new(
                    first_source.clone(),
                    source.clone(),
                    component.semantic_role(),
                    *first_target,
                    target,
                )
            })
        }
    }
}

fn failure_report(mut report: FitReport, failure: &CubicEqualityFailure) -> FitReport {
    match failure {
        CubicEqualityFailure::Backend(failure) => {
            report.attempts = public_attempts(kkt_failure_attempts(failure));
            report.execution_failure = public_kkt_failure(failure);
        }
        CubicEqualityFailure::Representation(failure) => {
            if let RepresentationFailure::AffineReproductionBackend(failure) = failure.as_ref() {
                report.attempts = public_attempts(kkt_failure_attempts(failure));
                report.execution_failure = public_kkt_failure(failure);
            }
        }
        CubicEqualityFailure::RecoveryVerification { evidence, backend } => {
            report.backend_fingerprint = Some(public_backend_fingerprint(&backend.backend));
            report.attempts = public_attempts(&backend.attempts);
            report.recovery_verification = Some(public_recovery_evidence(evidence));
        }
        CubicEqualityFailure::EmptyEqualitySet | CubicEqualityFailure::NonFiniteTarget { .. } => {}
    }
    report
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

fn public_recovery_evidence(
    evidence: &RecoveryVerificationFailureEvidence,
) -> RecoveryVerificationEvidence {
    RecoveryVerificationEvidence::new(
        evidence.reasons.clone(),
        evidence
            .hard_equality_violations
            .map(|envelope| (envelope.field_value, envelope.field_value_per_length)),
        evidence.polynomial_round_trip_error,
        evidence.field_coefficient_round_trip_error,
        evidence.field_energy_round_trip_error,
        evidence.tolerance_round_trip_error,
        evidence.no_model_produced,
    )
}

fn diagnose(failure: &CubicEqualityFailure) -> ProblemDiagnosis {
    match failure {
        CubicEqualityFailure::EmptyEqualitySet | CubicEqualityFailure::NonFiniteTarget { .. } => {
            ProblemDiagnosis::InvalidProblem
        }
        CubicEqualityFailure::Representation(failure) => diagnose_representation(failure),
        CubicEqualityFailure::Backend(failure) => diagnose_kkt(failure),
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
