//! Strict fit outcomes and physical-unit fit reports.

use std::error::Error;
use std::fmt;

use crate::cubic_equality::{
    CubicEqualityCore, CubicEqualityFailure, CubicEqualitySolution, HardEquality,
    RecoveryVerificationFailureEvidence, RecoveryVerificationFailureReason,
    ReducedPairingFailureClassification, RepresentationFailure,
};
use crate::diagnostics::{
    BackendFingerprint, ProblemDiagnosis, RecoveryVerificationEvidence, RecoveryVerificationReason,
    ResidualDimension, SolveAttemptRecord, SolveAttemptTermination,
};
use crate::functional::{
    CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, RelationId,
    ResidualId, SemanticRolePath, SourceId, UsageProvenance,
};
use crate::kernel::{FieldEnergyNormalization, KernelConfig};
use crate::kkt::{
    BackendFingerprint as InternalBackendFingerprint, KktAttemptRecord, KktFailure,
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
    observations: usize,
    scalar_hard_relations: usize,
}

impl ProblemSize {
    /// Returns the number of top-level domain observations.
    pub fn observations(self) -> usize {
        self.observations
    }

    /// Returns the number of scalar hard equality components after lowering.
    pub fn scalar_hard_relations(self) -> usize {
        self.scalar_hard_relations
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
}

pub(crate) fn fit_snapshot(snapshot: &ProblemSnapshot) -> Result<FitSuccess, FitFailure> {
    let equalities = lower_observations(&snapshot.inner.observations);
    let problem_size = ProblemSize {
        observations: snapshot.inner.observations.len(),
        scalar_hard_relations: equalities.len(),
    };
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
    };
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
    problem_size: ProblemSize,
    solution: &CubicEqualitySolution,
) -> FitReport {
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
    }
}

fn failure_report(mut report: FitReport, failure: &CubicEqualityFailure) -> FitReport {
    match failure {
        CubicEqualityFailure::Backend(failure) => {
            report.attempts = public_attempts(kkt_failure_attempts(failure));
        }
        CubicEqualityFailure::Representation(failure) => {
            if let RepresentationFailure::AffineReproductionBackend(failure) = failure.as_ref() {
                report.attempts = public_attempts(kkt_failure_attempts(failure));
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

fn public_attempts(attempts: &[KktAttemptRecord]) -> Vec<SolveAttemptRecord> {
    attempts
        .iter()
        .map(|attempt| {
            SolveAttemptRecord::new(
                attempt.sequence,
                match attempt.termination {
                    InternalTermination::AcceptedCandidate => {
                        SolveAttemptTermination::AcceptedCandidate
                    }
                    InternalTermination::RejectedCandidate => {
                        SolveAttemptTermination::RejectedCandidate
                    }
                    InternalTermination::NumericalError => SolveAttemptTermination::NumericalError,
                },
                attempt
                    .residual
                    .map(|residual| residual.normalized_backward_error),
                public_backend_fingerprint(&attempt.backend),
            )
        })
        .collect()
}

fn public_backend_fingerprint(backend: &InternalBackendFingerprint) -> BackendFingerprint {
    BackendFingerprint::new(
        backend.crate_name,
        backend.crate_version,
        backend.algorithm,
        backend.requested_threads,
        backend.actual_threads,
    )
}

fn public_recovery_evidence(
    evidence: &RecoveryVerificationFailureEvidence,
) -> RecoveryVerificationEvidence {
    RecoveryVerificationEvidence::new(
        evidence
            .reasons
            .iter()
            .copied()
            .map(|reason| match reason {
                RecoveryVerificationFailureReason::InvalidRecoveryMap => {
                    RecoveryVerificationReason::InvalidRecoveryMap
                }
                RecoveryVerificationFailureReason::ProvenanceMismatch => {
                    RecoveryVerificationReason::ProvenanceMismatch
                }
                RecoveryVerificationFailureReason::NonFiniteRecoveredQuantity => {
                    RecoveryVerificationReason::NonFiniteRecoveredQuantity
                }
                RecoveryVerificationFailureReason::SideConditionViolation => {
                    RecoveryVerificationReason::SideConditionViolation
                }
                RecoveryVerificationFailureReason::SideConditionRoundTripViolation => {
                    RecoveryVerificationReason::SideConditionRoundTripViolation
                }
                RecoveryVerificationFailureReason::HardEqualityViolation => {
                    RecoveryVerificationReason::HardEqualityViolation
                }
                RecoveryVerificationFailureReason::PolynomialRoundTripViolation => {
                    RecoveryVerificationReason::PolynomialRoundTripViolation
                }
                RecoveryVerificationFailureReason::FieldCoefficientRoundTripViolation => {
                    RecoveryVerificationReason::FieldCoefficientRoundTripViolation
                }
                RecoveryVerificationFailureReason::FieldEnergyRoundTripViolation => {
                    RecoveryVerificationReason::FieldEnergyRoundTripViolation
                }
                RecoveryVerificationFailureReason::ToleranceRoundTripViolation => {
                    RecoveryVerificationReason::ToleranceRoundTripViolation
                }
            })
            .collect(),
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
