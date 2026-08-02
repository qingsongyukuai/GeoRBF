//! Typed public fit diagnoses and backend-attempt evidence.

pub use crate::cubic_equality::RecoveryVerificationFailureReason as RecoveryVerificationReason;
pub use crate::functional::SemanticRolePath;
use crate::functional::SourceId;
pub use crate::numerical::NumericalPolicyId;

/// GeoRBF's semantic conclusion for a failed fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProblemDiagnosis {
    /// The supplied problem could not form a valid fitting problem.
    InvalidProblem,
    /// Two locally comparable hard inputs have incompatible exact targets.
    DirectInputConflict,
    /// The observations did not identify every observable field mode.
    UnidentifiedFieldMode,
    /// A numerical decision fell between the versioned accept/reject bands.
    NumericalDecisionGrayZone,
    /// The checked peak-memory plan exceeded the supported capacity.
    CapacityExceeded,
    /// A backend candidate violated its backend-standard-form contract.
    BackendContractViolation,
    /// Backend-standard form passed, but physical recovery verification failed.
    RecoveryVerificationFailure,
    /// Numerical execution failed without proving a stronger diagnosis.
    NumericalFailure,
}

/// Stable source and target evidence for a direct hard-input contradiction.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectInputConflictEvidence {
    first_source: SourceId,
    second_source: SourceId,
    semantic_role: SemanticRolePath,
    first_target: f64,
    second_target: f64,
}

impl DirectInputConflictEvidence {
    pub(crate) fn new(
        first_source: SourceId,
        second_source: SourceId,
        semantic_role: SemanticRolePath,
        first_target: f64,
        second_target: f64,
    ) -> Self {
        Self {
            first_source,
            second_source,
            semantic_role,
            first_target,
            second_target,
        }
    }

    /// Returns the first conflicting SourceId in stable order.
    pub fn first_source(&self) -> &SourceId {
        &self.first_source
    }

    /// Returns the second conflicting SourceId in stable order.
    pub fn second_source(&self) -> &SourceId {
        &self.second_source
    }

    /// Returns the conflicting scalar semantic component.
    pub fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }

    /// Returns the first exact hard target.
    pub fn first_target(&self) -> f64 {
        self.first_target
    }

    /// Returns the incompatible exact hard target.
    pub fn second_target(&self) -> f64 {
        self.second_target
    }
}

/// Physical dimension of one scalar hard-relation residual.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ResidualDimension {
    /// Scalar field-value units.
    FieldValue,
    /// Scalar field-value-per-length units.
    FieldValuePerLength,
}

/// Backend termination evidence, distinct from [`ProblemDiagnosis`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolveAttemptTermination {
    /// A candidate passed the backend-standard-form contract.
    AcceptedCandidate,
    /// A candidate was rejected by that contract.
    RejectedCandidate,
    /// The backend attempt stopped on a numerical error.
    NumericalError,
}

/// Auditable identity of the backend used by one fit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BackendFingerprint {
    schema_version: u32,
    crate_name: Box<str>,
    crate_version: Box<str>,
    features: [Box<str>; 2],
    algorithm: Box<str>,
    target_arch: Box<str>,
    target_os: Box<str>,
    requested_threads: usize,
    actual_threads: usize,
}

pub(crate) struct BackendFingerprintParts {
    pub(crate) schema_version: u32,
    pub(crate) crate_name: &'static str,
    pub(crate) crate_version: &'static str,
    pub(crate) features: [&'static str; 2],
    pub(crate) algorithm: &'static str,
    pub(crate) target_arch: &'static str,
    pub(crate) target_os: &'static str,
    pub(crate) requested_threads: usize,
    pub(crate) actual_threads: usize,
}

impl BackendFingerprint {
    pub(crate) fn new(parts: BackendFingerprintParts) -> Self {
        Self {
            schema_version: parts.schema_version,
            crate_name: parts.crate_name.into(),
            crate_version: parts.crate_version.into(),
            features: parts.features.map(Into::into),
            algorithm: parts.algorithm.into(),
            target_arch: parts.target_arch.into(),
            target_os: parts.target_os.into(),
            requested_threads: parts.requested_threads,
            actual_threads: parts.actual_threads,
        }
    }

    /// Returns the fingerprint schema version.
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Returns the backend crate name.
    pub fn crate_name(&self) -> &str {
        &self.crate_name
    }

    /// Returns the backend crate version.
    pub fn crate_version(&self) -> &str {
        &self.crate_version
    }

    /// Returns the exact enabled backend features.
    pub fn features(&self) -> [&str; 2] {
        self.features.each_ref().map(|feature| feature.as_ref())
    }

    /// Returns the resolved backend algorithm.
    pub fn algorithm(&self) -> &str {
        &self.algorithm
    }

    /// Returns the target architecture recorded by the adapter.
    pub fn target_arch(&self) -> &str {
        &self.target_arch
    }

    /// Returns the target operating system recorded by the adapter.
    pub fn target_os(&self) -> &str {
        &self.target_os
    }

    /// Returns the requested thread count recorded by the backend adapter.
    pub fn requested_threads(&self) -> usize {
        self.requested_threads
    }

    /// Returns the actual thread count recorded by the backend adapter.
    pub fn actual_threads(&self) -> usize {
        self.actual_threads
    }
}

/// Algorithm selected for one bounded backend attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolveAttemptKind {
    /// Symmetric-indefinite Bunch-Kaufman factorization with refinement.
    BunchKaufmanRefinement,
    /// Full-SVD rescue after the primary candidate was rejected.
    FullSvdRescue,
}

/// Resolved settings for one backend attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BackendAttemptSettings {
    kind: SolveAttemptKind,
    pivoting: Option<Box<str>>,
    block_size: Option<usize>,
    parallelism_threshold: Option<usize>,
    factor_workspace_source: Option<Box<str>>,
    maximum_refinement_steps: Option<usize>,
    settings_id: Option<Box<str>>,
    left_vectors: Option<Box<str>>,
    right_vectors: Option<Box<str>>,
}

impl BackendAttemptSettings {
    pub(crate) fn lblt(
        pivoting: impl Into<Box<str>>,
        block_size: usize,
        parallelism_threshold: usize,
        factor_workspace_source: impl Into<Box<str>>,
        maximum_refinement_steps: usize,
    ) -> Self {
        Self {
            kind: SolveAttemptKind::BunchKaufmanRefinement,
            pivoting: Some(pivoting.into()),
            block_size: Some(block_size),
            parallelism_threshold: Some(parallelism_threshold),
            factor_workspace_source: Some(factor_workspace_source.into()),
            maximum_refinement_steps: Some(maximum_refinement_steps),
            settings_id: None,
            left_vectors: None,
            right_vectors: None,
        }
    }

    pub(crate) fn full_svd(
        settings_id: impl Into<Box<str>>,
        left_vectors: impl Into<Box<str>>,
        right_vectors: impl Into<Box<str>>,
    ) -> Self {
        Self {
            kind: SolveAttemptKind::FullSvdRescue,
            pivoting: None,
            block_size: None,
            parallelism_threshold: None,
            factor_workspace_source: None,
            maximum_refinement_steps: None,
            settings_id: Some(settings_id.into()),
            left_vectors: Some(left_vectors.into()),
            right_vectors: Some(right_vectors.into()),
        }
    }

    /// Returns the algorithm family these settings configure.
    pub fn kind(&self) -> SolveAttemptKind {
        self.kind
    }

    /// Returns the pivoting strategy for a factorization attempt.
    pub fn pivoting(&self) -> Option<&str> {
        self.pivoting.as_deref()
    }

    /// Returns the factorization block size when applicable.
    pub fn block_size(&self) -> Option<usize> {
        self.block_size
    }

    /// Returns the backend parallelism threshold when applicable.
    pub fn parallelism_threshold(&self) -> Option<usize> {
        self.parallelism_threshold
    }

    /// Returns how factorization workspace was provisioned.
    pub fn factor_workspace_source(&self) -> Option<&str> {
        self.factor_workspace_source.as_deref()
    }

    /// Returns the maximum refinement-step budget when applicable.
    pub fn maximum_refinement_steps(&self) -> Option<usize> {
        self.maximum_refinement_steps
    }

    /// Returns the backend settings identity for an SVD rescue.
    pub fn settings_id(&self) -> Option<&str> {
        self.settings_id.as_deref()
    }

    /// Returns the requested left-singular-vector mode when applicable.
    pub fn left_vectors(&self) -> Option<&str> {
        self.left_vectors.as_deref()
    }

    /// Returns the requested right-singular-vector mode when applicable.
    pub fn right_vectors(&self) -> Option<&str> {
        self.right_vectors.as_deref()
    }
}

/// Ruiz scaling summary applied to one backend attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalingSummary {
    method: Box<str>,
    rounds: usize,
    saturated_outside_target: usize,
}

impl ScalingSummary {
    pub(crate) fn new(
        method: impl Into<Box<str>>,
        rounds: usize,
        saturated_outside_target: usize,
    ) -> Self {
        Self {
            method: method.into(),
            rounds,
            saturated_outside_target,
        }
    }

    /// Returns the stable scaling-method identity.
    pub fn method(&self) -> &str {
        &self.method
    }

    /// Returns the number of completed fixed scaling rounds.
    pub fn rounds(&self) -> usize {
        self.rounds
    }

    /// Returns how many rows saturated outside the target norm band.
    pub fn saturated_outside_target(&self) -> usize {
        self.saturated_outside_target
    }
}

/// Complete dimensionless residual evidence for a backend candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearResidualEvidence {
    infinity_norm: f64,
    matrix_infinity_norm: f64,
    solution_infinity_norm: f64,
    rhs_infinity_norm: f64,
    normalized_backward_error: f64,
}

impl LinearResidualEvidence {
    pub(crate) fn new(values: [f64; 5]) -> Self {
        Self {
            infinity_norm: values[0],
            matrix_infinity_norm: values[1],
            solution_infinity_norm: values[2],
            rhs_infinity_norm: values[3],
            normalized_backward_error: values[4],
        }
    }

    /// Returns the candidate residual infinity norm.
    pub fn infinity_norm(self) -> f64 {
        self.infinity_norm
    }

    /// Returns the backend-standard-form matrix infinity norm.
    pub fn matrix_infinity_norm(self) -> f64 {
        self.matrix_infinity_norm
    }

    /// Returns the backend candidate infinity norm.
    pub fn solution_infinity_norm(self) -> f64 {
        self.solution_infinity_norm
    }

    /// Returns the right-hand-side infinity norm.
    pub fn rhs_infinity_norm(self) -> f64 {
        self.rhs_infinity_norm
    }

    /// Returns the dimensionless normalized backward error.
    pub fn normalized_backward_error(self) -> f64 {
        self.normalized_backward_error
    }
}

/// Category of a rejected backend attempt or execution path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AttemptFailureCategory {
    /// The backend returned a non-finite candidate.
    NonFiniteCandidate,
    /// The candidate exceeded the normalized backward-error limit.
    BackwardErrorExceeded,
    /// Scaling recovery exceeded its round-trip limit.
    ScalingRoundTripExceeded,
    /// The backend decomposition itself failed numerically.
    BackendDecompositionFailure,
}

/// Structured reason why an attempt or execution path was rejected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AttemptFailureEvidence {
    category: AttemptFailureCategory,
    observed: Option<f64>,
    limit: Option<f64>,
}

impl AttemptFailureEvidence {
    pub(crate) fn new(
        category: AttemptFailureCategory,
        observed: Option<f64>,
        limit: Option<f64>,
    ) -> Self {
        Self {
            category,
            observed,
            limit,
        }
    }

    /// Returns the stable failure category.
    pub fn category(self) -> AttemptFailureCategory {
        self.category
    }

    /// Returns the observed value for thresholded failures.
    pub fn observed(self) -> Option<f64> {
        self.observed
    }

    /// Returns the applied acceptance limit for thresholded failures.
    pub fn limit(self) -> Option<f64> {
        self.limit
    }
}

/// Evidence retained for one bounded backend attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveAttemptRecord {
    sequence: usize,
    kind: SolveAttemptKind,
    termination: SolveAttemptTermination,
    settings: BackendAttemptSettings,
    scaling: ScalingSummary,
    refinement_steps: usize,
    residual: Option<LinearResidualEvidence>,
    certificate_present: bool,
    failure_reason: Option<AttemptFailureEvidence>,
    backend_fingerprint: BackendFingerprint,
}

pub(crate) struct SolveAttemptRecordParts {
    pub(crate) sequence: usize,
    pub(crate) kind: SolveAttemptKind,
    pub(crate) termination: SolveAttemptTermination,
    pub(crate) settings: BackendAttemptSettings,
    pub(crate) scaling: ScalingSummary,
    pub(crate) refinement_steps: usize,
    pub(crate) residual: Option<LinearResidualEvidence>,
    pub(crate) certificate_present: bool,
    pub(crate) failure_reason: Option<AttemptFailureEvidence>,
    pub(crate) backend_fingerprint: BackendFingerprint,
}

impl SolveAttemptRecord {
    pub(crate) fn new(parts: SolveAttemptRecordParts) -> Self {
        Self {
            sequence: parts.sequence,
            kind: parts.kind,
            termination: parts.termination,
            settings: parts.settings,
            scaling: parts.scaling,
            refinement_steps: parts.refinement_steps,
            residual: parts.residual,
            certificate_present: parts.certificate_present,
            failure_reason: parts.failure_reason,
            backend_fingerprint: parts.backend_fingerprint,
        }
    }

    /// Returns the deterministic attempt sequence number.
    pub fn sequence(&self) -> usize {
        self.sequence
    }

    /// Returns the algorithm family used for this attempt.
    pub fn kind(&self) -> SolveAttemptKind {
        self.kind
    }

    /// Returns backend termination evidence for this attempt.
    pub fn termination(&self) -> SolveAttemptTermination {
        self.termination
    }

    /// Returns the complete resolved attempt settings.
    pub fn settings(&self) -> &BackendAttemptSettings {
        &self.settings
    }

    /// Returns the scaling summary applied to this attempt.
    pub fn scaling(&self) -> &ScalingSummary {
        &self.scaling
    }

    /// Returns the number of completed iterative-refinement steps.
    pub fn refinement_steps(&self) -> usize {
        self.refinement_steps
    }

    /// Returns complete residual evidence when a candidate was produced.
    pub fn residual(&self) -> Option<LinearResidualEvidence> {
        self.residual
    }

    /// Returns normalized backward error when the attempt produced a candidate.
    pub fn normalized_backward_error(&self) -> Option<f64> {
        self.residual
            .map(|evidence| evidence.normalized_backward_error())
    }

    /// Reports whether this attempt retained a backend certificate.
    pub fn certificate_present(&self) -> bool {
        self.certificate_present
    }

    /// Returns structured rejection evidence for a failed attempt.
    pub fn failure_reason(&self) -> Option<AttemptFailureEvidence> {
        self.failure_reason
    }

    /// Returns the adapter-recorded backend identity and settings.
    pub fn backend_fingerprint(&self) -> &BackendFingerprint {
        &self.backend_fingerprint
    }
}

/// Structured physical evidence retained when recovery rejects a candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct RecoveryVerificationEvidence {
    reasons: Vec<RecoveryVerificationReason>,
    field_value_hard_residual_max: Option<f64>,
    field_value_per_length_hard_residual_max: Option<f64>,
    polynomial_round_trip_error: Option<f64>,
    field_coefficient_round_trip_error: Option<f64>,
    field_energy_round_trip_error: Option<f64>,
    tolerance_round_trip_error: Option<f64>,
    no_model_produced: bool,
}

impl RecoveryVerificationEvidence {
    pub(crate) fn new(
        reasons: Vec<RecoveryVerificationReason>,
        hard_residual_maxima: Option<(f64, f64)>,
        polynomial_round_trip_error: Option<f64>,
        field_coefficient_round_trip_error: Option<f64>,
        field_energy_round_trip_error: Option<f64>,
        tolerance_round_trip_error: Option<f64>,
        no_model_produced: bool,
    ) -> Self {
        Self {
            reasons,
            field_value_hard_residual_max: hard_residual_maxima.map(|maxima| maxima.0),
            field_value_per_length_hard_residual_max: hard_residual_maxima.map(|maxima| maxima.1),
            polynomial_round_trip_error,
            field_coefficient_round_trip_error,
            field_energy_round_trip_error,
            tolerance_round_trip_error,
            no_model_produced,
        }
    }

    /// Returns every recovery rejection reason in deterministic order.
    pub fn reasons(&self) -> &[RecoveryVerificationReason] {
        &self.reasons
    }

    /// Returns the maximum hard field-value residual when recovery reached it.
    pub fn field_value_hard_residual_max(&self) -> Option<f64> {
        self.field_value_hard_residual_max
    }

    /// Returns the maximum hard derivative residual when recovery reached it.
    pub fn field_value_per_length_hard_residual_max(&self) -> Option<f64> {
        self.field_value_per_length_hard_residual_max
    }

    /// Returns polynomial round-trip error when available.
    pub fn polynomial_round_trip_error(&self) -> Option<f64> {
        self.polynomial_round_trip_error
    }

    /// Returns field-coefficient round-trip error when available.
    pub fn field_coefficient_round_trip_error(&self) -> Option<f64> {
        self.field_coefficient_round_trip_error
    }

    /// Returns FieldEnergy round-trip error when available.
    pub fn field_energy_round_trip_error(&self) -> Option<f64> {
        self.field_energy_round_trip_error
    }

    /// Returns relation-tolerance round-trip error when available.
    pub fn tolerance_round_trip_error(&self) -> Option<f64> {
        self.tolerance_round_trip_error
    }

    /// Confirms that rejection produced no public model.
    pub fn no_model_produced(&self) -> bool {
        self.no_model_produced
    }
}
