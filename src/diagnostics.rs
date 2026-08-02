//! Typed public fit diagnoses and backend-attempt evidence.

pub use crate::functional::SemanticRolePath;
pub use crate::numerical::NumericalPolicyId;

/// GeoRBF's semantic conclusion for a failed fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProblemDiagnosis {
    /// The supplied problem could not form a valid fitting problem.
    InvalidProblem,
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
    crate_name: Box<str>,
    crate_version: Box<str>,
    algorithm: Box<str>,
    requested_threads: usize,
    actual_threads: usize,
}

impl BackendFingerprint {
    pub(crate) fn new(
        crate_name: impl Into<Box<str>>,
        crate_version: impl Into<Box<str>>,
        algorithm: impl Into<Box<str>>,
        requested_threads: usize,
        actual_threads: usize,
    ) -> Self {
        Self {
            crate_name: crate_name.into(),
            crate_version: crate_version.into(),
            algorithm: algorithm.into(),
            requested_threads,
            actual_threads,
        }
    }

    /// Returns the backend crate name.
    pub fn crate_name(&self) -> &str {
        &self.crate_name
    }

    /// Returns the backend crate version.
    pub fn crate_version(&self) -> &str {
        &self.crate_version
    }

    /// Returns the resolved backend algorithm.
    pub fn algorithm(&self) -> &str {
        &self.algorithm
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

/// Evidence retained for one bounded backend attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveAttemptRecord {
    sequence: usize,
    termination: SolveAttemptTermination,
    normalized_backward_error: Option<f64>,
    backend_fingerprint: BackendFingerprint,
}

impl SolveAttemptRecord {
    pub(crate) fn new(
        sequence: usize,
        termination: SolveAttemptTermination,
        normalized_backward_error: Option<f64>,
        backend_fingerprint: BackendFingerprint,
    ) -> Self {
        Self {
            sequence,
            termination,
            normalized_backward_error,
            backend_fingerprint,
        }
    }

    /// Returns the deterministic attempt sequence number.
    pub fn sequence(&self) -> usize {
        self.sequence
    }

    /// Returns backend termination evidence for this attempt.
    pub fn termination(&self) -> SolveAttemptTermination {
        self.termination
    }

    /// Returns normalized backward error when the attempt produced a candidate.
    pub fn normalized_backward_error(&self) -> Option<f64> {
        self.normalized_backward_error
    }

    /// Returns the adapter-recorded backend identity and settings.
    pub fn backend_fingerprint(&self) -> &BackendFingerprint {
        &self.backend_fingerprint
    }
}

/// A physical Recover-and-Verify rejection reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RecoveryVerificationReason {
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
    /// Relation-tolerance recovery exceeded its round-trip limit.
    ToleranceRoundTripViolation,
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
