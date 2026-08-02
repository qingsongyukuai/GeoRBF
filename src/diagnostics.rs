//! Typed public fit diagnoses and backend-attempt evidence.

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
}

impl SolveAttemptRecord {
    pub(crate) fn new(
        sequence: usize,
        termination: SolveAttemptTermination,
        normalized_backward_error: Option<f64>,
    ) -> Self {
        Self {
            sequence,
            termination,
            normalized_backward_error,
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
}
