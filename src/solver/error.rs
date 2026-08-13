//! Typed failures for the pure-Rust partial-pivot LU path.

use std::fmt;

use crate::{DenseVector, Error};

use super::{LuFactorizationEvidence, LuResidualEvidence, LuValidation};

/// Precise failure stage for one LU request.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum LuSolveErrorKind {
    InvalidStorageLength,
    EmptySystem,
    NonSquareMatrix,
    DimensionMismatch,
    NonFiniteMatrix,
    NonFiniteRightHandSide,
    SingularSystem,
    NonFiniteSolution,
    NonFiniteResidual,
    ResidualTooLarge,
}

impl fmt::Display for LuSolveErrorKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::InvalidStorageLength => "matrix storage length does not match its dimensions",
            Self::EmptySystem => "the linear system is empty",
            Self::NonSquareMatrix => "partial-pivot LU requires a square matrix",
            Self::DimensionMismatch => "right-hand-side length does not match matrix rows",
            Self::NonFiniteMatrix => "the interpolation matrix contains a non-finite value",
            Self::NonFiniteRightHandSide => {
                "the linear right-hand side contains a non-finite value"
            }
            Self::SingularSystem => "an exact zero pivot produced a non-finite solution",
            Self::NonFiniteSolution => "LU produced a non-finite solution",
            Self::NonFiniteResidual => "the solved system produced a non-finite residual",
            Self::ResidualTooLarge => "the solved system failed its backward-error check",
        })
    }
}

/// Failure plus evidence showing whether frozen-style factorization was tried.
#[derive(Clone, Debug, PartialEq)]
pub struct LuSolveError {
    pub(super) kind: LuSolveErrorKind,
    pub(super) attempted: bool,
    pub(super) validation: LuValidation,
    pub(super) factorization: Option<Box<LuFactorizationEvidence>>,
    pub(super) candidate_weights: Option<Box<DenseVector>>,
    pub(super) residual: Option<LuResidualEvidence>,
}

impl LuSolveError {
    pub const fn kind(&self) -> LuSolveErrorKind {
        self.kind
    }

    pub const fn attempted(&self) -> bool {
        self.attempted
    }

    pub const fn validation(&self) -> LuValidation {
        self.validation
    }

    pub fn factorization(&self) -> Option<&LuFactorizationEvidence> {
        self.factorization.as_deref()
    }

    pub fn candidate_weights(&self) -> Option<&DenseVector> {
        self.candidate_weights.as_deref()
    }

    pub fn residual(&self) -> Option<&LuResidualEvidence> {
        self.residual.as_ref()
    }

    /// Frozen model call sites map any `solve() == false` to this category.
    pub const fn surfe_error(&self) -> Error {
        Error::LinearSolverFailure
    }
}

impl fmt::Display for LuSolveError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.kind.fmt(formatter)
    }
}

impl std::error::Error for LuSolveError {}
