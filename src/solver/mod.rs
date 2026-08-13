//! Pure-Rust solver paths corresponding to frozen Surfe.

mod error;
mod lu;
mod qp_loqo;
mod qp_predictor_corrector;

pub use error::{LuSolveError, LuSolveErrorKind};
pub use lu::{
    solve_dense_partial_pivot_lu, solve_partial_pivot_lu, validate_lu_system,
    LuFactorizationEvidence, LuResidualEvidence, LuSolution, LuValidation,
};
pub use qp_loqo::{
    loqo_step_divisor, solve_loqo_qp, solve_loqo_qp_with_options, validate_loqo_qp,
    LoqoIterationEvidence, LoqoKktFailure, LoqoKktSolveEvidence, LoqoKktStage, LoqoOptions,
    LoqoResidualEvidence, LoqoSolution, LoqoSolveError, LoqoSolveErrorKind, LoqoStepError,
    LoqoStopReason, LoqoValidation,
};
pub use qp_predictor_corrector::{
    predictor_corrector_step_length, solve_predictor_corrector_qp,
    solve_predictor_corrector_qp_with_options, validate_predictor_corrector_qp,
    QpIterationEvidence, QpKktFailure, QpKktSolveEvidence, QpKktStage, QpOptions,
    QpResidualEvidence, QpSolution, QpSolveError, QpSolveErrorKind, QpStepLengthError,
    QpStopReason, QpValidation,
};
