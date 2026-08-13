//! Pure-Rust solver paths corresponding to frozen Surfe.

mod error;
mod lu;

pub use error::{LuSolveError, LuSolveErrorKind};
pub use lu::{
    solve_dense_partial_pivot_lu, solve_partial_pivot_lu, validate_lu_system,
    LuFactorizationEvidence, LuResidualEvidence, LuSolution, LuValidation,
};
