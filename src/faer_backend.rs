use faer::dyn_stack::StackReq;
use faer::linalg::cholesky::lblt::{
    factor::{self, LbltParams, PivotingStrategy},
    solve,
};
use faer::{Par, Spec};

pub(crate) const CRATE_NAME: &str = "faer";
pub(crate) const CRATE_VERSION: &str = "0.24.4";
pub(crate) const FEATURES: [&str; 2] = ["linalg", "std"];
pub(crate) const ALGORITHM: &str = "LBLT Bunch-Kaufman";
pub(crate) const PIVOTING: &str = "PartialDiag";
pub(crate) const BLOCK_SIZE: usize = 64;
pub(crate) const PARALLELISM_THRESHOLD: usize = 128 * 128;
pub(crate) const REQUESTED_THREADS: usize = 1;
pub(crate) const FACTOR_WORKSPACE_SOURCE: &str =
    "faer::linalg::cholesky::lblt::factor::cholesky_in_place_scratch";

pub(crate) fn parallelism() -> Par {
    Par::Seq
}

pub(crate) fn lblt_params() -> Spec<LbltParams, f64> {
    // faer's non-exhaustive marker prevents direct construction. Start from
    // its exact-version value, then overwrite every behavioral field so no
    // backend default defines GeoRBF's policy.
    let mut params = Spec::<LbltParams, f64>::default();
    params.pivoting = PivotingStrategy::PartialDiag;
    params.block_size = BLOCK_SIZE;
    params.par_threshold = PARALLELISM_THRESHOLD;
    params
}

pub(crate) fn factor_workspace_requirement(dimension: usize) -> StackReq {
    factor::cholesky_in_place_scratch::<usize, f64>(dimension, parallelism(), lblt_params())
}

pub(crate) fn solve_workspace_requirement(dimension: usize) -> StackReq {
    solve::solve_in_place_scratch::<usize, f64>(dimension, 1, parallelism())
}
