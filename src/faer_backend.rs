use faer::dyn_stack::StackReq;
use faer::linalg::cholesky::lblt::{
    factor::{self, LbltParams, PivotingStrategy},
    solve,
};
use faer::linalg::evd::{self, ComputeEigenvectors};
use faer::linalg::qr::col_pivoting;
use faer::linalg::svd::{self, ComputeSvdVectors};
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

pub(crate) fn rrqr_block_size(dimension: usize) -> usize {
    faer::linalg::qr::no_pivoting::factor::recommended_block_size::<f64>(dimension, dimension)
}

pub(crate) fn rrqr_workspace_requirement(dimension: usize) -> StackReq {
    col_pivoting::factor::qr_in_place_scratch::<usize, f64>(
        dimension,
        dimension,
        rrqr_block_size(dimension),
        parallelism(),
        Default::default(),
    )
}

pub(crate) fn singular_values_workspace_requirement(dimension: usize) -> StackReq {
    svd::svd_scratch::<f64>(
        dimension,
        dimension,
        ComputeSvdVectors::No,
        ComputeSvdVectors::No,
        parallelism(),
        Default::default(),
    )
}

pub(crate) fn inertia_workspace_requirement(dimension: usize) -> StackReq {
    evd::self_adjoint_evd_scratch::<f64>(
        dimension,
        ComputeEigenvectors::No,
        parallelism(),
        Default::default(),
    )
}

pub(crate) fn svd_rescue_workspace_requirement(dimension: usize) -> StackReq {
    svd::svd_scratch::<f64>(
        dimension,
        dimension,
        ComputeSvdVectors::Full,
        ComputeSvdVectors::Full,
        parallelism(),
        Default::default(),
    )
}
