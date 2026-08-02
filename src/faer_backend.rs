use faer::diag::DiagMut;
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::linalg::cholesky::lblt::{
    factor::{self, LbltParams, PivotingStrategy},
    solve,
};
use faer::linalg::cholesky::llt::factor as llt;
use faer::linalg::evd::{self, ComputeEigenvectors};
use faer::linalg::qr::{col_pivoting, no_pivoting};
use faer::linalg::svd::{self, ComputeSvdVectors};
use faer::{Mat, MatRef, Par, Spec};

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
pub(crate) const RRQR_ALGORITHM: &str = "faer column-pivoted Householder QR";
pub(crate) const RRQR_SETTINGS_ID: &str =
    "georbf-faer-rrqr-v1:block-size=256,blocking=2304,parallel=49152";
pub(crate) const SVD_ALGORITHM: &str = "faer bidiagonal SVD";
pub(crate) const SVD_SETTINGS_ID: &str = "georbf-faer-svd-v1:recursion=128,qr-ratio=11/6,bidiag-parallel=49152,qr-blocking=2304,qr-parallel=49152";
pub(crate) const EVD_ALGORITHM: &str = "faer self-adjoint tridiagonal EVD";
pub(crate) const EVD_SETTINGS_ID: &str = "georbf-faer-evd-v1:recursion=128,tridiag-parallel=49152";
pub(crate) const HOUSEHOLDER_QR_SETTINGS_ID: &str =
    "georbf-faer-qr-v1:block-size=256,blocking=2304,parallel=49152";
pub(crate) const LLT_SETTINGS_ID: &str =
    "georbf-faer-llt-v1:recursion=64,block-size=128,regularization=disabled";

const BLOCKING_THRESHOLD: usize = 48 * 48;
const ANALYSIS_PARALLELISM_THRESHOLD: usize = 192 * 256;
const SVD_RECURSION_THRESHOLD: usize = 128;
const SVD_QR_RATIO_THRESHOLD: f64 = 11.0 / 6.0;
const EVD_RECURSION_THRESHOLD: usize = 128;
const QR_BLOCK_SIZE: usize = 256;
const LLT_RECURSION_THRESHOLD: usize = 64;
const LLT_BLOCK_SIZE: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct WorkspaceAllocationFailure {
    pub(crate) bytes: u64,
    pub(crate) alignment: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DecompositionFailure {
    WorkspaceAllocation(WorkspaceAllocationFailure),
    NumericalError,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CholeskyFailure {
    WorkspaceAllocation(WorkspaceAllocationFailure),
    NonPositivePivot,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BunchKaufmanInertia {
    pub(crate) positive: usize,
    pub(crate) negative: usize,
    pub(crate) zero: usize,
}

pub(crate) struct HouseholderQrFactors {
    pub(crate) basis: Mat<f64>,
    pub(crate) coefficients: Mat<f64>,
}

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

pub(crate) fn rrqr_params() -> Spec<col_pivoting::factor::ColPivQrParams, f64> {
    let mut params = Spec::<col_pivoting::factor::ColPivQrParams, f64>::default();
    params.blocking_threshold = BLOCKING_THRESHOLD;
    params.par_threshold = ANALYSIS_PARALLELISM_THRESHOLD;
    params
}

pub(crate) fn qr_params() -> Spec<no_pivoting::factor::QrParams, f64> {
    let mut params = Spec::<no_pivoting::factor::QrParams, f64>::default();
    params.blocking_threshold = BLOCKING_THRESHOLD;
    params.par_threshold = ANALYSIS_PARALLELISM_THRESHOLD;
    params
}

pub(crate) fn svd_params() -> Spec<svd::SvdParams, f64> {
    let mut params = Spec::<svd::SvdParams, f64>::default();
    params.bidiag.par_threshold = ANALYSIS_PARALLELISM_THRESHOLD;
    params.qr.blocking_threshold = BLOCKING_THRESHOLD;
    params.qr.par_threshold = ANALYSIS_PARALLELISM_THRESHOLD;
    params.recursion_threshold = SVD_RECURSION_THRESHOLD;
    params.qr_ratio_threshold = SVD_QR_RATIO_THRESHOLD;
    params
}

pub(crate) fn evd_params() -> Spec<evd::SelfAdjointEvdParams, f64> {
    let mut params = Spec::<evd::SelfAdjointEvdParams, f64>::default();
    params.tridiag.par_threshold = ANALYSIS_PARALLELISM_THRESHOLD;
    params.recursion_threshold = EVD_RECURSION_THRESHOLD;
    params
}

pub(crate) fn llt_params() -> Spec<llt::LltParams, f64> {
    let mut params = Spec::<llt::LltParams, f64>::default();
    params.recursion_threshold = LLT_RECURSION_THRESHOLD;
    params.block_size = LLT_BLOCK_SIZE;
    params
}

pub(crate) fn factor_workspace_requirement(dimension: usize) -> StackReq {
    factor::cholesky_in_place_scratch::<usize, f64>(dimension, parallelism(), lblt_params())
}

pub(crate) fn solve_workspace_requirement(dimension: usize) -> StackReq {
    solve::solve_in_place_scratch::<usize, f64>(dimension, 1, parallelism())
}

pub(crate) fn qr_block_size() -> usize {
    QR_BLOCK_SIZE
}

pub(crate) fn rrqr_workspace_requirement(dimension: usize) -> StackReq {
    col_pivoting::factor::qr_in_place_scratch::<usize, f64>(
        dimension,
        dimension,
        qr_block_size(),
        parallelism(),
        rrqr_params(),
    )
}

pub(crate) fn qr_workspace_requirement(dimension: usize) -> StackReq {
    no_pivoting::factor::qr_in_place_scratch::<f64>(
        dimension,
        dimension,
        qr_block_size(),
        parallelism(),
        qr_params(),
    )
}

pub(crate) fn llt_workspace_requirement(dimension: usize) -> StackReq {
    llt::cholesky_in_place_scratch::<f64>(dimension, parallelism(), llt_params())
}

pub(crate) fn singular_values_workspace_requirement(dimension: usize) -> StackReq {
    svd::svd_scratch::<f64>(
        dimension,
        dimension,
        ComputeSvdVectors::No,
        ComputeSvdVectors::No,
        parallelism(),
        svd_params(),
    )
}

pub(crate) fn inertia_workspace_requirement(dimension: usize) -> StackReq {
    evd::self_adjoint_evd_scratch::<f64>(
        dimension,
        ComputeEigenvectors::No,
        parallelism(),
        evd_params(),
    )
}

pub(crate) fn svd_rescue_workspace_requirement(dimension: usize) -> StackReq {
    svd::svd_scratch::<f64>(
        dimension,
        dimension,
        ComputeSvdVectors::Full,
        ComputeSvdVectors::Full,
        parallelism(),
        svd_params(),
    )
}

pub(crate) fn rrqr_diagonal(
    matrix: MatRef<'_, f64>,
) -> Result<Vec<f64>, WorkspaceAllocationFailure> {
    let (nrows, ncols) = matrix.shape();
    let size = nrows.min(ncols);
    let block_size = qr_block_size();
    let mut factors = Mat::from_fn(nrows, ncols, |row, column| matrix[(row, column)]);
    let mut householder = Mat::<f64>::zeros(block_size, size);
    let mut permutation = vec![0usize; ncols];
    let mut inverse_permutation = vec![0usize; ncols];
    let requirement = col_pivoting::factor::qr_in_place_scratch::<usize, f64>(
        nrows,
        ncols,
        block_size,
        parallelism(),
        rrqr_params(),
    );
    let mut memory = allocate(requirement)?;
    col_pivoting::factor::qr_in_place(
        factors.as_mut(),
        householder.as_mut(),
        &mut permutation,
        &mut inverse_permutation,
        parallelism(),
        MemStack::new(&mut memory),
        rrqr_params(),
    );
    Ok((0..size)
        .map(|index| factors[(index, index)].abs())
        .collect())
}

pub(crate) fn householder_qr(
    matrix: MatRef<'_, f64>,
) -> Result<HouseholderQrFactors, WorkspaceAllocationFailure> {
    let (nrows, ncols) = matrix.shape();
    let size = nrows.min(ncols);
    let block_size = qr_block_size();
    let mut basis = Mat::from_fn(nrows, ncols, |row, column| matrix[(row, column)]);
    let mut coefficients = Mat::<f64>::zeros(block_size, size);
    let requirement = no_pivoting::factor::qr_in_place_scratch::<f64>(
        nrows,
        ncols,
        block_size,
        parallelism(),
        qr_params(),
    );
    let mut memory = allocate(requirement)?;
    no_pivoting::factor::qr_in_place(
        basis.as_mut(),
        coefficients.as_mut(),
        parallelism(),
        MemStack::new(&mut memory),
        qr_params(),
    );
    Ok(HouseholderQrFactors {
        basis,
        coefficients,
    })
}

pub(crate) fn cholesky_minimum_diagonal(matrix: MatRef<'_, f64>) -> Result<f64, CholeskyFailure> {
    let dimension = matrix.nrows();
    let mut factors = Mat::from_fn(dimension, dimension, |row, column| matrix[(row, column)]);
    let requirement = llt_workspace_requirement(dimension);
    let mut memory = allocate(requirement).map_err(CholeskyFailure::WorkspaceAllocation)?;
    llt::cholesky_in_place(
        factors.as_mut(),
        llt::LltRegularization {
            dynamic_regularization_delta: 0.0,
            dynamic_regularization_epsilon: 0.0,
        },
        parallelism(),
        MemStack::new(&mut memory),
        llt_params(),
    )
    .map_err(|_| CholeskyFailure::NonPositivePivot)?;
    Ok((0..dimension)
        .map(|index| factors[(index, index)])
        .fold(f64::INFINITY, f64::min))
}

pub(crate) fn bunch_kaufman_inertia(
    matrix: MatRef<'_, f64>,
) -> Result<BunchKaufmanInertia, DecompositionFailure> {
    let dimension = matrix.nrows();
    let mut factors = Mat::from_fn(dimension, dimension, |row, column| matrix[(row, column)]);
    let mut subdiagonal = vec![0.0; dimension];
    let mut permutation = vec![0usize; dimension];
    let mut inverse_permutation = vec![0usize; dimension];
    let requirement = factor_workspace_requirement(dimension);
    let mut memory = allocate(requirement).map_err(DecompositionFailure::WorkspaceAllocation)?;
    factor::cholesky_in_place(
        factors.as_mut(),
        DiagMut::from_slice_mut(&mut subdiagonal),
        &mut permutation,
        &mut inverse_permutation,
        parallelism(),
        MemStack::new(&mut memory),
        lblt_params(),
    );
    let mut inertia = BunchKaufmanInertia {
        positive: 0,
        negative: 0,
        zero: 0,
    };
    let mut index = 0;
    while index < dimension {
        if index + 1 < dimension && subdiagonal[index] != 0.0 {
            let diagonal = [factors[(index, index)], factors[(index + 1, index + 1)]];
            let off_diagonal = subdiagonal[index];
            if diagonal
                .into_iter()
                .chain([off_diagonal])
                .any(|value| !value.is_finite())
            {
                return Err(DecompositionFailure::NumericalError);
            }
            let scale = diagonal[0]
                .abs()
                .max(diagonal[1].abs())
                .max(off_diagonal.abs());
            if scale == 0.0 {
                inertia.zero += 2;
            } else {
                let center = 0.5 * (diagonal[0] / scale + diagonal[1] / scale);
                let radius =
                    (0.5 * (diagonal[0] / scale - diagonal[1] / scale)).hypot(off_diagonal / scale);
                record_sign(center - radius, &mut inertia);
                record_sign(center + radius, &mut inertia);
            }
            index += 2;
        } else {
            let value = factors[(index, index)];
            if !value.is_finite() {
                return Err(DecompositionFailure::NumericalError);
            }
            record_sign(value, &mut inertia);
            index += 1;
        }
    }
    Ok(inertia)
}

fn record_sign(value: f64, inertia: &mut BunchKaufmanInertia) {
    if value > 0.0 {
        inertia.positive += 1;
    } else if value < 0.0 {
        inertia.negative += 1;
    } else {
        inertia.zero += 1;
    }
}

pub(crate) fn singular_values(matrix: MatRef<'_, f64>) -> Result<Vec<f64>, DecompositionFailure> {
    let (nrows, ncols) = matrix.shape();
    let mut values = vec![0.0; nrows.min(ncols)];
    let requirement = svd::svd_scratch::<f64>(
        nrows,
        ncols,
        ComputeSvdVectors::No,
        ComputeSvdVectors::No,
        parallelism(),
        svd_params(),
    );
    let mut memory = allocate(requirement).map_err(DecompositionFailure::WorkspaceAllocation)?;
    svd::svd(
        matrix,
        DiagMut::from_slice_mut(&mut values),
        None,
        None,
        parallelism(),
        MemStack::new(&mut memory),
        svd_params(),
    )
    .map_err(|_| DecompositionFailure::NumericalError)?;
    Ok(values)
}

pub(crate) fn self_adjoint_eigenvalues(
    matrix: MatRef<'_, f64>,
) -> Result<Vec<f64>, DecompositionFailure> {
    let dimension = matrix.nrows();
    let mut values = vec![0.0; dimension];
    let requirement = evd::self_adjoint_evd_scratch::<f64>(
        dimension,
        ComputeEigenvectors::No,
        parallelism(),
        evd_params(),
    );
    let mut memory = allocate(requirement).map_err(DecompositionFailure::WorkspaceAllocation)?;
    evd::self_adjoint_evd(
        matrix,
        DiagMut::from_slice_mut(&mut values),
        None,
        parallelism(),
        MemStack::new(&mut memory),
        evd_params(),
    )
    .map_err(|_| DecompositionFailure::NumericalError)?;
    Ok(values)
}

pub(crate) fn solve_with_full_svd(
    matrix: MatRef<'_, f64>,
    rhs: &[f64],
) -> Result<Vec<f64>, DecompositionFailure> {
    let dimension = matrix.nrows();
    let mut singular_values = vec![0.0; dimension];
    let mut left = Mat::<f64>::zeros(dimension, dimension);
    let mut right = Mat::<f64>::zeros(dimension, dimension);
    let requirement = svd_rescue_workspace_requirement(dimension);
    let mut memory = allocate(requirement).map_err(DecompositionFailure::WorkspaceAllocation)?;
    svd::svd(
        matrix,
        DiagMut::from_slice_mut(&mut singular_values),
        Some(left.as_mut()),
        Some(right.as_mut()),
        parallelism(),
        MemStack::new(&mut memory),
        svd_params(),
    )
    .map_err(|_| DecompositionFailure::NumericalError)?;
    if singular_values
        .iter()
        .any(|value| !value.is_finite() || *value == 0.0)
    {
        return Err(DecompositionFailure::NumericalError);
    }
    let projected = (0..dimension)
        .map(|column| {
            (0..dimension)
                .map(|row| left[(row, column)] * rhs[row])
                .sum::<f64>()
                / singular_values[column]
        })
        .collect::<Vec<_>>();
    Ok((0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| right[(row, column)] * projected[column])
                .sum()
        })
        .collect())
}

fn allocate(requirement: StackReq) -> Result<MemBuffer, WorkspaceAllocationFailure> {
    MemBuffer::try_new(requirement).map_err(|_| WorkspaceAllocationFailure {
        bytes: u64::try_from(requirement.size_bytes()).unwrap_or(u64::MAX),
        alignment: requirement.align_bytes(),
    })
}
