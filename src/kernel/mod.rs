//! Pure-Rust radial kernels used by GeoRBF's later functional and model layers.

mod derivatives;
mod isotropic;

pub use isotropic::{IsotropicKernel, KernelError, KernelEvaluation};
