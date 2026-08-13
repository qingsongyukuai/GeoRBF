#![forbid(unsafe_code)]
//! GeoRBF is a pure-Rust port of the non-visual core of Surfe.
//!
//! The crate currently exposes Surfe-compatible parameters and stable error
//! categories. Mathematical and modelling modules are added in the dependency
//! order fixed by the port plan.

mod error;
mod parameters;

pub use error::Error;
pub use parameters::{
    Axis, DerivativePoint, FirstDerivative, InputParameters, InternalParameters, ModelType,
    Parameters, RbfKernel, SecondDerivative, SolverType, DEGREES_TO_RADIANS, POSITION_EPSILON,
    RADIANS_TO_DEGREES,
};
