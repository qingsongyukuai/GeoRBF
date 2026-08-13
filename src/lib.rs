#![forbid(unsafe_code)]
//! GeoRBF is a pure-Rust port of the non-visual core of Surfe.
//!
//! The crate currently exposes Surfe-compatible parameters, stable error
//! categories, and safe geological constraint values. Mathematical and
//! modelling modules are added in the dependency order fixed by the port plan.

mod constraints;
mod error;
mod geometry;
mod parameters;

pub use constraints::{Constraints, Inequality, Interface, Planar, Polarity, Tangent};
pub use error::Error;
pub use geometry::{ConstraintError, Point};
pub use parameters::{
    Axis, DerivativePoint, FirstDerivative, InputParameters, InternalParameters, ModelType,
    Parameters, RbfKernel, SecondDerivative, SolverType, DEGREES_TO_RADIANS, POSITION_EPSILON,
    RADIANS_TO_DEGREES,
};
