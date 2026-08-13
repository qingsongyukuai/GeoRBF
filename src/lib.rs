#![forbid(unsafe_code)]
//! GeoRBF is a pure-Rust port of the non-visual core of Surfe.
//!
//! The crate currently exposes Surfe-compatible parameters, stable error
//! categories, safe geological constraint values, and deterministic constraint
//! cleaning/grouping. Mathematical and modelling modules are added in the
//! dependency order fixed by the port plan.

mod constraints;
mod error;
mod geometry;
mod ordering;
mod parameters;

pub use constraints::{
    CollocationRemoval, Constraints, Inequality, Interface, InterfaceGrouping, Planar, Polarity,
    Tangent,
};
pub use error::Error;
pub use geometry::{ConstraintError, Point};
pub use ordering::{collocated, compare_points, sort_values_with_indices};
pub use parameters::{
    Axis, DerivativePoint, FirstDerivative, InputParameters, InternalParameters, ModelType,
    Parameters, RbfKernel, SecondDerivative, SolverType, DEGREES_TO_RADIANS, POSITION_EPSILON,
    RADIANS_TO_DEGREES,
};
