#![forbid(unsafe_code)]
//! GeoRBF is a pure-Rust port of the non-visual core of Surfe.
//!
//! The crate currently exposes Surfe-compatible parameters, stable error
//! categories, safe geological constraint values, and deterministic constraint
//! cleaning/grouping, spatial helpers, fixed-order polynomial bases, and
//! isotropic radial kernels, frozen global anisotropy, modified kernels, and
//! model-independent linear functionals plus deterministic five-model
//! row/column layouts.
//! Mathematical and modelling modules are added in the dependency order fixed
//! by the port plan.

mod constraints;
mod error;
mod functional;
mod geometry;
mod kernel;
mod layout;
mod model;
mod ordering;
mod parameters;
mod polynomial;
mod spatial;

pub use constraints::{
    CollocationRemoval, Constraints, Inequality, Interface, InterfaceGrouping, Planar, Polarity,
    Tangent,
};
pub use error::Error;
pub use functional::{
    DofLabel, FunctionalKernel, FunctionalPrimitive, FunctionalTerm, LinearFunctional,
};
pub use geometry::{ConstraintError, Point};
pub use kernel::{
    AnisotropicKernel, AnisotropyError, IsotropicKernel, KernelError, KernelEvaluation,
    ModifiedKernel,
};
pub use layout::{
    constraint_layout, ConstraintLayout, DifferenceKind, IndexRange, LayoutDof, LayoutPartitions,
    LayoutPointRef, LayoutRole, LayoutSection, LayoutSectionKind, SourceConstraintCounts,
};
pub use ordering::{collocated, compare_points, sort_values_with_indices};
pub use parameters::{
    Axis, DerivativePoint, FirstDerivative, InputParameters, InternalParameters, ModelType,
    Parameters, RbfKernel, SecondDerivative, SolverType, DEGREES_TO_RADIANS, POSITION_EPSILON,
    RADIANS_TO_DEGREES,
};
pub use polynomial::{LagrangianPolynomialBasis, PolynomialBasis, PolynomialOrder};
pub use spatial::{
    average_nearest_neighbour_distance, bounds, closest_to_distance_index, constraints_to_points,
    distance_between_points, extremal_point_indices, farthest_from_other_set_index,
    farthest_neighbour_index, farthest_pair_indices, largest_distance_between_points,
    maximal_axial_variability_order, nearest_neighbour_index, nearest_neighbour_indices,
    spatial_metrics, ConstraintAverageNearestNeighbourDistances, SpatialError, SpatialParameters,
};
