#![forbid(unsafe_code)]
//! GeoRBF is a pure-Rust port of the non-visual core of Surfe.
//!
//! The crate currently exposes Surfe-compatible parameters, stable error
//! categories, safe geological constraint values, and deterministic constraint
//! cleaning/grouping, spatial helpers, fixed-order polynomial bases, and
//! isotropic radial kernels, frozen global anisotropy, modified kernels, and
//! model-independent linear functionals plus deterministic five-model
//! row/column layouts, dense system assembly, and partial-pivot LU solving.
//! The ordinary and restricted-range LOQO-style predictor-corrector
//! quadratic-programming paths are also available with iteration and
//! feasibility evidence.
//! Mathematical and modelling modules are added in the dependency order fixed
//! by the port plan.

mod assembly;
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
mod solver;
mod spatial;

pub use assembly::{
    assemble_system, AssembledSystem, AssemblyConstraints, AssemblyError, BoundedConstraintSystem,
    ConstraintSystem, DenseMatrix, DenseVector,
};
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
pub use solver::{
    loqo_step_divisor, predictor_corrector_step_length, solve_dense_partial_pivot_lu,
    solve_loqo_qp, solve_loqo_qp_with_options, solve_partial_pivot_lu,
    solve_predictor_corrector_qp, solve_predictor_corrector_qp_with_options, validate_loqo_qp,
    validate_lu_system, validate_predictor_corrector_qp, LoqoIterationEvidence, LoqoKktFailure,
    LoqoKktSolveEvidence, LoqoKktStage, LoqoOptions, LoqoResidualEvidence, LoqoSolution,
    LoqoSolveError, LoqoSolveErrorKind, LoqoStepError, LoqoStopReason, LoqoValidation,
    LuFactorizationEvidence, LuResidualEvidence, LuSolution, LuSolveError, LuSolveErrorKind,
    LuValidation, QpIterationEvidence, QpKktFailure, QpKktSolveEvidence, QpKktStage, QpOptions,
    QpResidualEvidence, QpSolution, QpSolveError, QpSolveErrorKind, QpStepLengthError,
    QpStopReason, QpValidation,
};
pub use spatial::{
    average_nearest_neighbour_distance, bounds, closest_to_distance_index, constraints_to_points,
    distance_between_points, extremal_point_indices, farthest_from_other_set_index,
    farthest_neighbour_index, farthest_pair_indices, largest_distance_between_points,
    maximal_axial_variability_order, nearest_neighbour_index, nearest_neighbour_indices,
    spatial_metrics, ConstraintAverageNearestNeighbourDistances, SpatialError, SpatialParameters,
};
