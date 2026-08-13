#![forbid(unsafe_code)]
//! GeoRBF is a pure-Rust port of the non-visual core of Surfe.
//!
//! [`Builder`] is the safe public fitting entry point. It owns parameters and
//! geological constraints, validates dynamic compatibility tables, and returns
//! an immutable [`FittedModel`] snapshot. Fitted models dispatch all five Surfe
//! model kinds and support scalar and gradient evaluation for one point or an
//! ordered batch. Because evaluation has no shared mutable kernel state, a
//! fitted model can be read concurrently through ordinary Rust `Send + Sync`
//! ownership.
//!
//! Public failures expose their frozen Surfe exception class through
//! [`BuildError::surfe_category`] and [`EvaluationError::surfe_category`]. A
//! `None` category identifies a deliberate safe-Rust rejection that had no
//! stable C++ exception, so callers never need to classify diagnostic text.
//!
//! Lower-level model, assembly, solver, kernel, polynomial, and spatial APIs
//! remain available for parity inspection and advanced use.

mod assembly;
mod builder;
mod constraints;
mod error;
mod functional;
mod geometry;
mod greedy;
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
pub use builder::{BuildError, Builder};
pub use constraints::{
    CollocationRemoval, Constraints, Inequality, Interface, InterfaceGrouping, Planar, Polarity,
    Tangent,
};
pub use error::Error;
pub use functional::{
    DofLabel, FunctionalKernel, FunctionalPrimitive, FunctionalTerm, LinearFunctional,
};
pub use geometry::{ConstraintError, Point};
pub use greedy::{
    GreedyHookBody, GreedyModelAudit, GreedyRoundEvidence, GreedyStopReason, GreedyTrace,
    GREEDY_MODEL_AUDIT,
};
pub use kernel::{
    AnisotropicKernel, AnisotropyError, IsotropicKernel, KernelError, KernelEvaluation,
    ModifiedKernel,
};
pub use layout::{
    constraint_layout, ConstraintLayout, DifferenceKind, IndexRange, LayoutDof, LayoutPartitions,
    LayoutPointRef, LayoutRole, LayoutSection, LayoutSectionKind, SourceConstraintCounts,
};
pub use model::continuous_property::{
    fit_continuous_property, ContinuousPropertyError, ContinuousPropertyModel,
};
pub use model::fitted::{EvaluationError, FitBranch, FittedModel};
pub use model::lajaunie::{
    fit_lajaunie_linear, fit_lajaunie_restricted, fit_lajaunie_restricted_with_options,
    LajaunieIsoValueEvidence, LajaunieLinearError, LajaunieLinearModel,
    LajaunieRestrictedBoundEvidence, LajaunieRestrictedError, LajaunieRestrictedModel,
};
pub use model::reconstruct::{
    reconstruct_from_qp_weights, solve_and_reconstruct, ReconstructionAssemblyError,
    ReconstructionDofMapping, ReconstructionError, ReconstructionPredictionWitness,
    ReconstructionResult, ReconstructionSourceSolution, ReconstructionStage,
};
pub use model::single_surface::{
    fit_single_surface_inequality, fit_single_surface_inequality_with_options,
    fit_single_surface_linear, fit_single_surface_restricted,
    fit_single_surface_restricted_with_options, SingleSurfaceInequalityError,
    SingleSurfaceInequalityEvidence, SingleSurfaceInequalityModel, SingleSurfaceLinearError,
    SingleSurfaceLinearModel, SingleSurfaceRestrictedBoundEvidence, SingleSurfaceRestrictedError,
    SingleSurfaceRestrictedModel,
};
pub use model::stratigraphic::{
    fit_stratigraphic, fit_stratigraphic_restricted, fit_stratigraphic_restricted_with_options,
    fit_stratigraphic_with_options, StratigraphicError, StratigraphicIsoValueEvidence,
    StratigraphicLayerRelationEvidence, StratigraphicModel, StratigraphicRestrictedBoundEvidence,
    StratigraphicRestrictedError, StratigraphicRestrictedModel,
};
pub use model::vector_field::{fit_vector_field, VectorFieldError, VectorFieldModel};
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
