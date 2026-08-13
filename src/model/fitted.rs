//! Immutable public fitted-model dispatch over the five migrated model kinds.
//!
//! Source: `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`ComputeInterpolant`, getters, and scalar/vector single/batch evaluation).

use std::fmt;

use crate::{
    builder::{
        continuous_property_category, inequality_constraint_matrix, interface_constraint_matrix,
        lajaunie_linear_category, lajaunie_restricted_category, planar_constraint_matrix,
        single_surface_inequality_category, single_surface_linear_category,
        single_surface_restricted_category, stratigraphic_category,
        stratigraphic_restricted_category, tangent_constraint_matrix, vector_field_category,
    },
    constraints_to_points, fit_continuous_property, fit_lajaunie_linear, fit_lajaunie_restricted,
    fit_single_surface_inequality, fit_single_surface_linear, fit_single_surface_restricted,
    fit_stratigraphic, fit_stratigraphic_restricted, fit_vector_field, spatial_metrics, BuildError,
    ConstraintError, Constraints, ContinuousPropertyError, ContinuousPropertyModel, DenseMatrix,
    Error, GreedyTrace, LajaunieLinearError, LajaunieLinearModel, LajaunieRestrictedError,
    LajaunieRestrictedModel, ModelType, Parameters, Point, SingleSurfaceInequalityError,
    SingleSurfaceInequalityModel, SingleSurfaceLinearError, SingleSurfaceLinearModel,
    SingleSurfaceRestrictedError, SingleSurfaceRestrictedModel, SpatialError, SpatialParameters,
    StratigraphicError, StratigraphicModel, StratigraphicRestrictedError,
    StratigraphicRestrictedModel, VectorFieldError, VectorFieldModel,
};

/// Solver branch selected by the public fitted-model dispatcher.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum FitBranch {
    Linear,
    OrdinaryQuadratic,
    RestrictedRange,
}

/// A scalar/gradient input or model evaluation failure.
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub enum EvaluationError {
    Constraint(ConstraintError),
    IncorrectArrayDimensions {
        rows: usize,
        columns: usize,
        expected_columns: usize,
    },
    SingleSurfaceLinear(SingleSurfaceLinearError),
    SingleSurfaceInequality(SingleSurfaceInequalityError),
    SingleSurfaceRestricted(SingleSurfaceRestrictedError),
    LajaunieLinear(LajaunieLinearError),
    LajaunieRestricted(LajaunieRestrictedError),
    Stratigraphic(StratigraphicError),
    StratigraphicRestricted(StratigraphicRestrictedError),
    ContinuousProperty(ContinuousPropertyError),
    VectorField(VectorFieldError),
}

impl fmt::Display for EvaluationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Constraint(error) => error.fmt(formatter),
            Self::IncorrectArrayDimensions {
                rows,
                columns,
                expected_columns,
            } => write!(
                formatter,
                "evaluation array must be non-empty with {expected_columns} columns; got {rows}x{columns}"
            ),
            Self::SingleSurfaceLinear(error) => error.fmt(formatter),
            Self::SingleSurfaceInequality(error) => error.fmt(formatter),
            Self::SingleSurfaceRestricted(error) => error.fmt(formatter),
            Self::LajaunieLinear(error) => error.fmt(formatter),
            Self::LajaunieRestricted(error) => error.fmt(formatter),
            Self::Stratigraphic(error) => error.fmt(formatter),
            Self::StratigraphicRestricted(error) => error.fmt(formatter),
            Self::ContinuousProperty(error) => error.fmt(formatter),
            Self::VectorField(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for EvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Constraint(error) => Some(error),
            Self::SingleSurfaceLinear(error) => Some(error),
            Self::SingleSurfaceInequality(error) => Some(error),
            Self::SingleSurfaceRestricted(error) => Some(error),
            Self::LajaunieLinear(error) => Some(error),
            Self::LajaunieRestricted(error) => Some(error),
            Self::Stratigraphic(error) => Some(error),
            Self::StratigraphicRestricted(error) => Some(error),
            Self::ContinuousProperty(error) => Some(error),
            Self::VectorField(error) => Some(error),
            Self::IncorrectArrayDimensions { .. } => None,
        }
    }
}

impl EvaluationError {
    /// Return the frozen public exception category for this evaluation error.
    ///
    /// A fitted Rust model cannot express Surfe's `missinginterpolant` or
    /// `interpolantneedsupdate` states. Those categories remain available in
    /// [`Error`], while this method reports only failures reachable after a
    /// successful immutable fit.
    pub const fn surfe_category(&self) -> Option<Error> {
        match self {
            Self::Constraint(_) => None,
            Self::IncorrectArrayDimensions { .. } => Some(Error::IncorrectArrayDimensions),
            Self::SingleSurfaceLinear(error) => single_surface_linear_category(error),
            Self::SingleSurfaceInequality(error) => single_surface_inequality_category(error),
            Self::SingleSurfaceRestricted(error) => single_surface_restricted_category(error),
            Self::LajaunieLinear(error) => lajaunie_linear_category(error),
            Self::LajaunieRestricted(error) => lajaunie_restricted_category(error),
            Self::Stratigraphic(error) => stratigraphic_category(error),
            Self::StratigraphicRestricted(error) => stratigraphic_restricted_category(error),
            Self::ContinuousProperty(error) => continuous_property_category(error),
            Self::VectorField(error) => vector_field_category(error),
        }
    }
}

/// Immutable fitted state for one of the five Surfe model kinds.
///
/// The variants expose completed model-specific evidence without permitting
/// mutation of matrices, weights, kernels, or source constraints. Restricted
/// and Stratigraphic public evaluation follows frozen `Surfe_API`: it evaluates
/// the fitted Modified-Kernel field, because frozen `ComputeInterpolant` never
/// calls the separately migrated explicit reconstruction operation.
#[derive(Debug)]
#[non_exhaustive]
pub enum FittedModel {
    SingleSurfaceLinear(SingleSurfaceLinearModel),
    SingleSurfaceInequality(SingleSurfaceInequalityModel),
    SingleSurfaceRestricted(SingleSurfaceRestrictedModel),
    LajaunieLinear(LajaunieLinearModel),
    LajaunieRestricted(LajaunieRestrictedModel),
    Stratigraphic(StratigraphicModel),
    StratigraphicRestricted(StratigraphicRestrictedModel),
    ContinuousProperty(ContinuousPropertyModel),
    VectorField(VectorFieldModel),
}

impl FittedModel {
    pub(crate) fn fit(
        constraints: &Constraints,
        parameters: &Parameters,
    ) -> Result<Self, BuildError> {
        match parameters.model_type {
            ModelType::SingleSurface if parameters.use_restricted_range => {
                fit_single_surface_restricted(constraints, parameters)
                    .map(Self::SingleSurfaceRestricted)
                    .map_err(BuildError::SingleSurfaceRestricted)
            }
            ModelType::SingleSurface if !constraints.inequalities.is_empty() => {
                fit_single_surface_inequality(constraints, parameters)
                    .map(Self::SingleSurfaceInequality)
                    .map_err(BuildError::SingleSurfaceInequality)
            }
            ModelType::SingleSurface => fit_single_surface_linear(constraints, parameters)
                .map(Self::SingleSurfaceLinear)
                .map_err(BuildError::SingleSurfaceLinear),
            ModelType::LajaunieApproach if parameters.use_restricted_range => {
                fit_lajaunie_restricted(constraints, parameters)
                    .map(Self::LajaunieRestricted)
                    .map_err(BuildError::LajaunieRestricted)
            }
            ModelType::LajaunieApproach => fit_lajaunie_linear(constraints, parameters)
                .map(Self::LajaunieLinear)
                .map_err(BuildError::LajaunieLinear),
            ModelType::StratigraphicHorizons if parameters.use_restricted_range => {
                fit_stratigraphic_restricted(constraints, parameters)
                    .map(Self::StratigraphicRestricted)
                    .map_err(BuildError::StratigraphicRestricted)
            }
            ModelType::StratigraphicHorizons => fit_stratigraphic(constraints, parameters)
                .map(Self::Stratigraphic)
                .map_err(BuildError::Stratigraphic),
            ModelType::ContinuousProperty => fit_continuous_property(constraints, parameters)
                .map(Self::ContinuousProperty)
                .map_err(BuildError::ContinuousProperty),
            ModelType::VectorField => fit_vector_field(constraints, parameters)
                .map(Self::VectorField)
                .map_err(BuildError::VectorField),
        }
    }

    pub const fn model_type(&self) -> ModelType {
        match self {
            Self::SingleSurfaceLinear(_)
            | Self::SingleSurfaceInequality(_)
            | Self::SingleSurfaceRestricted(_) => ModelType::SingleSurface,
            Self::LajaunieLinear(_) | Self::LajaunieRestricted(_) => ModelType::LajaunieApproach,
            Self::Stratigraphic(_) | Self::StratigraphicRestricted(_) => {
                ModelType::StratigraphicHorizons
            }
            Self::ContinuousProperty(_) => ModelType::ContinuousProperty,
            Self::VectorField(_) => ModelType::VectorField,
        }
    }

    pub const fn fit_branch(&self) -> FitBranch {
        match self {
            Self::SingleSurfaceLinear(_)
            | Self::LajaunieLinear(_)
            | Self::ContinuousProperty(_)
            | Self::VectorField(_) => FitBranch::Linear,
            Self::SingleSurfaceInequality(_) | Self::Stratigraphic(_) => {
                FitBranch::OrdinaryQuadratic
            }
            Self::SingleSurfaceRestricted(_)
            | Self::LajaunieRestricted(_)
            | Self::StratigraphicRestricted(_) => FitBranch::RestrictedRange,
        }
    }

    pub const fn parameters(&self) -> &Parameters {
        match self {
            Self::SingleSurfaceLinear(model) => model.parameters(),
            Self::SingleSurfaceInequality(model) => model.parameters(),
            Self::SingleSurfaceRestricted(model) => model.parameters(),
            Self::LajaunieLinear(model) => model.parameters(),
            Self::LajaunieRestricted(model) => model.parameters(),
            Self::Stratigraphic(model) => model.parameters(),
            Self::StratigraphicRestricted(model) => model.parameters(),
            Self::ContinuousProperty(model) => model.parameters(),
            Self::VectorField(model) => model.parameters(),
        }
    }

    pub const fn constraints(&self) -> &Constraints {
        match self {
            Self::SingleSurfaceLinear(model) => model.constraints(),
            Self::SingleSurfaceInequality(model) => model.constraints(),
            Self::SingleSurfaceRestricted(model) => model.constraints(),
            Self::LajaunieLinear(model) => model.constraints(),
            Self::LajaunieRestricted(model) => model.constraints(),
            Self::Stratigraphic(model) => model.constraints(),
            Self::StratigraphicRestricted(model) => model.constraints(),
            Self::ContinuousProperty(model) => model.constraints(),
            Self::VectorField(model) => model.constraints(),
        }
    }

    /// Return the deterministic Greedy evidence for this public fit.
    ///
    /// Frozen `Surfe_API::ComputeInterpolant` never calls the source-only
    /// Greedy loop, so the trace always has zero rounds even when the frozen
    /// setter stored `use_greedy = true`.
    pub fn greedy_trace(&self) -> GreedyTrace {
        GreedyTrace::public_fit(self.parameters())
    }

    pub fn interface_reference_points(&self) -> Vec<[f64; 3]> {
        let grouping = match self {
            Self::SingleSurfaceLinear(model) => Some(model.interface_grouping()),
            Self::SingleSurfaceInequality(model) => Some(model.interface_grouping()),
            Self::SingleSurfaceRestricted(model) => Some(model.interface_grouping()),
            Self::LajaunieLinear(model) => Some(model.interface_grouping()),
            Self::LajaunieRestricted(model) => Some(model.interface_grouping()),
            Self::Stratigraphic(model) => Some(model.interface_grouping()),
            Self::StratigraphicRestricted(model) => Some(model.interface_grouping()),
            Self::ContinuousProperty(_) | Self::VectorField(_) => None,
        };
        grouping
            .into_iter()
            .flat_map(|value| value.reference_indices())
            .map(|index| self.constraints().interfaces[*index].point().position())
            .collect()
    }

    pub fn number_of_interfaces(&self) -> usize {
        self.interface_reference_points().len()
    }

    pub fn interface_constraint_matrix(&self) -> DenseMatrix {
        interface_constraint_matrix(self.constraints())
    }

    pub fn inequality_constraint_matrix(&self) -> DenseMatrix {
        inequality_constraint_matrix(self.constraints())
    }

    pub fn planar_constraint_matrix(&self) -> DenseMatrix {
        planar_constraint_matrix(self.constraints())
    }

    pub fn tangent_constraint_matrix(&self) -> DenseMatrix {
        tangent_constraint_matrix(self.constraints())
    }

    pub fn data_bounds_and_resolution(&self) -> Result<SpatialParameters, SpatialError> {
        spatial_metrics(&constraints_to_points(self.constraints()))
    }

    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, EvaluationError> {
        match self {
            Self::SingleSurfaceLinear(model) => model
                .evaluate_scalar(point)
                .map_err(EvaluationError::SingleSurfaceLinear),
            Self::SingleSurfaceInequality(model) => model
                .evaluate_scalar(point)
                .map_err(EvaluationError::SingleSurfaceInequality),
            Self::SingleSurfaceRestricted(model) => model
                .evaluate_modified_scalar(point)
                .map_err(EvaluationError::SingleSurfaceRestricted),
            Self::LajaunieLinear(model) => model
                .evaluate_scalar(point)
                .map_err(EvaluationError::LajaunieLinear),
            Self::LajaunieRestricted(model) => model
                .evaluate_modified_scalar(point)
                .map_err(EvaluationError::LajaunieRestricted),
            Self::Stratigraphic(model) => model
                .evaluate_modified_scalar(point)
                .map_err(EvaluationError::Stratigraphic),
            Self::StratigraphicRestricted(model) => model
                .evaluate_modified_scalar(point)
                .map_err(EvaluationError::StratigraphicRestricted),
            Self::ContinuousProperty(model) => model
                .evaluate_scalar(point)
                .map_err(EvaluationError::ContinuousProperty),
            Self::VectorField(model) => model
                .evaluate_potential(point)
                .map_err(EvaluationError::VectorField),
        }
    }

    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], EvaluationError> {
        match self {
            Self::SingleSurfaceLinear(model) => model
                .evaluate_gradient(point)
                .map_err(EvaluationError::SingleSurfaceLinear),
            Self::SingleSurfaceInequality(model) => model
                .evaluate_gradient(point)
                .map_err(EvaluationError::SingleSurfaceInequality),
            Self::SingleSurfaceRestricted(model) => model
                .evaluate_modified_gradient(point)
                .map_err(EvaluationError::SingleSurfaceRestricted),
            Self::LajaunieLinear(model) => model
                .evaluate_gradient(point)
                .map_err(EvaluationError::LajaunieLinear),
            Self::LajaunieRestricted(model) => model
                .evaluate_modified_gradient(point)
                .map_err(EvaluationError::LajaunieRestricted),
            Self::Stratigraphic(model) => model
                .evaluate_modified_gradient(point)
                .map_err(EvaluationError::Stratigraphic),
            Self::StratigraphicRestricted(model) => model
                .evaluate_modified_gradient(point)
                .map_err(EvaluationError::StratigraphicRestricted),
            Self::ContinuousProperty(model) => model
                .evaluate_gradient(point)
                .map_err(EvaluationError::ContinuousProperty),
            Self::VectorField(model) => model
                .evaluate_gradient(point)
                .map_err(EvaluationError::VectorField),
        }
    }

    pub fn evaluate_scalars(&self, points: &[Point]) -> Result<Vec<f64>, EvaluationError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }

    pub fn evaluate_gradients(&self, points: &[Point]) -> Result<Vec<[f64; 3]>, EvaluationError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }

    pub fn evaluate_scalar_matrix(
        &self,
        locations: &DenseMatrix,
    ) -> Result<Vec<f64>, EvaluationError> {
        let points = points_from_matrix(locations)?;
        self.evaluate_scalars(&points)
    }

    pub fn evaluate_gradient_matrix(
        &self,
        locations: &DenseMatrix,
    ) -> Result<Vec<[f64; 3]>, EvaluationError> {
        let points = points_from_matrix(locations)?;
        self.evaluate_gradients(&points)
    }
}

fn points_from_matrix(locations: &DenseMatrix) -> Result<Vec<Point>, EvaluationError> {
    if locations.rows() == 0 || locations.cols() != 3 {
        return Err(EvaluationError::IncorrectArrayDimensions {
            rows: locations.rows(),
            columns: locations.cols(),
            expected_columns: 3,
        });
    }
    locations
        .data()
        .chunks_exact(3)
        .map(|row| Point::new(row[0], row[1], row[2]).map_err(EvaluationError::Constraint))
        .collect()
}
