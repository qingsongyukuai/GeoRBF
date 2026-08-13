//! Owned, validated configuration for the public GeoRBF fitting lifecycle.
//!
//! The compatibility surface follows
//! `surfe_lib/surfe_api.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`.
//! Unlike the mutable C++ API, fitting snapshots this builder into an immutable
//! [`crate::FittedModel`]. Later builder changes cannot stale an existing fit.

use std::fmt;

use crate::{
    constraints_to_points, spatial_metrics, ConstraintError, Constraints, ContinuousPropertyError,
    DenseMatrix, Error, FittedModel, Inequality, Interface, LajaunieLinearError,
    LajaunieRestrictedError, ModelType, Parameters, Planar, Polarity, RbfKernel,
    SingleSurfaceInequalityError, SingleSurfaceLinearError, SingleSurfaceRestrictedError,
    SpatialError, SpatialParameters, StratigraphicError, StratigraphicRestrictedError, Tangent,
    VectorFieldError,
};

/// A configuration, input, or model-fit failure from [`Builder`].
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub enum BuildError {
    Constraint(ConstraintError),
    Surfe(Error),
    IncorrectArrayDimensions {
        rows: usize,
        columns: usize,
        expected_columns: usize,
    },
    NonFiniteParameter(&'static str),
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

impl fmt::Display for BuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Constraint(error) => error.fmt(formatter),
            Self::Surfe(error) => error.fmt(formatter),
            Self::IncorrectArrayDimensions {
                rows,
                columns,
                expected_columns,
            } => write!(
                formatter,
                "constraint array must be non-empty with {expected_columns} columns; got {rows}x{columns}"
            ),
            Self::NonFiniteParameter(name) => {
                write!(formatter, "parameter {name} must be finite")
            }
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

impl std::error::Error for BuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Constraint(error) => Some(error),
            Self::Surfe(error) => Some(error),
            Self::SingleSurfaceLinear(error) => Some(error),
            Self::SingleSurfaceInequality(error) => Some(error),
            Self::SingleSurfaceRestricted(error) => Some(error),
            Self::LajaunieLinear(error) => Some(error),
            Self::LajaunieRestricted(error) => Some(error),
            Self::Stratigraphic(error) => Some(error),
            Self::StratigraphicRestricted(error) => Some(error),
            Self::ContinuousProperty(error) => Some(error),
            Self::VectorField(error) => Some(error),
            Self::IncorrectArrayDimensions { .. } | Self::NonFiniteParameter(_) => None,
        }
    }
}

impl From<ConstraintError> for BuildError {
    fn from(error: ConstraintError) -> Self {
        Self::Constraint(error)
    }
}

impl From<Error> for BuildError {
    fn from(error: Error) -> Self {
        Self::Surfe(error)
    }
}

impl BuildError {
    /// Return the frozen Surfe exception category visible at the public API.
    ///
    /// `None` marks a Rust safety rejection for source behavior that had no
    /// stable C++ exception (for example non-finite input accepted into invalid
    /// state, or the Continuous Property equality-vector out-of-bounds bug).
    /// Callers should match this value rather than parsing `to_string()`.
    pub const fn surfe_category(&self) -> Option<Error> {
        match self {
            Self::Constraint(_) | Self::NonFiniteParameter(_) => None,
            Self::Surfe(error) => Some(*error),
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

pub(crate) const fn assembly_category(error: &crate::AssemblyError) -> Option<Error> {
    match error {
        crate::AssemblyError::Surfe(error) => Some(*error),
        crate::AssemblyError::Kernel(_)
        | crate::AssemblyError::Constraint(_)
        | crate::AssemblyError::KernelLayoutMismatch => None,
    }
}

pub(crate) const fn reconstruction_assembly_category(
    error: &crate::ReconstructionAssemblyError,
) -> Option<Error> {
    match error {
        crate::ReconstructionAssemblyError::Assembly(error) => assembly_category(error),
        crate::ReconstructionAssemblyError::UnsupportedModel
        | crate::ReconstructionAssemblyError::NotQuadratic
        | crate::ReconstructionAssemblyError::SourceKernelNotModified
        | crate::ReconstructionAssemblyError::OrdinaryKernelIsModified
        | crate::ReconstructionAssemblyError::SourceLayoutMismatch
        | crate::ReconstructionAssemblyError::SourceWeightLengthMismatch
        | crate::ReconstructionAssemblyError::NonFinitePrediction => None,
    }
}

const fn reconstruction_category(error: &crate::ReconstructionError) -> Option<Error> {
    match error {
        crate::ReconstructionError::SourceAssembly(error) => assembly_category(error),
        crate::ReconstructionError::PredictorCorrector(error) => Some(error.surfe_error()),
        crate::ReconstructionError::Loqo(error) => Some(error.surfe_error()),
        crate::ReconstructionError::Reassembly(error) => reconstruction_assembly_category(error),
        crate::ReconstructionError::Lu(error) => Some(error.surfe_error()),
    }
}

pub(crate) const fn single_surface_linear_category(
    error: &SingleSurfaceLinearError,
) -> Option<Error> {
    match error {
        SingleSurfaceLinearError::WrongModel => Some(Error::UnknownModel),
        SingleSurfaceLinearError::InequalityBranchNotAvailable
        | SingleSurfaceLinearError::RestrictedRangeBranchNotAvailable
        | SingleSurfaceLinearError::Evaluation(_) => None,
        SingleSurfaceLinearError::Surfe(error) => Some(*error),
        SingleSurfaceLinearError::Anisotropy(_) => Some(Error::BasisFunctionSetupFailure),
        SingleSurfaceLinearError::Assembly(error) => assembly_category(error),
        SingleSurfaceLinearError::Lu(error) => Some(error.surfe_error()),
    }
}

pub(crate) const fn single_surface_inequality_category(
    error: &SingleSurfaceInequalityError,
) -> Option<Error> {
    match error {
        SingleSurfaceInequalityError::WrongModel => Some(Error::UnknownModel),
        SingleSurfaceInequalityError::NoInequalities
        | SingleSurfaceInequalityError::RestrictedRangeBranchNotAvailable
        | SingleSurfaceInequalityError::Evaluation(_) => None,
        SingleSurfaceInequalityError::Surfe(error) => Some(*error),
        SingleSurfaceInequalityError::Anisotropy(_) | SingleSurfaceInequalityError::Basis(_) => {
            Some(Error::BasisFunctionSetupFailure)
        }
        SingleSurfaceInequalityError::Assembly(error) => assembly_category(error),
        SingleSurfaceInequalityError::Qp(error) => Some(error.surfe_error()),
    }
}

pub(crate) const fn single_surface_restricted_category(
    error: &SingleSurfaceRestrictedError,
) -> Option<Error> {
    match error {
        SingleSurfaceRestrictedError::WrongModel => Some(Error::UnknownModel),
        SingleSurfaceRestrictedError::RestrictedRangeRequired => None,
        SingleSurfaceRestrictedError::Surfe(error) => Some(*error),
        SingleSurfaceRestrictedError::Anisotropy(_) | SingleSurfaceRestrictedError::Basis(_) => {
            Some(Error::BasisFunctionSetupFailure)
        }
        SingleSurfaceRestrictedError::SourceAssembly(error) => assembly_category(error),
        SingleSurfaceRestrictedError::Loqo(error) => Some(error.surfe_error()),
        SingleSurfaceRestrictedError::Reconstruction(error) => reconstruction_category(error),
        SingleSurfaceRestrictedError::Evaluation(error) => reconstruction_assembly_category(error),
    }
}

pub(crate) const fn lajaunie_linear_category(error: &LajaunieLinearError) -> Option<Error> {
    match error {
        LajaunieLinearError::WrongModel => Some(Error::UnknownModel),
        LajaunieLinearError::RestrictedRangeBranchNotAvailable
        | LajaunieLinearError::Evaluation(_) => None,
        LajaunieLinearError::Surfe(error) => Some(*error),
        LajaunieLinearError::Anisotropy(_) => Some(Error::BasisFunctionSetupFailure),
        LajaunieLinearError::Assembly(error) => assembly_category(error),
        LajaunieLinearError::Lu(error) => Some(error.surfe_error()),
    }
}

pub(crate) const fn lajaunie_restricted_category(error: &LajaunieRestrictedError) -> Option<Error> {
    match error {
        LajaunieRestrictedError::WrongModel => Some(Error::UnknownModel),
        LajaunieRestrictedError::RestrictedRangeRequired => None,
        LajaunieRestrictedError::Surfe(error) => Some(*error),
        LajaunieRestrictedError::Anisotropy(_) | LajaunieRestrictedError::Basis(_) => {
            Some(Error::BasisFunctionSetupFailure)
        }
        LajaunieRestrictedError::SourceAssembly(error) => assembly_category(error),
        // Lajaunie's frozen restricted branch constructs LOQO but maps
        // `solve() == false` to `pcquadratricsolverfailure`.
        LajaunieRestrictedError::Loqo(_) => Some(Error::PredictorCorrectorSolverFailure),
        LajaunieRestrictedError::Reconstruction(error) => reconstruction_category(error),
        LajaunieRestrictedError::Evaluation(error) => reconstruction_assembly_category(error),
    }
}

pub(crate) const fn stratigraphic_category(error: &StratigraphicError) -> Option<Error> {
    match error {
        StratigraphicError::WrongModel => Some(Error::UnknownModel),
        StratigraphicError::RestrictedRangeBranchNotAvailable
        | StratigraphicError::Evaluation(_) => None,
        StratigraphicError::Surfe(error) => Some(*error),
        StratigraphicError::Anisotropy(_) | StratigraphicError::Basis(_) => {
            Some(Error::BasisFunctionSetupFailure)
        }
        StratigraphicError::SourceAssembly(error) => assembly_category(error),
        StratigraphicError::Qp(error) => Some(error.surfe_error()),
        StratigraphicError::Reconstruction(error) => reconstruction_category(error),
    }
}

pub(crate) const fn stratigraphic_restricted_category(
    error: &StratigraphicRestrictedError,
) -> Option<Error> {
    match error {
        StratigraphicRestrictedError::WrongModel => Some(Error::UnknownModel),
        StratigraphicRestrictedError::RestrictedRangeRequired => None,
        StratigraphicRestrictedError::Surfe(error) => Some(*error),
        StratigraphicRestrictedError::Anisotropy(_) | StratigraphicRestrictedError::Basis(_) => {
            Some(Error::BasisFunctionSetupFailure)
        }
        StratigraphicRestrictedError::SourceAssembly(error) => assembly_category(error),
        StratigraphicRestrictedError::Loqo(error) => Some(error.surfe_error()),
        StratigraphicRestrictedError::Reconstruction(error) => reconstruction_category(error),
        StratigraphicRestrictedError::Evaluation(error) => reconstruction_assembly_category(error),
    }
}

pub(crate) const fn continuous_property_category(error: &ContinuousPropertyError) -> Option<Error> {
    match error {
        ContinuousPropertyError::WrongModel => Some(Error::UnknownModel),
        ContinuousPropertyError::EqualityVectorOutOfBounds { .. }
        | ContinuousPropertyError::Evaluation(_) => None,
        ContinuousPropertyError::Surfe(error) => Some(*error),
        ContinuousPropertyError::Anisotropy(_) => Some(Error::BasisFunctionSetupFailure),
        ContinuousPropertyError::Assembly(error) => assembly_category(error),
        ContinuousPropertyError::Lu(error) => Some(error.surfe_error()),
    }
}

pub(crate) const fn vector_field_category(error: &VectorFieldError) -> Option<Error> {
    match error {
        VectorFieldError::WrongModel => Some(Error::UnknownModel),
        VectorFieldError::Anisotropy(_) => Some(Error::BasisFunctionSetupFailure),
        VectorFieldError::Assembly(error) => assembly_category(error),
        VectorFieldError::Lu(error) => Some(error.surfe_error()),
        VectorFieldError::Evaluation(_) => None,
    }
}

/// Mutable, owned input configuration that produces immutable fitted models.
///
/// Constraint values and parameters are owned by the builder. [`Self::fit`]
/// snapshots them into a new [`FittedModel`], so the builder may be edited and
/// fitted again without changing any earlier model.
///
/// ```
/// use georbf::{Builder, ModelType, Point, RbfKernel};
///
/// let mut builder = Builder::new(ModelType::SingleSurface);
/// builder.set_rbf_kernel(RbfKernel::Cubic).set_polynomial_order(1);
/// for [x, y, z] in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]] {
///     builder.add_interface_xyz(x, y, z, 0.0)?;
/// }
/// builder.add_planar_normal(0.5, 0.5, 0.0, 0.0, 0.0, 1.0)?;
/// let fitted = builder.fit()?;
/// let value = fitted.evaluate_scalar(&Point::new(0.25, 0.75, 0.5)?)?;
/// assert_eq!(value, 0.5);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// Evaluation is intentionally unavailable before fitting:
///
/// ```compile_fail
/// use georbf::{Builder, ModelType, Point};
///
/// let builder = Builder::new(ModelType::SingleSurface);
/// let point = Point::new(0.0, 0.0, 0.0)?;
/// let _ = builder.evaluate_scalar(&point);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Builder {
    parameters: Parameters,
    constraints: Constraints,
}

impl Builder {
    pub fn new(model_type: ModelType) -> Self {
        let parameters = Parameters {
            model_type,
            ..Parameters::default()
        };
        Self {
            parameters,
            constraints: Constraints::default(),
        }
    }

    /// Construct from frozen `Surfe_API(int)` model codes 1 through 5.
    pub fn from_surfe_model_code(code: i32) -> Result<Self, BuildError> {
        Ok(Self::new(ModelType::try_from(code)?))
    }

    pub fn from_parameters(parameters: Parameters) -> Self {
        Self {
            parameters,
            constraints: Constraints::default(),
        }
    }

    pub const fn parameters(&self) -> &Parameters {
        &self.parameters
    }

    pub const fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub fn set_constraints(&mut self, constraints: Constraints) -> &mut Self {
        self.constraints = constraints;
        self.sync_constraint_flags();
        self
    }

    pub fn set_model_type(&mut self, model_type: ModelType) -> &mut Self {
        self.parameters.model_type = model_type;
        self
    }

    pub fn set_rbf_kernel(&mut self, kernel: RbfKernel) -> &mut Self {
        self.parameters.set_rbf_kernel(kernel);
        self
    }

    pub fn set_rbf_kernel_name(&mut self, name: &str) -> Result<&mut Self, BuildError> {
        self.parameters.set_rbf_kernel_name(name)?;
        Ok(self)
    }

    pub fn set_rbf_shape_parameter(&mut self, shape_parameter: f64) -> &mut Self {
        self.parameters.set_rbf_shape_parameter(shape_parameter);
        self
    }

    pub fn set_polynomial_order(&mut self, polynomial_order: i32) -> &mut Self {
        self.parameters.set_polynomial_order(polynomial_order);
        self
    }

    pub fn set_global_anisotropy(&mut self, enabled: bool) -> &mut Self {
        self.parameters.set_global_anisotropy(enabled);
        self
    }

    pub fn set_min_stratigraphic_thickness(&mut self, thickness: f64) -> &mut Self {
        self.parameters.min_stratigraphic_thickness = thickness;
        self
    }

    pub fn set_regression_smoothing(&mut self, enabled: bool, amount: f64) -> &mut Self {
        self.parameters.set_regression_smoothing(enabled, amount);
        self
    }

    /// Retain frozen Surfe's flag-writing setter.
    ///
    /// The frozen public fit has no call edge to its source-only Greedy loop.
    /// [`FittedModel::greedy_trace`](crate::FittedModel::greedy_trace) records
    /// that compatible zero-round behavior.
    pub fn set_greedy_algorithm(
        &mut self,
        enabled: bool,
        interface_uncertainty: f64,
        angular_uncertainty: f64,
    ) -> &mut Self {
        self.parameters
            .set_greedy_algorithm(enabled, interface_uncertainty, angular_uncertainty);
        self
    }

    pub fn set_restricted_range(
        &mut self,
        enabled: bool,
        interface_uncertainty: f64,
        angular_uncertainty: f64,
    ) -> &mut Self {
        self.parameters
            .set_restricted_range(enabled, interface_uncertainty, angular_uncertainty);
        self
    }

    pub fn add_interface(&mut self, constraint: Interface) -> &mut Self {
        self.constraints.interfaces.push(constraint);
        self.parameters.use_interface = true;
        self
    }

    pub fn add_interface_xyz(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        level: f64,
    ) -> Result<&mut Self, BuildError> {
        let constraint = Interface::new(x, y, z, level)?;
        Ok(self.add_interface(constraint))
    }

    pub fn add_inequality(&mut self, constraint: Inequality) -> &mut Self {
        self.constraints.inequalities.push(constraint);
        self.parameters.use_inequality = true;
        self
    }

    pub fn add_inequality_xyz(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        level: f64,
    ) -> Result<&mut Self, BuildError> {
        let constraint = Inequality::new(x, y, z, level)?;
        Ok(self.add_inequality(constraint))
    }

    pub fn add_planar(&mut self, constraint: Planar) -> &mut Self {
        self.constraints.planars.push(constraint);
        self.parameters.use_planar = true;
        self
    }

    #[allow(clippy::too_many_arguments)]
    pub fn add_planar_normal(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        nx: f64,
        ny: f64,
        nz: f64,
    ) -> Result<&mut Self, BuildError> {
        let constraint = Planar::from_normal(x, y, z, nx, ny, nz)?;
        Ok(self.add_planar(constraint))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn add_planar_strike_dip_polarity(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        strike: f64,
        dip: f64,
        polarity: Polarity,
    ) -> Result<&mut Self, BuildError> {
        let constraint = Planar::from_strike_dip_polarity(x, y, z, strike, dip, polarity)?;
        Ok(self.add_planar(constraint))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn add_planar_azimuth_dip_polarity(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        azimuth: f64,
        dip: f64,
        polarity: Polarity,
    ) -> Result<&mut Self, BuildError> {
        let constraint = Planar::from_azimuth_dip_polarity(x, y, z, azimuth, dip, polarity)?;
        Ok(self.add_planar(constraint))
    }

    pub fn add_tangent(&mut self, constraint: Tangent) -> &mut Self {
        self.constraints.tangents.push(constraint);
        self.parameters.use_tangent = true;
        self
    }

    #[allow(clippy::too_many_arguments)]
    pub fn add_tangent_xyz(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        tx: f64,
        ty: f64,
        tz: f64,
    ) -> Result<&mut Self, BuildError> {
        let constraint = Tangent::new(x, y, z, tx, ty, tz)?;
        Ok(self.add_tangent(constraint))
    }

    pub fn replace_interfaces(&mut self, values: Vec<Interface>) -> &mut Self {
        self.constraints.interfaces = values;
        self.parameters.use_interface = !self.constraints.interfaces.is_empty();
        self
    }

    pub fn replace_inequalities(&mut self, values: Vec<Inequality>) -> &mut Self {
        self.constraints.inequalities = values;
        self.parameters.use_inequality = !self.constraints.inequalities.is_empty();
        self
    }

    pub fn replace_planars(&mut self, values: Vec<Planar>) -> &mut Self {
        self.constraints.planars = values;
        self.parameters.use_planar = !self.constraints.planars.is_empty();
        self
    }

    pub fn replace_tangents(&mut self, values: Vec<Tangent>) -> &mut Self {
        self.constraints.tangents = values;
        self.parameters.use_tangent = !self.constraints.tangents.is_empty();
        self
    }

    pub fn set_interface_constraint_matrix(
        &mut self,
        matrix: &DenseMatrix,
    ) -> Result<&mut Self, BuildError> {
        validate_table(matrix, 4)?;
        let values = matrix
            .data()
            .chunks_exact(4)
            .map(|row| Interface::new(row[0], row[1], row[2], row[3]))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self.replace_interfaces(values))
    }

    pub fn set_inequality_constraint_matrix(
        &mut self,
        matrix: &DenseMatrix,
    ) -> Result<&mut Self, BuildError> {
        validate_table(matrix, 4)?;
        let values = matrix
            .data()
            .chunks_exact(4)
            .map(|row| Inequality::new(row[0], row[1], row[2], row[3]))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self.replace_inequalities(values))
    }

    pub fn set_planar_constraint_matrix(
        &mut self,
        matrix: &DenseMatrix,
    ) -> Result<&mut Self, BuildError> {
        validate_table(matrix, 6)?;
        let values = matrix
            .data()
            .chunks_exact(6)
            .map(|row| Planar::from_normal(row[0], row[1], row[2], row[3], row[4], row[5]))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self.replace_planars(values))
    }

    /// Replace tangents from Surfe's `n x 6` compatibility layout.
    ///
    /// Frozen `SetTangentConstraints` accidentally routed rows to the planar
    /// adder. GeoRBF uses the documented tangent meaning and performs an
    /// atomic validated replacement instead of copying that source defect.
    pub fn set_tangent_constraint_matrix(
        &mut self,
        matrix: &DenseMatrix,
    ) -> Result<&mut Self, BuildError> {
        validate_table(matrix, 6)?;
        let values = matrix
            .data()
            .chunks_exact(6)
            .map(|row| Tangent::new(row[0], row[1], row[2], row[3], row[4], row[5]))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self.replace_tangents(values))
    }

    pub fn interface_constraint_matrix(&self) -> DenseMatrix {
        interface_constraint_matrix(&self.constraints)
    }

    pub fn inequality_constraint_matrix(&self) -> DenseMatrix {
        inequality_constraint_matrix(&self.constraints)
    }

    pub fn planar_constraint_matrix(&self) -> DenseMatrix {
        planar_constraint_matrix(&self.constraints)
    }

    pub fn tangent_constraint_matrix(&self) -> DenseMatrix {
        tangent_constraint_matrix(&self.constraints)
    }

    pub fn data_bounds_and_resolution(&self) -> Result<SpatialParameters, SpatialError> {
        spatial_metrics(&constraints_to_points(&self.constraints))
    }

    /// Fit a fully owned, immutable snapshot of this configuration.
    pub fn fit(&self) -> Result<FittedModel, BuildError> {
        validate_parameters(&self.parameters)?;
        FittedModel::fit(&self.constraints, &self.parameters)
    }

    fn sync_constraint_flags(&mut self) {
        self.parameters.use_interface = !self.constraints.interfaces.is_empty();
        self.parameters.use_inequality = !self.constraints.inequalities.is_empty();
        self.parameters.use_planar = !self.constraints.planars.is_empty();
        self.parameters.use_tangent = !self.constraints.tangents.is_empty();
    }
}

fn validate_table(matrix: &DenseMatrix, expected_columns: usize) -> Result<(), BuildError> {
    if matrix.rows() == 0 || matrix.cols() != expected_columns {
        Err(BuildError::IncorrectArrayDimensions {
            rows: matrix.rows(),
            columns: matrix.cols(),
            expected_columns,
        })
    } else {
        Ok(())
    }
}

fn validate_parameters(parameters: &Parameters) -> Result<(), BuildError> {
    for (name, value) in [
        (
            "min_stratigraphic_thickness",
            parameters.min_stratigraphic_thickness,
        ),
        ("shape_parameter", parameters.shape_parameter),
        ("smoothing_amount", parameters.smoothing_amount),
        ("interface_uncertainty", parameters.interface_uncertainty),
        ("angular_uncertainty", parameters.angular_uncertainty),
    ] {
        if !value.is_finite() {
            return Err(BuildError::NonFiniteParameter(name));
        }
    }
    Ok(())
}

pub(crate) fn interface_constraint_matrix(constraints: &Constraints) -> DenseMatrix {
    let data = constraints
        .interfaces
        .iter()
        .flat_map(|value| {
            let point = value.point();
            [point.x(), point.y(), point.z(), value.level()]
        })
        .collect();
    DenseMatrix::from_row_major(constraints.interfaces.len(), 4, data)
        .expect("four values are emitted for every interface")
}

pub(crate) fn inequality_constraint_matrix(constraints: &Constraints) -> DenseMatrix {
    let data = constraints
        .inequalities
        .iter()
        .flat_map(|value| {
            let point = value.point();
            [point.x(), point.y(), point.z(), value.level()]
        })
        .collect();
    DenseMatrix::from_row_major(constraints.inequalities.len(), 4, data)
        .expect("four values are emitted for every inequality")
}

pub(crate) fn planar_constraint_matrix(constraints: &Constraints) -> DenseMatrix {
    let data = constraints
        .planars
        .iter()
        .flat_map(|value| {
            let point = value.point();
            [
                point.x(),
                point.y(),
                point.z(),
                value.nx(),
                value.ny(),
                value.nz(),
            ]
        })
        .collect();
    DenseMatrix::from_row_major(constraints.planars.len(), 6, data)
        .expect("six values are emitted for every planar")
}

pub(crate) fn tangent_constraint_matrix(constraints: &Constraints) -> DenseMatrix {
    let data = constraints
        .tangents
        .iter()
        .flat_map(|value| {
            let point = value.point();
            [
                point.x(),
                point.y(),
                point.z(),
                value.tx(),
                value.ty(),
                value.tz(),
            ]
        })
        .collect();
    DenseMatrix::from_row_major(constraints.tangents.len(), 6, data)
        .expect("six values are emitted for every tangent")
}
