//! Surfe-compatible discrete parameters and defaults.
//!
//! Sources:
//! - `surfe_lib/modelling_parameters.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`

use std::{fmt, str::FromStr};

use crate::Error;

/// Degrees-to-radians factor used by frozen Surfe (`D2R`).
pub const DEGREES_TO_RADIANS: f64 = 0.017_453_292_519_943_295;
/// Radians-to-degrees factor used by frozen Surfe (`R2D`).
pub const RADIANS_TO_DEGREES: f64 = 57.295_779_513_082_32;
/// Position-comparison threshold used by frozen Surfe (`Epilson`).
pub const POSITION_EPSILON: f64 = 1E-3;

/// Which kernel argument a first derivative acts on (`Parameter_Types::DWRT`).
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum DerivativePoint {
    First,
    Second,
}

/// Mixed second-derivative component in Surfe declaration order.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum SecondDerivative {
    DxDx,
    DxDy,
    DxDz,
    DyDx,
    DyDy,
    DyDz,
    DzDx,
    DzDy,
    DzDz,
}

/// First-derivative component in Surfe declaration order.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum FirstDerivative {
    Dx,
    Dy,
    Dz,
}

/// Radial basis function choice in frozen Surfe declaration order.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum RbfKernel {
    Cubic,
    Gaussian,
    Multiquadric,
    MultiquadricCubic,
    InverseMultiquadric,
    ThinPlateSpline,
    Linear,
    WendlandC2,
    MaternC4,
}

impl RbfKernel {
    /// Every frozen public kernel choice in C++ declaration order.
    ///
    /// This table lets compatibility tests prove that no public name branch is
    /// silently omitted when the enum grows or is refactored.
    pub const ALL: [Self; 9] = [
        Self::Cubic,
        Self::Gaussian,
        Self::Multiquadric,
        Self::MultiquadricCubic,
        Self::InverseMultiquadric,
        Self::ThinPlateSpline,
        Self::Linear,
        Self::WendlandC2,
        Self::MaternC4,
    ];

    /// The exact spelling accepted by `Surfe_API::SetRBFKernel(const char*)`.
    pub const fn surfe_name(self) -> &'static str {
        match self {
            Self::Cubic => "r3",
            Self::Gaussian => "Gaussian",
            Self::Multiquadric => "Multiquadratics",
            Self::MultiquadricCubic => "Multiquadratics3",
            Self::InverseMultiquadric => "Inverse Multiquadratics",
            Self::ThinPlateSpline => "Thin Plate Spline",
            Self::Linear => "r",
            Self::WendlandC2 => "WendlandC2",
            Self::MaternC4 => "MaternC4",
        }
    }
}

impl FromStr for RbfKernel {
    type Err = Error;

    fn from_str(name: &str) -> Result<Self, Self::Err> {
        match name {
            "r3" => Ok(Self::Cubic),
            "WendlandC2" => Ok(Self::WendlandC2),
            "r" => Ok(Self::Linear),
            "Gaussian" => Ok(Self::Gaussian),
            "Multiquadratics" => Ok(Self::Multiquadric),
            "Multiquadratics3" => Ok(Self::MultiquadricCubic),
            "Thin Plate Spline" => Ok(Self::ThinPlateSpline),
            "Inverse Multiquadratics" => Ok(Self::InverseMultiquadric),
            "MaternC4" => Ok(Self::MaternC4),
            _ => Err(Error::UnknownRbf),
        }
    }
}

impl fmt::Display for RbfKernel {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.surfe_name())
    }
}

/// Solver branch selected by model setup.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum SolverType {
    Linear,
    Quadratic,
}

/// Geological model kind in frozen `ModelType` declaration order.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum ModelType {
    SingleSurface,
    LajaunieApproach,
    StratigraphicHorizons,
    ContinuousProperty,
    VectorField,
}

impl ModelType {
    /// Every frozen model kind in C++ `ModelType` declaration order.
    pub const ALL: [Self; 5] = [
        Self::SingleSurface,
        Self::LajaunieApproach,
        Self::StratigraphicHorizons,
        Self::ContinuousProperty,
        Self::VectorField,
    ];

    /// The frozen C++ enum identifier, used as GeoRBF's exact text form.
    pub const fn surfe_enum_name(self) -> &'static str {
        match self {
            Self::SingleSurface => "Single_surface",
            Self::LajaunieApproach => "Lajaunie_approach",
            Self::StratigraphicHorizons => "Stratigraphic_horizons",
            Self::ContinuousProperty => "Continuous_property",
            Self::VectorField => "Vector_field",
        }
    }

    /// Integer accepted by frozen `Surfe_API(int)` for this model.
    ///
    /// This deliberately differs from the C++ enum discriminant order for
    /// Vector Field, Stratigraphic Horizons, and Continuous Property.
    pub const fn surfe_api_code(self) -> i32 {
        match self {
            Self::SingleSurface => 1,
            Self::LajaunieApproach => 2,
            Self::VectorField => 3,
            Self::StratigraphicHorizons => 4,
            Self::ContinuousProperty => 5,
        }
    }
}

impl TryFrom<i32> for ModelType {
    type Error = Error;

    fn try_from(code: i32) -> Result<Self, Self::Error> {
        match code {
            1 => Ok(Self::SingleSurface),
            2 => Ok(Self::LajaunieApproach),
            3 => Ok(Self::VectorField),
            4 => Ok(Self::StratigraphicHorizons),
            5 => Ok(Self::ContinuousProperty),
            _ => Err(Error::UnknownModel),
        }
    }
}

impl FromStr for ModelType {
    type Err = Error;

    fn from_str(name: &str) -> Result<Self, Self::Err> {
        match name {
            "Single_surface" => Ok(Self::SingleSurface),
            "Lajaunie_approach" => Ok(Self::LajaunieApproach),
            "Stratigraphic_horizons" => Ok(Self::StratigraphicHorizons),
            "Continuous_property" => Ok(Self::ContinuousProperty),
            "Vector_field" => Ok(Self::VectorField),
            _ => Err(Error::UnknownModel),
        }
    }
}

impl fmt::Display for ModelType {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.surfe_enum_name())
    }
}

/// Cartesian axis in frozen Surfe declaration order.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(u8)]
pub enum Axis {
    X,
    Y,
    Z,
}

/// User-visible modelling parameters with frozen Surfe defaults.
#[derive(Clone, Debug, PartialEq)]
pub struct Parameters {
    pub model_type: ModelType,
    pub min_stratigraphic_thickness: f64,
    pub use_interface: bool,
    pub use_planar: bool,
    pub use_tangent: bool,
    pub use_inequality: bool,
    pub basis_type: RbfKernel,
    pub shape_parameter: f64,
    pub polynomial_order: i32,
    pub advanced_parameters: bool,
    pub model_global_anisotropy: bool,
    pub use_greedy: bool,
    pub use_restricted_range: bool,
    pub smoothing_amount: f64,
    pub use_regression_smoothing: bool,
    pub interface_uncertainty: f64,
    pub angular_uncertainty: f64,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            model_type: ModelType::SingleSurface,
            min_stratigraphic_thickness: 0.0,
            use_interface: false,
            use_planar: false,
            use_tangent: false,
            use_inequality: false,
            basis_type: RbfKernel::Cubic,
            shape_parameter: 100.0,
            polynomial_order: 1,
            advanced_parameters: false,
            model_global_anisotropy: false,
            use_greedy: false,
            use_restricted_range: false,
            smoothing_amount: 0.0,
            use_regression_smoothing: false,
            interface_uncertainty: 0.0,
            angular_uncertainty: 0.0,
        }
    }
}

impl Parameters {
    pub fn set_rbf_kernel(&mut self, kernel: RbfKernel) {
        self.basis_type = kernel;
    }

    pub fn set_rbf_kernel_name(&mut self, name: &str) -> Result<(), Error> {
        let kernel = name.parse()?;
        self.set_rbf_kernel(kernel);
        Ok(())
    }

    pub fn set_rbf_shape_parameter(&mut self, shape_parameter: f64) {
        self.shape_parameter = shape_parameter;
    }

    pub fn set_polynomial_order(&mut self, polynomial_order: i32) {
        self.polynomial_order = polynomial_order;
    }

    pub fn set_global_anisotropy(&mut self, enabled: bool) {
        self.model_global_anisotropy = enabled;
    }

    pub fn set_restricted_range(
        &mut self,
        enabled: bool,
        interface_uncertainty: f64,
        angular_uncertainty: f64,
    ) {
        self.use_restricted_range = enabled;
        self.interface_uncertainty = interface_uncertainty;
        self.angular_uncertainty = angular_uncertainty;
    }

    /// Match frozen `SetRegressionSmoothing`, which ignores its boolean input.
    pub fn set_regression_smoothing(&mut self, _enabled: bool, amount: f64) {
        self.use_regression_smoothing = true;
        self.smoothing_amount = amount;
    }

    /// Match frozen `SetGreedyAlgorithm`, which ignores its boolean input.
    pub fn set_greedy_algorithm(
        &mut self,
        _enabled: bool,
        interface_uncertainty: f64,
        angular_uncertainty: f64,
    ) {
        self.use_greedy = true;
        self.interface_uncertainty = interface_uncertainty;
        self.angular_uncertainty = angular_uncertainty;
    }
}

/// Derived model parameters initialized before layout and assembly.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InternalParameters {
    pub n_interface: usize,
    pub n_planar: usize,
    pub n_inequality: usize,
    pub n_tangent: usize,
    pub n_constraints: usize,
    pub n_equality: usize,
    pub modified_basis: bool,
    pub poly_term: bool,
    pub n_poly_terms: usize,
    pub problem_type: SolverType,
    pub restricted_range: bool,
}

impl Default for InternalParameters {
    fn default() -> Self {
        Self {
            n_interface: 0,
            n_planar: 0,
            n_inequality: 0,
            n_tangent: 0,
            n_constraints: 0,
            n_equality: 0,
            modified_basis: false,
            poly_term: true,
            n_poly_terms: 4,
            problem_type: SolverType::Linear,
            restricted_range: false,
        }
    }
}

/// Frozen file-oriented parameter aggregate, represented with owned Rust data.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InputParameters {
    pub parameters: Parameters,
    pub interface_file: String,
    pub planar_file: String,
    pub tangent_file: String,
    pub inequality_file: String,
}
