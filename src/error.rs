//! Stable error categories corresponding to Surfe exceptions.
//!
//! Source: `surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`.

use std::fmt;

/// A stable GeoRBF category for one frozen Surfe exception type.
///
/// The category, rather than the historical English message, is intended for
/// programmatic matching. [`Error::message`] and
/// [`Error::surfe_exception_name`] retain source compatibility evidence.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[non_exhaustive]
pub enum Error {
    NoInterfaceData,
    NoInterfaceIncrementPairs,
    NoPlanarData,
    InvalidInputData,
    GlobalAnisotropyFailure,
    AnisotropicKernelCreationFailure,
    BasisFunctionSetupFailure,
    ModifiedKernelCreationFailure,
    LagrangianBasisCreationFailure,
    LinearSolverFailure,
    PredictorCorrectorSolverFailure,
    LoqoSolverFailure,
    InterpolationMatrixFailure,
    EqualityVectorFailure,
    InequalityVectorFailure,
    InterfaceIsoValueUpdateFailure,
    InterpolantComputationFailure,
    MissingInterpolant,
    UnknownRbf,
    InterpolantNeedsUpdate,
    UnknownModel,
    SpatialParametersFailure,
    IncorrectArrayDimensions,
}

impl Error {
    /// Every exception category declared by frozen Surfe, in declaration order.
    pub const ALL: [Self; 23] = [
        Self::NoInterfaceData,
        Self::NoInterfaceIncrementPairs,
        Self::NoPlanarData,
        Self::InvalidInputData,
        Self::GlobalAnisotropyFailure,
        Self::AnisotropicKernelCreationFailure,
        Self::BasisFunctionSetupFailure,
        Self::ModifiedKernelCreationFailure,
        Self::LagrangianBasisCreationFailure,
        Self::LinearSolverFailure,
        Self::PredictorCorrectorSolverFailure,
        Self::LoqoSolverFailure,
        Self::InterpolationMatrixFailure,
        Self::EqualityVectorFailure,
        Self::InequalityVectorFailure,
        Self::InterfaceIsoValueUpdateFailure,
        Self::InterpolantComputationFailure,
        Self::MissingInterpolant,
        Self::UnknownRbf,
        Self::InterpolantNeedsUpdate,
        Self::UnknownModel,
        Self::SpatialParametersFailure,
        Self::IncorrectArrayDimensions,
    ];

    /// The exact C++ exception class mapped to this category.
    pub const fn surfe_exception_name(self) -> &'static str {
        match self {
            Self::NoInterfaceData => "nointerfacedata",
            Self::NoInterfaceIncrementPairs => "nointerfaceincrementpairs",
            Self::NoPlanarData => "noplanardata",
            Self::InvalidInputData => "invalidinputdata",
            Self::GlobalAnisotropyFailure => "failurecomputingglobalanisotropy",
            Self::AnisotropicKernelCreationFailure => "failurecreatinganisotropickernel",
            Self::BasisFunctionSetupFailure => "failuresettingupbasisfunctions",
            Self::ModifiedKernelCreationFailure => "failurecreatingmodifiedkernel",
            Self::LagrangianBasisCreationFailure => "failurecreatinglagrangianpolynomialbasis",
            Self::LinearSolverFailure => "linearsolverfailure",
            Self::PredictorCorrectorSolverFailure => "pcquadratricsolverfailure",
            Self::LoqoSolverFailure => "loqoquadratricsolverfailure",
            Self::InterpolationMatrixFailure => "errorcomputinginterpolationmatrix",
            Self::EqualityVectorFailure => "errorcomputingequalityvector",
            Self::InequalityVectorFailure => "errorcomputinginequalityvector",
            Self::InterfaceIsoValueUpdateFailure => "errorupdatinginterfaceisovalues",
            Self::InterpolantComputationFailure => "errorcomputinginterpolant",
            Self::MissingInterpolant => "missinginterpolant",
            Self::UnknownRbf => "unknownrbf",
            Self::InterpolantNeedsUpdate => "interpolantneedsupdate",
            Self::UnknownModel => "unknownmodellingmode",
            Self::SpatialParametersFailure => "problemcomputingspatialparameters",
            Self::IncorrectArrayDimensions => "arrayhasincorrectdimensions",
        }
    }

    /// The exact `what()` text from the mapped frozen Surfe exception.
    pub const fn message(self) -> &'static str {
        match self {
            Self::NoInterfaceData => "No interface data",
            Self::NoInterfaceIncrementPairs => "There are no interface increment pairs",
            Self::NoPlanarData => "No planar data",
            Self::InvalidInputData => {
                "Invalid input data as determined by check_input_data()"
            }
            Self::GlobalAnisotropyFailure => {
                "Failure computing global anisotropy because there are less than 2 planar constraints"
            }
            Self::AnisotropicKernelCreationFailure => {
                "Failure creating an anisotropic kernel"
            }
            Self::BasisFunctionSetupFailure => "Failure setting up basis functions",
            Self::ModifiedKernelCreationFailure => "Failure creating modified kernel",
            Self::LagrangianBasisCreationFailure => {
                "Failure creating Lagrangian Polynomial basis"
            }
            Self::LinearSolverFailure => "Eigen's linear solver failed",
            Self::PredictorCorrectorSolverFailure => {
                "Predictor-Corrector Quadratic Solver failure"
            }
            Self::LoqoSolverFailure => "LOQO Quadratic Solver failure",
            Self::InterpolationMatrixFailure => "Error computing interpolation matrix",
            Self::EqualityVectorFailure => "Error computing equality vector",
            Self::InequalityVectorFailure => "Error computing inequality vector",
            Self::InterfaceIsoValueUpdateFailure => "Error updating interface iso values",
            Self::InterpolantComputationFailure => "Error computing Interpolant",
            Self::MissingInterpolant => "Interpolant has not yet been computed",
            Self::UnknownRbf => "Entered RBF kernel name is unknown",
            Self::InterpolantNeedsUpdate => {
                "Constraints or Parameters have changed please recompute/update interpolant"
            }
            Self::UnknownModel => "Modelling mode code unknown; choose 1 - 5",
            Self::SpatialParametersFailure => "Problem computing spatial parameters",
            Self::IncorrectArrayDimensions => "Input array has incorrect dimensions!",
        }
    }

    /// Reproduce the message format of Surfe's `SurfeExceptions` wrapper.
    ///
    /// This is diagnostic compatibility only; callers should match the typed
    /// [`Error`] values instead of parsing the returned text.
    pub fn format_surfe_exception_chain(chain: &[Self]) -> String {
        let messages = chain
            .iter()
            .map(|error| error.message())
            .collect::<Vec<_>>()
            .join(", ");
        format!("Exceptions thrown: {messages}")
    }
}

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.message())
    }
}

impl std::error::Error for Error {}
