use std::str::FromStr;

use georbf::{
    Axis, DerivativePoint, Error, FirstDerivative, InputParameters, InternalParameters, ModelType,
    Parameters, RbfKernel, SecondDerivative, SolverType, DEGREES_TO_RADIANS, POSITION_EPSILON,
    RADIANS_TO_DEGREES,
};

#[test]
fn parameter_enum_discriminants_match_surfe_declaration_order() {
    assert_eq!(DerivativePoint::First as u8, 0);
    assert_eq!(DerivativePoint::Second as u8, 1);

    assert_eq!(FirstDerivative::Dx as u8, 0);
    assert_eq!(FirstDerivative::Dy as u8, 1);
    assert_eq!(FirstDerivative::Dz as u8, 2);

    let second_derivatives = [
        SecondDerivative::DxDx,
        SecondDerivative::DxDy,
        SecondDerivative::DxDz,
        SecondDerivative::DyDx,
        SecondDerivative::DyDy,
        SecondDerivative::DyDz,
        SecondDerivative::DzDx,
        SecondDerivative::DzDy,
        SecondDerivative::DzDz,
    ];
    for (index, derivative) in second_derivatives.into_iter().enumerate() {
        assert_eq!(derivative as usize, index);
    }

    let kernels = [
        RbfKernel::Cubic,
        RbfKernel::Gaussian,
        RbfKernel::Multiquadric,
        RbfKernel::MultiquadricCubic,
        RbfKernel::InverseMultiquadric,
        RbfKernel::ThinPlateSpline,
        RbfKernel::Linear,
        RbfKernel::WendlandC2,
        RbfKernel::MaternC4,
    ];
    for (index, kernel) in kernels.into_iter().enumerate() {
        assert_eq!(kernel as usize, index);
    }

    assert_eq!(SolverType::Linear as u8, 0);
    assert_eq!(SolverType::Quadratic as u8, 1);

    let models = [
        ModelType::SingleSurface,
        ModelType::LajaunieApproach,
        ModelType::StratigraphicHorizons,
        ModelType::ContinuousProperty,
        ModelType::VectorField,
    ];
    for (index, model) in models.into_iter().enumerate() {
        assert_eq!(model as usize, index);
    }

    assert_eq!(Axis::X as u8, 0);
    assert_eq!(Axis::Y as u8, 1);
    assert_eq!(Axis::Z as u8, 2);
}

#[test]
fn parameter_defaults_match_frozen_surfe() {
    let parameters = Parameters::default();
    assert_eq!(parameters.model_type, ModelType::SingleSurface);
    assert_eq!(parameters.min_stratigraphic_thickness, 0.0);
    assert!(!parameters.use_interface);
    assert!(!parameters.use_planar);
    assert!(!parameters.use_tangent);
    assert!(!parameters.use_inequality);
    assert_eq!(parameters.basis_type, RbfKernel::Cubic);
    assert_eq!(parameters.shape_parameter, 100.0);
    assert_eq!(parameters.polynomial_order, 1);
    assert!(!parameters.advanced_parameters);
    assert!(!parameters.model_global_anisotropy);
    assert!(!parameters.use_greedy);
    assert!(!parameters.use_restricted_range);
    assert_eq!(parameters.smoothing_amount, 0.0);
    assert!(!parameters.use_regression_smoothing);
    assert_eq!(parameters.interface_uncertainty, 0.0);
    assert_eq!(parameters.angular_uncertainty, 0.0);

    let internal = InternalParameters::default();
    assert_eq!(internal.n_interface, 0);
    assert_eq!(internal.n_planar, 0);
    assert_eq!(internal.n_inequality, 0);
    assert_eq!(internal.n_tangent, 0);
    assert_eq!(internal.n_constraints, 0);
    assert_eq!(internal.n_equality, 0);
    assert!(!internal.modified_basis);
    assert!(internal.poly_term);
    assert_eq!(internal.n_poly_terms, 4);
    assert_eq!(internal.problem_type, SolverType::Linear);
    assert!(!internal.restricted_range);

    let input = InputParameters::default();
    assert_eq!(input.parameters, parameters);
    assert!(input.interface_file.is_empty());
    assert!(input.planar_file.is_empty());
    assert!(input.tangent_file.is_empty());
    assert!(input.inequality_file.is_empty());
}

#[test]
fn constants_match_frozen_surfe_binary64_values() {
    assert_eq!(DEGREES_TO_RADIANS.to_bits(), 0x3f91_df46_a252_9d39);
    assert_eq!(RADIANS_TO_DEGREES.to_bits(), 0x404c_a5dc_1a63_c1f8);
    assert_eq!(POSITION_EPSILON.to_bits(), 0x3f50_624d_d2f1_a9fc);
}

#[test]
fn rbf_names_are_exact_and_round_trip() {
    let accepted = [
        ("r3", RbfKernel::Cubic),
        ("WendlandC2", RbfKernel::WendlandC2),
        ("r", RbfKernel::Linear),
        ("Gaussian", RbfKernel::Gaussian),
        ("Multiquadratics", RbfKernel::Multiquadric),
        ("Multiquadratics3", RbfKernel::MultiquadricCubic),
        ("Thin Plate Spline", RbfKernel::ThinPlateSpline),
        ("Inverse Multiquadratics", RbfKernel::InverseMultiquadric),
        ("MaternC4", RbfKernel::MaternC4),
    ];
    for (name, expected) in accepted {
        assert_eq!(RbfKernel::from_str(name), Ok(expected));
        assert_eq!(expected.to_string(), name);
    }

    for rejected in [
        "",
        "Cubic",
        "cubic",
        "R3",
        "wendlandc2",
        "Gaussian ",
        " Gaussian",
        "MQ",
        "TPS",
        "InverseMultiquadric",
        "Matern C4",
    ] {
        assert_eq!(RbfKernel::from_str(rejected), Err(Error::UnknownRbf));
    }
}

#[test]
fn model_enum_names_and_public_integer_codes_are_distinct_exact_mappings() {
    let names = [
        ("Single_surface", ModelType::SingleSurface),
        ("Lajaunie_approach", ModelType::LajaunieApproach),
        ("Stratigraphic_horizons", ModelType::StratigraphicHorizons),
        ("Continuous_property", ModelType::ContinuousProperty),
        ("Vector_field", ModelType::VectorField),
    ];
    for (name, expected) in names {
        assert_eq!(ModelType::from_str(name), Ok(expected));
        assert_eq!(expected.to_string(), name);
    }
    for rejected in [
        "",
        "single_surface",
        "Single Surface",
        "VectorField",
        "vector_field",
        "1",
    ] {
        assert_eq!(ModelType::from_str(rejected), Err(Error::UnknownModel));
    }

    let legacy_codes = [
        (1, ModelType::SingleSurface),
        (2, ModelType::LajaunieApproach),
        (3, ModelType::VectorField),
        (4, ModelType::StratigraphicHorizons),
        (5, ModelType::ContinuousProperty),
    ];
    for (code, expected) in legacy_codes {
        assert_eq!(ModelType::try_from(code), Ok(expected));
        assert_eq!(expected.surfe_api_code(), code);
    }
    for rejected in [i32::MIN, -1, 0, 6, i32::MAX] {
        assert_eq!(ModelType::try_from(rejected), Err(Error::UnknownModel));
    }
}

#[test]
fn parameter_setters_preserve_surfe_mutation_semantics() {
    let mut parameters = Parameters::default();
    parameters.set_rbf_kernel_name("MaternC4").unwrap();
    parameters.set_rbf_shape_parameter(7.5);
    parameters.set_polynomial_order(2);
    parameters.set_global_anisotropy(true);
    parameters.set_restricted_range(true, 0.25, 4.0);
    assert_eq!(parameters.basis_type, RbfKernel::MaternC4);
    assert_eq!(parameters.shape_parameter, 7.5);
    assert_eq!(parameters.polynomial_order, 2);
    assert!(parameters.model_global_anisotropy);
    assert!(parameters.use_restricted_range);
    assert_eq!(parameters.interface_uncertainty, 0.25);
    assert_eq!(parameters.angular_uncertainty, 4.0);

    let before = parameters.clone();
    assert_eq!(
        parameters.set_rbf_kernel_name("maternc4"),
        Err(Error::UnknownRbf)
    );
    assert_eq!(parameters, before);

    parameters.set_regression_smoothing(false, 3.0);
    assert!(parameters.use_regression_smoothing);
    assert_eq!(parameters.smoothing_amount, 3.0);
    parameters.set_greedy_algorithm(false, 0.5, 6.0);
    assert!(parameters.use_greedy);
    assert_eq!(parameters.interface_uncertainty, 0.5);
    assert_eq!(parameters.angular_uncertainty, 6.0);
}

#[test]
fn every_surfe_exception_has_one_stable_typed_mapping() {
    let expected = [
        (
            Error::NoInterfaceData,
            "nointerfacedata",
            "No interface data",
        ),
        (
            Error::NoInterfaceIncrementPairs,
            "nointerfaceincrementpairs",
            "There are no interface increment pairs",
        ),
        (Error::NoPlanarData, "noplanardata", "No planar data"),
        (
            Error::InvalidInputData,
            "invalidinputdata",
            "Invalid input data as determined by check_input_data()",
        ),
        (
            Error::GlobalAnisotropyFailure,
            "failurecomputingglobalanisotropy",
            "Failure computing global anisotropy because there are less than 2 planar constraints",
        ),
        (
            Error::AnisotropicKernelCreationFailure,
            "failurecreatinganisotropickernel",
            "Failure creating an anisotropic kernel",
        ),
        (
            Error::BasisFunctionSetupFailure,
            "failuresettingupbasisfunctions",
            "Failure setting up basis functions",
        ),
        (
            Error::ModifiedKernelCreationFailure,
            "failurecreatingmodifiedkernel",
            "Failure creating modified kernel",
        ),
        (
            Error::LagrangianBasisCreationFailure,
            "failurecreatinglagrangianpolynomialbasis",
            "Failure creating Lagrangian Polynomial basis",
        ),
        (
            Error::LinearSolverFailure,
            "linearsolverfailure",
            "Eigen's linear solver failed",
        ),
        (
            Error::PredictorCorrectorSolverFailure,
            "pcquadratricsolverfailure",
            "Predictor-Corrector Quadratic Solver failure",
        ),
        (
            Error::LoqoSolverFailure,
            "loqoquadratricsolverfailure",
            "LOQO Quadratic Solver failure",
        ),
        (
            Error::InterpolationMatrixFailure,
            "errorcomputinginterpolationmatrix",
            "Error computing interpolation matrix",
        ),
        (
            Error::EqualityVectorFailure,
            "errorcomputingequalityvector",
            "Error computing equality vector",
        ),
        (
            Error::InequalityVectorFailure,
            "errorcomputinginequalityvector",
            "Error computing inequality vector",
        ),
        (
            Error::InterfaceIsoValueUpdateFailure,
            "errorupdatinginterfaceisovalues",
            "Error updating interface iso values",
        ),
        (
            Error::InterpolantComputationFailure,
            "errorcomputinginterpolant",
            "Error computing Interpolant",
        ),
        (
            Error::MissingInterpolant,
            "missinginterpolant",
            "Interpolant has not yet been computed",
        ),
        (
            Error::UnknownRbf,
            "unknownrbf",
            "Entered RBF kernel name is unknown",
        ),
        (
            Error::InterpolantNeedsUpdate,
            "interpolantneedsupdate",
            "Constraints or Parameters have changed please recompute/update interpolant",
        ),
        (
            Error::UnknownModel,
            "unknownmodellingmode",
            "Modelling mode code unknown; choose 1 - 5",
        ),
        (
            Error::SpatialParametersFailure,
            "problemcomputingspatialparameters",
            "Problem computing spatial parameters",
        ),
        (
            Error::IncorrectArrayDimensions,
            "arrayhasincorrectdimensions",
            "Input array has incorrect dimensions!",
        ),
    ];

    assert_eq!(Error::ALL.len(), expected.len());
    for (index, (error, exception_name, message)) in expected.into_iter().enumerate() {
        assert_eq!(Error::ALL[index], error);
        assert_eq!(error.surfe_exception_name(), exception_name);
        assert_eq!(error.message(), message);
        assert_eq!(error.to_string(), message);
    }
}

#[test]
fn surfe_exception_wrapper_format_is_reproducible_without_using_text_as_category() {
    let chain = Error::format_surfe_exception_chain(&[
        Error::BasisFunctionSetupFailure,
        Error::AnisotropicKernelCreationFailure,
    ]);
    assert_eq!(
        chain,
        "Exceptions thrown: Failure setting up basis functions, Failure creating an anisotropic kernel"
    );
    assert_eq!(
        Error::format_surfe_exception_chain(&[]),
        "Exceptions thrown: "
    );
}
