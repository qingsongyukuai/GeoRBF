use std::str::FromStr;

use georbf::{
    Builder, Constraints, Error, Inequality, Interface, ModelType, Parameters, Planar, Point,
    RbfKernel,
};

const SOURCE_COMMIT: &str = "290dbe0ab344f4258a4935f05cad0f153f0f69a4";

fn single_surface_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 0.0).unwrap(),
            Interface::new(1.0, 0.0, 0.0, 0.0).unwrap(),
            Interface::new(0.0, 1.0, 0.0, 0.0).unwrap(),
            Interface::new(1.0, 1.0, 0.0, 0.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(0.5, 0.5, 0.0, 0.0, 0.0, 1.0).unwrap()],
        ..Constraints::default()
    }
}

fn lajaunie_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 30.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        ..Constraints::default()
    }
}

fn continuous_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 1.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, -0.5).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 2.25).unwrap(),
            Interface::new(0.0, 0.0, 4.0, -1.75).unwrap(),
            Interface::new(1.5, -2.0, 2.5, 0.75).unwrap(),
        ],
        ..Constraints::default()
    }
}

fn assert_oracle_close(actual: f64, expected_bits: u64, absolute: f64, relative: f64) {
    let expected = f64::from_bits(expected_bits);
    let tolerance = absolute.max(relative * expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn stratigraphic_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 30.0).unwrap(),
            Interface::new(5.0, 1.0, 1.0, 20.0).unwrap(),
            Interface::new(7.0, -1.0, 2.0, 10.0).unwrap(),
        ],
        inequalities: vec![
            Inequality::new(3.0, -2.0, 2.0, 25.0).unwrap(),
            Inequality::new(4.0, 2.0, -1.0, 5.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        ..Constraints::default()
    }
}

fn vector_constraints() -> Constraints {
    Constraints {
        planars: vec![
            Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap(),
            Planar::from_normal(-2.0, 1.0, 0.5, 0.0, 0.8, 0.6).unwrap(),
            Planar::from_normal(0.5, -1.0, 2.0, 0.36, -0.48, 0.8).unwrap(),
        ],
        ..Constraints::default()
    }
}

fn constant_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 3.0).unwrap(),
            Interface::new(1.0, 0.0, 0.0, 3.0).unwrap(),
            Interface::new(0.0, 1.0, 0.0, 3.0).unwrap(),
            Interface::new(0.0, 0.0, 1.0, 3.0).unwrap(),
        ],
        ..Constraints::default()
    }
}

fn quadratic_constraints() -> Constraints {
    let value = |x: f64, y: f64, z: f64| {
        x * x
            + 2.0 * y * y
            + 3.0 * z * z
            + 4.0 * x * y
            + 5.0 * x * z
            + 6.0 * y * z
            + 7.0 * x
            + 8.0 * y
            + 9.0 * z
            + 10.0
    };
    Constraints {
        interfaces: [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
        .into_iter()
        .map(|[x, y, z]| Interface::new(x, y, z, value(x, y, z)).unwrap())
        .collect(),
        ..Constraints::default()
    }
}

#[test]
fn public_kernel_and_model_compatibility_tables_are_exhaustive_and_exact() {
    assert_eq!(SOURCE_COMMIT.len(), 40);

    let kernel_names = [
        "r3",
        "Gaussian",
        "Multiquadratics",
        "Multiquadratics3",
        "Inverse Multiquadratics",
        "Thin Plate Spline",
        "r",
        "WendlandC2",
        "MaternC4",
    ];
    assert_eq!(RbfKernel::ALL.len(), kernel_names.len());
    for (kernel, name) in RbfKernel::ALL.into_iter().zip(kernel_names) {
        assert_eq!(kernel.surfe_name(), name);
        assert_eq!(RbfKernel::from_str(name), Ok(kernel));
        let mut builder = Builder::new(ModelType::SingleSurface);
        builder.set_rbf_kernel_name(name).unwrap();
        assert_eq!(builder.parameters().basis_type, kernel);
    }

    for rejected in [
        "",
        "Cubic",
        "R3",
        "gaussian",
        "Gaussian ",
        "Wendland C2",
        "MQ",
        "TPS",
        "Matern C4",
    ] {
        let mut builder = Builder::new(ModelType::SingleSurface);
        let before = builder.parameters().basis_type;
        assert!(builder.set_rbf_kernel_name(rejected).is_err());
        assert_eq!(builder.parameters().basis_type, before);
    }

    assert_eq!(ModelType::ALL.len(), 5);
    for model in ModelType::ALL {
        assert_eq!(ModelType::from_str(model.surfe_enum_name()), Ok(model));
        assert_eq!(ModelType::try_from(model.surfe_api_code()), Ok(model));
        assert_eq!(
            Builder::from_surfe_model_code(model.surfe_api_code())
                .unwrap()
                .parameters()
                .model_type,
            model
        );
    }
}

#[test]
fn every_public_model_has_a_default_configuration_snapshot() {
    let defaults = Parameters::default();
    assert_eq!(defaults.model_type, ModelType::SingleSurface);
    assert_eq!(defaults.min_stratigraphic_thickness, 0.0);
    assert_eq!(defaults.basis_type, RbfKernel::Cubic);
    assert_eq!(defaults.shape_parameter, 100.0);
    assert_eq!(defaults.polynomial_order, 1);
    assert!(!defaults.advanced_parameters);
    assert!(!defaults.model_global_anisotropy);
    assert!(!defaults.use_greedy);
    assert!(!defaults.use_restricted_range);
    assert_eq!(defaults.smoothing_amount, 0.0);
    assert!(!defaults.use_regression_smoothing);
    assert_eq!(defaults.interface_uncertainty, 0.0);
    assert_eq!(defaults.angular_uncertainty, 0.0);
    assert!(!defaults.use_interface);
    assert!(!defaults.use_planar);
    assert!(!defaults.use_tangent);
    assert!(!defaults.use_inequality);

    for model in ModelType::ALL {
        let builder = Builder::new(model);
        let actual = builder.parameters();
        assert_eq!(actual.model_type, model);
        assert_eq!(
            actual.min_stratigraphic_thickness,
            defaults.min_stratigraphic_thickness
        );
        assert_eq!(actual.basis_type, defaults.basis_type);
        assert_eq!(actual.shape_parameter, defaults.shape_parameter);
        assert_eq!(actual.polynomial_order, defaults.polynomial_order);
        assert_eq!(actual.advanced_parameters, defaults.advanced_parameters);
        assert_eq!(
            actual.model_global_anisotropy,
            defaults.model_global_anisotropy
        );
        assert_eq!(actual.use_greedy, defaults.use_greedy);
        assert_eq!(actual.use_restricted_range, defaults.use_restricted_range);
        assert_eq!(actual.smoothing_amount, defaults.smoothing_amount);
        assert_eq!(
            actual.use_regression_smoothing,
            defaults.use_regression_smoothing
        );
        assert_eq!(actual.interface_uncertainty, defaults.interface_uncertainty);
        assert_eq!(actual.angular_uncertainty, defaults.angular_uncertainty);
        assert!(!actual.use_interface);
        assert!(!actual.use_planar);
        assert!(!actual.use_tangent);
        assert!(!actual.use_inequality);
    }
}

#[test]
fn default_public_builds_match_the_frozen_oracle() {
    let cases = [
        (
            ModelType::SingleSurface,
            0x3fe0_0000_0000_0000,
            [
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x3ff0_0000_0000_0000,
            ],
        ),
        (
            ModelType::LajaunieApproach,
            0x3fec_10f9_d2e2_9c4a,
            [
                0x3fa6_5f27_d256_bc1c,
                0xbfaa_f183_7d94_5da0,
                0xbfcf_9639_7e7e_aa36,
            ],
        ),
        (
            ModelType::VectorField,
            0xbfc2_6789_da3f_9380,
            [
                0x3fcd_628d_517a_9eaf,
                0x3fca_ca92_78ec_a8f8,
                0x3fe6_0baa_7247_f0f4,
            ],
        ),
        (
            ModelType::StratigraphicHorizons,
            0xbfc2_25a5_e3c2_1329,
            [
                0x3fba_f566_40fd_4ebd,
                0xbfa9_238f_043f_e56e,
                0xbfd0_454c_51d3_9c86,
            ],
        ),
        (
            ModelType::ContinuousProperty,
            0x3ff2_1391_a86d_8a8d,
            [
                0xbfe7_f9ec_d70b_dd0f,
                0x3fd9_8d71_533f_4e3e,
                0x3fba_5e42_1840_a4e0,
            ],
        ),
    ];
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();

    for (model_type, scalar, gradient) in cases {
        let constraints = match model_type {
            ModelType::SingleSurface => single_surface_constraints(),
            ModelType::LajaunieApproach => lajaunie_constraints(),
            ModelType::VectorField => vector_constraints(),
            ModelType::StratigraphicHorizons => stratigraphic_constraints(),
            ModelType::ContinuousProperty => continuous_constraints(),
        };
        let fitted = Builder::new(model_type)
            .set_constraints(constraints)
            .fit()
            .unwrap();
        assert_oracle_close(
            fitted.evaluate_scalar(&witness).unwrap(),
            scalar,
            1.0e-9,
            1.0e-8,
        );
        for (actual, expected) in fitted
            .evaluate_gradient(&witness)
            .unwrap()
            .into_iter()
            .zip(gradient)
        {
            assert_oracle_close(actual, expected, 1.0e-9, 1.0e-8);
        }
    }
}

#[test]
fn all_public_kernel_names_match_the_frozen_oracle() {
    let expected_scalars = [
        0x3ff2_1391_a86d_8a8d,
        0x3f9e_ec0a_0e17_f202,
        0x3fe6_e775_983f_ce58,
        0x3ff1_d7d4_61e5_1c37,
        0x3fed_1d74_b8dc_8b20,
        0x3fe7_445b_91e8_a342,
        0x3fe7_8d83_58f3_09c4,
        0x3fcd_7ff4_6cdb_6e64,
        0x3fe8_1666_65d7_05b4,
    ];
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();

    for (kernel, expected) in RbfKernel::ALL.into_iter().zip(expected_scalars) {
        let mut builder = Builder::new(ModelType::ContinuousProperty);
        builder
            .set_constraints(continuous_constraints())
            .set_rbf_kernel(kernel)
            .set_rbf_shape_parameter(2.0);
        let actual = builder.fit().unwrap().evaluate_scalar(&witness).unwrap();
        assert_oracle_close(actual, expected, 2.0e-11, 2.0e-10);
    }
}

#[test]
fn legal_public_configuration_matrix_reaches_all_model_and_kernel_branches() {
    for kernel in RbfKernel::ALL {
        let mut builder = Builder::new(ModelType::ContinuousProperty);
        builder.set_constraints(continuous_constraints());
        builder.set_rbf_kernel(kernel).set_rbf_shape_parameter(2.0);
        let model = builder.fit().unwrap();
        assert_eq!(model.model_type(), ModelType::ContinuousProperty);
    }

    for model in ModelType::ALL {
        let constraints = match model {
            ModelType::SingleSurface => single_surface_constraints(),
            ModelType::LajaunieApproach => lajaunie_constraints(),
            ModelType::StratigraphicHorizons => stratigraphic_constraints(),
            ModelType::ContinuousProperty => continuous_constraints(),
            ModelType::VectorField => vector_constraints(),
        };
        let fitted = Builder::new(model)
            .set_constraints(constraints)
            .fit()
            .unwrap();
        assert_eq!(fitted.model_type(), model);
    }

    for (order, constraints) in [
        (0, constant_constraints()),
        (1, single_surface_constraints()),
        (2, quadratic_constraints()),
    ] {
        let mut builder = Builder::new(ModelType::SingleSurface);
        builder
            .set_constraints(constraints)
            .set_polynomial_order(order);
        assert_eq!(builder.fit().unwrap().parameters().polynomial_order, order);
    }
}

#[test]
fn anisotropic_kernel_acceptance_matrix_matches_the_frozen_factory() {
    let supported = [
        RbfKernel::Cubic,
        RbfKernel::Gaussian,
        RbfKernel::Multiquadric,
        RbfKernel::InverseMultiquadric,
        RbfKernel::ThinPlateSpline,
        RbfKernel::Linear,
    ];
    let unsupported = [
        RbfKernel::MultiquadricCubic,
        RbfKernel::WendlandC2,
        RbfKernel::MaternC4,
    ];
    for kernel in supported {
        let mut builder = Builder::new(ModelType::SingleSurface);
        let mut constraints = single_surface_constraints();
        constraints
            .planars
            .push(Planar::from_normal(0.25, 0.75, 0.0, 0.0, 0.0, 1.0).unwrap());
        builder
            .set_constraints(constraints)
            .set_rbf_kernel(kernel)
            .set_rbf_shape_parameter(2.0)
            .set_global_anisotropy(true);
        if let Err(error) = builder.fit() {
            assert_ne!(
                error.surfe_category(),
                Some(Error::BasisFunctionSetupFailure),
                "{kernel} should pass frozen anisotropic factory construction"
            );
        }
    }
    for kernel in unsupported {
        let mut builder = Builder::new(ModelType::SingleSurface);
        let mut constraints = single_surface_constraints();
        constraints
            .planars
            .push(Planar::from_normal(0.25, 0.75, 0.0, 0.0, 0.0, 1.0).unwrap());
        builder
            .set_constraints(constraints)
            .set_rbf_kernel(kernel)
            .set_global_anisotropy(true);
        assert_eq!(
            builder.fit().unwrap_err().surfe_category(),
            Some(Error::BasisFunctionSetupFailure),
            "{kernel} should be rejected by the frozen anisotropic factory"
        );
    }
}
