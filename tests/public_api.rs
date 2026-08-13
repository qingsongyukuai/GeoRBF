use std::{sync::Arc, thread};

use georbf::{
    fit_continuous_property, fit_lajaunie_linear, fit_stratigraphic, fit_vector_field, BuildError,
    Builder, Constraints, DenseMatrix, EvaluationError, FitBranch, FittedModel, Inequality,
    Interface, ModelType, Parameters, Planar, Point, RbfKernel, Tangent,
};

fn single_surface_builder() -> Builder {
    let mut builder = Builder::new(ModelType::SingleSurface);
    builder
        .set_rbf_kernel(RbfKernel::Cubic)
        .set_rbf_shape_parameter(2.0)
        .set_polynomial_order(1);
    for [x, y, z] in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ] {
        builder.add_interface_xyz(x, y, z, 0.0).unwrap();
    }
    builder
        .add_planar_normal(0.5, 0.5, 0.0, 0.0, 0.0, 1.0)
        .unwrap();
    builder
}

fn lajaunie_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 30.0).unwrap(),
            Interface::new(5.0, 1.0, 1.0, 20.0).unwrap(),
            Interface::new(6.0, 2.0, 1.0, 20.0).unwrap(),
            Interface::new(7.0, -1.0, 2.0, 10.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        tangents: vec![Tangent::new(-1.0, 2.0, 1.0, 0.2, 0.7, 0.1).unwrap()],
        ..Constraints::default()
    }
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

fn parameters(model_type: ModelType) -> Parameters {
    Parameters {
        model_type,
        basis_type: RbfKernel::Cubic,
        shape_parameter: 2.0,
        polynomial_order: 1,
        ..Parameters::default()
    }
}

fn assert_close(actual: f64, expected_bits: u64, tolerance: f64) {
    let expected = f64::from_bits(expected_bits);
    let scale = 1.0_f64.max(actual.abs()).max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn assert_gradient_close(actual: [f64; 3], expected: [u64; 3]) {
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert_close(actual, expected, 4.0e-11);
    }
}

fn assert_batch_matches_single(model: &FittedModel, points: &[Point]) {
    let scalars = model.evaluate_scalars(points).unwrap();
    let gradients = model.evaluate_gradients(points).unwrap();
    for index in 0..points.len() {
        assert_eq!(
            scalars[index],
            model.evaluate_scalar(&points[index]).unwrap()
        );
        assert_eq!(
            gradients[index],
            model.evaluate_gradient(&points[index]).unwrap()
        );
    }

    let locations = DenseMatrix::from_row_major(
        points.len(),
        3,
        points.iter().flat_map(|point| point.position()).collect(),
    )
    .unwrap();
    assert_eq!(model.evaluate_scalar_matrix(&locations).unwrap(), scalars);
    assert_eq!(
        model.evaluate_gradient_matrix(&locations).unwrap(),
        gradients
    );
}

#[test]
fn public_single_surface_matches_frozen_smoke_and_batch_order() {
    let mut builder = single_surface_builder();
    assert_eq!(builder.constraints().interfaces.len(), 4);
    assert_eq!(builder.interface_constraint_matrix().dimensions(), (4, 4));
    assert_eq!(builder.planar_constraint_matrix().dimensions(), (1, 6));

    let model = builder.fit().unwrap();
    builder.add_interface_xyz(9.0, 9.0, 9.0, 0.0).unwrap();
    assert_eq!(builder.constraints().interfaces.len(), 5);
    assert_eq!(model.constraints().interfaces.len(), 4);
    assert_eq!(model.model_type(), ModelType::SingleSurface);
    assert_eq!(model.fit_branch(), FitBranch::Linear);
    assert_eq!(model.number_of_interfaces(), 1);
    assert_eq!(model.interface_reference_points(), vec![[0.0, 0.0, 0.0]]);

    let points = [
        Point::new(0.25, 0.75, 0.5).unwrap(),
        Point::new(0.75, 0.25, -0.5).unwrap(),
    ];
    assert_eq!(model.evaluate_scalar(&points[0]).unwrap(), 0.5);
    assert_eq!(
        model.evaluate_gradient(&points[0]).unwrap(),
        [0.0, 0.0, 1.0]
    );
    assert_batch_matches_single(&model, &points);
    assert!(model.evaluate_scalars(&[]).unwrap().is_empty());
    assert!(model.evaluate_gradients(&[]).unwrap().is_empty());

    let spatial = model.data_bounds_and_resolution().unwrap();
    assert_eq!(spatial.bounds(), [0.0, 1.0, 0.0, 1.0, 0.0, 0.0]);
    assert_eq!(spatial.resolution.to_bits(), 0x3fd6_a09e_667f_3bcd);
}

#[test]
fn five_model_dispatch_matches_the_completed_internal_models() {
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    let batch = [witness.clone(), Point::new(-1.25, 2.5, -0.75).unwrap()];

    let lajaunie_parameters = parameters(ModelType::LajaunieApproach);
    let direct_lajaunie =
        fit_lajaunie_linear(&lajaunie_constraints(), &lajaunie_parameters).unwrap();
    let public_lajaunie = Builder::from_parameters(lajaunie_parameters)
        .set_constraints(lajaunie_constraints())
        .fit()
        .unwrap();
    assert_eq!(public_lajaunie.fit_branch(), FitBranch::Linear);
    assert_eq!(
        public_lajaunie.evaluate_scalar(&witness).unwrap(),
        direct_lajaunie.evaluate_scalar(&witness).unwrap()
    );
    assert_eq!(public_lajaunie.number_of_interfaces(), 3);
    assert_close(
        public_lajaunie.evaluate_scalar(&witness).unwrap(),
        0x3fb1_7133_d7c0_5e86,
        2.0e-12,
    );
    assert_gradient_close(
        public_lajaunie.evaluate_gradient(&witness).unwrap(),
        [
            0x3fa9_cef5_3cb2_2c2b,
            0xbfa6_0093_8a9a_2f5c,
            0xbfce_52af_dd98_707a,
        ],
    );
    assert_batch_matches_single(&public_lajaunie, &batch);

    let stratigraphic_parameters = parameters(ModelType::StratigraphicHorizons);
    let direct_stratigraphic =
        fit_stratigraphic(&stratigraphic_constraints(), &stratigraphic_parameters).unwrap();
    let public_stratigraphic = Builder::from_parameters(stratigraphic_parameters)
        .set_constraints(stratigraphic_constraints())
        .fit()
        .unwrap();
    assert_eq!(
        public_stratigraphic.fit_branch(),
        FitBranch::OrdinaryQuadratic
    );
    assert_eq!(
        public_stratigraphic.evaluate_scalar(&witness).unwrap(),
        direct_stratigraphic
            .evaluate_modified_scalar(&witness)
            .unwrap()
    );
    assert_eq!(public_stratigraphic.number_of_interfaces(), 3);
    assert_close(
        public_stratigraphic.evaluate_scalar(&witness).unwrap(),
        0xbfc2_25a5_e3c2_1329,
        2.0e-12,
    );
    assert_gradient_close(
        public_stratigraphic.evaluate_gradient(&witness).unwrap(),
        [
            0x3fba_f566_40fd_4ebd,
            0xbfa9_238f_043f_e56e,
            0xbfd0_454c_51d3_9c86,
        ],
    );
    assert_batch_matches_single(&public_stratigraphic, &batch);

    let continuous_parameters = parameters(ModelType::ContinuousProperty);
    let direct_continuous =
        fit_continuous_property(&continuous_constraints(), &continuous_parameters).unwrap();
    let public_continuous = Builder::from_parameters(continuous_parameters)
        .set_constraints(continuous_constraints())
        .fit()
        .unwrap();
    assert_eq!(public_continuous.fit_branch(), FitBranch::Linear);
    assert_eq!(
        public_continuous.evaluate_scalar(&witness).unwrap(),
        direct_continuous.evaluate_scalar(&witness).unwrap()
    );
    assert_eq!(public_continuous.number_of_interfaces(), 0);
    assert_close(
        public_continuous.evaluate_scalar(&witness).unwrap(),
        0x3ff2_1391_a86d_8a8d,
        2.0e-12,
    );
    assert_gradient_close(
        public_continuous.evaluate_gradient(&witness).unwrap(),
        [
            0xbfe7_f9ec_d70b_dd0f,
            0x3fd9_8d71_533f_4e3e,
            0x3fba_5e42_1840_a4e0,
        ],
    );
    assert_batch_matches_single(&public_continuous, &batch);

    let vector_parameters = parameters(ModelType::VectorField);
    let direct_vector = fit_vector_field(&vector_constraints(), &vector_parameters).unwrap();
    let public_vector = Builder::from_parameters(vector_parameters)
        .set_constraints(vector_constraints())
        .fit()
        .unwrap();
    assert_eq!(public_vector.fit_branch(), FitBranch::Linear);
    assert_eq!(
        public_vector.evaluate_scalar(&witness).unwrap(),
        direct_vector.evaluate_potential(&witness).unwrap()
    );
    assert_eq!(
        public_vector.evaluate_gradient(&witness).unwrap(),
        direct_vector.evaluate_gradient(&witness).unwrap()
    );
    assert_eq!(public_vector.number_of_interfaces(), 0);
    assert_close(
        public_vector.evaluate_scalar(&witness).unwrap(),
        0xbfc2_6789_da3f_9380,
        2.0e-12,
    );
    assert_gradient_close(
        public_vector.evaluate_gradient(&witness).unwrap(),
        [
            0x3fcd_628d_517a_9eaf,
            0x3fca_ca92_78ec_a8f8,
            0x3fe6_0baa_7247_f0f4,
        ],
    );
    assert_batch_matches_single(&public_vector, &batch);
}

#[test]
fn builder_dispatches_single_surface_inequality_and_owns_configuration() {
    let mut builder = Builder::from_parameters(parameters(ModelType::SingleSurface));
    builder.set_constraints(Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 30.0).unwrap(),
            Interface::new(5.0, 1.0, 1.0, 20.0).unwrap(),
            Interface::new(7.0, -1.0, 2.0, 10.0).unwrap(),
        ],
        inequalities: vec![
            Inequality::new(3.0, -2.0, 2.0, 1.0).unwrap(),
            Inequality::new(4.0, 2.0, -1.0, -1.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        ..Constraints::default()
    });
    builder.set_greedy_algorithm(false, 0.25, 5.0);
    let source_copy = builder.constraints().clone();
    let fitted = builder.fit().unwrap();
    assert_eq!(fitted.fit_branch(), FitBranch::OrdinaryQuadratic);
    assert_eq!(fitted.constraints().inequalities.len(), 2);
    assert_eq!(source_copy.inequalities.len(), 2);
    assert!(fitted.parameters().use_greedy);
}

#[test]
fn dynamic_constraint_tables_are_dimension_checked_and_atomic() {
    let mut builder = Builder::new(ModelType::SingleSurface);
    builder.add_interface_xyz(9.0, 8.0, 7.0, 6.0).unwrap();
    let wrong = DenseMatrix::from_row_major(1, 3, vec![0.0; 3]).unwrap();
    assert_eq!(
        builder.set_interface_constraint_matrix(&wrong).unwrap_err(),
        BuildError::IncorrectArrayDimensions {
            rows: 1,
            columns: 3,
            expected_columns: 4,
        }
    );
    assert_eq!(builder.constraints().interfaces.len(), 1);

    let empty = DenseMatrix::from_row_major(0, 4, Vec::new()).unwrap();
    assert_eq!(
        builder.set_interface_constraint_matrix(&empty).unwrap_err(),
        BuildError::IncorrectArrayDimensions {
            rows: 0,
            columns: 4,
            expected_columns: 4,
        }
    );
    assert_eq!(builder.constraints().interfaces.len(), 1);

    let tangents = DenseMatrix::from_row_major(
        2,
        6,
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0],
    )
    .unwrap();
    builder.set_tangent_constraint_matrix(&tangents).unwrap();
    assert_eq!(builder.constraints().tangents.len(), 2);
    assert!(builder.constraints().planars.is_empty());
    assert_eq!(builder.tangent_constraint_matrix(), tangents);
}

#[test]
fn fitted_model_is_send_sync_and_supports_concurrent_read_only_evaluation() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<FittedModel>();

    let model = Arc::new(single_surface_builder().fit().unwrap());
    let expected = model
        .evaluate_scalar(&Point::new(0.25, 0.75, 0.5).unwrap())
        .unwrap();
    let handles = (0..8)
        .map(|_| {
            let model = Arc::clone(&model);
            thread::spawn(move || {
                let point = Point::new(0.25, 0.75, 0.5).unwrap();
                (
                    model.evaluate_scalar(&point).unwrap(),
                    model.evaluate_gradient(&point).unwrap(),
                )
            })
        })
        .collect::<Vec<_>>();
    for handle in handles {
        assert_eq!(handle.join().unwrap(), (expected, [0.0, 0.0, 1.0]));
    }
}

#[test]
fn invalid_parameter_and_evaluation_dimensions_fail_before_numeric_work() {
    let mut invalid = single_surface_builder();
    invalid.set_rbf_shape_parameter(f64::NAN);
    assert_eq!(
        invalid.fit().unwrap_err(),
        BuildError::NonFiniteParameter("shape_parameter")
    );

    let fitted = single_surface_builder().fit().unwrap();
    let empty = DenseMatrix::from_row_major(0, 3, Vec::new()).unwrap();
    assert_eq!(
        fitted.evaluate_scalar_matrix(&empty).unwrap_err(),
        EvaluationError::IncorrectArrayDimensions {
            rows: 0,
            columns: 3,
            expected_columns: 3,
        }
    );
    let wrong = DenseMatrix::from_row_major(1, 2, vec![0.0, 0.0]).unwrap();
    assert_eq!(
        fitted.evaluate_gradient_matrix(&wrong).unwrap_err(),
        EvaluationError::IncorrectArrayDimensions {
            rows: 1,
            columns: 2,
            expected_columns: 3,
        }
    );
}

#[test]
fn all_restricted_model_dispatches_are_publicly_reachable() {
    let source = Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 30.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 30.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        ..Constraints::default()
    };
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();

    for model_type in [
        ModelType::SingleSurface,
        ModelType::LajaunieApproach,
        ModelType::StratigraphicHorizons,
    ] {
        let mut configuration = parameters(model_type);
        configuration.set_restricted_range(true, 0.25, 8.0);
        let fitted = Builder::from_parameters(configuration)
            .set_constraints(source.clone())
            .fit()
            .unwrap();
        assert_eq!(fitted.model_type(), model_type);
        assert_eq!(fitted.fit_branch(), FitBranch::RestrictedRange);
        assert!(fitted.evaluate_scalar(&witness).unwrap().is_finite());
        assert!(fitted
            .evaluate_gradient(&witness)
            .unwrap()
            .into_iter()
            .all(f64::is_finite));
    }
}

#[test]
fn surfe_model_codes_and_exact_kernel_names_share_the_typed_builder() {
    let expected = [
        ModelType::SingleSurface,
        ModelType::LajaunieApproach,
        ModelType::VectorField,
        ModelType::StratigraphicHorizons,
        ModelType::ContinuousProperty,
    ];
    for (code, model_type) in (1..=5).zip(expected) {
        assert_eq!(
            Builder::from_surfe_model_code(code)
                .unwrap()
                .parameters()
                .model_type,
            model_type
        );
    }
    assert_eq!(
        Builder::from_surfe_model_code(0).unwrap_err(),
        BuildError::Surfe(georbf::Error::UnknownModel)
    );

    let mut builder = Builder::new(ModelType::SingleSurface);
    builder.set_rbf_kernel_name("Thin Plate Spline").unwrap();
    assert_eq!(builder.parameters().basis_type, RbfKernel::ThinPlateSpline);
    assert_eq!(
        builder
            .set_rbf_kernel_name("thin plate spline")
            .unwrap_err(),
        BuildError::Surfe(georbf::Error::UnknownRbf)
    );
    assert_eq!(builder.parameters().basis_type, RbfKernel::ThinPlateSpline);
}
