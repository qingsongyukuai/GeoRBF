use georbf::{
    fit_continuous_property, AnisotropyError, Constraints, ContinuousPropertyError, Error,
    Inequality, Interface, KernelError, ModelType, Parameters, Planar, Point, RbfKernel, Tangent,
};

fn parameters() -> Parameters {
    Parameters {
        model_type: ModelType::ContinuousProperty,
        basis_type: RbfKernel::Cubic,
        shape_parameter: 1.7,
        ..Parameters::default()
    }
}

fn interfaces() -> Vec<Interface> {
    vec![
        Interface::new(0.0, 0.0, 0.0, 1.0).unwrap(),
        Interface::new(2.0, 0.0, 0.0, -0.5).unwrap(),
        Interface::new(0.0, 3.0, 0.0, 2.25).unwrap(),
        Interface::new(0.0, 0.0, 4.0, -1.75).unwrap(),
        Interface::new(1.5, -2.0, 2.5, 0.75).unwrap(),
    ]
}

fn constraints() -> Constraints {
    Constraints {
        interfaces: interfaces(),
        ..Constraints::default()
    }
}

fn query_points() -> Vec<Point> {
    vec![
        Point::new(0.5, 0.75, 1.25).unwrap(),
        Point::new(-1.0, 2.0, 0.5).unwrap(),
    ]
}

fn hash_values(values: &[f64]) -> u64 {
    let mut hash = 1_469_598_103_934_665_603_u64;
    mix_u64(&mut hash, values.len() as u64);
    for value in values {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn hash_matrix(rows: usize, columns: usize, values: &[f64]) -> u64 {
    let mut hash = 1_469_598_103_934_665_603_u64;
    mix_u64(&mut hash, rows as u64);
    mix_u64(&mut hash, columns as u64);
    for value in values {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn mix_u64(hash: &mut u64, value: u64) {
    for byte in 0..8 {
        *hash ^= (value >> (8 * byte)) & 0xff;
        *hash = hash.wrapping_mul(1_099_511_628_211);
    }
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let scale = 1.0_f64.max(actual.abs()).max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn public_reachable_cubic_fit_matrix_rhs_lu_and_field_match_frozen_surfe() {
    let model = fit_continuous_property(&constraints(), &parameters()).unwrap();

    assert_eq!(model.layout().matrix_size(), 5);
    assert_eq!(
        model.layout().internal_parameters().n_interface,
        model.constraints().interfaces.len()
    );
    assert_eq!(model.layout().internal_parameters().n_planar, 0);
    assert_eq!(model.layout().internal_parameters().n_tangent, 0);
    assert_eq!(model.layout().internal_parameters().n_inequality, 0);
    assert!(!model.layout().internal_parameters().poly_term);
    assert!(!model.layout().internal_parameters().modified_basis);
    assert_eq!(
        hash_matrix(
            model.interpolation_matrix().rows(),
            model.interpolation_matrix().cols(),
            model.interpolation_matrix().data(),
        ),
        0x6253_1c8a_e9a7_fb27
    );
    assert_eq!(
        hash_values(model.right_hand_side().values()),
        0xcbfe_51c8_e0af_6568
    );
    assert_eq!(
        hash_values(model.lu_solution().weights().values()),
        0x793c_4439_db53_b5ac
    );

    let residual = model.lu_solution().residual();
    assert!(residual.accepted());
    assert!(residual.relative_l2().is_finite());
    let points = query_points();
    let scalars = model.evaluate_scalars(&points).unwrap();
    let gradients = model.evaluate_gradients(&points).unwrap();
    let expected_scalars = [
        f64::from_bits(0x3fee_9b3f_e3a6_471d),
        f64::from_bits(0x4007_34c5_e50d_e036),
    ];
    let expected_gradients = [
        [
            f64::from_bits(0xbfe6_f7c2_762a_9b60),
            f64::from_bits(0x3fad_dc02_75bd_ce3c),
            f64::from_bits(0xbfc1_52a3_664d_66af),
        ],
        [
            f64::from_bits(0xbfea_9a6b_73ca_320e),
            f64::from_bits(0x3fdc_846e_71cd_cfb8),
            f64::from_bits(0xbfd9_8886_a0ec_e730),
        ],
    ];
    for index in 0..2 {
        assert_close(scalars[index], expected_scalars[index], 2.0e-12);
        for component in 0..3 {
            assert_close(
                gradients[index][component],
                expected_gradients[index][component],
                4.0e-11,
            );
        }
    }
}

#[test]
fn public_batch_evaluation_matches_single_point_order() {
    let model = fit_continuous_property(&constraints(), &parameters()).unwrap();
    let points = query_points();
    let scalar_batch = model.evaluate_scalars(&points).unwrap();
    let gradient_batch = model.evaluate_gradients(&points).unwrap();
    for (index, point) in points.iter().enumerate() {
        assert_eq!(scalar_batch[index], model.evaluate_scalar(point).unwrap());
        assert_eq!(
            gradient_batch[index],
            model.evaluate_gradient(point).unwrap()
        );
    }
}

#[test]
fn gaussian_value_matrix_and_field_match_the_frozen_reachable_branch() {
    let mut gaussian = parameters();
    gaussian.basis_type = RbfKernel::Gaussian;
    let model = fit_continuous_property(&constraints(), &gaussian).unwrap();
    assert_eq!(
        hash_matrix(
            model.interpolation_matrix().rows(),
            model.interpolation_matrix().cols(),
            model.interpolation_matrix().data(),
        ),
        0xb418_df9e_6ac2_9bc2
    );
    assert_eq!(
        hash_values(model.lu_solution().weights().values()),
        0x36e0_1b04_2cef_a9b4
    );
    let points = query_points();
    let scalars = model.evaluate_scalars(&points).unwrap();
    let gradients = model.evaluate_gradients(&points).unwrap();
    let expected_scalars = [
        f64::from_bits(0x3f51_1845_47c7_5ccb),
        f64::from_bits(0x3f6b_a4f1_614b_b58e),
    ];
    let expected_gradients = [
        [
            f64::from_bits(0xbf68_dad6_0c6e_944b),
            f64::from_bits(0xbf72_86af_e8b6_13d9),
            f64::from_bits(0xbf7e_e0a3_df1e_a1d3),
        ],
        [
            f64::from_bits(0x3f93_f917_5b6a_8a0f),
            f64::from_bits(0x3f93_f7eb_c98c_502d),
            f64::from_bits(0xbf83_f917_5b6a_8a36),
        ],
    ];
    for index in 0..2 {
        assert_close(scalars[index], expected_scalars[index], 2.0e-12);
        for component in 0..3 {
            assert_close(
                gradients[index][component],
                expected_gradients[index][component],
                4.0e-11,
            );
        }
    }
}

#[test]
fn inequality_and_inactive_parameter_families_are_reachable_but_ignored() {
    let baseline = fit_continuous_property(&constraints(), &parameters()).unwrap();
    let mut extra = constraints();
    extra.inequalities = vec![
        Inequality::new(-3.0, 1.0, 2.0, 100.0).unwrap(),
        Inequality::new(6.0, -2.0, 1.0, -100.0).unwrap(),
    ];
    let mut ignored = parameters();
    ignored.polynomial_order = 2;
    ignored.use_regression_smoothing = true;
    ignored.smoothing_amount = 9.5;
    ignored.use_restricted_range = true;
    ignored.interface_uncertainty = 0.25;
    ignored.angular_uncertainty = 12.0;
    ignored.use_greedy = true;

    let model = fit_continuous_property(&extra, &ignored).unwrap();
    assert_eq!(model.constraints().inequalities.len(), 2);
    assert_eq!(model.layout().source_counts().inequalities, 2);
    assert_eq!(model.layout().internal_parameters().n_inequality, 0);
    assert_eq!(model.layout().internal_parameters().n_poly_terms, 0);
    assert!(!model.layout().internal_parameters().restricted_range);
    assert_eq!(
        model.interpolation_matrix(),
        baseline.interpolation_matrix()
    );
    assert_eq!(model.right_hand_side(), baseline.right_hand_side());
    assert_eq!(
        model.lu_solution().weights(),
        baseline.lu_solution().weights()
    );
    for point in query_points() {
        assert_eq!(
            model.evaluate_scalar(&point).unwrap(),
            baseline.evaluate_scalar(&point).unwrap()
        );
        assert_eq!(
            model.evaluate_gradient(&point).unwrap(),
            baseline.evaluate_gradient(&point).unwrap()
        );
    }
}

#[test]
fn exact_collocation_cleaning_precedes_the_interface_value_fit() {
    let mut duplicated = constraints();
    duplicated
        .interfaces
        .push(Interface::new(0.0, 0.0, 0.0, 99.0).unwrap());
    let model = fit_continuous_property(&duplicated, &parameters()).unwrap();
    assert_eq!(model.collocation_removal().interfaces, 1);
    assert_eq!(model.constraints().interfaces.len(), 5);
    assert_eq!(model.right_hand_side().values()[0], 1.0);
}

#[test]
fn linear_value_path_fits_but_gradient_preserves_the_frozen_sentinel() {
    let mut linear = parameters();
    linear.basis_type = RbfKernel::Linear;
    let model = fit_continuous_property(&constraints(), &linear).unwrap();
    assert!(model
        .evaluate_scalar(&query_points()[0])
        .unwrap()
        .is_finite());
    assert_eq!(
        model.evaluate_gradient(&query_points()[0]),
        Err(ContinuousPropertyError::Evaluation(
            KernelError::LinearDerivativeUnavailable
        ))
    );
}

#[test]
fn source_ub_and_unavailable_anisotropy_are_safe_typed_failures() {
    let mut with_planar = constraints();
    with_planar
        .planars
        .push(Planar::from_normal(1.0, 1.0, 1.0, 0.0, 0.0, 1.0).unwrap());
    assert_eq!(
        fit_continuous_property(&with_planar, &parameters()).unwrap_err(),
        ContinuousPropertyError::EqualityVectorOutOfBounds {
            planar_count: 1,
            tangent_count: 0,
        }
    );

    let mut with_tangent = constraints();
    with_tangent
        .tangents
        .push(Tangent::new(1.0, 1.0, 1.0, 1.0, 0.0, 0.0).unwrap());
    assert_eq!(
        fit_continuous_property(&with_tangent, &parameters()).unwrap_err(),
        ContinuousPropertyError::EqualityVectorOutOfBounds {
            planar_count: 0,
            tangent_count: 1,
        }
    );

    let mut anisotropic = parameters();
    anisotropic.model_global_anisotropy = true;
    assert_eq!(
        fit_continuous_property(&constraints(), &anisotropic).unwrap_err(),
        ContinuousPropertyError::Anisotropy(AnisotropyError::InsufficientPlanars)
    );
}

#[test]
fn wrong_model_and_missing_interface_errors_are_explicit() {
    let mut wrong = parameters();
    wrong.model_type = ModelType::VectorField;
    assert_eq!(
        fit_continuous_property(&constraints(), &wrong).unwrap_err(),
        ContinuousPropertyError::WrongModel
    );
    assert_eq!(
        fit_continuous_property(&Constraints::default(), &parameters()).unwrap_err(),
        ContinuousPropertyError::Surfe(Error::NoInterfaceData)
    );
}
