use georbf::{
    fit_vector_field, AnisotropyError, AssemblyError, Constraints, Inequality, Interface,
    KernelError, LuSolveErrorKind, ModelType, Parameters, Planar, Point, RbfKernel, Tangent,
    VectorFieldError,
};

fn parameters(kernel: RbfKernel) -> Parameters {
    Parameters {
        model_type: ModelType::VectorField,
        basis_type: kernel,
        shape_parameter: 1.7,
        ..Parameters::default()
    }
}

fn planar(x: f64, y: f64, z: f64, normal: [f64; 3]) -> Planar {
    Planar::from_normal(x, y, z, normal[0], normal[1], normal[2]).unwrap()
}

fn ordinary_planars() -> Vec<Planar> {
    vec![
        planar(1.0, 2.0, 3.0, [0.6, 0.0, 0.8]),
        planar(-2.0, 1.0, 0.5, [0.0, 0.8, 0.6]),
        planar(0.5, -1.0, 2.0, [0.36, -0.48, 0.8]),
    ]
}

fn anisotropic_planars() -> Vec<Planar> {
    vec![
        planar(1.0, 2.0, 3.0, [1.0, 0.0, 0.0]),
        planar(-2.0, 1.0, 0.5, [0.0, 1.0, 0.0]),
        planar(0.5, -1.0, 2.0, [0.0, 0.0, 1.0]),
        planar(3.0, 0.5, -2.0, [0.36, -0.48, 0.8]),
        planar(-1.5, -2.5, 1.0, [-0.48, 0.64, 0.6]),
    ]
}

fn ordinary_constraints() -> Constraints {
    Constraints {
        planars: ordinary_planars(),
        ..Constraints::default()
    }
}

fn query_points() -> Vec<Point> {
    vec![
        Point::new(0.25, -0.5, 1.25).unwrap(),
        Point::new(-1.25, 2.5, -0.75).unwrap(),
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

fn assert_matrix_symmetric(model: &georbf::VectorFieldModel) {
    let matrix = model.interpolation_matrix();
    for row in 0..matrix.rows() {
        for column in 0..matrix.cols() {
            assert_close(
                matrix.get(row, column).unwrap(),
                matrix.get(column, row).unwrap(),
                2.0e-12,
            );
        }
    }
}

#[test]
fn cubic_hessian_system_rhs_lu_potential_and_gradient_match_frozen_surfe() {
    let model = fit_vector_field(&ordinary_constraints(), &parameters(RbfKernel::Cubic)).unwrap();

    assert_eq!(model.layout().matrix_size(), 9);
    assert_eq!(model.layout().internal_parameters().n_planar, 3);
    assert_eq!(model.layout().internal_parameters().n_interface, 0);
    assert_eq!(model.layout().internal_parameters().n_inequality, 0);
    assert_eq!(model.layout().internal_parameters().n_tangent, 0);
    assert!(!model.layout().internal_parameters().poly_term);
    assert_eq!(model.layout().internal_parameters().n_poly_terms, 0);
    assert_eq!(
        hash_matrix(
            model.interpolation_matrix().rows(),
            model.interpolation_matrix().cols(),
            model.interpolation_matrix().data(),
        ),
        0xeb24_aafe_28b2_dd7b
    );
    assert_eq!(
        hash_values(model.right_hand_side().values()),
        0x963b_9b6a_4ac7_7002
    );
    assert!(model.lu_solution().unwrap().residual().accepted());
    assert_matrix_symmetric(&model);

    let expected_scalars = [
        f64::from_bits(0x3fd7_7112_10c5_0ee5),
        f64::from_bits(0xbfe0_a2d4_1a86_2aca),
    ];
    let expected_gradients = [
        [
            f64::from_bits(0x3fd7_c005_7bd3_fac5),
            f64::from_bits(0xbfc8_dc24_161c_51ee),
            f64::from_bits(0x3fe8_f91d_7cf5_dc75),
        ],
        [
            f64::from_bits(0xbfb1_50da_1d9b_c70a),
            f64::from_bits(0x3fd5_ff0f_1c97_6b02),
            f64::from_bits(0x3ff0_7e44_97a6_18e3),
        ],
    ];
    let points = query_points();
    let scalars = model.evaluate_potentials(&points).unwrap();
    let gradients = model.evaluate_gradients(&points).unwrap();
    for query in 0..points.len() {
        assert_close(scalars[query], expected_scalars[query], 2.0e-12);
        for component in 0..3 {
            assert_close(
                gradients[query][component],
                expected_gradients[query][component],
                4.0e-11,
            );
        }
    }

    for planar in &model.constraints().planars {
        let gradient = model.evaluate_gradient(planar.point()).unwrap();
        for (component, value) in gradient.iter().enumerate() {
            assert_close(*value, planar.normal()[component], 4.0e-11);
        }
    }
}

#[test]
fn gaussian_field_matches_oracle_and_gradient_is_the_potential_derivative() {
    let model =
        fit_vector_field(&ordinary_constraints(), &parameters(RbfKernel::Gaussian)).unwrap();
    assert_eq!(
        hash_matrix(
            model.interpolation_matrix().rows(),
            model.interpolation_matrix().cols(),
            model.interpolation_matrix().data(),
        ),
        0x30b9_abf6_7cf1_e90f
    );
    assert!(model.lu_solution().unwrap().residual().accepted());
    let point = &query_points()[0];
    let scalar = model.evaluate_potential(point).unwrap();
    let gradient = model.evaluate_gradient(point).unwrap();
    assert_close(scalar, f64::from_bits(0xbfb2_fd2b_3b76_2720), 2.0e-12);
    let expected = [
        f64::from_bits(0xbfb4_16a5_69c2_c060),
        f64::from_bits(0x3fc6_89e4_d00d_d75d),
        f64::from_bits(0xbfd0_7ee1_0c50_bd40),
    ];
    for component in 0..3 {
        assert_close(gradient[component], expected[component], 4.0e-11);
    }

    let step = 1.0e-6;
    for component in 0..3 {
        let mut positive = [point.x(), point.y(), point.z()];
        let mut negative = positive;
        positive[component] += step;
        negative[component] -= step;
        let positive = Point::new(positive[0], positive[1], positive[2]).unwrap();
        let negative = Point::new(negative[0], negative[1], negative[2]).unwrap();
        let numerical = (model.evaluate_potential(&positive).unwrap()
            - model.evaluate_potential(&negative).unwrap())
            / (2.0 * step);
        assert_close(numerical, gradient[component], 2.0e-7);
    }
}

#[test]
fn global_anisotropy_hessian_system_and_field_match_frozen_surfe() {
    let constraints = Constraints {
        planars: anisotropic_planars(),
        ..Constraints::default()
    };
    let mut anisotropic = parameters(RbfKernel::Cubic);
    anisotropic.model_global_anisotropy = true;
    let model = fit_vector_field(&constraints, &anisotropic).unwrap();
    assert_eq!(model.layout().matrix_size(), 15);
    assert_eq!(
        hash_matrix(
            model.interpolation_matrix().rows(),
            model.interpolation_matrix().cols(),
            model.interpolation_matrix().data(),
        ),
        0x6a6f_4483_455c_56df
    );
    assert_eq!(
        hash_values(model.right_hand_side().values()),
        0x563b_3fc1_bf08_9ce3
    );
    assert!(model.lu_solution().unwrap().residual().accepted());
    let scalars = model.evaluate_potentials(&query_points()).unwrap();
    let gradients = model.evaluate_gradients(&query_points()).unwrap();
    let expected_scalars = [
        f64::from_bits(0xbfd7_6d52_da18_8a54),
        f64::from_bits(0x3fec_1d31_a634_39e6),
    ];
    let expected_gradients = [
        [
            f64::from_bits(0xbfb3_7373_06bc_ad7e),
            f64::from_bits(0x3fd3_dc0c_813a_fbf4),
            f64::from_bits(0x3fe3_5064_b9f2_698e),
        ],
        [
            f64::from_bits(0xbfc9_d29f_6791_1478),
            f64::from_bits(0x3fe1_5ad9_c372_ed08),
            f64::from_bits(0x3fe4_eb8e_abba_ef11),
        ],
    ];
    for query in 0..2 {
        assert_close(scalars[query], expected_scalars[query], 2.0e-12);
        for component in 0..3 {
            assert_close(
                gradients[query][component],
                expected_gradients[query][component],
                4.0e-11,
            );
        }
    }
}

#[test]
fn all_eight_differentiable_isotropic_kernels_fit_the_frozen_hessian_path() {
    let cases = [
        (RbfKernel::Cubic, 0xeb24_aafe_28b2_dd7b),
        (RbfKernel::Gaussian, 0x30b9_abf6_7cf1_e90f),
        (RbfKernel::Multiquadric, 0x44ec_803d_0366_3d05),
        (RbfKernel::MultiquadricCubic, 0x9138_bbde_93b3_83f5),
        (RbfKernel::ThinPlateSpline, 0x7e4b_3d8a_8fcf_7ca3),
        (RbfKernel::InverseMultiquadric, 0xcccc_4628_9fc4_f36f),
        (RbfKernel::WendlandC2, 0xd83d_22b3_7a43_e063),
        (RbfKernel::MaternC4, 0xca73_57bb_4626_39b7),
    ];
    for (kernel, expected_matrix_hash) in cases {
        let model = fit_vector_field(&ordinary_constraints(), &parameters(kernel)).unwrap();
        assert_eq!(
            hash_matrix(
                model.interpolation_matrix().rows(),
                model.interpolation_matrix().cols(),
                model.interpolation_matrix().data(),
            ),
            expected_matrix_hash,
            "{kernel:?}"
        );
        assert!(model.lu_solution().unwrap().residual().accepted());
        for planar in &model.constraints().planars {
            let gradient = model.evaluate_gradient(planar.point()).unwrap();
            for (component, value) in gradient.iter().enumerate() {
                assert_close(*value, planar.normal()[component], 4.0e-10);
            }
        }
    }
}

#[test]
fn inactive_categories_and_parameters_are_exactly_ignored_after_cleaning() {
    let baseline =
        fit_vector_field(&ordinary_constraints(), &parameters(RbfKernel::Cubic)).unwrap();
    let mut constraints = ordinary_constraints();
    constraints.inequalities = vec![
        Inequality::new(8.0, 1.0, -2.0, 50.0).unwrap(),
        Inequality::new(-4.0, 2.0, 1.0, -30.0).unwrap(),
    ];
    constraints.interfaces = vec![
        Interface::new(7.0, -3.0, 2.0, 100.0).unwrap(),
        Interface::new(-5.0, 4.0, 0.0, -100.0).unwrap(),
    ];
    constraints.tangents = vec![Tangent::new(2.0, -4.0, 1.0, 0.2, 0.7, 0.1).unwrap()];
    constraints
        .planars
        .push(planar(1.0, 2.0, 3.0, [-1.0, 0.0, 0.0]));
    let mut ignored = parameters(RbfKernel::Cubic);
    ignored.polynomial_order = 2;
    ignored.use_regression_smoothing = true;
    ignored.smoothing_amount = 9.5;
    ignored.use_restricted_range = true;
    ignored.interface_uncertainty = 0.25;
    ignored.angular_uncertainty = 12.0;
    ignored.use_greedy = true;

    let model = fit_vector_field(&constraints, &ignored).unwrap();
    assert_eq!(model.collocation_removal().planars, 1);
    assert_eq!(model.constraints().planars.len(), 3);
    assert_eq!(model.constraints().inequalities.len(), 2);
    assert_eq!(model.constraints().interfaces.len(), 2);
    assert_eq!(model.constraints().tangents.len(), 1);
    assert_eq!(model.layout().source_counts().inequalities, 2);
    assert_eq!(model.layout().source_counts().interfaces, 2);
    assert_eq!(model.layout().source_counts().tangents, 1);
    assert_eq!(model.layout().internal_parameters().n_inequality, 0);
    assert_eq!(model.layout().internal_parameters().n_interface, 0);
    assert_eq!(model.layout().internal_parameters().n_tangent, 0);
    assert_eq!(model.layout().internal_parameters().n_poly_terms, 0);
    assert_eq!(
        model.interpolation_matrix(),
        baseline.interpolation_matrix()
    );
    assert_eq!(model.right_hand_side(), baseline.right_hand_side());
    assert_eq!(
        model.lu_solution().unwrap().weights(),
        baseline.lu_solution().unwrap().weights()
    );
    let points = query_points();
    assert_eq!(
        model.evaluate_potentials(&points).unwrap(),
        baseline.evaluate_potentials(&points).unwrap()
    );
    assert_eq!(
        model.evaluate_gradients(&points).unwrap(),
        baseline.evaluate_gradients(&points).unwrap()
    );
    for (index, point) in points.iter().enumerate() {
        assert_eq!(
            model.evaluate_potentials(&points).unwrap()[index],
            model.evaluate_potential(point).unwrap()
        );
        assert_eq!(
            model.evaluate_gradients(&points).unwrap()[index],
            model.evaluate_gradient(point).unwrap()
        );
    }
}

#[test]
fn frozen_empty_planar_fit_is_a_zero_potential_and_zero_gradient_success() {
    let model = fit_vector_field(&Constraints::default(), &parameters(RbfKernel::Cubic)).unwrap();
    assert_eq!(model.layout().matrix_size(), 0);
    assert!(model.lu_solution().is_none());
    assert!(model.interpolation_matrix().data().is_empty());
    assert!(model.right_hand_side().is_empty());
    for point in query_points() {
        assert_eq!(model.evaluate_potential(&point).unwrap().to_bits(), 0);
        assert_eq!(model.evaluate_gradient(&point).unwrap(), [0.0; 3]);
    }
}

#[test]
fn degenerate_and_unsupported_paths_have_stable_source_compatible_errors() {
    let wrong = Parameters {
        model_type: ModelType::ContinuousProperty,
        ..parameters(RbfKernel::Cubic)
    };
    assert!(matches!(
        fit_vector_field(&ordinary_constraints(), &wrong),
        Err(VectorFieldError::WrongModel)
    ));

    let singleton = Constraints {
        planars: vec![planar(1.0, 2.0, 3.0, [0.6, 0.0, 0.8])],
        ..Constraints::default()
    };
    match fit_vector_field(&singleton, &parameters(RbfKernel::Cubic)) {
        Err(VectorFieldError::Lu(error)) => {
            assert_eq!(error.kind(), LuSolveErrorKind::SingularSystem);
            assert!(error.attempted());
            assert_eq!(error.surfe_error(), georbf::Error::LinearSolverFailure);
        }
        result => panic!("expected frozen linear-solver failure, got {result:?}"),
    }

    match fit_vector_field(&ordinary_constraints(), &parameters(RbfKernel::Linear)) {
        Err(VectorFieldError::Assembly(AssemblyError::Kernel(
            KernelError::LinearDerivativeUnavailable,
        ))) => {}
        result => panic!("expected frozen -666 derivative classification, got {result:?}"),
    }

    let mut anisotropic = parameters(RbfKernel::Cubic);
    anisotropic.model_global_anisotropy = true;
    assert!(matches!(
        fit_vector_field(&singleton, &anisotropic),
        Err(VectorFieldError::Anisotropy(
            AnisotropyError::InsufficientPlanars
        ))
    ));
}
