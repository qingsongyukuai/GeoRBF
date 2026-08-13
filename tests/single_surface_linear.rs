use georbf::{
    fit_single_surface_linear, AnisotropyError, Constraints, DenseMatrix, DenseVector, Error,
    Interface, ModelType, Parameters, Planar, Point, RbfKernel, SingleSurfaceLinearError,
};

fn smoke_constraints() -> Constraints {
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

fn parameters() -> Parameters {
    Parameters {
        model_type: ModelType::SingleSurface,
        use_interface: true,
        use_planar: true,
        basis_type: RbfKernel::Cubic,
        shape_parameter: 2.0,
        polynomial_order: 1,
        ..Parameters::default()
    }
}

fn hash_matrix(matrix: &DenseMatrix) -> u64 {
    let mut hash = 1_469_598_103_934_665_603_u64;
    mix_u64(&mut hash, matrix.rows() as u64);
    mix_u64(&mut hash, matrix.cols() as u64);
    for value in matrix.data() {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn hash_vector(vector: &DenseVector) -> u64 {
    let mut hash = 1_469_598_103_934_665_603_u64;
    mix_u64(&mut hash, vector.len() as u64);
    for value in vector.values() {
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

fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    assert!(
        actual.is_finite() && delta <= absolute + relative * scale,
        "actual={actual:?} expected={expected:?} delta={delta:?}"
    );
}

#[test]
fn frozen_t03_public_smoke_matches_layout_matrix_rhs_weights_and_field() {
    let model = fit_single_surface_linear(&smoke_constraints(), &parameters()).unwrap();

    assert_eq!(model.layout().matrix_size(), 11);
    assert_eq!(
        hash_matrix(model.interpolation_matrix()),
        0x908a_282d_5b21_e173
    );
    assert_eq!(hash_vector(model.right_hand_side()), 0xe3c8_2ac7_f225_5b15);
    assert_eq!(
        hash_vector(model.lu_solution().weights()),
        0xadba_05c4_7e42_f835
    );
    assert_eq!(model.lu_solution().residual().l2(), 0.0);

    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    assert_eq!(
        model.evaluate_scalar(&witness).unwrap().to_bits(),
        0.5_f64.to_bits()
    );
    assert_eq!(model.evaluate_gradient(&witness).unwrap(), [0.0, 0.0, 1.0]);
}

#[test]
fn every_differentiable_isotropic_kernel_matches_the_frozen_system_and_field() {
    let cases = [
        (RbfKernel::Cubic, 0x908a_282d_5b21_e173),
        (RbfKernel::Gaussian, 0x8b36_850d_616d_7d17),
        (RbfKernel::Multiquadric, 0x7357_f30a_f79d_8bee),
        (RbfKernel::MultiquadricCubic, 0xa8a3_6fd2_c507_9f20),
        (RbfKernel::ThinPlateSpline, 0xc9a8_1f8c_075e_40cf),
        (RbfKernel::InverseMultiquadric, 0x67ce_b249_02a9_0ace),
        (RbfKernel::WendlandC2, 0x935c_fbbc_5f34_d80f),
        (RbfKernel::MaternC4, 0x4652_7884_323e_d557),
    ];
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    for (kernel, matrix_hash) in cases {
        let mut parameters = parameters();
        parameters.basis_type = kernel;
        let model = fit_single_surface_linear(&smoke_constraints(), &parameters).unwrap();
        let actual_matrix_hash = hash_matrix(model.interpolation_matrix());
        assert_eq!(actual_matrix_hash, matrix_hash, "{kernel:?}");
        assert_eq!(hash_vector(model.right_hand_side()), 0xe3c8_2ac7_f225_5b15);
        assert_eq!(
            hash_vector(model.lu_solution().weights()),
            0xadba_05c4_7e42_f835
        );
        assert_eq!(model.evaluate_scalar(&witness).unwrap(), 0.5);
        assert_eq!(model.evaluate_gradient(&witness).unwrap(), [0.0, 0.0, 1.0]);
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

fn quadratic_value(x: f64, y: f64, z: f64) -> f64 {
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
}

fn quadratic_constraints() -> Constraints {
    let points = [
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
    ];
    Constraints {
        interfaces: points
            .into_iter()
            .map(|[x, y, z]| Interface::new(x, y, z, quadratic_value(x, y, z)).unwrap())
            .collect(),
        ..Constraints::default()
    }
}

#[test]
fn zero_first_and_second_order_polynomial_paths_match_frozen_surfe() {
    let witness = Point::new(0.25, 0.5, 0.75).unwrap();

    let mut zero = parameters();
    zero.use_planar = false;
    zero.polynomial_order = 0;
    let zero_model = fit_single_surface_linear(&constant_constraints(), &zero).unwrap();
    assert_eq!(
        hash_matrix(zero_model.interpolation_matrix()),
        0x95dc_9d5c_8204_c457
    );
    assert_eq!(
        hash_vector(zero_model.right_hand_side()),
        0x9343_eb7f_8ed2_4706
    );
    assert_eq!(
        hash_vector(zero_model.lu_solution().weights()),
        0x6b35_ba6e_4487_248e
    );
    assert_eq!(zero_model.evaluate_scalar(&witness).unwrap(), 3.0);
    assert_eq!(zero_model.evaluate_gradient(&witness).unwrap(), [0.0; 3]);

    let first_model = fit_single_surface_linear(&smoke_constraints(), &parameters()).unwrap();
    assert_eq!(first_model.evaluate_scalar(&witness).unwrap(), 0.75);
    assert_eq!(
        first_model.evaluate_gradient(&witness).unwrap(),
        [0.0, 0.0, 1.0]
    );

    let mut second = parameters();
    second.use_planar = false;
    second.polynomial_order = 2;
    let second_model = fit_single_surface_linear(&quadratic_constraints(), &second).unwrap();
    assert_eq!(
        hash_matrix(second_model.interpolation_matrix()),
        0x454b_3e95_77f7_904b
    );
    assert_eq!(
        hash_vector(second_model.right_hand_side()),
        0xd6da_1c34_9f96_ef50
    );
    assert!(second_model
        .lu_solution()
        .weights()
        .values()
        .iter()
        .all(|value| value.is_finite()));
    assert!(second_model.lu_solution().residual().accepted());
    assert_close(
        second_model.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0x403c_6fff_ffff_ffff),
        1.0e-9,
        1.0e-8,
    );
    let actual_gradient = second_model.evaluate_gradient(&witness).unwrap();
    let expected_gradient = [
        f64::from_bits(0x402a_8000_0000_0002),
        f64::from_bits(0x402e_ffff_ffff_ffff),
        f64::from_bits(0x4031_c000_0000_0000),
    ];
    for component in 0..3 {
        assert_close(
            actual_gradient[component],
            expected_gradient[component],
            1.0e-8,
            1.0e-7,
        );
    }
}

#[test]
fn linear_radius_value_path_and_derivative_error_match_frozen_surfe() {
    let mut parameters = parameters();
    parameters.use_planar = false;
    parameters.polynomial_order = 0;
    parameters.basis_type = RbfKernel::Linear;
    let model = fit_single_surface_linear(&constant_constraints(), &parameters).unwrap();
    assert_eq!(
        hash_matrix(model.interpolation_matrix()),
        0x64f8_be52_0f23_cc03
    );
    assert_eq!(
        hash_vector(model.lu_solution().weights()),
        0x6b35_ba6e_4487_248e
    );
    let witness = Point::new(0.25, 0.5, 0.75).unwrap();
    assert_eq!(model.evaluate_scalar(&witness).unwrap(), 3.0);
    assert!(matches!(
        model.evaluate_gradient(&witness),
        Err(SingleSurfaceLinearError::Evaluation(_))
    ));
}

#[test]
fn anisotropy_and_regression_smoothing_match_frozen_end_to_end_evidence() {
    let mut anisotropic_constraints = smoke_constraints();
    anisotropic_constraints
        .planars
        .push(Planar::from_normal(0.25, 0.75, 0.0, 0.0, 0.0, 1.0).unwrap());
    let mut anisotropic = parameters();
    anisotropic.model_global_anisotropy = true;
    let anisotropic_model =
        fit_single_surface_linear(&anisotropic_constraints, &anisotropic).unwrap();
    assert_eq!(
        hash_matrix(anisotropic_model.interpolation_matrix()),
        0xaaa5_ebef_44f6_0df7
    );
    assert_eq!(
        hash_vector(anisotropic_model.lu_solution().weights()),
        0x33af_162c_9d1e_f1d0
    );
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    assert_eq!(anisotropic_model.evaluate_scalar(&witness).unwrap(), 0.5);
    assert_eq!(
        anisotropic_model.evaluate_gradient(&witness).unwrap(),
        [0.0, 0.0, 1.0]
    );

    let constraints = Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, -1.0).unwrap(),
            Interface::new(2.0, 1.0, 3.0, 2.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        tangents: vec![georbf::Tangent::new(-1.0, 2.0, 1.0, 0.2, 0.7, 0.1).unwrap()],
        ..Constraints::default()
    };
    let mut smoothing = parameters();
    smoothing.use_tangent = true;
    smoothing.set_regression_smoothing(true, 0.75);
    let smoothing_model = fit_single_surface_linear(&constraints, &smoothing).unwrap();
    assert_eq!(
        hash_matrix(smoothing_model.interpolation_matrix()),
        0x99a1_5a49_5ec3_44c7
    );
    assert!(smoothing_model
        .lu_solution()
        .weights()
        .values()
        .iter()
        .all(|value| value.is_finite()));
    assert!(smoothing_model.lu_solution().residual().accepted());
    assert_close(
        smoothing_model.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0xbfe5_e308_80d3_ff4f),
        1.0e-9,
        1.0e-8,
    );
    let actual_gradient = smoothing_model.evaluate_gradient(&witness).unwrap();
    let expected_gradient = [
        f64::from_bits(0x3fde_7399_f232_4d80),
        f64::from_bits(0xbfc6_8f72_924e_78e2),
        f64::from_bits(0x3fe5_813c_5c39_90dd),
    ];
    for component in 0..3 {
        assert_close(
            actual_gradient[component],
            expected_gradient[component],
            1.0e-8,
            1.0e-7,
        );
    }
}

#[test]
fn vertical_slice_batch_evaluation_matches_repeated_single_point_evaluation() {
    let model = fit_single_surface_linear(&smoke_constraints(), &parameters()).unwrap();
    let points = [-1.0, 0.0, 0.5, 2.0]
        .into_iter()
        .map(|z| Point::new(0.25, 0.75, z).unwrap())
        .collect::<Vec<_>>();

    let scalars = model.evaluate_scalars(&points).unwrap();
    let gradients = model.evaluate_gradients(&points).unwrap();
    assert_eq!(scalars, vec![-1.0, 0.0, 0.5, 2.0]);
    assert_eq!(gradients, vec![[0.0, 0.0, 1.0]; 4]);
    for (index, point) in points.iter().enumerate() {
        assert_eq!(scalars[index], model.evaluate_scalar(point).unwrap());
        assert_eq!(gradients[index], model.evaluate_gradient(point).unwrap());
    }
}

#[test]
fn fitting_owns_the_cleaned_constraints_and_interface_grouping_evidence() {
    let mut constraints = smoke_constraints();
    constraints
        .interfaces
        .push(Interface::new(0.0, 0.0, 0.0, 99.0).unwrap());
    constraints
        .planars
        .push(Planar::from_normal(0.5, 0.5, 0.0, 1.0, 0.0, 0.0).unwrap());

    let model = fit_single_surface_linear(&constraints, &parameters()).unwrap();
    assert_eq!(model.collocation_removal().interfaces, 1);
    assert_eq!(model.collocation_removal().planars, 1);
    assert_eq!(model.constraints().interfaces.len(), 4);
    assert_eq!(model.constraints().planars.len(), 1);
    assert_eq!(model.constraints().interfaces[0].level(), 0.0);
    assert_eq!(model.interface_grouping().levels_descending(), &[0.0]);
    assert_eq!(model.interface_grouping().reference_indices(), &[0]);
    assert_eq!(
        model.interface_grouping().multi_point_groups(),
        &[vec![0, 1, 2, 3]]
    );
}

#[test]
fn t22_rejects_branches_owned_by_later_single_surface_tasks() {
    let mut with_inequality = smoke_constraints();
    with_inequality
        .inequalities
        .push(georbf::Inequality::new(0.5, 0.5, 1.0, 1.0).unwrap());
    assert!(matches!(
        fit_single_surface_linear(&with_inequality, &parameters()),
        Err(SingleSurfaceLinearError::InequalityBranchNotAvailable)
    ));

    let mut restricted = parameters();
    restricted.use_restricted_range = true;
    assert!(matches!(
        fit_single_surface_linear(&smoke_constraints(), &restricted),
        Err(SingleSurfaceLinearError::RestrictedRangeBranchNotAvailable)
    ));
}

#[test]
fn input_kernel_and_solver_failures_keep_their_frozen_stage() {
    let mut wrong_model = parameters();
    wrong_model.model_type = ModelType::LajaunieApproach;
    assert!(matches!(
        fit_single_surface_linear(&smoke_constraints(), &wrong_model),
        Err(SingleSurfaceLinearError::WrongModel)
    ));

    assert!(matches!(
        fit_single_surface_linear(&Constraints::default(), &parameters()),
        Err(SingleSurfaceLinearError::Surfe(Error::NoInterfaceData))
    ));
    let mut restricted_without_interfaces = parameters();
    restricted_without_interfaces.use_restricted_range = true;
    assert!(matches!(
        fit_single_surface_linear(&Constraints::default(), &restricted_without_interfaces),
        Err(SingleSurfaceLinearError::Surfe(Error::NoInterfaceData))
    ));

    let mut anisotropic = parameters();
    anisotropic.model_global_anisotropy = true;
    let no_planars = Constraints {
        interfaces: smoke_constraints().interfaces,
        ..Constraints::default()
    };
    assert!(matches!(
        fit_single_surface_linear(&no_planars, &anisotropic),
        Err(SingleSurfaceLinearError::Anisotropy(
            AnisotropyError::InsufficientPlanars
        ))
    ));

    let mut linear_radius = parameters();
    linear_radius.basis_type = RbfKernel::Linear;
    assert!(matches!(
        fit_single_surface_linear(&smoke_constraints(), &linear_radius),
        Err(SingleSurfaceLinearError::Assembly(_))
    ));
}
