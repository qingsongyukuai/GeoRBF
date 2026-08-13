use georbf::{
    fit_single_surface_inequality, fit_single_surface_inequality_with_options, AssemblyError,
    Constraints, DenseMatrix, DenseVector, Error, Inequality, Interface, KernelError, ModelType,
    Parameters, Planar, Point, QpOptions, QpSolveErrorKind, RbfKernel,
    SingleSurfaceInequalityError, Tangent,
};

fn constraints() -> Constraints {
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
            Inequality::new(3.0, -2.0, 2.0, 1.0).unwrap(),
            Inequality::new(4.0, 2.0, -1.0, -1.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        ..Constraints::default()
    }
}

fn parameters() -> Parameters {
    Parameters {
        model_type: ModelType::SingleSurface,
        basis_type: RbfKernel::Cubic,
        shape_parameter: 2.0,
        polynomial_order: 1,
        use_interface: true,
        use_inequality: true,
        use_planar: true,
        ..Parameters::default()
    }
}

fn anisotropic_constraints() -> Constraints {
    let mut constraints = constraints();
    constraints.planars = vec![
        Planar::from_normal(1.0, 2.0, 3.0, 1.0, 0.0, 0.0).unwrap(),
        Planar::from_normal(-2.0, 1.0, 0.5, 0.0, 1.0, 0.0).unwrap(),
        Planar::from_normal(0.5, -1.0, 2.0, 0.0, 0.0, 1.0).unwrap(),
        Planar::from_normal(3.0, 0.5, -2.0, 0.36, -0.48, 0.8).unwrap(),
        Planar::from_normal(-1.5, -2.5, 1.0, -0.48, 0.64, 0.6).unwrap(),
    ];
    constraints
}

fn tangent_constraints() -> Constraints {
    let mut constraints = constraints();
    constraints.tangents = vec![Tangent::new(-1.0, 2.0, 1.0, 0.2, 0.7, 0.1).unwrap()];
    constraints
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
fn ordinary_single_surface_qp_exposes_solver_and_inequality_evidence() {
    let model = fit_single_surface_inequality(&constraints(), &parameters()).unwrap();
    assert_eq!(model.layout().matrix_size(), 11);
    assert_eq!(model.layout().internal_parameters().n_inequality, 2);
    assert_eq!(model.layout().internal_parameters().n_equality, 9);
    assert_eq!(
        hash_matrix(model.interpolation_matrix()),
        0x63ea_e652_a446_a5e0
    );
    assert_eq!(
        hash_matrix(model.equality_system().matrix()),
        0xadd3_d205_8def_61bf
    );
    assert_eq!(
        hash_vector(model.equality_system().values()),
        0xf86c_694c_5801_cd63
    );
    assert_eq!(
        hash_matrix(model.inequality_system().matrix()),
        0xd8bd_5b3f_c812_2beb
    );
    assert_eq!(
        hash_vector(model.inequality_system().values()),
        0x8b03_8a41_009b_3de1
    );
    // T04 treats weights as diagnostic-only. Lock the deterministic Rust
    // vector while using objective, feasibility, and fields for parity.
    assert_eq!(
        hash_vector(model.qp_solution().weights()),
        0x4ee4_0259_a92e_68fc
    );
    assert_eq!(model.qp_solution().trace().len(), 9);
    assert_close(
        model.qp_solution().objective(),
        f64::from_bits(0x40ac_3130_0c85_d85f),
        1.0e-10,
        1.0e-8,
    );
    assert_eq!(model.inequality_evidence().len(), 2);
    assert!(model.qp_solution().residual().accepted());
    assert!(model
        .inequality_evidence()
        .iter()
        .all(|evidence| evidence.transformed_value() >= -1.0e-8));

    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    assert_close(
        model.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0x4040_6850_5733_409f),
        1.0e-9,
        1.0e-8,
    );
    for (actual, expected) in model.evaluate_gradient(&witness).unwrap().into_iter().zip([
        0x3fe9_ecc8_119e_9754,
        0x3fe3_cc95_a65d_5675,
        0x400b_0675_0743_28a5,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-8, 1.0e-7);
    }
}

#[test]
fn inequality_signs_follow_the_frozen_positive_and_nonpositive_split() {
    let model = fit_single_surface_inequality(&constraints(), &parameters()).unwrap();
    assert_eq!(model.inequality_evidence()[0].row_sign(), 1.0);
    assert_eq!(model.inequality_evidence()[1].row_sign(), -1.0);
    assert_eq!(model.inequality_evidence()[0].source_level(), 1.0);
    assert_eq!(model.inequality_evidence()[1].source_level(), -1.0);
    assert_close(
        model.inequality_evidence()[0].scalar_field(),
        f64::from_bits(0x403f_d7fe_28e3_0fc1),
        1.0e-9,
        1.0e-8,
    );
    assert_close(
        model.inequality_evidence()[0].matrix_slack(),
        f64::from_bits(0x403f_d7fe_28e3_0fc0),
        1.0e-10,
        1.0e-8,
    );
    assert_close(
        model.inequality_evidence()[1].scalar_field(),
        f64::from_bits(0xbd87_fa00_0000_0000),
        1.0e-9,
        1.0e-8,
    );
    assert_close(
        model.inequality_evidence()[1].matrix_slack(),
        f64::from_bits(0x3d88_1e00_0000_0000),
        1.0e-10,
        1.0e-8,
    );
    assert!(!model.inequality_evidence()[0].active_within_solver_tolerance());
    assert!(model.inequality_evidence()[1].active_within_solver_tolerance());

    let mut zero_is_nonpositive = constraints();
    zero_is_nonpositive.inequalities[1] = Inequality::new(4.0, 2.0, -1.0, 0.0).unwrap();
    let zero_is_nonpositive =
        fit_single_surface_inequality(&zero_is_nonpositive, &parameters()).unwrap();
    assert_eq!(
        zero_is_nonpositive.inequality_evidence()[1].row_sign(),
        -1.0
    );
}

#[test]
fn anisotropic_and_smoothing_qp_branches_match_frozen_model_goldens() {
    let mut anisotropic_parameters = parameters();
    anisotropic_parameters.model_global_anisotropy = true;
    let anisotropic =
        fit_single_surface_inequality(&anisotropic_constraints(), &anisotropic_parameters).unwrap();
    assert_eq!(anisotropic.layout().matrix_size(), 23);
    assert_eq!(
        hash_matrix(anisotropic.interpolation_matrix()),
        0x17e3_6493_b87b_ec03
    );
    assert_eq!(
        hash_matrix(anisotropic.equality_system().matrix()),
        0x4797_40fa_1b41_3f07
    );
    assert_eq!(
        hash_vector(anisotropic.equality_system().values()),
        0xaf6f_ceaa_4a71_8ee5
    );
    assert_eq!(
        hash_matrix(anisotropic.inequality_system().matrix()),
        0x9725_9701_20fa_2ec8
    );
    assert_eq!(
        hash_vector(anisotropic.qp_solution().weights()),
        0xa2a3_b38e_0483_d0cd
    );
    assert_eq!(anisotropic.qp_solution().trace().len(), 9);
    assert_close(
        anisotropic.qp_solution().objective(),
        f64::from_bits(0x40ac_4302_0da9_478c),
        1.0e-10,
        1.0e-8,
    );
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    assert_close(
        anisotropic.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0x403f_8eb3_7122_9a62),
        1.0e-9,
        1.0e-8,
    );
    for (actual, expected) in anisotropic
        .evaluate_gradient(&witness)
        .unwrap()
        .into_iter()
        .zip([
            0xbfd8_6314_bfb9_c418,
            0xbff5_b5f0_feab_6b40,
            0x4011_8909_da1a_073c,
        ])
    {
        assert_close(actual, f64::from_bits(expected), 1.0e-8, 1.0e-7);
    }

    let mut smoothing_parameters = parameters();
    smoothing_parameters.use_regression_smoothing = true;
    smoothing_parameters.smoothing_amount = 0.75;
    let smoothing = fit_single_surface_inequality(&constraints(), &smoothing_parameters).unwrap();
    assert_eq!(
        hash_matrix(smoothing.interpolation_matrix()),
        0xbecd_d5d0_87d7_7fb4
    );
    assert_eq!(
        hash_matrix(smoothing.equality_system().matrix()),
        0xc143_7fbc_c584_114e
    );
    assert_eq!(
        hash_matrix(smoothing.inequality_system().matrix()),
        0xdd6c_c21d_1d62_fe02
    );
    assert_eq!(
        hash_vector(smoothing.qp_solution().weights()),
        0x5123_ae32_54fc_9a7d
    );
    assert_eq!(smoothing.qp_solution().trace().len(), 12);
    assert_close(
        smoothing.qp_solution().objective(),
        f64::from_bits(0x40b1_5447_4619_df36),
        1.0e-10,
        1.0e-8,
    );
    assert_close(
        smoothing.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0x4043_28c0_f4ee_4de6),
        1.0e-9,
        1.0e-8,
    );
    for (actual, expected) in smoothing
        .evaluate_gradient(&witness)
        .unwrap()
        .into_iter()
        .zip([
            0x4010_8736_70e4_ba32,
            0x3fc4_23d2_03f1_eced,
            0xbfab_7cf1_9d9a_24b0,
        ])
    {
        assert_close(actual, f64::from_bits(expected), 1.0e-8, 1.0e-7);
    }
}

#[test]
fn tangent_equality_remains_in_the_qp_layout_and_field_sum() {
    let model = fit_single_surface_inequality(&tangent_constraints(), &parameters()).unwrap();
    assert_eq!(model.layout().matrix_size(), 12);
    assert_eq!(model.layout().internal_parameters().n_equality, 10);
    assert_eq!(
        hash_matrix(model.interpolation_matrix()),
        0xf21a_39c5_5811_8048
    );
    assert_eq!(
        hash_matrix(model.equality_system().matrix()),
        0x4041_0411_ecce_6d7e
    );
    assert_eq!(
        hash_vector(model.equality_system().values()),
        0x3428_af86_0e45_6914
    );
    assert_eq!(
        hash_matrix(model.inequality_system().matrix()),
        0x8f8e_7533_a3de_45bd
    );
    assert_close(
        model.qp_solution().objective(),
        f64::from_bits(0x40ac_3206_b735_5940),
        1.0e-10,
        1.0e-8,
    );
    assert_eq!(model.qp_solution().trace().len(), 9);
    assert!(model.qp_solution().residual().accepted());

    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    assert_close(
        model.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0x4040_7297_bc5b_f3ee),
        1.0e-9,
        1.0e-8,
    );
    for (actual, expected) in model.evaluate_gradient(&witness).unwrap().into_iter().zip([
        0x3fd4_08fa_854c_1d98,
        0x3fe9_89a6_86e1_b998,
        0x400a_4a2a_58f9_c94d,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-8, 1.0e-7);
    }
}

#[test]
fn cleaning_and_batch_evaluation_preserve_frozen_order() {
    let baseline = fit_single_surface_inequality(&constraints(), &parameters()).unwrap();
    let mut duplicated = constraints();
    duplicated
        .inequalities
        .push(Inequality::new(3.0, -2.0, 2.0, -1.0).unwrap());
    duplicated
        .interfaces
        .push(Interface::new(0.0, 0.0, 0.0, 99.0).unwrap());
    let cleaned = fit_single_surface_inequality(&duplicated, &parameters()).unwrap();
    assert_eq!(cleaned.collocation_removal().inequalities, 1);
    assert_eq!(cleaned.collocation_removal().interfaces, 1);
    assert_eq!(
        hash_matrix(cleaned.interpolation_matrix()),
        hash_matrix(baseline.interpolation_matrix())
    );
    assert_eq!(
        hash_vector(cleaned.qp_solution().weights()),
        hash_vector(baseline.qp_solution().weights())
    );

    let points = vec![
        Point::new(0.25, 0.75, 0.5).unwrap(),
        Point::new(3.0, -2.0, 2.0).unwrap(),
        Point::new(4.0, 2.0, -1.0).unwrap(),
    ];
    assert_eq!(
        cleaned.evaluate_scalars(&points).unwrap(),
        points
            .iter()
            .map(|point| cleaned.evaluate_scalar(point).unwrap())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        cleaned.evaluate_gradients(&points).unwrap(),
        points
            .iter()
            .map(|point| cleaned.evaluate_gradient(point).unwrap())
            .collect::<Vec<_>>()
    );
}

#[test]
fn basis_assembly_and_qp_failures_keep_their_model_stage() {
    let mut too_few = constraints();
    too_few.interfaces.truncate(3);
    assert_eq!(
        fit_single_surface_inequality(&too_few, &parameters()).unwrap_err(),
        SingleSurfaceInequalityError::Basis(Error::ModifiedKernelCreationFailure)
    );

    let mut linear_kernel = parameters();
    linear_kernel.basis_type = RbfKernel::Linear;
    assert_eq!(
        fit_single_surface_inequality(&constraints(), &linear_kernel).unwrap_err(),
        SingleSurfaceInequalityError::Assembly(AssemblyError::Kernel(
            KernelError::LinearDerivativeUnavailable
        ))
    );

    let error = fit_single_surface_inequality_with_options(
        &constraints(),
        &parameters(),
        QpOptions { max_iterations: 0 },
    )
    .unwrap_err();
    let SingleSurfaceInequalityError::Qp(error) = error else {
        panic!("expected model-level QP error")
    };
    assert_eq!(error.kind(), QpSolveErrorKind::IterationLimit);
    assert!(error.attempted());
    assert_eq!(error.surfe_error(), Error::PredictorCorrectorSolverFailure);

    let mut infeasible = constraints();
    infeasible.inequalities = vec![Inequality::new(0.0, 0.0, 0.0, -1.0).unwrap()];
    let infeasible_error = fit_single_surface_inequality_with_options(
        &infeasible,
        &parameters(),
        QpOptions { max_iterations: 8 },
    )
    .unwrap_err();
    let SingleSurfaceInequalityError::Qp(infeasible_error) = infeasible_error else {
        panic!("expected infeasible model-level QP error")
    };
    assert_eq!(infeasible_error.kind(), QpSolveErrorKind::KktSolveFailure);
    assert!(infeasible_error.attempted());
    assert_eq!(
        infeasible_error.surfe_error(),
        Error::PredictorCorrectorSolverFailure
    );
}

#[test]
fn specialized_path_rejects_wrong_branch_before_solving() {
    let mut wrong_model = parameters();
    wrong_model.model_type = ModelType::LajaunieApproach;
    assert_eq!(
        fit_single_surface_inequality(&constraints(), &wrong_model).unwrap_err(),
        SingleSurfaceInequalityError::WrongModel
    );

    let mut no_inequality = constraints();
    no_inequality.inequalities.clear();
    assert_eq!(
        fit_single_surface_inequality(&no_inequality, &parameters()).unwrap_err(),
        SingleSurfaceInequalityError::NoInequalities
    );

    let mut restricted = parameters();
    restricted.use_restricted_range = true;
    assert_eq!(
        fit_single_surface_inequality(&constraints(), &restricted).unwrap_err(),
        SingleSurfaceInequalityError::RestrictedRangeBranchNotAvailable
    );

    let mut no_interface = constraints();
    no_interface.interfaces.clear();
    assert_eq!(
        fit_single_surface_inequality(&no_interface, &parameters()).unwrap_err(),
        SingleSurfaceInequalityError::Surfe(Error::NoInterfaceData)
    );
}
