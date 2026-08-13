use georbf::{
    fit_single_surface_restricted, fit_single_surface_restricted_with_options, Constraints,
    DenseMatrix, DenseVector, Inequality, Interface, LayoutDof, LoqoOptions, LoqoSolveErrorKind,
    ModelType, Parameters, Planar, Point, RbfKernel, ReconstructionStage,
    SingleSurfaceRestrictedError,
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
        tangents: Vec::new(),
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
        use_tangent: true,
        use_restricted_range: true,
        interface_uncertainty: 0.25,
        angular_uncertainty: 8.0,
        ..Parameters::default()
    }
}

fn basic_constraints() -> Constraints {
    let mut constraints = constraints();
    constraints.interfaces.truncate(4);
    constraints.inequalities.clear();
    constraints
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

fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    assert!(
        actual.is_finite() && delta <= absolute + relative * scale,
        "actual={actual:?} expected={expected:?} delta={delta:?}"
    );
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

#[test]
fn restricted_path_exposes_bounds_loqo_reconstruction_and_final_field() {
    let model = fit_single_surface_restricted(&basic_constraints(), &parameters()).unwrap();
    assert!(model.layout().internal_parameters().restricted_range);
    assert!(model.layout().internal_parameters().modified_basis);
    assert_eq!(model.layout().matrix_size(), 7);
    assert_eq!(model.bound_evidence().len(), 7);
    assert_eq!(model.bounded_system().lower().len(), 7);
    assert_eq!(model.bounded_system().range().len(), 7);
    assert_eq!(
        hash_matrix(model.modified_interpolation_matrix()),
        0xf814_3785_32f8_20b8
    );
    assert_eq!(
        hash_vector(model.bounded_system().lower()),
        0x0e82_94ea_6866_a5af
    );
    assert_eq!(
        hash_vector(model.bounded_system().range()),
        0xcf41_b808_7e30_4b08
    );
    assert!(model
        .bound_evidence()
        .iter()
        .all(|bound| bound.upper() == bound.lower() + bound.range()));
    assert!(matches!(
        model.bound_evidence()[0].dof(),
        LayoutDof::InterfaceValue { index: 0 }
    ));
    assert_eq!(model.bound_evidence()[0].lower(), -0.25);
    assert_eq!(model.bound_evidence()[0].upper(), 0.25);
    assert_eq!(model.bound_evidence()[1].range(), 0.5);
    assert!(model.loqo_solution().residual().accepted());
    assert_eq!(model.loqo_solution().trace().len(), 11);
    assert_close(
        model.loqo_solution().objective(),
        f64::from_bits(0x3fad_eeea_0809_2c2c),
        1.0e-10,
        1.0e-8,
    );

    let reconstruction = model.reconstruction();
    assert_eq!(reconstruction.mappings().len(), 7);
    assert_eq!(
        reconstruction
            .reconstructed_constraints()
            .inequalities
            .len(),
        0
    );
    assert_eq!(
        reconstruction.reconstructed_constraints().interfaces.len(),
        4
    );
    assert_eq!(reconstruction.layout().constraint_dof_count(), 7);
    assert_eq!(reconstruction.layout().polynomial_dof_count(), 4);
    assert!(reconstruction.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(reconstruction.interpolation_matrix()),
        0xef66_ae9b_c107_e6ff
    );

    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    let scalar = model.evaluate_scalar(&witness).unwrap();
    let gradient = model.evaluate_gradient(&witness).unwrap();
    assert_close(
        model.evaluate_modified_scalar(&witness).unwrap(),
        f64::from_bits(0xbfc3_da1f_b38b_e001),
        1.0e-10,
        1.0e-8,
    );
    for (actual, expected) in model
        .evaluate_modified_gradient(&witness)
        .unwrap()
        .into_iter()
        .zip([
            0x3faa_a5a0_0c4c_c82a,
            0xbfa5_362f_a47d_f96a,
            0xbfca_a335_20c9_ae36,
        ])
    {
        assert_close(actual, f64::from_bits(expected), 1.0e-9, 1.0e-8);
    }
    assert!(scalar.is_finite());
    assert!(gradient.into_iter().all(f64::is_finite));
    assert_close(
        scalar,
        f64::from_bits(0xbfc3_da1f_b38b_e008),
        1.0e-10,
        1.0e-8,
    );
    for (actual, expected) in gradient.into_iter().zip([
        0x3faa_a5a0_0c4c_c830,
        0xbfa5_362f_a47d_f958,
        0xbfca_a335_20c9_ae36,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-9, 1.0e-8);
    }
    assert_eq!(
        model.evaluate_scalars(&[witness.clone()]).unwrap(),
        vec![scalar]
    );
    assert_eq!(
        model.evaluate_gradients(&[witness]).unwrap(),
        vec![gradient]
    );
}

#[test]
fn anisotropic_restricted_path_uses_the_same_full_reconstruction() {
    let mut parameters = parameters();
    parameters.model_global_anisotropy = true;
    let model = fit_single_surface_restricted(&anisotropic_constraints(), &parameters).unwrap();
    assert_eq!(model.layout().matrix_size(), 23);
    assert_eq!(
        hash_matrix(model.modified_interpolation_matrix()),
        0x17e3_6493_b87b_ec03
    );
    assert_eq!(
        hash_vector(model.bounded_system().lower()),
        0x57e0_f4e0_8bf2_5864
    );
    assert_eq!(
        hash_vector(model.bounded_system().range()),
        0xc68e_6723_8bb3_62b6
    );
    assert!(model.loqo_solution().residual().accepted());
    assert_eq!(model.loqo_solution().trace().len(), 13);
    assert_close(
        model.loqo_solution().objective(),
        f64::from_bits(0x3fd0_c737_9aaa_4fc7),
        1.0e-10,
        1.0e-8,
    );
    assert!(model.reconstruction().lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(model.reconstruction().interpolation_matrix()),
        0x4018_0e3c_11ea_40cb
    );
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    assert_close(
        model.evaluate_modified_scalar(&witness).unwrap(),
        f64::from_bits(0x3fc1_531f_ad3e_8396),
        1.0e-9,
        1.0e-8,
    );
    for (actual, expected) in model
        .evaluate_modified_gradient(&witness)
        .unwrap()
        .into_iter()
        .zip([
            0x3fbf_4464_ba28_4abd,
            0x3fbc_3b14_aa27_ef60,
            0x3f8a_a336_f350_fa00,
        ])
    {
        assert_close(actual, f64::from_bits(expected), 1.0e-8, 1.0e-7);
    }
    assert_close(
        model.evaluate_scalar(&witness).unwrap(),
        f64::from_bits(0x3fc1_531f_ad3e_83b4),
        1.0e-9,
        1.0e-8,
    );
    for (actual, expected) in model.evaluate_gradient(&witness).unwrap().into_iter().zip([
        0x3fbf_4464_ba28_4aa2,
        0x3fbc_3b14_aa27_ef73,
        0x3f8a_a336_f350_f850,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-8, 1.0e-7);
    }
}

#[test]
fn restricted_path_keeps_failure_stages_distinct() {
    let mut not_restricted = parameters();
    not_restricted.use_restricted_range = false;
    assert_eq!(
        fit_single_surface_restricted(&constraints(), &not_restricted).unwrap_err(),
        SingleSurfaceRestrictedError::RestrictedRangeRequired
    );

    let error = fit_single_surface_restricted_with_options(
        &constraints(),
        &parameters(),
        LoqoOptions { max_iterations: 0 },
    )
    .unwrap_err();
    assert_eq!(error.stage(), Some(ReconstructionStage::Qp));
    assert!(matches!(
        error,
        SingleSurfaceRestrictedError::Loqo(ref source)
            if source.kind() == LoqoSolveErrorKind::IterationLimit
    ));

    let mut invalid_bounds = parameters();
    invalid_bounds.interface_uncertainty = -0.25;
    let error = fit_single_surface_restricted(&constraints(), &invalid_bounds).unwrap_err();
    assert_eq!(error.stage(), Some(ReconstructionStage::Qp));

    let mut unsupported_polynomial = parameters();
    unsupported_polynomial.polynomial_order = 3;
    let error =
        fit_single_surface_restricted(&basic_constraints(), &unsupported_polynomial).unwrap_err();
    assert_eq!(error.stage(), Some(ReconstructionStage::Reassembly));
    assert!(matches!(
        error,
        SingleSurfaceRestrictedError::Reconstruction(_)
    ));
}

#[test]
fn final_field_matches_reconstruction_witnesses() {
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    let model = fit_single_surface_restricted(&basic_constraints(), &parameters()).unwrap();
    let scalar = model.evaluate_scalar(&witness).unwrap();
    let gradient = model.evaluate_gradient(&witness).unwrap();
    let source = model
        .evaluate_modified_scalar(&witness)
        .expect("source Modified-Kernel evaluation");
    let source_gradient = model
        .evaluate_modified_gradient(&witness)
        .expect("source Modified-Kernel gradient evaluation");
    assert_close(source, scalar, 2.0e-8, 1.0e-7);
    for axis in 0..3 {
        assert_close(source_gradient[axis], gradient[axis], 2.0e-8, 1.0e-7);
    }
}
