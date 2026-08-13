use georbf::{
    fit_lajaunie_linear, fit_lajaunie_restricted, fit_lajaunie_restricted_with_options,
    Constraints, DenseMatrix, DenseVector, Inequality, Interface, LajaunieLinearError,
    LajaunieRestrictedError, LayoutDof, LayoutPointRef, LoqoOptions, LoqoSolveErrorKind, ModelType,
    Parameters, Planar, Point, RbfKernel, ReconstructionStage, Tangent,
};

fn parameters(restricted: bool) -> Parameters {
    Parameters {
        model_type: ModelType::LajaunieApproach,
        basis_type: RbfKernel::Cubic,
        polynomial_order: 1,
        shape_parameter: 2.0,
        use_interface: true,
        use_inequality: true,
        use_planar: true,
        use_tangent: true,
        use_restricted_range: restricted,
        interface_uncertainty: 0.25,
        angular_uncertainty: 8.0,
        ..Parameters::default()
    }
}

fn constraints() -> Constraints {
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
        inequalities: vec![Inequality::new(9.0, 8.0, 7.0, 30.0).unwrap()],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        tangents: vec![Tangent::new(-1.0, 2.0, 1.0, 0.2, 0.7, 0.1).unwrap()],
    }
}

fn restricted_constraints() -> Constraints {
    let mut values = constraints();
    values.interfaces.truncate(4);
    values.inequalities.clear();
    values.tangents.clear();
    values
}

fn five_planar_constraints() -> Constraints {
    let mut values = constraints();
    values.planars = vec![
        Planar::from_normal(1.0, 2.0, 3.0, 1.0, 0.0, 0.0).unwrap(),
        Planar::from_normal(-2.0, 1.0, 0.5, 0.0, 1.0, 0.0).unwrap(),
        Planar::from_normal(0.5, -1.0, 2.0, 0.0, 0.0, 1.0).unwrap(),
        Planar::from_normal(3.0, 0.5, -2.0, 0.36, -0.48, 0.8).unwrap(),
        Planar::from_normal(-1.5, -2.5, 1.0, -0.48, 0.64, 0.6).unwrap(),
    ];
    values
}

fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    assert!(
        actual.is_finite() && expected.is_finite() && delta <= absolute + relative * scale,
        "actual={actual:.17e} expected={expected:.17e} delta={delta:.3e}"
    );
}

fn mix_u64(hash: &mut u64, value: u64) {
    for byte in 0..8 {
        *hash ^= (value >> (8 * byte)) & 0xff;
        *hash = hash.wrapping_mul(1_099_511_628_211);
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

#[test]
fn ordinary_fit_preserves_reference_increment_and_iso_value_semantics() {
    let model = fit_lajaunie_linear(&constraints(), &parameters(false)).unwrap();
    assert_eq!(model.layout().constraint_dof_count(), 8);
    assert_eq!(model.layout().polynomial_dof_count(), 3);
    assert_eq!(model.layout().matrix_size(), 11);
    assert_eq!(
        model.interface_grouping().levels_descending(),
        &[30.0, 20.0, 10.0]
    );
    assert_eq!(model.interface_grouping().reference_indices(), &[0, 4, 6]);
    assert_eq!(
        &model.layout().dofs()[..4],
        &[
            LayoutDof::Difference {
                kind: georbf::DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(1),
            },
            LayoutDof::Difference {
                kind: georbf::DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(2),
            },
            LayoutDof::Difference {
                kind: georbf::DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(3),
            },
            LayoutDof::Difference {
                kind: georbf::DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(4),
                negative: LayoutPointRef::Interface(5),
            },
        ]
    );
    assert!(model.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(model.interpolation_matrix()),
        0xd5dc_0647_b570_080c
    );
    assert_eq!(hash_vector(model.right_hand_side()), 0x2d8a_0134_4965_3d59);
    assert_eq!(model.interface_iso_value_evidence().len(), 3);
    for ((evidence, source_level), reference_index) in model
        .interface_iso_value_evidence()
        .iter()
        .zip([30.0, 20.0, 10.0])
        .zip([0, 4, 6])
    {
        assert_eq!(evidence.source_level(), source_level);
        assert_eq!(evidence.reference_index(), reference_index);
        assert!(evidence.iso_value().is_finite());
    }
    assert_eq!(
        model.interface_iso_values(),
        model
            .interface_iso_value_evidence()
            .iter()
            .map(|value| value.iso_value())
            .collect::<Vec<_>>()
    );

    let points = [
        Point::new(0.25, 0.75, 0.5).unwrap(),
        Point::new(4.5, 1.25, 1.0).unwrap(),
    ];
    let scalars = model.evaluate_scalars(&points).unwrap();
    let gradients = model.evaluate_gradients(&points).unwrap();
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
    for (actual, expected) in model.interface_iso_values().iter().copied().zip([
        0x3fcb_2ebe_a777_9e50,
        0x3fe8_b7f8_2473_f158,
        0x3ff5_2c6e_42f7_9da4,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-9, 1.0e-8);
    }
    let field = [
        model.evaluate_scalar(&points[0]).unwrap(),
        model.evaluate_gradient(&points[0]).unwrap()[0],
        model.evaluate_gradient(&points[0]).unwrap()[1],
        model.evaluate_gradient(&points[0]).unwrap()[2],
    ];
    for (actual, expected) in field.into_iter().zip([
        0x3fb1_7133_d7c0_5e86,
        0x3fa9_cef5_3cb2_2c2b,
        0xbfa6_0093_8a9a_2f5c,
        0xbfce_52af_dd98_707a,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-9, 1.0e-8);
    }
}

#[test]
fn inequalities_are_cleaned_but_do_not_create_a_lajaunie_qp_branch() {
    let with_inequality = fit_lajaunie_linear(&constraints(), &parameters(false)).unwrap();
    let mut without = constraints();
    without.inequalities.clear();
    let without_inequality = fit_lajaunie_linear(&without, &parameters(false)).unwrap();

    assert_eq!(with_inequality.constraints().inequalities.len(), 1);
    assert_eq!(
        with_inequality.layout().internal_parameters().n_inequality,
        0
    );
    assert_eq!(
        with_inequality.interpolation_matrix().data(),
        without_inequality.interpolation_matrix().data()
    );
    assert_eq!(
        with_inequality.lu_solution().weights().values(),
        without_inequality.lu_solution().weights().values()
    );
    assert_eq!(
        with_inequality.interface_iso_values(),
        without_inequality.interface_iso_values()
    );
}

#[test]
fn restricted_fit_exposes_modified_loqo_reconstruction_and_both_iso_updates() {
    let model = fit_lajaunie_restricted(&restricted_constraints(), &parameters(true)).unwrap();
    assert!(model.layout().internal_parameters().restricted_range);
    assert!(model.layout().internal_parameters().modified_basis);
    assert_eq!(model.layout().constraint_dof_count(), 6);
    assert_eq!(model.bound_evidence().len(), 6);
    assert_eq!(model.bound_evidence()[0].lower(), -0.25);
    assert_eq!(model.bound_evidence()[0].upper(), 0.25);
    assert!(model.loqo_solution().residual().accepted());
    assert_eq!(
        hash_matrix(model.modified_interpolation_matrix()),
        0x2ad2_3a37_599c_f2ec
    );
    assert_eq!(
        hash_vector(model.bounded_system().lower()),
        0x7de9_47f2_0460_df13
    );
    assert_eq!(
        hash_vector(model.bounded_system().range()),
        0xcbb2_6e3a_a0b2_ff84
    );
    assert_eq!(model.loqo_solution().trace().len(), 11);
    assert_eq!(model.source_interface_iso_value_evidence().len(), 1);
    assert_close(
        model.source_interface_iso_values()[0],
        f64::from_bits(0xbf93_49f8_526b_4120),
        1.0e-9,
        1.0e-8,
    );

    let reconstruction = model.reconstruction();
    assert_eq!(reconstruction.layout().constraint_dof_count(), 6);
    assert_eq!(reconstruction.layout().polynomial_dof_count(), 3);
    assert!(reconstruction.lu_solution().residual().accepted());
    assert_eq!(reconstruction.mappings().len(), 6);
    assert_eq!(model.interface_iso_value_evidence().len(), 1);
    assert_eq!(
        hash_matrix(reconstruction.interpolation_matrix()),
        0x5404_f3ec_faa0_b5e0
    );
    assert_close(
        model.interface_iso_values()[0],
        f64::from_bits(0x3fef_f5bc_3bbb_0afa),
        1.0e-9,
        1.0e-8,
    );

    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    let source_scalar = model.evaluate_modified_scalar(&witness).unwrap();
    let final_scalar = model.evaluate_scalar(&witness).unwrap();
    assert!(source_scalar.is_finite());
    assert!(final_scalar.is_finite());
    for value in model.evaluate_modified_gradient(&witness).unwrap() {
        assert!(value.is_finite());
    }
    for value in model.evaluate_gradient(&witness).unwrap() {
        assert!(value.is_finite());
    }
    let source_field = [
        source_scalar,
        model.evaluate_modified_gradient(&witness).unwrap()[0],
        model.evaluate_modified_gradient(&witness).unwrap()[1],
        model.evaluate_modified_gradient(&witness).unwrap()[2],
    ];
    let final_field = [
        final_scalar,
        model.evaluate_gradient(&witness).unwrap()[0],
        model.evaluate_gradient(&witness).unwrap()[1],
        model.evaluate_gradient(&witness).unwrap()[2],
    ];
    for (actual, expected) in source_field.into_iter().zip([
        0xbfc3_da1f_8f9c_5ee9,
        0x3faa_a59f_4d27_03d5,
        0xbfa5_3630_4179_f652,
        0xbfca_a335_3c85_ffe5,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-9, 1.0e-8);
    }
    for (actual, expected) in final_field.into_iter().zip([
        0x3feb_9984_1a67_4d4b,
        0x3faa_a59f_4d27_03dc,
        0xbfa5_3630_4179_f65c,
        0xbfca_a335_3c85_ffe8,
    ]) {
        assert_close(actual, f64::from_bits(expected), 1.0e-9, 1.0e-8);
    }
}

#[test]
fn exact_levels_singletons_cleaning_anisotropy_and_smoothing_remain_deterministic() {
    let mut values = constraints();
    values.interfaces.push(values.interfaces[0].clone());
    let mut smooth = parameters(false);
    smooth.set_regression_smoothing(false, 0.75);
    let smoothed = fit_lajaunie_linear(&values, &smooth).unwrap();
    assert_eq!(smoothed.collocation_removal().interfaces, 1);
    assert_eq!(smoothed.interface_grouping().levels_descending().len(), 3);
    let smoothing = smoothed.assembled_system().smoothing_value().unwrap();
    for index in 0..4 {
        assert_eq!(
            smoothed.interpolation_matrix().get(index, index),
            Some(smoothing)
        );
    }
    assert_eq!(
        hash_matrix(smoothed.interpolation_matrix()),
        0x4c05_3d29_43c1_3d0d
    );

    let anisotropic_values = five_planar_constraints();
    let mut anisotropic = parameters(false);
    anisotropic.model_global_anisotropy = true;
    let model = fit_lajaunie_linear(&anisotropic_values, &anisotropic).unwrap();
    assert!(model.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(model.interpolation_matrix()),
        0xc940_bdbf_6714_c172
    );
    assert_eq!(hash_vector(model.right_hand_side()), 0xab1f_6e19_f912_456b);
    let witness = Point::new(0.5, -0.25, 1.5).unwrap();
    assert!(model.evaluate_scalar(&witness).unwrap().is_finite());
}

#[test]
fn ordinary_zero_and_second_order_truncated_polynomials_use_source_dimensions() {
    let mut zero = parameters(false);
    zero.polynomial_order = 0;
    zero.basis_type = RbfKernel::Gaussian;
    let zero_model = fit_lajaunie_linear(&constraints(), &zero).unwrap();
    assert_eq!(zero_model.layout().polynomial_dof_count(), 0);
    assert_eq!(zero_model.layout().matrix_size(), 8);
    assert!(zero_model.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(zero_model.interpolation_matrix()),
        0x8bfe_0624_4bbb_fdcf
    );
    assert_eq!(
        hash_vector(zero_model.right_hand_side()),
        0xa233_da71_99bc_8726
    );

    let mut second = parameters(false);
    second.polynomial_order = 2;
    let second_model = fit_lajaunie_linear(&five_planar_constraints(), &second).unwrap();
    assert_eq!(second_model.layout().polynomial_dof_count(), 9);
    assert_eq!(second_model.layout().matrix_size(), 29);
    assert!(second_model.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(second_model.interpolation_matrix()),
        0xa160_2b98_397e_73bc
    );
    assert_eq!(
        hash_vector(second_model.right_hand_side()),
        0x7d6f_2071_5a7b_2fbd
    );

    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    for model in [&zero_model, &second_model] {
        assert!(model.evaluate_scalar(&witness).unwrap().is_finite());
        assert!(model
            .evaluate_gradient(&witness)
            .unwrap()
            .into_iter()
            .all(f64::is_finite));
    }
}

#[test]
fn failures_preserve_model_input_basis_qp_reassembly_and_iteration_stages() {
    let mut wrong = parameters(false);
    wrong.model_type = ModelType::SingleSurface;
    assert!(matches!(
        fit_lajaunie_linear(&constraints(), &wrong),
        Err(LajaunieLinearError::WrongModel)
    ));

    let empty = Constraints::default();
    assert!(matches!(
        fit_lajaunie_linear(&empty, &parameters(false)),
        Err(LajaunieLinearError::Surfe(georbf::Error::NoInterfaceData))
    ));
    let singleton_levels = Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 2.0).unwrap(),
            Interface::new(1.0, 0.0, 0.0, 1.0).unwrap(),
        ],
        ..Constraints::default()
    };
    assert!(matches!(
        fit_lajaunie_linear(&singleton_levels, &parameters(false)),
        Err(LajaunieLinearError::Surfe(
            georbf::Error::NoInterfaceIncrementPairs
        ))
    ));
    assert!(matches!(
        fit_lajaunie_linear(&constraints(), &parameters(true)),
        Err(LajaunieLinearError::RestrictedRangeBranchNotAvailable)
    ));
    assert!(matches!(
        fit_lajaunie_restricted(&constraints(), &parameters(false)),
        Err(LajaunieRestrictedError::RestrictedRangeRequired)
    ));

    let iteration = fit_lajaunie_restricted_with_options(
        &restricted_constraints(),
        &parameters(true),
        LoqoOptions { max_iterations: 0 },
    )
    .unwrap_err();
    assert_eq!(iteration.stage(), Some(ReconstructionStage::Qp));
    assert!(matches!(
        iteration,
        LajaunieRestrictedError::Loqo(error)
            if error.kind() == LoqoSolveErrorKind::IterationLimit
    ));

    let mut second_order = parameters(true);
    second_order.polynomial_order = 2;
    let reconstruction =
        fit_lajaunie_restricted(&restricted_constraints(), &second_order).unwrap_err();
    assert_eq!(
        reconstruction.stage(),
        Some(ReconstructionStage::Reassembly)
    );
}

#[test]
fn iso_values_equal_direct_reference_point_predictions() {
    let model = fit_lajaunie_linear(&constraints(), &parameters(false)).unwrap();
    for evidence in model.interface_iso_value_evidence() {
        let point = model.constraints().interfaces[evidence.reference_index()].point();
        assert_close(
            evidence.iso_value(),
            model.evaluate_scalar(point).unwrap(),
            0.0,
            0.0,
        );
    }
}
