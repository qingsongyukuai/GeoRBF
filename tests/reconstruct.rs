use georbf::{
    assemble_system, reconstruct_from_qp_weights, solve_and_reconstruct, ConstraintLayout,
    Constraints, DenseVector, FunctionalKernel, Inequality, Interface, IsotropicKernel, LayoutDof,
    ModelType, ModifiedKernel, Parameters, Planar, Point, RbfKernel, ReconstructionError,
    ReconstructionStage,
};

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let scale = 1.0_f64.max(actual.abs()).max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn assert_tolerance(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && delta <= absolute + relative * scale,
        "actual={actual:.17e}, expected={expected:.17e}, delta={delta:.3e}, absolute={absolute:.3e}, relative={relative:.3e}"
    );
}

fn assert_bits_values(actual: &[f64], expected_bits: &[u64], absolute: f64, relative: f64) {
    assert_eq!(actual.len(), expected_bits.len());
    for (actual, bits) in actual.iter().zip(expected_bits) {
        assert_tolerance(*actual, f64::from_bits(*bits), absolute, relative);
    }
}

fn vector_from_bits(bits: &[u64]) -> DenseVector {
    DenseVector::from_values(bits.iter().copied().map(f64::from_bits).collect())
}

fn assert_witness_golden(
    witness: &georbf::ReconstructionPredictionWitness,
    source_bits: [u64; 4],
    reconstructed_bits: [u64; 4],
) {
    let source = witness.source_gradient();
    assert_tolerance(
        witness.source_scalar(),
        f64::from_bits(source_bits[0]),
        1.0e-9,
        1.0e-8,
    );
    for axis in 0..3 {
        assert_tolerance(
            source[axis],
            f64::from_bits(source_bits[axis + 1]),
            1.0e-8,
            1.0e-7,
        );
    }
    let reconstructed = witness.reconstructed_gradient();
    assert_tolerance(
        witness.reconstructed_scalar(),
        f64::from_bits(reconstructed_bits[0]),
        1.0e-9,
        1.0e-8,
    );
    for axis in 0..3 {
        assert_tolerance(
            reconstructed[axis],
            f64::from_bits(reconstructed_bits[axis + 1]),
            1.0e-8,
            1.0e-7,
        );
    }
}

fn mix_u64(hash: &mut u64, value: u64) {
    for byte in 0..8 {
        *hash ^= (value >> (8 * byte)) & 0xff;
        *hash = hash.wrapping_mul(1_099_511_628_211);
    }
}

fn hash_matrix(matrix: &georbf::DenseMatrix) -> u64 {
    let mut hash = 1_469_598_103_934_665_603;
    mix_u64(&mut hash, matrix.rows() as u64);
    mix_u64(&mut hash, matrix.cols() as u64);
    for value in matrix.data() {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn single_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 0.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 0.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 0.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 0.0).unwrap(),
        ],
        inequalities: vec![Inequality::new(1.0, 1.0, 2.0, 1.0).unwrap()],
        planars: vec![Planar::from_normal(0.5, 0.5, 1.0, 0.0, 0.0, 1.0).unwrap()],
        ..Constraints::default()
    }
}

fn single_parameters() -> Parameters {
    Parameters {
        model_type: ModelType::SingleSurface,
        basis_type: RbfKernel::Cubic,
        polynomial_order: 1,
        shape_parameter: 2.0,
        ..Parameters::default()
    }
}

fn lajaunie_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 20.0).unwrap(),
            Interface::new(2.0, 0.0, 0.0, 20.0).unwrap(),
            Interface::new(0.0, 3.0, 0.0, 20.0).unwrap(),
            Interface::new(0.0, 0.0, 4.0, 20.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(0.75, 0.5, 1.25, 0.0, 0.0, 1.0).unwrap()],
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
        inequalities: vec![Inequality::new(-2.0, 1.0, 3.0, 35.0).unwrap()],
        planars: vec![Planar::from_normal(1.0, 1.5, 2.0, 0.0, 0.0, 1.0).unwrap()],
        ..Constraints::default()
    }
}

fn assert_linear_layout(layout: &ConstraintLayout, constraints: usize, polynomial: usize) {
    assert!(!layout.internal_parameters().modified_basis);
    assert!(!layout.internal_parameters().restricted_range);
    assert_eq!(layout.constraint_dof_count(), constraints);
    assert_eq!(layout.polynomial_dof_count(), polynomial);
    assert_eq!(layout.partitions().equality().len(), constraints);
    assert!(layout.partitions().inequality().is_empty());
    assert!(layout.partitions().bounded().is_empty());
}

#[test]
fn single_qp_weights_reconstruct_every_constraint_then_resolve_with_lu() {
    let constraints = single_constraints();
    let parameters = single_parameters();
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces.clone()]).unwrap();
    let qp_system =
        assemble_system(&constraints, &parameters, FunctionalKernel::from(&modified)).unwrap();
    let weights = DenseVector::from_values(vec![0.125, -0.25, 0.5, 0.75, -0.375, 0.2, -0.1, 0.3]);
    assert_eq!(weights.len(), qp_system.layout().constraint_dof_count());
    let witnesses = [
        Point::new(0.25, 0.5, 0.75).unwrap(),
        Point::new(1.25, -0.5, 2.0).unwrap(),
    ];

    let result = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &qp_system,
        &weights,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        &witnesses,
    )
    .unwrap();
    // Frozen Single Surface prepends every inequality as a reconstructed
    // interface, then appends every original interface. It does not filter by
    // an active-set tolerance.
    assert_eq!(result.reconstructed_constraints().inequalities.len(), 0);
    assert_eq!(result.reconstructed_constraints().interfaces.len(), 5);
    assert_linear_layout(result.layout(), 8, 4);
    assert_eq!(result.mappings().len(), 8);
    assert_eq!(result.mappings()[0].source_index(), 0);
    assert_eq!(result.mappings()[0].reconstructed_index(), 0);
    assert!(matches!(
        result.mappings()[0].source_dof(),
        LayoutDof::InequalityValue { index: 0 }
    ));
    assert!(matches!(
        result.mappings()[0].reconstructed_dof(),
        LayoutDof::InterfaceValue { index: 0 }
    ));
    assert!(result.lu_solution().residual().accepted());

    for mapping in result.mappings() {
        assert_close(
            result
                .right_hand_side()
                .get(mapping.reconstructed_index())
                .unwrap(),
            mapping.target_value(),
            0.0,
        );
    }
    assert_eq!(result.prediction_witnesses().len(), witnesses.len());
    for witness in result.prediction_witnesses() {
        assert_close(
            witness.source_scalar(),
            witness.reconstructed_scalar(),
            2.0e-10,
        );
        for axis in 0..3 {
            assert_close(
                witness.source_gradient()[axis],
                witness.reconstructed_gradient()[axis],
                2.0e-10,
            );
        }
    }
}

#[test]
fn lajaunie_reconstruction_retains_pairs_and_forces_three_polynomial_terms() {
    let constraints = lajaunie_constraints();
    let parameters = Parameters {
        model_type: ModelType::LajaunieApproach,
        use_restricted_range: true,
        interface_uncertainty: 0.25,
        angular_uncertainty: 8.0,
        polynomial_order: 1,
        ..Parameters::default()
    };
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces.clone()]).unwrap();
    let source =
        assemble_system(&constraints, &parameters, FunctionalKernel::from(&modified)).unwrap();
    assert_eq!(source.layout().constraint_dof_count(), 6);
    let weights = DenseVector::from_values(vec![0.4, -0.3, 0.2, 0.1, -0.05, 0.025]);
    let witnesses = [
        Point::new(0.4, -0.25, 1.5).unwrap(),
        Point::new(1.5, 0.75, -0.5).unwrap(),
    ];
    let result = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &source,
        &weights,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        &witnesses,
    )
    .unwrap();

    assert_linear_layout(result.layout(), 6, 3);
    assert_eq!(result.layout().dofs()[..6], source.layout().dofs()[..6]);
    assert_eq!(result.mappings().len(), 6);
    assert!(result.lu_solution().residual().accepted());
    let scalar_offsets = result
        .prediction_witnesses()
        .iter()
        .map(|prediction| prediction.source_scalar() - prediction.reconstructed_scalar())
        .collect::<Vec<_>>();
    assert_close(scalar_offsets[0], scalar_offsets[1], 2.0e-10);
    for prediction in result.prediction_witnesses() {
        for axis in 0..3 {
            assert_close(
                prediction.source_gradient()[axis],
                prediction.reconstructed_gradient()[axis],
                2.0e-10,
            );
        }
    }
}

#[test]
fn stratigraphic_reconstruction_converts_all_pair_kinds_to_equalities() {
    let constraints = stratigraphic_constraints();
    let parameters = Parameters {
        model_type: ModelType::StratigraphicHorizons,
        min_stratigraphic_thickness: 0.5,
        polynomial_order: 2,
        ..Parameters::default()
    };
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces[..4].to_vec()]).unwrap();
    let source =
        assemble_system(&constraints, &parameters, FunctionalKernel::from(&modified)).unwrap();
    assert_eq!(source.layout().constraint_dof_count(), 9);
    let weights = DenseVector::from_values(vec![0.1, -0.2, 0.3, 0.4, -0.1, 0.05, 0.2, -0.15, 0.25]);
    let witnesses = [
        Point::new(0.75, 0.25, 1.0).unwrap(),
        Point::new(4.0, -0.5, 2.5).unwrap(),
    ];
    let result = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &source,
        &weights,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        &witnesses,
    )
    .unwrap();

    assert_linear_layout(result.layout(), 9, 3);
    assert_eq!(result.layout().internal_parameters().n_inequality, 3);
    assert_eq!(result.layout().dofs()[..9], source.layout().dofs()[..9]);
    assert_eq!(result.mappings().len(), 9);
    assert!(result.lu_solution().residual().accepted());
    let scalar_offsets = result
        .prediction_witnesses()
        .iter()
        .map(|prediction| prediction.source_scalar() - prediction.reconstructed_scalar())
        .collect::<Vec<_>>();
    assert_close(scalar_offsets[0], scalar_offsets[1], 5.0e-10);
    for prediction in result.prediction_witnesses() {
        for axis in 0..3 {
            assert_close(
                prediction.source_gradient()[axis],
                prediction.reconstructed_gradient()[axis],
                5.0e-10,
            );
        }
    }
}

#[test]
fn reconstruction_failures_preserve_source_qp_reassembly_and_lu_stages() {
    let constraints = single_constraints();
    let parameters = single_parameters();
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces.clone()]).unwrap();
    let source =
        assemble_system(&constraints, &parameters, FunctionalKernel::from(&modified)).unwrap();
    let weights = DenseVector::from_values(vec![0.0; source.layout().constraint_dof_count()]);
    let witness = Point::new(0.25, 0.5, 0.75).unwrap();

    let wrong_weights = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &source,
        &DenseVector::zeros(weights.len() - 1),
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        std::slice::from_ref(&witness),
    )
    .unwrap_err();
    assert_eq!(wrong_weights.stage(), ReconstructionStage::Reassembly);

    let modified_target = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &source,
        &weights,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&modified),
        std::slice::from_ref(&witness),
    )
    .unwrap_err();
    assert_eq!(modified_target.stage(), ReconstructionStage::Reassembly);

    let nonfinite_ordinary = IsotropicKernel::new(RbfKernel::Gaussian, f64::NAN);
    let lu = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &source,
        &weights,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&nonfinite_ordinary),
        std::slice::from_ref(&witness),
    )
    .unwrap_err();
    assert_eq!(lu.stage(), ReconstructionStage::Lu);

    let source_assembly = solve_and_reconstruct(
        &constraints,
        &parameters,
        FunctionalKernel::from(&ordinary),
        FunctionalKernel::from(&ordinary),
        std::slice::from_ref(&witness),
    )
    .unwrap_err();
    assert!(matches!(
        source_assembly,
        ReconstructionError::SourceAssembly(_)
    ));
    assert_eq!(source_assembly.stage(), ReconstructionStage::SourceAssembly);

    let lajaunie_constraints = lajaunie_constraints();
    let invalid_bounds = Parameters {
        model_type: ModelType::LajaunieApproach,
        use_restricted_range: true,
        interface_uncertainty: -0.25,
        angular_uncertainty: 8.0,
        polynomial_order: 1,
        ..Parameters::default()
    };
    let lajaunie_modified =
        ModifiedKernel::from_isotropic(ordinary, &[lajaunie_constraints.interfaces.clone()])
            .unwrap();
    let qp = solve_and_reconstruct(
        &lajaunie_constraints,
        &invalid_bounds,
        FunctionalKernel::from(&lajaunie_modified),
        FunctionalKernel::from(&ordinary),
        &[],
    )
    .unwrap_err();
    assert_eq!(qp.stage(), ReconstructionStage::Qp);
}

#[test]
fn high_level_single_path_keeps_the_qp_solution_as_evidence() {
    let constraints = single_constraints();
    let parameters = single_parameters();
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces.clone()]).unwrap();
    let result = solve_and_reconstruct(
        &constraints,
        &parameters,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        &[Point::new(0.25, 0.5, 0.75).unwrap()],
    )
    .unwrap();
    let georbf::ReconstructionSourceSolution::PredictorCorrector(source) =
        result.source_solution().unwrap()
    else {
        panic!("Single ordinary inequalities must use predictor-corrector");
    };
    assert!(source.slack().get(0).unwrap() > 1.0);
    assert!(source.dual_inequality().get(0).unwrap().abs() < 1.0e-10);
    assert!(result.lu_solution().residual().accepted());
}

#[test]
fn single_conversion_keeps_active_and_inactive_inequalities_in_source_order() {
    let mut constraints = single_constraints();
    constraints
        .inequalities
        .push(Inequality::new(1.5, -0.5, 0.25, -1.0).unwrap());
    constraints
        .inequalities
        .push(Inequality::new(-0.5, 1.5, 0.25, -1.0).unwrap());
    let parameters = single_parameters();
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces.clone()]).unwrap();
    let result = solve_and_reconstruct(
        &constraints,
        &parameters,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        &[],
    )
    .unwrap();
    let georbf::ReconstructionSourceSolution::PredictorCorrector(source) =
        result.source_solution().unwrap()
    else {
        panic!("Single ordinary inequalities must use predictor-corrector");
    };
    assert!(source.slack().get(0).unwrap() > 1.0);
    assert!(source.dual_inequality().get(0).unwrap().abs() < 1.0e-10);
    for index in 1..3 {
        assert!(source.slack().get(index).unwrap().abs() < 1.0e-8);
        assert!(source.dual_inequality().get(index).unwrap() > 0.05);
    }
    assert_eq!(result.mappings().len(), 10);
    for index in 0..3 {
        assert!(matches!(
            result.mappings()[index].source_dof(),
            LayoutDof::InequalityValue { index: mapped } if *mapped == index
        ));
        assert!(matches!(
            result.mappings()[index].reconstructed_dof(),
            LayoutDof::InterfaceValue { index: mapped } if *mapped == index
        ));
    }
}

#[test]
fn lajaunie_non_first_order_records_the_frozen_reassembly_failure() {
    let constraints = lajaunie_constraints();
    let parameters = Parameters {
        model_type: ModelType::LajaunieApproach,
        use_restricted_range: true,
        interface_uncertainty: 0.25,
        angular_uncertainty: 8.0,
        polynomial_order: 2,
        ..Parameters::default()
    };
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let modified =
        ModifiedKernel::from_isotropic(ordinary, &[constraints.interfaces.clone()]).unwrap();
    let source =
        assemble_system(&constraints, &parameters, FunctionalKernel::from(&modified)).unwrap();
    let weights = DenseVector::from_values(vec![0.1; source.layout().constraint_dof_count()]);
    let error = reconstruct_from_qp_weights(
        &constraints,
        &parameters,
        &source,
        &weights,
        FunctionalKernel::from(&modified),
        FunctionalKernel::from(&ordinary),
        &[],
    )
    .unwrap_err();
    assert_eq!(error.stage(), ReconstructionStage::Reassembly);
}

#[test]
fn frozen_cpp_reconstruction_matrix_rhs_lu_and_field_golden() {
    let witnesses = [
        Point::new(0.25, 0.5, 0.75).unwrap(),
        Point::new(1.25, -0.5, 2.0).unwrap(),
    ];
    let ordinary = IsotropicKernel::new(RbfKernel::Cubic, 2.0);

    let single_constraints = single_constraints();
    let single_parameters = single_parameters();
    let single_modified =
        ModifiedKernel::from_isotropic(ordinary, &[single_constraints.interfaces.clone()]).unwrap();
    let single_source = assemble_system(
        &single_constraints,
        &single_parameters,
        FunctionalKernel::from(&single_modified),
    )
    .unwrap();
    let single_source_weights = vector_from_bits(&[
        0x3cb1_5df6_0000_0000,
        0x3f99_b0a4_88ed_909e,
        0x3f74_123f_9f8b_a303,
        0x3f66_ca39_2688_e86c,
        0xbfa0_c73d_cad0_cb59,
        0xbf84_123f_9f8b_a38f,
        0xbf81_17aa_dce6_aed3,
        0x3fc0_c73d_cad0_cb48,
    ]);
    let single = reconstruct_from_qp_weights(
        &single_constraints,
        &single_parameters,
        &single_source,
        &single_source_weights,
        FunctionalKernel::from(&single_modified),
        FunctionalKernel::from(&ordinary),
        &witnesses,
    )
    .unwrap();
    assert_eq!(
        hash_matrix(single.interpolation_matrix()),
        0xa06d_69b3_396a_05a3
    );
    assert_bits_values(
        single.right_hand_side().values(),
        &[
            0x3ff4_d310_ae6d_925a,
            0xbc50_0000_0000_0000,
            0xbc30_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x3c64_0000_0000_0000,
            0x3c64_0000_0000_0000,
            0x3fef_ffff_ffff_ffff,
            0,
            0,
            0,
            0,
        ],
        1.0e-11,
        1.0e-10,
    );
    assert_witness_golden(
        &single.prediction_witnesses()[0],
        [
            0x3fe2_3701_5616_2ae3,
            0xbf9b_f4f8_4366_fdea,
            0xbf83_69d1_0598_831f,
            0x3fec_77f6_aa54_02a8,
        ],
        [
            0x3fe2_3701_5616_2ae2,
            0xbf9b_f4f8_4366_fdb0,
            0xbf83_69d1_0598_82f0,
            0x3fec_77f6_aa54_02a7,
        ],
    );
    assert_witness_golden(
        &single.prediction_witnesses()[1],
        [
            0x3ff0_b8de_726f_8787,
            0xbfab_1ba0_655c_e20b,
            0x3fd4_c937_c709_4357,
            0x3fa6_2b6e_f18b_948e,
        ],
        [
            0x3ff0_b8de_726f_8787,
            0xbfab_1ba0_655c_e20c,
            0x3fd4_c937_c709_435d,
            0x3fa6_2b6e_f18b_949c,
        ],
    );

    let lajaunie_constraints = lajaunie_constraints();
    let lajaunie_parameters = Parameters {
        model_type: ModelType::LajaunieApproach,
        use_restricted_range: true,
        min_stratigraphic_thickness: 0.5,
        interface_uncertainty: 0.25,
        angular_uncertainty: 8.0,
        polynomial_order: 1,
        shape_parameter: 2.0,
        ..Parameters::default()
    };
    let lajaunie_modified =
        ModifiedKernel::from_isotropic(ordinary, &[lajaunie_constraints.interfaces.clone()])
            .unwrap();
    let lajaunie_source = assemble_system(
        &lajaunie_constraints,
        &lajaunie_parameters,
        FunctionalKernel::from(&lajaunie_modified),
    )
    .unwrap();
    let lajaunie_source_weights = vector_from_bits(&[
        0x3e60_94a0_8ba2_90e0,
        0x3e5e_42ee_b8e4_0080,
        0x3e71_02cf_d158_ce40,
        0xbf84_28b6_df2d_922b,
        0xbf81_38c3_d451_c23b,
        0x3fc2_a83f_bedd_1cba,
    ]);
    let lajaunie = reconstruct_from_qp_weights(
        &lajaunie_constraints,
        &lajaunie_parameters,
        &lajaunie_source,
        &lajaunie_source_weights,
        FunctionalKernel::from(&lajaunie_modified),
        FunctionalKernel::from(&ordinary),
        &witnesses,
    )
    .unwrap();
    assert_eq!(
        hash_matrix(lajaunie.interpolation_matrix()),
        0x02b8_296f_5a23_584c
    );
    assert_bits_values(
        lajaunie.right_hand_side().values(),
        &[
            0xbf98_5d4e_760a_cd30,
            0xbf9a_88b1_0e2d_b740,
            0xbfb0_adfe_5f5c_2213,
            0x3fb1_aea4_384e_f462,
            0x3f83_dd9b_44d7_9061,
            0x3fef_ec0b_7170_feef,
            0,
            0,
            0,
        ],
        1.0e-11,
        1.0e-10,
    );
    assert_witness_golden(
        &lajaunie.prediction_witnesses()[0],
        [
            0x3fd9_b14a_d7ed_5953,
            0xbf90_91e1_ff5e_cb43,
            0xbf8c_892c_8430_3177,
            0x3fe6_dbd4_2d58_5155,
        ],
        [
            0xbfef_f0c9_0147_6fde,
            0xbf90_91e1_ff5e_cae0,
            0xbf8c_892c_8430_3240,
            0x3fe6_dbd4_2d58_5154,
        ],
    );
    assert_witness_golden(
        &lajaunie.prediction_witnesses()[1],
        [
            0x3ff1_8bc1_e25e_835a,
            0x3fb6_3416_9668_ce66,
            0x3fd4_95bf_9000_35a5,
            0x3fc9_6229_054d_43f2,
        ],
        [
            0xbfd3_63d5_5102_2b96,
            0x3fb6_3416_9668_ce90,
            0x3fd4_95bf_9000_3595,
            0x3fc9_6229_054d_43f7,
        ],
    );

    let stratigraphic_constraints = stratigraphic_constraints();
    let stratigraphic_parameters = Parameters {
        model_type: ModelType::StratigraphicHorizons,
        min_stratigraphic_thickness: 0.5,
        interface_uncertainty: 0.25,
        angular_uncertainty: 8.0,
        polynomial_order: 1,
        shape_parameter: 2.0,
        ..Parameters::default()
    };
    let stratigraphic_modified = ModifiedKernel::from_isotropic(
        ordinary,
        &[stratigraphic_constraints.interfaces[..4].to_vec()],
    )
    .unwrap();
    let stratigraphic_source = assemble_system(
        &stratigraphic_constraints,
        &stratigraphic_parameters,
        FunctionalKernel::from(&stratigraphic_modified),
    )
    .unwrap();
    let stratigraphic_source_weights = vector_from_bits(&[
        0x3f77_028d_663e_c5eb,
        0x3f4f_4919_77ef_0968,
        0x3f4f_3939_8706_fca9,
        0xbf94_a2fa_1866_d564,
        0xbf88_e2e7_2f2f_8940,
        0x3fa1_e7f8_e083_a3aa,
        0xbf81_353c_821f_a5f2,
        0xbfa1_4109_6877_6270,
        0x3fc2_61aa_b898_73ee,
    ]);
    let stratigraphic = reconstruct_from_qp_weights(
        &stratigraphic_constraints,
        &stratigraphic_parameters,
        &stratigraphic_source,
        &stratigraphic_source_weights,
        FunctionalKernel::from(&stratigraphic_modified),
        FunctionalKernel::from(&ordinary),
        &witnesses,
    )
    .unwrap();
    assert_eq!(
        hash_matrix(stratigraphic.interpolation_matrix()),
        0x7918_b629_f414_a98c
    );
    assert_bits_values(
        stratigraphic.right_hand_side().values(),
        &[
            0x3fe0_0000_0000_0004,
            0x3fdf_ffff_ffff_ffe8,
            0x3da3_ea53_8000_0000,
            0xbc6c_0000_0000_0000,
            0xbc58_0000_0000_0000,
            0xbc50_0000_0000_0000,
            0xbc98_0000_0000_0000,
            0x3c90_0000_0000_0000,
            0x3ff0_0000_0000_0001,
            0,
            0,
            0,
        ],
        1.0e-11,
        1.0e-10,
    );
    assert_witness_golden(
        &stratigraphic.prediction_witnesses()[0],
        [
            0xbf9b_9122_3f27_9028,
            0xbf80_6eda_fd33_98fe,
            0xbfb7_b283_3227_6b89,
            0x3fc0_2d7c_0a02_efd2,
        ],
        [
            0xbfed_4936_8ee7_4617,
            0xbf80_6eda_fd33_9910,
            0xbfb7_b283_3227_6b90,
            0x3fc0_2d7c_0a02_efd3,
        ],
    );
    assert_witness_golden(
        &stratigraphic.prediction_witnesses()[1],
        [
            0x3fcc_ff2f_f6fe_f0f6,
            0x3fc1_1bbf_ccc2_e96c,
            0x3fd0_fbc6_0a00_89e9,
            0x3fbc_1305_dfda_1762,
        ],
        [
            0xbfe5_2ce1_7f2e_4d56,
            0x3fc1_1bbf_ccc2_e973,
            0x3fd0_fbc6_0a00_89e3,
            0x3fbc_1305_dfda_1786,
        ],
    );
}
