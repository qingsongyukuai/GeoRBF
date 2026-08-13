use georbf::{
    constraint_layout, fit_stratigraphic, fit_stratigraphic_restricted,
    fit_stratigraphic_restricted_with_options, fit_stratigraphic_with_options, Constraints,
    DenseMatrix, DenseVector, DifferenceKind, Error, Inequality, Interface, LayoutDof,
    LayoutPointRef, LoqoOptions, LoqoSolveErrorKind, ModelType, Parameters, Planar, Point,
    QpOptions, QpSolveErrorKind, RbfKernel, ReconstructionStage, StratigraphicError,
    StratigraphicRestrictedError,
};

fn interface(x: f64, y: f64, z: f64, level: f64) -> Interface {
    Interface::new(x, y, z, level).unwrap()
}

fn constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            interface(0.0, 0.0, 0.0, 30.0),
            interface(2.0, 0.0, 0.0, 30.0),
            interface(0.0, 3.0, 0.0, 30.0),
            interface(0.0, 0.0, 4.0, 30.0),
            interface(5.0, 1.0, 1.0, 20.0),
            interface(7.0, -1.0, 2.0, 10.0),
        ],
        inequalities: vec![
            Inequality::new(3.0, -2.0, 2.0, 25.0).unwrap(),
            Inequality::new(4.0, 2.0, -1.0, 5.0).unwrap(),
        ],
        planars: vec![Planar::from_normal(1.0, 2.0, 3.0, 0.6, 0.0, 0.8).unwrap()],
        ..Constraints::default()
    }
}

fn parameters(restricted: bool) -> Parameters {
    Parameters {
        model_type: ModelType::StratigraphicHorizons,
        basis_type: RbfKernel::Cubic,
        polynomial_order: 1,
        shape_parameter: 2.0,
        min_stratigraphic_thickness: if restricted { 0.05 } else { 0.0 },
        interface_uncertainty: 0.2,
        angular_uncertainty: 8.0,
        use_restricted_range: restricted,
        ..Parameters::default()
    }
}

fn restricted_constraints() -> Constraints {
    let mut values = constraints();
    values.interfaces.truncate(4);
    values.inequalities.clear();
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

fn assert_bits_values(actual: &[f64], expected: &[u64], absolute: f64, relative: f64) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_close(*actual, f64::from_bits(*expected), absolute, relative);
    }
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
fn ordinary_and_restricted_paths_are_model_scoped() {
    let mut wrong = parameters(false);
    wrong.model_type = ModelType::SingleSurface;
    assert_eq!(
        fit_stratigraphic(&constraints(), &wrong).unwrap_err(),
        StratigraphicError::WrongModel
    );

    let mut wrong_restricted = wrong;
    wrong_restricted.use_restricted_range = true;
    assert_eq!(
        fit_stratigraphic_restricted(&constraints(), &wrong_restricted).unwrap_err(),
        StratigraphicRestrictedError::WrongModel
    );
    assert_eq!(
        fit_stratigraphic(&constraints(), &parameters(true)).unwrap_err(),
        StratigraphicError::RestrictedRangeBranchNotAvailable
    );
    assert_eq!(
        fit_stratigraphic_restricted(&constraints(), &parameters(false)).unwrap_err(),
        StratigraphicRestrictedError::RestrictedRangeRequired
    );
}

#[test]
fn boundary_lithology_points_use_only_strict_nearest_horizons_in_source_order() {
    let mut values = constraints();
    values.inequalities = vec![
        Inequality::new(2.0, 5.0, -3.0, 35.0).unwrap(),
        Inequality::new(3.0, -2.0, 2.0, 25.0).unwrap(),
        Inequality::new(4.0, 2.0, -1.0, 5.0).unwrap(),
    ];
    values.remove_collocated();
    let layout = constraint_layout(
        ModelType::StratigraphicHorizons,
        &values,
        &parameters(false),
    )
    .unwrap();
    assert_eq!(
        &layout.dofs()[..6],
        &[
            LayoutDof::Difference {
                kind: DifferenceKind::SequencedInterfaces,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(4),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::SequencedInterfaces,
                positive: LayoutPointRef::Interface(4),
                negative: LayoutPointRef::Interface(5),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityAboveLowerInterface,
                positive: LayoutPointRef::Inequality(0),
                negative: LayoutPointRef::Interface(0),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityBelowUpperInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Inequality(1),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityAboveLowerInterface,
                positive: LayoutPointRef::Inequality(1),
                negative: LayoutPointRef::Interface(4),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityBelowUpperInterface,
                positive: LayoutPointRef::Interface(5),
                negative: LayoutPointRef::Inequality(2),
            },
        ]
    );

    values.inequalities = vec![Inequality::new(9.0, 9.0, 9.0, 20.0).unwrap()];
    assert_eq!(
        constraint_layout(
            ModelType::StratigraphicHorizons,
            &values,
            &parameters(false),
        )
        .unwrap_err(),
        Error::InvalidInputData
    );
}

#[test]
fn three_level_ordinary_fit_matches_frozen_source_and_conversion_goldens() {
    let model = fit_stratigraphic(&constraints(), &parameters(false)).unwrap();
    assert_eq!(
        model.interface_grouping().levels_descending(),
        &[30.0, 20.0, 10.0]
    );
    assert_eq!(model.interface_grouping().reference_indices(), &[0, 4, 5]);
    assert_eq!(model.layout().partitions().inequality().len(), 5);
    assert_eq!(model.layout().partitions().equality().len(), 6);
    assert_eq!(model.layout().constraint_dof_count(), 11);
    assert_eq!(model.reconstruction().layout().constraint_dof_count(), 11);
    assert_eq!(model.reconstruction().layout().polynomial_dof_count(), 3);
    assert_eq!(model.reconstruction().mappings().len(), 11);

    assert_eq!(
        &model.layout().dofs()[..8],
        &[
            LayoutDof::Difference {
                kind: DifferenceKind::SequencedInterfaces,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(4),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::SequencedInterfaces,
                positive: LayoutPointRef::Interface(4),
                negative: LayoutPointRef::Interface(5),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityBelowUpperInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Inequality(0),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityAboveLowerInterface,
                positive: LayoutPointRef::Inequality(0),
                negative: LayoutPointRef::Interface(4),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::InequalityBelowUpperInterface,
                positive: LayoutPointRef::Interface(5),
                negative: LayoutPointRef::Inequality(1),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(1),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(2),
            },
            LayoutDof::Difference {
                kind: DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(0),
                negative: LayoutPointRef::Interface(3),
            },
        ]
    );

    assert_eq!(
        hash_matrix(model.modified_interpolation_matrix()),
        0xf96e_82f5_7bfd_f2f3
    );
    assert_eq!(
        hash_matrix(model.inequality_system().matrix()),
        0x5fed_4c46_742e_e72a
    );
    assert_eq!(
        hash_matrix(model.equality_system().matrix()),
        0x0bc4_cd2e_3e88_6bfd
    );
    assert_eq!(
        hash_vector(model.inequality_system().values()),
        0x6b1a_ca6e_4470_7806
    );
    assert_eq!(
        hash_vector(model.equality_system().values()),
        0x53ba_b84d_92de_6d90
    );
    assert!(model.qp_solution().residual().accepted());
    assert!(model.qp_solution().attempted());
    assert!(!model.qp_solution().trace().is_empty());

    assert_eq!(model.layer_relation_evidence().len(), 5);
    for relation in model.layer_relation_evidence() {
        assert!(relation.accepted());
        assert_close(
            relation.increment(),
            relation.matrix_value(),
            1.0e-10,
            1.0e-9,
        );
        assert_eq!(relation.maximum(), None);
    }
    assert_eq!(model.layer_relation_evidence()[0].minimum(), 0.0);
    assert_eq!(model.layer_relation_evidence()[1].minimum(), 0.0);
    assert_eq!(model.layer_relation_evidence()[2].minimum(), 0.0);

    assert_bits_values(
        model.source_interface_iso_values(),
        &[
            0x3ce9_4800_0000_0000,
            0xbdcd_4500_0000_0000,
            0xbec6_bdf1_af80_0000,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_eq!(model.source_interface_iso_value_evidence().len(), 3);
    for ((evidence, source_level), reference_index) in model
        .source_interface_iso_value_evidence()
        .iter()
        .zip([30.0, 20.0, 10.0])
        .zip([0, 4, 5])
    {
        assert_eq!(evidence.source_level(), source_level);
        assert_eq!(evidence.reference_index(), reference_index);
    }

    let reconstruction = model.reconstruction();
    assert!(reconstruction.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(reconstruction.interpolation_matrix()),
        0x6edb_fd0f_3806_992b
    );
    assert_bits_values(
        reconstruction.right_hand_side().values(),
        &[
            0x3dcd_4565_2000_0000,
            0x3ec6_bdd4_6a80_0000,
            0xbdf1_db62_dc00_0000,
            0x3df5_840f_8000_0000,
            0x3fc3_f930_9f24_4d90,
            0xbd04_5c00_0000_0000,
            0x3d11_a980_0000_0000,
            0x3cf1_a400_0000_0000,
            0x3fe3_3333_3333_359a,
            0xbd39_cf00_0000_0000,
            0x3fe9_9999_9999_9b15,
            0,
            0,
            0,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_bits_values(
        model.interface_iso_values(),
        &[
            0x3fda_9994_f5e0_7128,
            0x3fda_9994_f5d1_ce9a,
            0x3fda_9989_96e7_9970,
        ],
        1.0e-9,
        1.0e-8,
    );

    let point = Point::new(1.5, 0.5, 1.25).unwrap();
    let source_field = [
        model.evaluate_modified_scalar(&point).unwrap(),
        model.evaluate_modified_gradient(&point).unwrap()[0],
        model.evaluate_modified_gradient(&point).unwrap()[1],
        model.evaluate_modified_gradient(&point).unwrap()[2],
    ];
    let final_field = [
        model.evaluate_scalar(&point).unwrap(),
        model.evaluate_gradient(&point).unwrap()[0],
        model.evaluate_gradient(&point).unwrap()[1],
        model.evaluate_gradient(&point).unwrap()[2],
    ];
    assert_bits_values(
        &source_field,
        &[
            0xbfa2_c57d_6515_d91c,
            0x3fc5_d8dd_8d68_826e,
            0xbf82_ba9a_6aab_b7b8,
            0x3fb3_d6af_0c52_ace6,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_bits_values(
        &final_field,
        &[
            0x3fd8_40e5_493d_b279,
            0x3fc5_d8dd_8d68_823c,
            0xbf82_ba9a_6aab_c880,
            0x3fb3_d6af_0c52_a812,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_eq!(
        model.evaluate_scalars(&[point.clone()]).unwrap(),
        vec![final_field[0]]
    );
    assert_eq!(
        model.evaluate_gradients(&[point]).unwrap(),
        vec![[final_field[1], final_field[2], final_field[3]]]
    );
}

#[test]
fn restricted_fit_matches_all_bounds_loqo_conversion_and_final_field() {
    let model = fit_stratigraphic_restricted(&restricted_constraints(), &parameters(true)).unwrap();
    assert_eq!(model.layout().constraint_dof_count(), 6);
    assert_eq!(model.bound_evidence().len(), 6);
    assert!(model.layer_relation_evidence().is_empty());
    assert_eq!(model.bound_evidence()[0].lower(), -0.2);
    assert_eq!(model.bound_evidence()[0].range(), 0.4);
    assert_eq!(model.bound_evidence()[0].upper(), 0.2);
    assert_eq!(
        hash_matrix(model.modified_interpolation_matrix()),
        0x2ad2_3a37_599c_f2ec
    );
    assert_eq!(
        hash_vector(model.bounded_system().lower()),
        0x272d_6ebf_d118_5f8b
    );
    assert_eq!(
        hash_vector(model.bounded_system().range()),
        0xbdcb_d578_54cd_5a8c
    );
    assert!(model.loqo_solution().attempted());
    assert!(model.loqo_solution().residual().accepted());
    assert_eq!(model.loqo_solution().trace().len(), 11);
    assert_bits_values(
        model.source_interface_iso_values(),
        &[0xbf93_49f4_850a_61bc],
        1.0e-9,
        1.0e-8,
    );

    let reconstruction = model.reconstruction();
    assert_eq!(reconstruction.layout().constraint_dof_count(), 6);
    assert_eq!(reconstruction.layout().polynomial_dof_count(), 3);
    assert!(reconstruction.lu_solution().residual().accepted());
    assert_eq!(
        hash_matrix(reconstruction.interpolation_matrix()),
        0x5404_f3ec_faa0_b5e0
    );
    assert_bits_values(
        reconstruction.right_hand_side().values(),
        &[
            0xbfa2_10ab_41d0_f001,
            0xbf86_e05c_268b_3502,
            0xbf9d_964d_7d42_0c6e,
            0x3fe1_32cf_af44_ed5f,
            0x3fb7_5018_663e_6122,
            0x3fe8_32c5_0816_85cd,
            0,
            0,
            0,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_bits_values(
        model.interface_iso_values(),
        &[0x3fef_f5bc_50ef_b6c1],
        1.0e-9,
        1.0e-8,
    );

    let point = Point::new(1.5, 0.5, 1.25).unwrap();
    let source_field = [
        model.evaluate_modified_scalar(&point).unwrap(),
        model.evaluate_modified_gradient(&point).unwrap()[0],
        model.evaluate_modified_gradient(&point).unwrap()[1],
        model.evaluate_modified_gradient(&point).unwrap()[2],
    ];
    let final_field = [
        model.evaluate_scalar(&point).unwrap(),
        model.evaluate_gradient(&point).unwrap()[0],
        model.evaluate_gradient(&point).unwrap()[1],
        model.evaluate_gradient(&point).unwrap()[2],
    ];
    assert_bits_values(
        &source_field,
        &[
            0xbfa6_0666_28ae_ff7d,
            0x3fcd_1d03_4981_b598,
            0x3f8e_0ce4_aed6_eb34,
            0x3fb9_1c91_fb44_a5bd,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_bits_values(
        &final_field,
        &[
            0x3fef_2fa5_928d_19d8,
            0x3fcd_1d03_4981_b59b,
            0x3f8e_0ce4_aed6_eac0,
            0x3fb9_1c91_fb44_a5b8,
        ],
        1.0e-9,
        1.0e-8,
    );
    assert_eq!(
        model.evaluate_scalars(&[point.clone()]).unwrap(),
        vec![final_field[0]]
    );
    assert_eq!(
        model.evaluate_gradients(&[point]).unwrap(),
        vec![[final_field[1], final_field[2], final_field[3]]]
    );
}

#[test]
fn multilayer_restricted_frozen_candidate_is_preserved_as_a_safe_failure() {
    let error = fit_stratigraphic_restricted(&constraints(), &parameters(true)).unwrap_err();
    assert_eq!(error.stage(), Some(ReconstructionStage::Qp));
    let StratigraphicRestrictedError::Loqo(error) = error else {
        panic!("expected terminal LOQO evidence");
    };
    assert_eq!(error.kind(), LoqoSolveErrorKind::ResidualTooLarge);
    assert!(error.attempted());
    assert_eq!(error.trace().len(), 15);
    assert!(error.candidate_weights().is_some());
    let residual = error.residual().unwrap();
    assert!(residual.is_finite());
    assert!(residual.minimum_lower_slack() >= -residual.feasibility_limit());
    assert!(residual.minimum_upper_slack() >= -residual.feasibility_limit());
    assert!(!residual.accepted());
}

#[test]
fn invalid_inputs_iteration_caps_and_positive_separation_fail_with_stage_evidence() {
    let empty = Constraints::default();
    assert_eq!(
        fit_stratigraphic(&empty, &parameters(false)).unwrap_err(),
        StratigraphicError::Surfe(Error::NoInterfaceData)
    );

    let mut exact = constraints();
    exact.inequalities = vec![Inequality::new(9.0, 9.0, 9.0, 20.0).unwrap()];
    assert_eq!(
        fit_stratigraphic(&exact, &parameters(false)).unwrap_err(),
        StratigraphicError::Surfe(Error::InvalidInputData)
    );

    let error = fit_stratigraphic_with_options(
        &constraints(),
        &parameters(false),
        QpOptions { max_iterations: 0 },
    )
    .unwrap_err();
    assert_eq!(error.stage(), Some(ReconstructionStage::Qp));
    let StratigraphicError::Qp(error) = error else {
        panic!("expected ordinary QP evidence");
    };
    assert_eq!(error.kind(), QpSolveErrorKind::IterationLimit);
    assert!(error.attempted());

    let error = fit_stratigraphic_restricted_with_options(
        &restricted_constraints(),
        &parameters(true),
        LoqoOptions { max_iterations: 0 },
    )
    .unwrap_err();
    let StratigraphicRestrictedError::Loqo(error) = error else {
        panic!("expected restricted QP evidence");
    };
    assert_eq!(error.kind(), LoqoSolveErrorKind::IterationLimit);
    assert!(error.attempted());

    let mut separated = parameters(false);
    separated.min_stratigraphic_thickness = 0.5;
    let StratigraphicError::Qp(error) = fit_stratigraphic(&constraints(), &separated).unwrap_err()
    else {
        panic!("expected safe rejection of the frozen terminal candidate");
    };
    assert_eq!(error.kind(), QpSolveErrorKind::InfeasibleSolution);
    assert!(error.attempted());
    assert!(error.candidate_weights().is_some());
    assert!(error
        .residual()
        .is_some_and(|residual| !residual.accepted()));
}
