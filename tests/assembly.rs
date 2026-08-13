use georbf::{
    assemble_system, AnisotropicKernel, AssemblyConstraints, AssemblyError, ConstraintSystem,
    Constraints, DenseMatrix, DenseVector, Error, FunctionalKernel, Inequality, Interface,
    IsotropicKernel, KernelError, ModelType, ModifiedKernel, Parameters, Planar, RbfKernel,
    Tangent,
};

fn parameters(model: ModelType) -> Parameters {
    Parameters {
        model_type: model,
        basis_type: RbfKernel::Cubic,
        polynomial_order: 1,
        shape_parameter: 2.0,
        min_stratigraphic_thickness: 0.75,
        interface_uncertainty: 0.2,
        angular_uncertainty: 8.0,
        ..Parameters::default()
    }
}

fn interfaces() -> Vec<Interface> {
    vec![
        Interface::new(0.0, 0.0, 0.0, 30.0).unwrap(),
        Interface::new(2.0, 0.0, 0.0, 30.0).unwrap(),
        Interface::new(0.0, 3.0, 0.0, 30.0).unwrap(),
        Interface::new(0.0, 0.0, 4.0, 30.0).unwrap(),
        Interface::new(5.0, 1.0, 1.0, 20.0).unwrap(),
        Interface::new(7.0, -1.0, 2.0, 10.0).unwrap(),
    ]
}

fn planar(x: f64, y: f64, z: f64, normal: [f64; 3]) -> Planar {
    Planar::from_normal(x, y, z, normal[0], normal[1], normal[2]).unwrap()
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

fn base_constraints() -> Constraints {
    Constraints {
        interfaces: interfaces(),
        planars: vec![planar(1.0, 2.0, 3.0, [0.6, 0.0, 0.8])],
        tangents: vec![Tangent::new(-1.0, 2.0, 1.0, 0.2, 0.7, 0.1).unwrap()],
        ..Constraints::default()
    }
}

fn single_linear_constraints() -> Constraints {
    Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, -1.0).unwrap(),
            Interface::new(2.0, 1.0, 3.0, 2.0).unwrap(),
        ],
        planars: vec![planar(1.0, 2.0, 3.0, [0.6, 0.0, 0.8])],
        tangents: vec![Tangent::new(-1.0, 2.0, 1.0, 0.2, 0.7, 0.1).unwrap()],
        ..Constraints::default()
    }
}

fn single_quadratic_constraints() -> Constraints {
    let mut constraints = base_constraints();
    constraints.inequalities = vec![
        Inequality::new(3.0, -2.0, 2.0, 1.0).unwrap(),
        Inequality::new(4.0, 2.0, -1.0, -1.0).unwrap(),
    ];
    constraints
}

fn stratigraphic_constraints() -> Constraints {
    let mut constraints = base_constraints();
    constraints.inequalities = vec![
        Inequality::new(3.0, -2.0, 2.0, 25.0).unwrap(),
        Inequality::new(4.0, 2.0, -1.0, 5.0).unwrap(),
    ];
    constraints
}

fn modified_kernel(constraints: &Constraints) -> ModifiedKernel {
    let lists = interface_lists(constraints);
    ModifiedKernel::from_isotropic(IsotropicKernel::new(RbfKernel::Cubic, 2.0), &lists).unwrap()
}

fn interface_lists(constraints: &Constraints) -> Vec<Vec<Interface>> {
    let grouping = constraints.interface_grouping().unwrap();
    grouping
        .multi_point_groups()
        .iter()
        .map(|indices| {
            indices
                .iter()
                .map(|index| constraints.interfaces[*index].clone())
                .collect::<Vec<_>>()
        })
        .collect()
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

fn assert_symmetric(matrix: &DenseMatrix) {
    assert_eq!(matrix.rows(), matrix.cols());
    for row in 0..matrix.rows() {
        for column in 0..matrix.cols() {
            let left = matrix.get(row, column).unwrap();
            let right = matrix.get(column, row).unwrap();
            let scale = 1.0_f64.max(left.abs()).max(right.abs());
            assert!(
                (left - right).abs() <= 2.0e-12 * scale,
                "asymmetric ({row},{column}): {left} != {right}"
            );
        }
    }
}

fn quadratic_parts(constraints: &AssemblyConstraints) -> (&ConstraintSystem, &ConstraintSystem) {
    match constraints {
        AssemblyConstraints::Quadratic {
            equality,
            inequality,
        } => (equality, inequality),
        branch => panic!("expected ordinary quadratic branch, got {branch:?}"),
    }
}

#[test]
fn ordinary_five_model_matrices_and_rhs_match_frozen_full_matrix_hashes() {
    let cubic = IsotropicKernel::new(RbfKernel::Cubic, 2.0);

    let single_constraints = single_linear_constraints();
    let single = assemble_system(
        &single_constraints,
        &parameters(ModelType::SingleSurface),
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(single.interpolation_matrix().dimensions(), (10, 10));
    assert_eq!(
        hash_matrix(single.interpolation_matrix()),
        0xd85c_edb5_3723_977b
    );
    let single_rhs = single.constraints().linear_rhs().unwrap();
    assert_eq!(hash_vector(single_rhs), 0x9906_1e44_c270_4865);
    assert_eq!(single.interpolation_matrix().get(9, 9), Some(0.0));
    assert_symmetric(single.interpolation_matrix());
    // Complete first-order P/PT blocks retain `[x,y,z,1]`, derivative,
    // tangent-contraction, transpose, and zero polynomial suffix semantics.
    assert_eq!(
        [
            single.interpolation_matrix().get(6, 1).unwrap(),
            single.interpolation_matrix().get(7, 1).unwrap(),
            single.interpolation_matrix().get(8, 1).unwrap(),
            single.interpolation_matrix().get(9, 1).unwrap(),
        ],
        [2.0, 1.0, 3.0, 1.0]
    );
    assert_eq!(single.interpolation_matrix().get(6, 2), Some(1.0));
    assert_eq!(single.interpolation_matrix().get(7, 3), Some(1.0));
    assert_eq!(single.interpolation_matrix().get(8, 4), Some(1.0));
    assert_eq!(single.interpolation_matrix().get(6, 5), Some(0.2));
    assert_eq!(single.interpolation_matrix().get(7, 5), Some(0.7));
    assert_eq!(single.interpolation_matrix().get(8, 5), Some(0.1));

    let lajaunie_constraints = base_constraints();
    let lajaunie = assemble_system(
        &lajaunie_constraints,
        &parameters(ModelType::LajaunieApproach),
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(lajaunie.interpolation_matrix().dimensions(), (10, 10));
    assert_eq!(
        hash_matrix(lajaunie.interpolation_matrix()),
        0xcf1b_ba14_0abe_1850
    );
    assert_eq!(
        hash_vector(lajaunie.constraints().linear_rhs().unwrap()),
        0xdb39_7c52_5510_1064
    );
    assert_eq!(lajaunie.interpolation_matrix().get(0, 0), Some(-16.0));
    assert_symmetric(lajaunie.interpolation_matrix());
    assert_eq!(lajaunie.interpolation_matrix().get(7, 0), Some(-2.0));
    assert_eq!(lajaunie.interpolation_matrix().get(8, 0), Some(0.0));
    assert_eq!(lajaunie.interpolation_matrix().get(9, 0), Some(0.0));

    let continuous_constraints = Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, -1.0).unwrap(),
            Interface::new(2.0, 1.0, 3.0, 2.0).unwrap(),
            Interface::new(-1.0, 4.0, 2.0, 0.5).unwrap(),
        ],
        ..Constraints::default()
    };
    let continuous = assemble_system(
        &continuous_constraints,
        &parameters(ModelType::ContinuousProperty),
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(
        hash_matrix(continuous.interpolation_matrix()),
        0xe070_2c9a_0160_f9f7
    );
    assert_eq!(
        hash_vector(continuous.constraints().linear_rhs().unwrap()),
        0x2586_bdc1_1e58_2110
    );
    assert_symmetric(continuous.interpolation_matrix());

    let vector_constraints = Constraints {
        planars: vec![
            planar(1.0, 2.0, 3.0, [0.6, 0.0, 0.8]),
            planar(-2.0, 1.0, 0.5, [0.0, 0.8, 0.6]),
        ],
        ..Constraints::default()
    };
    let vector = assemble_system(
        &vector_constraints,
        &parameters(ModelType::VectorField),
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(
        hash_matrix(vector.interpolation_matrix()),
        0xbcb0_20bb_51aa_1d37
    );
    assert_eq!(
        hash_vector(vector.constraints().linear_rhs().unwrap()),
        0x33b7_90ee_7056_b051
    );
    assert_symmetric(vector.interpolation_matrix());

    let strat_constraints = stratigraphic_constraints();
    let strat_kernel = modified_kernel(&strat_constraints);
    let stratigraphic = assemble_system(
        &strat_constraints,
        &parameters(ModelType::StratigraphicHorizons),
        FunctionalKernel::from(&strat_kernel),
    )
    .unwrap();
    assert_eq!(
        hash_matrix(stratigraphic.interpolation_matrix()),
        0x8246_9ab1_9f37_3d0f
    );
    let (equality, inequality) = quadratic_parts(stratigraphic.constraints());
    assert_eq!(hash_matrix(equality.matrix()), 0x4bf2_609b_580f_fef9);
    assert_eq!(hash_vector(equality.values()), 0x5a16_d363_7084_4715);
    assert_eq!(hash_matrix(inequality.matrix()), 0x14d6_384c_8fd5_25eb);
    assert_eq!(hash_vector(inequality.values()), 0x4318_9583_5a08_35a6);
    assert_symmetric(stratigraphic.interpolation_matrix());
}

#[test]
fn modified_single_quadratic_signs_and_partition_slices_match_frozen_surfe() {
    let constraints = single_quadratic_constraints();
    let modified = modified_kernel(&constraints);
    let assembly = assemble_system(
        &constraints,
        &parameters(ModelType::SingleSurface),
        FunctionalKernel::from(&modified),
    )
    .unwrap();
    assert_eq!(
        hash_matrix(assembly.interpolation_matrix()),
        0x5e69_78a2_e0b2_58a4
    );
    let (equality, inequality) = quadratic_parts(assembly.constraints());
    assert_eq!(hash_matrix(equality.matrix()), 0xae71_d005_5165_d6ea);
    assert_eq!(hash_vector(equality.values()), 0x3428_af86_0e45_6914);
    assert_eq!(hash_matrix(inequality.matrix()), 0xa998_bf6c_29a4_660d);
    assert_eq!(hash_vector(inequality.values()), 0x8b03_8a41_009b_3de1);
    assert_symmetric(assembly.interpolation_matrix());
    for column in 0..assembly.interpolation_matrix().cols() {
        assert_eq!(
            inequality.matrix().get(0, column),
            assembly.interpolation_matrix().get(0, column)
        );
        assert_eq!(
            inequality.matrix().get(1, column),
            assembly
                .interpolation_matrix()
                .get(1, column)
                .map(|value| -value)
        );
    }
}

#[test]
fn restricted_range_lower_and_range_vectors_match_all_three_frozen_models() {
    let single_constraints = single_quadratic_constraints();
    let single_kernel = modified_kernel(&single_constraints);
    let mut single_parameters = parameters(ModelType::SingleSurface);
    single_parameters.use_restricted_range = true;
    let single = assemble_system(
        &single_constraints,
        &single_parameters,
        FunctionalKernel::from(&single_kernel),
    )
    .unwrap();
    let single_bounded = single.constraints().bounded().unwrap();
    assert_eq!(hash_matrix(single_bounded.matrix()), 0x5e69_78a2_e0b2_58a4);
    assert_eq!(hash_vector(single_bounded.lower()), 0xd5a2_dabe_a748_ef82);
    assert_eq!(hash_vector(single_bounded.range()), 0x6f71_4714_ffc2_3c20);

    let lajaunie_constraints = base_constraints();
    let lajaunie_kernel = modified_kernel(&lajaunie_constraints);
    let mut lajaunie_parameters = parameters(ModelType::LajaunieApproach);
    lajaunie_parameters.use_restricted_range = true;
    let lajaunie = assemble_system(
        &lajaunie_constraints,
        &lajaunie_parameters,
        FunctionalKernel::from(&lajaunie_kernel),
    )
    .unwrap();
    let lajaunie_bounded = lajaunie.constraints().bounded().unwrap();
    assert_eq!(
        hash_matrix(lajaunie_bounded.matrix()),
        0x305e_dc68_d5e3_1884
    );
    assert_eq!(hash_vector(lajaunie_bounded.lower()), 0x96c2_4168_3b58_c576);
    assert_eq!(hash_vector(lajaunie_bounded.range()), 0x7027_0be6_0e1e_0d1a);

    let strat_constraints = stratigraphic_constraints();
    let strat_kernel = modified_kernel(&strat_constraints);
    let mut strat_parameters = parameters(ModelType::StratigraphicHorizons);
    strat_parameters.use_restricted_range = true;
    let stratigraphic = assemble_system(
        &strat_constraints,
        &strat_parameters,
        FunctionalKernel::from(&strat_kernel),
    )
    .unwrap();
    let strat_bounded = stratigraphic.constraints().bounded().unwrap();
    assert_eq!(hash_matrix(strat_bounded.matrix()), 0x8246_9ab1_9f37_3d0f);
    assert_eq!(hash_vector(strat_bounded.lower()), 0x129a_5171_1516_2db9);
    assert_eq!(hash_vector(strat_bounded.range()), 0x9f09_c941_db1f_5678);
}

#[test]
fn ordinary_and_modified_anisotropic_assembly_match_frozen_full_matrix_hashes() {
    let vector_constraints = Constraints {
        planars: anisotropic_planars(),
        ..Constraints::default()
    };
    let vector_kernel =
        AnisotropicKernel::new(RbfKernel::Cubic, 2.0, &vector_constraints.planars).unwrap();
    let mut vector_parameters = parameters(ModelType::VectorField);
    vector_parameters.model_global_anisotropy = true;
    let vector = assemble_system(
        &vector_constraints,
        &vector_parameters,
        FunctionalKernel::from(&vector_kernel),
    )
    .unwrap();
    assert_eq!(
        hash_matrix(vector.interpolation_matrix()),
        0xe509_0944_f1c4_674f
    );
    assert_eq!(
        hash_vector(vector.constraints().linear_rhs().unwrap()),
        0xf17e_8a9b_20e5_fb6b
    );

    let mut single_constraints = single_quadratic_constraints();
    single_constraints.planars = anisotropic_planars();
    let radial =
        AnisotropicKernel::new(RbfKernel::Cubic, 2.0, &single_constraints.planars).unwrap();
    let modified =
        ModifiedKernel::from_anisotropic(radial, &interface_lists(&single_constraints)).unwrap();
    let mut single_parameters = parameters(ModelType::SingleSurface);
    single_parameters.model_global_anisotropy = true;
    let single = assemble_system(
        &single_constraints,
        &single_parameters,
        FunctionalKernel::from(&modified),
    )
    .unwrap();
    assert_eq!(single.interpolation_matrix().dimensions(), (24, 24));
    assert_eq!(
        hash_matrix(single.interpolation_matrix()),
        0x609f_48bc_3638_f808
    );
    let (equality, inequality) = quadratic_parts(single.constraints());
    assert_eq!(hash_matrix(equality.matrix()), 0xc6d4_dadc_661a_207c);
    assert_eq!(hash_vector(equality.values()), 0x85df_7582_a2c2_1f5e);
    assert_eq!(hash_matrix(inequality.matrix()), 0x9f4c_37d9_fa4a_2553);
    assert_eq!(hash_vector(inequality.values()), 0x8b03_8a41_009b_3de1);
    assert_symmetric(single.interpolation_matrix());
}

#[test]
fn smoothing_replaces_only_the_source_selected_diagonal_entries() {
    let cubic = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let mut single_parameters = parameters(ModelType::SingleSurface);
    single_parameters.set_regression_smoothing(false, 0.75);
    let single = assemble_system(
        &single_linear_constraints(),
        &single_parameters,
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(single.smoothing_value(), Some(0.421875));
    assert_eq!(
        hash_matrix(single.interpolation_matrix()),
        0x99a1_5a49_5ec3_44c7
    );
    assert_eq!(single.interpolation_matrix().get(0, 0), Some(0.421875));
    assert_eq!(single.interpolation_matrix().get(1, 1), Some(0.421875));
    assert_ne!(single.interpolation_matrix().get(2, 2), Some(0.421875));

    let mut lajaunie_parameters = parameters(ModelType::LajaunieApproach);
    lajaunie_parameters.set_regression_smoothing(false, 0.75);
    let lajaunie = assemble_system(
        &base_constraints(),
        &lajaunie_parameters,
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(
        hash_matrix(lajaunie.interpolation_matrix()),
        0x131c_91c3_ac6e_dfdd
    );
    assert_eq!(lajaunie.interpolation_matrix().get(0, 0), Some(0.421875));
    assert_eq!(lajaunie.interpolation_matrix().get(1, 1), Some(0.421875));
    assert_eq!(lajaunie.interpolation_matrix().get(2, 2), Some(0.421875));
    assert_ne!(lajaunie.interpolation_matrix().get(3, 3), Some(0.421875));

    let mut zero_parameters = parameters(ModelType::LajaunieApproach);
    zero_parameters.set_regression_smoothing(false, 0.0);
    let zero = assemble_system(
        &base_constraints(),
        &zero_parameters,
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(zero.smoothing_value(), Some(0.0));
    for index in 0..3 {
        assert_eq!(zero.interpolation_matrix().get(index, index), Some(0.0));
    }
}

#[test]
fn matrix_snapshot_and_assembly_failures_are_explicit() {
    let cubic = IsotropicKernel::new(RbfKernel::Cubic, 2.0);
    let continuous_constraints = Constraints {
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, -1.0).unwrap(),
            Interface::new(1.0, 0.0, 0.0, 2.0).unwrap(),
        ],
        ..Constraints::default()
    };
    let assembly = assemble_system(
        &continuous_constraints,
        &parameters(ModelType::ContinuousProperty),
        FunctionalKernel::from(&cubic),
    )
    .unwrap();
    assert_eq!(
        assembly.interpolation_matrix().debug_snapshot(),
        "2x2 row-major bits\n0x0000000000000000 0x3ff0000000000000\n\
0x3ff0000000000000 0x0000000000000000"
    );
    assert!(assembly.debug_snapshot().contains("branch=linear"));

    let mut unsupported = parameters(ModelType::SingleSurface);
    unsupported.polynomial_order = 3;
    assert_eq!(
        assemble_system(
            &single_linear_constraints(),
            &unsupported,
            FunctionalKernel::from(&cubic)
        ),
        Err(AssemblyError::Surfe(Error::InterpolationMatrixFailure))
    );

    let linear = IsotropicKernel::new(RbfKernel::Linear, 2.0);
    assert_eq!(
        assemble_system(
            &single_linear_constraints(),
            &parameters(ModelType::SingleSurface),
            FunctionalKernel::from(&linear)
        ),
        Err(AssemblyError::Kernel(
            KernelError::LinearDerivativeUnavailable
        ))
    );

    let modified = modified_kernel(&base_constraints());
    assert_eq!(
        assemble_system(
            &single_linear_constraints(),
            &parameters(ModelType::SingleSurface),
            FunctionalKernel::from(&modified)
        ),
        Err(AssemblyError::KernelLayoutMismatch)
    );
}
