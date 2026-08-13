use georbf::{
    Builder, Constraints, DenseMatrix, FittedModel, GreedyHookBody, GreedyStopReason, Inequality,
    Interface, ModelType, Planar, Point, Tangent, GREEDY_MODEL_AUDIT,
};

fn single_surface_constraints() -> Constraints {
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

fn constraints(model: ModelType) -> Constraints {
    match model {
        ModelType::SingleSurface => single_surface_constraints(),
        ModelType::LajaunieApproach => lajaunie_constraints(),
        ModelType::StratigraphicHorizons => stratigraphic_constraints(),
        ModelType::ContinuousProperty => continuous_constraints(),
        ModelType::VectorField => vector_constraints(),
    }
}

fn interpolation_matrix(model: &FittedModel) -> &DenseMatrix {
    match model {
        FittedModel::SingleSurfaceLinear(model) => model.interpolation_matrix(),
        FittedModel::LajaunieLinear(model) => model.interpolation_matrix(),
        FittedModel::Stratigraphic(model) => model.modified_interpolation_matrix(),
        FittedModel::ContinuousProperty(model) => model.interpolation_matrix(),
        FittedModel::VectorField(model) => model.interpolation_matrix(),
        unexpected => panic!("unexpected public branch: {unexpected:?}"),
    }
}

fn matrix_bits(model: &FittedModel) -> Vec<u64> {
    interpolation_matrix(model)
        .data()
        .iter()
        .map(|value| value.to_bits())
        .collect()
}

fn mix_u64(hash: &mut u64, value: u64) {
    for byte in value.to_le_bytes() {
        *hash ^= u64::from(byte);
        *hash = hash.wrapping_mul(1_099_511_628_211);
    }
}

fn matrix_hash(model: &FittedModel) -> u64 {
    let matrix = interpolation_matrix(model);
    let mut hash = 1_469_598_103_934_665_603;
    mix_u64(&mut hash, matrix.rows() as u64);
    mix_u64(&mut hash, matrix.cols() as u64);
    for value in matrix.data() {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn assert_close(actual: f64, expected_bits: u64, tolerance: f64) {
    let expected = f64::from_bits(expected_bits);
    let scale = 1.0_f64.max(actual.abs()).max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn frozen_greedy_hook_table_is_complete_and_honest() {
    assert_eq!(GREEDY_MODEL_AUDIT.len(), ModelType::ALL.len());
    assert_eq!(GREEDY_MODEL_AUDIT.map(|entry| entry.model), ModelType::ALL);
    assert!(GREEDY_MODEL_AUDIT
        .iter()
        .all(|entry| !entry.reachable_from_surfe_api));

    assert_eq!(GREEDY_MODEL_AUDIT[0].minimal, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[0].residual, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[0].append, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[1].minimal, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[1].residual, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[1].append, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[2].minimal, GreedyHookBody::TodoStub);
    assert_eq!(GREEDY_MODEL_AUDIT[2].residual, GreedyHookBody::TodoStub);
    assert_eq!(GREEDY_MODEL_AUDIT[2].append, GreedyHookBody::TodoStub);
    assert_eq!(GREEDY_MODEL_AUDIT[3].minimal, GreedyHookBody::TodoStub);
    assert_eq!(GREEDY_MODEL_AUDIT[3].residual, GreedyHookBody::Implemented);
    assert_eq!(GREEDY_MODEL_AUDIT[3].append, GreedyHookBody::Implemented);
    assert_eq!(
        GREEDY_MODEL_AUDIT[4].greedy_factory_model,
        ModelType::ContinuousProperty
    );
    assert_eq!(GREEDY_MODEL_AUDIT[4].minimal, GreedyHookBody::TodoStub);
    assert_eq!(GREEDY_MODEL_AUDIT[4].residual, GreedyHookBody::TodoStub);
    assert_eq!(GREEDY_MODEL_AUDIT[4].append, GreedyHookBody::TodoStub);
}

#[test]
fn public_greedy_request_runs_zero_rounds_for_all_five_models() {
    let witness = Point::new(0.25, 0.75, 0.5).unwrap();
    let frozen = [
        (
            0x908a_282d_5b21_e173,
            0x3fe0_0000_0000_0000,
            [0, 0, 0x3ff0_0000_0000_0000],
        ),
        (
            0xd5dc_0647_b570_080c,
            0x3fb1_7133_d7c0_5e86,
            [
                0x3fa9_cef5_3cb2_2c2b,
                0xbfa6_0093_8a9a_2f5c,
                0xbfce_52af_dd98_707a,
            ],
        ),
        (
            0xf96e_82f5_7bfd_f2f3,
            0xbfc2_25a5_e3c2_1329,
            [
                0x3fba_f566_40fd_4ebd,
                0xbfa9_238f_043f_e56e,
                0xbfd0_454c_51d3_9c86,
            ],
        ),
        (
            0x6253_1c8a_e9a7_fb27,
            0x3ff2_1391_a86d_8a8d,
            [
                0xbfe7_f9ec_d70b_dd0f,
                0x3fd9_8d71_533f_4e3e,
                0x3fba_5e42_1840_a4e0,
            ],
        ),
        (
            0xeb24_aafe_28b2_dd7b,
            0xbfc2_6789_da3f_9380,
            [
                0x3fcd_628d_517a_9eaf,
                0x3fca_ca92_78ec_a8f8,
                0x3fe6_0baa_7247_f0f4,
            ],
        ),
    ];

    for (model_type, (expected_matrix, expected_scalar, expected_gradient)) in
        ModelType::ALL.into_iter().zip(frozen)
    {
        let baseline = Builder::new(model_type)
            .set_constraints(constraints(model_type))
            .fit()
            .unwrap();
        let mut requested_builder = Builder::new(model_type);
        requested_builder
            .set_constraints(constraints(model_type))
            .set_greedy_algorithm(false, 0.125, 4.5);
        let requested = requested_builder.fit().unwrap();

        let trace = requested.greedy_trace();
        assert!(trace.stored_use_greedy);
        assert_eq!(trace.interface_uncertainty.to_bits(), 0x3fc0_0000_0000_0000);
        assert_eq!(trace.angular_uncertainty.to_bits(), 0x4012_0000_0000_0000);
        assert!(trace.rounds.is_empty());
        assert_eq!(trace.stop_reason, GreedyStopReason::NotCalledBySurfeApi);

        assert_eq!(matrix_bits(&requested), matrix_bits(&baseline));
        assert_eq!(matrix_hash(&requested), expected_matrix);
        assert_eq!(
            requested.evaluate_scalar(&witness).unwrap().to_bits(),
            baseline.evaluate_scalar(&witness).unwrap().to_bits()
        );
        assert_close(
            requested.evaluate_scalar(&witness).unwrap(),
            expected_scalar,
            2.0e-12,
        );
        assert_eq!(
            requested
                .evaluate_gradient(&witness)
                .unwrap()
                .map(f64::to_bits),
            baseline
                .evaluate_gradient(&witness)
                .unwrap()
                .map(f64::to_bits)
        );
        for (actual, expected) in requested
            .evaluate_gradient(&witness)
            .unwrap()
            .into_iter()
            .zip(expected_gradient)
        {
            assert_close(actual, expected, 4.0e-11);
        }
    }
}
