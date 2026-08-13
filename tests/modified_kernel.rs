use georbf::{
    AnisotropicKernel, Axis, DerivativePoint, Error, FirstDerivative, Interface, IsotropicKernel,
    KernelError, ModifiedKernel, Planar, Point, RbfKernel, SecondDerivative, Tangent,
};

const AXES: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];
const DIFFERENTIABLE_ISOTROPIC: [RbfKernel; 8] = [
    RbfKernel::Cubic,
    RbfKernel::Gaussian,
    RbfKernel::Multiquadric,
    RbfKernel::MultiquadricCubic,
    RbfKernel::ThinPlateSpline,
    RbfKernel::InverseMultiquadric,
    RbfKernel::WendlandC2,
    RbfKernel::MaternC4,
];
const DIFFERENTIABLE_ANISOTROPIC: [RbfKernel; 5] = [
    RbfKernel::Cubic,
    RbfKernel::Gaussian,
    RbfKernel::Multiquadric,
    RbfKernel::ThinPlateSpline,
    RbfKernel::InverseMultiquadric,
];
const SECOND_DERIVATIVES: [SecondDerivative; 9] = [
    SecondDerivative::DxDx,
    SecondDerivative::DxDy,
    SecondDerivative::DxDz,
    SecondDerivative::DyDx,
    SecondDerivative::DyDy,
    SecondDerivative::DyDz,
    SecondDerivative::DzDx,
    SecondDerivative::DzDy,
    SecondDerivative::DzDz,
];

fn point(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z).unwrap()
}

fn interface(x: f64, y: f64, z: f64) -> Interface {
    Interface::new(x, y, z, 20.0).unwrap()
}

fn groups() -> Vec<Vec<Interface>> {
    vec![vec![
        interface(0.0, 0.0, 0.0),
        interface(10.0, 1.0, 1.0),
        interface(2.0, -3.0, 0.0),
        interface(3.0, 4.0, 2.0),
        interface(4.0, 2.0, 0.5),
        interface(5.0, -1.0, 1.5),
    ]]
}

fn oblique_planars() -> Vec<Planar> {
    vec![
        Planar::from_normal(0.0, 0.0, 0.0, 1.0, 0.0, 0.0).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, 0.0, 1.0, 0.0).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, 0.36, -0.48, 0.8).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, -0.48, 0.64, 0.6).unwrap(),
    ]
}

const fn axis_index(axis: Axis) -> usize {
    match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
    }
}

const fn first_derivative(axis: Axis) -> FirstDerivative {
    match axis {
        Axis::X => FirstDerivative::Dx,
        Axis::Y => FirstDerivative::Dy,
        Axis::Z => FirstDerivative::Dz,
    }
}

fn displaced(value: &Point, axis: Axis, amount: f64) -> Point {
    let mut position = value.position();
    position[axis_index(axis)] += amount;
    Point::with_c(position[0], position[1], position[2], value.c()).unwrap()
}

fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let limit = absolute + relative * actual.abs().max(expected.abs());
    assert!(
        delta <= limit,
        "actual={actual:.17e}, expected={expected:.17e}, delta={delta:.3e}, limit={limit:.3e}"
    );
}

fn evaluation_values(kernel: &ModifiedKernel, first: &Point, second: &Point) -> [f64; 16] {
    let evaluation = kernel.evaluate(first, second).unwrap();
    let first_gradient = evaluation.first_at_first();
    let second_gradient = evaluation.first_at_second();
    let hessian = evaluation.mixed_hessian();
    [
        evaluation.basis(),
        first_gradient[0],
        first_gradient[1],
        first_gradient[2],
        second_gradient[0],
        second_gradient[1],
        second_gradient[2],
        hessian[0][0],
        hessian[0][1],
        hessian[0][2],
        hessian[1][0],
        hessian[1][1],
        hessian[1][2],
        hessian[2][0],
        hessian[2][1],
        hessian[2][2],
    ]
}

fn tangent_values(
    kernel: &ModifiedKernel,
    first: &Point,
    second: &Point,
    first_tangent: &Tangent,
    second_tangent: &Tangent,
) -> [f64; 9] {
    [
        kernel.basis_pt_tangent(first, second_tangent).unwrap(),
        kernel.basis_tangent_pt(first_tangent, second).unwrap(),
        kernel
            .basis_tangent_tangent(first_tangent, second_tangent)
            .unwrap(),
        kernel
            .basis_planar_tangent(first, second_tangent, FirstDerivative::Dx)
            .unwrap(),
        kernel
            .basis_planar_tangent(first, second_tangent, FirstDerivative::Dy)
            .unwrap(),
        kernel
            .basis_planar_tangent(first, second_tangent, FirstDerivative::Dz)
            .unwrap(),
        kernel
            .basis_tangent_planar(first_tangent, second, FirstDerivative::Dx)
            .unwrap(),
        kernel
            .basis_tangent_planar(first_tangent, second, FirstDerivative::Dy)
            .unwrap(),
        kernel
            .basis_tangent_planar(first_tangent, second, FirstDerivative::Dz)
            .unwrap(),
    ]
}

#[test]
fn every_combination_matches_frozen_cpp_hexfloat_golden() {
    const ISOTROPIC_EVALUATION: [u64; 16] = [
        0x4081_3660_9881_6242,
        0xc04a_4f64_fe7f_fdb2,
        0xc051_fcba_cd92_7580,
        0x406f_7178_6eaa_8fbb,
        0xc069_ff9a_8c7a_23dc,
        0xc062_edf6_108e_a854,
        0x4085_1250_e9a5_ce9c,
        0x4045_01b9_6ebe_829e,
        0x4037_0d80_bc4b_03c7,
        0xc054_bdd5_6fe5_5c80,
        0x4037_a463_df07_288a,
        0x4039_0230_709c_7320,
        0xc054_d1fb_d4e4_2390,
        0xc054_d25c_c6f4_b28a,
        0xc050_d4b4_3b8f_61f8,
        0x4072_b3ae_ad1a_ad7b,
    ];
    const ISOTROPIC_TANGENT: [u64; 9] = [
        0x407a_6bd0_8363_258a,
        0x4061_5886_cd5c_2c57,
        0x4058_11a9_ba0f_40d2,
        0xc04a_52dc_7a50_b416,
        0xc048_028f_382a_c680,
        0x4067_372e_2792_cdbc,
        0xc043_3f85_ebf0_2d10,
        0xc042_5fea_9c3d_053a,
        0x4063_c134_3a19_19eb,
    ];
    const ISOTROPIC_ZERO: [u64; 16] = [
        0x4098_2ab3_edd7_c9f9,
        0xc06c_9eab_f3e9_c958,
        0xc065_b94c_4268_d058,
        0x4084_e5d1_6cd3_fe39,
        0xc06c_9eab_f3e9_c958,
        0xc065_b94c_4268_d056,
        0x4084_e5d1_6cd3_fe39,
        0x4050_63f0_9918_e5a1,
        0x4034_713b_f61f_2cb2,
        0xc055_f073_f9ed_aa28,
        0x4034_713b_f61f_2cae,
        0x4044_1708_1b0f_8575,
        0xc054_0ca9_7f46_a7ee,
        0xc055_f073_f9ed_aa28,
        0xc054_0ca9_7f46_a7f0,
        0x4073_0a48_b093_2f6e,
    ];
    const ANISOTROPIC_EVALUATION: [u64; 16] = [
        0x408c_4479_384b_9013,
        0xc050_aec6_a637_2d8d,
        0xc062_73be_4d74_f594,
        0x407a_673b_295e_0326,
        0xc073_0fa1_5cb5_b8bd,
        0xc06c_a25a_0e30_4e7b,
        0x4090_f384_06a2_0006,
        0x404d_7b67_e0e1_99c4,
        0x4038_bf84_e583_9e56,
        0xc05c_7331_4a1b_8690,
        0x403b_88ae_9833_5832,
        0x4045_c150_08b9_881a,
        0xc063_4ac5_197e_151a,
        0xc05e_4c8e_79ef_04d1,
        0xc05b_4a33_ec00_3974,
        0x407e_9c0f_e3d2_0677,
    ];
    const ANISOTROPIC_TANGENT: [u64; 9] = [
        0x4085_57bb_76d4_2056,
        0x406f_481c_af6a_eb59,
        0x4064_58ed_6828_2594,
        0xc053_26a5_c189_4905,
        0xc055_dbeb_0f15_05ec,
        0x4072_d924_0ef8_88f3,
        0xc04a_f612_3b4f_3520,
        0xc050_23cc_ccee_8c12,
        0x4071_078b_ab0c_f6f6,
    ];
    const ANISOTROPIC_ZERO: [u64; 16] = [
        0x40a3_0ef3_ff6d_15b2,
        0xc074_c54e_1b4d_bef6,
        0xc071_be8e_e742_9de6,
        0x4090_a50e_fc5a_a20b,
        0xc074_c54e_1b4d_bef5,
        0xc071_be8e_e742_9de7,
        0x4090_a50e_fc5a_a20b,
        0x4059_1a6b_4a39_f886,
        0x401b_dba0_6ec0_f34e,
        0xc05f_177b_29b8_988a,
        0x401b_dba0_6ec0_f35a,
        0x4054_ddbb_6a81_72ff,
        0xc061_d936_e13d_1cae,
        0xc05f_177b_29b8_988a,
        0xc061_d936_e13d_1cab,
        0x407f_7340_9cea_1160,
    ];

    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let first_tangent = Tangent::new(-1.25, 0.5, 2.0, 0.3, -0.4, 0.5).unwrap();
    let second_tangent = Tangent::new(0.75, -0.25, 1.0, -0.2, 0.6, 0.7).unwrap();
    let isotropic =
        ModifiedKernel::from_isotropic(IsotropicKernel::new(RbfKernel::Cubic, 0.7), &groups())
            .unwrap();
    let anisotropic = ModifiedKernel::from_anisotropic(
        AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &oblique_planars()).unwrap(),
        &groups(),
    )
    .unwrap();

    assert_eq!(
        evaluation_values(&isotropic, &first, &second).map(f64::to_bits),
        ISOTROPIC_EVALUATION
    );
    assert_eq!(
        tangent_values(&isotropic, &first, &second, &first_tangent, &second_tangent,)
            .map(f64::to_bits),
        ISOTROPIC_TANGENT
    );
    assert_eq!(
        evaluation_values(&isotropic, &first, &first).map(f64::to_bits),
        ISOTROPIC_ZERO
    );
    assert_eq!(
        evaluation_values(&anisotropic, &first, &second).map(f64::to_bits),
        ANISOTROPIC_EVALUATION
    );
    assert_eq!(
        tangent_values(
            &anisotropic,
            &first,
            &second,
            &first_tangent,
            &second_tangent,
        )
        .map(f64::to_bits),
        ANISOTROPIC_TANGENT
    );
    assert_eq!(
        evaluation_values(&anisotropic, &first, &first).map(f64::to_bits),
        ANISOTROPIC_ZERO
    );
}

#[test]
fn construction_maps_lagrangian_failure_to_modified_kernel_category() {
    let base = IsotropicKernel::new(RbfKernel::Cubic, 0.7);
    assert_eq!(
        ModifiedKernel::from_isotropic(base, &[]).unwrap_err(),
        Error::ModifiedKernelCreationFailure
    );
    assert_eq!(
        ModifiedKernel::from_isotropic(
            base,
            &[vec![
                interface(0.0, 0.0, 0.0),
                interface(1.0, 0.0, 0.0),
                interface(0.0, 1.0, 0.0),
            ]],
        )
        .unwrap_err(),
        Error::ModifiedKernelCreationFailure
    );
}

#[test]
fn every_value_planar_tangent_combination_matches_evaluation_components() {
    let kernel =
        ModifiedKernel::from_isotropic(IsotropicKernel::new(RbfKernel::Cubic, 0.7), &groups())
            .unwrap();
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let first_tangent = Tangent::new(first.x(), first.y(), first.z(), 0.3, -0.4, 0.5).unwrap();
    let second_tangent = Tangent::new(second.x(), second.y(), second.z(), -0.2, 0.6, 0.7).unwrap();
    let evaluation = kernel.evaluate(&first, &second).unwrap();

    assert_eq!(kernel.basis_pt_pt(&first, &second), evaluation.basis());
    for axis in AXES {
        let index = axis_index(axis);
        assert_eq!(
            kernel.basis_planar_pt(&first, &second, axis).unwrap(),
            evaluation.first_at_first()[index]
        );
        assert_eq!(
            kernel.basis_pt_planar(&first, &second, axis).unwrap(),
            evaluation.first_at_second()[index]
        );
    }
    for (index, component) in SECOND_DERIVATIVES.into_iter().enumerate() {
        assert_eq!(
            kernel
                .basis_planar_planar(&first, &second, component)
                .unwrap(),
            evaluation.mixed_hessian()[index / 3][index % 3]
        );
    }

    let first_vector = first_tangent.vector();
    let second_vector = second_tangent.vector();
    let expected_pt_tangent = (0..3)
        .map(|index| evaluation.first_at_second()[index] * second_vector[index])
        .sum::<f64>();
    let expected_tangent_pt = (0..3)
        .map(|index| evaluation.first_at_first()[index] * first_vector[index])
        .sum::<f64>();
    assert_eq!(
        kernel.basis_pt_tangent(&first, &second_tangent).unwrap(),
        expected_pt_tangent
    );
    assert_eq!(
        kernel.basis_tangent_pt(&first_tangent, &second).unwrap(),
        expected_tangent_pt
    );

    let hessian = evaluation.mixed_hessian();
    let expected_tangent_tangent = (0..3)
        .flat_map(|row| (0..3).map(move |column| (row, column)))
        .map(|(row, column)| first_vector[row] * second_vector[column] * hessian[row][column])
        .sum::<f64>();
    assert_eq!(
        kernel
            .basis_tangent_tangent(&first_tangent, &second_tangent)
            .unwrap(),
        expected_tangent_tangent
    );
    for axis in AXES {
        let row = axis_index(axis);
        let expected_planar_tangent = (0..3)
            .map(|column| second_vector[column] * hessian[row][column])
            .sum::<f64>();
        let expected_tangent_planar = (0..3)
            .map(|first_axis| first_vector[first_axis] * hessian[first_axis][row])
            .sum::<f64>();
        assert_eq!(
            kernel
                .basis_planar_tangent(&first, &second_tangent, first_derivative(axis))
                .unwrap(),
            expected_planar_tangent
        );
        assert_eq!(
            kernel
                .basis_tangent_planar(&first_tangent, &second, first_derivative(axis))
                .unwrap(),
            expected_tangent_planar
        );
    }
}

#[test]
fn gaussian_unit_diagonal_is_eliminated_at_every_unisolvent_point() {
    let kernel =
        ModifiedKernel::from_isotropic(IsotropicKernel::new(RbfKernel::Gaussian, 1.0), &groups())
            .unwrap();
    let sample = point(1.25, -0.75, 2.5);
    let mut actual_bits = Vec::with_capacity(8);
    for unisolvent in kernel.lagrangian_basis().unisolvent_points() {
        let forward = kernel.basis_pt_pt(unisolvent, &sample);
        let reverse = kernel.basis_pt_pt(&sample, unisolvent);
        actual_bits.extend([forward.to_bits(), reverse.to_bits()]);
        assert_close(forward, 0.0, 2.0e-14, 2.0e-14);
        assert_close(reverse, 0.0, 2.0e-14, 2.0e-14);
    }
    assert_eq!(
        actual_bits,
        [
            0xbc93_d7a0_0000_0000,
            0xbc93_d7a0_0000_0000,
            0x3cc0_0000_0006_158f,
            0x3cc0_0000_0006_158f,
            0x3cc7_2b42_0000_0000,
            0x3cc7_2b42_0000_0000,
            0x3c93_c7ec_1bb8_2000,
            0x3c93_c7ec_1bb8_2000,
        ]
    );
}

#[test]
fn exchange_symmetry_and_finite_difference_hold_for_isotropic_and_anisotropic() {
    let isotropic =
        ModifiedKernel::from_isotropic(IsotropicKernel::new(RbfKernel::Gaussian, 0.7), &groups())
            .unwrap();
    let anisotropic = ModifiedKernel::from_anisotropic(
        AnisotropicKernel::new(RbfKernel::Gaussian, 0.7, &oblique_planars()).unwrap(),
        &groups(),
    )
    .unwrap();
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let step = 2.0e-5;

    for kernel in [isotropic, anisotropic] {
        let forward = kernel.evaluate(&first, &second).unwrap();
        let reverse = kernel.evaluate(&second, &first).unwrap();
        assert_close(forward.basis(), reverse.basis(), 2.0e-13, 2.0e-13);
        for row in 0..3 {
            assert_close(
                forward.first_at_first()[row],
                reverse.first_at_second()[row],
                2.0e-13,
                2.0e-13,
            );
            for column in 0..3 {
                assert_close(
                    forward.mixed_hessian()[row][column],
                    reverse.mixed_hessian()[column][row],
                    5.0e-12,
                    5.0e-12,
                );
            }
        }

        for axis in AXES {
            let plus = displaced(&first, axis, step);
            let minus = displaced(&first, axis, -step);
            let numerical = (kernel.basis_pt_pt(&plus, &second)
                - kernel.basis_pt_pt(&minus, &second))
                / (2.0 * step);
            assert_close(
                forward.first_at_first()[axis_index(axis)],
                numerical,
                2.0e-8,
                2.0e-8,
            );
            for second_axis in AXES {
                let plus_second = displaced(&second, second_axis, step);
                let minus_second = displaced(&second, second_axis, -step);
                let plus_gradient = kernel.basis_planar_pt(&first, &plus_second, axis).unwrap();
                let minus_gradient = kernel.basis_planar_pt(&first, &minus_second, axis).unwrap();
                let numerical_hessian = (plus_gradient - minus_gradient) / (2.0 * step);
                assert_close(
                    forward.mixed_hessian()[axis_index(axis)][axis_index(second_axis)],
                    numerical_hessian,
                    2.0e-7,
                    2.0e-7,
                );
            }
        }
    }
}

#[test]
fn every_supported_kernel_and_zero_distance_path_is_finite() {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    for kind in DIFFERENTIABLE_ISOTROPIC {
        let parameter = if kind == RbfKernel::WendlandC2 {
            10.0
        } else {
            0.7
        };
        let kernel =
            ModifiedKernel::from_isotropic(IsotropicKernel::new(kind, parameter), &groups())
                .unwrap();
        assert!(kernel.evaluate(&first, &second).unwrap().is_finite());
        assert!(kernel.evaluate(&first, &first).unwrap().is_finite());
    }
    for kind in DIFFERENTIABLE_ANISOTROPIC {
        let kernel = ModifiedKernel::from_anisotropic(
            AnisotropicKernel::new(kind, 0.7, &oblique_planars()).unwrap(),
            &groups(),
        )
        .unwrap();
        assert!(kernel.evaluate(&first, &second).unwrap().is_finite());
        assert!(kernel.evaluate(&first, &first).unwrap().is_finite());
    }
}

#[test]
fn linear_value_remains_available_and_every_derivative_returns_typed_error() {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let isotropic =
        ModifiedKernel::from_isotropic(IsotropicKernel::new(RbfKernel::Linear, 0.7), &groups())
            .unwrap();
    let anisotropic = ModifiedKernel::from_anisotropic(
        AnisotropicKernel::new(RbfKernel::Linear, 0.7, &oblique_planars()).unwrap(),
        &groups(),
    )
    .unwrap();
    assert_eq!(
        isotropic.basis_pt_pt(&first, &second).to_bits(),
        0xc020_0a91_80e9_3007
    );
    for kernel in [isotropic, anisotropic] {
        assert!(kernel.basis_pt_pt(&first, &second).is_finite());
        assert_eq!(
            kernel
                .basis_planar_pt(&first, &second, Axis::X)
                .unwrap_err(),
            KernelError::LinearDerivativeUnavailable
        );
        assert_eq!(
            kernel
                .basis_pt_planar(&first, &second, Axis::X)
                .unwrap_err(),
            KernelError::LinearDerivativeUnavailable
        );
        assert_eq!(
            kernel
                .basis_planar_planar(&first, &second, SecondDerivative::DxDx)
                .unwrap_err(),
            KernelError::LinearDerivativeUnavailable
        );
        assert_eq!(
            kernel.evaluate(&first, &second).unwrap_err(),
            KernelError::LinearDerivativeUnavailable
        );
    }
}

#[test]
fn modified_first_and_hessian_are_the_derivatives_of_the_modified_value() {
    let kernel = ModifiedKernel::from_isotropic(
        IsotropicKernel::new(RbfKernel::Multiquadric, 0.7),
        &groups(),
    )
    .unwrap();
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let evaluation = kernel.evaluate(&first, &second).unwrap();
    let step = 2.0e-5;
    for first_axis in AXES {
        let plus = displaced(&first, first_axis, step);
        let minus = displaced(&first, first_axis, -step);
        let numerical = (kernel.basis_pt_pt(&plus, &second) - kernel.basis_pt_pt(&minus, &second))
            / (2.0 * step);
        assert_close(
            evaluation.first_at_first()[axis_index(first_axis)],
            numerical,
            2.0e-9,
            2.0e-9,
        );
        for second_axis in AXES {
            let plus_second = displaced(&second, second_axis, step);
            let minus_second = displaced(&second, second_axis, -step);
            let plus_gradient = kernel
                .basis_planar_pt(&first, &plus_second, first_axis)
                .unwrap();
            let minus_gradient = kernel
                .basis_planar_pt(&first, &minus_second, first_axis)
                .unwrap();
            let numerical_hessian = (plus_gradient - minus_gradient) / (2.0 * step);
            assert_close(
                evaluation.mixed_hessian()[axis_index(first_axis)][axis_index(second_axis)],
                numerical_hessian,
                2.0e-8,
                2.0e-8,
            );
        }
    }

    for axis in AXES {
        assert_close(
            kernel.basis_planar_pt(&first, &second, axis).unwrap(),
            kernel
                .first_derivative(&first, &second, DerivativePoint::First, axis)
                .unwrap(),
            0.0,
            0.0,
        );
    }
}
