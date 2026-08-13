use georbf::{
    Axis, DerivativePoint, IsotropicKernel, KernelError, Point, RbfKernel, SecondDerivative,
};

const KERNELS: [RbfKernel; 9] = [
    RbfKernel::Cubic,
    RbfKernel::Gaussian,
    RbfKernel::Multiquadric,
    RbfKernel::MultiquadricCubic,
    RbfKernel::ThinPlateSpline,
    RbfKernel::InverseMultiquadric,
    RbfKernel::Linear,
    RbfKernel::WendlandC2,
    RbfKernel::MaternC4,
];

const DIFFERENTIABLE_KERNELS: [RbfKernel; 8] = [
    RbfKernel::Cubic,
    RbfKernel::Gaussian,
    RbfKernel::Multiquadric,
    RbfKernel::MultiquadricCubic,
    RbfKernel::ThinPlateSpline,
    RbfKernel::InverseMultiquadric,
    RbfKernel::WendlandC2,
    RbfKernel::MaternC4,
];

const AXES: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];
const SECOND_DERIVATIVES: [(SecondDerivative, Axis, Axis); 9] = [
    (SecondDerivative::DxDx, Axis::X, Axis::X),
    (SecondDerivative::DxDy, Axis::X, Axis::Y),
    (SecondDerivative::DxDz, Axis::X, Axis::Z),
    (SecondDerivative::DyDx, Axis::Y, Axis::X),
    (SecondDerivative::DyDy, Axis::Y, Axis::Y),
    (SecondDerivative::DyDz, Axis::Y, Axis::Z),
    (SecondDerivative::DzDx, Axis::Z, Axis::X),
    (SecondDerivative::DzDy, Axis::Z, Axis::Y),
    (SecondDerivative::DzDz, Axis::Z, Axis::Z),
];

fn point(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z).unwrap()
}

fn point_with_c(x: f64, y: f64, z: f64, c: f64) -> Point {
    Point::with_c(x, y, z, c).unwrap()
}

fn displaced(point: &Point, axis: Axis, amount: f64) -> Point {
    let mut position = point.position();
    position[axis_index(axis)] += amount;
    point_with_c(position[0], position[1], position[2], point.c())
}

const fn axis_index(axis: Axis) -> usize {
    match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
    }
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.1e}"
    );
}

#[test]
fn separated_case_matches_frozen_cpp_hexfloat_golden() {
    const GOLDEN: [(RbfKernel, [u64; 16]); 8] = [
        (
            RbfKernel::Cubic,
            [
                0x402a_3cfe_9be1_0415,
                0xc02c_4d4c_2993_cae2,
                0x4015_39f9_1f2e_d82a,
                0x401c_4d4c_2993_cae2,
                0x402c_4d4c_2993_cae2,
                0xc015_39f9_1f2e_d82a,
                0xc01c_4d4c_2993_cae2,
                0xc028_53b3_155b_1480,
                0x3ffe_8727_01b3_8d2e,
                0x4004_5a1a_0122_5e1e,
                0x3ffe_8727_01b3_8d2e,
                0xc01f_29f7_d1bc_a01e,
                0xbfee_8727_01b3_8d2e,
                0x4004_5a1a_0122_5e1e,
                0xbfee_8727_01b3_8d2e,
                0xc020_b1e9_54ee_3135,
            ],
        ),
        (
            RbfKernel::Gaussian,
            [
                0x3fb0_c4f3_b6c6_82a7,
                0x3fc0_6f17_ccb8_4cd6,
                0xbfa8_a6a3_b314_7341,
                0xbfb0_6f17_ccb8_4cd6,
                0xbfc0_6f17_ccb8_4cd6,
                0x3fa8_a6a3_b314_7341,
                0x3fb0_6f17_ccb8_4cd6,
                0xbfc7_fe5b_114a_84a9,
                0x3fb8_286d_39bd_004f,
                0x3fc0_1af3_7bd3_558a,
                0x3fb8_286d_39bd_004f,
                0x3f9d_7fbb_5c45_b2e3,
                0xbfa8_286d_39bd_004f,
                0x3fc0_1af3_7bd3_558a,
                0xbfa8_286d_39bd_004f,
                0x3f55_0914_393d_d311,
            ],
        ),
        (
            RbfKernel::Multiquadric,
            [
                0x4004_051e_10b7_24a2,
                0xbfe9_930e_640b_2b39,
                0x3fd3_2e4a_cb08_606b,
                0x3fd9_930e_640b_2b39,
                0x3fe9_930e_640b_2b39,
                0xbfd3_2e4a_cb08_606b,
                0xbfd9_930e_640b_2b39,
                0xbfc2_7a9b_8412_4a26,
                0xbfb8_80a0_f303_0939,
                0xbfc0_55c0_a202_0626,
                0xbfb8_80a0_f303_0939,
                0xbfd7_46ff_4d42_e25c,
                0x3fa8_80a0_f303_0939,
                0xbfc0_55c0_a202_0626,
                0x3fa8_80a0_f303_0939,
                0xbfd5_7d9e_3b8a_a9b0,
            ],
        ),
        (
            RbfKernel::MultiquadricCubic,
            [
                0x402f_5803_122b_888d,
                0xc02e_07ad_1912_b6f3,
                0x4016_85c1_d2ce_0936,
                0x401e_07ad_1912_b6f3,
                0x402e_07ad_1912_b6f3,
                0xc016_85c1_d2ce_0936,
                0xc01e_07ad_1912_b6f3,
                0xc028_9afb_f20d_8baf,
                0x3ffc_c570_308c_90a0,
                0x4003_2e4a_cb08_606b,
                0x3ffc_c570_308c_90a0,
                0xc020_5d17_cecf_f241,
                0xbfec_c570_308c_90a0,
                0x4003_2e4a_cb08_606b,
                0xbfec_c570_308c_90a0,
                0xc021_699f_e5ea_6787,
            ],
        ),
        (
            RbfKernel::ThinPlateSpline,
            [
                0x403a_8c68_1e02_2d8a,
                0xc048_a74d_bf49_79eb,
                0x4032_7d7a_4f77_1b71,
                0x4038_a74d_bf49_79eb,
                0x4048_a74d_bf49_79eb,
                0xc032_7d7a_4f77_1b71,
                0xc038_a74d_bf49_79eb,
                0xc053_070f_0803_cc79,
                0x4033_4bd9_644a_24fe,
                0x4039_ba77_3062_dbfd,
                0x4033_4bd9_644a_24fe,
                0xc03f_e3bf_44e5_47ca,
                0xc023_4bd9_644a_24fe,
                0x4039_ba77_3062_dbfd,
                0xc023_4bd9_644a_24fe,
                0xc042_c244_abbd_73f4,
            ],
        ),
        (
            RbfKernel::InverseMultiquadric,
            [
                0x3fd9_930e_640b_2b39,
                0x3fc0_55c0_a202_0626,
                0xbfa8_80a0_f303_0939,
                0xbfb0_55c0_a202_0626,
                0xbfc0_55c0_a202_0626,
                0x3fa8_80a0_f303_0939,
                0x3fb0_55c0_a202_0626,
                0xbfad_ee5f_6f5e_26de,
                0x3fa7_79b4_4344_d32f,
                0x3faf_4cf0_59b1_1995,
                0x3fa7_79b4_4344_d32f,
                0x3fa7_dddd_aaca_3d1a,
                0xbf97_79b4_4344_d32f,
                0x3faf_4cf0_59b1_1995,
                0xbf97_79b4_4344_d32f,
                0x3fa1_0509_172b_7f82,
            ],
        ),
        (RbfKernel::WendlandC2, [0; 16]),
        (
            RbfKernel::MaternC4,
            [
                0x4000_640f_7e76_eef7,
                0x3fdf_e6c3_6917_c07c,
                0xbfc7_ed12_8ed1_d05d,
                0xbfcf_e6c3_6917_c07c,
                0xbfdf_e6c3_6917_c07c,
                0x3fc7_ed12_8ed1_d05d,
                0x3fcf_e6c3_6917_c07c,
                0x3fb0_a131_d999_d74c,
                0x3fb1_b09f_dd38_1fa0,
                0x3fb7_962a_7c4a_d4d6,
                0x3fb1_b09f_dd38_1fa0,
                0x3fcc_95a5_6f9d_3a8f,
                0xbfa1_b09f_dd38_1fa0,
                0x3fb7_962a_7c4a_d4d6,
                0xbfa1_b09f_dd38_1fa0,
                0x3fca_0138_ca05_0b47,
            ],
        ),
    ];

    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    for (kind, expected) in GOLDEN {
        let evaluation = IsotropicKernel::new(kind, 0.7)
            .evaluate(&first, &second)
            .unwrap();
        let mut actual = Vec::with_capacity(16);
        actual.push(evaluation.basis());
        actual.extend(evaluation.first_at_first());
        actual.extend(evaluation.first_at_second());
        actual.extend(evaluation.mixed_hessian().into_iter().flatten());
        assert_eq!(
            actual.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
            expected,
            "frozen hexfloat mismatch for {kind:?}"
        );
    }

    assert_eq!(
        IsotropicKernel::new(RbfKernel::Linear, 0.7)
            .basis(&first, &second)
            .to_bits(),
        0x4002_de32_c662_8741
    );
}

#[test]
fn nine_kernel_values_and_eight_complete_derivative_surfaces_are_available() {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);

    for kind in KERNELS {
        let kernel = IsotropicKernel::new(kind, 0.7);
        assert!(kernel.basis(&first, &second).is_finite());
        if kind == RbfKernel::Linear {
            continue;
        }
        for axis in AXES {
            assert!(kernel
                .first_derivative(&first, &second, DerivativePoint::First, axis)
                .unwrap()
                .is_finite());
            assert!(kernel
                .first_derivative(&first, &second, DerivativePoint::Second, axis)
                .unwrap()
                .is_finite());
        }
        for first_axis in AXES {
            for second_axis in AXES {
                assert!(kernel
                    .mixed_second_derivative(&first, &second, first_axis, second_axis)
                    .unwrap()
                    .is_finite());
            }
        }
        for (component, first_axis, second_axis) in SECOND_DERIVATIVES {
            assert_eq!(
                kernel
                    .mixed_second_derivative_component(&first, &second, component)
                    .unwrap(),
                kernel
                    .mixed_second_derivative(&first, &second, first_axis, second_axis)
                    .unwrap()
            );
        }
        assert!(kernel.evaluate(&first, &second).unwrap().is_finite());
    }
}

#[test]
fn radial_exchange_sign_and_hessian_symmetry_hold() {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);

    for kind in DIFFERENTIABLE_KERNELS {
        let kernel = IsotropicKernel::new(kind, 0.7);
        assert_close(
            kernel.basis(&first, &second),
            kernel.basis(&second, &first),
            2.0e-14,
        );
        for axis in AXES {
            let first_derivative = kernel
                .first_derivative(&first, &second, DerivativePoint::First, axis)
                .unwrap();
            let second_derivative = kernel
                .first_derivative(&first, &second, DerivativePoint::Second, axis)
                .unwrap();
            let swapped_first = kernel
                .first_derivative(&second, &first, DerivativePoint::First, axis)
                .unwrap();
            assert_close(first_derivative, -second_derivative, 2.0e-14);
            assert_close(second_derivative, swapped_first, 2.0e-14);
        }
        for first_axis in AXES {
            for second_axis in AXES {
                let mixed = kernel
                    .mixed_second_derivative(&first, &second, first_axis, second_axis)
                    .unwrap();
                let transposed = kernel
                    .mixed_second_derivative(&first, &second, second_axis, first_axis)
                    .unwrap();
                let swapped = kernel
                    .mixed_second_derivative(&second, &first, first_axis, second_axis)
                    .unwrap();
                assert_close(mixed, transposed, 2.0e-14);
                assert_close(mixed, swapped, 2.0e-14);
            }
        }
    }
}

#[test]
fn analytic_derivatives_match_finite_differences_at_both_points() {
    let step = 1.0e-5;

    for kind in DIFFERENTIABLE_KERNELS {
        let (first, second) = if kind == RbfKernel::WendlandC2 {
            (point(0.2, -0.15, 0.1), point(0.0, 0.0, 0.0))
        } else {
            (point(-1.25, 0.5, 2.0), point(0.75, -0.25, 1.0))
        };
        let kernel = IsotropicKernel::new(kind, 0.7);
        for axis in AXES {
            let first_lower = displaced(&first, axis, -step);
            let first_upper = displaced(&first, axis, step);
            let numeric_first = (kernel.basis(&first_upper, &second)
                - kernel.basis(&first_lower, &second))
                / (2.0 * step);
            let analytic_first = kernel
                .first_derivative(&first, &second, DerivativePoint::First, axis)
                .unwrap();
            assert_close(analytic_first, numeric_first, 3.0e-8);

            let second_lower = displaced(&second, axis, -step);
            let second_upper = displaced(&second, axis, step);
            let numeric_second = (kernel.basis(&first, &second_upper)
                - kernel.basis(&first, &second_lower))
                / (2.0 * step);
            let analytic_second = kernel
                .first_derivative(&first, &second, DerivativePoint::Second, axis)
                .unwrap();
            assert_close(analytic_second, numeric_second, 3.0e-8);
        }

        for first_axis in AXES {
            for second_axis in AXES {
                let second_lower = displaced(&second, second_axis, -step);
                let second_upper = displaced(&second, second_axis, step);
                let lower = kernel
                    .first_derivative(&first, &second_lower, DerivativePoint::First, first_axis)
                    .unwrap();
                let upper = kernel
                    .first_derivative(&first, &second_upper, DerivativePoint::First, first_axis)
                    .unwrap();
                let numeric = (upper - lower) / (2.0 * step);
                let analytic = kernel
                    .mixed_second_derivative(&first, &second, first_axis, second_axis)
                    .unwrap();
                assert_close(analytic, numeric, 5.0e-8);
            }
        }
    }
}

#[test]
fn zero_distance_branches_match_frozen_surfe() {
    let same = point(2.0, -3.0, 4.0);
    let shape = 0.7_f64;
    let expected = [
        (RbfKernel::Cubic, 0.0, 0.0),
        (RbfKernel::Gaussian, 1.0, 2.0 * shape * shape),
        (
            RbfKernel::Multiquadric,
            shape.powf(0.5),
            -1.0 / shape.powf(0.5),
        ),
        (
            RbfKernel::MultiquadricCubic,
            shape.powf(1.5),
            -3.0 * shape.powf(0.5),
        ),
        (RbfKernel::ThinPlateSpline, 0.0, 0.0),
        (
            RbfKernel::InverseMultiquadric,
            1.0 / shape.powf(0.5),
            1.0 / shape.powf(1.5),
        ),
        (RbfKernel::WendlandC2, 1.0, 20.0 / (shape * shape)),
        (RbfKernel::MaternC4, 3.0, shape * shape),
    ];

    for (kind, expected_basis, expected_diagonal) in expected {
        let kernel = IsotropicKernel::new(kind, shape);
        assert_close(kernel.basis(&same, &same), expected_basis, 2.0e-14);
        for axis in AXES {
            assert_eq!(
                kernel
                    .first_derivative(&same, &same, DerivativePoint::First, axis)
                    .unwrap(),
                0.0
            );
            assert_close(
                kernel
                    .mixed_second_derivative(&same, &same, axis, axis)
                    .unwrap(),
                expected_diagonal,
                2.0e-14,
            );
        }
        for (first_axis, second_axis) in
            [(Axis::X, Axis::Y), (Axis::X, Axis::Z), (Axis::Y, Axis::Z)]
        {
            assert_eq!(
                kernel
                    .mixed_second_derivative(&same, &same, first_axis, second_axis)
                    .unwrap(),
                0.0
            );
        }
    }
}

#[test]
fn wendland_support_boundary_uses_strict_frozen_branch() {
    let center = point(0.0, 0.0, 0.0);
    let support = 2.0;
    let kernel = IsotropicKernel::new(RbfKernel::WendlandC2, support);
    let just_inside = point(support * (1.0 - 2.0_f64.powi(-20)), 0.0, 0.0);
    let boundary = point(support, 0.0, 0.0);
    let just_outside = point(support * (1.0 + 2.0_f64.powi(-20)), 0.0, 0.0);

    assert!(kernel.basis(&just_inside, &center) > 0.0);
    assert_eq!(kernel.basis(&boundary, &center), 0.0);
    assert_eq!(kernel.basis(&just_outside, &center), 0.0);
    for point in [&boundary, &just_outside] {
        assert_eq!(
            kernel
                .first_derivative(point, &center, DerivativePoint::First, Axis::X)
                .unwrap(),
            0.0
        );
        assert_eq!(
            kernel
                .mixed_second_derivative(point, &center, Axis::X, Axis::X)
                .unwrap(),
            0.0
        );
    }

    let interior = point(0.2, -0.15, 0.1);
    let evaluation = IsotropicKernel::new(RbfKernel::WendlandC2, 0.7)
        .evaluate(&interior, &center)
        .unwrap();
    let mut actual = Vec::with_capacity(16);
    actual.push(evaluation.basis());
    actual.extend(evaluation.first_at_first());
    actual.extend(evaluation.first_at_second());
    actual.extend(evaluation.mixed_hessian().into_iter().flatten());
    assert_eq!(
        actual.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
        [
            0x3fd7_4b63_9057_bae5,
            0xbffe_6ec6_7332_607d,
            0x3ff6_d314_d665_c85e,
            0xbfee_6ec6_7332_607d,
            0x3ffe_6ec6_7332_607d,
            0xbff6_d314_d665_c85e,
            0x3fee_6ec6_7332_607d,
            0xbfd5_17ad_2a6c_01a9,
            0x401d_84f6_29fc_4a88,
            0xc013_adf9_7152_dc5c,
            0x401d_84f6_29fc_4a88,
            0x400f_cd7e_e103_8171,
            0x400d_84f6_29fc_4a88,
            0xc013_adf9_7152_dc5c,
            0x400d_84f6_29fc_4a88,
            0x401c_337b_5755_8a71,
        ]
    );
}

#[test]
fn fourth_coordinate_participates_in_radius_and_preserves_source_hessian_order() {
    let first = point_with_c(-1.25, 0.5, 2.0, 1.5);
    let second = point_with_c(0.75, -0.25, 1.0, -0.25);
    let without_c_first = point(-1.25, 0.5, 2.0);
    let without_c_second = point(0.75, -0.25, 1.0);

    for kind in KERNELS {
        let kernel = IsotropicKernel::new(kind, 10.0);
        assert_ne!(
            kernel.basis(&first, &second),
            kernel.basis(&without_c_first, &without_c_second)
        );
    }
    let tps = IsotropicKernel::new(RbfKernel::ThinPlateSpline, 10.0);
    let evaluation = tps.evaluate(&first, &second).unwrap();
    let mut actual = Vec::with_capacity(16);
    actual.push(evaluation.basis());
    actual.extend(evaluation.first_at_first());
    actual.extend(evaluation.first_at_second());
    actual.extend(evaluation.mixed_hessian().into_iter().flatten());
    assert_eq!(
        actual.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
        [
            0x4054_092e_10d3_2129,
            0xc056_e580_0f9b_04ca,
            0x4041_2c20_0bb4_4397,
            0x4046_e580_0f9b_04ca,
            0x4056_e580_0f9b_04ca,
            0xc041_2c20_0bb4_4397,
            0xc046_e580_0f9b_04ca,
            0xc059_4d20_85c3_c766,
            0x4035_ed90_bcf1_6782,
            0x403d_3cc0_fbec_8a02,
            0x4035_ed90_bcf1_6782,
            0xc049_7a0b_3308_4832,
            0xc025_ed90_bcf1_6782,
            0x403d_3cc0_fbec_8a02,
            0xc025_ed90_bcf1_6782,
            0xc04c_acb0_4e96_274a,
        ]
    );
}

#[test]
fn finite_multiscale_inputs_remain_finite() {
    for scale in [1.0e-3, 1.0, 1.0e3] {
        let first = point(-1.25 * scale, 0.5 * scale, 2.0 * scale);
        let second = point(0.75 * scale, -0.25 * scale, 1.0 * scale);
        for kind in DIFFERENTIABLE_KERNELS {
            let evaluation = IsotropicKernel::new(kind, 0.7)
                .evaluate(&first, &second)
                .unwrap();
            assert!(evaluation.is_finite(), "non-finite {kind:?} at {scale}");
        }
    }
}

#[test]
fn linear_kernel_derivative_sentinel_is_a_typed_error() {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let kernel = IsotropicKernel::new(RbfKernel::Linear, 0.7);

    assert!(kernel.basis(&first, &second).is_finite());
    assert_eq!(
        kernel.first_derivative(&first, &second, DerivativePoint::First, Axis::X),
        Err(KernelError::LinearDerivativeUnavailable)
    );
    assert_eq!(
        kernel.mixed_second_derivative(&first, &second, Axis::X, Axis::X),
        Err(KernelError::LinearDerivativeUnavailable)
    );
    assert_eq!(
        kernel.evaluate(&first, &second),
        Err(KernelError::LinearDerivativeUnavailable)
    );
}
