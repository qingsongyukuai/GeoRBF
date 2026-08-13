use georbf::{
    AnisotropicKernel, AnisotropyError, Axis, DerivativePoint, KernelError, Planar, Point,
    RbfKernel,
};

const KERNELS: [RbfKernel; 6] = [
    RbfKernel::Cubic,
    RbfKernel::Gaussian,
    RbfKernel::Multiquadric,
    RbfKernel::ThinPlateSpline,
    RbfKernel::InverseMultiquadric,
    RbfKernel::Linear,
];

const DIFFERENTIABLE_KERNELS: [RbfKernel; 5] = [
    RbfKernel::Cubic,
    RbfKernel::Gaussian,
    RbfKernel::Multiquadric,
    RbfKernel::ThinPlateSpline,
    RbfKernel::InverseMultiquadric,
];

const AXES: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

fn point(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z).unwrap()
}

fn planar(nx: f64, ny: f64, nz: f64) -> Planar {
    Planar::from_normal(0.0, 0.0, 0.0, nx, ny, nz).unwrap()
}

fn oblique_planars() -> Vec<Planar> {
    vec![
        planar(1.0, 0.0, 0.0),
        planar(0.0, 1.0, 0.0),
        planar(0.0, 0.0, 1.0),
        planar(0.36, -0.48, 0.8),
        planar(-0.48, 0.64, 0.6),
    ]
}

fn isotropic_limit_planars() -> Vec<Planar> {
    vec![
        planar(1.0, 0.0, 0.0),
        planar(0.0, 1.0, 0.0),
        planar(0.0, 0.0, 1.0),
    ]
}

const fn axis_index(axis: Axis) -> usize {
    match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
    }
}

fn displaced(point: &Point, axis: Axis, amount: f64) -> Point {
    let mut position = point.position();
    position[axis_index(axis)] += amount;
    Point::with_c(position[0], position[1], position[2], point.c()).unwrap()
}

fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let limit = absolute + relative * actual.abs().max(expected.abs());
    assert!(
        delta <= limit,
        "actual={actual:.17e}, expected={expected:.17e}, delta={delta:.3e}, limit={limit:.3e}"
    );
}

#[test]
fn construction_exposes_f32_anisotropy_evidence_and_rejects_fewer_than_two_planars() {
    assert_eq!(
        AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &[]),
        Err(AnisotropyError::InsufficientPlanars)
    );
    assert_eq!(
        AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &[planar(1.0, 0.0, 0.0)]),
        Err(AnisotropyError::InsufficientPlanars)
    );

    let kernel = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &oblique_planars()).unwrap();
    let eigenvalues = kernel.eigenvalues();
    let transform = kernel.transform();
    let plunge = kernel.global_plunge();
    assert_eq!(
        eigenvalues.map(f32::to_bits),
        [0x3f80_0000, 0x4000_0000, 0x4000_0000]
    );
    assert_eq!(
        transform.map(|row| row.map(f32::to_bits)),
        [
            [0x3f93_1644, 0xbe4b_9820, 0x0000_0000],
            [0xbe4b_9824, 0x3fa1_eeb0, 0x0000_0000],
            [0x0000_0000, 0x0000_0000, 0x3fb5_04f3],
        ]
    );
    assert_eq!(
        plunge.map(f32::to_bits),
        [0xbf4c_cccd, 0xbf19_999a, 0x8000_0000]
    );

    assert!(eigenvalues[0] <= eigenvalues[1] && eigenvalues[1] <= eigenvalues[2]);
    assert!(eigenvalues.into_iter().all(f32::is_finite));
    assert!(transform.into_iter().flatten().all(f32::is_finite));
    assert!(plunge.into_iter().all(f32::is_finite));
    for row in 0..3 {
        for column in 0..3 {
            assert!((transform[row][column] - transform[column][row]).abs() <= 2.0e-6);
        }
    }
    let plunge_norm = plunge
        .into_iter()
        .map(|component| component * component)
        .sum::<f32>()
        .sqrt();
    assert!((plunge_norm - 1.0).abs() <= 2.0e-6);
}

#[test]
fn isotropic_covariance_limit_produces_identity_support_matrix() {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let planars = isotropic_limit_planars();

    for kind in KERNELS {
        let kernel = AnisotropicKernel::new(kind, 0.7, &planars).unwrap();
        assert_eq!(
            kernel.transform(),
            [[1.0_f32, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_close(
            kernel.basis(&first, &second),
            georbf::IsotropicKernel::new(kind, 0.7).basis(&first, &second),
            1.0e-12,
            1.0e-11,
        );
    }
}

#[test]
fn eigenvalue_floor_is_applied_before_support_scaling() {
    let planars = vec![planar(1.0, 0.0, 0.0), planar(1.0, 0.0, 0.0)];
    let kernel = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &planars).unwrap();

    assert_eq!(kernel.eigenvalues(), [0.0, 0.0, 2.0]);
    assert_eq!(
        kernel.clamped_eigenvalues()[0].to_bits(),
        0.0001_f32.to_bits()
    );
    assert_eq!(
        kernel.clamped_eigenvalues()[1].to_bits(),
        0.0001_f32.to_bits()
    );
    assert_eq!(kernel.clamped_eigenvalues()[2], 2.0);
    assert!(kernel.transform().into_iter().flatten().all(f32::is_finite));
}

#[test]
fn eigenvalue_floor_boundary_and_uniform_normal_scale_are_deterministic() {
    let below = vec![
        planar(0.009_999, 0.0, 0.0),
        planar(0.0, 1.0, 0.0),
        planar(0.0, 0.0, 1.0),
    ];
    let above = vec![
        planar(0.010_001, 0.0, 0.0),
        planar(0.0, 1.0, 0.0),
        planar(0.0, 0.0, 1.0),
    ];
    let below_kernel = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &below).unwrap();
    let above_kernel = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &above).unwrap();
    assert!(below_kernel.eigenvalues()[0] < 0.0001);
    assert_eq!(below_kernel.clamped_eigenvalues()[0], 0.0001);
    assert!(above_kernel.eigenvalues()[0] > 0.0001);
    assert_eq!(
        above_kernel.clamped_eigenvalues()[0],
        above_kernel.eigenvalues()[0]
    );

    let scaled = oblique_planars()
        .into_iter()
        .map(|planar| {
            let [nx, ny, nz] = planar.normal();
            Planar::from_normal(0.0, 0.0, 0.0, nx * 0.5, ny * 0.5, nz * 0.5).unwrap()
        })
        .collect::<Vec<_>>();
    let original = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &oblique_planars()).unwrap();
    let scaled = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &scaled).unwrap();
    for row in 0..3 {
        for column in 0..3 {
            assert_close(
                scaled.transform()[row][column] as f64,
                original.transform()[row][column] as f64,
                2.0e-6,
                2.0e-5,
            );
        }
    }
}

#[test]
fn oblique_kernel_values_and_derivatives_match_frozen_cpp_golden() {
    let planars = oblique_planars();
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);

    // Filled from the ignored T13 frozen C++ probe. Each row is basis,
    // point-1 xyz, point-2 xyz, then the row-major mixed Hessian.
    const GOLDEN: [(RbfKernel, [u64; 16]); 5] = [
        (
            RbfKernel::Cubic,
            [
                0x403e_b0cb_0e38_df50,
                0xc03c_edf1_e075_1366,
                0x4034_91ea_1f3f_ab2c,
                0x4032_c90d_d250_b6c1,
                0x403c_edf1_e075_1366,
                0xc034_91ea_1f3f_ab2c,
                0xc032_c90d_d250_b6c1,
                0xc035_dd25_8a3a_2eaf,
                0x4025_f181_8009_60c2,
                0x4017_9c20_5e10_b7d8,
                0x4025_f181_8009_60c2,
                0xc033_ffe5_1594_225e,
                0xc010_c99b_d6ba_3394,
                0x4017_9c20_5e10_b7d8,
                0xc010_c99b_d6ba_3394,
                0xc036_9e3e_3556_bf5b,
            ],
        ),
        (
            RbfKernel::Gaussian,
            [
                0x3f80_ccdc_c761_bb1c,
                0x3f99_5adc_e17d_7134,
                0xbf92_0741_7396_e9ef,
                0xbf90_76d8_538a_14b0,
                0xbf99_5adc_e17d_7134,
                0x3f92_0741_7396_e9ef,
                0x3f90_76d8_538a_14b0,
                0xbfb0_557b_7286_0407,
                0x3fa9_3b90_4eb9_8f99,
                0x3fa8_d90b_9b82_d774,
                0x3fa9_3b90_4eb9_8f99,
                0xbf99_3116_d344_b965,
                0xbfa1_aaf3_4d88_25c7,
                0x3fa8_d90b_9b82_d774,
                0xbfa1_aaf3_4d88_25c7,
                0xbf8f_9c80_895a_1b30,
            ],
        ),
        (
            RbfKernel::Multiquadric,
            [
                0x4009_ed11_9f7b_8e19,
                0xbfee_69a4_1889_defd,
                0x3fe5_9fe4_ad3c_0f1c,
                0x3fe3_bf9c_2199_42d6,
                0x3fee_69a4_1889_defd,
                0xbfe5_9fe4_ad3c_0f1c,
                0xbfe3_bf9c_2199_42d6,
                0xbfc2_0a48_9a89_1e18,
                0xbfa9_a210_6144_e548,
                0xbfc7_2a81_0943_f23a,
                0xbfa9_a210_6144_e548,
                0xbfd7_5e7d_effe_47be,
                0x3fc0_78d1_9db5_42ba,
                0xbfc7_2a81_0943_f23a,
                0x3fc0_78d1_9db5_42ba,
                0xbfdff9bb4c19c4f3,
            ],
        ),
        (
            RbfKernel::ThinPlateSpline,
            [
                0x405b_6ac4_08ee_170c,
                0xc065_00ce_333e_d640,
                0x405d_de3c_3744_98c2,
                0x405b_46dc_d70c_5483,
                0x4065_00ce_333e_d640,
                0xc05d_de3c_3744_98c2,
                0xc05b_46dc_d70c_5483,
                0xc06b_3747_3051_6d0a,
                0x4060_07d1_bc4b_661a,
                0x4057_4d11_b0c6_4db3,
                0x4060_07d1_bc4b_661a,
                0xc064_4123_619e_64d2,
                0xc050_9165_5a64_bf23,
                0x4057_4d11_b0c6_4db3,
                0xc050_9165_5a64_bf23,
                0xc065_3424_3c66_3f19,
            ],
        ),
        (
            RbfKernel::InverseMultiquadric,
            [
                0x3fd3_bf9c_2cf0_7e6f,
                0x3fb7_2a81_1691_986e,
                0xbfb0_78d1_a72a_cf41,
                0xbfae_15f3_aa66_e9cc,
                0xbfb7_2a81_1691_986e,
                0x3fb0_78d1_a72a_cf41,
                0x3fae_15f3_aa66_e9cc,
                0xbfa4_4de0_16cd_d49a,
                0x3fa5_c35d_d477_34ea,
                0x3faa_781a_d2e2_48f4,
                0x3fa5_c35d_d477_34ea,
                0x3f80_3f3c_dbaa_9724,
                0xbfa2_d212_e3e9_c1d6,
                0x3faa_781a_d2e2_48f4,
                0xbfa2_d212_e3e9_c1d6,
                0x3f99_cbc4_e872_9bf6,
            ],
        ),
    ];

    for (kind, expected) in GOLDEN {
        let evaluation = AnisotropicKernel::new(kind, 0.7, &planars)
            .unwrap()
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
}

#[test]
fn analytic_anisotropic_derivatives_match_finite_differences() {
    let planars = oblique_planars();
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let step = 1.0e-5;

    for kind in DIFFERENTIABLE_KERNELS {
        let kernel = AnisotropicKernel::new(kind, 0.7, &planars).unwrap();
        for axis in AXES {
            let lower = displaced(&first, axis, -step);
            let upper = displaced(&first, axis, step);
            let numeric =
                (kernel.basis(&upper, &second) - kernel.basis(&lower, &second)) / (2.0 * step);
            let analytic = kernel
                .first_derivative(&first, &second, DerivativePoint::First, axis)
                .unwrap();
            assert_close(analytic, numeric, 1.0e-8, 3.0e-7);
        }

        for first_axis in AXES {
            for second_axis in AXES {
                let lower = displaced(&second, second_axis, -step);
                let upper = displaced(&second, second_axis, step);
                let lower_value = kernel
                    .first_derivative(&first, &lower, DerivativePoint::First, first_axis)
                    .unwrap();
                let upper_value = kernel
                    .first_derivative(&first, &upper, DerivativePoint::First, first_axis)
                    .unwrap();
                let numeric = (upper_value - lower_value) / (2.0 * step);
                let analytic = kernel
                    .mixed_second_derivative(&first, &second, first_axis, second_axis)
                    .unwrap();
                assert_close(analytic, numeric, 1.0e-7, 2.0e-6);
            }
        }
    }
}

#[test]
fn linear_anisotropic_kernel_retains_value_and_typed_derivative_error() {
    let kernel = AnisotropicKernel::new(RbfKernel::Linear, 0.7, &oblique_planars()).unwrap();
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);

    assert_eq!(
        kernel.basis(&first, &second).to_bits(),
        0x4009_0c12_7c23_2568
    );
    assert_eq!(
        kernel.first_derivative(&first, &second, DerivativePoint::First, Axis::X),
        Err(KernelError::LinearDerivativeUnavailable)
    );
}

#[test]
fn anisotropic_scaled_radius_ignores_fourth_coordinate() {
    let planars = oblique_planars();
    let kernel = AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &planars).unwrap();
    let first = Point::with_c(-1.25, 0.5, 2.0, 20.0).unwrap();
    let second = Point::with_c(0.75, -0.25, 1.0, -30.0).unwrap();
    let without_c_first = point(-1.25, 0.5, 2.0);
    let without_c_second = point(0.75, -0.25, 1.0);

    assert_eq!(
        kernel.basis(&first, &second).to_bits(),
        kernel.basis(&without_c_first, &without_c_second).to_bits()
    );
}

#[test]
fn unsupported_anisotropic_variants_are_rejected_without_fabrication() {
    for kind in [
        RbfKernel::MultiquadricCubic,
        RbfKernel::WendlandC2,
        RbfKernel::MaternC4,
    ] {
        assert_eq!(
            AnisotropicKernel::new(kind, 0.7, &oblique_planars()),
            Err(AnisotropyError::UnsupportedKernel(kind))
        );
    }

    let huge = planar(f64::MAX, 0.0, 0.0);
    assert_eq!(
        AnisotropicKernel::new(RbfKernel::Cubic, 0.7, &[huge.clone(), huge]),
        Err(AnisotropyError::NonFiniteComputation)
    );
}
