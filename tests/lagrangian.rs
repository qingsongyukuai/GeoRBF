use georbf::{Error, Interface, LagrangianPolynomialBasis, Point, POSITION_EPSILON};

const TOLERANCE: f64 = 2.0e-14;

fn interface(x: f64, y: f64, z: f64, level: f64) -> Interface {
    Interface::new(x, y, z, level).unwrap()
}

fn point(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z).unwrap()
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.1e}"
    );
}

fn selected_horizon() -> Vec<Interface> {
    vec![
        interface(0.0, 0.0, 0.0, 20.0),
        interface(10.0, 1.0, 1.0, 20.0),
        interface(2.0, -3.0, 0.0, 20.0),
        interface(3.0, 4.0, 2.0, 20.0),
        interface(4.0, 2.0, 0.5, 20.0),
        interface(5.0, -1.0, 1.5, 20.0),
    ]
}

#[test]
fn largest_horizon_selection_and_coefficients_match_frozen_surfe() {
    let tied_group = selected_horizon()[..4].to_vec();
    let tied = LagrangianPolynomialBasis::new(&[tied_group.clone(), tied_group]).unwrap();
    assert_eq!(tied.selected_horizon_index(), 0);

    let groups = vec![
        vec![
            interface(-1.0, -1.0, -1.0, 30.0),
            interface(1.0, -1.0, -1.0, 30.0),
            interface(-1.0, 1.0, -1.0, 30.0),
            interface(-1.0, -1.0, 1.0, 30.0),
        ],
        selected_horizon(),
    ];
    let basis = LagrangianPolynomialBasis::new(&groups).unwrap();

    assert_eq!(basis.selected_horizon_index(), 1);
    assert_eq!(basis.selected_source_indices(), [0, 1, 2, 3]);
    assert_eq!(
        std::array::from_fn(|index| basis.unisolvent_points()[index].position()),
        [
            [0.0, 0.0, 0.0],
            [10.0, 1.0, 1.0],
            [2.0, -3.0, 0.0],
            [3.0, 4.0, 2.0],
        ]
    );

    let denominator = 47.0;
    let expected = [
        [
            1.0,
            -1.0 / denominator,
            15.0 / denominator,
            -52.0 / denominator,
        ],
        [
            0.0,
            6.0 / denominator,
            4.0 / denominator,
            -17.0 / denominator,
        ],
        [
            0.0,
            -2.0 / denominator,
            -17.0 / denominator,
            37.0 / denominator,
        ],
        [
            0.0,
            -3.0 / denominator,
            -2.0 / denominator,
            32.0 / denominator,
        ],
    ];
    for (actual, expected) in basis
        .coefficients()
        .iter()
        .flatten()
        .zip(expected.iter().flatten())
    {
        assert_close(*actual, *expected, TOLERANCE);
    }

    let sample = point(1.25, -0.75, 2.5);
    let expected_values = [
        1.0 + (-1.25 - 11.25 - 130.0) / denominator,
        (7.5 - 3.0 - 42.5) / denominator,
        (-2.5 + 12.75 + 92.5) / denominator,
        (-3.75 + 1.5 + 80.0) / denominator,
    ];
    for (actual, expected) in basis.values(&sample).into_iter().zip(expected_values) {
        assert_close(actual, expected, TOLERANCE);
    }

    // Exact binary64 output from the frozen T11 C++ probe.
    assert_eq!(
        basis.values(&sample).map(f64::to_bits),
        [
            0xc000_415c_9882_b932,
            0xbfe9_df51_b3be_a368,
            0x4001_7d46_cefa_8d9e,
            0x3ffa_77d4_6cef_a8da,
        ]
    );
    assert_eq!(
        basis.dx(&sample).map(f64::to_bits),
        [
            0xbf95_c988_2b93_1057,
            0x3fc0_5726_20ae_4c41,
            0xbfa5_c988_2b93_1057,
            0xbfb0_5726_20ae_4c41,
        ]
    );
    assert_eq!(
        basis.dy(&sample).map(f64::to_bits),
        [
            0x3fd4_6cef_a8d9_df52,
            0x3fb5_c988_2b93_1057,
            0xbfd7_2620_ae4c_415d,
            0xbfa5_c988_2b93_1057,
        ]
    );
    assert_eq!(
        basis.dz(&sample).map(f64::to_bits),
        [
            0xbff1_b3be_a367_7d47,
            0xbfd7_2620_ae4c_415d,
            0x3fe9_3105_7262_0ae5,
            0x3fe5_c988_2b93_1057,
        ]
    );
}

#[test]
fn selected_points_satisfy_kronecker_identity_and_derivatives() {
    let basis = LagrangianPolynomialBasis::new(&[selected_horizon()]).unwrap();
    for (selected_index, selected) in basis.unisolvent_points().iter().enumerate() {
        for (basis_index, value) in basis.values(selected).into_iter().enumerate() {
            assert_close(
                value,
                usize::from(selected_index == basis_index) as f64,
                TOLERANCE,
            );
        }
    }

    let sample = [1.25, -0.75, 2.5];
    let center = point(sample[0], sample[1], sample[2]);
    let analytic = [basis.dx(&center), basis.dy(&center), basis.dz(&center)];
    let step = 1.0e-6;
    for axis in 0..3 {
        let mut lower = sample;
        let mut upper = sample;
        lower[axis] -= step;
        upper[axis] += step;
        let lower = basis.values(&point(lower[0], lower[1], lower[2]));
        let upper = basis.values(&point(upper[0], upper[1], upper[2]));
        for term in 0..4 {
            assert_close(
                analytic[axis][term],
                (upper[term] - lower[term]) / (2.0 * step),
                2.0e-9,
            );
        }
    }
}

#[test]
fn point_permutation_preserves_geometry_and_basis_with_remapped_indices() {
    let forward = LagrangianPolynomialBasis::new(&[selected_horizon()]).unwrap();
    let mut reversed_horizon = selected_horizon();
    reversed_horizon.reverse();
    let reversed = LagrangianPolynomialBasis::new(&[reversed_horizon]).unwrap();

    assert_eq!(reversed.selected_source_indices(), [5, 4, 3, 2]);
    assert_eq!(
        std::array::from_fn::<_, 4, _>(|index| { forward.unisolvent_points()[index].position() }),
        std::array::from_fn::<_, 4, _>(|index| { reversed.unisolvent_points()[index].position() })
    );
    let sample = point(-4.0, 2.25, 0.125);
    assert_eq!(forward.values(&sample), reversed.values(&sample));
    assert_eq!(forward.dx(&sample), reversed.dx(&sample));
    assert_eq!(forward.dy(&sample), reversed.dy(&sample));
    assert_eq!(forward.dz(&sample), reversed.dz(&sample));
}

#[test]
fn axis_aligned_plane_uses_frozen_epsilon_adjustment() {
    let plane = vec![
        interface(2.0, 0.0, 0.0, 7.0),
        interface(2.0, 10.0, 1.0, 7.0),
        interface(2.0, 2.0, -3.0, 7.0),
        interface(2.0, 4.0, 4.0, 7.0),
    ];
    let basis = LagrangianPolynomialBasis::new(&[plane]).unwrap();

    assert_eq!(basis.selected_source_indices(), [0, 1, 2, 3]);
    assert_eq!(basis.unisolvent_points()[0].x(), 2.0 + POSITION_EPSILON);
    assert_eq!(basis.unisolvent_points()[1].x(), 2.0);
    for (selected_index, selected) in basis.unisolvent_points().iter().enumerate() {
        for (basis_index, value) in basis.values(selected).into_iter().enumerate() {
            assert_close(
                value,
                usize::from(selected_index == basis_index) as f64,
                5.0e-12,
            );
        }
    }
}

#[test]
fn empty_insufficient_tied_axes_and_degenerate_subsets_fail_safely() {
    assert_eq!(
        LagrangianPolynomialBasis::new(&[]).unwrap_err(),
        Error::LagrangianBasisCreationFailure
    );
    assert_eq!(
        LagrangianPolynomialBasis::new(&[vec![
            interface(0.0, 0.0, 0.0, 1.0),
            interface(1.0, 0.0, 0.0, 1.0),
            interface(0.0, 1.0, 0.0, 1.0),
        ]])
        .unwrap_err(),
        Error::LagrangianBasisCreationFailure
    );

    let equal_axis_ranges = vec![
        interface(0.0, 0.0, 0.0, 1.0),
        interface(1.0, 0.0, 0.0, 1.0),
        interface(0.0, 1.0, 0.0, 1.0),
        interface(0.0, 0.0, 1.0, 1.0),
    ];
    assert_eq!(
        LagrangianPolynomialBasis::new(&[equal_axis_ranges]).unwrap_err(),
        Error::LagrangianBasisCreationFailure
    );

    let oblique_coplanar = vec![
        interface(0.0, 0.0, 0.0, 1.0),
        interface(10.0, 1.0, 11.0, 1.0),
        interface(2.0, -3.0, -1.0, 1.0),
        interface(3.0, 4.0, 7.0, 1.0),
        interface(4.0, 2.0, 6.0, 1.0),
    ];
    assert_eq!(
        LagrangianPolynomialBasis::new(&[oblique_coplanar]).unwrap_err(),
        Error::LagrangianBasisCreationFailure
    );

    let collinear = vec![
        interface(0.0, 0.0, 0.0, 1.0),
        interface(12.0, 6.0, 3.0, 1.0),
        interface(2.0, 1.0, 0.5, 1.0),
        interface(4.0, 2.0, 1.0, 1.0),
    ];
    assert_eq!(
        LagrangianPolynomialBasis::new(&[collinear]).unwrap_err(),
        Error::LagrangianBasisCreationFailure
    );

    let repeated = vec![
        interface(0.0, 0.0, 0.0, 1.0),
        interface(10.0, 1.0, 1.0, 1.0),
        interface(2.0, -3.0, 0.0, 1.0),
        interface(2.0, -3.0, 0.0, 1.0),
    ];
    assert_eq!(
        LagrangianPolynomialBasis::new(&[repeated]).unwrap_err(),
        Error::LagrangianBasisCreationFailure
    );
}

#[test]
fn nonzero_near_degenerate_determinant_is_attempted_not_pre_rejected() {
    let almost_coplanar = vec![
        interface(0.0, 0.0, 0.0, 1.0),
        interface(10.0, 1.0, 0.0, 1.0),
        interface(2.0, -3.0, 0.0, 1.0),
        interface(3.0, 4.0, 1.0e-12, 1.0),
    ];
    let basis = LagrangianPolynomialBasis::new(&[almost_coplanar]).unwrap();
    assert!(basis
        .coefficients()
        .iter()
        .flatten()
        .all(|value| value.is_finite()));
    for (selected_index, selected) in basis.unisolvent_points().iter().enumerate() {
        for (basis_index, value) in basis.values(selected).into_iter().enumerate() {
            assert_close(
                value,
                usize::from(selected_index == basis_index) as f64,
                1.0e-8,
            );
        }
    }

    let sample = point(1.25, -0.75, 2.5);
    assert_eq!(
        basis.values(&sample).map(f64::to_bits),
        [
            0xc28d_8efe_f487_fac0,
            0xc273_53a6_b393_fee0,
            0x4285_0835_6912_0230,
            0x4282_309c_e540_0000,
        ]
    );
    assert_eq!(
        basis.dz(&sample).map(f64::to_bits),
        [
            0xc277_a598_c3a0_0000,
            0xc25e_ec3d_ec20_0000,
            0x4270_d35d_eda8_0000,
            0x426d_1a94_a200_0000,
        ]
    );
}
