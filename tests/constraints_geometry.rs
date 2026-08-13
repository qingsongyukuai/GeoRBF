use std::f64::consts::FRAC_1_SQRT_2;

use georbf::{
    ConstraintError, Constraints, Inequality, Interface, Planar, Point, Polarity, Tangent,
};

const EPSILON: f64 = 2.0e-14;

fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= EPSILON * scale,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn assert_vector_close(actual: [f64; 3], expected: [f64; 3]) {
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert_close(actual, expected);
    }
}

#[test]
fn point_state_is_owned_initialized_and_finite() {
    let mut point = Point::with_c(1.0, -2.0, 3.5, 4.25).unwrap();
    assert_eq!(point.position(), [1.0, -2.0, 3.5]);
    assert_eq!(
        (point.x(), point.y(), point.z(), point.c()),
        (1.0, -2.0, 3.5, 4.25)
    );
    assert_eq!(point.scalar_field(), 0.0);
    assert_eq!(point.vector_field(), [0.0; 3]);

    point.set_c(-7.0).unwrap();
    point.set_scalar_field(8.5).unwrap();
    point.set_vector_field(-1.0, 2.0, -3.0).unwrap();
    assert_eq!(point.c(), -7.0);
    assert_eq!(point.scalar_field(), 8.5);
    assert_eq!(point.vector_field(), [-1.0, 2.0, -3.0]);
    assert_eq!(
        (point.nx_interp(), point.ny_interp(), point.nz_interp()),
        (-1.0, 2.0, -3.0)
    );

    assert_eq!(
        Point::new(f64::NAN, 0.0, 0.0).unwrap_err(),
        ConstraintError::NonFiniteInput
    );
    assert_eq!(
        point.set_c(f64::INFINITY),
        Err(ConstraintError::NonFiniteInput)
    );
    assert_eq!(
        point.set_scalar_field(f64::NEG_INFINITY),
        Err(ConstraintError::NonFiniteInput)
    );
    assert_eq!(
        point.set_vector_field(0.0, f64::NAN, 0.0),
        Err(ConstraintError::NonFiniteInput)
    );
}

#[test]
fn scalar_constraint_types_match_surfe_constructor_and_bound_semantics() {
    let mut interface = Interface::with_c(1.0, 2.0, 3.0, -4.0, 5.0).unwrap();
    assert_eq!(interface.point().position(), [1.0, 2.0, 3.0]);
    assert_eq!(interface.point().c(), 5.0);
    assert_eq!(interface.level(), -4.0);
    assert_eq!(interface.level_bounds()[0].to_bits(), 0);
    assert_eq!(interface.level_bounds()[1].to_bits(), 0);
    interface.set_level_bounds(0.0).unwrap();
    assert_eq!(
        interface.level_lower_bound().to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(interface.level_upper_bound().to_bits(), 0);
    interface.set_level_bounds(0.75).unwrap();
    assert_eq!(interface.level_bounds(), [-0.75, 0.75]);

    let inequality = Inequality::with_c(-1.0, -2.0, -3.0, 9.0, 4.0).unwrap();
    assert_eq!(inequality.point().position(), [-1.0, -2.0, -3.0]);
    assert_eq!(inequality.point().c(), 4.0);
    assert_eq!(inequality.level(), 9.0);

    assert_eq!(
        Interface::new(0.0, 0.0, 0.0, f64::NAN).unwrap_err(),
        ConstraintError::NonFiniteInput
    );
    assert_eq!(
        interface.set_level_bounds(f64::INFINITY),
        Err(ConstraintError::NonFiniteInput)
    );
    assert_eq!(
        Inequality::new(0.0, 0.0, 0.0, f64::NEG_INFINITY).unwrap_err(),
        ConstraintError::NonFiniteInput
    );
}

#[test]
fn normal_to_geological_angles_matches_frozen_surfe_quadrants() {
    let cases = [
        ([1.0, 0.0, 0.0], 90.0, 360.0, Polarity::Upright),
        ([0.0, 1.0, 0.0], 90.0, 270.0, Polarity::Upright),
        ([-1.0, 0.0, 0.0], 90.0, 180.0, Polarity::Upright),
        ([0.0, -1.0, 0.0], 90.0, 90.0, Polarity::Upright),
        ([0.0, 0.0, 1.0], 0.0, 360.0, Polarity::Upright),
        ([0.0, 0.0, -1.0], 180.0, 360.0, Polarity::Overturned),
    ];

    for (normal, expected_dip, expected_strike, expected_polarity) in cases {
        let planar =
            Planar::from_normal(10.0, 20.0, 30.0, normal[0], normal[1], normal[2]).unwrap();
        assert_eq!(planar.normal(), normal);
        assert_eq!([planar.nx(), planar.ny(), planar.nz()], normal);
        assert_close(planar.dip(), expected_dip);
        assert_close(planar.strike(), expected_strike);
        assert_eq!(planar.polarity(), expected_polarity);
    }

    let not_normalized = Planar::from_normal(0.0, 0.0, 0.0, 2.0, 0.0, 0.0).unwrap();
    assert_eq!(not_normalized.normal(), [2.0, 0.0, 0.0]);
    assert_close(not_normalized.dip(), 90.0);
    assert_close(not_normalized.strike(), 360.0);
}

#[test]
fn strike_dip_polarity_and_azimuth_constructors_match_surfe() {
    let upright =
        Planar::from_strike_dip_polarity(1.0, 2.0, 3.0, 30.0, 40.0, Polarity::Upright).unwrap();
    assert_vector_close(
        upright.normal(),
        [
            0.556_670_399_226_419_4,
            -0.321_393_804_843_269_6,
            0.766_044_443_118_978,
        ],
    );
    assert_close(upright.strike(), 30.0);
    assert_close(upright.dip(), 40.0);

    let overturned =
        Planar::from_strike_dip_polarity(1.0, 2.0, 3.0, 30.0, 40.0, Polarity::Overturned).unwrap();
    assert_vector_close(
        overturned.normal(),
        [
            -0.556_670_399_226_419_4,
            0.321_393_804_843_269_6,
            -0.766_044_443_118_978,
        ],
    );

    let from_azimuth =
        Planar::from_azimuth_dip_polarity(1.0, 2.0, 3.0, 120.0, 40.0, Polarity::Upright).unwrap();
    assert_close(from_azimuth.strike(), 30.0);
    assert_vector_close(from_azimuth.normal(), upright.normal());

    let wrapped_branch =
        Planar::from_azimuth_dip_polarity(0.0, 0.0, 0.0, 0.0, 45.0, Polarity::Upright).unwrap();
    assert_close(wrapped_branch.strike(), 270.0);
    let branch_boundary =
        Planar::from_azimuth_dip_polarity(0.0, 0.0, 0.0, 90.0, 45.0, Polarity::Upright).unwrap();
    assert_close(branch_boundary.strike(), 0.0);
    let below_boundary =
        Planar::from_azimuth_dip_polarity(0.0, 0.0, 0.0, 89.999, 45.0, Polarity::Upright).unwrap();
    assert_close(below_boundary.strike(), 359.999);

    for code in [i32::MIN, -1, 2, i32::MAX] {
        assert_eq!(
            Polarity::try_from(code),
            Err(ConstraintError::InvalidPolarity)
        );
    }
    assert_eq!(Polarity::try_from(0), Ok(Polarity::Upright));
    assert_eq!(Polarity::try_from(1), Ok(Polarity::Overturned));
}

#[test]
fn planar_round_trip_and_direction_vectors_satisfy_geometry_identities() {
    for (strike, dip, polarity) in [
        (0.0, 0.0, Polarity::Upright),
        (30.0, 40.0, Polarity::Upright),
        (30.0, 40.0, Polarity::Overturned),
        (270.0, 89.0, Polarity::Upright),
    ] {
        let from_angles =
            Planar::from_strike_dip_polarity(0.0, 0.0, 0.0, strike, dip, polarity).unwrap();
        let normal = from_angles.normal();
        let from_normal =
            Planar::from_normal(0.0, 0.0, 0.0, normal[0], normal[1], normal[2]).unwrap();
        let reconstructed = Planar::from_strike_dip_polarity(
            0.0,
            0.0,
            0.0,
            from_normal.strike(),
            from_normal.dip(),
            from_normal.polarity(),
        )
        .unwrap();
        assert_vector_close(reconstructed.normal(), normal);
    }

    let planar =
        Planar::from_strike_dip_polarity(0.0, 0.0, 0.0, 30.0, 40.0, Polarity::Upright).unwrap();
    let dip = planar.dip_vector();
    let strike = planar.strike_vector();
    assert_vector_close(
        dip,
        [
            0.663_413_948_168_938_4,
            -0.383_022_221_559_489,
            -0.642_787_609_686_539_3,
        ],
    );
    assert_vector_close(strike, [0.5, 0.866_025_403_784_438_6, 0.0]);
    assert_close(dip.into_iter().map(|value| value * value).sum::<f64>(), 1.0);
    assert_close(
        strike.into_iter().map(|value| value * value).sum::<f64>(),
        1.0,
    );
    assert_close(
        dip.into_iter()
            .zip(strike)
            .map(|(left, right)| left * right)
            .sum::<f64>(),
        0.0,
    );
}

#[test]
fn planar_normal_bounds_match_frozen_surfe_corner_envelope() {
    let mut planar =
        Planar::from_strike_dip_polarity(0.0, 0.0, 0.0, 30.0, 40.0, Polarity::Upright).unwrap();
    assert_eq!(planar.normal_bounds(), None);
    assert_eq!(planar.nx_bounds(), None);
    assert_eq!(planar.ny_bounds(), None);
    assert_eq!(planar.nz_bounds(), None);
    planar.set_normal_bounds(10.0, 5.0).unwrap();
    let bounds = planar.normal_bounds().unwrap();
    assert_eq!(planar.nx_bounds(), Some(bounds[0]));
    assert_eq!(planar.ny_bounds(), Some(bounds[1]));
    assert_eq!(planar.nz_bounds(), Some(bounds[2]));
    assert_vector_close(
        [bounds[0][0], bounds[1][0], bounds[2][0]],
        [
            0.439_385_041_770_705,
            -0.454_519_477_672_043_6,
            FRAC_1_SQRT_2,
        ],
    );
    assert_vector_close(
        [bounds[0][1], bounds[1][1], bounds[2][1]],
        [
            0.664_463_024_388_674_7,
            -0.196_174_694_969_011_05,
            0.819_152_044_288_991_8,
        ],
    );
}

#[test]
fn tangent_defaults_bounds_and_validation_match_surfe() {
    let mut tangent = Tangent::with_c(1.0, 2.0, 3.0, 1.0, -2.0, 3.0, 4.0).unwrap();
    assert_eq!(tangent.point().position(), [1.0, 2.0, 3.0]);
    assert_eq!(tangent.point().c(), 4.0);
    assert_eq!(tangent.vector(), [1.0, -2.0, 3.0]);
    assert_eq!(tangent.inner_product_constraint(), 0.0);
    assert_eq!(tangent.angle_bounds(), None);
    assert_eq!(tangent.angle_lower_bound(), None);
    assert_eq!(tangent.angle_upper_bound(), None);

    tangent.set_angle_bounds(10.0).unwrap();
    let bounds = tangent.angle_bounds().unwrap();
    assert_eq!(tangent.angle_lower_bound(), Some(bounds[0]));
    assert_eq!(tangent.angle_upper_bound(), Some(bounds[1]));
    assert_close(bounds[0], 0.0);
    assert_close(bounds[1], 0.347_296_355_333_860_83);

    tangent.set_angle_bounds(-10.0).unwrap();
    let negative_bounds = tangent.angle_bounds().unwrap();
    assert_close(negative_bounds[0], -0.347_296_355_333_860_6);
    assert_close(negative_bounds[1], 0.0);

    assert_eq!(
        Tangent::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap_err(),
        ConstraintError::ZeroTangent
    );
    assert_eq!(
        Tangent::new(0.0, 0.0, 0.0, f64::INFINITY, 0.0, 0.0).unwrap_err(),
        ConstraintError::NonFiniteInput
    );
    assert_eq!(
        tangent.set_angle_bounds(f64::NAN),
        Err(ConstraintError::NonFiniteInput)
    );
}

#[test]
fn invalid_planar_directions_fail_safely_instead_of_copying_cpp_nan_or_ub() {
    assert_eq!(
        Planar::from_normal(0.0, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap_err(),
        ConstraintError::ZeroNormal
    );
    assert_eq!(
        Planar::from_normal(0.0, 0.0, 0.0, f64::NAN, 0.0, 1.0).unwrap_err(),
        ConstraintError::NonFiniteInput
    );
    assert_eq!(
        Planar::from_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.000_000_000_000_000_2).unwrap_err(),
        ConstraintError::NormalZOutOfRange
    );
    assert_eq!(
        Planar::from_strike_dip_polarity(0.0, 0.0, 0.0, f64::INFINITY, 45.0, Polarity::Upright,)
            .unwrap_err(),
        ConstraintError::NonFiniteInput
    );
    assert_eq!(
        Planar::from_azimuth_dip_polarity(0.0, 0.0, 0.0, f64::NAN, 45.0, Polarity::Upright,)
            .unwrap_err(),
        ConstraintError::NonFiniteInput
    );
}

#[test]
fn constraints_container_preserves_category_and_insertion_order() {
    let mut constraints = Constraints::default();
    constraints
        .interfaces
        .push(Interface::new(1.0, 0.0, 0.0, 1.0).unwrap());
    constraints
        .interfaces
        .push(Interface::new(2.0, 0.0, 0.0, 2.0).unwrap());
    constraints
        .inequalities
        .push(Inequality::new(3.0, 0.0, 0.0, 3.0).unwrap());
    constraints
        .planars
        .push(Planar::from_normal(4.0, 0.0, 0.0, 0.0, 0.0, 1.0).unwrap());
    constraints
        .tangents
        .push(Tangent::new(5.0, 0.0, 0.0, 1.0, 0.0, 0.0).unwrap());

    assert_eq!(constraints.interfaces.len(), 2);
    assert_eq!(constraints.interfaces[0].level(), 1.0);
    assert_eq!(constraints.interfaces[1].level(), 2.0);
    assert_eq!(constraints.inequalities.len(), 1);
    assert_eq!(constraints.planars.len(), 1);
    assert_eq!(constraints.tangents.len(), 1);
}
