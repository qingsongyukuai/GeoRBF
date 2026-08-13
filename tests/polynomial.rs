use georbf::{Point, PolynomialBasis, PolynomialOrder};

const FINITE_DIFFERENCE_STEP: f64 = 1.0e-6;
const FINITE_DIFFERENCE_TOLERANCE: f64 = 2.0e-9;

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

fn assert_vector_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert_close(actual, expected, tolerance);
    }
}

#[test]
fn complete_basis_term_order_and_derivatives_match_frozen_surfe() {
    let sample = point(-2.0, 3.0, -4.0);

    let zero = PolynomialBasis::complete(PolynomialOrder::Zero);
    assert_eq!(zero.term_count(), 1);
    assert_eq!(zero.values(&sample), [1.0]);
    assert_eq!(zero.dx(&sample), [0.0]);
    assert_eq!(zero.dy(&sample), [0.0]);
    assert_eq!(zero.dz(&sample), [0.0]);

    let first = PolynomialBasis::complete(PolynomialOrder::First);
    assert_eq!(first.term_count(), 4);
    assert_eq!(first.values(&sample), [-2.0, 3.0, -4.0, 1.0]);
    assert_eq!(first.dx(&sample), [1.0, 0.0, 0.0, 0.0]);
    assert_eq!(first.dy(&sample), [0.0, 1.0, 0.0, 0.0]);
    assert_eq!(first.dz(&sample), [0.0, 0.0, 1.0, 0.0]);

    let second = PolynomialBasis::complete(PolynomialOrder::Second);
    assert_eq!(second.term_count(), 10);
    assert_eq!(
        second.values(&sample),
        [4.0, 9.0, 16.0, -6.0, 8.0, -12.0, -2.0, 3.0, -4.0, 1.0]
    );
    assert_eq!(
        second.dx(&sample),
        [-4.0, 0.0, 0.0, 3.0, -4.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        second.dy(&sample),
        [0.0, 6.0, 0.0, -2.0, 0.0, -4.0, 0.0, 1.0, 0.0, 0.0]
    );
    assert_eq!(
        second.dz(&sample),
        [0.0, 0.0, -8.0, 0.0, -2.0, 3.0, 0.0, 0.0, 1.0, 0.0]
    );
}

#[test]
fn truncation_removes_only_the_constant_term_in_surfe_order() {
    let sample = point(-2.0, 3.0, -4.0);

    let zero = PolynomialBasis::truncated(PolynomialOrder::Zero);
    assert!(zero.is_truncated());
    assert_eq!(zero.term_count(), 0);
    assert!(zero.values(&sample).is_empty());
    assert!(zero.dx(&sample).is_empty());
    assert!(zero.dy(&sample).is_empty());
    assert!(zero.dz(&sample).is_empty());

    let first = PolynomialBasis::truncated(PolynomialOrder::First);
    assert_eq!(first.term_count(), 3);
    assert_eq!(first.values(&sample), [-2.0, 3.0, -4.0]);
    assert_eq!(first.dx(&sample), [1.0, 0.0, 0.0]);
    assert_eq!(first.dy(&sample), [0.0, 1.0, 0.0]);
    assert_eq!(first.dz(&sample), [0.0, 0.0, 1.0]);

    let second = PolynomialBasis::truncated(PolynomialOrder::Second);
    assert_eq!(second.term_count(), 9);
    assert_eq!(
        second.values(&sample),
        [4.0, 9.0, 16.0, -6.0, 8.0, -12.0, -2.0, 3.0, -4.0]
    );
    assert_eq!(
        second.dx(&sample),
        [-4.0, 0.0, 0.0, 3.0, -4.0, 0.0, 1.0, 0.0, 0.0]
    );
    assert_eq!(
        second.dy(&sample),
        [0.0, 6.0, 0.0, -2.0, 0.0, -4.0, 0.0, 1.0, 0.0]
    );
    assert_eq!(
        second.dz(&sample),
        [0.0, 0.0, -8.0, 0.0, -2.0, 3.0, 0.0, 0.0, 1.0]
    );
}

#[test]
fn origin_multiscale_and_fourth_coordinate_do_not_change_term_layout() {
    let second = PolynomialBasis::complete(PolynomialOrder::Second);
    assert_eq!(
        second.values(&point(0.0, -0.0, 0.0)),
        [0.0, 0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, 1.0]
    );

    let multiscale = Point::with_c(1.0e6, -2.0e-6, 3.0, 99.0).unwrap();
    let without_c = point(1.0e6, -2.0e-6, 3.0);
    assert_eq!(second.values(&multiscale), second.values(&without_c));
    assert_eq!(
        second.values(&multiscale),
        [1.0e12, 4.0e-12, 9.0, -2.0, 3.0e6, -6.0e-6, 1.0e6, -2.0e-6, 3.0, 1.0,]
    );
}

#[test]
fn analytic_derivatives_match_central_finite_differences() {
    let sample = [0.75, -1.25, 2.5];
    for order in [
        PolynomialOrder::Zero,
        PolynomialOrder::First,
        PolynomialOrder::Second,
    ] {
        for basis in [
            PolynomialBasis::complete(order),
            PolynomialBasis::truncated(order),
        ] {
            let center = point(sample[0], sample[1], sample[2]);
            let analytic = [basis.dx(&center), basis.dy(&center), basis.dz(&center)];

            for axis in 0..3 {
                let mut lower = sample;
                let mut upper = sample;
                lower[axis] -= FINITE_DIFFERENCE_STEP;
                upper[axis] += FINITE_DIFFERENCE_STEP;
                let lower_values = basis.values(&point(lower[0], lower[1], lower[2]));
                let upper_values = basis.values(&point(upper[0], upper[1], upper[2]));
                let finite_difference = upper_values
                    .iter()
                    .zip(lower_values)
                    .map(|(upper, lower)| (upper - lower) / (2.0 * FINITE_DIFFERENCE_STEP))
                    .collect::<Vec<_>>();
                assert_vector_close(
                    &analytic[axis],
                    &finite_difference,
                    FINITE_DIFFERENCE_TOLERANCE,
                );
            }
        }
    }
}
