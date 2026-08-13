use crate::{Axis, Point, SecondDerivative};

pub(crate) const fn axis_index(axis: Axis) -> usize {
    match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
    }
}

pub(crate) const fn second_derivative_axes(component: SecondDerivative) -> (Axis, Axis) {
    match component {
        SecondDerivative::DxDx => (Axis::X, Axis::X),
        SecondDerivative::DxDy => (Axis::X, Axis::Y),
        SecondDerivative::DxDz => (Axis::X, Axis::Z),
        SecondDerivative::DyDx => (Axis::Y, Axis::X),
        SecondDerivative::DyDy => (Axis::Y, Axis::Y),
        SecondDerivative::DyDz => (Axis::Y, Axis::Z),
        SecondDerivative::DzDx => (Axis::Z, Axis::X),
        SecondDerivative::DzDy => (Axis::Z, Axis::Y),
        SecondDerivative::DzDz => (Axis::Z, Axis::Z),
    }
}

pub(crate) fn radius_and_deltas(first: &Point, second: &Point) -> (f64, [f64; 3]) {
    let x_delta = first.x() - second.x();
    let y_delta = first.y() - second.y();
    let z_delta = first.z() - second.z();
    let c_delta = first.c() - second.c();
    let radius =
        (x_delta * x_delta + y_delta * y_delta + z_delta * z_delta + c_delta * c_delta).sqrt();
    (radius, [x_delta, y_delta, z_delta])
}
