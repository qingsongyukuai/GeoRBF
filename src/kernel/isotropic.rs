//! Frozen Surfe isotropic radial kernels and their spatial derivatives.
//!
//! Source: `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`.
//! The fourth point coordinate participates in the radius exactly as it does
//! in Surfe, while public derivatives remain the three spatial coordinates.

use std::fmt;

use crate::{Axis, DerivativePoint, Point, RbfKernel, SecondDerivative};

use super::derivatives::{axis_index, radius_and_deltas, second_derivative_axes};

/// A stable replacement for frozen `R`'s integer `-666` derivative sentinel.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum KernelError {
    /// Frozen Surfe implements the linear-radius value but not its derivatives.
    LinearDerivativeUnavailable,
}

impl fmt::Display for KernelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::LinearDerivativeUnavailable => {
                "the frozen linear-radius kernel does not implement derivatives"
            }
        })
    }
}

impl std::error::Error for KernelError {}

/// All values exposed by a differentiable frozen isotropic kernel at two points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KernelEvaluation {
    basis: f64,
    first_at_first: [f64; 3],
    first_at_second: [f64; 3],
    mixed_hessian: [[f64; 3]; 3],
}

impl KernelEvaluation {
    pub(super) const fn from_components(
        basis: f64,
        first_at_first: [f64; 3],
        first_at_second: [f64; 3],
        mixed_hessian: [[f64; 3]; 3],
    ) -> Self {
        Self {
            basis,
            first_at_first,
            first_at_second,
            mixed_hessian,
        }
    }

    pub const fn basis(self) -> f64 {
        self.basis
    }

    pub const fn first_at_first(self) -> [f64; 3] {
        self.first_at_first
    }

    pub const fn first_at_second(self) -> [f64; 3] {
        self.first_at_second
    }

    /// Rows act on the first point; columns act on the second point.
    pub const fn mixed_hessian(self) -> [[f64; 3]; 3] {
        self.mixed_hessian
    }

    pub fn is_finite(self) -> bool {
        self.basis.is_finite()
            && self.first_at_first.into_iter().all(f64::is_finite)
            && self.first_at_second.into_iter().all(f64::is_finite)
            && self.mixed_hessian.into_iter().flatten().all(f64::is_finite)
    }
}

/// One of the nine isotropic kernels selectable in frozen Surfe.
///
/// `parameter` is used as the shape parameter by Gaussian, MQ, MQ3, IMQ, and
/// MaternC4, and as the compact support cutoff by WendlandC2. It is retained
/// but ignored by Cubic, TPS, and Linear, matching the frozen factory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsotropicKernel {
    kind: RbfKernel,
    parameter: f64,
}

impl IsotropicKernel {
    pub const fn new(kind: RbfKernel, parameter: f64) -> Self {
        Self { kind, parameter }
    }

    pub const fn kind(self) -> RbfKernel {
        self.kind
    }

    pub const fn parameter(self) -> f64 {
        self.parameter
    }

    /// Frozen `basis_pt_pt`, including the fourth coordinate in the radius.
    pub fn basis(self, first: &Point, second: &Point) -> f64 {
        let (radius, _) = radius_and_deltas(first, second);
        match self.kind {
            RbfKernel::Cubic => radius * radius * radius,
            RbfKernel::Gaussian => (-(self.parameter * self.parameter * radius * radius)).exp(),
            RbfKernel::Multiquadric => (self.parameter + radius * radius).powf(0.5),
            RbfKernel::MultiquadricCubic => (self.parameter + radius * radius).powf(1.5),
            RbfKernel::ThinPlateSpline => {
                if radius != 0.0 {
                    radius.powf(4.0) * radius.ln()
                } else {
                    0.0
                }
            }
            RbfKernel::InverseMultiquadric => 1.0 / (self.parameter + radius * radius).powf(0.5),
            RbfKernel::Linear => radius,
            RbfKernel::WendlandC2 => {
                if radius > self.parameter {
                    0.0
                } else {
                    let scaled_radius = radius / self.parameter;
                    (1.0 - scaled_radius).powf(4.0) * (1.0 + 4.0 * scaled_radius)
                }
            }
            RbfKernel::MaternC4 => {
                let scaled_radius = self.parameter * radius;
                (-scaled_radius).exp() * (3.0 + 3.0 * scaled_radius + scaled_radius * scaled_radius)
            }
        }
    }

    /// A spatial first derivative with respect to either kernel point.
    pub fn first_derivative(
        self,
        first: &Point,
        second: &Point,
        with_respect_to: DerivativePoint,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        if self.kind == RbfKernel::Linear {
            return Err(KernelError::LinearDerivativeUnavailable);
        }
        let (radius, deltas) = radius_and_deltas(first, second);
        if self.kind == RbfKernel::WendlandC2 && radius > self.parameter {
            return Ok(0.0);
        }
        let delta = deltas[axis_index(axis)];
        if self.kind == RbfKernel::ThinPlateSpline {
            if radius == 0.0 {
                return Ok(0.0);
            }
            return Ok(match with_respect_to {
                DerivativePoint::First => {
                    delta * radius * radius + 4.0 * delta * radius * radius * radius.ln()
                }
                // Keep the two source subtractions instead of negating the
                // point-one expression: they observably produce +0 for a
                // zero component when the logarithmic term is -0.
                DerivativePoint::Second => {
                    -delta * radius * radius - 4.0 * delta * radius * radius * radius.ln()
                }
            });
        }
        let derivative_at_first = match self.kind {
            RbfKernel::Cubic => 3.0 * radius * delta,
            RbfKernel::Gaussian => {
                -2.0 * self.parameter
                    * self.parameter
                    * delta
                    * (-(self.parameter * self.parameter * radius * radius)).exp()
            }
            RbfKernel::Multiquadric => delta / (self.parameter + radius * radius).powf(0.5),
            RbfKernel::MultiquadricCubic => {
                3.0 * delta * (self.parameter + radius * radius).powf(0.5)
            }
            RbfKernel::ThinPlateSpline => unreachable!("TPS derivatives return above"),
            RbfKernel::InverseMultiquadric => -delta / (self.parameter + radius * radius).powf(1.5),
            RbfKernel::WendlandC2 => {
                20.0 * delta * (radius - self.parameter).powf(3.0) / self.parameter.powf(5.0)
            }
            RbfKernel::MaternC4 => {
                let scaled_radius = self.parameter * radius;
                -(-scaled_radius).exp()
                    * self.parameter
                    * self.parameter
                    * (1.0 + scaled_radius)
                    * delta
            }
            RbfKernel::Linear => unreachable!("linear derivatives return above"),
        };

        Ok(match with_respect_to {
            DerivativePoint::First => derivative_at_first,
            DerivativePoint::Second => -derivative_at_first,
        })
    }

    /// A mixed spatial derivative, first by point 1 and then by point 2.
    pub fn mixed_second_derivative(
        self,
        first: &Point,
        second: &Point,
        first_axis: Axis,
        second_axis: Axis,
    ) -> Result<f64, KernelError> {
        if self.kind == RbfKernel::Linear {
            return Err(KernelError::LinearDerivativeUnavailable);
        }
        if axis_index(first_axis) > axis_index(second_axis) {
            // Frozen `dyx`, `dzx`, and `dzy` call `dxy`, `dxz`, and `dyz`
            // rather than recomputing the symmetric expression with reversed
            // multiplication order. That order is observable at one ULP in
            // complete Vector Field Hessian matrices.
            return self.mixed_second_derivative(first, second, second_axis, first_axis);
        }

        let (radius, deltas) = radius_and_deltas(first, second);
        let first_index = axis_index(first_axis);
        let second_index = axis_index(second_axis);
        let first_delta = deltas[first_index];
        let second_delta = deltas[second_index];
        let diagonal = first_index == second_index;

        let value = match self.kind {
            RbfKernel::Cubic => {
                if radius == 0.0 {
                    0.0
                } else if diagonal {
                    -3.0 * (((first_delta * first_delta) / radius) + radius)
                } else {
                    -3.0 * ((first_delta * second_delta) / radius)
                }
            }
            RbfKernel::Gaussian => {
                let exponential = (-(self.parameter * self.parameter * radius * radius)).exp();
                if diagonal {
                    (2.0 * self.parameter * self.parameter
                        - 4.0 * self.parameter.powf(4.0) * first_delta * first_delta)
                        * exponential
                } else {
                    (-4.0 * self.parameter.powf(4.0) * first_delta * second_delta) * exponential
                }
            }
            RbfKernel::Multiquadric => {
                let shifted_radius = self.parameter + radius * radius;
                if diagonal {
                    (first_delta * first_delta) / shifted_radius.powf(1.5)
                        - (1.0 / shifted_radius.powf(0.5))
                } else {
                    (first_delta * second_delta) / shifted_radius.powf(1.5)
                }
            }
            RbfKernel::MultiquadricCubic => {
                let shifted_radius = self.parameter + radius * radius;
                if diagonal {
                    (-3.0 * first_delta * first_delta * shifted_radius.powf(-0.5))
                        - (3.0 * shifted_radius.powf(0.5))
                } else {
                    -3.0 * first_delta * second_delta * shifted_radius.powf(-0.5)
                }
            }
            RbfKernel::ThinPlateSpline => {
                if radius == 0.0 {
                    0.0
                } else if diagonal {
                    self.thin_plate_spline_diagonal(axis_index(first_axis), radius, deltas)
                } else {
                    -6.0 * first_delta * second_delta
                        - 8.0 * first_delta * second_delta * radius.ln()
                }
            }
            RbfKernel::InverseMultiquadric => {
                let shifted_radius = self.parameter + radius * radius;
                if diagonal {
                    (-3.0 * first_delta * first_delta / shifted_radius.powf(2.5))
                        + 1.0 / shifted_radius.powf(1.5)
                } else {
                    -3.0 * first_delta * second_delta / shifted_radius.powf(2.5)
                }
            }
            RbfKernel::WendlandC2 => {
                if radius > self.parameter {
                    0.0
                } else if radius == 0.0 {
                    if diagonal {
                        20.0 / (self.parameter * self.parameter)
                    } else {
                        0.0
                    }
                } else if diagonal {
                    self.wendland_diagonal(axis_index(first_axis), radius, deltas)
                } else {
                    let a =
                        -60.0 * first_delta * second_delta / (self.parameter.powf(5.0) * radius);
                    let b = (self.parameter - radius).powf(2.0);
                    a * b
                }
            }
            RbfKernel::MaternC4 => {
                let scaled_radius = self.parameter * radius;
                if diagonal {
                    (-scaled_radius).exp()
                        * (self.parameter * self.parameter
                            + radius * self.parameter * self.parameter * self.parameter
                            - self.parameter
                                * self.parameter
                                * self.parameter
                                * self.parameter
                                * first_delta
                                * first_delta)
                } else {
                    -(-scaled_radius).exp()
                        * self.parameter
                        * self.parameter
                        * self.parameter
                        * self.parameter
                        * first_delta
                        * second_delta
                }
            }
            RbfKernel::Linear => unreachable!("linear derivatives return above"),
        };
        Ok(value)
    }

    /// Frozen `basis_planar_planar` enum dispatch.
    pub fn mixed_second_derivative_component(
        self,
        first: &Point,
        second: &Point,
        component: SecondDerivative,
    ) -> Result<f64, KernelError> {
        let (first_axis, second_axis) = second_derivative_axes(component);
        self.mixed_second_derivative(first, second, first_axis, second_axis)
    }

    pub fn evaluate(self, first: &Point, second: &Point) -> Result<KernelEvaluation, KernelError> {
        let mut first_at_first = [0.0; 3];
        let mut first_at_second = [0.0; 3];
        let mut mixed_hessian = [[0.0; 3]; 3];
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let index = axis_index(axis);
            first_at_first[index] =
                self.first_derivative(first, second, DerivativePoint::First, axis)?;
            first_at_second[index] =
                self.first_derivative(first, second, DerivativePoint::Second, axis)?;
            for second_axis in [Axis::X, Axis::Y, Axis::Z] {
                mixed_hessian[index][axis_index(second_axis)] =
                    self.mixed_second_derivative(first, second, axis, second_axis)?;
            }
        }
        Ok(KernelEvaluation::from_components(
            self.basis(first, second),
            first_at_first,
            first_at_second,
            mixed_hessian,
        ))
    }

    fn thin_plate_spline_diagonal(self, axis: usize, radius: f64, deltas: [f64; 3]) -> f64 {
        let logarithm = radius.ln();
        match axis {
            0 => {
                -7.0 * deltas[0] * deltas[0]
                    - deltas[1] * deltas[1]
                    - deltas[2] * deltas[2]
                    - 8.0 * deltas[0] * deltas[0] * logarithm
                    - 4.0 * radius * radius * logarithm
            }
            1 => {
                -7.0 * deltas[1] * deltas[1]
                    - deltas[0] * deltas[0]
                    - deltas[2] * deltas[2]
                    - 8.0 * deltas[1] * deltas[1] * logarithm
                    - 4.0 * radius * radius * logarithm
            }
            2 => {
                -7.0 * deltas[2] * deltas[2]
                    - deltas[1] * deltas[1]
                    - deltas[0] * deltas[0]
                    - 8.0 * deltas[2] * deltas[2] * logarithm
                    - 4.0 * radius * radius * logarithm
            }
            _ => unreachable!("axis index is always 0..3"),
        }
    }

    fn wendland_diagonal(self, axis: usize, radius: f64, deltas: [f64; 3]) -> f64 {
        let a = -20.0 / (self.parameter.powf(5.0) * radius * radius);
        let b = (self.parameter - radius).powf(2.0);
        let c = match axis {
            0 => {
                -self.parameter * radius * radius
                    + radius
                        * (4.0 * deltas[0] * deltas[0]
                            + deltas[1] * deltas[1]
                            + deltas[2] * deltas[2])
            }
            1 => {
                -self.parameter * radius * radius
                    + radius
                        * (deltas[0] * deltas[0]
                            + 4.0 * deltas[1] * deltas[1]
                            + deltas[2] * deltas[2])
            }
            2 => {
                -self.parameter * radius * radius
                    + radius
                        * (deltas[0] * deltas[0]
                            + deltas[1] * deltas[1]
                            + 4.0 * deltas[2] * deltas[2])
            }
            _ => unreachable!("axis index is always 0..3"),
        };
        a * b * c
    }
}
