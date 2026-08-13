//! Frozen Surfe's Lagrangian-polynomial modified radial kernel.
//!
//! Source: `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`,
//! `Modified_Kernel::{Modified_Kernel,basis_pt_pt,basis_pt_planar_x,
//! basis_planar_x_pt,basis_pt_planar_y,basis_planar_y_pt,
//! basis_pt_planar_z,basis_planar_z_pt,basis_pt_tangent,basis_tangent_pt,
//! basis_planar_planar,basis_tangent_tangent,basis_planar_tangent,
//! basis_tangent_planar}`.

use crate::{
    Axis, DerivativePoint, Error, FirstDerivative, Interface, LagrangianPolynomialBasis, Point,
    SecondDerivative, Tangent,
};

use super::{
    derivatives::{axis_index, second_derivative_axes},
    AnisotropicKernel, IsotropicKernel, KernelError, KernelEvaluation,
};

#[derive(Clone, Copy, Debug, PartialEq)]
enum RadialKernel {
    Isotropic(IsotropicKernel),
    Anisotropic(AnisotropicKernel),
}

impl RadialKernel {
    fn basis(self, first: &Point, second: &Point) -> f64 {
        match self {
            Self::Isotropic(kernel) => kernel.basis(first, second),
            Self::Anisotropic(kernel) => kernel.basis(first, second),
        }
    }

    fn first_derivative(
        self,
        first: &Point,
        second: &Point,
        with_respect_to: DerivativePoint,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        match self {
            Self::Isotropic(kernel) => {
                kernel.first_derivative(first, second, with_respect_to, axis)
            }
            Self::Anisotropic(kernel) => {
                kernel.first_derivative(first, second, with_respect_to, axis)
            }
        }
    }

    fn mixed_second_derivative(
        self,
        first: &Point,
        second: &Point,
        first_axis: Axis,
        second_axis: Axis,
    ) -> Result<f64, KernelError> {
        match self {
            Self::Isotropic(kernel) => {
                kernel.mixed_second_derivative(first, second, first_axis, second_axis)
            }
            Self::Anisotropic(kernel) => {
                kernel.mixed_second_derivative(first, second, first_axis, second_axis)
            }
        }
    }
}

/// A pure-Rust equivalent of frozen `Modified_Kernel`.
///
/// Construction owns both the underlying radial kernel and the T11
/// Lagrangian basis. Invalid interface groups are reported at the same outer
/// category used by Surfe's `setup_basis_functions` call chain.
#[derive(Clone, Debug)]
pub struct ModifiedKernel {
    radial: RadialKernel,
    lagrangian: LagrangianPolynomialBasis,
}

impl ModifiedKernel {
    /// Construct a modified kernel over any frozen isotropic radial kernel.
    pub fn from_isotropic(
        radial: IsotropicKernel,
        interface_point_lists: &[Vec<Interface>],
    ) -> Result<Self, Error> {
        Self::new(RadialKernel::Isotropic(radial), interface_point_lists)
    }

    /// Construct a modified kernel over a previously validated anisotropic
    /// radial kernel.
    pub fn from_anisotropic(
        radial: AnisotropicKernel,
        interface_point_lists: &[Vec<Interface>],
    ) -> Result<Self, Error> {
        Self::new(RadialKernel::Anisotropic(radial), interface_point_lists)
    }

    fn new(radial: RadialKernel, interface_point_lists: &[Vec<Interface>]) -> Result<Self, Error> {
        let lagrangian = LagrangianPolynomialBasis::new(interface_point_lists)
            .map_err(|_| Error::ModifiedKernelCreationFailure)?;
        Ok(Self { radial, lagrangian })
    }

    /// The selected points and coefficients used by the modification.
    pub const fn lagrangian_basis(&self) -> &LagrangianPolynomialBasis {
        &self.lagrangian
    }

    /// Frozen `basis_pt_pt`.
    pub fn basis_pt_pt(&self, first: &Point, second: &Point) -> f64 {
        let first_polynomial = self.lagrangian.values(first);
        let second_polynomial = self.lagrangian.values(second);
        self.modified_value(first, second, first_polynomial, second_polynomial)
    }

    /// A spatial first derivative of the modified value at either point.
    pub fn first_derivative(
        &self,
        first: &Point,
        second: &Point,
        with_respect_to: DerivativePoint,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        match with_respect_to {
            DerivativePoint::First => self.modified_first_at_first(first, second, axis),
            DerivativePoint::Second => self.modified_first_at_second(first, second, axis),
        }
    }

    /// Frozen `basis_pt_planar_{x,y,z}`.
    pub fn basis_pt_planar(
        &self,
        first: &Point,
        second: &Point,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        self.modified_first_at_second(first, second, axis)
    }

    /// Frozen `basis_planar_{x,y,z}_pt`.
    pub fn basis_planar_pt(
        &self,
        first: &Point,
        second: &Point,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        self.modified_first_at_first(first, second, axis)
    }

    /// A row-by-column mixed derivative of the modified value.
    pub fn mixed_second_derivative(
        &self,
        first: &Point,
        second: &Point,
        first_axis: Axis,
        second_axis: Axis,
    ) -> Result<f64, KernelError> {
        let first_polynomial = self.polynomial_derivative(first, first_axis);
        let second_polynomial = self.polynomial_derivative(second, second_axis);
        let unisolvent = self.lagrangian.unisolvent_points();
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        let mut t3 = 0.0;
        let mut t4 = 0.0;

        for j in 0..4 {
            let b1 = self.radial.first_derivative(
                &unisolvent[j],
                second,
                DerivativePoint::Second,
                second_axis,
            )?;
            let b2 = self.radial.first_derivative(
                first,
                &unisolvent[j],
                DerivativePoint::First,
                first_axis,
            )?;
            t1 += first_polynomial[j] * b1;
            t2 += second_polynomial[j] * b2;
            t3 += first_polynomial[j] * second_polynomial[j];
            for k in 0..4 {
                if k != j {
                    let b3 = self.radial.basis(&unisolvent[j], &unisolvent[k]);
                    t4 += second_polynomial[k] * b3 * first_polynomial[j];
                }
            }
        }

        let base = self
            .radial
            .mixed_second_derivative(first, second, first_axis, second_axis)?;
        Ok(base - t1 - t2 + t3 + t4)
    }

    /// Frozen `basis_planar_planar` enum dispatch.
    pub fn basis_planar_planar(
        &self,
        first: &Point,
        second: &Point,
        component: SecondDerivative,
    ) -> Result<f64, KernelError> {
        let (first_axis, second_axis) = second_derivative_axes(component);
        self.mixed_second_derivative(first, second, first_axis, second_axis)
    }

    /// Frozen `basis_pt_tangent`.
    pub fn basis_pt_tangent(&self, first: &Point, second: &Tangent) -> Result<f64, KernelError> {
        let dx = self.modified_first_at_second(first, second.point(), Axis::X)?;
        let dy = self.modified_first_at_second(first, second.point(), Axis::Y)?;
        let dz = self.modified_first_at_second(first, second.point(), Axis::Z)?;
        Ok(dx * second.tx() + dy * second.ty() + dz * second.tz())
    }

    /// Frozen `basis_tangent_pt`.
    pub fn basis_tangent_pt(&self, first: &Tangent, second: &Point) -> Result<f64, KernelError> {
        let dx = self.modified_first_at_first(first.point(), second, Axis::X)?;
        let dy = self.modified_first_at_first(first.point(), second, Axis::Y)?;
        let dz = self.modified_first_at_first(first.point(), second, Axis::Z)?;
        Ok(dx * first.tx() + dy * first.ty() + dz * first.tz())
    }

    /// Frozen `basis_tangent_tangent`.
    pub fn basis_tangent_tangent(
        &self,
        first: &Tangent,
        second: &Tangent,
    ) -> Result<f64, KernelError> {
        let hessian = self.modified_hessian(first.point(), second.point())?;
        Ok(first.tx() * second.tx() * hessian[0][0]
            + first.tx() * second.ty() * hessian[0][1]
            + first.tx() * second.tz() * hessian[0][2]
            + first.ty() * second.tx() * hessian[1][0]
            + first.ty() * second.ty() * hessian[1][1]
            + first.ty() * second.tz() * hessian[1][2]
            + first.tz() * second.tx() * hessian[2][0]
            + first.tz() * second.ty() * hessian[2][1]
            + first.tz() * second.tz() * hessian[2][2])
    }

    /// Frozen `basis_planar_tangent`.
    pub fn basis_planar_tangent(
        &self,
        first: &Point,
        second: &Tangent,
        component: FirstDerivative,
    ) -> Result<f64, KernelError> {
        let row = axis_index(first_derivative_axis(component));
        let hessian = self.modified_hessian(first, second.point())?;
        Ok(second.tx() * hessian[row][0]
            + second.ty() * hessian[row][1]
            + second.tz() * hessian[row][2])
    }

    /// Frozen `basis_tangent_planar`.
    pub fn basis_tangent_planar(
        &self,
        first: &Tangent,
        second: &Point,
        component: FirstDerivative,
    ) -> Result<f64, KernelError> {
        let column = axis_index(first_derivative_axis(component));
        let hessian = self.modified_hessian(first.point(), second)?;
        Ok(first.tx() * hessian[0][column]
            + first.ty() * hessian[1][column]
            + first.tz() * hessian[2][column])
    }

    /// Value, both spatial gradients, and the row-by-column mixed Hessian.
    pub fn evaluate(&self, first: &Point, second: &Point) -> Result<KernelEvaluation, KernelError> {
        let mut first_at_first = [0.0; 3];
        let mut first_at_second = [0.0; 3];
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let index = axis_index(axis);
            first_at_first[index] = self.modified_first_at_first(first, second, axis)?;
            first_at_second[index] = self.modified_first_at_second(first, second, axis)?;
        }
        Ok(KernelEvaluation::from_components(
            self.basis_pt_pt(first, second),
            first_at_first,
            first_at_second,
            self.modified_hessian(first, second)?,
        ))
    }

    fn modified_value(
        &self,
        first: &Point,
        second: &Point,
        first_polynomial: [f64; 4],
        second_polynomial: [f64; 4],
    ) -> f64 {
        let unisolvent = self.lagrangian.unisolvent_points();
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        let mut t3 = 0.0;
        let mut t4 = 0.0;

        for j in 0..4 {
            let b1 = self.radial.basis(&unisolvent[j], second);
            let b2 = self.radial.basis(first, &unisolvent[j]);
            t1 += first_polynomial[j] * b1;
            t2 += second_polynomial[j] * b2;
            t3 += first_polynomial[j] * second_polynomial[j];
            for k in 0..4 {
                if k != j {
                    let b3 = self.radial.basis(&unisolvent[j], &unisolvent[k]);
                    t4 += first_polynomial[j] * second_polynomial[k] * b3;
                }
            }
        }

        self.radial.basis(first, second) - t1 - t2 + t3 + t4
    }

    fn modified_first_at_first(
        &self,
        first: &Point,
        second: &Point,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        let first_polynomial = self.polynomial_derivative(first, axis);
        let second_polynomial = self.lagrangian.values(second);
        let unisolvent = self.lagrangian.unisolvent_points();
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        let mut t3 = 0.0;
        let mut t4 = 0.0;

        for j in 0..4 {
            let b1 = self.radial.basis(&unisolvent[j], second);
            let b2 = self.radial.first_derivative(
                first,
                &unisolvent[j],
                DerivativePoint::First,
                axis,
            )?;
            t1 += first_polynomial[j] * b1;
            t2 += second_polynomial[j] * b2;
            t3 += first_polynomial[j] * second_polynomial[j];
            for k in 0..4 {
                if k != j {
                    let b3 = self.radial.basis(&unisolvent[j], &unisolvent[k]);
                    t4 += first_polynomial[j] * second_polynomial[k] * b3;
                }
            }
        }

        let base = self
            .radial
            .first_derivative(first, second, DerivativePoint::First, axis)?;
        Ok(base - t1 - t2 + t3 + t4)
    }

    fn modified_first_at_second(
        &self,
        first: &Point,
        second: &Point,
        axis: Axis,
    ) -> Result<f64, KernelError> {
        let first_polynomial = self.lagrangian.values(first);
        let second_polynomial = self.polynomial_derivative(second, axis);
        let unisolvent = self.lagrangian.unisolvent_points();
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        let mut t3 = 0.0;
        let mut t4 = 0.0;

        for j in 0..4 {
            let b1 = self.radial.first_derivative(
                &unisolvent[j],
                second,
                DerivativePoint::Second,
                axis,
            )?;
            let b2 = self.radial.basis(first, &unisolvent[j]);
            t1 += first_polynomial[j] * b1;
            t2 += second_polynomial[j] * b2;
            t3 += first_polynomial[j] * second_polynomial[j];
            for k in 0..4 {
                if k != j {
                    let b3 = self.radial.basis(&unisolvent[j], &unisolvent[k]);
                    t4 += first_polynomial[j] * second_polynomial[k] * b3;
                }
            }
        }

        let base = self
            .radial
            .first_derivative(first, second, DerivativePoint::Second, axis)?;
        Ok(base - t1 - t2 + t3 + t4)
    }

    fn modified_hessian(
        &self,
        first: &Point,
        second: &Point,
    ) -> Result<[[f64; 3]; 3], KernelError> {
        let mut hessian = [[0.0; 3]; 3];
        for first_axis in [Axis::X, Axis::Y, Axis::Z] {
            for second_axis in [Axis::X, Axis::Y, Axis::Z] {
                hessian[axis_index(first_axis)][axis_index(second_axis)] =
                    self.mixed_second_derivative(first, second, first_axis, second_axis)?;
            }
        }
        Ok(hessian)
    }

    fn polynomial_derivative(&self, point: &Point, axis: Axis) -> [f64; 4] {
        match axis {
            Axis::X => self.lagrangian.dx(point),
            Axis::Y => self.lagrangian.dy(point),
            Axis::Z => self.lagrangian.dz(point),
        }
    }
}

const fn first_derivative_axis(component: FirstDerivative) -> Axis {
    match component {
        FirstDerivative::Dx => Axis::X,
        FirstDerivative::Dy => Axis::Y,
        FirstDerivative::Dz => Axis::Z,
    }
}
