//! Frozen Surfe global-anisotropy support and anisotropic radial kernels.
//!
//! Source: `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`,
//! especially `RBFKernel::{get_global_anisotropy,scaled_radius}` and the
//! `ACubic/AGaussian/AMQ/ATPS/AIMQ/AR` implementations. The covariance,
//! eigensystem, eigenvalue floor, and support matrix deliberately retain
//! binary32 intermediate arithmetic.

use std::fmt;

use crate::{Axis, DerivativePoint, Planar, Point, RbfKernel, SecondDerivative};

use super::{
    derivatives::{axis_index, second_derivative_axes},
    KernelError, KernelEvaluation,
};

const EIGENVALUE_FLOOR: f32 = 0.0001;

/// Failures while constructing frozen Surfe's global anisotropy support.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum AnisotropyError {
    /// `get_global_anisotropy` requires at least two planar constraints.
    InsufficientPlanars,
    /// Frozen Surfe has no anisotropic class for this otherwise valid kernel.
    UnsupportedKernel(RbfKernel),
    /// Binary32 covariance or eigensystem arithmetic did not stay finite.
    NonFiniteComputation,
    /// The symmetric binary32 eigensolver did not converge.
    EigenSolverFailure,
}

impl fmt::Display for AnisotropyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InsufficientPlanars => {
                formatter.write_str("global anisotropy requires at least two planar constraints")
            }
            Self::UnsupportedKernel(kernel) => {
                write!(formatter, "frozen Surfe has no anisotropic {kernel} kernel")
            }
            Self::NonFiniteComputation => {
                formatter.write_str("global anisotropy produced a non-finite binary32 value")
            }
            Self::EigenSolverFailure => {
                formatter.write_str("global anisotropy eigensolver did not converge")
            }
        }
    }
}

impl std::error::Error for AnisotropyError {}

/// One of the six anisotropic kernel classes available in frozen Surfe.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnisotropicKernel {
    kind: RbfKernel,
    parameter: f64,
    eigenvalues: [f32; 3],
    clamped_eigenvalues: [f32; 3],
    transform: [[f32; 3]; 3],
    global_plunge: [f32; 3],
}

impl AnisotropicKernel {
    /// Construct the exact anisotropic branch selected by the frozen factory.
    pub fn new(
        kind: RbfKernel,
        parameter: f64,
        planars: &[Planar],
    ) -> Result<Self, AnisotropyError> {
        if !matches!(
            kind,
            RbfKernel::Cubic
                | RbfKernel::Gaussian
                | RbfKernel::Multiquadric
                | RbfKernel::ThinPlateSpline
                | RbfKernel::InverseMultiquadric
                | RbfKernel::Linear
        ) {
            return Err(AnisotropyError::UnsupportedKernel(kind));
        }
        if planars.len() < 2 {
            return Err(AnisotropyError::InsufficientPlanars);
        }

        let covariance = covariance_matrix(planars)?;
        let (eigenvalues, eigenvectors) = self_adjoint_eigen(covariance)?;
        let global_plunge = [eigenvectors[0][0], eigenvectors[1][0], eigenvectors[2][0]];
        let mut clamped_eigenvalues = eigenvalues;
        if clamped_eigenvalues[0] < EIGENVALUE_FLOOR {
            clamped_eigenvalues[0] = EIGENVALUE_FLOOR;
        }
        if clamped_eigenvalues[1] < EIGENVALUE_FLOOR {
            clamped_eigenvalues[1] = EIGENVALUE_FLOOR;
        }

        let scale = [
            1.0_f32,
            (clamped_eigenvalues[1] / clamped_eigenvalues[0]).sqrt(),
            (clamped_eigenvalues[2] / clamped_eigenvalues[0]).sqrt(),
        ];
        let mut scaled_eigenvectors = [[0.0_f32; 3]; 3];
        for (row, scaled_row) in scaled_eigenvectors.iter_mut().enumerate() {
            for (column, output) in scaled_row.iter_mut().enumerate() {
                *output = eigenvectors[row][column] * scale[column];
            }
        }
        let mut transform = [[0.0_f32; 3]; 3];
        for (row, transform_row) in transform.iter_mut().enumerate() {
            for (column, output) in transform_row.iter_mut().enumerate() {
                *output = scaled_eigenvectors[row][0] * eigenvectors[column][0]
                    + scaled_eigenvectors[row][1] * eigenvectors[column][1]
                    + scaled_eigenvectors[row][2] * eigenvectors[column][2];
            }
        }

        if !parameter.is_finite()
            || !eigenvalues.into_iter().all(f32::is_finite)
            || !clamped_eigenvalues.into_iter().all(f32::is_finite)
            || !global_plunge.into_iter().all(f32::is_finite)
            || !transform.into_iter().flatten().all(f32::is_finite)
        {
            return Err(AnisotropyError::NonFiniteComputation);
        }

        Ok(Self {
            kind,
            parameter,
            eigenvalues,
            clamped_eigenvalues,
            transform,
            global_plunge,
        })
    }

    pub const fn kind(self) -> RbfKernel {
        self.kind
    }

    pub const fn parameter(self) -> f64 {
        self.parameter
    }

    /// Ascending binary32 eigenvalues before Surfe's floor is applied.
    pub const fn eigenvalues(self) -> [f32; 3] {
        self.eigenvalues
    }

    /// Eigenvalues after flooring only entries zero and one at `1e-4`.
    pub const fn clamped_eigenvalues(self) -> [f32; 3] {
        self.clamped_eigenvalues
    }

    /// Surfe's binary32 symmetric support matrix, in logical row-major form.
    pub const fn transform(self) -> [[f32; 3]; 3] {
        self.transform
    }

    /// First eigenvector written by frozen Surfe; it is otherwise unused there.
    pub const fn global_plunge(self) -> [f32; 3] {
        self.global_plunge
    }

    /// Frozen `scaled_radius`, which intentionally ignores the point `c` fields.
    pub fn scaled_radius(self, first: &Point, second: &Point) -> f64 {
        self.transformed_deltas(first, second).0
    }

    pub fn basis(self, first: &Point, second: &Point) -> f64 {
        let (radius, _) = self.transformed_deltas(first, second);
        match self.kind {
            RbfKernel::Cubic => radius * radius * radius,
            RbfKernel::Gaussian => (-(self.parameter * self.parameter * radius * radius)).exp(),
            RbfKernel::Multiquadric => (self.parameter + radius * radius).powf(0.5),
            RbfKernel::ThinPlateSpline => {
                if radius != 0.0 {
                    radius.powf(4.0) * radius.ln()
                } else {
                    0.0
                }
            }
            RbfKernel::InverseMultiquadric => 1.0 / (self.parameter + radius * radius).powf(0.5),
            RbfKernel::Linear => radius,
            _ => unreachable!("constructor rejects unsupported anisotropic kernels"),
        }
    }

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
        let (radius, transformed) = self.transformed_deltas(first, second);
        let projected = self.projected_deltas(transformed);
        let q = projected[axis_index(axis)];
        let derivative_at_first = match self.kind {
            RbfKernel::Cubic => 3.0 * radius * q,
            RbfKernel::Gaussian => {
                -2.0 * self.parameter
                    * self.parameter
                    * q
                    * (-(self.parameter * self.parameter * radius * radius)).exp()
            }
            RbfKernel::Multiquadric => q / (self.parameter + radius * radius).powf(0.5),
            RbfKernel::ThinPlateSpline => {
                if radius != 0.0 {
                    let a = q * radius * radius;
                    a + 4.0 * a * radius.ln()
                } else {
                    0.0
                }
            }
            RbfKernel::InverseMultiquadric => -q / (self.parameter + radius * radius).powf(1.5),
            RbfKernel::Linear => unreachable!("linear derivatives return above"),
            _ => unreachable!("constructor rejects unsupported anisotropic kernels"),
        };
        Ok(match with_respect_to {
            DerivativePoint::First => derivative_at_first,
            DerivativePoint::Second => -derivative_at_first,
        })
    }

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
            // Frozen `dyx`, `dzx`, and `dzy` are aliases to `dxy`, `dxz`,
            // and `dyz`; preserve their exact recomputation path.
            return self.mixed_second_derivative(first, second, second_axis, first_axis);
        }
        let (radius, transformed) = self.transformed_deltas(first, second);
        let projected = self.projected_deltas(transformed);
        let first_index = axis_index(first_axis);
        let second_index = axis_index(second_axis);
        let a = projected[first_index];
        let b = projected[second_index];
        let metric = self.metric(first_index, second_index);

        let value = match self.kind {
            RbfKernel::Cubic => {
                if radius == 0.0 {
                    0.0
                } else {
                    -3.0 * (((a * b) / radius) + metric * radius)
                }
            }
            RbfKernel::Gaussian => {
                let c = self.parameter
                    * self.parameter
                    * (-(self.parameter * self.parameter * radius * radius)).exp();
                2.0 * c * (metric - 2.0 * self.parameter * self.parameter * a * b)
            }
            RbfKernel::Multiquadric => {
                let shifted = self.parameter + radius * radius;
                (a * b) / shifted.powf(1.5) - metric / shifted.powf(0.5)
            }
            RbfKernel::ThinPlateSpline => {
                if radius == 0.0 {
                    0.0
                } else {
                    -6.0 * a * b
                        - metric * radius * radius
                        - 8.0 * a * b * radius.ln()
                        - 4.0 * metric * radius * radius * radius.ln()
                }
            }
            RbfKernel::InverseMultiquadric => {
                let shifted = self.parameter + radius * radius;
                -3.0 * a * b / shifted.powf(2.5) + metric / shifted.powf(1.5)
            }
            RbfKernel::Linear => unreachable!("linear derivatives return above"),
            _ => unreachable!("constructor rejects unsupported anisotropic kernels"),
        };
        Ok(value)
    }

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

    fn transformed_deltas(self, first: &Point, second: &Point) -> (f64, [f64; 3]) {
        let dx = first.x() - second.x();
        let dy = first.y() - second.y();
        let dz = first.z() - second.z();
        let transformed = [
            self.transform[0][0] as f64 * dx
                + self.transform[0][1] as f64 * dy
                + self.transform[0][2] as f64 * dz,
            self.transform[1][0] as f64 * dx
                + self.transform[1][1] as f64 * dy
                + self.transform[1][2] as f64 * dz,
            self.transform[2][0] as f64 * dx
                + self.transform[2][1] as f64 * dy
                + self.transform[2][2] as f64 * dz,
        ];
        let radius = (transformed[0] * transformed[0]
            + transformed[1] * transformed[1]
            + transformed[2] * transformed[2])
            .sqrt();
        (radius, transformed)
    }

    fn projected_deltas(self, transformed: [f64; 3]) -> [f64; 3] {
        [
            self.transform[0][0] as f64 * transformed[0]
                + self.transform[1][0] as f64 * transformed[1]
                + self.transform[2][0] as f64 * transformed[2],
            self.transform[0][1] as f64 * transformed[0]
                + self.transform[1][1] as f64 * transformed[1]
                + self.transform[2][1] as f64 * transformed[2],
            self.transform[0][2] as f64 * transformed[0]
                + self.transform[1][2] as f64 * transformed[1]
                + self.transform[2][2] as f64 * transformed[2],
        ]
    }

    fn metric(self, first_axis: usize, second_axis: usize) -> f64 {
        // Both operands are `Matrix3f` coefficients in the frozen source, so
        // every product and the sum round in binary32 before assignment to a
        // C++ `double` local.
        (self.transform[0][first_axis] * self.transform[0][second_axis]
            + self.transform[1][first_axis] * self.transform[1][second_axis]
            + self.transform[2][first_axis] * self.transform[2][second_axis]) as f64
    }
}

fn covariance_matrix(planars: &[Planar]) -> Result<[[f32; 3]; 3], AnisotropyError> {
    let mut sum_xx = 0.0_f64;
    let mut sum_xy = 0.0_f64;
    let mut sum_xz = 0.0_f64;
    let mut sum_yy = 0.0_f64;
    let mut sum_yz = 0.0_f64;
    let mut sum_zz = 0.0_f64;
    for planar in planars {
        let [nx, ny, nz] = planar.normal();
        sum_xx += nx * nx;
        sum_yy += ny * ny;
        sum_zz += nz * nz;
        sum_xy += nx * ny;
        sum_xz += nx * nz;
        sum_yz += ny * nz;
    }
    let matrix = [
        [sum_xx as f32, sum_xy as f32, sum_xz as f32],
        [sum_xy as f32, sum_yy as f32, sum_yz as f32],
        [sum_xz as f32, sum_yz as f32, sum_zz as f32],
    ];
    if matrix.into_iter().flatten().all(f32::is_finite) {
        Ok(matrix)
    } else {
        Err(AnisotropyError::NonFiniteComputation)
    }
}

/// Pure-Rust equivalent of the frozen Eigen `SelfAdjointEigenSolver<Matrix3f>`
/// path: lower-triangle scaling, the 3x3 tridiagonal specialization, implicit
/// QR iterations, ascending selection sort, and binary32 rescaling.
fn self_adjoint_eigen(matrix: [[f32; 3]; 3]) -> Result<([f32; 3], [[f32; 3]; 3]), AnisotropyError> {
    let mut scale = 0.0_f32;
    for (row, values) in matrix.iter().enumerate() {
        for value in values.iter().take(row + 1) {
            scale = scale.max(value.abs());
        }
    }
    if scale == 0.0 {
        scale = 1.0;
    }
    let mut lower = [[0.0_f32; 3]; 3];
    for (row, lower_row) in lower.iter_mut().enumerate() {
        for (column, output) in lower_row.iter_mut().enumerate().take(row + 1) {
            *output = matrix[row][column] / scale;
        }
    }

    let mut diagonal = [0.0_f32; 3];
    let mut subdiagonal = [0.0_f32; 2];
    let mut eigenvectors;
    diagonal[0] = lower[0][0];
    let v1_norm_squared = lower[2][0] * lower[2][0];
    if v1_norm_squared <= f32::MIN_POSITIVE {
        diagonal[1] = lower[1][1];
        diagonal[2] = lower[2][2];
        subdiagonal[0] = lower[1][0];
        subdiagonal[1] = lower[2][1];
        eigenvectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    } else {
        let beta = (lower[1][0] * lower[1][0] + v1_norm_squared).sqrt();
        let inverse_beta = 1.0 / beta;
        let m01 = lower[1][0] * inverse_beta;
        let m02 = lower[2][0] * inverse_beta;
        let q = 2.0 * m01 * lower[2][1] + m02 * (lower[2][2] - lower[1][1]);
        diagonal[1] = lower[1][1] + m02 * q;
        diagonal[2] = lower[2][2] - m02 * q;
        subdiagonal[0] = beta;
        subdiagonal[1] = lower[2][1] - m01 * q;
        eigenvectors = [[1.0, 0.0, 0.0], [0.0, m01, m02], [0.0, m02, -m01]];
    }

    let precision = 2.0 * f32::EPSILON;
    let mut end = 2_usize;
    let mut start = 0_usize;
    let mut iterations = 0_usize;
    while end > 0 {
        for index in start..end {
            if subdiagonal[index].abs()
                <= (diagonal[index].abs() + diagonal[index + 1].abs()) * precision
                || subdiagonal[index].abs() <= f32::MIN_POSITIVE
            {
                subdiagonal[index] = 0.0;
            }
        }
        while end > 0 && subdiagonal[end - 1] == 0.0 {
            end -= 1;
        }
        if end == 0 {
            break;
        }
        iterations += 1;
        if iterations > 90 {
            return Err(AnisotropyError::EigenSolverFailure);
        }
        start = end - 1;
        while start > 0 && subdiagonal[start - 1] != 0.0 {
            start -= 1;
        }
        tridiagonal_qr_step(
            &mut diagonal,
            &mut subdiagonal,
            start,
            end,
            &mut eigenvectors,
        );
    }

    for index in 0..2 {
        let mut minimum = index;
        for candidate in (index + 1)..3 {
            if diagonal[candidate] < diagonal[minimum] {
                minimum = candidate;
            }
        }
        if minimum != index {
            diagonal.swap(index, minimum);
            for row in &mut eigenvectors {
                row.swap(index, minimum);
            }
        }
    }
    for value in &mut diagonal {
        *value *= scale;
    }
    if diagonal.into_iter().all(f32::is_finite)
        && eigenvectors.into_iter().flatten().all(f32::is_finite)
    {
        Ok((diagonal, eigenvectors))
    } else {
        Err(AnisotropyError::NonFiniteComputation)
    }
}

fn tridiagonal_qr_step(
    diagonal: &mut [f32; 3],
    subdiagonal: &mut [f32; 2],
    start: usize,
    end: usize,
    eigenvectors: &mut [[f32; 3]; 3],
) {
    let td = (diagonal[end - 1] - diagonal[end]) * 0.5;
    let e = subdiagonal[end - 1];
    let mut shift = diagonal[end];
    if td == 0.0 {
        shift -= e.abs();
    } else {
        let e_squared = e * e;
        let hypotenuse = td.hypot(e);
        if e_squared == 0.0 {
            let sign = if td > 0.0 { 1.0 } else { -1.0 };
            shift -= (e / (td + sign)) * (e / hypotenuse);
        } else {
            shift -= e_squared / (td + if td > 0.0 { hypotenuse } else { -hypotenuse });
        }
    }

    let mut x = diagonal[start] - shift;
    let mut z = subdiagonal[start];
    for index in start..end {
        let (cosine, sine) = make_givens(x, z);
        let sdk = sine * diagonal[index] + cosine * subdiagonal[index];
        let dkp1 = sine * subdiagonal[index] + cosine * diagonal[index + 1];
        diagonal[index] = cosine * (cosine * diagonal[index] - sine * subdiagonal[index])
            - sine * (cosine * subdiagonal[index] - sine * diagonal[index + 1]);
        diagonal[index + 1] = sine * sdk + cosine * dkp1;
        subdiagonal[index] = cosine * sdk - sine * dkp1;
        if index > start {
            subdiagonal[index - 1] = cosine * subdiagonal[index - 1] - sine * z;
        }
        x = subdiagonal[index];
        if index < end - 1 {
            z = -sine * subdiagonal[index + 1];
            subdiagonal[index + 1] *= cosine;
        }
        for row in &mut *eigenvectors {
            let first = row[index];
            let second = row[index + 1];
            row[index] = cosine * first - sine * second;
            row[index + 1] = sine * first + cosine * second;
        }
    }
}

fn make_givens(first: f32, second: f32) -> (f32, f32) {
    if second == 0.0 {
        (if first < 0.0 { -1.0 } else { 1.0 }, 0.0)
    } else if first == 0.0 {
        (0.0, if second < 0.0 { 1.0 } else { -1.0 })
    } else if first.abs() > second.abs() {
        let ratio = second / first;
        let mut scale = (1.0 + ratio * ratio).sqrt();
        if first < 0.0 {
            scale = -scale;
        }
        let cosine = 1.0 / scale;
        (cosine, -ratio * cosine)
    } else {
        let ratio = first / second;
        let mut scale = (1.0 + ratio * ratio).sqrt();
        if second < 0.0 {
            scale = -scale;
        }
        let sine = -1.0 / scale;
        (-ratio * sine, sine)
    }
}
