//! Model-independent linear functionals over frozen Surfe kernels.
//!
//! Sources:
//! - `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`Kernel::basis_*` and both `RBFKernel`/`Modified_Kernel` implementations)
//! - `surfe_lib/{single_surface,lajaunie,stratigraphic_surfaces,
//!   continuous_property,vector_field}.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`get_interpolation_matrix` and increment-pair difference blocks)

use crate::{
    AnisotropicKernel, Axis, DerivativePoint, IsotropicKernel, KernelError, ModifiedKernel, Point,
    Tangent,
};

/// A model-independent label for one scalar degree of freedom.
///
/// Planar constraints contribute one `Derivative` label per Cartesian axis;
/// the later layout layer is responsible for assigning matrix indices.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DofLabel {
    Value,
    Derivative(Axis),
    Tangent,
    Difference,
}

/// A frozen Surfe linear functional with all of its geometric support.
///
/// `Difference` means value at `positive` minus value at `negative`, matching
/// Lajaunie and Stratigraphic increment pairs. It is deliberately not reduced
/// to a model-specific row or column index here.
#[derive(Clone, Debug)]
pub enum LinearFunctional {
    Value(Point),
    Derivative { point: Point, axis: Axis },
    Tangent(Tangent),
    Difference { positive: Point, negative: Point },
}

impl LinearFunctional {
    pub const fn value(point: Point) -> Self {
        Self::Value(point)
    }

    pub const fn derivative(point: Point, axis: Axis) -> Self {
        Self::Derivative { point, axis }
    }

    pub const fn tangent(tangent: Tangent) -> Self {
        Self::Tangent(tangent)
    }

    /// Construct `Value(positive) - Value(negative)`.
    pub const fn difference(positive: Point, negative: Point) -> Self {
        Self::Difference { positive, negative }
    }

    pub const fn label(&self) -> DofLabel {
        match self {
            Self::Value(_) => DofLabel::Value,
            Self::Derivative { axis, .. } => DofLabel::Derivative(*axis),
            Self::Tangent(_) => DofLabel::Tangent,
            Self::Difference { .. } => DofLabel::Difference,
        }
    }

    /// Expand into primitive terms without assigning a model or matrix index.
    ///
    /// Application retains the frozen nested subtraction order independently
    /// of this inspection representation.
    pub fn expansion(&self) -> Vec<FunctionalTerm<'_>> {
        match self {
            Self::Value(point) => vec![FunctionalTerm::new(1.0, FunctionalPrimitive::Value(point))],
            Self::Derivative { point, axis } => vec![FunctionalTerm::new(
                1.0,
                FunctionalPrimitive::Derivative { point, axis: *axis },
            )],
            Self::Tangent(tangent) => vec![FunctionalTerm::new(
                1.0,
                FunctionalPrimitive::Tangent(tangent),
            )],
            Self::Difference { positive, negative } => vec![
                FunctionalTerm::new(1.0, FunctionalPrimitive::Value(positive)),
                FunctionalTerm::new(-1.0, FunctionalPrimitive::Value(negative)),
            ],
        }
    }

    const fn primitive(&self) -> Option<FunctionalPrimitive<'_>> {
        match self {
            Self::Value(point) => Some(FunctionalPrimitive::Value(point)),
            Self::Derivative { point, axis } => {
                Some(FunctionalPrimitive::Derivative { point, axis: *axis })
            }
            Self::Tangent(tangent) => Some(FunctionalPrimitive::Tangent(tangent)),
            Self::Difference { .. } => None,
        }
    }
}

/// A borrowed primitive used by [`LinearFunctional::expansion`].
#[derive(Clone, Copy, Debug)]
pub enum FunctionalPrimitive<'a> {
    Value(&'a Point),
    Derivative { point: &'a Point, axis: Axis },
    Tangent(&'a Tangent),
}

impl FunctionalPrimitive<'_> {
    pub const fn label(self) -> DofLabel {
        match self {
            Self::Value(_) => DofLabel::Value,
            Self::Derivative { axis, .. } => DofLabel::Derivative(axis),
            Self::Tangent(_) => DofLabel::Tangent,
        }
    }
}

/// One coefficient and primitive in a functional expansion.
#[derive(Clone, Copy, Debug)]
pub struct FunctionalTerm<'a> {
    coefficient: f64,
    primitive: FunctionalPrimitive<'a>,
}

impl<'a> FunctionalTerm<'a> {
    const fn new(coefficient: f64, primitive: FunctionalPrimitive<'a>) -> Self {
        Self {
            coefficient,
            primitive,
        }
    }

    pub const fn coefficient(self) -> f64 {
        self.coefficient
    }

    pub const fn primitive(self) -> FunctionalPrimitive<'a> {
        self.primitive
    }
}

/// A borrowed, unified call layer for ordinary and Modified kernels.
#[derive(Clone, Copy, Debug)]
pub enum FunctionalKernel<'a> {
    Isotropic(&'a IsotropicKernel),
    Anisotropic(&'a AnisotropicKernel),
    Modified(&'a ModifiedKernel),
}

impl<'a> From<&'a IsotropicKernel> for FunctionalKernel<'a> {
    fn from(kernel: &'a IsotropicKernel) -> Self {
        Self::Isotropic(kernel)
    }
}

impl<'a> From<&'a AnisotropicKernel> for FunctionalKernel<'a> {
    fn from(kernel: &'a AnisotropicKernel) -> Self {
        Self::Anisotropic(kernel)
    }
}

impl<'a> From<&'a ModifiedKernel> for FunctionalKernel<'a> {
    fn from(kernel: &'a ModifiedKernel) -> Self {
        Self::Modified(kernel)
    }
}

impl FunctionalKernel<'_> {
    /// Apply a row functional to the first kernel parameter and a column
    /// functional to the second parameter.
    pub fn apply(
        self,
        first: &LinearFunctional,
        second: &LinearFunctional,
    ) -> Result<f64, KernelError> {
        match first {
            LinearFunctional::Difference { positive, negative } => {
                let positive_value = self
                    .apply_primitive_to_functional(FunctionalPrimitive::Value(positive), second)?;
                let negative_value = self
                    .apply_primitive_to_functional(FunctionalPrimitive::Value(negative), second)?;
                Ok(positive_value - negative_value)
            }
            _ => self.apply_primitive_to_functional(
                first
                    .primitive()
                    .expect("non-difference functional must be primitive"),
                second,
            ),
        }
    }

    fn apply_primitive_to_functional(
        self,
        first: FunctionalPrimitive<'_>,
        second: &LinearFunctional,
    ) -> Result<f64, KernelError> {
        match second {
            LinearFunctional::Difference { positive, negative } => {
                let positive_value =
                    self.apply_primitives(first, FunctionalPrimitive::Value(positive))?;
                let negative_value =
                    self.apply_primitives(first, FunctionalPrimitive::Value(negative))?;
                Ok(positive_value - negative_value)
            }
            _ => self.apply_primitives(
                first,
                second
                    .primitive()
                    .expect("non-difference functional must be primitive"),
            ),
        }
    }

    fn apply_primitives(
        self,
        first: FunctionalPrimitive<'_>,
        second: FunctionalPrimitive<'_>,
    ) -> Result<f64, KernelError> {
        match (first, second) {
            (FunctionalPrimitive::Value(first), FunctionalPrimitive::Value(second)) => {
                Ok(self.basis(first, second))
            }
            (
                FunctionalPrimitive::Value(first),
                FunctionalPrimitive::Derivative {
                    point: second,
                    axis,
                },
            ) => self.first_derivative(first, second, DerivativePoint::Second, axis),
            (
                FunctionalPrimitive::Derivative { point: first, axis },
                FunctionalPrimitive::Value(second),
            ) => self.first_derivative(first, second, DerivativePoint::First, axis),
            (
                FunctionalPrimitive::Derivative {
                    point: first,
                    axis: first_axis,
                },
                FunctionalPrimitive::Derivative {
                    point: second,
                    axis: second_axis,
                },
            ) => self.mixed_second_derivative(first, second, first_axis, second_axis),
            (FunctionalPrimitive::Value(first), FunctionalPrimitive::Tangent(second)) => {
                let direction = second.vector();
                let dx =
                    self.first_derivative(first, second.point(), DerivativePoint::Second, Axis::X)?;
                let dy =
                    self.first_derivative(first, second.point(), DerivativePoint::Second, Axis::Y)?;
                let dz =
                    self.first_derivative(first, second.point(), DerivativePoint::Second, Axis::Z)?;
                Ok(dx * direction[0] + dy * direction[1] + dz * direction[2])
            }
            (FunctionalPrimitive::Tangent(first), FunctionalPrimitive::Value(second)) => {
                let direction = first.vector();
                let dx =
                    self.first_derivative(first.point(), second, DerivativePoint::First, Axis::X)?;
                let dy =
                    self.first_derivative(first.point(), second, DerivativePoint::First, Axis::Y)?;
                let dz =
                    self.first_derivative(first.point(), second, DerivativePoint::First, Axis::Z)?;
                Ok(dx * direction[0] + dy * direction[1] + dz * direction[2])
            }
            (
                FunctionalPrimitive::Derivative { point: first, axis },
                FunctionalPrimitive::Tangent(second),
            ) => {
                let direction = second.vector();
                let dx = self.mixed_second_derivative(first, second.point(), axis, Axis::X)?;
                let dy = self.mixed_second_derivative(first, second.point(), axis, Axis::Y)?;
                let dz = self.mixed_second_derivative(first, second.point(), axis, Axis::Z)?;
                Ok(direction[0] * dx + direction[1] * dy + direction[2] * dz)
            }
            (
                FunctionalPrimitive::Tangent(first),
                FunctionalPrimitive::Derivative {
                    point: second,
                    axis,
                },
            ) => {
                let direction = first.vector();
                let dx = self.mixed_second_derivative(first.point(), second, Axis::X, axis)?;
                let dy = self.mixed_second_derivative(first.point(), second, Axis::Y, axis)?;
                let dz = self.mixed_second_derivative(first.point(), second, Axis::Z, axis)?;
                Ok(direction[0] * dx + direction[1] * dy + direction[2] * dz)
            }
            (FunctionalPrimitive::Tangent(first), FunctionalPrimitive::Tangent(second)) => {
                let first_direction = first.vector();
                let second_direction = second.vector();
                let dxx =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::X, Axis::X)?;
                let dxy =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::X, Axis::Y)?;
                let dxz =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::X, Axis::Z)?;
                let dyx =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::Y, Axis::X)?;
                let dyy =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::Y, Axis::Y)?;
                let dyz =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::Y, Axis::Z)?;
                let dzx =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::Z, Axis::X)?;
                let dzy =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::Z, Axis::Y)?;
                let dzz =
                    self.mixed_second_derivative(first.point(), second.point(), Axis::Z, Axis::Z)?;
                Ok(first_direction[0] * second_direction[0] * dxx
                    + first_direction[0] * second_direction[1] * dxy
                    + first_direction[0] * second_direction[2] * dxz
                    + first_direction[1] * second_direction[0] * dyx
                    + first_direction[1] * second_direction[1] * dyy
                    + first_direction[1] * second_direction[2] * dyz
                    + first_direction[2] * second_direction[0] * dzx
                    + first_direction[2] * second_direction[1] * dzy
                    + first_direction[2] * second_direction[2] * dzz)
            }
        }
    }

    fn basis(self, first: &Point, second: &Point) -> f64 {
        match self {
            Self::Isotropic(kernel) => kernel.basis(first, second),
            Self::Anisotropic(kernel) => kernel.basis(first, second),
            Self::Modified(kernel) => kernel.basis_pt_pt(first, second),
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
            Self::Modified(kernel) => kernel.first_derivative(first, second, with_respect_to, axis),
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
            Self::Modified(kernel) => {
                kernel.mixed_second_derivative(first, second, first_axis, second_axis)
            }
        }
    }
}
