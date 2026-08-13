//! Dense matrix, right-hand-side, and regression-smoothing assembly.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - the five model `get_interpolation_matrix` and `get_*values` bodies;
//! - the Single/Lajaunie polynomial block helpers;
//! - `GRBF_Modelling_Methods::get_equality_matrix`;
//! - the Single/Lajaunie regression-smoothing branches.

use std::fmt::{self, Write};

use crate::{
    constraint_layout, model, Axis, ConstraintError, ConstraintLayout, Constraints, Error,
    FunctionalKernel, KernelError, LayoutDof, LayoutPointRef, LinearFunctional, ModelType,
    Parameters, Point, PolynomialBasis, PolynomialOrder,
};

/// Row-major, owned, pure-Rust dense matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct DenseMatrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl DenseMatrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows.saturating_mul(cols)],
        }
    }

    pub const fn rows(&self) -> usize {
        self.rows
    }

    pub const fn cols(&self) -> usize {
        self.cols
    }

    pub const fn dimensions(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.index(row, column).map(|index| self.data[index])
    }

    pub fn row(&self, row: usize) -> Option<&[f64]> {
        if row >= self.rows {
            return None;
        }
        let start = row * self.cols;
        Some(&self.data[start..start + self.cols])
    }

    /// Stable full-matrix evidence using exact binary64 bit patterns.
    pub fn debug_snapshot(&self) -> String {
        let mut snapshot = format!("{}x{} row-major bits", self.rows, self.cols);
        for row in 0..self.rows {
            snapshot.push('\n');
            for column in 0..self.cols {
                if column != 0 {
                    snapshot.push(' ');
                }
                write!(
                    snapshot,
                    "0x{:016x}",
                    self.get(row, column)
                        .expect("snapshot indices are in bounds")
                        .to_bits()
                )
                .expect("writing to String cannot fail");
            }
        }
        snapshot
    }

    pub(crate) fn set(&mut self, row: usize, column: usize, value: f64) {
        let index = self
            .index(row, column)
            .expect("assembly indices must match the T16 layout");
        self.data[index] = value;
    }

    pub(crate) fn rows_slice(&self, start: usize, end: usize) -> Self {
        debug_assert!(start <= end && end <= self.rows);
        let mut matrix = Self::zeros(end - start, self.cols);
        for row in start..end {
            let source = self.row(row).expect("validated row range");
            let destination_start = (row - start) * self.cols;
            matrix.data[destination_start..destination_start + self.cols].copy_from_slice(source);
        }
        matrix
    }

    fn index(&self, row: usize, column: usize) -> Option<usize> {
        (row < self.rows && column < self.cols).then_some(row * self.cols + column)
    }
}

/// Owned dense vector used by every T17 right-hand-side branch.
#[derive(Clone, Debug, PartialEq)]
pub struct DenseVector {
    values: Vec<f64>,
}

impl DenseVector {
    pub fn zeros(len: usize) -> Self {
        Self {
            values: vec![0.0; len],
        }
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn get(&self, index: usize) -> Option<f64> {
        self.values.get(index).copied()
    }

    pub(crate) fn from_values(values: Vec<f64>) -> Self {
        Self { values }
    }
}

/// One ordinary equality or inequality matrix and its value vector.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstraintSystem {
    matrix: DenseMatrix,
    values: DenseVector,
}

impl ConstraintSystem {
    pub(crate) const fn new(matrix: DenseMatrix, values: DenseVector) -> Self {
        Self { matrix, values }
    }

    pub const fn matrix(&self) -> &DenseMatrix {
        &self.matrix
    }

    pub const fn values(&self) -> &DenseVector {
        &self.values
    }
}

/// Frozen LOQO-style `lower <= A*x <= lower + range` inputs.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundedConstraintSystem {
    matrix: DenseMatrix,
    lower: DenseVector,
    range: DenseVector,
}

impl BoundedConstraintSystem {
    pub(crate) const fn new(matrix: DenseMatrix, lower: DenseVector, range: DenseVector) -> Self {
        Self {
            matrix,
            lower,
            range,
        }
    }

    pub const fn matrix(&self) -> &DenseMatrix {
        &self.matrix
    }

    pub const fn lower(&self) -> &DenseVector {
        &self.lower
    }

    pub const fn range(&self) -> &DenseVector {
        &self.range
    }
}

/// Solver-facing branch selected by the frozen model parameters.
#[derive(Clone, Debug, PartialEq)]
pub enum AssemblyConstraints {
    Linear {
        right_hand_side: DenseVector,
    },
    Quadratic {
        equality: ConstraintSystem,
        inequality: ConstraintSystem,
    },
    Bounded {
        system: BoundedConstraintSystem,
    },
}

impl AssemblyConstraints {
    pub const fn linear_rhs(&self) -> Option<&DenseVector> {
        match self {
            Self::Linear { right_hand_side } => Some(right_hand_side),
            _ => None,
        }
    }

    pub const fn bounded(&self) -> Option<&BoundedConstraintSystem> {
        match self {
            Self::Bounded { system } => Some(system),
            _ => None,
        }
    }

    const fn snapshot_name(&self) -> &'static str {
        match self {
            Self::Linear { .. } => "linear",
            Self::Quadratic { .. } => "quadratic",
            Self::Bounded { .. } => "bounded",
        }
    }
}

/// Complete T17 output without solving the system.
#[derive(Clone, Debug, PartialEq)]
pub struct AssembledSystem {
    layout: ConstraintLayout,
    interpolation_matrix: DenseMatrix,
    constraints: AssemblyConstraints,
    smoothing_value: Option<f64>,
}

impl AssembledSystem {
    pub const fn layout(&self) -> &ConstraintLayout {
        &self.layout
    }

    pub const fn interpolation_matrix(&self) -> &DenseMatrix {
        &self.interpolation_matrix
    }

    pub const fn constraints(&self) -> &AssemblyConstraints {
        &self.constraints
    }

    pub const fn smoothing_value(&self) -> Option<f64> {
        self.smoothing_value
    }

    /// Full layout, branch, smoothing, and binary64 matrix evidence.
    pub fn debug_snapshot(&self) -> String {
        let mut snapshot = self.layout.debug_snapshot();
        writeln!(snapshot).expect("writing to String cannot fail");
        writeln!(
            snapshot,
            "branch={} smoothing={}",
            self.constraints.snapshot_name(),
            match self.smoothing_value {
                Some(value) => format!("0x{:016x}", value.to_bits()),
                None => "none".to_owned(),
            }
        )
        .expect("writing to String cannot fail");
        snapshot.push_str(&self.interpolation_matrix.debug_snapshot());
        snapshot
    }
}

/// T17 failures retain either the frozen Surfe stage or a safe lower-level
/// reason that C++ represented with an integer sentinel or invalid state.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum AssemblyError {
    Surfe(Error),
    Kernel(KernelError),
    Constraint(ConstraintError),
    KernelLayoutMismatch,
}

impl fmt::Display for AssemblyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Surfe(error) => error.fmt(formatter),
            Self::Kernel(error) => error.fmt(formatter),
            Self::Constraint(error) => error.fmt(formatter),
            Self::KernelLayoutMismatch => formatter.write_str(
                "kernel kind does not match the modified-basis flag selected by the layout",
            ),
        }
    }
}

impl std::error::Error for AssemblyError {}

impl From<KernelError> for AssemblyError {
    fn from(error: KernelError) -> Self {
        Self::Kernel(error)
    }
}

impl From<ConstraintError> for AssemblyError {
    fn from(error: ConstraintError) -> Self {
        Self::Constraint(error)
    }
}

/// Assemble exactly one frozen model system without invoking LU or QP.
pub fn assemble_system(
    constraints: &Constraints,
    parameters: &Parameters,
    kernel: FunctionalKernel<'_>,
) -> Result<AssembledSystem, AssemblyError> {
    let layout = constraint_layout(parameters.model_type, constraints, parameters)
        .map_err(AssemblyError::Surfe)?;
    validate_kernel_kind(&layout, kernel)?;

    let constraint_count = layout.constraint_dof_count();
    let functionals = layout
        .dofs()
        .iter()
        .take(constraint_count)
        .map(|dof| functional_for_dof(dof, constraints))
        .collect::<Result<Vec<_>, _>>()?;
    let mut interpolation_matrix = DenseMatrix::zeros(layout.matrix_size(), layout.matrix_size());
    for (row, row_functional) in functionals.iter().enumerate() {
        for (column, column_functional) in functionals.iter().enumerate() {
            interpolation_matrix.set(
                row,
                column,
                kernel.apply(row_functional, column_functional)?,
            );
        }
    }
    insert_polynomial_blocks(
        &layout,
        constraints,
        parameters,
        &functionals,
        &mut interpolation_matrix,
    )?;
    let smoothing_value =
        apply_regression_smoothing(&layout, parameters, kernel, &mut interpolation_matrix)?;

    let assembled_constraints = match parameters.model_type {
        ModelType::SingleSurface => model::single_surface::assembly::build(
            &layout,
            constraints,
            parameters,
            &interpolation_matrix,
        )?,
        ModelType::LajaunieApproach => model::lajaunie::assembly::build(
            &layout,
            constraints,
            parameters,
            &interpolation_matrix,
        )?,
        ModelType::StratigraphicHorizons => model::stratigraphic::assembly::build(
            &layout,
            constraints,
            parameters,
            &interpolation_matrix,
        )?,
        ModelType::ContinuousProperty => {
            model::continuous_property::assembly::build(&layout, constraints)?
        }
        ModelType::VectorField => model::vector_field::assembly::build(&layout, constraints)?,
    };
    Ok(AssembledSystem {
        layout,
        interpolation_matrix,
        constraints: assembled_constraints,
        smoothing_value,
    })
}

fn validate_kernel_kind(
    layout: &ConstraintLayout,
    kernel: FunctionalKernel<'_>,
) -> Result<(), AssemblyError> {
    let is_modified = matches!(kernel, FunctionalKernel::Modified(_));
    if layout.internal_parameters().modified_basis == is_modified {
        Ok(())
    } else {
        Err(AssemblyError::KernelLayoutMismatch)
    }
}

fn functional_for_dof(
    dof: &LayoutDof,
    constraints: &Constraints,
) -> Result<LinearFunctional, AssemblyError> {
    let functional = match dof {
        LayoutDof::InequalityValue { index } => LinearFunctional::value(
            constraints
                .inequalities
                .get(*index)
                .ok_or(AssemblyError::Surfe(Error::InterpolationMatrixFailure))?
                .point()
                .clone(),
        ),
        LayoutDof::InterfaceValue { index } => LinearFunctional::value(
            constraints
                .interfaces
                .get(*index)
                .ok_or(AssemblyError::Surfe(Error::InterpolationMatrixFailure))?
                .point()
                .clone(),
        ),
        LayoutDof::Difference {
            positive, negative, ..
        } => LinearFunctional::difference(
            point_for_ref(*positive, constraints)?.clone(),
            point_for_ref(*negative, constraints)?.clone(),
        ),
        LayoutDof::PlanarDerivative { index, axis } => LinearFunctional::derivative(
            constraints
                .planars
                .get(*index)
                .ok_or(AssemblyError::Surfe(Error::InterpolationMatrixFailure))?
                .point()
                .clone(),
            *axis,
        ),
        LayoutDof::Tangent { index } => LinearFunctional::tangent(
            constraints
                .tangents
                .get(*index)
                .ok_or(AssemblyError::Surfe(Error::InterpolationMatrixFailure))?
                .clone(),
        ),
        LayoutDof::PolynomialTerm { .. } => {
            return Err(AssemblyError::Surfe(Error::InterpolationMatrixFailure));
        }
    };
    Ok(functional)
}

fn point_for_ref(
    reference: LayoutPointRef,
    constraints: &Constraints,
) -> Result<&Point, AssemblyError> {
    match reference {
        LayoutPointRef::Interface(index) => constraints
            .interfaces
            .get(index)
            .map(|value| value.point())
            .ok_or(AssemblyError::Surfe(Error::InterpolationMatrixFailure)),
        LayoutPointRef::Inequality(index) => constraints
            .inequalities
            .get(index)
            .map(|value| value.point())
            .ok_or(AssemblyError::Surfe(Error::InterpolationMatrixFailure)),
    }
}

fn insert_polynomial_blocks(
    layout: &ConstraintLayout,
    _constraints: &Constraints,
    parameters: &Parameters,
    functionals: &[LinearFunctional],
    matrix: &mut DenseMatrix,
) -> Result<(), AssemblyError> {
    let polynomial_count = layout.polynomial_dof_count();
    if polynomial_count == 0 {
        return Ok(());
    }
    let basis = polynomial_basis(parameters.model_type, parameters.polynomial_order);
    if basis.term_count() != polynomial_count {
        return Err(AssemblyError::Surfe(Error::InterpolationMatrixFailure));
    }
    let polynomial_start = layout.constraint_dof_count();
    for (column, functional) in functionals.iter().enumerate() {
        let values = polynomial_functional_values(basis, functional);
        if values.len() != polynomial_count {
            return Err(AssemblyError::Surfe(Error::InterpolationMatrixFailure));
        }
        for (term, value) in values.into_iter().enumerate() {
            matrix.set(polynomial_start + term, column, value);
            matrix.set(column, polynomial_start + term, value);
        }
    }
    Ok(())
}

fn polynomial_basis(model: ModelType, order: i32) -> PolynomialBasis {
    let order = match order {
        0 => PolynomialOrder::Zero,
        1 => PolynomialOrder::First,
        _ => PolynomialOrder::Second,
    };
    match model {
        ModelType::LajaunieApproach | ModelType::StratigraphicHorizons => {
            PolynomialBasis::truncated(order)
        }
        _ => PolynomialBasis::complete(order),
    }
}

fn polynomial_functional_values(basis: PolynomialBasis, functional: &LinearFunctional) -> Vec<f64> {
    match functional {
        LinearFunctional::Value(point) => basis.values(point),
        LinearFunctional::Derivative { point, axis } => polynomial_derivative(basis, point, *axis),
        LinearFunctional::Tangent(tangent) => {
            let direction = tangent.vector();
            let dx = basis.dx(tangent.point());
            let dy = basis.dy(tangent.point());
            let dz = basis.dz(tangent.point());
            dx.into_iter()
                .zip(dy)
                .zip(dz)
                .map(|((dx, dy), dz)| direction[0] * dx + direction[1] * dy + direction[2] * dz)
                .collect()
        }
        LinearFunctional::Difference { positive, negative } => basis
            .values(positive)
            .into_iter()
            .zip(basis.values(negative))
            .map(|(positive, negative)| positive - negative)
            .collect(),
    }
}

fn polynomial_derivative(basis: PolynomialBasis, point: &Point, axis: Axis) -> Vec<f64> {
    match axis {
        Axis::X => basis.dx(point),
        Axis::Y => basis.dy(point),
        Axis::Z => basis.dz(point),
    }
}

fn apply_regression_smoothing(
    layout: &ConstraintLayout,
    parameters: &Parameters,
    kernel: FunctionalKernel<'_>,
    matrix: &mut DenseMatrix,
) -> Result<Option<f64>, AssemblyError> {
    if !parameters.use_regression_smoothing {
        return Ok(None);
    }
    let diagonal_count = match layout.model() {
        ModelType::SingleSurface => {
            layout.internal_parameters().n_inequality + layout.internal_parameters().n_interface
        }
        ModelType::LajaunieApproach => layout
            .section(crate::LayoutSectionKind::SameLevelDifferences)
            .map_or(0, |range| range.len()),
        ModelType::StratigraphicHorizons
        | ModelType::ContinuousProperty
        | ModelType::VectorField => return Ok(None),
    };
    let origin = LinearFunctional::value(Point::new(0.0, 0.0, 0.0)?);
    let offset = LinearFunctional::value(Point::new(0.0, 0.0, parameters.smoothing_amount)?);
    let smoothing = kernel.apply(&origin, &offset)?;
    for index in 0..diagonal_count {
        matrix.set(index, index, smoothing);
    }
    Ok(Some(smoothing))
}

pub(crate) fn normal_component(normal: [f64; 3], axis: Axis) -> f64 {
    match axis {
        Axis::X => normal[0],
        Axis::Y => normal[1],
        Axis::Z => normal[2],
    }
}

pub(crate) fn rows_for_range(matrix: &DenseMatrix, range: crate::IndexRange) -> DenseMatrix {
    matrix.rows_slice(range.start(), range.end())
}
