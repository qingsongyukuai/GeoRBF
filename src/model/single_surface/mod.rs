//! Ordinary linear Single Surface fitting and evaluation.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/single_surface.{h,cpp}` (`process_input_data`,
//!   `get_method_parameters`, `get_interpolation_matrix`,
//!   `get_equality_values`, `setup_system_solver`, and `eval_*`);
//! - `surfe_lib/modeling_methods.{h,cpp}` (`get_interface_data`,
//!   `remove_collocated_constraints`, `setup_basis_functions`, and
//!   `check_interpolant`);
//! - `surfe_lib/surfe_api.cpp` (`Surfe_API::ComputeInterpolant`).

use std::fmt;

use crate::{
    assemble_system, solve_dense_partial_pivot_lu, AnisotropicKernel, AnisotropyError,
    AssembledSystem, AssemblyError, Axis, CollocationRemoval, Constraints, DenseMatrix,
    DenseVector, Error, FunctionalKernel, InterfaceGrouping, IsotropicKernel, KernelError,
    LinearFunctional, LuSolution, LuSolveError, ModelType, Parameters, Point,
};

pub(crate) mod assembly;
pub(crate) mod layout;

#[derive(Clone, Copy, Debug, PartialEq)]
enum OrdinaryKernel {
    Isotropic(IsotropicKernel),
    Anisotropic(AnisotropicKernel),
}

impl OrdinaryKernel {
    fn from_parameters(
        parameters: &Parameters,
        constraints: &Constraints,
    ) -> Result<Self, AnisotropyError> {
        if parameters.model_global_anisotropy {
            AnisotropicKernel::new(
                parameters.basis_type,
                parameters.shape_parameter,
                &constraints.planars,
            )
            .map(Self::Anisotropic)
        } else {
            Ok(Self::Isotropic(IsotropicKernel::new(
                parameters.basis_type,
                parameters.shape_parameter,
            )))
        }
    }

    const fn functional(&self) -> FunctionalKernel<'_> {
        match self {
            Self::Isotropic(kernel) => FunctionalKernel::Isotropic(kernel),
            Self::Anisotropic(kernel) => FunctionalKernel::Anisotropic(kernel),
        }
    }
}

/// Failure from the T22 ordinary-equality Single Surface path.
///
/// Inequality and restricted-range variants are scope guards: their complete
/// implementations belong to T23 and T24 respectively.
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub enum SingleSurfaceLinearError {
    WrongModel,
    InequalityBranchNotAvailable,
    RestrictedRangeBranchNotAvailable,
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Assembly(AssemblyError),
    Lu(LuSolveError),
    Evaluation(KernelError),
}

impl fmt::Display for SingleSurfaceLinearError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => formatter.write_str("parameters do not select Single Surface"),
            Self::InequalityBranchNotAvailable => {
                formatter.write_str("Single Surface inequalities require the T23 QP path")
            }
            Self::RestrictedRangeBranchNotAvailable => {
                formatter.write_str("Single Surface restricted range requires the T24 LOQO path")
            }
            Self::Surfe(error) => error.fmt(formatter),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::Assembly(error) => error.fmt(formatter),
            Self::Lu(error) => error.fmt(formatter),
            Self::Evaluation(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for SingleSurfaceLinearError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Surfe(error) => Some(error),
            Self::Anisotropy(error) => Some(error),
            Self::Assembly(error) => Some(error),
            Self::Lu(error) => Some(error),
            Self::Evaluation(error) => Some(error),
            Self::WrongModel
            | Self::InequalityBranchNotAvailable
            | Self::RestrictedRangeBranchNotAvailable => None,
        }
    }
}

impl From<AssemblyError> for SingleSurfaceLinearError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}

impl From<LuSolveError> for SingleSurfaceLinearError {
    fn from(error: LuSolveError) -> Self {
        Self::Lu(error)
    }
}

impl From<KernelError> for SingleSurfaceLinearError {
    fn from(error: KernelError) -> Self {
        Self::Evaluation(error)
    }
}

/// Immutable result of the frozen ordinary Single Surface fitting path.
#[derive(Debug)]
pub struct SingleSurfaceLinearModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    interface_grouping: InterfaceGrouping,
    kernel: OrdinaryKernel,
    assembled: AssembledSystem,
    solution: LuSolution,
}

impl SingleSurfaceLinearModel {
    pub const fn parameters(&self) -> &Parameters {
        &self.parameters
    }

    /// Constraints after the four independent frozen sort/dedup passes.
    pub const fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub const fn collocation_removal(&self) -> CollocationRemoval {
        self.collocation_removal
    }

    pub const fn interface_grouping(&self) -> &InterfaceGrouping {
        &self.interface_grouping
    }

    pub const fn assembled_system(&self) -> &AssembledSystem {
        &self.assembled
    }

    pub const fn layout(&self) -> &crate::ConstraintLayout {
        self.assembled.layout()
    }

    pub const fn interpolation_matrix(&self) -> &DenseMatrix {
        self.assembled.interpolation_matrix()
    }

    pub fn right_hand_side(&self) -> &DenseVector {
        self.assembled
            .constraints()
            .linear_rhs()
            .expect("T22 stores only the ordinary linear branch")
    }

    pub const fn lu_solution(&self) -> &LuSolution {
        &self.solution
    }

    /// Evaluate the scalar field with the frozen category-by-category sum
    /// order from `Single_Surface::eval_scalar_interpolant_at_point`.
    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, SingleSurfaceLinearError> {
        let kernel = self.kernel.functional();
        let weights = self.solution.weights().values();
        let query = LinearFunctional::value(point.clone());
        let mut interface_sum = 0.0;
        for (index, interface) in self.constraints.interfaces.iter().enumerate() {
            let functional = LinearFunctional::value(interface.point().clone());
            interface_sum += weights[index] * kernel.apply(&query, &functional)?;
        }

        let planar_offset = self.constraints.interfaces.len();
        let mut planar_sum = 0.0;
        for (index, planar) in self.constraints.planars.iter().enumerate() {
            for (component, axis) in [Axis::X, Axis::Y, Axis::Z].into_iter().enumerate() {
                let functional = LinearFunctional::derivative(planar.point().clone(), axis);
                planar_sum += weights[planar_offset + 3 * index + component]
                    * kernel.apply(&query, &functional)?;
            }
        }

        let tangent_offset = planar_offset + 3 * self.constraints.planars.len();
        let mut tangent_sum = 0.0;
        for (index, tangent) in self.constraints.tangents.iter().enumerate() {
            let functional = LinearFunctional::tangent(tangent.clone());
            tangent_sum += weights[tangent_offset + index] * kernel.apply(&query, &functional)?;
        }

        let polynomial_offset = tangent_offset + self.constraints.tangents.len();
        let basis = crate::assembly::polynomial_basis(
            ModelType::SingleSurface,
            self.parameters.polynomial_order,
        );
        let mut polynomial_sum = 0.0;
        for (index, value) in basis.values(point).into_iter().enumerate() {
            polynomial_sum += weights[polynomial_offset + index] * value;
        }

        Ok((((0.0 + interface_sum) + planar_sum) + tangent_sum) + polynomial_sum)
    }

    /// Evaluate the spatial gradient with the frozen component and summation
    /// order from `Single_Surface::eval_vector_interpolant_at_point`.
    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], SingleSurfaceLinearError> {
        let axes = [Axis::X, Axis::Y, Axis::Z];
        let kernel = self.kernel.functional();
        let weights = self.solution.weights().values();
        let queries = axes.map(|axis| LinearFunctional::derivative(point.clone(), axis));
        let mut interface_sums = [0.0; 3];
        for (index, interface) in self.constraints.interfaces.iter().enumerate() {
            let functional = LinearFunctional::value(interface.point().clone());
            for component in 0..3 {
                interface_sums[component] +=
                    weights[index] * kernel.apply(&queries[component], &functional)?;
            }
        }

        let planar_offset = self.constraints.interfaces.len();
        let mut planar_sums = [0.0; 3];
        for (index, planar) in self.constraints.planars.iter().enumerate() {
            for row_component in 0..3 {
                for (column_component, axis) in axes.into_iter().enumerate() {
                    let functional = LinearFunctional::derivative(planar.point().clone(), axis);
                    planar_sums[row_component] += weights
                        [planar_offset + 3 * index + column_component]
                        * kernel.apply(&queries[row_component], &functional)?;
                }
            }
        }

        let tangent_offset = planar_offset + 3 * self.constraints.planars.len();
        let mut tangent_sums = [0.0; 3];
        for (index, tangent) in self.constraints.tangents.iter().enumerate() {
            let functional = LinearFunctional::tangent(tangent.clone());
            for component in 0..3 {
                tangent_sums[component] += weights[tangent_offset + index]
                    * kernel.apply(&queries[component], &functional)?;
            }
        }

        let polynomial_offset = tangent_offset + self.constraints.tangents.len();
        let basis = crate::assembly::polynomial_basis(
            ModelType::SingleSurface,
            self.parameters.polynomial_order,
        );
        let derivatives = [basis.dx(point), basis.dy(point), basis.dz(point)];
        let mut polynomial_sums = [0.0; 3];
        for term in 0..derivatives[0].len() {
            for component in 0..3 {
                polynomial_sums[component] +=
                    derivatives[component][term] * weights[polynomial_offset + term];
            }
        }

        Ok(std::array::from_fn(|component| {
            (((0.0 + interface_sums[component]) + planar_sums[component]) + tangent_sums[component])
                + polynomial_sums[component]
        }))
    }

    pub fn evaluate_scalars(&self, points: &[Point]) -> Result<Vec<f64>, SingleSurfaceLinearError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }

    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, SingleSurfaceLinearError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

/// Fit the ordinary equality-only Single Surface path.
pub fn fit_single_surface_linear(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<SingleSurfaceLinearModel, SingleSurfaceLinearError> {
    if parameters.model_type != ModelType::SingleSurface {
        return Err(SingleSurfaceLinearError::WrongModel);
    }

    let mut constraints = constraints.clone();
    let collocation_removal = constraints.remove_collocated();
    let interface_grouping = constraints
        .interface_grouping()
        .ok_or(SingleSurfaceLinearError::Surfe(Error::NoInterfaceData))?;
    if parameters.use_restricted_range {
        return Err(SingleSurfaceLinearError::RestrictedRangeBranchNotAvailable);
    }
    if !constraints.inequalities.is_empty() {
        return Err(SingleSurfaceLinearError::InequalityBranchNotAvailable);
    }
    let kernel = OrdinaryKernel::from_parameters(parameters, &constraints)
        .map_err(SingleSurfaceLinearError::Anisotropy)?;
    let assembled = assemble_system(&constraints, parameters, kernel.functional())?;
    let right_hand_side = assembled
        .constraints()
        .linear_rhs()
        .expect("ordinary Single Surface layout must select the linear branch");
    let solution = solve_dense_partial_pivot_lu(assembled.interpolation_matrix(), right_hand_side)?;

    Ok(SingleSurfaceLinearModel {
        parameters: parameters.clone(),
        constraints,
        collocation_removal,
        interface_grouping,
        kernel,
        assembled,
        solution,
    })
}
