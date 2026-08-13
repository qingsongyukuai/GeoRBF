//! Frozen Vector Field Hessian fitting, potential, and gradient evaluation.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/vector_field.{h,cpp}` (constructors,
//!   `get_method_parameters`, Hessian matrix/RHS construction,
//!   `setup_system_solver`, scalar-potential evaluation, and vector-gradient
//!   evaluation);
//! - `surfe_lib/modeling_methods.{h,cpp}` (`remove_collocated_constraints`,
//!   `setup_basis_functions`, and ordinary/anisotropic kernel factories);
//! - `surfe_lib/matrix_solver.{h,cpp}` (partial-pivot LU consumption);
//! - `surfe_lib/surfe_api.{h,cpp}` (Vector Field factory and public fit and
//!   evaluation call order).
//!
//! Only cleaned planar constraints enter this frozen model. Each planar owns
//! three derivative degrees of freedom in x/y/z order; their mixed Hessian is
//! solved against the stored normal. Interface, inequality, and tangent data
//! remain observable input but do not enter the system. Polynomial,
//! smoothing, restricted-range, and Greedy settings are likewise inactive:
//! `get_method_parameters` always selects a non-modified linear system with no
//! polynomial terms. The frozen empty-planar 0×0 solve succeeds and evaluates
//! to the zero potential/gradient, while a singleton Cubic system attempts LU
//! and fails. GeoRBF preserves both observable branches without copying
//! shared mutable kernel state.

use std::fmt;

use crate::{
    assemble_system, solve_dense_partial_pivot_lu, AnisotropicKernel, AnisotropyError,
    AssembledSystem, AssemblyError, Axis, CollocationRemoval, Constraints, DenseMatrix,
    DenseVector, FunctionalKernel, IsotropicKernel, KernelError, LinearFunctional, LuSolution,
    LuSolveError, ModelType, Parameters, Point,
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

/// Failure from the frozen Vector Field fit or evaluation path.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum VectorFieldError {
    WrongModel,
    Anisotropy(AnisotropyError),
    Assembly(AssemblyError),
    Lu(LuSolveError),
    Evaluation(KernelError),
}

impl fmt::Display for VectorFieldError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => formatter.write_str("parameters do not select Vector Field"),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::Assembly(error) => error.fmt(formatter),
            Self::Lu(error) => error.fmt(formatter),
            Self::Evaluation(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for VectorFieldError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Anisotropy(error) => Some(error),
            Self::Assembly(error) => Some(error),
            Self::Lu(error) => Some(error),
            Self::Evaluation(error) => Some(error),
            Self::WrongModel => None,
        }
    }
}

impl From<AssemblyError> for VectorFieldError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}

impl From<LuSolveError> for VectorFieldError {
    fn from(error: LuSolveError) -> Self {
        Self::Lu(error)
    }
}

impl From<KernelError> for VectorFieldError {
    fn from(error: KernelError) -> Self {
        Self::Evaluation(error)
    }
}

/// Immutable result of the frozen Vector Field fit.
#[derive(Debug)]
pub struct VectorFieldModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    kernel: OrdinaryKernel,
    assembled: AssembledSystem,
    // Frozen Eigen treats the empty 0×0 system as a successful empty solve.
    // T18 intentionally rejects generic empty systems, so this model records
    // that source-specific success without fabricating LU pivots/residuals.
    solution: Option<LuSolution>,
}

impl VectorFieldModel {
    pub const fn parameters(&self) -> &Parameters {
        &self.parameters
    }

    /// Constraints after the four independent frozen cleaning passes.
    /// Only `planars` enter the active system.
    pub const fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub const fn collocation_removal(&self) -> CollocationRemoval {
        self.collocation_removal
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
            .expect("Vector Field stores only the linear branch")
    }

    /// LU evidence for a non-empty system. `None` is the frozen successful
    /// empty-planar branch, not an unattempted failed fit.
    pub const fn lu_solution(&self) -> Option<&LuSolution> {
        self.solution.as_ref()
    }

    fn weights(&self) -> &[f64] {
        self.solution
            .as_ref()
            .map_or(&[], |solution| solution.weights().values())
    }

    /// Evaluate the scalar potential in the exact source summation order from
    /// `Vector_Field::eval_scalar_interpolant_at_point`.
    pub fn evaluate_potential(&self, point: &Point) -> Result<f64, VectorFieldError> {
        let weights = self.weights();
        let kernel = self.kernel.functional();
        let query = LinearFunctional::value(point.clone());
        let mut sum = 0.0;
        for (index, planar) in self.constraints.planars.iter().enumerate() {
            let offset = 3 * index;
            let dx = LinearFunctional::derivative(planar.point().clone(), Axis::X);
            let dy = LinearFunctional::derivative(planar.point().clone(), Axis::Y);
            let dz = LinearFunctional::derivative(planar.point().clone(), Axis::Z);
            sum += weights[offset] * kernel.apply(&query, &dx)?;
            sum += weights[offset + 1] * kernel.apply(&query, &dy)?;
            sum += weights[offset + 2] * kernel.apply(&query, &dz)?;
        }
        Ok(sum)
    }

    /// Evaluate the gradient of the fitted potential with the frozen Hessian
    /// row/column and accumulation order.
    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], VectorFieldError> {
        let weights = self.weights();
        let kernel = self.kernel.functional();
        let queries = [Axis::X, Axis::Y, Axis::Z]
            .map(|axis| LinearFunctional::derivative(point.clone(), axis));
        let mut sums = [0.0; 3];
        for (index, planar) in self.constraints.planars.iter().enumerate() {
            let offset = 3 * index;
            let sources = [Axis::X, Axis::Y, Axis::Z]
                .map(|axis| LinearFunctional::derivative(planar.point().clone(), axis));
            for row in 0..3 {
                sums[row] += weights[offset] * kernel.apply(&queries[row], &sources[0])?;
                sums[row] += weights[offset + 1] * kernel.apply(&queries[row], &sources[1])?;
                sums[row] += weights[offset + 2] * kernel.apply(&queries[row], &sources[2])?;
            }
        }
        Ok(sums)
    }

    pub fn evaluate_potentials(&self, points: &[Point]) -> Result<Vec<f64>, VectorFieldError> {
        points
            .iter()
            .map(|point| self.evaluate_potential(point))
            .collect()
    }

    pub fn evaluate_gradients(&self, points: &[Point]) -> Result<Vec<[f64; 3]>, VectorFieldError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

/// Fit the Vector Field behavior reachable from
/// `Surfe_API::ComputeInterpolant`.
pub fn fit_vector_field(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<VectorFieldModel, VectorFieldError> {
    if parameters.model_type != ModelType::VectorField {
        return Err(VectorFieldError::WrongModel);
    }

    let mut constraints = constraints.clone();
    let collocation_removal = constraints.remove_collocated();
    let kernel = OrdinaryKernel::from_parameters(parameters, &constraints)
        .map_err(VectorFieldError::Anisotropy)?;
    let assembled = assemble_system(&constraints, parameters, kernel.functional())?;
    let right_hand_side = assembled
        .constraints()
        .linear_rhs()
        .expect("Vector Field layout must select the linear branch");
    let solution = if assembled.interpolation_matrix().rows() == 0 {
        None
    } else {
        Some(solve_dense_partial_pivot_lu(
            assembled.interpolation_matrix(),
            right_hand_side,
        )?)
    };

    Ok(VectorFieldModel {
        parameters: parameters.clone(),
        constraints,
        collocation_removal,
        kernel,
        assembled,
        solution,
    })
}
