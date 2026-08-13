//! Frozen Continuous Property interface-value fitting and evaluation.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/continuous_property.{h,cpp}` (`Continuous_Property`
//!   constructors/destructor, `create_polynomial_basis`,
//!   `get_method_parameters`, `process_input_data`, matrix/RHS construction,
//!   `setup_system_solver`, `eval_*`, `measure_residuals`,
//!   `append_greedy_input`, and the TODO conversion body);
//! - `surfe_lib/modeling_methods.{h,cpp}` (`remove_collocated_constraints`,
//!   `setup_basis_functions`, the Greedy call chain, and the model factory);
//! - `surfe_lib/surfe_api.cpp` (both Continuous Property factory branches,
//!   `ComputeInterpolant`, and scalar/vector evaluation entry points).
//!
//! The frozen reachable and defined fit is deliberately narrow: only cleaned
//! interface values enter an ordinary isotropic RBF/LU system. Inequalities,
//! polynomial settings, smoothing, restricted range, and the unused Greedy
//! flag do not enter this model's active counts or equations. Planar and
//! tangent containers are more severe: the frozen RHS is allocated for
//! interface values and then writes those extra containers out of bounds.
//! GeoRBF reports that source defect before assembly instead of copying UB.
//! Polynomial helpers and Modified-Kernel conversion are not exposed because
//! their active flags are always false or their source body is a TODO. Greedy
//! residual/append hooks remain reserved for T31 and are not claimed here.

use std::fmt;

use crate::{
    assemble_system, solve_dense_partial_pivot_lu, AnisotropicKernel, AnisotropyError,
    AssembledSystem, AssemblyError, Axis, CollocationRemoval, Constraints, DenseMatrix,
    DenseVector, Error, FunctionalKernel, IsotropicKernel, KernelError, LinearFunctional,
    LuSolution, LuSolveError, ModelType, Parameters, Point,
};

pub(crate) mod assembly;
pub(crate) mod layout;

/// Failure from the actually reachable frozen Continuous Property path.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum ContinuousPropertyError {
    WrongModel,
    /// Frozen `get_equality_values` writes these ignored categories beyond the
    /// interface-sized `VectorXd`; valid Rust never performs that write.
    EqualityVectorOutOfBounds {
        planar_count: usize,
        tangent_count: usize,
    },
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Assembly(AssemblyError),
    Lu(LuSolveError),
    Evaluation(KernelError),
}

impl fmt::Display for ContinuousPropertyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => formatter.write_str("parameters do not select Continuous Property"),
            Self::EqualityVectorOutOfBounds {
                planar_count,
                tangent_count,
            } => write!(
                formatter,
                "frozen Continuous Property would write {planar_count} planar and {tangent_count} tangent constraints beyond its interface-sized equality vector"
            ),
            Self::Surfe(error) => error.fmt(formatter),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::Assembly(error) => error.fmt(formatter),
            Self::Lu(error) => error.fmt(formatter),
            Self::Evaluation(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for ContinuousPropertyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Surfe(error) => Some(error),
            Self::Anisotropy(error) => Some(error),
            Self::Assembly(error) => Some(error),
            Self::Lu(error) => Some(error),
            Self::Evaluation(error) => Some(error),
            Self::WrongModel | Self::EqualityVectorOutOfBounds { .. } => None,
        }
    }
}

impl From<AssemblyError> for ContinuousPropertyError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}

impl From<LuSolveError> for ContinuousPropertyError {
    fn from(error: LuSolveError) -> Self {
        Self::Lu(error)
    }
}

impl From<KernelError> for ContinuousPropertyError {
    fn from(error: KernelError) -> Self {
        Self::Evaluation(error)
    }
}

/// Immutable result of the frozen interface-value Continuous Property fit.
#[derive(Debug)]
pub struct ContinuousPropertyModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    kernel: IsotropicKernel,
    assembled: AssembledSystem,
    solution: LuSolution,
}

impl ContinuousPropertyModel {
    pub const fn parameters(&self) -> &Parameters {
        &self.parameters
    }

    /// Constraints after the four independent frozen cleaning passes.
    /// Inequalities remain visible but do not enter the system.
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
            .expect("Continuous Property stores only the linear branch")
    }

    pub const fn lu_solution(&self) -> &LuSolution {
        &self.solution
    }

    /// Evaluate the scalar field in the exact active summation order from
    /// `Continuous_Property::eval_scalar_interpolant_at_point`.
    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, ContinuousPropertyError> {
        let kernel = FunctionalKernel::from(&self.kernel);
        let query = LinearFunctional::value(point.clone());
        let weights = self.solution.weights().values();
        let mut interface_sum = 0.0;
        for (index, interface) in self.constraints.interfaces.iter().enumerate() {
            let source = LinearFunctional::value(interface.point().clone());
            interface_sum += weights[index] * kernel.apply(&query, &source)?;
        }
        Ok((((0.0 + interface_sum) + 0.0) + 0.0) + 0.0)
    }

    /// Evaluate the spatial gradient serially without the frozen shared-kernel
    /// data race. The scalar mathematics and component order are unchanged.
    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], ContinuousPropertyError> {
        let kernel = FunctionalKernel::from(&self.kernel);
        let weights = self.solution.weights().values();
        let queries = [Axis::X, Axis::Y, Axis::Z]
            .map(|axis| LinearFunctional::derivative(point.clone(), axis));
        let mut interface_sums = [0.0; 3];
        for (index, interface) in self.constraints.interfaces.iter().enumerate() {
            let source = LinearFunctional::value(interface.point().clone());
            for component in 0..3 {
                interface_sums[component] +=
                    weights[index] * kernel.apply(&queries[component], &source)?;
            }
        }
        Ok(interface_sums.map(|sum| ((sum + 0.0) + 0.0) + 0.0))
    }

    pub fn evaluate_scalars(&self, points: &[Point]) -> Result<Vec<f64>, ContinuousPropertyError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }

    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, ContinuousPropertyError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

/// Fit the frozen Continuous Property behavior reachable from
/// `Surfe_API::ComputeInterpolant` without reproducing source UB.
pub fn fit_continuous_property(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<ContinuousPropertyModel, ContinuousPropertyError> {
    if parameters.model_type != ModelType::ContinuousProperty {
        return Err(ContinuousPropertyError::WrongModel);
    }

    let mut constraints = constraints.clone();
    let collocation_removal = constraints.remove_collocated();
    if constraints.interfaces.is_empty() {
        return Err(ContinuousPropertyError::Surfe(Error::NoInterfaceData));
    }

    // Frozen setup constructs anisotropy before allocating/filling the RHS.
    // Preserve that failure precedence even though no successful anisotropic
    // Continuous Property model can pass the later out-of-bounds RHS branch.
    if parameters.model_global_anisotropy {
        AnisotropicKernel::new(
            parameters.basis_type,
            parameters.shape_parameter,
            &constraints.planars,
        )
        .map_err(ContinuousPropertyError::Anisotropy)?;
    }

    if !constraints.planars.is_empty() || !constraints.tangents.is_empty() {
        return Err(ContinuousPropertyError::EqualityVectorOutOfBounds {
            planar_count: constraints.planars.len(),
            tangent_count: constraints.tangents.len(),
        });
    }

    let kernel = IsotropicKernel::new(parameters.basis_type, parameters.shape_parameter);
    let assembled = assemble_system(&constraints, parameters, FunctionalKernel::from(&kernel))?;
    let right_hand_side = assembled
        .constraints()
        .linear_rhs()
        .expect("Continuous Property layout must select the linear branch");
    let solution = solve_dense_partial_pivot_lu(assembled.interpolation_matrix(), right_hand_side)?;

    Ok(ContinuousPropertyModel {
        parameters: parameters.clone(),
        constraints,
        collocation_removal,
        kernel,
        assembled,
        solution,
    })
}
