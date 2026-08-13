//! Ordinary linear and inequality/QP Single Surface fitting and evaluation.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/single_surface.{h,cpp}` (`process_input_data`,
//!   `get_method_parameters`, `get_interpolation_matrix`,
//!   `get_equality_values`, `get_inequality_matrix`, both
//!   `get_inequality_values` overloads, `setup_system_solver`, and `eval_*`);
//! - `surfe_lib/modeling_methods.{h,cpp}` (`get_interface_data`,
//!   `remove_collocated_constraints`, `setup_basis_functions`, and
//!   `check_interpolant`);
//! - `surfe_lib/matrix_solver.{h,cpp}`
//!   (`Quadratic_Predictor_Corrector` call chain);
//! - `surfe_lib/surfe_api.cpp` (`Surfe_API::ComputeInterpolant`).

use std::fmt;

use crate::{
    assemble_system, solve_dense_partial_pivot_lu, solve_predictor_corrector_qp_with_options,
    AnisotropicKernel, AnisotropyError, AssembledSystem, AssemblyConstraints, AssemblyError, Axis,
    CollocationRemoval, ConstraintSystem, Constraints, DenseMatrix, DenseVector, Error,
    FunctionalKernel, Interface, InterfaceGrouping, IsotropicKernel, KernelError, LinearFunctional,
    LuSolution, LuSolveError, ModelType, ModifiedKernel, Parameters, Point, QpOptions, QpSolution,
    QpSolveError,
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

    fn modified(self, interface_point_lists: &[Vec<Interface>]) -> Result<ModifiedKernel, Error> {
        match self {
            Self::Isotropic(kernel) => {
                ModifiedKernel::from_isotropic(kernel, interface_point_lists)
            }
            Self::Anisotropic(kernel) => {
                ModifiedKernel::from_anisotropic(kernel, interface_point_lists)
            }
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

/// Failure from the frozen ordinary-inequality Single Surface QP path.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum SingleSurfaceInequalityError {
    WrongModel,
    NoInequalities,
    RestrictedRangeBranchNotAvailable,
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Basis(Error),
    Assembly(AssemblyError),
    Qp(QpSolveError),
    Evaluation(KernelError),
}

impl fmt::Display for SingleSurfaceInequalityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => formatter.write_str("parameters do not select Single Surface"),
            Self::NoInequalities => {
                formatter.write_str("the Single Surface inequality path needs inequality data")
            }
            Self::RestrictedRangeBranchNotAvailable => {
                formatter.write_str("Single Surface restricted range requires the T24 LOQO path")
            }
            Self::Surfe(error) | Self::Basis(error) => error.fmt(formatter),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::Assembly(error) => error.fmt(formatter),
            Self::Qp(error) => error.fmt(formatter),
            Self::Evaluation(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for SingleSurfaceInequalityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Surfe(error) | Self::Basis(error) => Some(error),
            Self::Anisotropy(error) => Some(error),
            Self::Assembly(error) => Some(error),
            Self::Qp(error) => Some(error),
            Self::Evaluation(error) => Some(error),
            Self::WrongModel | Self::NoInequalities | Self::RestrictedRangeBranchNotAvailable => {
                None
            }
        }
    }
}

impl From<AssemblyError> for SingleSurfaceInequalityError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}

impl From<QpSolveError> for SingleSurfaceInequalityError {
    fn from(error: QpSolveError) -> Self {
        Self::Qp(error)
    }
}

impl From<KernelError> for SingleSurfaceInequalityError {
    fn from(error: KernelError) -> Self {
        Self::Evaluation(error)
    }
}

/// Terminal field and primal-dual evidence for one source inequality row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SingleSurfaceInequalityEvidence {
    source_index: usize,
    source_level: f64,
    row_sign: f64,
    scalar_field: f64,
    transformed_value: f64,
    matrix_slack: f64,
    solver_slack: f64,
    dual_multiplier: f64,
    active_tolerance: f64,
}

impl SingleSurfaceInequalityEvidence {
    pub const fn source_index(self) -> usize {
        self.source_index
    }

    pub const fn source_level(self) -> f64 {
        self.source_level
    }

    /// `+1` for a positive level and `-1` for a non-positive level.
    pub const fn row_sign(self) -> f64 {
        self.row_sign
    }

    pub const fn scalar_field(self) -> f64 {
        self.scalar_field
    }

    /// The frozen signed row value, `row_sign * scalar_field`, constrained
    /// to be greater than or equal to zero.
    pub const fn transformed_value(self) -> f64 {
        self.transformed_value
    }

    /// `C*x-d` evaluated from the assembled signed inequality row.
    pub const fn matrix_slack(self) -> f64 {
        self.matrix_slack
    }

    /// The terminal primal slack carried by the predictor-corrector solver.
    pub const fn solver_slack(self) -> f64 {
        self.solver_slack
    }

    pub const fn dual_multiplier(self) -> f64 {
        self.dual_multiplier
    }

    pub const fn active_tolerance(self) -> f64 {
        self.active_tolerance
    }

    /// Diagnostic active-bound classification using the solver's own final
    /// inequality acceptance threshold. Frozen Surfe does not filter rows by
    /// this classification.
    pub fn active_within_solver_tolerance(self) -> bool {
        self.matrix_slack.abs() <= self.active_tolerance
    }
}

/// Immutable result of the frozen ordinary predictor-corrector path.
#[derive(Debug)]
pub struct SingleSurfaceInequalityModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    interface_grouping: InterfaceGrouping,
    kernel: ModifiedKernel,
    assembled: AssembledSystem,
    solution: QpSolution,
    inequality_evidence: Vec<SingleSurfaceInequalityEvidence>,
}

impl SingleSurfaceInequalityModel {
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

    pub fn equality_system(&self) -> &ConstraintSystem {
        match self.assembled.constraints() {
            AssemblyConstraints::Quadratic { equality, .. } => equality,
            _ => unreachable!("T23 stores only the ordinary quadratic branch"),
        }
    }

    pub fn inequality_system(&self) -> &ConstraintSystem {
        match self.assembled.constraints() {
            AssemblyConstraints::Quadratic { inequality, .. } => inequality,
            _ => unreachable!("T23 stores only the ordinary quadratic branch"),
        }
    }

    pub const fn qp_solution(&self) -> &QpSolution {
        &self.solution
    }

    pub fn inequality_evidence(&self) -> &[SingleSurfaceInequalityEvidence] {
        &self.inequality_evidence
    }

    /// Evaluate the scalar field in the exact source category order:
    /// inequality, interface, planar, then tangent.
    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, SingleSurfaceInequalityError> {
        evaluate_qp_scalar(
            &self.kernel,
            &self.constraints,
            self.solution.weights(),
            point,
        )
        .map_err(Into::into)
    }

    /// Evaluate the spatial gradient in the frozen source summation order.
    pub fn evaluate_gradient(
        &self,
        point: &Point,
    ) -> Result<[f64; 3], SingleSurfaceInequalityError> {
        evaluate_qp_gradient(
            &self.kernel,
            &self.constraints,
            self.solution.weights(),
            point,
        )
        .map_err(Into::into)
    }

    pub fn evaluate_scalars(
        &self,
        points: &[Point],
    ) -> Result<Vec<f64>, SingleSurfaceInequalityError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }

    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, SingleSurfaceInequalityError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

/// Fit the ordinary Single Surface inequality/QP branch.
pub fn fit_single_surface_inequality(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<SingleSurfaceInequalityModel, SingleSurfaceInequalityError> {
    fit_single_surface_inequality_with_options(constraints, parameters, QpOptions::default())
}

/// Fit with an explicit safety-only QP iteration cap.
pub fn fit_single_surface_inequality_with_options(
    constraints: &Constraints,
    parameters: &Parameters,
    options: QpOptions,
) -> Result<SingleSurfaceInequalityModel, SingleSurfaceInequalityError> {
    if parameters.model_type != ModelType::SingleSurface {
        return Err(SingleSurfaceInequalityError::WrongModel);
    }

    let mut constraints = constraints.clone();
    let collocation_removal = constraints.remove_collocated();
    let interface_grouping = constraints
        .interface_grouping()
        .ok_or(SingleSurfaceInequalityError::Surfe(Error::NoInterfaceData))?;
    if parameters.use_restricted_range {
        return Err(SingleSurfaceInequalityError::RestrictedRangeBranchNotAvailable);
    }
    if constraints.inequalities.is_empty() {
        return Err(SingleSurfaceInequalityError::NoInequalities);
    }

    let radial = OrdinaryKernel::from_parameters(parameters, &constraints)
        .map_err(SingleSurfaceInequalityError::Anisotropy)?;
    let interface_point_lists = interface_point_lists(&constraints, &interface_grouping);
    let kernel = radial
        .modified(&interface_point_lists)
        .map_err(SingleSurfaceInequalityError::Basis)?;
    let assembled = assemble_system(
        &constraints,
        parameters,
        FunctionalKernel::Modified(&kernel),
    )?;
    let (equality, inequality) = match assembled.constraints() {
        AssemblyConstraints::Quadratic {
            equality,
            inequality,
        } => (equality, inequality),
        _ => unreachable!("inequality Single Surface must assemble the quadratic branch"),
    };
    let solution = solve_predictor_corrector_qp_with_options(
        assembled.interpolation_matrix(),
        equality,
        inequality,
        options,
    )?;
    let inequality_evidence =
        build_inequality_evidence(&constraints, &kernel, inequality, &solution)?;

    Ok(SingleSurfaceInequalityModel {
        parameters: parameters.clone(),
        constraints,
        collocation_removal,
        interface_grouping,
        kernel,
        assembled,
        solution,
        inequality_evidence,
    })
}

fn interface_point_lists(
    constraints: &Constraints,
    grouping: &InterfaceGrouping,
) -> Vec<Vec<Interface>> {
    grouping
        .multi_point_groups()
        .iter()
        .map(|indices| {
            indices
                .iter()
                .map(|index| constraints.interfaces[*index].clone())
                .collect()
        })
        .collect()
}

fn build_inequality_evidence(
    constraints: &Constraints,
    kernel: &ModifiedKernel,
    inequality_system: &ConstraintSystem,
    solution: &QpSolution,
) -> Result<Vec<SingleSurfaceInequalityEvidence>, KernelError> {
    let tolerance = solution.residual().inequality_limit();
    constraints
        .inequalities
        .iter()
        .enumerate()
        .map(|(index, inequality)| {
            let scalar_field =
                evaluate_qp_scalar(kernel, constraints, solution.weights(), inequality.point())?;
            let row_sign = if inequality.level() > 0.0 { 1.0 } else { -1.0 };
            let matrix_slack = inequality_system
                .matrix()
                .row(index)
                .expect("T23 inequality rows align with source inequalities")
                .iter()
                .zip(solution.weights().values())
                .fold(
                    -inequality_system.values().get(index).unwrap_or(0.0),
                    |sum, (a, x)| sum + a * x,
                );
            Ok(SingleSurfaceInequalityEvidence {
                source_index: index,
                source_level: inequality.level(),
                row_sign,
                scalar_field,
                transformed_value: row_sign * scalar_field,
                matrix_slack,
                solver_slack: solution.slack().get(index).unwrap_or(f64::NAN),
                dual_multiplier: solution.dual_inequality().get(index).unwrap_or(f64::NAN),
                active_tolerance: tolerance,
            })
        })
        .collect()
}

fn evaluate_qp_scalar(
    kernel: &ModifiedKernel,
    constraints: &Constraints,
    weights: &DenseVector,
    point: &Point,
) -> Result<f64, KernelError> {
    let kernel = FunctionalKernel::Modified(kernel);
    let values = weights.values();
    let query = LinearFunctional::value(point.clone());
    let mut inequality_sum = 0.0;
    for (index, inequality) in constraints.inequalities.iter().enumerate() {
        let functional = LinearFunctional::value(inequality.point().clone());
        inequality_sum += values[index] * kernel.apply(&query, &functional)?;
    }

    let interface_offset = constraints.inequalities.len();
    let mut interface_sum = 0.0;
    for (index, interface) in constraints.interfaces.iter().enumerate() {
        let functional = LinearFunctional::value(interface.point().clone());
        interface_sum += values[interface_offset + index] * kernel.apply(&query, &functional)?;
    }

    let planar_offset = interface_offset + constraints.interfaces.len();
    let mut planar_sum = 0.0;
    for (index, planar) in constraints.planars.iter().enumerate() {
        for (component, axis) in [Axis::X, Axis::Y, Axis::Z].into_iter().enumerate() {
            let functional = LinearFunctional::derivative(planar.point().clone(), axis);
            planar_sum += values[planar_offset + 3 * index + component]
                * kernel.apply(&query, &functional)?;
        }
    }

    let tangent_offset = planar_offset + 3 * constraints.planars.len();
    let mut tangent_sum = 0.0;
    for (index, tangent) in constraints.tangents.iter().enumerate() {
        let functional = LinearFunctional::tangent(tangent.clone());
        tangent_sum += values[tangent_offset + index] * kernel.apply(&query, &functional)?;
    }

    Ok(inequality_sum + interface_sum + planar_sum + tangent_sum + 0.0)
}

fn evaluate_qp_gradient(
    kernel: &ModifiedKernel,
    constraints: &Constraints,
    weights: &DenseVector,
    point: &Point,
) -> Result<[f64; 3], KernelError> {
    let axes = [Axis::X, Axis::Y, Axis::Z];
    let kernel = FunctionalKernel::Modified(kernel);
    let values = weights.values();
    let queries = axes.map(|axis| LinearFunctional::derivative(point.clone(), axis));
    let mut point_sums = [0.0; 3];
    for (index, inequality) in constraints.inequalities.iter().enumerate() {
        let functional = LinearFunctional::value(inequality.point().clone());
        for component in 0..3 {
            point_sums[component] +=
                values[index] * kernel.apply(&queries[component], &functional)?;
        }
    }
    let interface_offset = constraints.inequalities.len();
    for (index, interface) in constraints.interfaces.iter().enumerate() {
        let functional = LinearFunctional::value(interface.point().clone());
        for component in 0..3 {
            point_sums[component] += values[interface_offset + index]
                * kernel.apply(&queries[component], &functional)?;
        }
    }

    let planar_offset = interface_offset + constraints.interfaces.len();
    let mut planar_sums = [0.0; 3];
    for (index, planar) in constraints.planars.iter().enumerate() {
        for row_component in 0..3 {
            for (column_component, axis) in axes.into_iter().enumerate() {
                let functional = LinearFunctional::derivative(planar.point().clone(), axis);
                planar_sums[row_component] += values[planar_offset + 3 * index + column_component]
                    * kernel.apply(&queries[row_component], &functional)?;
            }
        }
    }

    let tangent_offset = planar_offset + 3 * constraints.planars.len();
    let mut tangent_sums = [0.0; 3];
    for (index, tangent) in constraints.tangents.iter().enumerate() {
        let functional = LinearFunctional::tangent(tangent.clone());
        for component in 0..3 {
            tangent_sums[component] +=
                values[tangent_offset + index] * kernel.apply(&queries[component], &functional)?;
        }
    }

    Ok(std::array::from_fn(|component| {
        point_sums[component] + planar_sums[component] + tangent_sums[component] + 0.0
    }))
}
