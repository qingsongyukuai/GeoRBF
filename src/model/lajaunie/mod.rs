//! Ordinary and restricted-range Lajaunie fitting and evaluation.
//!
//! Frozen sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/lajaunie.{h,cpp}` (`_get_increment_pairs`,
//!   `process_input_data`, `get_method_parameters`, matrix/RHS/bounds,
//!   `setup_system_solver`, `convert_modified_kernel_to_rbf_kernel`, and
//!   `eval_*`);
//! - `surfe_lib/modeling_methods.{h,cpp}` (`get_interface_data`,
//!   `_update_interface_iso_values`, cleaning, and basis construction);
//! - `surfe_lib/surfe_api.cpp` (`ComputeInterpolant` call order).

use std::fmt;

use crate::{
    assemble_system, reconstruct_from_qp_weights, solve_dense_partial_pivot_lu,
    solve_loqo_qp_with_options, AnisotropicKernel, AnisotropyError, AssembledSystem,
    AssemblyConstraints, AssemblyError, BoundedConstraintSystem, CollocationRemoval, Constraints,
    DenseMatrix, DenseVector, Error, FunctionalKernel, Interface, InterfaceGrouping,
    IsotropicKernel, LayoutDof, LoqoOptions, LoqoSolution, LoqoSolveError, LuSolution,
    LuSolveError, ModelType, ModifiedKernel, Parameters, Point, ReconstructionAssemblyError,
    ReconstructionError, ReconstructionResult, ReconstructionStage,
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

/// Updated scalar field value for one exact input interface level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LajaunieIsoValueEvidence {
    source_level: f64,
    reference_index: usize,
    iso_value: f64,
}

impl LajaunieIsoValueEvidence {
    /// Exact source level before `_update_interface_iso_values` overwrites the
    /// output level vector.
    pub const fn source_level(self) -> f64 {
        self.source_level
    }

    /// First cleaned interface point at this exact level.
    pub const fn reference_index(self) -> usize {
        self.reference_index
    }

    /// Fitted scalar field evaluated at the reference point.
    pub const fn iso_value(self) -> f64 {
        self.iso_value
    }
}

/// Failure from the frozen ordinary linear Lajaunie path.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum LajaunieLinearError {
    WrongModel,
    RestrictedRangeBranchNotAvailable,
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Assembly(AssemblyError),
    Lu(LuSolveError),
    Evaluation(ReconstructionAssemblyError),
}

impl fmt::Display for LajaunieLinearError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => formatter.write_str("parameters do not select Lajaunie"),
            Self::RestrictedRangeBranchNotAvailable => {
                formatter.write_str("Lajaunie restricted range requires the LOQO path")
            }
            Self::Surfe(error) => error.fmt(formatter),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::Assembly(error) => error.fmt(formatter),
            Self::Lu(error) => error.fmt(formatter),
            Self::Evaluation(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for LajaunieLinearError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Surfe(error) => Some(error),
            Self::Anisotropy(error) => Some(error),
            Self::Assembly(error) => Some(error),
            Self::Lu(error) => Some(error),
            Self::Evaluation(error) => Some(error),
            Self::WrongModel | Self::RestrictedRangeBranchNotAvailable => None,
        }
    }
}

/// Immutable result of the ordinary linear Lajaunie fit.
#[derive(Debug)]
pub struct LajaunieLinearModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    interface_grouping: InterfaceGrouping,
    kernel: OrdinaryKernel,
    assembled: AssembledSystem,
    solution: LuSolution,
    interface_iso_value_evidence: Vec<LajaunieIsoValueEvidence>,
    interface_iso_values: Vec<f64>,
}

impl LajaunieLinearModel {
    pub const fn parameters(&self) -> &Parameters {
        &self.parameters
    }

    /// Constraints after the four independent frozen sort/dedup passes.
    /// Inequalities remain observable here but do not enter a Lajaunie row.
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
            .expect("ordinary Lajaunie stores only the linear branch")
    }

    pub const fn lu_solution(&self) -> &LuSolution {
        &self.solution
    }

    pub fn interface_iso_value_evidence(&self) -> &[LajaunieIsoValueEvidence] {
        &self.interface_iso_value_evidence
    }

    pub fn interface_iso_values(&self) -> &[f64] {
        &self.interface_iso_values
    }

    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, LajaunieLinearError> {
        evaluate_layout_field(
            self.layout(),
            &self.constraints,
            &self.parameters,
            self.solution.weights(),
            self.kernel.functional(),
            point,
        )
        .map(|field| field.scalar)
        .map_err(LajaunieLinearError::Evaluation)
    }

    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], LajaunieLinearError> {
        evaluate_layout_field(
            self.layout(),
            &self.constraints,
            &self.parameters,
            self.solution.weights(),
            self.kernel.functional(),
            point,
        )
        .map(|field| field.gradient)
        .map_err(LajaunieLinearError::Evaluation)
    }

    pub fn evaluate_scalars(&self, points: &[Point]) -> Result<Vec<f64>, LajaunieLinearError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }

    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, LajaunieLinearError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

/// Fit the ordinary equality-only Lajaunie path.
///
/// Frozen Lajaunie sets `n_inequality = 0`; supplied inequality values are
/// therefore retained in the cleaned input but do not create QP rows.
pub fn fit_lajaunie_linear(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<LajaunieLinearModel, LajaunieLinearError> {
    if parameters.model_type != ModelType::LajaunieApproach {
        return Err(LajaunieLinearError::WrongModel);
    }
    if parameters.use_restricted_range {
        return Err(LajaunieLinearError::RestrictedRangeBranchNotAvailable);
    }

    let mut constraints = constraints.clone();
    let collocation_removal = constraints.remove_collocated();
    let interface_grouping =
        grouping_with_increment_pairs(&constraints).map_err(LajaunieLinearError::Surfe)?;
    let kernel = OrdinaryKernel::from_parameters(parameters, &constraints)
        .map_err(LajaunieLinearError::Anisotropy)?;
    let assembled = assemble_system(&constraints, parameters, kernel.functional())
        .map_err(LajaunieLinearError::Assembly)?;
    let right_hand_side = assembled
        .constraints()
        .linear_rhs()
        .expect("ordinary Lajaunie layout must select the linear branch");
    let solution = solve_dense_partial_pivot_lu(assembled.interpolation_matrix(), right_hand_side)
        .map_err(LajaunieLinearError::Lu)?;
    let (interface_iso_value_evidence, interface_iso_values) = update_interface_iso_values(
        &interface_grouping,
        assembled.layout(),
        &constraints,
        parameters,
        solution.weights(),
        kernel.functional(),
    )
    .map_err(LajaunieLinearError::Evaluation)?;

    Ok(LajaunieLinearModel {
        parameters: parameters.clone(),
        constraints,
        collocation_removal,
        interface_grouping,
        kernel,
        assembled,
        solution,
        interface_iso_value_evidence,
        interface_iso_values,
    })
}

/// Exact lower/range/upper evidence for one bounded Lajaunie source row.
#[derive(Clone, Debug, PartialEq)]
pub struct LajaunieRestrictedBoundEvidence {
    source_index: usize,
    dof: LayoutDof,
    lower: f64,
    range: f64,
}

impl LajaunieRestrictedBoundEvidence {
    pub const fn source_index(&self) -> usize {
        self.source_index
    }

    pub const fn dof(&self) -> &LayoutDof {
        &self.dof
    }

    pub const fn lower(&self) -> f64 {
        self.lower
    }

    pub const fn range(&self) -> f64 {
        self.range
    }

    pub fn upper(&self) -> f64 {
        self.lower + self.range
    }
}

/// Failure from Lajaunie Modified-Kernel/LOQO and explicit reconstruction.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum LajaunieRestrictedError {
    WrongModel,
    RestrictedRangeRequired,
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Basis(Error),
    SourceAssembly(AssemblyError),
    Loqo(LoqoSolveError),
    Reconstruction(ReconstructionError),
    Evaluation(ReconstructionAssemblyError),
}

impl LajaunieRestrictedError {
    pub const fn stage(&self) -> Option<ReconstructionStage> {
        match self {
            Self::SourceAssembly(_) => Some(ReconstructionStage::SourceAssembly),
            Self::Loqo(_) => Some(ReconstructionStage::Qp),
            Self::Reconstruction(error) => Some(error.stage()),
            Self::WrongModel
            | Self::RestrictedRangeRequired
            | Self::Surfe(_)
            | Self::Anisotropy(_)
            | Self::Basis(_)
            | Self::Evaluation(_) => None,
        }
    }
}

impl fmt::Display for LajaunieRestrictedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => formatter.write_str("parameters do not select Lajaunie"),
            Self::RestrictedRangeRequired => {
                formatter.write_str("Lajaunie restricted range must be enabled")
            }
            Self::Surfe(error) | Self::Basis(error) => error.fmt(formatter),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::SourceAssembly(error) => write!(formatter, "source assembly failed: {error}"),
            Self::Loqo(error) => write!(formatter, "restricted-range QP failed: {error}"),
            Self::Reconstruction(error) => error.fmt(formatter),
            Self::Evaluation(error) => write!(formatter, "field evaluation failed: {error}"),
        }
    }
}

impl std::error::Error for LajaunieRestrictedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Surfe(error) | Self::Basis(error) => Some(error),
            Self::Anisotropy(error) => Some(error),
            Self::SourceAssembly(error) => Some(error),
            Self::Loqo(error) => Some(error),
            Self::Reconstruction(error) => Some(error),
            Self::Evaluation(error) => Some(error),
            Self::WrongModel | Self::RestrictedRangeRequired => None,
        }
    }
}

/// Complete bounded Lajaunie result, including both source and reconstructed
/// iso-value updates.
#[derive(Debug)]
pub struct LajaunieRestrictedModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    interface_grouping: InterfaceGrouping,
    ordinary_kernel: OrdinaryKernel,
    modified_kernel: ModifiedKernel,
    source_assembled: AssembledSystem,
    loqo_solution: LoqoSolution,
    bound_evidence: Vec<LajaunieRestrictedBoundEvidence>,
    source_interface_iso_value_evidence: Vec<LajaunieIsoValueEvidence>,
    source_interface_iso_values: Vec<f64>,
    reconstruction: ReconstructionResult,
    interface_iso_value_evidence: Vec<LajaunieIsoValueEvidence>,
    interface_iso_values: Vec<f64>,
}

impl LajaunieRestrictedModel {
    pub const fn parameters(&self) -> &Parameters {
        &self.parameters
    }

    pub const fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub const fn collocation_removal(&self) -> CollocationRemoval {
        self.collocation_removal
    }

    pub const fn interface_grouping(&self) -> &InterfaceGrouping {
        &self.interface_grouping
    }

    pub const fn source_assembled_system(&self) -> &AssembledSystem {
        &self.source_assembled
    }

    pub const fn layout(&self) -> &crate::ConstraintLayout {
        self.source_assembled.layout()
    }

    pub const fn modified_interpolation_matrix(&self) -> &DenseMatrix {
        self.source_assembled.interpolation_matrix()
    }

    pub fn bounded_system(&self) -> &BoundedConstraintSystem {
        match self.source_assembled.constraints() {
            AssemblyConstraints::Bounded { system } => system,
            _ => unreachable!("restricted Lajaunie stores only the bounded branch"),
        }
    }

    pub const fn loqo_solution(&self) -> &LoqoSolution {
        &self.loqo_solution
    }

    pub fn bound_evidence(&self) -> &[LajaunieRestrictedBoundEvidence] {
        &self.bound_evidence
    }

    pub fn source_interface_iso_value_evidence(&self) -> &[LajaunieIsoValueEvidence] {
        &self.source_interface_iso_value_evidence
    }

    pub fn source_interface_iso_values(&self) -> &[f64] {
        &self.source_interface_iso_values
    }

    pub const fn reconstruction(&self) -> &ReconstructionResult {
        &self.reconstruction
    }

    pub fn interface_iso_value_evidence(&self) -> &[LajaunieIsoValueEvidence] {
        &self.interface_iso_value_evidence
    }

    pub fn interface_iso_values(&self) -> &[f64] {
        &self.interface_iso_values
    }

    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, LajaunieRestrictedError> {
        self.evaluate_reconstructed_field(point)
            .map(|field| field.scalar)
    }

    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], LajaunieRestrictedError> {
        self.evaluate_reconstructed_field(point)
            .map(|field| field.gradient)
    }

    pub fn evaluate_scalars(&self, points: &[Point]) -> Result<Vec<f64>, LajaunieRestrictedError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }

    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, LajaunieRestrictedError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }

    pub fn evaluate_modified_scalar(&self, point: &Point) -> Result<f64, LajaunieRestrictedError> {
        self.evaluate_modified_field(point)
            .map(|field| field.scalar)
    }

    pub fn evaluate_modified_gradient(
        &self,
        point: &Point,
    ) -> Result<[f64; 3], LajaunieRestrictedError> {
        self.evaluate_modified_field(point)
            .map(|field| field.gradient)
    }

    fn evaluate_modified_field(
        &self,
        point: &Point,
    ) -> Result<crate::model::reconstruct::FieldValue, LajaunieRestrictedError> {
        evaluate_layout_field(
            self.source_assembled.layout(),
            &self.constraints,
            &self.parameters,
            self.loqo_solution.weights(),
            FunctionalKernel::Modified(&self.modified_kernel),
            point,
        )
        .map_err(LajaunieRestrictedError::Evaluation)
    }

    fn evaluate_reconstructed_field(
        &self,
        point: &Point,
    ) -> Result<crate::model::reconstruct::FieldValue, LajaunieRestrictedError> {
        let mut parameters = self.parameters.clone();
        parameters.use_restricted_range = false;
        evaluate_layout_field(
            self.reconstruction.layout(),
            self.reconstruction.reconstructed_constraints(),
            &parameters,
            self.reconstruction.lu_solution().weights(),
            self.ordinary_kernel.functional(),
            point,
        )
        .map_err(LajaunieRestrictedError::Evaluation)
    }
}

pub fn fit_lajaunie_restricted(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<LajaunieRestrictedModel, LajaunieRestrictedError> {
    fit_lajaunie_restricted_with_options(constraints, parameters, LoqoOptions::default())
}

/// Fit with an explicit safety-only LOQO iteration cap, then execute the
/// complete frozen conversion body as an explicit operation.
pub fn fit_lajaunie_restricted_with_options(
    constraints: &Constraints,
    parameters: &Parameters,
    options: LoqoOptions,
) -> Result<LajaunieRestrictedModel, LajaunieRestrictedError> {
    if parameters.model_type != ModelType::LajaunieApproach {
        return Err(LajaunieRestrictedError::WrongModel);
    }
    if !parameters.use_restricted_range {
        return Err(LajaunieRestrictedError::RestrictedRangeRequired);
    }

    let mut constraints = constraints.clone();
    let collocation_removal = constraints.remove_collocated();
    let interface_grouping =
        grouping_with_increment_pairs(&constraints).map_err(LajaunieRestrictedError::Surfe)?;
    let ordinary_kernel = OrdinaryKernel::from_parameters(parameters, &constraints)
        .map_err(LajaunieRestrictedError::Anisotropy)?;
    let interface_point_lists = interface_point_lists(&constraints, &interface_grouping);
    let modified_kernel = ordinary_kernel
        .modified(&interface_point_lists)
        .map_err(LajaunieRestrictedError::Basis)?;
    let source_assembled = assemble_system(
        &constraints,
        parameters,
        FunctionalKernel::Modified(&modified_kernel),
    )
    .map_err(LajaunieRestrictedError::SourceAssembly)?;
    let bounded = match source_assembled.constraints() {
        AssemblyConstraints::Bounded { system } => system,
        _ => unreachable!("restricted Lajaunie must assemble the bounded branch"),
    };
    let loqo_solution = solve_loqo_qp_with_options(
        source_assembled.interpolation_matrix(),
        bounded.matrix(),
        bounded.lower(),
        bounded.range(),
        options,
    )
    .map_err(LajaunieRestrictedError::Loqo)?;
    let bound_evidence = source_assembled
        .layout()
        .dofs()
        .iter()
        .take(source_assembled.layout().constraint_dof_count())
        .cloned()
        .enumerate()
        .map(|(source_index, dof)| LajaunieRestrictedBoundEvidence {
            source_index,
            dof,
            lower: bounded
                .lower()
                .get(source_index)
                .expect("bounded lower vector matches Lajaunie source layout"),
            range: bounded
                .range()
                .get(source_index)
                .expect("bounded range vector matches Lajaunie source layout"),
        })
        .collect();
    let (source_interface_iso_value_evidence, source_interface_iso_values) =
        update_interface_iso_values(
            &interface_grouping,
            source_assembled.layout(),
            &constraints,
            parameters,
            loqo_solution.weights(),
            FunctionalKernel::Modified(&modified_kernel),
        )
        .map_err(LajaunieRestrictedError::Evaluation)?;

    let witnesses = reference_points(&constraints, &interface_grouping)
        .map_err(LajaunieRestrictedError::Evaluation)?;
    let reconstruction = reconstruct_from_qp_weights(
        &constraints,
        parameters,
        &source_assembled,
        loqo_solution.weights(),
        FunctionalKernel::Modified(&modified_kernel),
        ordinary_kernel.functional(),
        &witnesses,
    )
    .map_err(LajaunieRestrictedError::Reconstruction)?;
    let mut reconstructed_parameters = parameters.clone();
    reconstructed_parameters.use_restricted_range = false;
    let (interface_iso_value_evidence, interface_iso_values) = update_interface_iso_values(
        &interface_grouping,
        reconstruction.layout(),
        reconstruction.reconstructed_constraints(),
        &reconstructed_parameters,
        reconstruction.lu_solution().weights(),
        ordinary_kernel.functional(),
    )
    .map_err(LajaunieRestrictedError::Evaluation)?;

    Ok(LajaunieRestrictedModel {
        parameters: parameters.clone(),
        constraints,
        collocation_removal,
        interface_grouping,
        ordinary_kernel,
        modified_kernel,
        source_assembled,
        loqo_solution,
        bound_evidence,
        source_interface_iso_value_evidence,
        source_interface_iso_values,
        reconstruction,
        interface_iso_value_evidence,
        interface_iso_values,
    })
}

fn grouping_with_increment_pairs(constraints: &Constraints) -> Result<InterfaceGrouping, Error> {
    let grouping = constraints
        .interface_grouping()
        .ok_or(Error::NoInterfaceData)?;
    if grouping.increment_pair_count() == 0 {
        return Err(Error::NoInterfaceIncrementPairs);
    }
    Ok(grouping)
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

fn reference_points(
    constraints: &Constraints,
    grouping: &InterfaceGrouping,
) -> Result<Vec<Point>, ReconstructionAssemblyError> {
    grouping
        .reference_indices()
        .iter()
        .map(|index| {
            constraints
                .interfaces
                .get(*index)
                .map(|interface| interface.point().clone())
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)
        })
        .collect()
}

fn update_interface_iso_values(
    grouping: &InterfaceGrouping,
    layout: &crate::ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
) -> Result<(Vec<LajaunieIsoValueEvidence>, Vec<f64>), ReconstructionAssemblyError> {
    if grouping.levels_descending().len() != grouping.reference_indices().len()
        || grouping.reference_indices().is_empty()
    {
        return Err(ReconstructionAssemblyError::SourceLayoutMismatch);
    }
    let functionals = layout
        .dofs()
        .iter()
        .take(layout.constraint_dof_count())
        .map(|dof| crate::assembly::functional_for_dof(dof, constraints))
        .collect::<Result<Vec<_>, _>>()
        .map_err(ReconstructionAssemblyError::from)?;
    let evidence = grouping
        .levels_descending()
        .iter()
        .copied()
        .zip(grouping.reference_indices().iter().copied())
        .map(|(source_level, reference_index)| {
            let point = constraints
                .interfaces
                .get(reference_index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?
                .point();
            let iso_value = crate::model::reconstruct::evaluate_field(
                layout,
                constraints,
                parameters,
                &functionals,
                weights,
                kernel,
                point,
            )?
            .scalar;
            Ok(LajaunieIsoValueEvidence {
                source_level,
                reference_index,
                iso_value,
            })
        })
        .collect::<Result<Vec<_>, ReconstructionAssemblyError>>()?;
    let values = evidence.iter().map(|value| value.iso_value).collect();
    Ok((evidence, values))
}

fn evaluate_layout_field(
    layout: &crate::ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
    point: &Point,
) -> Result<crate::model::reconstruct::FieldValue, ReconstructionAssemblyError> {
    let functionals = layout
        .dofs()
        .iter()
        .take(layout.constraint_dof_count())
        .map(|dof| crate::assembly::functional_for_dof(dof, constraints))
        .collect::<Result<Vec<_>, _>>()
        .map_err(ReconstructionAssemblyError::from)?;
    crate::model::reconstruct::evaluate_field(
        layout,
        constraints,
        parameters,
        &functionals,
        weights,
        kernel,
        point,
    )
}
