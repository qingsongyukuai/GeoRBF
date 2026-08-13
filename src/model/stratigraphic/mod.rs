//! Ordinary-QP and restricted-range Stratigraphic Horizons fitting.
//!
//! Frozen sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! `surfe_lib/stratigraphic_surfaces.{h,cpp}` (`_get_increment_pairs`,
//! lithostratigraphic neighbour selection, matrix/RHS/bounds, both solvers,
//! explicit Modified-Kernel conversion, iso-value update, and `eval_*`).

use std::fmt;

use crate::{
    assemble_system, constraint_layout, reconstruct_from_qp_weights, solve_loqo_qp_with_options,
    solve_predictor_corrector_qp_with_options, AnisotropicKernel, AnisotropyError, AssembledSystem,
    AssemblyConstraints, AssemblyError, BoundedConstraintSystem, CollocationRemoval,
    ConstraintLayout, ConstraintSystem, Constraints, DenseMatrix, DenseVector, Error,
    FunctionalKernel, IndexRange, Interface, InterfaceGrouping, IsotropicKernel, LayoutDof,
    LayoutPointRef, LayoutSectionKind, LoqoOptions, LoqoSolution, LoqoSolveError, ModelType,
    ModifiedKernel, Parameters, Point, QpOptions, QpSolution, QpSolveError,
    ReconstructionAssemblyError, ReconstructionError, ReconstructionResult, ReconstructionStage,
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

    fn modified(self, interface_lists: &[Vec<Interface>]) -> Result<ModifiedKernel, Error> {
        match self {
            Self::Isotropic(kernel) => ModifiedKernel::from_isotropic(kernel, interface_lists),
            Self::Anisotropic(kernel) => ModifiedKernel::from_anisotropic(kernel, interface_lists),
        }
    }
}

/// Fitted scalar at the first cleaned point of one exact horizon level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StratigraphicIsoValueEvidence {
    source_level: f64,
    reference_index: usize,
    iso_value: f64,
}

impl StratigraphicIsoValueEvidence {
    pub const fn source_level(self) -> f64 {
        self.source_level
    }

    pub const fn reference_index(self) -> usize {
        self.reference_index
    }

    pub const fn iso_value(self) -> f64 {
        self.iso_value
    }
}

/// Terminal fitted relation for a sequenced horizon or lithology row.
#[derive(Clone, Debug, PartialEq)]
pub struct StratigraphicLayerRelationEvidence {
    source_index: usize,
    dof: LayoutDof,
    minimum: f64,
    maximum: Option<f64>,
    increment: f64,
    matrix_value: f64,
    tolerance: f64,
}

impl StratigraphicLayerRelationEvidence {
    pub const fn source_index(&self) -> usize {
        self.source_index
    }

    pub const fn dof(&self) -> &LayoutDof {
        &self.dof
    }

    pub const fn minimum(&self) -> f64 {
        self.minimum
    }

    pub const fn maximum(&self) -> Option<f64> {
        self.maximum
    }

    /// Scalar field at the positive endpoint minus the negative endpoint.
    pub const fn increment(&self) -> f64 {
        self.increment
    }

    /// The corresponding assembled matrix row multiplied by the weights.
    pub const fn matrix_value(&self) -> f64 {
        self.matrix_value
    }

    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    pub fn accepted(&self) -> bool {
        self.increment + self.tolerance >= self.minimum
            && self
                .maximum
                .is_none_or(|maximum| self.increment - self.tolerance <= maximum)
    }
}

/// Failure from the ordinary Stratigraphic Horizons predictor-corrector path.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum StratigraphicError {
    WrongModel,
    RestrictedRangeBranchNotAvailable,
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Basis(Error),
    SourceAssembly(AssemblyError),
    Qp(QpSolveError),
    Reconstruction(ReconstructionError),
    Evaluation(ReconstructionAssemblyError),
}

impl StratigraphicError {
    pub const fn stage(&self) -> Option<ReconstructionStage> {
        match self {
            Self::SourceAssembly(_) => Some(ReconstructionStage::SourceAssembly),
            Self::Qp(_) => Some(ReconstructionStage::Qp),
            Self::Reconstruction(error) => Some(error.stage()),
            _ => None,
        }
    }
}

impl fmt::Display for StratigraphicError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => {
                formatter.write_str("parameters do not select Stratigraphic Horizons")
            }
            Self::RestrictedRangeBranchNotAvailable => {
                formatter.write_str("Stratigraphic restricted range requires the LOQO path")
            }
            Self::Surfe(error) | Self::Basis(error) => error.fmt(formatter),
            Self::Anisotropy(error) => error.fmt(formatter),
            Self::SourceAssembly(error) => write!(formatter, "source assembly failed: {error}"),
            Self::Qp(error) => write!(formatter, "ordinary QP failed: {error}"),
            Self::Reconstruction(error) => error.fmt(formatter),
            Self::Evaluation(error) => write!(formatter, "field evaluation failed: {error}"),
        }
    }
}

impl std::error::Error for StratigraphicError {}

/// Immutable ordinary-QP result plus the explicit frozen conversion body.
#[derive(Debug)]
pub struct StratigraphicModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    interface_grouping: InterfaceGrouping,
    ordinary_kernel: OrdinaryKernel,
    modified_kernel: ModifiedKernel,
    source_assembled: AssembledSystem,
    qp_solution: QpSolution,
    layer_relation_evidence: Vec<StratigraphicLayerRelationEvidence>,
    source_iso_value_evidence: Vec<StratigraphicIsoValueEvidence>,
    source_iso_values: Vec<f64>,
    reconstruction: ReconstructionResult,
    iso_value_evidence: Vec<StratigraphicIsoValueEvidence>,
    iso_values: Vec<f64>,
}

impl StratigraphicModel {
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
    pub const fn layout(&self) -> &ConstraintLayout {
        self.source_assembled.layout()
    }
    pub const fn modified_interpolation_matrix(&self) -> &DenseMatrix {
        self.source_assembled.interpolation_matrix()
    }
    pub const fn qp_solution(&self) -> &QpSolution {
        &self.qp_solution
    }
    pub fn equality_system(&self) -> &ConstraintSystem {
        match self.source_assembled.constraints() {
            AssemblyConstraints::Quadratic { equality, .. } => equality,
            _ => unreachable!("ordinary Stratigraphic stores the quadratic branch"),
        }
    }
    pub fn inequality_system(&self) -> &ConstraintSystem {
        match self.source_assembled.constraints() {
            AssemblyConstraints::Quadratic { inequality, .. } => inequality,
            _ => unreachable!("ordinary Stratigraphic stores the quadratic branch"),
        }
    }
    pub fn layer_relation_evidence(&self) -> &[StratigraphicLayerRelationEvidence] {
        &self.layer_relation_evidence
    }
    pub fn source_interface_iso_value_evidence(&self) -> &[StratigraphicIsoValueEvidence] {
        &self.source_iso_value_evidence
    }
    pub fn source_interface_iso_values(&self) -> &[f64] {
        &self.source_iso_values
    }
    pub const fn reconstruction(&self) -> &ReconstructionResult {
        &self.reconstruction
    }
    pub fn interface_iso_value_evidence(&self) -> &[StratigraphicIsoValueEvidence] {
        &self.iso_value_evidence
    }
    pub fn interface_iso_values(&self) -> &[f64] {
        &self.iso_values
    }

    pub fn evaluate_modified_scalar(&self, point: &Point) -> Result<f64, StratigraphicError> {
        evaluate_layout_field(
            self.layout(),
            &self.constraints,
            &self.parameters,
            self.qp_solution.weights(),
            FunctionalKernel::Modified(&self.modified_kernel),
            point,
        )
        .map(|field| field.scalar)
        .map_err(StratigraphicError::Evaluation)
    }
    pub fn evaluate_modified_gradient(
        &self,
        point: &Point,
    ) -> Result<[f64; 3], StratigraphicError> {
        evaluate_layout_field(
            self.layout(),
            &self.constraints,
            &self.parameters,
            self.qp_solution.weights(),
            FunctionalKernel::Modified(&self.modified_kernel),
            point,
        )
        .map(|field| field.gradient)
        .map_err(StratigraphicError::Evaluation)
    }
    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, StratigraphicError> {
        evaluate_reconstructed(
            self.reconstruction(),
            &self.parameters,
            self.ordinary_kernel.functional(),
            point,
        )
        .map(|field| field.scalar)
        .map_err(StratigraphicError::Evaluation)
    }
    pub fn evaluate_gradient(&self, point: &Point) -> Result<[f64; 3], StratigraphicError> {
        evaluate_reconstructed(
            self.reconstruction(),
            &self.parameters,
            self.ordinary_kernel.functional(),
            point,
        )
        .map(|field| field.gradient)
        .map_err(StratigraphicError::Evaluation)
    }
    pub fn evaluate_scalars(&self, points: &[Point]) -> Result<Vec<f64>, StratigraphicError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }
    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, StratigraphicError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

pub fn fit_stratigraphic(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<StratigraphicModel, StratigraphicError> {
    fit_stratigraphic_with_options(constraints, parameters, QpOptions::default())
}

/// Fit with an explicit safety-only ordinary-QP iteration cap.
pub fn fit_stratigraphic_with_options(
    constraints: &Constraints,
    parameters: &Parameters,
    options: QpOptions,
) -> Result<StratigraphicModel, StratigraphicError> {
    if parameters.model_type != ModelType::StratigraphicHorizons {
        return Err(StratigraphicError::WrongModel);
    }
    if parameters.use_restricted_range {
        return Err(StratigraphicError::RestrictedRangeBranchNotAvailable);
    }
    let prepared = prepare(constraints, parameters).map_err(map_prepare_ordinary)?;
    let source_assembled = assemble_system(
        &prepared.constraints,
        parameters,
        FunctionalKernel::Modified(&prepared.modified_kernel),
    )
    .map_err(StratigraphicError::SourceAssembly)?;
    let (equality, inequality) = match source_assembled.constraints() {
        AssemblyConstraints::Quadratic {
            equality,
            inequality,
        } => (equality, inequality),
        _ => unreachable!("ordinary Stratigraphic must assemble the quadratic branch"),
    };
    let qp_solution = solve_predictor_corrector_qp_with_options(
        source_assembled.interpolation_matrix(),
        equality,
        inequality,
        options,
    )
    .map_err(StratigraphicError::Qp)?;
    let source_kernel = FunctionalKernel::Modified(&prepared.modified_kernel);
    let layer_relation_evidence = ordinary_relation_evidence(
        source_assembled.layout(),
        &prepared.constraints,
        parameters,
        qp_solution.weights(),
        source_kernel,
        inequality,
        &qp_solution,
    )
    .map_err(StratigraphicError::Evaluation)?;
    let (source_iso_value_evidence, source_iso_values) = update_iso_values(
        &prepared.grouping,
        source_assembled.layout(),
        &prepared.constraints,
        parameters,
        qp_solution.weights(),
        source_kernel,
    )
    .map_err(StratigraphicError::Evaluation)?;
    let reconstruction = reconstruct(
        &prepared,
        parameters,
        &source_assembled,
        qp_solution.weights(),
        source_kernel,
    )
    .map_err(StratigraphicError::Reconstruction)?;
    let (iso_value_evidence, iso_values) =
        reconstructed_iso_values(&prepared, parameters, &reconstruction)
            .map_err(StratigraphicError::Evaluation)?;
    Ok(StratigraphicModel {
        parameters: parameters.clone(),
        constraints: prepared.constraints,
        collocation_removal: prepared.removal,
        interface_grouping: prepared.grouping,
        ordinary_kernel: prepared.ordinary_kernel,
        modified_kernel: prepared.modified_kernel,
        source_assembled,
        qp_solution,
        layer_relation_evidence,
        source_iso_value_evidence,
        source_iso_values,
        reconstruction,
        iso_value_evidence,
        iso_values,
    })
}

/// Exact lower/range/upper evidence for one bounded source row.
#[derive(Clone, Debug, PartialEq)]
pub struct StratigraphicRestrictedBoundEvidence {
    source_index: usize,
    dof: LayoutDof,
    lower: f64,
    range: f64,
}

impl StratigraphicRestrictedBoundEvidence {
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

#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum StratigraphicRestrictedError {
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

impl StratigraphicRestrictedError {
    pub const fn stage(&self) -> Option<ReconstructionStage> {
        match self {
            Self::SourceAssembly(_) => Some(ReconstructionStage::SourceAssembly),
            Self::Loqo(_) => Some(ReconstructionStage::Qp),
            Self::Reconstruction(error) => Some(error.stage()),
            _ => None,
        }
    }
}

impl fmt::Display for StratigraphicRestrictedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongModel => {
                formatter.write_str("parameters do not select Stratigraphic Horizons")
            }
            Self::RestrictedRangeRequired => {
                formatter.write_str("Stratigraphic restricted range must be enabled")
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

impl std::error::Error for StratigraphicRestrictedError {}

#[derive(Debug)]
pub struct StratigraphicRestrictedModel {
    parameters: Parameters,
    constraints: Constraints,
    collocation_removal: CollocationRemoval,
    interface_grouping: InterfaceGrouping,
    ordinary_kernel: OrdinaryKernel,
    modified_kernel: ModifiedKernel,
    source_assembled: AssembledSystem,
    loqo_solution: LoqoSolution,
    bound_evidence: Vec<StratigraphicRestrictedBoundEvidence>,
    layer_relation_evidence: Vec<StratigraphicLayerRelationEvidence>,
    source_iso_value_evidence: Vec<StratigraphicIsoValueEvidence>,
    source_iso_values: Vec<f64>,
    reconstruction: ReconstructionResult,
    iso_value_evidence: Vec<StratigraphicIsoValueEvidence>,
    iso_values: Vec<f64>,
}

impl StratigraphicRestrictedModel {
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
    pub const fn layout(&self) -> &ConstraintLayout {
        self.source_assembled.layout()
    }
    pub const fn modified_interpolation_matrix(&self) -> &DenseMatrix {
        self.source_assembled.interpolation_matrix()
    }
    pub fn bounded_system(&self) -> &BoundedConstraintSystem {
        match self.source_assembled.constraints() {
            AssemblyConstraints::Bounded { system } => system,
            _ => unreachable!("restricted Stratigraphic stores the bounded branch"),
        }
    }
    pub const fn loqo_solution(&self) -> &LoqoSolution {
        &self.loqo_solution
    }
    pub fn bound_evidence(&self) -> &[StratigraphicRestrictedBoundEvidence] {
        &self.bound_evidence
    }
    pub fn layer_relation_evidence(&self) -> &[StratigraphicLayerRelationEvidence] {
        &self.layer_relation_evidence
    }
    pub fn source_interface_iso_value_evidence(&self) -> &[StratigraphicIsoValueEvidence] {
        &self.source_iso_value_evidence
    }
    pub fn source_interface_iso_values(&self) -> &[f64] {
        &self.source_iso_values
    }
    pub const fn reconstruction(&self) -> &ReconstructionResult {
        &self.reconstruction
    }
    pub fn interface_iso_value_evidence(&self) -> &[StratigraphicIsoValueEvidence] {
        &self.iso_value_evidence
    }
    pub fn interface_iso_values(&self) -> &[f64] {
        &self.iso_values
    }
    pub fn evaluate_modified_scalar(
        &self,
        point: &Point,
    ) -> Result<f64, StratigraphicRestrictedError> {
        evaluate_layout_field(
            self.layout(),
            &self.constraints,
            &self.parameters,
            self.loqo_solution.weights(),
            FunctionalKernel::Modified(&self.modified_kernel),
            point,
        )
        .map(|field| field.scalar)
        .map_err(StratigraphicRestrictedError::Evaluation)
    }
    pub fn evaluate_modified_gradient(
        &self,
        point: &Point,
    ) -> Result<[f64; 3], StratigraphicRestrictedError> {
        evaluate_layout_field(
            self.layout(),
            &self.constraints,
            &self.parameters,
            self.loqo_solution.weights(),
            FunctionalKernel::Modified(&self.modified_kernel),
            point,
        )
        .map(|field| field.gradient)
        .map_err(StratigraphicRestrictedError::Evaluation)
    }
    pub fn evaluate_scalar(&self, point: &Point) -> Result<f64, StratigraphicRestrictedError> {
        evaluate_reconstructed(
            self.reconstruction(),
            &self.parameters,
            self.ordinary_kernel.functional(),
            point,
        )
        .map(|field| field.scalar)
        .map_err(StratigraphicRestrictedError::Evaluation)
    }
    pub fn evaluate_gradient(
        &self,
        point: &Point,
    ) -> Result<[f64; 3], StratigraphicRestrictedError> {
        evaluate_reconstructed(
            self.reconstruction(),
            &self.parameters,
            self.ordinary_kernel.functional(),
            point,
        )
        .map(|field| field.gradient)
        .map_err(StratigraphicRestrictedError::Evaluation)
    }
    pub fn evaluate_scalars(
        &self,
        points: &[Point],
    ) -> Result<Vec<f64>, StratigraphicRestrictedError> {
        points
            .iter()
            .map(|point| self.evaluate_scalar(point))
            .collect()
    }
    pub fn evaluate_gradients(
        &self,
        points: &[Point],
    ) -> Result<Vec<[f64; 3]>, StratigraphicRestrictedError> {
        points
            .iter()
            .map(|point| self.evaluate_gradient(point))
            .collect()
    }
}

pub fn fit_stratigraphic_restricted(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<StratigraphicRestrictedModel, StratigraphicRestrictedError> {
    fit_stratigraphic_restricted_with_options(constraints, parameters, LoqoOptions::default())
}

pub fn fit_stratigraphic_restricted_with_options(
    constraints: &Constraints,
    parameters: &Parameters,
    options: LoqoOptions,
) -> Result<StratigraphicRestrictedModel, StratigraphicRestrictedError> {
    if parameters.model_type != ModelType::StratigraphicHorizons {
        return Err(StratigraphicRestrictedError::WrongModel);
    }
    if !parameters.use_restricted_range {
        return Err(StratigraphicRestrictedError::RestrictedRangeRequired);
    }
    let prepared = prepare(constraints, parameters).map_err(map_prepare_restricted)?;
    let source_assembled = assemble_system(
        &prepared.constraints,
        parameters,
        FunctionalKernel::Modified(&prepared.modified_kernel),
    )
    .map_err(StratigraphicRestrictedError::SourceAssembly)?;
    let bounded = match source_assembled.constraints() {
        AssemblyConstraints::Bounded { system } => system,
        _ => unreachable!("restricted Stratigraphic must assemble the bounded branch"),
    };
    let loqo_solution = solve_loqo_qp_with_options(
        source_assembled.interpolation_matrix(),
        bounded.matrix(),
        bounded.lower(),
        bounded.range(),
        options,
    )
    .map_err(StratigraphicRestrictedError::Loqo)?;
    let bound_evidence = source_assembled
        .layout()
        .dofs()
        .iter()
        .take(source_assembled.layout().constraint_dof_count())
        .cloned()
        .enumerate()
        .map(|(source_index, dof)| StratigraphicRestrictedBoundEvidence {
            source_index,
            dof,
            lower: bounded
                .lower()
                .get(source_index)
                .expect("bounded lower matches layout"),
            range: bounded
                .range()
                .get(source_index)
                .expect("bounded range matches layout"),
        })
        .collect();
    let source_kernel = FunctionalKernel::Modified(&prepared.modified_kernel);
    let layer_relation_evidence = bounded_relation_evidence(
        source_assembled.layout(),
        &prepared.constraints,
        parameters,
        loqo_solution.weights(),
        source_kernel,
        bounded,
        &loqo_solution,
    )
    .map_err(StratigraphicRestrictedError::Evaluation)?;
    let (source_iso_value_evidence, source_iso_values) = update_iso_values(
        &prepared.grouping,
        source_assembled.layout(),
        &prepared.constraints,
        parameters,
        loqo_solution.weights(),
        source_kernel,
    )
    .map_err(StratigraphicRestrictedError::Evaluation)?;
    let reconstruction = reconstruct(
        &prepared,
        parameters,
        &source_assembled,
        loqo_solution.weights(),
        source_kernel,
    )
    .map_err(StratigraphicRestrictedError::Reconstruction)?;
    let (iso_value_evidence, iso_values) =
        reconstructed_iso_values(&prepared, parameters, &reconstruction)
            .map_err(StratigraphicRestrictedError::Evaluation)?;
    Ok(StratigraphicRestrictedModel {
        parameters: parameters.clone(),
        constraints: prepared.constraints,
        collocation_removal: prepared.removal,
        interface_grouping: prepared.grouping,
        ordinary_kernel: prepared.ordinary_kernel,
        modified_kernel: prepared.modified_kernel,
        source_assembled,
        loqo_solution,
        bound_evidence,
        layer_relation_evidence,
        source_iso_value_evidence,
        source_iso_values,
        reconstruction,
        iso_value_evidence,
        iso_values,
    })
}

struct Prepared {
    constraints: Constraints,
    removal: CollocationRemoval,
    grouping: InterfaceGrouping,
    ordinary_kernel: OrdinaryKernel,
    modified_kernel: ModifiedKernel,
}

#[derive(Clone, Debug)]
enum PrepareError {
    Surfe(Error),
    Anisotropy(AnisotropyError),
    Basis(Error),
}

fn prepare(constraints: &Constraints, parameters: &Parameters) -> Result<Prepared, PrepareError> {
    let mut constraints = constraints.clone();
    let removal = constraints.remove_collocated();
    let grouping = constraints
        .interface_grouping()
        .ok_or(PrepareError::Surfe(Error::NoInterfaceData))?;
    constraint_layout(ModelType::StratigraphicHorizons, &constraints, parameters)
        .map_err(PrepareError::Surfe)?;
    let ordinary_kernel = OrdinaryKernel::from_parameters(parameters, &constraints)
        .map_err(PrepareError::Anisotropy)?;
    let lists = interface_point_lists(&constraints, &grouping);
    let modified_kernel = ordinary_kernel
        .modified(&lists)
        .map_err(PrepareError::Basis)?;
    Ok(Prepared {
        constraints,
        removal,
        grouping,
        ordinary_kernel,
        modified_kernel,
    })
}

fn map_prepare_ordinary(error: PrepareError) -> StratigraphicError {
    match error {
        PrepareError::Surfe(e) => StratigraphicError::Surfe(e),
        PrepareError::Anisotropy(e) => StratigraphicError::Anisotropy(e),
        PrepareError::Basis(e) => StratigraphicError::Basis(e),
    }
}

fn map_prepare_restricted(error: PrepareError) -> StratigraphicRestrictedError {
    match error {
        PrepareError::Surfe(e) => StratigraphicRestrictedError::Surfe(e),
        PrepareError::Anisotropy(e) => StratigraphicRestrictedError::Anisotropy(e),
        PrepareError::Basis(e) => StratigraphicRestrictedError::Basis(e),
    }
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

fn reconstruct(
    prepared: &Prepared,
    parameters: &Parameters,
    source: &AssembledSystem,
    weights: &DenseVector,
    source_kernel: FunctionalKernel<'_>,
) -> Result<ReconstructionResult, ReconstructionError> {
    let witnesses = reference_points(&prepared.constraints, &prepared.grouping)
        .map_err(ReconstructionError::Reassembly)?;
    reconstruct_from_qp_weights(
        &prepared.constraints,
        parameters,
        source,
        weights,
        source_kernel,
        prepared.ordinary_kernel.functional(),
        &witnesses,
    )
}

fn reconstructed_iso_values(
    prepared: &Prepared,
    parameters: &Parameters,
    reconstruction: &ReconstructionResult,
) -> Result<(Vec<StratigraphicIsoValueEvidence>, Vec<f64>), ReconstructionAssemblyError> {
    let mut parameters = parameters.clone();
    parameters.use_restricted_range = false;
    parameters.polynomial_order = 1;
    update_iso_values(
        &prepared.grouping,
        reconstruction.layout(),
        reconstruction.reconstructed_constraints(),
        &parameters,
        reconstruction.lu_solution().weights(),
        prepared.ordinary_kernel.functional(),
    )
}

fn update_iso_values(
    grouping: &InterfaceGrouping,
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
) -> Result<(Vec<StratigraphicIsoValueEvidence>, Vec<f64>), ReconstructionAssemblyError> {
    if grouping.levels_descending().len() != grouping.reference_indices().len()
        || grouping.reference_indices().is_empty()
    {
        return Err(ReconstructionAssemblyError::SourceLayoutMismatch);
    }
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
            let iso_value =
                evaluate_layout_field(layout, constraints, parameters, weights, kernel, point)?
                    .scalar;
            Ok(StratigraphicIsoValueEvidence {
                source_level,
                reference_index,
                iso_value,
            })
        })
        .collect::<Result<Vec<_>, ReconstructionAssemblyError>>()?;
    let values = evidence.iter().map(|value| value.iso_value).collect();
    Ok((evidence, values))
}

fn relation_range(layout: &ConstraintLayout) -> IndexRange {
    let first = layout
        .section(LayoutSectionKind::SequencedInterfaceDifferences)
        .unwrap_or(IndexRange::new(0, 0));
    let second = layout
        .section(LayoutSectionKind::SequencedInequalityDifferences)
        .unwrap_or(IndexRange::new(first.end(), first.end()));
    IndexRange::new(first.start(), second.end())
}

fn ordinary_relation_evidence(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
    inequality: &ConstraintSystem,
    solution: &QpSolution,
) -> Result<Vec<StratigraphicLayerRelationEvidence>, ReconstructionAssemblyError> {
    relation_evidence(
        layout,
        constraints,
        parameters,
        weights,
        kernel,
        RelationSystem {
            matrix: inequality.matrix(),
            lower: inequality.values(),
            range: None,
            tolerance: solution.residual().inequality_limit(),
        },
    )
}

fn bounded_relation_evidence(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
    bounded: &BoundedConstraintSystem,
    solution: &LoqoSolution,
) -> Result<Vec<StratigraphicLayerRelationEvidence>, ReconstructionAssemblyError> {
    relation_evidence(
        layout,
        constraints,
        parameters,
        weights,
        kernel,
        RelationSystem {
            matrix: bounded.matrix(),
            lower: bounded.lower(),
            range: Some(bounded.range()),
            tolerance: solution.residual().feasibility_limit(),
        },
    )
}

struct RelationSystem<'a> {
    matrix: &'a DenseMatrix,
    lower: &'a DenseVector,
    range: Option<&'a DenseVector>,
    tolerance: f64,
}

fn relation_evidence(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
    system: RelationSystem<'_>,
) -> Result<Vec<StratigraphicLayerRelationEvidence>, ReconstructionAssemblyError> {
    let offset = relation_range(layout).start();
    layout.dofs()[relation_range(layout).start()..relation_range(layout).end()]
        .iter()
        .cloned()
        .enumerate()
        .map(|(local_index, dof)| {
            let source_index = offset + local_index;
            let (positive, negative) = match &dof {
                LayoutDof::Difference {
                    positive, negative, ..
                } => (*positive, *negative),
                _ => return Err(ReconstructionAssemblyError::SourceLayoutMismatch),
            };
            let positive = point_for_ref(positive, constraints)?;
            let negative = point_for_ref(negative, constraints)?;
            let increment =
                evaluate_layout_field(layout, constraints, parameters, weights, kernel, positive)?
                    .scalar
                    - evaluate_layout_field(
                        layout,
                        constraints,
                        parameters,
                        weights,
                        kernel,
                        negative,
                    )?
                    .scalar;
            let row_index = if system.range.is_some() {
                source_index
            } else {
                local_index
            };
            let matrix_value = system
                .matrix
                .row(row_index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?
                .iter()
                .zip(weights.values())
                .fold(0.0, |sum, (a, x)| sum + a * x);
            let minimum = system
                .lower
                .get(row_index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?;
            let maximum = system
                .range
                .map(|values| {
                    values
                        .get(row_index)
                        .map(|range| minimum + range)
                        .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)
                })
                .transpose()?;
            Ok(StratigraphicLayerRelationEvidence {
                source_index,
                dof,
                minimum,
                maximum,
                increment,
                matrix_value,
                tolerance: system.tolerance,
            })
        })
        .collect()
}

fn point_for_ref(
    reference: LayoutPointRef,
    constraints: &Constraints,
) -> Result<&Point, ReconstructionAssemblyError> {
    match reference {
        LayoutPointRef::Interface(index) => constraints.interfaces.get(index).map(Interface::point),
        LayoutPointRef::Inequality(index) => constraints
            .inequalities
            .get(index)
            .map(crate::Inequality::point),
    }
    .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)
}

fn evaluate_reconstructed(
    reconstruction: &ReconstructionResult,
    parameters: &Parameters,
    kernel: FunctionalKernel<'_>,
    point: &Point,
) -> Result<crate::model::reconstruct::FieldValue, ReconstructionAssemblyError> {
    let mut parameters = parameters.clone();
    parameters.use_restricted_range = false;
    parameters.polynomial_order = 1;
    evaluate_layout_field(
        reconstruction.layout(),
        reconstruction.reconstructed_constraints(),
        &parameters,
        reconstruction.lu_solution().weights(),
        kernel,
        point,
    )
}

fn evaluate_layout_field(
    layout: &ConstraintLayout,
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
