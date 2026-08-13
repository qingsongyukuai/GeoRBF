//! QP/Modified-Kernel to ordinary-RBF linear reconstruction.
//!
//! Frozen sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `single_surface.cpp::convert_modified_kernel_to_rbf_kernel`;
//! - `lajaunie.cpp::convert_modified_kernel_to_rbf_kernel`;
//! - `stratigraphic_surfaces.cpp::convert_modified_kernel_to_rbf_kernel`.
//!
//! These functions are not called by the frozen public fitting path, but their
//! complete observable bodies are retained here as an explicit internal
//! operation. In particular, they convert every QP degree of freedom; there is
//! no multiplier/slack-based active-set filter in the source.

use std::fmt;

use crate::{
    assembly::{
        assemble_matrix_for_layout, functional_for_dof, polynomial_basis, AssemblyConstraints,
    },
    layout::section,
    solve_dense_partial_pivot_lu, solve_loqo_qp, solve_predictor_corrector_qp, AssembledSystem,
    AssemblyError, ConstraintLayout, Constraints, DenseMatrix, DenseVector, FunctionalKernel,
    IndexRange, Interface, InternalParameters, LayoutDof, LayoutPartitions, LayoutPointRef,
    LayoutRole, LayoutSectionKind, LinearFunctional, LoqoSolution, LoqoSolveError, LuSolution,
    LuSolveError, ModelType, Parameters, Point, QpSolution, QpSolveError, SolverType,
};

/// The solver that produced the Modified-Kernel weights being reconstructed.
#[derive(Clone, Debug, PartialEq)]
pub enum ReconstructionSourceSolution {
    PredictorCorrector(QpSolution),
    Loqo(LoqoSolution),
}

impl ReconstructionSourceSolution {
    pub const fn weights(&self) -> &DenseVector {
        match self {
            Self::PredictorCorrector(solution) => solution.weights(),
            Self::Loqo(solution) => solution.weights(),
        }
    }
}

/// Stable failure stages for callers that must distinguish the original QP,
/// reconstruction assembly, and the final ordinary LU.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReconstructionStage {
    SourceAssembly,
    Qp,
    Reassembly,
    Lu,
}

/// Safe reasons the frozen reconstruction body cannot be assembled.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ReconstructionAssemblyError {
    UnsupportedModel,
    NotQuadratic,
    SourceKernelNotModified,
    OrdinaryKernelIsModified,
    SourceLayoutMismatch,
    SourceWeightLengthMismatch,
    NonFinitePrediction,
    Assembly(AssemblyError),
}

impl fmt::Display for ReconstructionAssemblyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedModel => formatter.write_str(
                "QP-to-RBF reconstruction only exists for Single Surface, Lajaunie, and Stratigraphic Horizons",
            ),
            Self::NotQuadratic => {
                formatter.write_str("the source system is not a quadratic/Modified-Kernel branch")
            }
            Self::SourceKernelNotModified => {
                formatter.write_str("the source reconstruction kernel is not Modified-Kernel")
            }
            Self::OrdinaryKernelIsModified => {
                formatter.write_str("the reconstruction target kernel must be an ordinary RBF")
            }
            Self::SourceLayoutMismatch => {
                formatter.write_str("source assembly, model parameters, and constraints disagree")
            }
            Self::SourceWeightLengthMismatch => {
                formatter.write_str("QP weights do not cover exactly the source constraint layout")
            }
            Self::NonFinitePrediction => {
                formatter.write_str("the QP interpolant produced a non-finite reconstruction target")
            }
            Self::Assembly(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for ReconstructionAssemblyError {}

impl From<AssemblyError> for ReconstructionAssemblyError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}

/// Full stage-preserving reconstruction error.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum ReconstructionError {
    SourceAssembly(AssemblyError),
    PredictorCorrector(QpSolveError),
    Loqo(LoqoSolveError),
    Reassembly(ReconstructionAssemblyError),
    Lu(LuSolveError),
}

impl ReconstructionError {
    pub const fn stage(&self) -> ReconstructionStage {
        match self {
            Self::SourceAssembly(_) => ReconstructionStage::SourceAssembly,
            Self::PredictorCorrector(_) | Self::Loqo(_) => ReconstructionStage::Qp,
            Self::Reassembly(_) => ReconstructionStage::Reassembly,
            Self::Lu(_) => ReconstructionStage::Lu,
        }
    }
}

impl fmt::Display for ReconstructionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceAssembly(error) => write!(formatter, "source assembly failed: {error}"),
            Self::PredictorCorrector(error) => write!(formatter, "ordinary QP failed: {error}"),
            Self::Loqo(error) => write!(formatter, "restricted-range QP failed: {error}"),
            Self::Reassembly(error) => write!(formatter, "RBF reconstruction failed: {error}"),
            Self::Lu(error) => write!(formatter, "reconstructed LU failed: {error}"),
        }
    }
}

impl std::error::Error for ReconstructionError {}

/// One exact source-QP index to reconstructed equality index mapping.
#[derive(Clone, Debug, PartialEq)]
pub struct ReconstructionDofMapping {
    source_index: usize,
    reconstructed_index: usize,
    source_dof: LayoutDof,
    reconstructed_dof: LayoutDof,
    target_value: f64,
}

impl ReconstructionDofMapping {
    pub const fn source_index(&self) -> usize {
        self.source_index
    }

    pub const fn reconstructed_index(&self) -> usize {
        self.reconstructed_index
    }

    pub const fn source_dof(&self) -> &LayoutDof {
        &self.source_dof
    }

    pub const fn reconstructed_dof(&self) -> &LayoutDof {
        &self.reconstructed_dof
    }

    pub const fn target_value(&self) -> f64 {
        self.target_value
    }
}

/// Source Modified-Kernel and reconstructed ordinary-RBF predictions at one
/// caller-supplied witness point.
#[derive(Clone, Debug, PartialEq)]
pub struct ReconstructionPredictionWitness {
    point: [f64; 4],
    source_scalar: f64,
    source_gradient: [f64; 3],
    reconstructed_scalar: f64,
    reconstructed_gradient: [f64; 3],
}

impl ReconstructionPredictionWitness {
    pub const fn point(&self) -> [f64; 4] {
        self.point
    }

    pub const fn source_scalar(&self) -> f64 {
        self.source_scalar
    }

    pub const fn source_gradient(&self) -> [f64; 3] {
        self.source_gradient
    }

    pub const fn reconstructed_scalar(&self) -> f64 {
        self.reconstructed_scalar
    }

    pub const fn reconstructed_gradient(&self) -> [f64; 3] {
        self.reconstructed_gradient
    }
}

/// Complete ordinary-RBF system, LU result, and conversion evidence.
#[derive(Clone, Debug)]
pub struct ReconstructionResult {
    source_solution: Option<ReconstructionSourceSolution>,
    reconstructed_constraints: Constraints,
    layout: ConstraintLayout,
    interpolation_matrix: DenseMatrix,
    right_hand_side: DenseVector,
    smoothing_value: Option<f64>,
    lu_solution: LuSolution,
    mappings: Vec<ReconstructionDofMapping>,
    prediction_witnesses: Vec<ReconstructionPredictionWitness>,
}

impl ReconstructionResult {
    pub const fn source_solution(&self) -> Option<&ReconstructionSourceSolution> {
        self.source_solution.as_ref()
    }

    pub const fn reconstructed_constraints(&self) -> &Constraints {
        &self.reconstructed_constraints
    }

    pub const fn layout(&self) -> &ConstraintLayout {
        &self.layout
    }

    pub const fn interpolation_matrix(&self) -> &DenseMatrix {
        &self.interpolation_matrix
    }

    pub const fn right_hand_side(&self) -> &DenseVector {
        &self.right_hand_side
    }

    pub const fn smoothing_value(&self) -> Option<f64> {
        self.smoothing_value
    }

    pub const fn lu_solution(&self) -> &LuSolution {
        &self.lu_solution
    }

    pub fn mappings(&self) -> &[ReconstructionDofMapping] {
        &self.mappings
    }

    pub fn prediction_witnesses(&self) -> &[ReconstructionPredictionWitness] {
        &self.prediction_witnesses
    }
}

/// Assemble and solve the frozen source QP branch, then execute its otherwise
/// unreachable conversion to an ordinary-RBF linear system.
pub fn solve_and_reconstruct(
    constraints: &Constraints,
    parameters: &Parameters,
    modified_kernel: FunctionalKernel<'_>,
    ordinary_kernel: FunctionalKernel<'_>,
    witness_points: &[Point],
) -> Result<ReconstructionResult, ReconstructionError> {
    let source = crate::assemble_system(constraints, parameters, modified_kernel)
        .map_err(ReconstructionError::SourceAssembly)?;
    let source_solution = match source.constraints() {
        AssemblyConstraints::Quadratic {
            equality,
            inequality,
        } => ReconstructionSourceSolution::PredictorCorrector(
            solve_predictor_corrector_qp(source.interpolation_matrix(), equality, inequality)
                .map_err(ReconstructionError::PredictorCorrector)?,
        ),
        AssemblyConstraints::Bounded { system } => ReconstructionSourceSolution::Loqo(
            solve_loqo_qp(
                source.interpolation_matrix(),
                system.matrix(),
                system.lower(),
                system.range(),
            )
            .map_err(ReconstructionError::Loqo)?,
        ),
        AssemblyConstraints::Linear { .. } => {
            return Err(ReconstructionError::Reassembly(
                ReconstructionAssemblyError::NotQuadratic,
            ));
        }
    };
    let weights = source_solution.weights().clone();
    let mut result = reconstruct_from_qp_weights(
        constraints,
        parameters,
        &source,
        &weights,
        modified_kernel,
        ordinary_kernel,
        witness_points,
    )?;
    result.source_solution = Some(source_solution);
    Ok(result)
}

/// Reconstruct from already-solved QP weights. This lower-level entry point is
/// also the deterministic oracle seam for testing conversion independently of
/// QP convergence.
pub fn reconstruct_from_qp_weights(
    constraints: &Constraints,
    parameters: &Parameters,
    source: &AssembledSystem,
    qp_weights: &DenseVector,
    modified_kernel: FunctionalKernel<'_>,
    ordinary_kernel: FunctionalKernel<'_>,
    witness_points: &[Point],
) -> Result<ReconstructionResult, ReconstructionError> {
    validate_source(
        constraints,
        parameters,
        source,
        qp_weights,
        modified_kernel,
        ordinary_kernel,
    )
    .map_err(ReconstructionError::Reassembly)?;

    let source_functionals = source
        .layout()
        .dofs()
        .iter()
        .take(source.layout().constraint_dof_count())
        .map(|dof| functional_for_dof(dof, constraints))
        .collect::<Result<Vec<_>, _>>()
        .map_err(ReconstructionAssemblyError::from)
        .map_err(ReconstructionError::Reassembly)?;

    let source_fields = witness_points
        .iter()
        .map(|point| {
            evaluate_field(
                source.layout(),
                constraints,
                parameters,
                &source_functionals,
                qp_weights,
                modified_kernel,
                point,
            )
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(ReconstructionError::Reassembly)?;

    let prepared = prepare_reconstruction(
        constraints,
        parameters,
        source.layout(),
        &source_functionals,
        qp_weights,
        modified_kernel,
    )
    .map_err(ReconstructionError::Reassembly)?;

    let (interpolation_matrix, smoothing_value) = assemble_matrix_for_layout(
        &prepared.layout,
        &prepared.constraints,
        &prepared.parameters,
        ordinary_kernel,
    )
    .map_err(ReconstructionAssemblyError::from)
    .map_err(ReconstructionError::Reassembly)?;
    let right_hand_side = DenseVector::from_values(prepared.right_hand_side);
    let lu_solution = solve_dense_partial_pivot_lu(&interpolation_matrix, &right_hand_side)
        .map_err(ReconstructionError::Lu)?;

    let reconstructed_functionals = prepared
        .layout
        .dofs()
        .iter()
        .take(prepared.layout.constraint_dof_count())
        .map(|dof| functional_for_dof(dof, &prepared.constraints))
        .collect::<Result<Vec<_>, _>>()
        .map_err(ReconstructionAssemblyError::from)
        .map_err(ReconstructionError::Reassembly)?;
    let reconstructed_fields = witness_points
        .iter()
        .map(|point| {
            evaluate_field(
                &prepared.layout,
                &prepared.constraints,
                &prepared.parameters,
                &reconstructed_functionals,
                lu_solution.weights(),
                ordinary_kernel,
                point,
            )
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(ReconstructionError::Reassembly)?;
    let prediction_witnesses = witness_points
        .iter()
        .zip(source_fields)
        .zip(reconstructed_fields)
        .map(
            |((point, source), reconstructed)| ReconstructionPredictionWitness {
                point: [point.x(), point.y(), point.z(), point.c()],
                source_scalar: source.scalar,
                source_gradient: source.gradient,
                reconstructed_scalar: reconstructed.scalar,
                reconstructed_gradient: reconstructed.gradient,
            },
        )
        .collect();

    Ok(ReconstructionResult {
        source_solution: None,
        reconstructed_constraints: prepared.constraints,
        layout: prepared.layout,
        interpolation_matrix,
        right_hand_side,
        smoothing_value,
        lu_solution,
        mappings: prepared.mappings,
        prediction_witnesses,
    })
}

fn validate_source(
    constraints: &Constraints,
    parameters: &Parameters,
    source: &AssembledSystem,
    qp_weights: &DenseVector,
    modified_kernel: FunctionalKernel<'_>,
    ordinary_kernel: FunctionalKernel<'_>,
) -> Result<(), ReconstructionAssemblyError> {
    if !matches!(
        parameters.model_type,
        ModelType::SingleSurface | ModelType::LajaunieApproach | ModelType::StratigraphicHorizons
    ) {
        return Err(ReconstructionAssemblyError::UnsupportedModel);
    }
    if !matches!(modified_kernel, FunctionalKernel::Modified(_)) {
        return Err(ReconstructionAssemblyError::SourceKernelNotModified);
    }
    if matches!(ordinary_kernel, FunctionalKernel::Modified(_)) {
        return Err(ReconstructionAssemblyError::OrdinaryKernelIsModified);
    }
    let layout = source.layout();
    if layout.model() != parameters.model_type
        || !layout.internal_parameters().modified_basis
        || layout.polynomial_dof_count() != 0
        || layout.source_counts().inequalities != constraints.inequalities.len()
        || layout.source_counts().interfaces != constraints.interfaces.len()
        || layout.source_counts().planars != constraints.planars.len()
        || layout.source_counts().tangents != constraints.tangents.len()
    {
        return Err(ReconstructionAssemblyError::SourceLayoutMismatch);
    }
    if qp_weights.len() != layout.constraint_dof_count() {
        return Err(ReconstructionAssemblyError::SourceWeightLengthMismatch);
    }
    if !matches!(
        source.constraints(),
        AssemblyConstraints::Quadratic { .. } | AssemblyConstraints::Bounded { .. }
    ) {
        return Err(ReconstructionAssemblyError::NotQuadratic);
    }
    Ok(())
}

struct PreparedReconstruction {
    constraints: Constraints,
    parameters: Parameters,
    layout: ConstraintLayout,
    right_hand_side: Vec<f64>,
    mappings: Vec<ReconstructionDofMapping>,
}

fn prepare_reconstruction(
    constraints: &Constraints,
    parameters: &Parameters,
    source_layout: &ConstraintLayout,
    source_functionals: &[LinearFunctional],
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
) -> Result<PreparedReconstruction, ReconstructionAssemblyError> {
    let model = parameters.model_type;
    let mut reconstructed = constraints.clone();
    let mut reconstructed_parameters = parameters.clone();
    reconstructed_parameters.use_restricted_range = false;
    if model == ModelType::StratigraphicHorizons {
        // Stratigraphic reconstruction explicitly constructs a first-order
        // truncated basis. Lajaunie instead keeps `parameters.polynomial_order`
        // while forcing `n_poly_terms = 3`, so non-first-order Lajaunie input
        // must fail during reconstruction assembly just as the source does.
        reconstructed_parameters.polynomial_order = 1;
    }

    let target_values = source_layout
        .dofs()
        .iter()
        .take(source_layout.constraint_dof_count())
        .map(|dof| {
            target_for_dof(
                dof,
                constraints,
                parameters,
                source_layout,
                source_functionals,
                weights,
                kernel,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;

    for (index, planar) in reconstructed.planars.iter_mut().enumerate() {
        let point = constraints
            .planars
            .get(index)
            .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?
            .point();
        let field = evaluate_field(
            source_layout,
            constraints,
            parameters,
            source_functionals,
            weights,
            kernel,
            point,
        )?;
        planar
            .set_reconstructed_normal(field.gradient)
            .map_err(AssemblyError::Constraint)?;
    }
    for (index, tangent) in reconstructed.tangents.iter_mut().enumerate() {
        let original = constraints
            .tangents
            .get(index)
            .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?;
        let field = evaluate_field(
            source_layout,
            constraints,
            parameters,
            source_functionals,
            weights,
            kernel,
            original.point(),
        )?;
        let direction = original.vector();
        let inner = field.gradient[0] * direction[0]
            + field.gradient[1] * direction[1]
            + field.gradient[2] * direction[2];
        tangent
            .set_reconstructed_inner_product(inner)
            .map_err(AssemblyError::Constraint)?;
    }

    let (layout, reconstructed_dofs) = match model {
        ModelType::SingleSurface => {
            let mut interfaces =
                Vec::with_capacity(constraints.inequalities.len() + constraints.interfaces.len());
            for inequality in &constraints.inequalities {
                let field = evaluate_field(
                    source_layout,
                    constraints,
                    parameters,
                    source_functionals,
                    weights,
                    kernel,
                    inequality.point(),
                )?;
                interfaces.push(
                    Interface::with_c(
                        inequality.point().x(),
                        inequality.point().y(),
                        inequality.point().z(),
                        field.scalar,
                        inequality.point().c(),
                    )
                    .map_err(AssemblyError::Constraint)?,
                );
            }
            for interface in &constraints.interfaces {
                let field = evaluate_field(
                    source_layout,
                    constraints,
                    parameters,
                    source_functionals,
                    weights,
                    kernel,
                    interface.point(),
                )?;
                interfaces.push(
                    Interface::with_c(
                        interface.point().x(),
                        interface.point().y(),
                        interface.point().z(),
                        field.scalar,
                        interface.point().c(),
                    )
                    .map_err(AssemblyError::Constraint)?,
                );
            }
            reconstructed.inequalities.clear();
            reconstructed.interfaces = interfaces;
            let layout = crate::constraint_layout(
                ModelType::SingleSurface,
                &reconstructed,
                &reconstructed_parameters,
            )
            .map_err(AssemblyError::Surfe)?;
            let dofs = source_layout
                .dofs()
                .iter()
                .take(source_layout.constraint_dof_count())
                .map(|dof| match dof {
                    LayoutDof::InequalityValue { index } => {
                        LayoutDof::InterfaceValue { index: *index }
                    }
                    LayoutDof::InterfaceValue { index } => LayoutDof::InterfaceValue {
                        index: constraints.inequalities.len() + *index,
                    },
                    other => other.clone(),
                })
                .collect::<Vec<_>>();
            (layout, dofs)
        }
        ModelType::LajaunieApproach | ModelType::StratigraphicHorizons => {
            let layout = retained_difference_layout(source_layout, model);
            let dofs = source_layout.dofs()[..source_layout.constraint_dof_count()].to_vec();
            (layout, dofs)
        }
        ModelType::ContinuousProperty | ModelType::VectorField => {
            return Err(ReconstructionAssemblyError::UnsupportedModel);
        }
    };

    if reconstructed_dofs.len() != target_values.len()
        || layout.constraint_dof_count() != target_values.len()
    {
        return Err(ReconstructionAssemblyError::SourceLayoutMismatch);
    }
    let mappings = source_layout.dofs()[..source_layout.constraint_dof_count()]
        .iter()
        .cloned()
        .zip(reconstructed_dofs)
        .zip(target_values.iter().copied())
        .enumerate()
        .map(
            |(index, ((source_dof, reconstructed_dof), target_value))| ReconstructionDofMapping {
                source_index: index,
                reconstructed_index: index,
                source_dof,
                reconstructed_dof,
                target_value,
            },
        )
        .collect::<Vec<_>>();
    if mappings.iter().any(|mapping| {
        layout.dof(mapping.reconstructed_index()) != Some(mapping.reconstructed_dof())
    }) {
        return Err(ReconstructionAssemblyError::SourceLayoutMismatch);
    }
    let mut right_hand_side = target_values;
    right_hand_side.resize(layout.matrix_size(), 0.0);
    if !right_hand_side.iter().all(|value| value.is_finite()) {
        return Err(ReconstructionAssemblyError::NonFinitePrediction);
    }
    Ok(PreparedReconstruction {
        constraints: reconstructed,
        parameters: reconstructed_parameters,
        layout,
        right_hand_side,
        mappings,
    })
}

fn retained_difference_layout(source: &ConstraintLayout, model: ModelType) -> ConstraintLayout {
    let constraint_count = source.constraint_dof_count();
    let polynomial_count = 3;
    let mut dofs = source.dofs()[..constraint_count].to_vec();
    dofs.extend((0..polynomial_count).map(|index| LayoutDof::PolynomialTerm { index }));
    let mut roles = vec![LayoutRole::Equality; constraint_count];
    roles.extend(std::iter::repeat_n(
        LayoutRole::Polynomial,
        polynomial_count,
    ));
    let mut sections = source.sections().to_vec();
    sections.push(section(
        LayoutSectionKind::Polynomial,
        constraint_count,
        constraint_count + polynomial_count,
    ));
    let source_internal = source.internal_parameters();
    let internal = InternalParameters {
        n_interface: source_internal.n_interface,
        n_planar: source_internal.n_planar,
        // Frozen Stratigraphic conversion clears the two sequenced-pair
        // counters but leaves this derived field stale. Preserve it as
        // inspection evidence; the rebuilt layout roles/partitions below are
        // nevertheless all equality, exactly as the matrix body uses them.
        n_inequality: if model == ModelType::StratigraphicHorizons {
            source_internal.n_inequality
        } else {
            0
        },
        n_tangent: source_internal.n_tangent,
        n_constraints: constraint_count,
        n_equality: constraint_count,
        modified_basis: false,
        poly_term: true,
        n_poly_terms: polynomial_count,
        problem_type: SolverType::Linear,
        restricted_range: false,
    };
    ConstraintLayout::new(
        model,
        source.source_counts(),
        internal,
        dofs,
        roles,
        sections,
        LayoutPartitions::new(
            IndexRange::new(0, constraint_count),
            IndexRange::new(0, 0),
            IndexRange::new(0, 0),
            IndexRange::new(constraint_count, constraint_count + polynomial_count),
        ),
    )
}

fn target_for_dof(
    dof: &LayoutDof,
    constraints: &Constraints,
    parameters: &Parameters,
    layout: &ConstraintLayout,
    functionals: &[LinearFunctional],
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
) -> Result<f64, ReconstructionAssemblyError> {
    let value = match dof {
        LayoutDof::InequalityValue { index } => {
            let point = constraints
                .inequalities
                .get(*index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?
                .point();
            evaluate_field(
                layout,
                constraints,
                parameters,
                functionals,
                weights,
                kernel,
                point,
            )?
            .scalar
        }
        LayoutDof::InterfaceValue { index } => {
            let point = constraints
                .interfaces
                .get(*index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?
                .point();
            evaluate_field(
                layout,
                constraints,
                parameters,
                functionals,
                weights,
                kernel,
                point,
            )?
            .scalar
        }
        LayoutDof::Difference {
            positive, negative, ..
        } => {
            let positive = point_for_ref(*positive, constraints)?;
            let negative = point_for_ref(*negative, constraints)?;
            let first = evaluate_field(
                layout,
                constraints,
                parameters,
                functionals,
                weights,
                kernel,
                positive,
            )?
            .scalar;
            let second = evaluate_field(
                layout,
                constraints,
                parameters,
                functionals,
                weights,
                kernel,
                negative,
            )?
            .scalar;
            first - second
        }
        LayoutDof::PlanarDerivative { index, axis } => {
            let point = constraints
                .planars
                .get(*index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?
                .point();
            let gradient = evaluate_field(
                layout,
                constraints,
                parameters,
                functionals,
                weights,
                kernel,
                point,
            )?
            .gradient;
            gradient[axis_index(*axis)]
        }
        LayoutDof::Tangent { index } => {
            let tangent = constraints
                .tangents
                .get(*index)
                .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?;
            let gradient = evaluate_field(
                layout,
                constraints,
                parameters,
                functionals,
                weights,
                kernel,
                tangent.point(),
            )?
            .gradient;
            let direction = tangent.vector();
            gradient[0] * direction[0] + gradient[1] * direction[1] + gradient[2] * direction[2]
        }
        LayoutDof::PolynomialTerm { .. } => {
            return Err(ReconstructionAssemblyError::SourceLayoutMismatch);
        }
    };
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ReconstructionAssemblyError::NonFinitePrediction)
    }
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
            .map(|value| value.point()),
    }
    .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)
}

#[derive(Clone, Copy)]
pub(crate) struct FieldValue {
    pub(crate) scalar: f64,
    pub(crate) gradient: [f64; 3],
}

pub(crate) fn evaluate_field(
    layout: &ConstraintLayout,
    _constraints: &Constraints,
    parameters: &Parameters,
    functionals: &[LinearFunctional],
    weights: &DenseVector,
    kernel: FunctionalKernel<'_>,
    point: &Point,
) -> Result<FieldValue, ReconstructionAssemblyError> {
    if functionals.len() != layout.constraint_dof_count() || weights.len() != layout.matrix_size() {
        return Err(ReconstructionAssemblyError::SourceWeightLengthMismatch);
    }
    let planar_start = layout
        .section(LayoutSectionKind::PlanarDerivatives)
        .map_or(layout.constraint_dof_count(), IndexRange::start);
    let tangent_range = layout
        .section(LayoutSectionKind::Tangents)
        .unwrap_or(IndexRange::new(planar_start, planar_start));
    let planar_range = layout
        .section(LayoutSectionKind::PlanarDerivatives)
        .unwrap_or(IndexRange::new(planar_start, planar_start));

    let scalar_row = LinearFunctional::value(point.clone());
    let scalar = match layout.model() {
        ModelType::SingleSurface => {
            let inequalities = layout
                .section(LayoutSectionKind::InequalityValues)
                .unwrap_or(IndexRange::new(0, 0));
            let interfaces = layout
                .section(LayoutSectionKind::InterfaceValues)
                .unwrap_or(IndexRange::new(inequalities.end(), inequalities.end()));
            weighted_kernel_sum(kernel, &scalar_row, functionals, weights, inequalities)?
                + weighted_kernel_sum(kernel, &scalar_row, functionals, weights, interfaces)?
                + weighted_kernel_sum(kernel, &scalar_row, functionals, weights, planar_range)?
                + weighted_kernel_sum(kernel, &scalar_row, functionals, weights, tangent_range)?
                + polynomial_scalar(layout, parameters, weights, point)
        }
        ModelType::LajaunieApproach | ModelType::StratigraphicHorizons => {
            weighted_kernel_sum(
                kernel,
                &scalar_row,
                functionals,
                weights,
                IndexRange::new(0, planar_start),
            )? + weighted_kernel_sum(kernel, &scalar_row, functionals, weights, planar_range)?
                + weighted_kernel_sum(kernel, &scalar_row, functionals, weights, tangent_range)?
                + polynomial_scalar(layout, parameters, weights, point)
        }
        ModelType::ContinuousProperty | ModelType::VectorField => {
            return Err(ReconstructionAssemblyError::UnsupportedModel);
        }
    };

    let mut gradient = [0.0; 3];
    for axis in [crate::Axis::X, crate::Axis::Y, crate::Axis::Z] {
        let row = LinearFunctional::derivative(point.clone(), axis);
        let value = match layout.model() {
            ModelType::SingleSurface
            | ModelType::LajaunieApproach
            | ModelType::StratigraphicHorizons => {
                weighted_kernel_sum(
                    kernel,
                    &row,
                    functionals,
                    weights,
                    IndexRange::new(0, planar_start),
                )? + weighted_kernel_sum(kernel, &row, functionals, weights, planar_range)?
                    + weighted_kernel_sum(kernel, &row, functionals, weights, tangent_range)?
                    + polynomial_gradient_component(layout, parameters, weights, point, axis)
            }
            ModelType::ContinuousProperty | ModelType::VectorField => {
                return Err(ReconstructionAssemblyError::UnsupportedModel);
            }
        };
        gradient[axis_index(axis)] = value;
    }
    if scalar.is_finite() && gradient.into_iter().all(f64::is_finite) {
        Ok(FieldValue { scalar, gradient })
    } else {
        Err(ReconstructionAssemblyError::NonFinitePrediction)
    }
}

fn weighted_kernel_sum(
    kernel: FunctionalKernel<'_>,
    row: &LinearFunctional,
    columns: &[LinearFunctional],
    weights: &DenseVector,
    range: IndexRange,
) -> Result<f64, ReconstructionAssemblyError> {
    let mut sum = 0.0;
    for index in range.start()..range.end() {
        let column = columns
            .get(index)
            .ok_or(ReconstructionAssemblyError::SourceLayoutMismatch)?;
        let weight = weights
            .get(index)
            .ok_or(ReconstructionAssemblyError::SourceWeightLengthMismatch)?;
        sum += weight * kernel.apply(row, column).map_err(AssemblyError::Kernel)?;
    }
    Ok(sum)
}

fn polynomial_scalar(
    layout: &ConstraintLayout,
    parameters: &Parameters,
    weights: &DenseVector,
    point: &Point,
) -> f64 {
    let basis = polynomial_basis(layout.model(), parameters.polynomial_order);
    polynomial_dot(layout, weights, basis.values(point))
}

fn polynomial_gradient_component(
    layout: &ConstraintLayout,
    parameters: &Parameters,
    weights: &DenseVector,
    point: &Point,
    axis: crate::Axis,
) -> f64 {
    let basis = polynomial_basis(layout.model(), parameters.polynomial_order);
    let values = match axis {
        crate::Axis::X => basis.dx(point),
        crate::Axis::Y => basis.dy(point),
        crate::Axis::Z => basis.dz(point),
    };
    polynomial_dot(layout, weights, values)
}

fn polynomial_dot(layout: &ConstraintLayout, weights: &DenseVector, values: Vec<f64>) -> f64 {
    let range = layout.partitions().polynomial();
    if range.is_empty() {
        return 0.0;
    }
    debug_assert_eq!(range.len(), values.len());
    let mut sum = 0.0;
    for (offset, value) in values.into_iter().enumerate() {
        sum += value
            * weights
                .get(range.start() + offset)
                .expect("validated reconstructed polynomial weights");
    }
    sum
}

const fn axis_index(axis: crate::Axis) -> usize {
    match axis {
        crate::Axis::X => 0,
        crate::Axis::Y => 1,
        crate::Axis::Z => 2,
    }
}
