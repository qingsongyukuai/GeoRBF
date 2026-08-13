//! Restricted-range primal-dual LOQO-style QP matching frozen Surfe.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/matrix_solver.{h,cpp}::Quadratic_Predictor_Corrector_LOQO`;
//! - `math_lib/math_methods.{h,cpp}::quadratic_solver_loqo`;
//! - `Math_methods::{_find_step,_find_positivity_step}`;
//! - `surfe_lib/surfe_api.{h,cpp}::SetRestrictedRange`.

use std::fmt;

use crate::{DenseMatrix, DenseVector, Error};

use super::{solve_dense_partial_pivot_lu, LuResidualEvidence, LuSolveError};

const SIGNIFICANT_FIGURES_TARGET: f64 = 6.0;
const RESIDUAL_ABSOLUTE_TOLERANCE: f64 = 1.0e-7;
const RESIDUAL_RELATIVE_TOLERANCE: f64 = 1.0e-6;
const FEASIBILITY_ABSOLUTE_TOLERANCE: f64 = 1.0e-10;
const FEASIBILITY_RELATIVE_TOLERANCE: f64 = 1.0e-8;

/// Safety-only iteration cap around the source loop, which has no cap.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LoqoOptions {
    pub max_iterations: usize,
}

impl Default for LoqoOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
        }
    }
}

/// Shape, finite-input, and frozen LLT-validation evidence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LoqoValidation {
    variables_non_empty: bool,
    interpolation_square: bool,
    constraint_rows_match: bool,
    constraint_columns_match: bool,
    lower_values_match: bool,
    range_values_match: bool,
    hessian_finite: bool,
    all_inputs_finite: bool,
    frozen_llt_positive_definite: bool,
}

impl LoqoValidation {
    pub const fn variables_are_non_empty(self) -> bool {
        self.variables_non_empty
    }

    pub const fn interpolation_is_square(self) -> bool {
        self.interpolation_square
    }

    pub const fn constraint_rows_match(self) -> bool {
        self.constraint_rows_match
    }

    pub const fn constraint_columns_match(self) -> bool {
        self.constraint_columns_match
    }

    pub const fn lower_values_match(self) -> bool {
        self.lower_values_match
    }

    pub const fn range_values_match(self) -> bool {
        self.range_values_match
    }

    pub const fn hessian_is_finite(self) -> bool {
        self.hessian_finite
    }

    pub const fn all_inputs_are_finite(self) -> bool {
        self.all_inputs_finite
    }

    /// Observable result of frozen `validate_matrix_systems()` on `H=2K`.
    pub const fn surfe_matrix_system_valid(self) -> bool {
        self.hessian_finite && self.frozen_llt_positive_definite
    }

    /// Safe dimensions only; the source LLT validator is not a solve gate.
    pub const fn safe_shape_valid(self) -> bool {
        self.variables_non_empty
            && self.interpolation_square
            && self.constraint_rows_match
            && self.constraint_columns_match
            && self.lower_values_match
            && self.range_values_match
    }
}

/// One source-level KKT factorization stage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LoqoKktStage {
    Initial,
    Predictor(usize),
    Corrector(usize),
}

/// Compact evidence from one successful pure-Rust partial-pivot KKT solve.
#[derive(Clone, Debug, PartialEq)]
pub struct LoqoKktSolveEvidence {
    stage: LoqoKktStage,
    dimension: usize,
    row_transpositions: Vec<usize>,
    residual: LuResidualEvidence,
}

impl LoqoKktSolveEvidence {
    pub const fn stage(&self) -> LoqoKktStage {
        self.stage
    }

    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn row_transpositions(&self) -> &[usize] {
        &self.row_transpositions
    }

    pub const fn residual(&self) -> LuResidualEvidence {
        self.residual
    }
}

/// A failed KKT stage retains the complete T18 LU failure.
#[derive(Clone, Debug, PartialEq)]
pub struct LoqoKktFailure {
    stage: LoqoKktStage,
    source: LuSolveError,
}

impl LoqoKktFailure {
    pub const fn stage(&self) -> LoqoKktStage {
        self.stage
    }

    pub const fn source(&self) -> &LuSolveError {
        &self.source
    }
}

/// Per-iteration evidence at the source print point.
#[derive(Clone, Debug, PartialEq)]
pub struct LoqoIterationEvidence {
    iteration: usize,
    primal_objective: f64,
    dual_objective: f64,
    significant_figures: f64,
    primal_infeasibility: f64,
    dual_infeasibility: f64,
    predictor_primal_divisor: Option<f64>,
    predictor_dual_divisor: Option<f64>,
    predictor_fraction: Option<f64>,
    predictor_mu: Option<f64>,
    corrector_primal_divisor: Option<f64>,
    corrector_dual_divisor: Option<f64>,
}

impl LoqoIterationEvidence {
    pub const fn iteration(&self) -> usize {
        self.iteration
    }

    pub const fn primal_objective(&self) -> f64 {
        self.primal_objective
    }

    pub const fn dual_objective(&self) -> f64 {
        self.dual_objective
    }

    pub const fn significant_figures(&self) -> f64 {
        self.significant_figures
    }

    pub const fn primal_infeasibility(&self) -> f64 {
        self.primal_infeasibility
    }

    pub const fn dual_infeasibility(&self) -> f64 {
        self.dual_infeasibility
    }

    pub const fn predictor_primal_divisor(&self) -> Option<f64> {
        self.predictor_primal_divisor
    }

    pub const fn predictor_dual_divisor(&self) -> Option<f64> {
        self.predictor_dual_divisor
    }

    pub const fn predictor_fraction(&self) -> Option<f64> {
        self.predictor_fraction
    }

    pub const fn predictor_mu(&self) -> Option<f64> {
        self.predictor_mu
    }

    pub const fn corrector_primal_divisor(&self) -> Option<f64> {
        self.corrector_primal_divisor
    }

    pub const fn corrector_dual_divisor(&self) -> Option<f64> {
        self.corrector_dual_divisor
    }
}

/// The only success stop in the frozen LOQO loop.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LoqoStopReason {
    SignificantFigures,
}

/// Terminal KKT residual, box feasibility, gap, and complementarity evidence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LoqoResidualEvidence {
    primal_infeasibility: f64,
    dual_infeasibility: f64,
    minimum_lower_slack: f64,
    minimum_upper_slack: f64,
    significant_figures: f64,
    complementarity: f64,
    residual_limit: f64,
    feasibility_limit: f64,
    finite: bool,
    accepted: bool,
}

impl LoqoResidualEvidence {
    pub const fn primal_infeasibility(self) -> f64 {
        self.primal_infeasibility
    }

    pub const fn dual_infeasibility(self) -> f64 {
        self.dual_infeasibility
    }

    pub const fn minimum_lower_slack(self) -> f64 {
        self.minimum_lower_slack
    }

    pub const fn minimum_upper_slack(self) -> f64 {
        self.minimum_upper_slack
    }

    pub const fn significant_figures(self) -> f64 {
        self.significant_figures
    }

    pub const fn complementarity(self) -> f64 {
        self.complementarity
    }

    pub const fn residual_limit(self) -> f64 {
        self.residual_limit
    }

    pub const fn feasibility_limit(self) -> f64 {
        self.feasibility_limit
    }

    pub const fn is_finite(self) -> bool {
        self.finite
    }

    pub const fn accepted(self) -> bool {
        self.accepted
    }
}

/// Successful restricted-range weights and complete solver evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct LoqoSolution {
    validation: LoqoValidation,
    weights: DenseVector,
    objective: f64,
    residual: LoqoResidualEvidence,
    stop_reason: LoqoStopReason,
    trace: Vec<LoqoIterationEvidence>,
    kkt_solves: Vec<LoqoKktSolveEvidence>,
}

impl LoqoSolution {
    pub const fn attempted(&self) -> bool {
        true
    }

    pub const fn validation(&self) -> LoqoValidation {
        self.validation
    }

    pub const fn weights(&self) -> &DenseVector {
        &self.weights
    }

    pub const fn objective(&self) -> f64 {
        self.objective
    }

    pub const fn residual(&self) -> LoqoResidualEvidence {
        self.residual
    }

    pub const fn stop_reason(&self) -> LoqoStopReason {
        self.stop_reason
    }

    pub fn trace(&self) -> &[LoqoIterationEvidence] {
        &self.trace
    }

    pub fn kkt_solves(&self) -> &[LoqoKktSolveEvidence] {
        &self.kkt_solves
    }
}

/// Stable failure classification for the restricted-range path.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum LoqoSolveErrorKind {
    EmptySystem,
    NonSquareInterpolation,
    ConstraintRowMismatch,
    ConstraintColumnMismatch,
    LowerValueMismatch,
    RangeValueMismatch,
    NonFiniteInput,
    KktSolveFailure,
    NonFiniteIterate,
    DualObjectiveAbovePrimal,
    SignificantFiguresDecreased,
    IterationLimit,
    InfeasibleSolution,
    ResidualTooLarge,
}

impl fmt::Display for LoqoSolveErrorKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::EmptySystem => "the restricted-range program has no variables",
            Self::NonSquareInterpolation => "the interpolation matrix is not square",
            Self::ConstraintRowMismatch => "frozen LOQO requires one bounded row per variable",
            Self::ConstraintColumnMismatch => "the bounded matrix has the wrong column count",
            Self::LowerValueMismatch => "the lower-bound count does not match bounded rows",
            Self::RangeValueMismatch => "the range count does not match bounded rows",
            Self::NonFiniteInput => "the restricted-range program contains a non-finite input",
            Self::KktSolveFailure => "a LOQO KKT solve failed",
            Self::NonFiniteIterate => "LOQO produced a non-finite iterate",
            Self::DualObjectiveAbovePrimal => {
                "the frozen LOQO dual objective exceeded its primal objective"
            }
            Self::SignificantFiguresDecreased => {
                "the frozen LOQO significant-figures measure decreased"
            }
            Self::IterationLimit => "LOQO reached its safety iteration limit",
            Self::InfeasibleSolution => "the terminal LOQO candidate violates its bounds",
            Self::ResidualTooLarge => "the terminal LOQO candidate failed residual checks",
        })
    }
}

/// Failure with source stop, KKT, trace, and terminal-candidate evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct LoqoSolveError {
    kind: LoqoSolveErrorKind,
    attempted: bool,
    validation: LoqoValidation,
    trace: Vec<LoqoIterationEvidence>,
    kkt_solves: Vec<LoqoKktSolveEvidence>,
    kkt_failure: Option<Box<LoqoKktFailure>>,
    candidate_weights: Option<DenseVector>,
    residual: Option<Box<LoqoResidualEvidence>>,
}

impl LoqoSolveError {
    pub const fn kind(&self) -> LoqoSolveErrorKind {
        self.kind
    }

    pub const fn attempted(&self) -> bool {
        self.attempted
    }

    pub const fn validation(&self) -> LoqoValidation {
        self.validation
    }

    pub fn trace(&self) -> &[LoqoIterationEvidence] {
        &self.trace
    }

    pub fn kkt_solves(&self) -> &[LoqoKktSolveEvidence] {
        &self.kkt_solves
    }

    pub fn kkt_failure(&self) -> Option<&LoqoKktFailure> {
        self.kkt_failure.as_deref()
    }

    pub const fn candidate_weights(&self) -> Option<&DenseVector> {
        self.candidate_weights.as_ref()
    }

    pub const fn residual(&self) -> Option<LoqoResidualEvidence> {
        match self.residual.as_ref() {
            Some(residual) => Some(**residual),
            None => None,
        }
    }

    pub const fn surfe_error(&self) -> Error {
        Error::LoqoSolverFailure
    }
}

impl fmt::Display for LoqoSolveError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.kind.fmt(formatter)
    }
}

impl std::error::Error for LoqoSolveError {}

/// Safe direct-use errors around frozen `_find_step`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LoqoStepError {
    EmptyVectors,
    DimensionMismatch,
    NonFiniteRatio,
}

impl fmt::Display for LoqoStepError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::EmptyVectors => "LOQO step vectors are empty",
            Self::DimensionMismatch => "LOQO step vectors have different lengths",
            Self::NonFiniteRatio => "LOQO step ratio is non-finite",
        })
    }
}

impl std::error::Error for LoqoStepError {}

/// Inspect inputs without using the frozen LLT result as a solve gate.
pub fn validate_loqo_qp(
    interpolation: &DenseMatrix,
    constraints: &DenseMatrix,
    lower: &DenseVector,
    range: &DenseVector,
) -> LoqoValidation {
    let n = interpolation.rows();
    let hessian = scaled_matrix(interpolation, 2.0);
    let hessian_finite = all_finite(hessian.data());
    let all_inputs_finite = hessian_finite
        && all_finite(constraints.data())
        && all_finite(lower.values())
        && all_finite(range.values());
    LoqoValidation {
        variables_non_empty: n != 0,
        interpolation_square: interpolation.cols() == n,
        constraint_rows_match: constraints.rows() == n,
        constraint_columns_match: constraints.cols() == n,
        lower_values_match: lower.len() == n,
        range_values_match: range.len() == n,
        hessian_finite,
        all_inputs_finite,
        frozen_llt_positive_definite: frozen_llt_positive_definite(&hessian),
    }
}

/// Frozen `_find_step`: the returned divisor is at least one and applies the
/// source's `0.95` fraction-to-boundary rule.
pub fn loqo_step_divisor(direction: &[f64], value: &[f64]) -> Result<f64, LoqoStepError> {
    if direction.is_empty() {
        return Err(LoqoStepError::EmptyVectors);
    }
    if direction.len() != value.len() {
        return Err(LoqoStepError::DimensionMismatch);
    }
    let divisor = step_divisor_unchecked(direction, value);
    if divisor.is_finite() {
        Ok(divisor)
    } else {
        Err(LoqoStepError::NonFiniteRatio)
    }
}

/// Solve with the default safety-only iteration cap.
pub fn solve_loqo_qp(
    interpolation: &DenseMatrix,
    constraints: &DenseMatrix,
    lower: &DenseVector,
    range: &DenseVector,
) -> Result<LoqoSolution, LoqoSolveError> {
    solve_loqo_qp_with_options(
        interpolation,
        constraints,
        lower,
        range,
        LoqoOptions::default(),
    )
}

/// Solve `min x' interpolation x` subject to `lower <= A*x <= lower+range`.
#[allow(clippy::too_many_lines)]
pub fn solve_loqo_qp_with_options(
    interpolation: &DenseMatrix,
    constraints: &DenseMatrix,
    lower: &DenseVector,
    range: &DenseVector,
    options: LoqoOptions,
) -> Result<LoqoSolution, LoqoSolveError> {
    let validation = validate_loqo_qp(interpolation, constraints, lower, range);
    if !validation.variables_non_empty {
        return Err(early_error(LoqoSolveErrorKind::EmptySystem, validation));
    }
    if !validation.interpolation_square {
        return Err(early_error(
            LoqoSolveErrorKind::NonSquareInterpolation,
            validation,
        ));
    }
    if !validation.constraint_rows_match {
        return Err(early_error(
            LoqoSolveErrorKind::ConstraintRowMismatch,
            validation,
        ));
    }
    if !validation.constraint_columns_match {
        return Err(early_error(
            LoqoSolveErrorKind::ConstraintColumnMismatch,
            validation,
        ));
    }
    if !validation.lower_values_match {
        return Err(early_error(
            LoqoSolveErrorKind::LowerValueMismatch,
            validation,
        ));
    }
    if !validation.range_values_match {
        return Err(early_error(
            LoqoSolveErrorKind::RangeValueMismatch,
            validation,
        ));
    }

    let n = interpolation.rows();
    let hessian = scaled_matrix(interpolation, 2.0);
    let mut trace = Vec::new();
    let mut kkt_solves = Vec::new();
    let initial_kkt = initial_kkt(&hessian, constraints);
    let mut initial_rhs = vec![0.0; 2 * n];
    initial_rhs[n..].copy_from_slice(lower.values());
    let initial = match solve_kkt_stage(
        &initial_kkt,
        &DenseVector::from_values(initial_rhs),
        LoqoKktStage::Initial,
        &mut kkt_solves,
    ) {
        Ok(solution) => solution,
        Err(failure) => match failure.source().candidate_weights() {
            Some(candidate) if candidate.len() == 2 * n => candidate.values().to_vec(),
            _ => {
                let kind = if validation.all_inputs_finite {
                    LoqoSolveErrorKind::KktSolveFailure
                } else {
                    LoqoSolveErrorKind::NonFiniteInput
                };
                return Err(attempted_error(
                    kind,
                    validation,
                    trace,
                    kkt_solves,
                    Some(failure),
                    None,
                    None,
                ));
            }
        },
    };

    let mut state = LoqoState::new(initial[..n].to_vec(), initial[n..].to_vec(), range.values());
    // The source reads uninitialized Eigen residual vectors before iteration
    // one. They only affect printed diagnostics. Deterministic zero evidence
    // matches the frozen probe without copying undefined memory behavior.
    let mut lagged = LoqoResidualVectors::zeros(n);
    let mut last_significant_figures = 0.0;

    for iteration in 1usize.. {
        let primal_objective = source_primal_objective(&hessian, &state.x);
        let dual_objective =
            dot(lower.values(), &state.y) - primal_objective - dot(range.values(), &state.q);
        let gap_ratio = (primal_objective - dual_objective).abs() / (primal_objective.abs() + 1.0);
        let significant_figures = cxx_max(-gap_ratio.log10(), 0.0);
        let (primal_infeasibility, dual_infeasibility) = infeasibility(&lagged, lower.values());
        trace.push(LoqoIterationEvidence {
            iteration,
            primal_objective,
            dual_objective,
            significant_figures,
            primal_infeasibility,
            dual_infeasibility,
            predictor_primal_divisor: None,
            predictor_dual_divisor: None,
            predictor_fraction: None,
            predictor_mu: None,
            corrector_primal_divisor: None,
            corrector_dual_divisor: None,
        });

        if significant_figures > SIGNIFICANT_FIGURES_TARGET {
            return finish_candidate(
                validation,
                interpolation,
                constraints,
                lower,
                range,
                state,
                trace,
                kkt_solves,
            );
        }
        if dual_objective > primal_objective {
            let residual = final_residual_evidence(&hessian, constraints, lower, range, &state);
            return Err(attempted_error(
                LoqoSolveErrorKind::DualObjectiveAbovePrimal,
                validation,
                trace,
                kkt_solves,
                None,
                Some(DenseVector::from_values(state.x)),
                Some(residual),
            ));
        }
        if significant_figures < last_significant_figures {
            let residual = final_residual_evidence(&hessian, constraints, lower, range, &state);
            return Err(attempted_error(
                LoqoSolveErrorKind::SignificantFiguresDecreased,
                validation,
                trace,
                kkt_solves,
                None,
                Some(DenseVector::from_values(state.x)),
                Some(residual),
            ));
        }
        if iteration > options.max_iterations {
            let residual = final_residual_evidence(&hessian, constraints, lower, range, &state);
            return Err(attempted_error(
                LoqoSolveErrorKind::IterationLimit,
                validation,
                trace,
                kkt_solves,
                None,
                Some(DenseVector::from_values(state.x)),
                Some(residual),
            ));
        }
        last_significant_figures = significant_figures;

        lagged = source_residuals(
            &hessian,
            constraints,
            lower.values(),
            range.values(),
            &state,
        );
        let diagonal = diagonal_blocks(&state);
        let predictor_rhs = reduced_rhs(&lagged, &state, &diagonal);
        let kkt = reduced_kkt(&hessian, constraints, &diagonal);
        let predictor = match solve_kkt_stage(
            &kkt,
            &DenseVector::from_values(predictor_rhs),
            LoqoKktStage::Predictor(iteration),
            &mut kkt_solves,
        ) {
            Ok(solution) => solution,
            Err(failure) => {
                let kind = if validation.all_inputs_finite {
                    LoqoSolveErrorKind::KktSolveFailure
                } else {
                    LoqoSolveErrorKind::NonFiniteInput
                };
                return Err(attempted_error(
                    kind,
                    validation,
                    trace,
                    kkt_solves,
                    Some(failure),
                    Some(DenseVector::from_values(state.x)),
                    None,
                ));
            }
        };
        let predictor_direction =
            directions(&state, &lagged, &diagonal, &predictor[..n], &predictor[n..]);
        let predictor_primal_divisor = positivity_divisor(
            &predictor_direction.dg,
            &state.g,
            &predictor_direction.dw,
            &state.w,
            &predictor_direction.dt,
            &state.t,
            &predictor_direction.dp,
            &state.p,
        );
        let predictor_dual_divisor = positivity_divisor(
            &predictor_direction.dz,
            &state.z,
            &predictor_direction.dv,
            &state.v,
            &predictor_direction.ds,
            &state.s,
            &predictor_direction.dq,
            &state.q,
        );
        let alpha = cxx_max(predictor_primal_divisor, predictor_dual_divisor);
        let fraction_base = (alpha - 1.0) / (alpha + 10.0);
        let fraction = fraction_base.powf(2.0);
        let predictor_mu = complementarity_sum(&state) * fraction / (4 * n) as f64;

        let corrected = corrected_residuals(&lagged, &state, &predictor_direction, predictor_mu);
        let corrected_rhs = reduced_rhs(&corrected, &state, &diagonal);
        let corrector = match solve_kkt_stage(
            &kkt,
            &DenseVector::from_values(corrected_rhs),
            LoqoKktStage::Corrector(iteration),
            &mut kkt_solves,
        ) {
            Ok(solution) => solution,
            Err(failure) => {
                return Err(attempted_error(
                    LoqoSolveErrorKind::KktSolveFailure,
                    validation,
                    trace,
                    kkt_solves,
                    Some(failure),
                    Some(DenseVector::from_values(state.x)),
                    None,
                ));
            }
        };
        let corrector_direction = directions(
            &state,
            &corrected,
            &diagonal,
            &corrector[..n],
            &corrector[n..],
        );
        let corrector_primal_divisor = positivity_divisor(
            &corrector_direction.dg,
            &state.g,
            &corrector_direction.dw,
            &state.w,
            &corrector_direction.dt,
            &state.t,
            &corrector_direction.dp,
            &state.p,
        );
        let corrector_dual_divisor = positivity_divisor(
            &corrector_direction.dz,
            &state.z,
            &corrector_direction.dv,
            &state.v,
            &corrector_direction.ds,
            &state.s,
            &corrector_direction.dq,
            &state.q,
        );
        let evidence = trace.last_mut().expect("current iteration was recorded");
        evidence.predictor_primal_divisor = Some(predictor_primal_divisor);
        evidence.predictor_dual_divisor = Some(predictor_dual_divisor);
        evidence.predictor_fraction = Some(fraction);
        evidence.predictor_mu = Some(predictor_mu);
        evidence.corrector_primal_divisor = Some(corrector_primal_divisor);
        evidence.corrector_dual_divisor = Some(corrector_dual_divisor);

        state.apply_direction(
            &corrector_direction,
            1.0 / corrector_primal_divisor,
            1.0 / corrector_dual_divisor,
        );
        if !state.is_finite() {
            return Err(attempted_error(
                LoqoSolveErrorKind::NonFiniteIterate,
                validation,
                trace,
                kkt_solves,
                None,
                Some(DenseVector::from_values(state.x)),
                None,
            ));
        }
    }
    unreachable!("the unbounded source loop returns through an explicit branch")
}

#[derive(Clone, Debug)]
struct LoqoState {
    x: Vec<f64>,
    y: Vec<f64>,
    g: Vec<f64>,
    z: Vec<f64>,
    t: Vec<f64>,
    s: Vec<f64>,
    v: Vec<f64>,
    w: Vec<f64>,
    p: Vec<f64>,
    q: Vec<f64>,
}

impl LoqoState {
    fn new(x: Vec<f64>, y: Vec<f64>, range: &[f64]) -> Self {
        let n = x.len();
        let mut g = vec![0.0; n];
        let mut z = vec![0.0; n];
        let mut t = vec![0.0; n];
        let mut s = vec![0.0; n];
        let mut v = vec![0.0; n];
        let mut w = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut q = vec![0.0; n];
        for index in 0..n {
            g[index] = cxx_max(x[index].abs(), 100.0);
            z[index] = g[index];
            t[index] = g[index];
            s[index] = g[index];
            v[index] = cxx_max(y[index].abs(), 100.0);
            w[index] = v[index];
            p[index] = cxx_max((range[index] - w[index]).abs(), 100.0);
            q[index] = v[index];
        }
        Self {
            x,
            y,
            g,
            z,
            t,
            s,
            v,
            w,
            p,
            q,
        }
    }

    fn is_finite(&self) -> bool {
        all_finite(&self.x)
            && all_finite(&self.y)
            && all_finite(&self.g)
            && all_finite(&self.z)
            && all_finite(&self.t)
            && all_finite(&self.s)
            && all_finite(&self.v)
            && all_finite(&self.w)
            && all_finite(&self.p)
            && all_finite(&self.q)
    }

    fn apply_direction(&mut self, direction: &LoqoDirection, primal: f64, dual: f64) {
        add_scaled(&mut self.x, &direction.dx, primal);
        add_scaled(&mut self.g, &direction.dg, primal);
        add_scaled(&mut self.w, &direction.dw, primal);
        add_scaled(&mut self.t, &direction.dt, primal);
        add_scaled(&mut self.p, &direction.dp, primal);
        add_scaled(&mut self.y, &direction.dy, dual);
        add_scaled(&mut self.z, &direction.dz, dual);
        add_scaled(&mut self.v, &direction.dv, dual);
        add_scaled(&mut self.s, &direction.ds, dual);
        add_scaled(&mut self.q, &direction.dq, dual);
    }
}

#[derive(Clone, Debug)]
struct LoqoResidualVectors {
    rho: Vec<f64>,
    nu: Vec<f64>,
    alpha: Vec<f64>,
    sigma: Vec<f64>,
    tau: Vec<f64>,
    beta: Vec<f64>,
    gamma_z: Vec<f64>,
    gamma_w: Vec<f64>,
    gamma_s: Vec<f64>,
    gamma_q: Vec<f64>,
}

impl LoqoResidualVectors {
    fn zeros(n: usize) -> Self {
        Self {
            rho: vec![0.0; n],
            nu: vec![0.0; n],
            alpha: vec![0.0; n],
            sigma: vec![0.0; n],
            tau: vec![0.0; n],
            beta: vec![0.0; n],
            gamma_z: vec![0.0; n],
            gamma_w: vec![0.0; n],
            gamma_s: vec![0.0; n],
            gamma_q: vec![0.0; n],
        }
    }
}

#[derive(Clone, Debug)]
struct DiagonalBlocks {
    d: Vec<f64>,
    e: Vec<f64>,
}

#[derive(Clone, Debug)]
struct LoqoDirection {
    dx: Vec<f64>,
    dy: Vec<f64>,
    dg: Vec<f64>,
    dz: Vec<f64>,
    dt: Vec<f64>,
    ds: Vec<f64>,
    dv: Vec<f64>,
    dw: Vec<f64>,
    dp: Vec<f64>,
    dq: Vec<f64>,
}

fn initial_kkt(hessian: &DenseMatrix, constraints: &DenseMatrix) -> DenseMatrix {
    let n = hessian.rows();
    let mut kkt = DenseMatrix::zeros(2 * n, 2 * n);
    for row in 0..n {
        for column in 0..n {
            let identity = if row == column { 1.0 } else { 0.0 };
            kkt.set(row, column, -(value(hessian, row, column) + identity));
            kkt.set(row, n + column, value(constraints, column, row));
            kkt.set(n + row, column, value(constraints, row, column));
            kkt.set(n + row, n + column, if row == column { 1.0 } else { 0.0 });
        }
    }
    kkt
}

fn diagonal_blocks(state: &LoqoState) -> DiagonalBlocks {
    let mut d = Vec::with_capacity(state.x.len());
    let mut e = Vec::with_capacity(state.x.len());
    for index in 0..state.x.len() {
        let s_inverse_t = (1.0 / state.s[index]) * state.t[index];
        let g_z_inverse = state.g[index] * (1.0 / state.z[index]);
        d.push(1.0 / (s_inverse_t + g_z_inverse));
        let v_w_inverse = state.v[index] * (1.0 / state.w[index]);
        let p_inverse_q = (1.0 / state.p[index]) * state.q[index];
        e.push(1.0 / (v_w_inverse + p_inverse_q));
    }
    DiagonalBlocks { d, e }
}

fn source_residuals(
    hessian: &DenseMatrix,
    constraints: &DenseMatrix,
    lower: &[f64],
    range: &[f64],
    state: &LoqoState,
) -> LoqoResidualVectors {
    let ax = matrix_vector(constraints, &state.x);
    let aty = transpose_matrix_vector(constraints, &state.y);
    let hx = matrix_vector(hessian, &state.x);
    let mut residuals = LoqoResidualVectors::zeros(state.x.len());
    for index in 0..state.x.len() {
        residuals.rho[index] = lower[index] - ax[index] + state.w[index];
        residuals.nu[index] = -state.x[index] + state.g[index] - state.t[index];
        residuals.alpha[index] = range[index] - state.w[index] - state.p[index];
        residuals.sigma[index] = -aty[index] - state.z[index] + hx[index];
        residuals.tau[index] = -state.z[index] - state.s[index];
        residuals.beta[index] = state.y[index] + state.q[index] - state.v[index];
        residuals.gamma_z[index] = -state.z[index];
        residuals.gamma_w[index] = -state.w[index];
        residuals.gamma_s[index] = -state.s[index];
        residuals.gamma_q[index] = -state.q[index];
    }
    residuals
}

fn corrected_residuals(
    residuals: &LoqoResidualVectors,
    state: &LoqoState,
    direction: &LoqoDirection,
    mu: f64,
) -> LoqoResidualVectors {
    let mut corrected = residuals.clone();
    for index in 0..state.x.len() {
        corrected.gamma_z[index] = mu * (1.0 / state.g[index])
            - state.z[index]
            - (1.0 / state.g[index]) * direction.dg[index] * direction.dz[index];
        corrected.gamma_w[index] = mu * (1.0 / state.v[index])
            - state.w[index]
            - (1.0 / state.v[index]) * direction.dv[index] * direction.dw[index];
        corrected.gamma_s[index] = mu * (1.0 / state.t[index])
            - state.s[index]
            - (1.0 / state.t[index]) * direction.dt[index] * direction.ds[index];
        corrected.gamma_q[index] = mu * (1.0 / state.p[index])
            - state.q[index]
            - (1.0 / state.p[index]) * direction.dp[index] * direction.dq[index];
    }
    corrected
}

fn reduced_rhs(
    residuals: &LoqoResidualVectors,
    state: &LoqoState,
    diagonal: &DiagonalBlocks,
) -> Vec<f64> {
    let n = state.x.len();
    let mut rhs = vec![0.0; 2 * n];
    for index in 0..n {
        let tau_hat = residuals.tau[index] - residuals.gamma_s[index];
        let v_w_inverse = state.v[index] * (1.0 / state.w[index]);
        let beta_hat = residuals.beta[index] - v_w_inverse * residuals.gamma_w[index];
        let p_inverse_q = (1.0 / state.p[index]) * state.q[index];
        let p_q_inverse = state.p[index] * (1.0 / state.q[index]);
        let alpha_hat = residuals.alpha[index] - p_q_inverse * residuals.gamma_q[index];
        let g_z_inverse = state.g[index] * (1.0 / state.z[index]);
        let nu_hat = residuals.nu[index] + g_z_inverse * residuals.gamma_z[index];
        let s_inverse_t = (1.0 / state.s[index]) * state.t[index];
        rhs[index] = residuals.sigma[index] - diagonal.d[index] * (nu_hat + s_inverse_t * tau_hat);
        rhs[n + index] =
            residuals.rho[index] - diagonal.e[index] * (beta_hat - p_inverse_q * alpha_hat);
    }
    rhs
}

fn reduced_kkt(
    hessian: &DenseMatrix,
    constraints: &DenseMatrix,
    diagonal: &DiagonalBlocks,
) -> DenseMatrix {
    let n = hessian.rows();
    let mut kkt = DenseMatrix::zeros(2 * n, 2 * n);
    for row in 0..n {
        for column in 0..n {
            let d = if row == column { diagonal.d[row] } else { 0.0 };
            kkt.set(row, column, -(value(hessian, row, column) + d));
            kkt.set(row, n + column, value(constraints, column, row));
            kkt.set(n + row, column, value(constraints, row, column));
            kkt.set(
                n + row,
                n + column,
                if row == column { diagonal.e[row] } else { 0.0 },
            );
        }
    }
    kkt
}

fn directions(
    state: &LoqoState,
    residuals: &LoqoResidualVectors,
    diagonal: &DiagonalBlocks,
    dx: &[f64],
    dy: &[f64],
) -> LoqoDirection {
    let n = state.x.len();
    let mut direction = LoqoDirection {
        dx: dx.to_vec(),
        dy: dy.to_vec(),
        dg: vec![0.0; n],
        dz: vec![0.0; n],
        dt: vec![0.0; n],
        ds: vec![0.0; n],
        dv: vec![0.0; n],
        dw: vec![0.0; n],
        dp: vec![0.0; n],
        dq: vec![0.0; n],
    };
    for index in 0..n {
        let tau_hat = residuals.tau[index] - residuals.gamma_s[index];
        let v_w_inverse = state.v[index] * (1.0 / state.w[index]);
        let beta_hat = residuals.beta[index] - v_w_inverse * residuals.gamma_w[index];
        let p_inverse_q = (1.0 / state.p[index]) * state.q[index];
        let p_q_inverse = state.p[index] * (1.0 / state.q[index]);
        let alpha_hat = residuals.alpha[index] - p_q_inverse * residuals.gamma_q[index];
        let g_z_inverse = state.g[index] * (1.0 / state.z[index]);
        let nu_hat = residuals.nu[index] + g_z_inverse * residuals.gamma_z[index];
        direction.dw[index] = -diagonal.e[index] * (beta_hat - p_inverse_q * alpha_hat + dy[index]);
        let d_s_inverse_t = (diagonal.d[index] * (1.0 / state.s[index])) * state.t[index];
        direction.dt[index] = -d_s_inverse_t * (g_z_inverse * tau_hat - nu_hat + dx[index]);
        direction.dz[index] =
            ((1.0 / state.g[index]) * state.z[index]) * (nu_hat - dx[index] - direction.dt[index]);
        direction.dq[index] = p_inverse_q * (direction.dw[index] - alpha_hat);
        direction.dv[index] = v_w_inverse * (residuals.gamma_w[index] - direction.dw[index]);
        direction.ds[index] = residuals.gamma_s[index]
            - (state.s[index] * (1.0 / state.t[index])) * direction.dt[index];
        direction.dp[index] = (state.p[index] * (1.0 / state.q[index]))
            * (residuals.gamma_q[index] - direction.dq[index]);
        direction.dg[index] = g_z_inverse * (residuals.gamma_z[index] - direction.dz[index]);
    }
    direction
}

#[allow(clippy::too_many_arguments)]
fn positivity_divisor(
    da: &[f64],
    a: &[f64],
    db: &[f64],
    b: &[f64],
    dc: &[f64],
    c: &[f64],
    dd: &[f64],
    d: &[f64],
) -> f64 {
    let mut maximum = step_divisor_unchecked(da, a);
    for divisor in [
        step_divisor_unchecked(db, b),
        step_divisor_unchecked(dc, c),
        step_divisor_unchecked(dd, d),
    ] {
        maximum = cxx_max(divisor, maximum);
    }
    maximum
}

fn step_divisor_unchecked(direction: &[f64], value: &[f64]) -> f64 {
    let mut maximum = (direction[0] / value[0]).abs();
    for index in 1..direction.len() {
        let ratio = (direction[index] / value[index]).abs();
        if ratio > maximum {
            maximum = ratio;
        }
    }
    cxx_max(maximum / 0.95, 1.0)
}

fn complementarity_sum(state: &LoqoState) -> f64 {
    dot(&state.z, &state.g)
        + dot(&state.v, &state.w)
        + dot(&state.s, &state.t)
        + dot(&state.p, &state.q)
}

fn infeasibility(residuals: &LoqoResidualVectors, lower: &[f64]) -> (f64, f64) {
    let primal_sum = dot(&residuals.rho, &residuals.rho)
        + dot(&residuals.tau, &residuals.tau)
        + dot(&residuals.alpha, &residuals.alpha)
        + dot(&residuals.nu, &residuals.nu);
    let primal = primal_sum.sqrt() / (dot(lower, lower).sqrt() + 1.0);
    let dual =
        (dot(&residuals.sigma, &residuals.sigma) + dot(&residuals.beta, &residuals.beta)).sqrt();
    (primal, dual)
}

fn solve_kkt_stage(
    kkt: &DenseMatrix,
    right_hand_side: &DenseVector,
    stage: LoqoKktStage,
    evidence: &mut Vec<LoqoKktSolveEvidence>,
) -> Result<Vec<f64>, LoqoKktFailure> {
    match solve_dense_partial_pivot_lu(kkt, right_hand_side) {
        Ok(solution) => {
            evidence.push(LoqoKktSolveEvidence {
                stage,
                dimension: kkt.rows(),
                row_transpositions: solution.factorization().row_transpositions().to_vec(),
                residual: solution.residual(),
            });
            Ok(solution.weights().values().to_vec())
        }
        Err(source) => {
            // Frozen LOQO checks only whether predictor/corrector weights are
            // finite. A stricter T18 residual verdict remains evidence, but
            // must not become an extra solve gate here.
            if let (Some(candidate), Some(factorization), Some(residual)) = (
                source.candidate_weights(),
                source.factorization(),
                source.residual(),
            ) {
                if candidate.len() == kkt.rows()
                    && candidate.values().iter().all(|value| value.is_finite())
                {
                    evidence.push(LoqoKktSolveEvidence {
                        stage,
                        dimension: kkt.rows(),
                        row_transpositions: factorization.row_transpositions().to_vec(),
                        residual: *residual,
                    });
                    return Ok(candidate.values().to_vec());
                }
            }
            Err(LoqoKktFailure { stage, source })
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn finish_candidate(
    validation: LoqoValidation,
    interpolation: &DenseMatrix,
    constraints: &DenseMatrix,
    lower: &DenseVector,
    range: &DenseVector,
    state: LoqoState,
    trace: Vec<LoqoIterationEvidence>,
    kkt_solves: Vec<LoqoKktSolveEvidence>,
) -> Result<LoqoSolution, LoqoSolveError> {
    let hessian = scaled_matrix(interpolation, 2.0);
    let residual = final_residual_evidence(&hessian, constraints, lower, range, &state);
    let objective = source_quadratic_objective(interpolation, &state.x);
    if residual.accepted {
        return Ok(LoqoSolution {
            validation,
            weights: DenseVector::from_values(state.x),
            objective,
            residual,
            stop_reason: LoqoStopReason::SignificantFigures,
            trace,
            kkt_solves,
        });
    }
    let kind = if !residual.finite {
        LoqoSolveErrorKind::NonFiniteIterate
    } else if residual.minimum_lower_slack < -residual.feasibility_limit
        || residual.minimum_upper_slack < -residual.feasibility_limit
    {
        LoqoSolveErrorKind::InfeasibleSolution
    } else {
        LoqoSolveErrorKind::ResidualTooLarge
    };
    Err(attempted_error(
        kind,
        validation,
        trace,
        kkt_solves,
        None,
        Some(DenseVector::from_values(state.x)),
        Some(residual),
    ))
}

fn final_residual_evidence(
    hessian: &DenseMatrix,
    constraints: &DenseMatrix,
    lower: &DenseVector,
    range: &DenseVector,
    state: &LoqoState,
) -> LoqoResidualEvidence {
    let residuals = source_residuals(hessian, constraints, lower.values(), range.values(), state);
    let (primal_infeasibility, dual_infeasibility) = infeasibility(&residuals, lower.values());
    let ax = matrix_vector(constraints, &state.x);
    let mut minimum_lower_slack = f64::INFINITY;
    let mut minimum_upper_slack = f64::INFINITY;
    for (index, value) in ax.iter().enumerate() {
        minimum_lower_slack = minimum_lower_slack.min(value - lower.values()[index]);
        minimum_upper_slack =
            minimum_upper_slack.min(lower.values()[index] + range.values()[index] - value);
    }
    let primal_objective = source_primal_objective(hessian, &state.x);
    let dual_objective =
        dot(lower.values(), &state.y) - primal_objective - dot(range.values(), &state.q);
    let significant_figures = cxx_max(
        -((primal_objective - dual_objective).abs() / (primal_objective.abs() + 1.0)).log10(),
        0.0,
    );
    let complementarity = complementarity_sum(state) / (4 * state.x.len()) as f64;
    let residual_scale = 1.0
        + primal_objective
            .abs()
            .max(dual_objective.abs())
            .max(linf(lower.values()))
            .max(linf(range.values()));
    let feasibility_scale = 1.0 + linf(lower.values()).max(linf(range.values()));
    let residual_limit = RESIDUAL_ABSOLUTE_TOLERANCE + RESIDUAL_RELATIVE_TOLERANCE * residual_scale;
    let feasibility_limit =
        FEASIBILITY_ABSOLUTE_TOLERANCE + FEASIBILITY_RELATIVE_TOLERANCE * feasibility_scale;
    let finite = state.is_finite()
        && primal_infeasibility.is_finite()
        && dual_infeasibility.is_finite()
        && minimum_lower_slack.is_finite()
        && minimum_upper_slack.is_finite()
        && significant_figures.is_finite()
        && complementarity.is_finite();
    let accepted = finite
        && significant_figures > SIGNIFICANT_FIGURES_TARGET
        && primal_infeasibility <= residual_limit
        && dual_infeasibility <= residual_limit
        && complementarity <= residual_limit
        && minimum_lower_slack >= -feasibility_limit
        && minimum_upper_slack >= -feasibility_limit;
    LoqoResidualEvidence {
        primal_infeasibility,
        dual_infeasibility,
        minimum_lower_slack,
        minimum_upper_slack,
        significant_figures,
        complementarity,
        residual_limit,
        feasibility_limit,
        finite,
        accepted,
    }
}

fn early_error(kind: LoqoSolveErrorKind, validation: LoqoValidation) -> LoqoSolveError {
    LoqoSolveError {
        kind,
        attempted: false,
        validation,
        trace: Vec::new(),
        kkt_solves: Vec::new(),
        kkt_failure: None,
        candidate_weights: None,
        residual: None,
    }
}

#[allow(clippy::too_many_arguments)]
fn attempted_error(
    kind: LoqoSolveErrorKind,
    validation: LoqoValidation,
    trace: Vec<LoqoIterationEvidence>,
    kkt_solves: Vec<LoqoKktSolveEvidence>,
    kkt_failure: Option<LoqoKktFailure>,
    candidate_weights: Option<DenseVector>,
    residual: Option<LoqoResidualEvidence>,
) -> LoqoSolveError {
    LoqoSolveError {
        kind,
        attempted: true,
        validation,
        trace,
        kkt_solves,
        kkt_failure: kkt_failure.map(Box::new),
        candidate_weights,
        residual: residual.map(Box::new),
    }
}

fn frozen_llt_positive_definite(matrix: &DenseMatrix) -> bool {
    if matrix.rows() == 0
        || matrix.rows() != matrix.cols()
        || matrix.data().iter().any(|value| !value.is_finite())
    {
        return false;
    }
    let n = matrix.rows();
    let mut lower = vec![0.0; n * n];
    for row in 0..n {
        for column in 0..=row {
            let mut sum = value(matrix, row, column);
            for index in 0..column {
                sum -= lower[row * n + index] * lower[column * n + index];
            }
            if row == column {
                if sum <= 0.0 || !sum.is_finite() {
                    return false;
                }
                lower[row * n + column] = sum.sqrt();
            } else {
                lower[row * n + column] = sum / lower[column * n + column];
            }
        }
    }
    true
}

fn scaled_matrix(matrix: &DenseMatrix, scale: f64) -> DenseMatrix {
    DenseMatrix::from_row_major(
        matrix.rows(),
        matrix.cols(),
        matrix.data().iter().map(|value| scale * value).collect(),
    )
    .expect("scaling preserves dense storage dimensions")
}

fn source_primal_objective(hessian: &DenseMatrix, weights: &[f64]) -> f64 {
    0.5 * source_quadratic_objective(hessian, weights)
}

fn source_quadratic_objective(matrix: &DenseMatrix, weights: &[f64]) -> f64 {
    let objective = dot(weights, &matrix_vector(matrix, weights));
    // The frozen Eigen expression preserves a negative zero for the observed
    // all-zero Hessian branch when its finite candidate is negative.
    if objective == 0.0
        && matrix.data().iter().all(|value| *value == 0.0)
        && weights.iter().any(|value| value.is_sign_negative())
    {
        -0.0
    } else {
        objective
    }
}

fn matrix_vector(matrix: &DenseMatrix, vector: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; matrix.rows()];
    for (row, output) in result.iter_mut().enumerate() {
        let mut sum = 0.0;
        for (column, vector_value) in vector.iter().enumerate().take(matrix.cols()) {
            sum += value(matrix, row, column) * vector_value;
        }
        *output = sum;
    }
    result
}

fn transpose_matrix_vector(matrix: &DenseMatrix, vector: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; matrix.cols()];
    for (column, output) in result.iter_mut().enumerate() {
        let mut sum = 0.0;
        for (row, vector_value) in vector.iter().enumerate().take(matrix.rows()) {
            sum += value(matrix, row, column) * vector_value;
        }
        *output = sum;
    }
    result
}

fn value(matrix: &DenseMatrix, row: usize, column: usize) -> f64 {
    matrix.get(row, column).expect("validated dense index")
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    let mut result = 0.0;
    for index in 0..left.len() {
        result += left[index] * right[index];
    }
    result
}

fn linf(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |maximum, value| maximum.max(value.abs()))
}

fn add_scaled(values: &mut [f64], direction: &[f64], scale: f64) {
    for index in 0..values.len() {
        values[index] += scale * direction[index];
    }
}

fn all_finite(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

fn cxx_max(first: f64, second: f64) -> f64 {
    if first < second {
        second
    } else {
        first
    }
}
