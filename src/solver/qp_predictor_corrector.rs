//! Ordinary primal-dual predictor-corrector QP matching frozen Surfe.
//!
//! Sources at `290dbe0ab344f4258a4935f05cad0f153f0f69a4`:
//! - `surfe_lib/matrix_solver.{h,cpp}::Quadratic_Predictor_Corrector`;
//! - `math_lib/math_methods.{h,cpp}::quadratic_solver`;
//! - `Math_methods::{_find_step_length,max_element_wrt_zero}`.

use std::fmt;

use crate::{ConstraintSystem, DenseMatrix, DenseVector, Error};

use super::{solve_dense_partial_pivot_lu, LuResidualEvidence, LuSolveError};

const COMPLEMENTARITY_TOLERANCE: f64 = 1.0e-8;
const FEASIBILITY_ABSOLUTE_TOLERANCE: f64 = 1.0e-10;
const FEASIBILITY_RELATIVE_TOLERANCE: f64 = 1.0e-8;

/// Safety-only iteration cap; the frozen convergent path normally terminates
/// through its complementarity tests long before this value.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QpOptions {
    pub max_iterations: usize,
}

impl Default for QpOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
        }
    }
}

/// Shape, finite-input, and frozen LLT-validation evidence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QpValidation {
    variables_non_empty: bool,
    interpolation_square: bool,
    equality_columns_match: bool,
    equality_values_match: bool,
    inequality_columns_match: bool,
    inequality_values_match: bool,
    has_inequalities: bool,
    interpolation_finite: bool,
    all_inputs_finite: bool,
    frozen_llt_positive_definite: bool,
}

impl QpValidation {
    pub const fn variables_are_non_empty(self) -> bool {
        self.variables_non_empty
    }

    pub const fn interpolation_is_square(self) -> bool {
        self.interpolation_square
    }

    pub const fn equality_columns_match(self) -> bool {
        self.equality_columns_match
    }

    pub const fn equality_values_match(self) -> bool {
        self.equality_values_match
    }

    pub const fn inequality_columns_match(self) -> bool {
        self.inequality_columns_match
    }

    pub const fn inequality_values_match(self) -> bool {
        self.inequality_values_match
    }

    pub const fn has_inequalities(self) -> bool {
        self.has_inequalities
    }

    pub const fn interpolation_is_finite(self) -> bool {
        self.interpolation_finite
    }

    pub const fn all_inputs_are_finite(self) -> bool {
        self.all_inputs_finite
    }

    pub const fn frozen_llt_is_positive_definite(self) -> bool {
        self.frozen_llt_positive_definite
    }

    /// Observable result of frozen `validate_matrix_systems()`.
    pub const fn surfe_matrix_system_valid(self) -> bool {
        self.interpolation_finite && self.frozen_llt_positive_definite
    }

    /// Safe dimension checks, deliberately excluding the frozen LLT result.
    /// The real frozen `solve()` path does not call that validator.
    pub const fn safe_shape_valid(self) -> bool {
        self.variables_non_empty
            && self.interpolation_square
            && self.equality_columns_match
            && self.equality_values_match
            && self.inequality_columns_match
            && self.inequality_values_match
            && self.has_inequalities
    }
}

/// One frozen KKT factorization stage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QpKktStage {
    Initial,
    Predictor(usize),
    Corrector(usize),
}

/// Compact evidence from one successful pure-Rust partial-pivot KKT solve.
#[derive(Clone, Debug, PartialEq)]
pub struct QpKktSolveEvidence {
    stage: QpKktStage,
    dimension: usize,
    row_transpositions: Vec<usize>,
    residual: LuResidualEvidence,
}

impl QpKktSolveEvidence {
    pub const fn stage(&self) -> QpKktStage {
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
pub struct QpKktFailure {
    stage: QpKktStage,
    source: LuSolveError,
}

impl QpKktFailure {
    pub const fn stage(&self) -> QpKktStage {
        self.stage
    }

    pub const fn source(&self) -> &LuSolveError {
        &self.source
    }
}

/// Per-iteration evidence at the same point where frozen Surfe prints `mu`.
#[derive(Clone, Debug, PartialEq)]
pub struct QpIterationEvidence {
    iteration: usize,
    mu: f64,
    objective: f64,
    stationarity_linf: f64,
    equality_linf: f64,
    inequality_residual_linf: f64,
    minimum_inequality_slack: f64,
    affine_step: Option<f64>,
    affine_mu: Option<f64>,
    centering_sigma: Option<f64>,
    corrector_step: Option<f64>,
}

impl QpIterationEvidence {
    pub const fn iteration(&self) -> usize {
        self.iteration
    }

    pub const fn mu(&self) -> f64 {
        self.mu
    }

    pub const fn objective(&self) -> f64 {
        self.objective
    }

    pub const fn stationarity_linf(&self) -> f64 {
        self.stationarity_linf
    }

    pub const fn equality_linf(&self) -> f64 {
        self.equality_linf
    }

    pub const fn inequality_residual_linf(&self) -> f64 {
        self.inequality_residual_linf
    }

    pub const fn minimum_inequality_slack(&self) -> f64 {
        self.minimum_inequality_slack
    }

    pub const fn affine_step(&self) -> Option<f64> {
        self.affine_step
    }

    pub const fn affine_mu(&self) -> Option<f64> {
        self.affine_mu
    }

    pub const fn centering_sigma(&self) -> Option<f64> {
        self.centering_sigma
    }

    pub const fn corrector_step(&self) -> Option<f64> {
        self.corrector_step
    }
}

/// The two source-level termination branches.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QpStopReason {
    ComplementarityTolerance,
    ComplementarityIncreased,
}

/// Final objective, residual, feasibility, and complementarity evidence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QpResidualEvidence {
    stationarity_linf: f64,
    equality_linf: f64,
    inequality_residual_linf: f64,
    minimum_inequality_slack: f64,
    complementarity: f64,
    stationarity_limit: f64,
    equality_limit: f64,
    inequality_limit: f64,
    complementarity_limit: f64,
    finite: bool,
    accepted: bool,
}

impl QpResidualEvidence {
    pub const fn stationarity_linf(self) -> f64 {
        self.stationarity_linf
    }

    pub const fn equality_linf(self) -> f64 {
        self.equality_linf
    }

    pub const fn inequality_residual_linf(self) -> f64 {
        self.inequality_residual_linf
    }

    pub const fn minimum_inequality_slack(self) -> f64 {
        self.minimum_inequality_slack
    }

    pub const fn complementarity(self) -> f64 {
        self.complementarity
    }

    pub const fn stationarity_limit(self) -> f64 {
        self.stationarity_limit
    }

    pub const fn equality_limit(self) -> f64 {
        self.equality_limit
    }

    pub const fn inequality_limit(self) -> f64 {
        self.inequality_limit
    }

    pub const fn complementarity_limit(self) -> f64 {
        self.complementarity_limit
    }

    pub const fn is_finite(self) -> bool {
        self.finite
    }

    pub const fn accepted(self) -> bool {
        self.accepted
    }
}

/// Successful ordinary-QP weights and all solver evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct QpSolution {
    validation: QpValidation,
    weights: DenseVector,
    dual_equality: DenseVector,
    dual_inequality: DenseVector,
    slack: DenseVector,
    objective: f64,
    residual: QpResidualEvidence,
    stop_reason: QpStopReason,
    trace: Vec<QpIterationEvidence>,
    kkt_solves: Vec<QpKktSolveEvidence>,
}

impl QpSolution {
    pub const fn attempted(&self) -> bool {
        true
    }

    pub const fn validation(&self) -> QpValidation {
        self.validation
    }

    pub const fn weights(&self) -> &DenseVector {
        &self.weights
    }

    pub const fn dual_equality(&self) -> &DenseVector {
        &self.dual_equality
    }

    pub const fn dual_inequality(&self) -> &DenseVector {
        &self.dual_inequality
    }

    pub const fn slack(&self) -> &DenseVector {
        &self.slack
    }

    pub const fn objective(&self) -> f64 {
        self.objective
    }

    pub const fn residual(&self) -> QpResidualEvidence {
        self.residual
    }

    pub const fn stop_reason(&self) -> QpStopReason {
        self.stop_reason
    }

    pub fn trace(&self) -> &[QpIterationEvidence] {
        &self.trace
    }

    pub fn kkt_solves(&self) -> &[QpKktSolveEvidence] {
        &self.kkt_solves
    }
}

/// Stable failure classification for the ordinary predictor-corrector path.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum QpSolveErrorKind {
    EmptySystem,
    NonSquareInterpolation,
    EqualityColumnMismatch,
    EqualityValueMismatch,
    InequalityColumnMismatch,
    InequalityValueMismatch,
    MissingInequalities,
    NonFiniteInput,
    KktSolveFailure,
    NonFiniteIterate,
    IterationLimit,
    InfeasibleSolution,
    ResidualTooLarge,
}

impl fmt::Display for QpSolveErrorKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::EmptySystem => "the quadratic program has no variables",
            Self::NonSquareInterpolation => "the interpolation matrix is not square",
            Self::EqualityColumnMismatch => "the equality matrix has the wrong column count",
            Self::EqualityValueMismatch => "the equality value count does not match its rows",
            Self::InequalityColumnMismatch => "the inequality matrix has the wrong column count",
            Self::InequalityValueMismatch => "the inequality value count does not match its rows",
            Self::MissingInequalities => {
                "frozen ordinary predictor-corrector requires an inequality row"
            }
            Self::NonFiniteInput => "the quadratic program contains a non-finite input",
            Self::KktSolveFailure => "a predictor-corrector KKT solve failed",
            Self::NonFiniteIterate => "predictor-corrector produced a non-finite iterate",
            Self::IterationLimit => "predictor-corrector reached its safety iteration limit",
            Self::InfeasibleSolution => "the terminal QP candidate violates its constraints",
            Self::ResidualTooLarge => "the terminal QP candidate failed residual checks",
        })
    }
}

/// Failure with attempted-state, trace, KKT, and terminal-candidate evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct QpSolveError {
    kind: QpSolveErrorKind,
    attempted: bool,
    validation: QpValidation,
    trace: Vec<QpIterationEvidence>,
    kkt_solves: Vec<QpKktSolveEvidence>,
    kkt_failure: Option<Box<QpKktFailure>>,
    candidate_weights: Option<DenseVector>,
    residual: Option<Box<QpResidualEvidence>>,
    stop_reason: Option<QpStopReason>,
}

impl QpSolveError {
    pub const fn kind(&self) -> QpSolveErrorKind {
        self.kind
    }

    pub const fn attempted(&self) -> bool {
        self.attempted
    }

    pub const fn validation(&self) -> QpValidation {
        self.validation
    }

    pub fn trace(&self) -> &[QpIterationEvidence] {
        &self.trace
    }

    pub fn kkt_solves(&self) -> &[QpKktSolveEvidence] {
        &self.kkt_solves
    }

    pub fn kkt_failure(&self) -> Option<&QpKktFailure> {
        self.kkt_failure.as_deref()
    }

    pub const fn candidate_weights(&self) -> Option<&DenseVector> {
        self.candidate_weights.as_ref()
    }

    pub const fn residual(&self) -> Option<QpResidualEvidence> {
        match self.residual.as_ref() {
            Some(residual) => Some(**residual),
            None => None,
        }
    }

    pub const fn stop_reason(&self) -> Option<QpStopReason> {
        self.stop_reason
    }

    pub const fn surfe_error(&self) -> Error {
        Error::PredictorCorrectorSolverFailure
    }
}

impl fmt::Display for QpSolveError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.kind.fmt(formatter)
    }
}

impl std::error::Error for QpSolveError {}

/// Safe error for direct use of the frozen step-length helper.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QpStepLengthError {
    EmptyVectors,
    DimensionMismatch,
}

impl fmt::Display for QpStepLengthError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::EmptyVectors => "step-length vectors are empty",
            Self::DimensionMismatch => "step-length vectors have different lengths",
        })
    }
}

impl std::error::Error for QpStepLengthError {}

/// Inspect the QP inputs without using the frozen LLT result as a solve gate.
pub fn validate_predictor_corrector_qp(
    interpolation: &DenseMatrix,
    equality: &ConstraintSystem,
    inequality: &ConstraintSystem,
) -> QpValidation {
    let variables = interpolation.rows();
    let interpolation_finite = all_finite(interpolation.data());
    let all_inputs_finite = interpolation_finite
        && all_finite(equality.matrix().data())
        && all_finite(equality.values().values())
        && all_finite(inequality.matrix().data())
        && all_finite(inequality.values().values());
    QpValidation {
        variables_non_empty: variables != 0,
        interpolation_square: interpolation.cols() == variables,
        equality_columns_match: equality.matrix().cols() == variables,
        equality_values_match: equality.values().len() == equality.matrix().rows(),
        inequality_columns_match: inequality.matrix().cols() == variables,
        inequality_values_match: inequality.values().len() == inequality.matrix().rows(),
        has_inequalities: inequality.matrix().rows() != 0,
        interpolation_finite,
        all_inputs_finite,
        frozen_llt_positive_definite: frozen_llt_positive_definite(interpolation),
    }
}

/// Frozen `_find_step_length`, with safe dimension checks around its exact
/// sign, `1e-14`, and unit-cap behavior.
pub fn predictor_corrector_step_length(
    slack: &[f64],
    slack_direction: &[f64],
    dual: &[f64],
    dual_direction: &[f64],
) -> Result<f64, QpStepLengthError> {
    if slack.is_empty() {
        return Err(QpStepLengthError::EmptyVectors);
    }
    if slack.len() != slack_direction.len()
        || slack.len() != dual.len()
        || slack.len() != dual_direction.len()
    {
        return Err(QpStepLengthError::DimensionMismatch);
    }
    Ok(step_length_unchecked(
        slack,
        slack_direction,
        dual,
        dual_direction,
    ))
}

/// Solve with the default safety-only iteration cap.
pub fn solve_predictor_corrector_qp(
    interpolation: &DenseMatrix,
    equality: &ConstraintSystem,
    inequality: &ConstraintSystem,
) -> Result<QpSolution, QpSolveError> {
    solve_predictor_corrector_qp_with_options(
        interpolation,
        equality,
        inequality,
        QpOptions::default(),
    )
}

/// Solve `min x' interpolation x` subject to `A*x=b` and `C*x>=d`.
///
/// This mirrors the frozen wrapper's `H = 2 * interpolation` before entering
/// `Math_methods::quadratic_solver`.
pub fn solve_predictor_corrector_qp_with_options(
    interpolation: &DenseMatrix,
    equality: &ConstraintSystem,
    inequality: &ConstraintSystem,
    options: QpOptions,
) -> Result<QpSolution, QpSolveError> {
    let validation = validate_predictor_corrector_qp(interpolation, equality, inequality);
    if !validation.variables_non_empty {
        return Err(early_error(QpSolveErrorKind::EmptySystem, validation));
    }
    if !validation.interpolation_square {
        return Err(early_error(
            QpSolveErrorKind::NonSquareInterpolation,
            validation,
        ));
    }
    if !validation.equality_columns_match {
        return Err(early_error(
            QpSolveErrorKind::EqualityColumnMismatch,
            validation,
        ));
    }
    if !validation.equality_values_match {
        return Err(early_error(
            QpSolveErrorKind::EqualityValueMismatch,
            validation,
        ));
    }
    if !validation.inequality_columns_match {
        return Err(early_error(
            QpSolveErrorKind::InequalityColumnMismatch,
            validation,
        ));
    }
    if !validation.inequality_values_match {
        return Err(early_error(
            QpSolveErrorKind::InequalityValueMismatch,
            validation,
        ));
    }
    if !validation.has_inequalities {
        return Err(early_error(
            QpSolveErrorKind::MissingInequalities,
            validation,
        ));
    }

    let n = interpolation.rows();
    let na = equality.matrix().rows();
    let nc = inequality.matrix().rows();
    let hessian = scaled_matrix(interpolation, 2.0);

    let mut x = vec![0.0; n];
    let mut y = vec![0.0; na];
    let mut z = vec![0.0; nc];
    let mut slack = vec![0.0; nc];
    let datanorm = max_coefficient(&hessian).sqrt();
    for index in 0..nc {
        z[index] = datanorm;
        slack[index] = datanorm;
    }

    let mut kkt = build_kkt(&hessian, equality.matrix(), inequality.matrix(), &slack, &z);
    let mut trace = Vec::new();
    let mut kkt_solves = Vec::new();

    let initial_residuals = residuals(&hessian, equality, inequality, &x, &y, &z, &slack);
    let initial_direction = match solve_kkt_stage(
        &kkt,
        &initial_residuals.right_hand_side,
        QpKktStage::Initial,
        &mut kkt_solves,
    ) {
        Ok(direction) => direction,
        Err(failure) => {
            let kind = if validation.all_inputs_finite {
                QpSolveErrorKind::KktSolveFailure
            } else {
                QpSolveErrorKind::NonFiniteInput
            };
            return Err(attempted_error(
                kind,
                validation,
                trace,
                kkt_solves,
                Some(failure),
                Some(DenseVector::from_values(x)),
                None,
                None,
            ));
        }
    };

    let mut dx = initial_direction[0..n].to_vec();
    let mut dy = initial_direction[n..n + na].to_vec();
    let mut dz = initial_direction[n + na..n + na + nc].to_vec();
    let mut ds = initial_direction[n + na + nc..].to_vec();
    add_full_step(&mut x, &dx);
    add_full_step(&mut y, &dy);
    add_full_step(&mut z, &dz);
    add_full_step(&mut slack, &ds);

    let mut max_violation_list = Vec::with_capacity(nc);
    for index in 0..nc {
        max_violation_list.push(max_element_wrt_zero(-z[index], -slack[index]));
    }
    max_violation_list.sort_by(|left, right| {
        left.partial_cmp(right)
            .expect("finite KKT result has ordered violations")
    });
    let max_violation = max_violation_list[nc - 1];
    let shift = 1000.0 + 2.0 * max_violation;
    for index in 0..nc {
        z[index] += shift;
        slack[index] += shift;
    }

    let mut previous_mu = None;
    let mut iteration = 0usize;
    loop {
        update_complementarity_blocks(&mut kkt, n, na, nc, &slack, &z);
        let current = residuals(&hessian, equality, inequality, &x, &y, &z, &slack);
        let mu = current.complementarity_sum / nc as f64;
        trace.push(iteration_evidence(
            iteration,
            mu,
            interpolation,
            inequality,
            &x,
            &current,
        ));

        let stop_reason = if iteration > 5 && previous_mu.is_some_and(|value| mu > value) {
            Some(QpStopReason::ComplementarityIncreased)
        } else {
            previous_mu = Some(mu);
            (mu < COMPLEMENTARITY_TOLERANCE).then_some(QpStopReason::ComplementarityTolerance)
        };
        if let Some(stop_reason) = stop_reason {
            return finish_candidate(
                validation,
                interpolation,
                equality,
                inequality,
                x,
                y,
                z,
                slack,
                mu,
                stop_reason,
                trace,
                kkt_solves,
            );
        }
        if iteration >= options.max_iterations {
            let residual =
                final_residual_evidence(&hessian, equality, inequality, &x, &y, &z, &slack, mu);
            return Err(attempted_error(
                QpSolveErrorKind::IterationLimit,
                validation,
                trace,
                kkt_solves,
                None,
                Some(DenseVector::from_values(x)),
                Some(residual),
                None,
            ));
        }

        let predictor = match solve_kkt_stage(
            &kkt,
            &current.right_hand_side,
            QpKktStage::Predictor(iteration),
            &mut kkt_solves,
        ) {
            Ok(direction) => direction,
            Err(failure) => {
                let kind = if validation.all_inputs_finite {
                    QpSolveErrorKind::KktSolveFailure
                } else {
                    QpSolveErrorKind::NonFiniteInput
                };
                return Err(attempted_error(
                    kind,
                    validation,
                    trace,
                    kkt_solves,
                    Some(failure),
                    Some(DenseVector::from_values(x)),
                    None,
                    None,
                ));
            }
        };
        let dz_aff = predictor[n + na..n + na + nc].to_vec();
        let ds_aff = predictor[n + na + nc..].to_vec();
        let affine_step = step_length_unchecked(&slack, &ds_aff, &z, &dz_aff);
        let mut affine_complementarity_sum = 0.0;
        for index in 0..nc {
            affine_complementarity_sum += (z[index] + affine_step * dz_aff[index])
                * (slack[index] + affine_step * ds_aff[index]);
        }
        let affine_mu = affine_complementarity_sum / nc as f64;
        let ratio = affine_mu / mu;
        let sigma = ratio * ratio * ratio;

        let mut corrected_values = current.right_hand_side.values().to_vec();
        for index in 0..nc {
            let corrected =
                current.complementarity[index] - sigma * mu + dz_aff[index] * ds_aff[index];
            corrected_values[n + na + nc + index] = -corrected;
        }
        let corrected_right_hand_side = DenseVector::from_values(corrected_values);
        let corrector = match solve_kkt_stage(
            &kkt,
            &corrected_right_hand_side,
            QpKktStage::Corrector(iteration),
            &mut kkt_solves,
        ) {
            Ok(direction) => direction,
            Err(failure) => {
                return Err(attempted_error(
                    QpSolveErrorKind::KktSolveFailure,
                    validation,
                    trace,
                    kkt_solves,
                    Some(failure),
                    Some(DenseVector::from_values(x)),
                    None,
                    None,
                ));
            }
        };

        // Preserve the frozen extraction loop, including its `j < n` guards
        // for equality and inequality directions.
        for index in 0..n {
            dx[index] = corrector[index];
            if index < na {
                dy[index] = corrector[n + index];
            }
            if index < nc {
                dz[index] = corrector[n + na + index];
                ds[index] = corrector[n + na + nc + index];
            }
        }
        let corrector_step = step_length_unchecked(&slack, &ds, &z, &dz);
        let evidence = trace.last_mut().expect("current iteration was recorded");
        evidence.affine_step = Some(affine_step);
        evidence.affine_mu = Some(affine_mu);
        evidence.centering_sigma = Some(sigma);
        evidence.corrector_step = Some(corrector_step);

        add_scaled_step(&mut x, &dx, corrector_step);
        add_scaled_step(&mut y, &dy, corrector_step);
        add_scaled_step(&mut z, &dz, corrector_step);
        add_scaled_step(&mut slack, &ds, corrector_step);
        if !all_finite(&x) || !all_finite(&y) || !all_finite(&z) || !all_finite(&slack) {
            return Err(attempted_error(
                QpSolveErrorKind::NonFiniteIterate,
                validation,
                trace,
                kkt_solves,
                None,
                Some(DenseVector::from_values(x)),
                None,
                None,
            ));
        }
        iteration += 1;
    }
}

#[derive(Clone, Debug)]
struct ResidualVectors {
    stationarity: Vec<f64>,
    equality: Vec<f64>,
    inequality: Vec<f64>,
    complementarity: Vec<f64>,
    complementarity_sum: f64,
    right_hand_side: DenseVector,
}

fn residuals(
    hessian: &DenseMatrix,
    equality: &ConstraintSystem,
    inequality: &ConstraintSystem,
    x: &[f64],
    y: &[f64],
    z: &[f64],
    slack: &[f64],
) -> ResidualVectors {
    let hx = matrix_vector(hessian, x);
    let aty = transpose_matrix_vector(equality.matrix(), y);
    let ctz = transpose_matrix_vector(inequality.matrix(), z);
    let ax = matrix_vector(equality.matrix(), x);
    let cx = matrix_vector(inequality.matrix(), x);
    let mut stationarity = Vec::with_capacity(x.len());
    for index in 0..x.len() {
        stationarity.push(hx[index] - aty[index] - ctz[index]);
    }
    let mut equality_residual = Vec::with_capacity(y.len());
    for (value, expected) in ax.iter().zip(equality.values().values()) {
        equality_residual.push(value - expected);
    }
    let mut inequality_residual = Vec::with_capacity(z.len());
    let mut complementarity = Vec::with_capacity(z.len());
    let mut complementarity_sum = 0.0;
    for index in 0..z.len() {
        inequality_residual.push(cx[index] - slack[index] - inequality.values().values()[index]);
        let product = slack[index] * z[index];
        complementarity.push(product);
        complementarity_sum += product;
    }
    let mut right_hand_side = Vec::with_capacity(x.len() + y.len() + 2 * z.len());
    right_hand_side.extend(stationarity.iter().map(|value| -*value));
    right_hand_side.extend(equality_residual.iter().map(|value| -*value));
    right_hand_side.extend(inequality_residual.iter().map(|value| -*value));
    right_hand_side.extend(complementarity.iter().map(|value| -*value));
    ResidualVectors {
        stationarity,
        equality: equality_residual,
        inequality: inequality_residual,
        complementarity,
        complementarity_sum,
        right_hand_side: DenseVector::from_values(right_hand_side),
    }
}

fn solve_kkt_stage(
    kkt: &DenseMatrix,
    right_hand_side: &DenseVector,
    stage: QpKktStage,
    evidence: &mut Vec<QpKktSolveEvidence>,
) -> Result<Vec<f64>, QpKktFailure> {
    match solve_dense_partial_pivot_lu(kkt, right_hand_side) {
        Ok(solution) => {
            evidence.push(QpKktSolveEvidence {
                stage,
                dimension: kkt.rows(),
                row_transpositions: solution.factorization().row_transpositions().to_vec(),
                residual: solution.residual(),
            });
            Ok(solution.weights().values().to_vec())
        }
        Err(source) => Err(QpKktFailure { stage, source }),
    }
}

#[allow(clippy::too_many_arguments)]
fn finish_candidate(
    validation: QpValidation,
    interpolation: &DenseMatrix,
    equality: &ConstraintSystem,
    inequality: &ConstraintSystem,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    slack: Vec<f64>,
    mu: f64,
    stop_reason: QpStopReason,
    trace: Vec<QpIterationEvidence>,
    kkt_solves: Vec<QpKktSolveEvidence>,
) -> Result<QpSolution, QpSolveError> {
    let hessian = scaled_matrix(interpolation, 2.0);
    let objective = dot(&x, &matrix_vector(interpolation, &x));
    let residual = final_residual_evidence(&hessian, equality, inequality, &x, &y, &z, &slack, mu);
    if residual.accepted {
        return Ok(QpSolution {
            validation,
            weights: DenseVector::from_values(x),
            dual_equality: DenseVector::from_values(y),
            dual_inequality: DenseVector::from_values(z),
            slack: DenseVector::from_values(slack),
            objective,
            residual,
            stop_reason,
            trace,
            kkt_solves,
        });
    }
    let kind = if !residual.finite {
        QpSolveErrorKind::NonFiniteIterate
    } else if residual.equality_linf > residual.equality_limit
        || residual.minimum_inequality_slack < -residual.inequality_limit
    {
        QpSolveErrorKind::InfeasibleSolution
    } else {
        QpSolveErrorKind::ResidualTooLarge
    };
    Err(attempted_error(
        kind,
        validation,
        trace,
        kkt_solves,
        None,
        Some(DenseVector::from_values(x)),
        Some(residual),
        Some(stop_reason),
    ))
}

#[allow(clippy::too_many_arguments)]
fn final_residual_evidence(
    hessian: &DenseMatrix,
    equality: &ConstraintSystem,
    inequality: &ConstraintSystem,
    x: &[f64],
    y: &[f64],
    z: &[f64],
    slack: &[f64],
    mu: f64,
) -> QpResidualEvidence {
    let current = residuals(hessian, equality, inequality, x, y, z, slack);
    let actual_slack = matrix_vector(inequality.matrix(), x)
        .into_iter()
        .zip(inequality.values().values())
        .map(|(value, lower)| value - lower)
        .collect::<Vec<_>>();
    let stationarity_linf = linf(&current.stationarity);
    let equality_linf = linf(&current.equality);
    let inequality_residual_linf = linf(&current.inequality);
    let minimum_inequality_slack = actual_slack.iter().copied().fold(f64::INFINITY, f64::min);
    let stationarity_scale = 1.0
        + linf(&matrix_vector(hessian, x))
            .max(linf(&transpose_matrix_vector(equality.matrix(), y)))
            .max(linf(&transpose_matrix_vector(inequality.matrix(), z)));
    let equality_scale = 1.0 + linf(equality.values().values());
    let inequality_scale = 1.0 + linf(inequality.values().values());
    let stationarity_limit = layered_limit(stationarity_scale);
    let equality_limit = layered_limit(equality_scale);
    let inequality_limit = layered_limit(inequality_scale);
    let complementarity_limit = layered_limit(1.0);
    let finite = all_finite(x)
        && all_finite(y)
        && all_finite(z)
        && all_finite(slack)
        && stationarity_linf.is_finite()
        && equality_linf.is_finite()
        && inequality_residual_linf.is_finite()
        && minimum_inequality_slack.is_finite()
        && mu.is_finite();
    let accepted = finite
        && stationarity_linf <= stationarity_limit
        && equality_linf <= equality_limit
        && inequality_residual_linf <= inequality_limit
        && minimum_inequality_slack >= -inequality_limit
        && mu.abs() <= complementarity_limit;
    QpResidualEvidence {
        stationarity_linf,
        equality_linf,
        inequality_residual_linf,
        minimum_inequality_slack,
        complementarity: mu,
        stationarity_limit,
        equality_limit,
        inequality_limit,
        complementarity_limit,
        finite,
        accepted,
    }
}

fn iteration_evidence(
    iteration: usize,
    mu: f64,
    interpolation: &DenseMatrix,
    inequality: &ConstraintSystem,
    x: &[f64],
    residuals: &ResidualVectors,
) -> QpIterationEvidence {
    let actual_slack = matrix_vector(inequality.matrix(), x)
        .into_iter()
        .zip(inequality.values().values())
        .map(|(value, lower)| value - lower)
        .collect::<Vec<_>>();
    QpIterationEvidence {
        iteration,
        mu,
        objective: dot(x, &matrix_vector(interpolation, x)),
        stationarity_linf: linf(&residuals.stationarity),
        equality_linf: linf(&residuals.equality),
        inequality_residual_linf: linf(&residuals.inequality),
        minimum_inequality_slack: actual_slack.into_iter().fold(f64::INFINITY, f64::min),
        affine_step: None,
        affine_mu: None,
        centering_sigma: None,
        corrector_step: None,
    }
}

fn build_kkt(
    hessian: &DenseMatrix,
    equality: &DenseMatrix,
    inequality: &DenseMatrix,
    slack: &[f64],
    dual: &[f64],
) -> DenseMatrix {
    let n = hessian.rows();
    let na = equality.rows();
    let nc = inequality.rows();
    let dimension = n + na + 2 * nc;
    let mut kkt = DenseMatrix::zeros(dimension, dimension);
    for row in 0..n {
        for column in 0..n {
            kkt.set(row, column, value(hessian, row, column));
        }
        for column in 0..na {
            kkt.set(row, n + column, -value(equality, column, row));
        }
        for column in 0..nc {
            kkt.set(row, n + na + column, -value(inequality, column, row));
        }
    }
    for row in 0..na {
        for column in 0..n {
            kkt.set(n + row, column, value(equality, row, column));
        }
    }
    for row in 0..nc {
        for column in 0..n {
            kkt.set(n + na + row, column, value(inequality, row, column));
        }
        kkt.set(n + na + row, n + na + nc + row, -1.0);
    }
    update_complementarity_blocks(&mut kkt, n, na, nc, slack, dual);
    kkt
}

fn update_complementarity_blocks(
    kkt: &mut DenseMatrix,
    n: usize,
    na: usize,
    nc: usize,
    slack: &[f64],
    dual: &[f64],
) {
    for row in 0..nc {
        for column in 0..nc {
            kkt.set(
                n + na + nc + row,
                n + na + column,
                if row == column { slack[row] } else { 0.0 },
            );
            kkt.set(
                n + na + nc + row,
                n + na + nc + column,
                if row == column { dual[row] } else { 0.0 },
            );
        }
    }
}

fn step_length_unchecked(a: &[f64], da: &[f64], b: &[f64], db: &[f64]) -> f64 {
    let mut min_alpha_a = 100.0;
    let mut min_alpha_b = 100.0;
    for index in 0..a.len() {
        let alpha_b = if b[index] > 0.0 {
            b[index] / db[index]
        } else {
            -b[index] / db[index]
        };
        let alpha_a = if a[index] > 0.0 {
            a[index] / da[index]
        } else {
            -a[index] / da[index]
        };
        if alpha_b < min_alpha_b && alpha_b > 1.0e-14 {
            min_alpha_b = alpha_b;
        }
        if alpha_a < min_alpha_a && alpha_a > 1.0e-14 {
            min_alpha_a = alpha_a;
        }
    }
    let mut alpha = if min_alpha_b < min_alpha_a {
        min_alpha_b
    } else {
        min_alpha_a
    };
    if alpha > 1.0 || alpha == 0.0 {
        alpha = 1.0;
    }
    alpha
}

fn max_element_wrt_zero(first: f64, second: f64) -> f64 {
    let mut maximum = first;
    if second > maximum {
        maximum = second;
    }
    if 0.0 > maximum {
        maximum = 0.0;
    }
    maximum
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

fn max_coefficient(matrix: &DenseMatrix) -> f64 {
    let mut maximum = matrix.data()[0];
    for value in &matrix.data()[1..] {
        if *value > maximum {
            maximum = *value;
        }
    }
    maximum
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

fn layered_limit(scale: f64) -> f64 {
    FEASIBILITY_ABSOLUTE_TOLERANCE + FEASIBILITY_RELATIVE_TOLERANCE * scale
}

fn add_full_step(values: &mut [f64], direction: &[f64]) {
    for index in 0..values.len() {
        values[index] += direction[index];
    }
}

fn add_scaled_step(values: &mut [f64], direction: &[f64], step: f64) {
    for index in 0..values.len() {
        values[index] += step * direction[index];
    }
}

fn all_finite(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

fn early_error(kind: QpSolveErrorKind, validation: QpValidation) -> QpSolveError {
    QpSolveError {
        kind,
        attempted: false,
        validation,
        trace: Vec::new(),
        kkt_solves: Vec::new(),
        kkt_failure: None,
        candidate_weights: None,
        residual: None,
        stop_reason: None,
    }
}

#[allow(clippy::too_many_arguments)]
fn attempted_error(
    kind: QpSolveErrorKind,
    validation: QpValidation,
    trace: Vec<QpIterationEvidence>,
    kkt_solves: Vec<QpKktSolveEvidence>,
    kkt_failure: Option<QpKktFailure>,
    candidate_weights: Option<DenseVector>,
    residual: Option<QpResidualEvidence>,
    stop_reason: Option<QpStopReason>,
) -> QpSolveError {
    QpSolveError {
        kind,
        attempted: true,
        validation,
        trace,
        kkt_solves,
        kkt_failure: kkt_failure.map(Box::new),
        candidate_weights,
        residual: residual.map(Box::new),
        stop_reason,
    }
}
