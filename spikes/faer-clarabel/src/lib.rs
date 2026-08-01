use std::error::Error;
use std::fmt::{Display, Formatter};

use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSettings, DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SolverStatus,
    ZeroConeT,
};
use faer::Side;
use faer::linalg::solvers::Solve;
use faer::prelude::*;

pub const FAER_VERSION: &str = "0.24.4";
pub const CLARABEL_VERSION: &str = "0.11.1";
pub const CAPACITY_LIMIT_BYTES: usize = 8 * 1024 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProbeError {
    pub backend: &'static str,
    pub version: &'static str,
    pub problem_class: &'static str,
    pub reason: AttemptFailureReason,
    pub detail: String,
}

impl Display for ProbeError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "backend={} version={} class={} reason={:?} detail={}",
            self.backend, self.version, self.problem_class, self.reason, self.detail
        )
    }
}

impl Error for ProbeError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapabilityState {
    Proven,
    Unavailable,
    Ambiguous,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveAttemptTermination {
    Unsolved,
    Solved,
    PrimalInfeasible,
    DualInfeasible,
    AlmostSolved,
    AlmostPrimalInfeasible,
    AlmostDualInfeasible,
    MaxIterations,
    MaxTime,
    NumericalError,
    InsufficientProgress,
    CallbackTerminated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttemptFailureReason {
    BackendSetup,
    FactorizationFailure,
    UnexpectedTermination,
    NonFiniteCandidate,
    CertificateVerification,
}

fn faer_probe_error(
    problem_class: &'static str,
    reason: AttemptFailureReason,
    detail: impl Into<String>,
) -> ProbeError {
    ProbeError {
        backend: "faer",
        version: FAER_VERSION,
        problem_class,
        reason,
        detail: detail.into(),
    }
}

fn clarabel_probe_error(
    problem_class: &'static str,
    reason: AttemptFailureReason,
    detail: impl Into<String>,
) -> ProbeError {
    ProbeError {
        backend: "clarabel",
        version: CLARABEL_VERSION,
        problem_class,
        reason,
        detail: detail.into(),
    }
}

impl Display for SolveAttemptTermination {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl From<SolverStatus> for SolveAttemptTermination {
    fn from(status: SolverStatus) -> Self {
        match status {
            SolverStatus::Unsolved => Self::Unsolved,
            SolverStatus::Solved => Self::Solved,
            SolverStatus::PrimalInfeasible => Self::PrimalInfeasible,
            SolverStatus::DualInfeasible => Self::DualInfeasible,
            SolverStatus::AlmostSolved => Self::AlmostSolved,
            SolverStatus::AlmostPrimalInfeasible => Self::AlmostPrimalInfeasible,
            SolverStatus::AlmostDualInfeasible => Self::AlmostDualInfeasible,
            SolverStatus::MaxIterations => Self::MaxIterations,
            SolverStatus::MaxTime => Self::MaxTime,
            SolverStatus::NumericalError => Self::NumericalError,
            SolverStatus::InsufficientProgress => Self::InsufficientProgress,
            SolverStatus::CallbackTerminated => Self::CallbackTerminated,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapacityOutcome {
    Accepted { bytes: usize },
    RejectedBeforeAllocation,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CapacityEvidence {
    pub state: CapabilityState,
    pub first_oversize_square_dimension: usize,
    pub first_oversize_square_bytes: usize,
    pub representable_oversize: CapacityOutcome,
    pub arithmetic_overflow: CapacityOutcome,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThreadEvidence {
    pub requested: usize,
    pub actual: usize,
    pub process_global_state_modified: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InertiaEvidence {
    pub positive: usize,
    pub negative: usize,
    pub zero: usize,
}

#[derive(Debug, Clone)]
pub struct KktEvidence {
    pub solution: Vec<f64>,
    pub normalized_backward_error: f64,
    pub inertia: InertiaEvidence,
    pub adapter_scaling_applied: bool,
}

#[derive(Debug, Clone)]
pub struct FactorizationEvidence {
    pub cholesky_spd_succeeded: bool,
    pub cholesky_indefinite_rejected: bool,
    pub col_pivoted_qr_rank: usize,
    pub svd_rank: usize,
    pub singular_values: Vec<f64>,
    pub capacity: CapacityEvidence,
}

#[derive(Debug, Clone)]
pub struct FaerEvidence {
    pub version: &'static str,
    pub features: Vec<&'static str>,
    pub kkt: KktEvidence,
    pub factorizations: FactorizationEvidence,
    pub threads: ThreadEvidence,
    pub failure_reason: Option<AttemptFailureReason>,
}

#[derive(Debug, Clone)]
pub struct ClarabelSettingsEvidence {
    pub max_iterations: u32,
    pub max_threads: u32,
    pub equilibration: bool,
    pub iterative_refinement: bool,
    pub static_regularization: bool,
    pub dynamic_regularization: bool,
    pub direct_solve_method: String,
    pub feasibility_tolerance: f64,
}

#[derive(Debug, Clone)]
pub struct ClarabelScalingEvidence {
    pub variable: Vec<f64>,
    pub inverse_variable: Vec<f64>,
    pub constraint: Vec<f64>,
    pub inverse_constraint: Vec<f64>,
    pub objective: f64,
}

#[derive(Debug, Clone)]
pub struct ClarabelEvidence {
    pub version: &'static str,
    pub features: Vec<&'static str>,
    pub problem_class: &'static str,
    pub cones: Vec<String>,
    pub primal: Vec<f64>,
    pub dual: Vec<f64>,
    pub slack: Vec<f64>,
    pub termination: SolveAttemptTermination,
    pub primal_residual: f64,
    pub dual_residual: f64,
    pub primal_infeasibility_residual: f64,
    pub dual_infeasibility_residual: f64,
    pub absolute_gap: f64,
    pub relative_gap: f64,
    pub iterations: u32,
    pub settings: ClarabelSettingsEvidence,
    pub scaling: ClarabelScalingEvidence,
    pub threads: ThreadEvidence,
    pub linear_solver: String,
    pub failure_reason: Option<AttemptFailureReason>,
}

#[derive(Debug, Clone)]
pub struct PrimalInfeasibilityEvidence {
    pub attempt: ClarabelEvidence,
    pub certificate: Vec<f64>,
    pub certificate_residual: f64,
    pub cone_violation: f64,
    pub separation_margin: f64,
}

#[derive(Debug, Clone)]
pub struct UnboundednessEvidence {
    pub attempt: ClarabelEvidence,
    pub certificate: Vec<f64>,
    pub certificate_residual: f64,
    pub cone_violation: f64,
    pub descent_margin: f64,
}

pub fn probe_faer() -> Result<FaerEvidence, ProbeError> {
    let parallelism_before = faer::get_global_parallelism().degree();
    let kkt: Mat<f64> = faer::mat![[2.0, 0.0, 1.0], [0.0, 2.0, 1.0], [1.0, 1.0, 0.0],];
    let rhs: Mat<f64> = faer::mat![[2.0], [0.0], [0.0]];
    let factor = kkt.lblt(Side::Lower);
    let solved = factor.solve(&rhs);
    let solution = (0..solved.nrows())
        .map(|row| solved[(row, 0)])
        .collect::<Vec<_>>();

    if solution.iter().any(|value| !value.is_finite()) {
        return Err(faer_probe_error(
            "symmetric_indefinite_kkt",
            AttemptFailureReason::NonFiniteCandidate,
            "non-finite KKT candidate",
        ));
    }

    let diagonal = (0..factor.B_diag().dim())
        .map(|index| factor.B_diag()[index])
        .collect::<Vec<_>>();
    let subdiagonal = (0..factor.B_subdiag().dim())
        .map(|index| factor.B_subdiag()[index])
        .collect::<Vec<_>>();
    let inertia = inertia_from_lblt_blocks(&diagonal, &subdiagonal, 1.0e-12);

    let rank_deficient: Mat<f64> = faer::mat![[1.0, 2.0], [2.0, 4.0], [3.0, 6.0],];
    let qr = rank_deficient.col_piv_qr();
    let qr_diagonal = (0..2)
        .map(|index| qr.R()[(index, index)].abs())
        .collect::<Vec<_>>();
    let qr_threshold = qr_diagonal[0] * 1.0e-12;
    let col_pivoted_qr_rank = qr_diagonal
        .iter()
        .filter(|value| **value > qr_threshold)
        .count();

    let singular_values = rank_deficient.singular_values().map_err(|error| {
        faer_probe_error(
            "rank_svd",
            AttemptFailureReason::FactorizationFailure,
            format!("{error:?}"),
        )
    })?;
    let svd_threshold = singular_values[0] * 1.0e-12;
    let svd_rank = singular_values
        .iter()
        .filter(|value| **value > svd_threshold)
        .count();

    let spd: Mat<f64> = faer::mat![[4.0, 1.0], [1.0, 3.0]];
    let indefinite: Mat<f64> = faer::mat![[1.0, 2.0], [2.0, 1.0]];

    let parallelism_after = faer::get_global_parallelism().degree();
    Ok(FaerEvidence {
        version: FAER_VERSION,
        features: vec!["linalg", "std"],
        kkt: KktEvidence {
            normalized_backward_error: normalized_backward_error(&kkt, &solution, &rhs),
            solution,
            inertia,
            adapter_scaling_applied: false,
        },
        factorizations: FactorizationEvidence {
            cholesky_spd_succeeded: spd.llt(Side::Lower).is_ok(),
            cholesky_indefinite_rejected: indefinite.llt(Side::Lower).is_err(),
            col_pivoted_qr_rank,
            svd_rank,
            singular_values,
            capacity: capacity_evidence(),
        },
        threads: ThreadEvidence {
            requested: 1,
            actual: parallelism_after,
            process_global_state_modified: parallelism_before != parallelism_after,
        },
        failure_reason: None,
    })
}

pub fn checked_dense_capacity(rows: usize, columns: usize) -> CapacityOutcome {
    let bytes = rows
        .checked_mul(columns)
        .and_then(|elements| elements.checked_mul(std::mem::size_of::<f64>()));
    match bytes {
        Some(bytes) if bytes <= CAPACITY_LIMIT_BYTES => CapacityOutcome::Accepted { bytes },
        _ => CapacityOutcome::RejectedBeforeAllocation,
    }
}

fn capacity_evidence() -> CapacityEvidence {
    let first_oversize_square_dimension = 32_769;
    let first_oversize_square_bytes = first_oversize_square_dimension
        * first_oversize_square_dimension
        * std::mem::size_of::<f64>();
    CapacityEvidence {
        state: CapabilityState::Ambiguous,
        first_oversize_square_dimension,
        first_oversize_square_bytes,
        representable_oversize: checked_dense_capacity(
            first_oversize_square_dimension,
            first_oversize_square_dimension,
        ),
        arithmetic_overflow: checked_dense_capacity(usize::MAX, usize::MAX),
    }
}

fn inertia_from_lblt_blocks(
    diagonal: &[f64],
    subdiagonal: &[f64],
    tolerance: f64,
) -> InertiaEvidence {
    let mut eigenvalues = Vec::with_capacity(diagonal.len());
    let mut index = 0;
    while index < diagonal.len() {
        if index + 1 < diagonal.len() && subdiagonal[index] != 0.0 {
            let a = diagonal[index];
            let b = subdiagonal[index];
            let c = diagonal[index + 1];
            let radius = (a - c).hypot(2.0 * b);
            eigenvalues.push(0.5 * (a + c + radius));
            eigenvalues.push(0.5 * (a + c - radius));
            index += 2;
        } else {
            eigenvalues.push(diagonal[index]);
            index += 1;
        }
    }

    InertiaEvidence {
        positive: eigenvalues
            .iter()
            .filter(|value| **value > tolerance)
            .count(),
        negative: eigenvalues
            .iter()
            .filter(|value| **value < -tolerance)
            .count(),
        zero: eigenvalues
            .iter()
            .filter(|value| value.abs() <= tolerance)
            .count(),
    }
}

fn normalized_backward_error(matrix: &Mat<f64>, solution: &[f64], rhs: &Mat<f64>) -> f64 {
    let residual_norm = (0..matrix.nrows())
        .map(|row| {
            let product = (0..matrix.ncols())
                .map(|column| matrix[(row, column)] * solution[column])
                .sum::<f64>();
            (product - rhs[(row, 0)]).abs()
        })
        .fold(0.0_f64, f64::max);
    let matrix_norm = (0..matrix.nrows())
        .map(|row| {
            (0..matrix.ncols())
                .map(|column| matrix[(row, column)].abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let solution_norm = solution
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let rhs_norm = (0..rhs.nrows())
        .map(|row| rhs[(row, 0)].abs())
        .fold(0.0_f64, f64::max);
    residual_norm / (matrix_norm * solution_norm + rhs_norm)
}

fn clarabel_settings() -> DefaultSettings<f64> {
    DefaultSettings {
        verbose: false,
        max_threads: 1,
        direct_solve_method: "qdldl".into(),
        ..DefaultSettings::default()
    }
}

fn settings_evidence(settings: &DefaultSettings<f64>) -> ClarabelSettingsEvidence {
    ClarabelSettingsEvidence {
        max_iterations: settings.max_iter,
        max_threads: settings.max_threads,
        equilibration: settings.equilibrate_enable,
        iterative_refinement: settings.iterative_refinement_enable,
        static_regularization: settings.static_regularization_enable,
        dynamic_regularization: settings.dynamic_regularization_enable,
        direct_solve_method: settings.direct_solve_method.clone(),
        feasibility_tolerance: settings.tol_feas,
    }
}

pub fn probe_clarabel_qp() -> Result<ClarabelEvidence, ProbeError> {
    let parallelism_before = faer::get_global_parallelism().degree();
    let p = CscMatrix::from(&[[6.0, 0.0], [0.0, 4.0]]);
    let q = vec![-1.0, -4.0];
    let a = CscMatrix::from(&[
        [1.0, -2.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [-1.0, 0.0],
        [0.0, -1.0],
    ]);
    let b = vec![0.0, 1.0, 1.0, 1.0, 1.0];
    let cones = [ZeroConeT(1), NonnegativeConeT(4)];
    let settings = clarabel_settings();
    let settings_evidence = settings_evidence(&settings);
    let mut solver = DefaultSolver::new(&p, &q, &a, &b, &cones, settings).map_err(|error| {
        clarabel_probe_error(
            "convex_qp",
            AttemptFailureReason::BackendSetup,
            error.to_string(),
        )
    })?;
    solver.solve();
    clarabel_evidence(
        "convex_qp",
        vec!["zero(1)".into(), "nonnegative(4)".into()],
        solver,
        settings_evidence,
        parallelism_before,
    )
}

pub fn probe_clarabel_socp() -> Result<ClarabelEvidence, ProbeError> {
    let parallelism_before = faer::get_global_parallelism().degree();
    let p = CscMatrix::from(&[[0.0, 0.0], [0.0, 2.0]]);
    let q = vec![0.0, 0.0];
    let a = CscMatrix::from(&[[0.0, 0.0], [-2.0, 0.0], [0.0, -1.0]]);
    let b = vec![1.0, -2.0, -2.0];
    let cones = [SecondOrderConeT(3)];
    let settings = clarabel_settings();
    let settings_evidence = settings_evidence(&settings);
    let mut solver = DefaultSolver::new(&p, &q, &a, &b, &cones, settings).map_err(|error| {
        clarabel_probe_error(
            "convex_socp",
            AttemptFailureReason::BackendSetup,
            error.to_string(),
        )
    })?;
    solver.solve();
    clarabel_evidence(
        "convex_socp",
        vec!["second_order(3)".into()],
        solver,
        settings_evidence,
        parallelism_before,
    )
}

fn clarabel_evidence(
    problem_class: &'static str,
    cones: Vec<String>,
    solver: DefaultSolver<f64>,
    settings: ClarabelSettingsEvidence,
    parallelism_before: usize,
) -> Result<ClarabelEvidence, ProbeError> {
    if solver.solution.x.iter().any(|value| !value.is_finite()) {
        return Err(clarabel_probe_error(
            problem_class,
            AttemptFailureReason::NonFiniteCandidate,
            "non-finite primal candidate",
        ));
    }

    let parallelism_after = faer::get_global_parallelism().degree();
    let scaling = ClarabelScalingEvidence {
        variable: solver.data.equilibration.d.clone(),
        inverse_variable: solver.data.equilibration.dinv.clone(),
        constraint: solver.data.equilibration.e.clone(),
        inverse_constraint: solver.data.equilibration.einv.clone(),
        objective: solver.data.equilibration.c,
    };
    Ok(ClarabelEvidence {
        version: CLARABEL_VERSION,
        features: vec!["serde"],
        problem_class,
        cones,
        primal: solver.solution.x.clone(),
        dual: solver.solution.z.clone(),
        slack: solver.solution.s.clone(),
        termination: solver.solution.status.into(),
        primal_residual: solver.info.res_primal,
        dual_residual: solver.info.res_dual,
        primal_infeasibility_residual: solver.info.res_primal_inf,
        dual_infeasibility_residual: solver.info.res_dual_inf,
        absolute_gap: solver.info.gap_abs,
        relative_gap: solver.info.gap_rel,
        iterations: solver.info.iterations,
        threads: ThreadEvidence {
            requested: settings.max_threads as usize,
            actual: solver.info.linsolver.threads,
            process_global_state_modified: parallelism_before != parallelism_after,
        },
        linear_solver: solver.info.linsolver.name.clone(),
        settings,
        scaling,
        failure_reason: None,
    })
}

pub fn probe_clarabel_primal_infeasible() -> Result<PrimalInfeasibilityEvidence, ProbeError> {
    let parallelism_before = faer::get_global_parallelism().degree();
    let p = CscMatrix::from(&[[1.0]]);
    let q = vec![0.0];
    let a = CscMatrix::from(&[[1.0], [-1.0]]);
    let b = vec![0.0, -1.0];
    let cones = [NonnegativeConeT(2)];
    let settings = clarabel_settings();
    let settings_evidence = settings_evidence(&settings);
    let mut solver = DefaultSolver::new(&p, &q, &a, &b, &cones, settings).map_err(|error| {
        clarabel_probe_error(
            "infeasible_qp",
            AttemptFailureReason::BackendSetup,
            error.to_string(),
        )
    })?;
    solver.solve();
    require_status(
        solver.solution.status,
        SolverStatus::PrimalInfeasible,
        "primal infeasibility",
    )?;
    let attempt = clarabel_evidence(
        "infeasible_qp",
        vec!["nonnegative(2)".into()],
        solver,
        settings_evidence,
        parallelism_before,
    )?;
    let certificate = attempt.dual.clone();
    let certificate_residual = (certificate[0] - certificate[1]).abs();
    let cone_violation = certificate
        .iter()
        .map(|value| (-value).max(0.0))
        .fold(0.0_f64, f64::max);
    let separation_margin = certificate[1];

    Ok(PrimalInfeasibilityEvidence {
        attempt,
        certificate,
        certificate_residual,
        cone_violation,
        separation_margin,
    })
}

pub fn probe_clarabel_unbounded() -> Result<UnboundednessEvidence, ProbeError> {
    let parallelism_before = faer::get_global_parallelism().degree();
    let p = CscMatrix::from(&[[0.0]]);
    let q = vec![-1.0];
    let a = CscMatrix::from(&[[-1.0]]);
    let b = vec![0.0];
    let cones = [NonnegativeConeT(1)];
    let settings = clarabel_settings();
    let settings_evidence = settings_evidence(&settings);
    let mut solver = DefaultSolver::new(&p, &q, &a, &b, &cones, settings).map_err(|error| {
        clarabel_probe_error(
            "unbounded_qp",
            AttemptFailureReason::BackendSetup,
            error.to_string(),
        )
    })?;
    solver.solve();
    require_status(
        solver.solution.status,
        SolverStatus::DualInfeasible,
        "dual infeasibility",
    )?;
    let attempt = clarabel_evidence(
        "unbounded_qp",
        vec!["nonnegative(1)".into()],
        solver,
        settings_evidence,
        parallelism_before,
    )?;
    let certificate = attempt.primal.clone();
    let certificate_residual = (-certificate[0] + attempt.slack[0]).abs();
    let cone_violation = (-attempt.slack[0]).max(0.0);
    let descent_margin = certificate[0];

    Ok(UnboundednessEvidence {
        attempt,
        certificate,
        certificate_residual,
        cone_violation,
        descent_margin,
    })
}

fn require_status(
    actual: SolverStatus,
    expected: SolverStatus,
    probe: &'static str,
) -> Result<(), ProbeError> {
    if actual == expected {
        Ok(())
    } else {
        Err(clarabel_probe_error(
            probe,
            AttemptFailureReason::UnexpectedTermination,
            format!("returned {actual}, expected {expected}"),
        ))
    }
}
