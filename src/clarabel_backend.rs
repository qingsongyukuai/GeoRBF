use clarabel::algebra::CscMatrix;
use clarabel::solver::{DefaultSettings, DefaultSolver, IPSolver, SolverStatus, SupportedConeT};

pub(crate) const CRATE_NAME: &str = "clarabel";
pub(crate) const CRATE_VERSION: &str = "0.11.1";
pub(crate) const FEATURES: [&str; 1] = ["serde"];
pub(crate) const DIRECT_SOLVER: &str = "qdldl";
pub(crate) const REQUESTED_THREADS: usize = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ClarabelAttemptProfile {
    Standard,
    Robust,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClarabelSettingsFingerprint {
    pub(crate) profile: ClarabelAttemptProfile,
    pub(crate) max_iterations: u32,
    pub(crate) max_threads: u32,
    pub(crate) direct_solver: &'static str,
    pub(crate) feasibility_tolerance: f64,
    pub(crate) absolute_gap_tolerance: f64,
    pub(crate) relative_gap_tolerance: f64,
    pub(crate) reduced_feasibility_tolerance: f64,
    pub(crate) equilibration_enabled: bool,
    pub(crate) equilibration_rounds: u32,
    pub(crate) static_regularization_enabled: bool,
    pub(crate) dynamic_regularization_enabled: bool,
    pub(crate) iterative_refinement_enabled: bool,
    pub(crate) iterative_refinement_max_iterations: u32,
    pub(crate) presolve_enabled: bool,
    pub(crate) all_settings: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ClarabelBackendFingerprint {
    pub(crate) crate_name: &'static str,
    pub(crate) crate_version: &'static str,
    pub(crate) features: [&'static str; FEATURES.len()],
    pub(crate) direct_solver: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ClarabelCertificateEvidence {
    pub(crate) primal_infeasibility_residual: f64,
    pub(crate) dual_infeasibility_residual: f64,
    pub(crate) kappa_tau_ratio: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClarabelInternalScalingEvidence {
    pub(crate) variable: Vec<f64>,
    pub(crate) inverse_variable: Vec<f64>,
    pub(crate) constraint: Vec<f64>,
    pub(crate) inverse_constraint: Vec<f64>,
    pub(crate) objective: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ClarabelTermination {
    Solved,
    AlmostSolved,
    PrimalInfeasible,
    DualInfeasible,
    AlmostPrimalInfeasible,
    AlmostDualInfeasible,
    IterationLimit,
    TimeLimit,
    NumericalError,
    InsufficientProgress,
    CallbackTerminated,
    Unsolved,
}

impl From<SolverStatus> for ClarabelTermination {
    fn from(status: SolverStatus) -> Self {
        match status {
            SolverStatus::Solved => Self::Solved,
            SolverStatus::AlmostSolved => Self::AlmostSolved,
            SolverStatus::PrimalInfeasible => Self::PrimalInfeasible,
            SolverStatus::DualInfeasible => Self::DualInfeasible,
            SolverStatus::AlmostPrimalInfeasible => Self::AlmostPrimalInfeasible,
            SolverStatus::AlmostDualInfeasible => Self::AlmostDualInfeasible,
            SolverStatus::MaxIterations => Self::IterationLimit,
            SolverStatus::MaxTime => Self::TimeLimit,
            SolverStatus::NumericalError => Self::NumericalError,
            SolverStatus::InsufficientProgress => Self::InsufficientProgress,
            SolverStatus::CallbackTerminated => Self::CallbackTerminated,
            SolverStatus::Unsolved => Self::Unsolved,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClarabelAttemptEvidence {
    pub(crate) sequence: usize,
    pub(crate) profile: ClarabelAttemptProfile,
    pub(crate) backend: ClarabelBackendFingerprint,
    pub(crate) termination: ClarabelTermination,
    pub(crate) settings: ClarabelSettingsFingerprint,
    pub(crate) requested_threads: usize,
    pub(crate) actual_threads: usize,
    pub(crate) iterations: u32,
    pub(crate) reported_primal_residual: f64,
    pub(crate) reported_dual_residual: f64,
    pub(crate) reported_absolute_gap: f64,
    pub(crate) reported_relative_gap: f64,
    pub(crate) certificate: ClarabelCertificateEvidence,
    pub(crate) linear_solver: String,
    pub(crate) internal_scaling: ClarabelInternalScalingEvidence,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClarabelCandidateEnvelope {
    pub(crate) primal: Vec<f64>,
    pub(crate) dual: Vec<f64>,
    pub(crate) slack: Vec<f64>,
    pub(crate) attempt: ClarabelAttemptEvidence,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ClarabelQpInput<'a> {
    pub(crate) variables: usize,
    pub(crate) constraints: usize,
    pub(crate) hessian: &'a [f64],
    pub(crate) linear_objective: &'a [f64],
    pub(crate) constraint_matrix: &'a [f64],
    pub(crate) constraint_rhs: &'a [f64],
    pub(crate) equality_constraints: usize,
    pub(crate) inequality_constraints: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ClarabelAdapterFailureKind {
    InvalidShape,
    Setup,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ClarabelAdapterFailure {
    pub(crate) kind: ClarabelAdapterFailureKind,
    pub(crate) profile: ClarabelAttemptProfile,
    pub(crate) detail: String,
    pub(crate) backend_invoked: bool,
}

pub(crate) fn solve_qp(
    input: ClarabelQpInput<'_>,
    profile: ClarabelAttemptProfile,
    sequence: usize,
) -> Result<ClarabelCandidateEnvelope, ClarabelAdapterFailure> {
    if input.hessian.len() != input.variables.saturating_mul(input.variables)
        || input.linear_objective.len() != input.variables
        || input.constraint_matrix.len() != input.constraints.saturating_mul(input.variables)
        || input.constraint_rhs.len() != input.constraints
        || input
            .equality_constraints
            .saturating_add(input.inequality_constraints)
            != input.constraints
    {
        return Err(ClarabelAdapterFailure {
            kind: ClarabelAdapterFailureKind::InvalidShape,
            profile,
            detail: "solver-independent QP form has inconsistent dimensions".into(),
            backend_invoked: false,
        });
    }

    let hessian = upper_csc(input.variables, input.hessian);
    let constraints = dense_csc(input.constraints, input.variables, input.constraint_matrix);
    let mut cones = Vec::<SupportedConeT<f64>>::new();
    if input.equality_constraints > 0 {
        cones.push(SupportedConeT::ZeroConeT(input.equality_constraints));
    }
    if input.inequality_constraints > 0 {
        cones.push(SupportedConeT::NonnegativeConeT(
            input.inequality_constraints,
        ));
    }
    let settings = resolved_settings(profile);
    let settings_fingerprint = settings_fingerprint(profile, &settings);
    let mut solver = DefaultSolver::new(
        &hessian,
        input.linear_objective,
        &constraints,
        input.constraint_rhs,
        &cones,
        settings,
    )
    .map_err(|error| ClarabelAdapterFailure {
        kind: ClarabelAdapterFailureKind::Setup,
        profile,
        detail: error.to_string(),
        backend_invoked: true,
    })?;
    solver.solve();

    Ok(ClarabelCandidateEnvelope {
        primal: solver.solution.x.clone(),
        dual: solver.solution.z.clone(),
        slack: solver.solution.s.clone(),
        attempt: ClarabelAttemptEvidence {
            sequence,
            profile,
            backend: ClarabelBackendFingerprint {
                crate_name: CRATE_NAME,
                crate_version: CRATE_VERSION,
                features: FEATURES,
                direct_solver: DIRECT_SOLVER,
            },
            termination: solver.solution.status.into(),
            settings: settings_fingerprint,
            requested_threads: REQUESTED_THREADS,
            actual_threads: solver.info.linsolver.threads,
            iterations: solver.info.iterations,
            reported_primal_residual: solver.info.res_primal,
            reported_dual_residual: solver.info.res_dual,
            reported_absolute_gap: solver.info.gap_abs,
            reported_relative_gap: solver.info.gap_rel,
            certificate: ClarabelCertificateEvidence {
                primal_infeasibility_residual: solver.info.res_primal_inf,
                dual_infeasibility_residual: solver.info.res_dual_inf,
                kappa_tau_ratio: solver.info.ktratio,
            },
            linear_solver: solver.info.linsolver.name.clone(),
            internal_scaling: ClarabelInternalScalingEvidence {
                variable: solver.data.equilibration.d.clone(),
                inverse_variable: solver.data.equilibration.dinv.clone(),
                constraint: solver.data.equilibration.e.clone(),
                inverse_constraint: solver.data.equilibration.einv.clone(),
                objective: solver.data.equilibration.c,
            },
        },
    })
}

fn upper_csc(dimension: usize, row_major: &[f64]) -> CscMatrix<f64> {
    let entries = dimension * (dimension + 1) / 2;
    let mut column_pointers = Vec::with_capacity(dimension + 1);
    let mut row_indices = Vec::with_capacity(entries);
    let mut values = Vec::with_capacity(entries);
    column_pointers.push(0);
    for column in 0..dimension {
        for row in 0..=column {
            row_indices.push(row);
            values.push(row_major[row * dimension + column]);
        }
        column_pointers.push(values.len());
    }
    CscMatrix::new(dimension, dimension, column_pointers, row_indices, values)
}

fn dense_csc(rows: usize, columns: usize, row_major: &[f64]) -> CscMatrix<f64> {
    let mut column_pointers = Vec::with_capacity(columns + 1);
    let mut row_indices = Vec::with_capacity(rows * columns);
    let mut values = Vec::with_capacity(rows * columns);
    column_pointers.push(0);
    for column in 0..columns {
        for row in 0..rows {
            row_indices.push(row);
            values.push(row_major[row * columns + column]);
        }
        column_pointers.push(values.len());
    }
    CscMatrix::new(rows, columns, column_pointers, row_indices, values)
}

fn resolved_settings(profile: ClarabelAttemptProfile) -> DefaultSettings<f64> {
    let (max_iter, tolerance, refinement_iterations) = match profile {
        // Backend stopping targets stay one decade inside GeoRBF's verified
        // convex residual envelope so two active sides of a closed interval do
        // not consume the complete public acceptance budget in backend noise.
        ClarabelAttemptProfile::Standard => (200, 1.0e-9, 10),
        ClarabelAttemptProfile::Robust => (400, 1.0e-10, 20),
    };
    DefaultSettings {
        max_iter,
        time_limit: f64::INFINITY,
        verbose: false,
        max_step_fraction: 0.99,
        tol_gap_abs: tolerance,
        tol_gap_rel: tolerance,
        tol_feas: tolerance,
        tol_infeas_abs: 1.0e-8,
        tol_infeas_rel: 1.0e-8,
        tol_ktratio: 1.0e-6,
        reduced_tol_gap_abs: 5.0e-5,
        reduced_tol_gap_rel: 5.0e-5,
        reduced_tol_feas: 1.0e-4,
        reduced_tol_infeas_abs: 5.0e-12,
        reduced_tol_infeas_rel: 5.0e-5,
        reduced_tol_ktratio: 1.0e-4,
        equilibrate_enable: true,
        equilibrate_max_iter: 10,
        equilibrate_min_scaling: 1.0e-4,
        equilibrate_max_scaling: 1.0e4,
        linesearch_backtrack_step: 0.8,
        min_switch_step_length: 0.1,
        min_terminate_step_length: 1.0e-4,
        max_threads: REQUESTED_THREADS as u32,
        direct_kkt_solver: true,
        direct_solve_method: DIRECT_SOLVER.into(),
        static_regularization_enable: true,
        static_regularization_constant: 1.0e-8,
        static_regularization_proportional: f64::EPSILON * f64::EPSILON,
        dynamic_regularization_enable: true,
        dynamic_regularization_eps: 1.0e-13,
        dynamic_regularization_delta: 2.0e-7,
        iterative_refinement_enable: true,
        iterative_refinement_reltol: 1.0e-13,
        iterative_refinement_abstol: 1.0e-12,
        iterative_refinement_max_iter: refinement_iterations,
        iterative_refinement_stop_ratio: 5.0,
        presolve_enable: true,
        input_sparse_dropzeros: false,
    }
}

pub(crate) fn expected_settings(profile: ClarabelAttemptProfile) -> ClarabelSettingsFingerprint {
    let settings = resolved_settings(profile);
    settings_fingerprint(profile, &settings)
}

fn settings_fingerprint(
    profile: ClarabelAttemptProfile,
    settings: &DefaultSettings<f64>,
) -> ClarabelSettingsFingerprint {
    ClarabelSettingsFingerprint {
        profile,
        max_iterations: settings.max_iter,
        max_threads: settings.max_threads,
        direct_solver: DIRECT_SOLVER,
        feasibility_tolerance: settings.tol_feas,
        absolute_gap_tolerance: settings.tol_gap_abs,
        relative_gap_tolerance: settings.tol_gap_rel,
        reduced_feasibility_tolerance: settings.reduced_tol_feas,
        equilibration_enabled: settings.equilibrate_enable,
        equilibration_rounds: settings.equilibrate_max_iter,
        static_regularization_enabled: settings.static_regularization_enable,
        dynamic_regularization_enabled: settings.dynamic_regularization_enable,
        iterative_refinement_enabled: settings.iterative_refinement_enable,
        iterative_refinement_max_iterations: settings.iterative_refinement_max_iter,
        presolve_enabled: settings.presolve_enable,
        all_settings: format!("{settings:?}"),
    }
}
