use crate::capacity::{CapacityExceededEvidence, ConvexQpCapacityPlan, plan_convex_qp_capacity};
use crate::clarabel_backend::{
    self, ClarabelAdapterFailure, ClarabelAttemptEvidence, ClarabelAttemptProfile,
    ClarabelCandidateEnvelope, ClarabelQpInput, ClarabelTermination,
};
use crate::cubic::GlobalAnisotropyMetric;
use crate::cubic_equality::{
    CanonicalEqualityParticipation, CanonicalInequalitySense, CanonicalRelationToleranceEvidence,
    CanonicalSoftResidualBlockKind, CanonicalViolationLoss, CpdEvidence, CubicCanonicalProblem,
    CubicEqualityCore, CubicEqualityFailure, CubicEqualitySolution, CubicRepresentation,
    FunctionalViolationEnvelope, POLYNOMIAL_DIMENSION, PhysicalSideConditionEvidence,
    RecoveredCubicField, RecoveredHardEquality, RecoveredSemanticLatent, RecoveredSoftEquality,
    RecoveredSoftObjective, RepresentationFailure, SemanticLatentCoefficient,
    canonical_characteristic_field_scale, canonical_fitting_uses, canonical_gauge_offset,
    dense_matrix_vector_product, dot_product, relative_slice_error,
};
use crate::functional::{
    DerivedBlockId, DerivedColumnId, DerivedRowId, FunctionalDimension, GroupId, ResidualId,
    UsageProvenance,
};
use crate::numerical::{EQUALITY_KKT_POLICY_V1, NumericalPolicyId};
#[cfg(test)]
use std::cell::RefCell;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CubicFormKind {
    SymmetricKkt,
    ConvexQp,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CubicAlgebraicPlan {
    pub(crate) form: CubicFormKind,
    pub(crate) hard_equalities: usize,
    pub(crate) affine_inequalities: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveredAffineInequality {
    pub(crate) provenance: UsageProvenance,
    pub(crate) dimension: FunctionalDimension,
    pub(crate) sense: CanonicalInequalitySense,
    pub(crate) bound: f64,
    pub(crate) value: f64,
    pub(crate) slack: f64,
    pub(crate) tolerance: f64,
    pub(crate) violation: f64,
    pub(crate) recovered_violation_channel: Option<f64>,
    pub(crate) violation_loss: Option<CanonicalViolationLoss>,
    pub(crate) objective_contribution: Option<f64>,
    pub(crate) backend_slack: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct ConvexResidualEvidence {
    pub(crate) primal: f64,
    pub(crate) dual: f64,
    pub(crate) stationarity: f64,
    pub(crate) complementarity: f64,
    pub(crate) relative_gap: f64,
}

impl ConvexResidualEvidence {
    fn maximum(self) -> f64 {
        self.primal
            .max(self.dual)
            .max(self.stationarity)
            .max(self.complementarity)
            .max(self.relative_gap)
    }

    fn is_finite(self) -> bool {
        [
            self.primal,
            self.dual,
            self.stationarity,
            self.complementarity,
            self.relative_gap,
        ]
        .into_iter()
        .all(f64::is_finite)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum QpAttemptFailureReason {
    NonFiniteCandidate,
    BackendResidualExceeded,
    ThreadContractViolation,
    BackendFingerprintMismatch,
    UnverifiedTermination,
    InvalidInfeasibilityCertificate,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ValidatedInfeasibilityEvidence {
    pub(crate) finite: bool,
    pub(crate) normalized_ray_norm: f64,
    pub(crate) stationarity_residual: f64,
    pub(crate) dual_cone_violation: f64,
    pub(crate) separation_margin: f64,
    pub(crate) residual_limit: f64,
    pub(crate) separation_limit: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct QpAttemptRecord {
    pub(crate) backend: ClarabelAttemptEvidence,
    pub(crate) georbf_scaling: QpScalingEvidence,
    pub(crate) georbf_scaling_round_trip_error: f64,
    pub(crate) residuals: Option<ConvexResidualEvidence>,
    pub(crate) infeasibility_certificate: Option<ValidatedInfeasibilityEvidence>,
    pub(crate) failure_reason: Option<QpAttemptFailureReason>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct QpAttemptPlan {
    pub(crate) numerical_policy: NumericalPolicyId,
    pub(crate) profiles: [ClarabelAttemptProfile; 2],
    pub(crate) maximum_attempts: usize,
    pub(crate) canonical_tolerance_is_immutable: bool,
    pub(crate) form_family_is_immutable: bool,
}

impl QpAttemptPlan {
    fn georbf_v1() -> Self {
        Self {
            numerical_policy: EQUALITY_KKT_POLICY_V1.id,
            profiles: [
                ClarabelAttemptProfile::Standard,
                ClarabelAttemptProfile::Robust,
            ],
            maximum_attempts: 2,
            canonical_tolerance_is_immutable: true,
            form_family_is_immutable: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct QpScalingRoundEvidence {
    pub(crate) variable_exponents: Vec<i32>,
    pub(crate) constraint_exponents: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct QpScalingEvidence {
    pub(crate) rounds: Vec<QpScalingRoundEvidence>,
    pub(crate) cumulative_variable_exponents: Vec<i32>,
    pub(crate) cumulative_constraint_exponents: Vec<i32>,
    pub(crate) saturated_outside_target: usize,
}

impl QpScalingEvidence {
    fn variable_factors(&self) -> Vec<f64> {
        self.cumulative_variable_exponents
            .iter()
            .map(|exponent| 2.0_f64.powi(*exponent))
            .collect()
    }

    fn constraint_factors(&self) -> Vec<f64> {
        self.cumulative_constraint_exponents
            .iter()
            .map(|exponent| 2.0_f64.powi(*exponent))
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CubicQpEvidence {
    pub(crate) capacity: ConvexQpCapacityPlan,
    pub(crate) scaling: QpScalingEvidence,
    pub(crate) scaling_round_trip_error: f64,
    pub(crate) backend_residuals: ConvexResidualEvidence,
    pub(crate) backend_internal_scaling_round_trip_error: f64,
    pub(crate) attempts: Vec<QpAttemptRecord>,
    pub(crate) attempt_plan: QpAttemptPlan,
    pub(crate) accepted_attempt: usize,
    pub(crate) physical_standard_form_violation: f64,
    pub(crate) reduction_round_trip_error: f64,
    pub(crate) hard_relation_tolerances: Vec<CanonicalRelationToleranceEvidence>,
    pub(crate) affine_relation_tolerances: Vec<CanonicalRelationToleranceEvidence>,
    pub(crate) polynomial_round_trip_error: f64,
    pub(crate) field_coefficient_round_trip_error: f64,
    pub(crate) field_energy_round_trip_error: f64,
    pub(crate) whitening_round_trip_error: f64,
    pub(crate) objective_round_trip_error: f64,
    pub(crate) provenance_verified: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct CubicExecutionSolution {
    pub(crate) plan: CubicAlgebraicPlan,
    pub(crate) representation: CpdEvidence,
    pub(crate) field: RecoveredCubicField,
    pub(crate) semantic_latents: Vec<RecoveredSemanticLatent>,
    pub(crate) hard_equalities: Vec<RecoveredHardEquality>,
    pub(crate) affine_inequalities: Vec<RecoveredAffineInequality>,
    pub(crate) soft_equalities: Vec<RecoveredSoftEquality>,
    pub(crate) soft_objectives: Vec<RecoveredSoftObjective>,
    pub(crate) hard_equality_violations: FunctionalViolationEnvelope,
    pub(crate) side_condition: PhysicalSideConditionEvidence,
    pub(crate) field_energy: f64,
    pub(crate) total_objective: f64,
    pub(crate) backend_standard_form_verified: bool,
    pub(crate) canonical_acceptance_verified: bool,
    pub(crate) qp: Option<CubicQpEvidence>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum QpAssemblyFailureReason {
    InvalidCanonicalProblem,
    InvalidShape,
    ZeroScalingNorm,
    NonFiniteScalingNorm,
    ScalingLimitExceeded,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum QpRecoveryFailureReason {
    InvalidRecoveryMap,
    ProvenanceMismatch,
    NonFiniteRecoveredQuantity,
    SideConditionViolation,
    SideConditionRoundTripViolation,
    HardEqualityViolation,
    AffineInequalityViolation,
    BackendSlackMismatch,
    ReductionRoundTripViolation,
    ScalingRoundTripViolation,
    PolynomialRoundTripViolation,
    FieldCoefficientRoundTripViolation,
    FieldEnergyRoundTripViolation,
    WhiteningRoundTripViolation,
    ObjectiveRoundTripViolation,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct QpRecoveryFailureEvidence {
    pub(crate) reasons: Vec<QpRecoveryFailureReason>,
    pub(crate) backend_standard_form_violation: f64,
    pub(crate) reduction_round_trip_error: f64,
    pub(crate) scaling_round_trip_error: f64,
    pub(crate) no_model_produced: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CubicExecutionFailure {
    Equality(Box<CubicEqualityFailure>),
    Capacity(Box<CapacityExceededEvidence>),
    Representation(Box<RepresentationFailure>),
    Assembly(QpAssemblyFailureReason),
    BackendAdapter(Box<ClarabelAdapterFailure>),
    BackendContract {
        attempts: Vec<QpAttemptRecord>,
        observed: f64,
        limit: f64,
    },
    AttemptsExhausted {
        attempts: Vec<QpAttemptRecord>,
    },
    ValidatedInfeasible {
        evidence: ValidatedInfeasibilityEvidence,
        attempts: Vec<QpAttemptRecord>,
    },
    RecoveryVerification {
        evidence: Box<QpRecoveryFailureEvidence>,
        attempts: Vec<QpAttemptRecord>,
    },
}

pub(crate) struct CubicExecutionCore;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum QpFaultInjection {
    Provenance,
    ScalingMap,
    RecoveryMap,
    BackendResidual,
    BackendStandardForm,
    Objective,
    StandardRetry,
    InfeasibilityCertificate,
}

#[cfg(test)]
thread_local! {
    static INJECTED_QP_FAULT: RefCell<Option<QpFaultInjection>> = const { RefCell::new(None) };
}

#[cfg(test)]
pub(crate) fn inject_qp_fault_once(fault: QpFaultInjection) {
    INJECTED_QP_FAULT.with(|slot| {
        assert!(
            slot.replace(Some(fault)).is_none(),
            "only one QP fault may be injected per execution"
        );
    });
}

#[cfg(test)]
fn take_qp_fault() -> Option<QpFaultInjection> {
    INJECTED_QP_FAULT.with(|slot| slot.borrow_mut().take())
}

impl CubicExecutionCore {
    pub(crate) fn plan(problem: &CubicCanonicalProblem) -> CubicAlgebraicPlan {
        CubicAlgebraicPlan {
            form: if problem.affine_inequalities.is_empty() {
                CubicFormKind::SymmetricKkt
            } else {
                CubicFormKind::ConvexQp
            },
            hard_equalities: problem.equalities.len(),
            affine_inequalities: problem.affine_inequalities.len(),
        }
    }

    pub(crate) fn solve(
        problem: CubicCanonicalProblem,
        metric: GlobalAnisotropyMetric,
    ) -> Result<CubicExecutionSolution, CubicExecutionFailure> {
        let plan = Self::plan(&problem);
        match plan.form {
            CubicFormKind::SymmetricKkt => CubicEqualityCore::solve_canonical(problem, metric)
                .map(|solution| CubicExecutionSolution {
                    plan,
                    representation: solution.representation,
                    field: solution.field,
                    semantic_latents: solution.semantic_latents,
                    hard_equalities: solution.hard_equalities,
                    affine_inequalities: Vec::new(),
                    soft_equalities: solution.soft_equalities,
                    soft_objectives: solution.soft_objectives,
                    hard_equality_violations: solution.hard_equality_violations,
                    side_condition: solution.side_condition,
                    field_energy: solution.field_energy,
                    total_objective: solution.total_objective,
                    backend_standard_form_verified: true,
                    canonical_acceptance_verified: solution.recovery_finite
                        && solution.provenance_verified
                        && solution.objective_verified,
                    qp: None,
                })
                .map_err(|failure| CubicExecutionFailure::Equality(Box::new(failure))),
            CubicFormKind::ConvexQp => solve_convex_qp(plan, problem, metric),
        }
    }

    pub(crate) fn solve_equality_production(
        problem: CubicCanonicalProblem,
        metric: GlobalAnisotropyMetric,
    ) -> Result<CubicEqualitySolution, CubicEqualityFailure> {
        if Self::plan(&problem).form != CubicFormKind::SymmetricKkt {
            return Err(CubicEqualityFailure::AffineInequalityRequiresConvexQp);
        }
        CubicEqualityCore::solve_canonical(problem, metric)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct QpCanonicalRow {
    canonical_index: usize,
    provenance: UsageProvenance,
    derived_block: DerivedBlockId,
    derived_row: DerivedRowId,
    derived_column: Option<DerivedColumnId>,
    coefficients: Vec<f64>,
    rhs: f64,
    violation_variable: Option<usize>,
    provenance_edges: Vec<QpProvenanceEdge>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QpConeRole {
    Zero,
    Nonnegative,
}

#[derive(Debug, Clone, PartialEq)]
struct QpProvenanceEdge {
    provenance: UsageProvenance,
    derived_block: DerivedBlockId,
    derived_row: DerivedRowId,
    derived_column: Option<DerivedColumnId>,
    backend_row: usize,
    backend_column: Option<usize>,
    cone: QpConeRole,
}

#[derive(Debug, Clone, PartialEq)]
struct QpViolationNonnegativeRow {
    coefficients: Vec<f64>,
    rhs: f64,
    violation_variable: usize,
    provenance_edge: QpProvenanceEdge,
}

#[derive(Debug, Clone, PartialEq)]
struct QpSoftObjectiveBlock {
    objective_index: usize,
    canonical_indices: Vec<usize>,
    provenances: Vec<UsageProvenance>,
    residuals: Vec<ResidualId>,
    derived_blocks: Vec<DerivedBlockId>,
    derived_rows: Vec<DerivedRowId>,
    derived_columns: Vec<DerivedColumnId>,
    rows: Vec<Vec<f64>>,
    targets: Vec<f64>,
    precision: Vec<f64>,
    whitening: Vec<f64>,
    inverse_whitening: Vec<f64>,
    covariance_group: Option<GroupId>,
    block_kind: CanonicalSoftResidualBlockKind,
}

#[derive(Debug, Clone, PartialEq)]
struct ConvexQpForm {
    variables: usize,
    reduced_field_variables: usize,
    polynomial_variables: usize,
    semantic_latents: usize,
    soft_violation_variables: usize,
    hessian: Vec<f64>,
    linear_objective: Vec<f64>,
    constraints: Vec<f64>,
    constraint_rhs: Vec<f64>,
    equality_constraints: usize,
    inequality_constraints: usize,
    hard_equality_rows: Vec<QpCanonicalRow>,
    affine_inequality_rows: Vec<QpCanonicalRow>,
    violation_nonnegative_rows: Vec<QpViolationNonnegativeRow>,
    soft_objective_blocks: Vec<QpSoftObjectiveBlock>,
    capacity: ConvexQpCapacityPlan,
}

#[derive(Debug, Clone, PartialEq)]
struct ScaledQpForm {
    hessian: Vec<f64>,
    linear_objective: Vec<f64>,
    constraints: Vec<f64>,
    constraint_rhs: Vec<f64>,
    equality_constraints: usize,
    inequality_constraints: usize,
    scaling: QpScalingEvidence,
    round_trip_error: f64,
}

fn solve_convex_qp(
    plan: CubicAlgebraicPlan,
    problem: CubicCanonicalProblem,
    metric: GlobalAnisotropyMetric,
) -> Result<CubicExecutionSolution, CubicExecutionFailure> {
    #[cfg(test)]
    let fault = take_qp_fault();
    #[cfg(not(test))]
    let fault: Option<QpFaultInjection> = None;
    validate_qp_problem(&problem)?;
    let fitting_uses = canonical_fitting_uses(
        &problem.equalities,
        &problem.soft_equalities,
        &problem.affine_inequalities,
    );
    let variables = fitting_uses
        .len()
        .checked_add(problem.semantic_latents.len())
        .and_then(|count| {
            count.checked_add(
                problem
                    .affine_inequalities
                    .iter()
                    .filter(|relation| relation.violation_channel().is_some())
                    .count(),
            )
        })
        .ok_or(QpAssemblyFailureReason::InvalidShape)
        .map_err(CubicExecutionFailure::Assembly)?;
    let equality_constraints = problem
        .equalities
        .iter()
        .filter(|equality| {
            equality.participation() == CanonicalEqualityParticipation::SolverConstraint
        })
        .count();
    let constraints = equality_constraints
        .checked_add(problem.affine_inequalities.len())
        .and_then(|count| {
            count.checked_add(
                problem
                    .affine_inequalities
                    .iter()
                    .filter(|relation| relation.violation_channel().is_some())
                    .count(),
            )
        })
        .ok_or(QpAssemblyFailureReason::InvalidShape)
        .map_err(CubicExecutionFailure::Assembly)?;
    let canonical_relations = problem
        .equalities
        .len()
        .checked_add(problem.soft_equalities.len())
        .and_then(|count| count.checked_add(problem.affine_inequalities.len()))
        .ok_or(QpAssemblyFailureReason::InvalidShape)
        .map_err(CubicExecutionFailure::Assembly)?;
    let capacity = plan_convex_qp_capacity(variables, constraints, canonical_relations)
        .map_err(|failure| CubicExecutionFailure::Capacity(Box::new(failure)))?;
    let mut representation = CubicRepresentation::new(fitting_uses, metric)
        .map_err(|failure| CubicExecutionFailure::Representation(Box::new(failure)))?;
    representation.set_field_energy_normalization(problem.field_energy_normalization);
    let mut form = assemble_qp_form(&representation, &problem, capacity)?;
    if matches!(fault, Some(QpFaultInjection::Provenance)) {
        form.hard_equality_rows[0].derived_column = Some(DerivedColumnId::from_residual(
            problem.equalities[1].provenance().residual(),
        ));
    }
    let mut scaled = scale_qp_form(&form)?;
    if matches!(fault, Some(QpFaultInjection::ScalingMap)) {
        scaled.round_trip_error = 1.0e-2;
    }
    let (candidate, residuals, attempts, accepted_attempt, internal_scaling_error) =
        execute_qp_attempts(&form, &scaled, fault)?;
    recover_and_verify_qp(
        plan,
        representation,
        problem,
        form,
        scaled,
        candidate,
        residuals,
        attempts,
        accepted_attempt,
        internal_scaling_error,
        fault,
    )
}

fn validate_qp_problem(problem: &CubicCanonicalProblem) -> Result<(), CubicExecutionFailure> {
    if problem.affine_inequalities.is_empty()
        || (problem.equalities.is_empty() && problem.soft_equalities.is_empty())
    {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::InvalidCanonicalProblem,
        ));
    }
    if problem.equalities.iter().any(|equality| {
        !equality.target().is_finite()
            || equality.latent_coefficients().iter().any(|term| {
                !term.coefficient.is_finite() || term.latent >= problem.semantic_latents.len()
            })
    }) || problem.affine_inequalities.iter().any(|inequality| {
        !inequality.bound().is_finite()
            || inequality.latent_coefficients().iter().any(|term| {
                !term.coefficient.is_finite() || term.latent >= problem.semantic_latents.len()
            })
            || inequality
                .violation_channel()
                .is_some_and(|channel| !channel.loss().is_valid())
    }) {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::InvalidCanonicalProblem,
        ));
    }
    let objective_residuals = problem
        .soft_objectives
        .iter()
        .flat_map(|objective| objective.residuals())
        .collect::<Vec<_>>();
    if problem
        .soft_equalities
        .iter()
        .any(|relation| !relation.target().is_finite())
        || objective_residuals.len() != problem.soft_equalities.len()
        || problem
            .soft_equalities
            .iter()
            .zip(objective_residuals)
            .any(|(relation, residual)| relation.provenance().residual() != residual)
        || problem.soft_objectives.iter().any(|objective| {
            !objective.loss().is_valid(objective.residuals().len())
                || !objective
                    .block_kind()
                    .is_valid(objective.residuals().len(), objective.covariance_group())
        })
    {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::InvalidCanonicalProblem,
        ));
    }
    Ok(())
}

fn assemble_qp_form(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    capacity: ConvexQpCapacityPlan,
) -> Result<ConvexQpForm, CubicExecutionFailure> {
    let reduced_field_variables = representation.null_space().reduced_dimension();
    let soft_violation_variables = problem
        .affine_inequalities
        .iter()
        .filter(|relation| relation.violation_channel().is_some())
        .count();
    let variables = reduced_field_variables
        .checked_add(POLYNOMIAL_DIMENSION)
        .and_then(|count| count.checked_add(problem.semantic_latents.len()))
        .and_then(|count| count.checked_add(soft_violation_variables))
        .ok_or(QpAssemblyFailureReason::InvalidShape)
        .map_err(CubicExecutionFailure::Assembly)?;
    if variables != capacity.variables {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::InvalidShape,
        ));
    }
    let hard_equalities = problem
        .equalities
        .iter()
        .enumerate()
        .filter(|(_, equality)| {
            equality.participation() == CanonicalEqualityParticipation::SolverConstraint
        })
        .collect::<Vec<_>>();
    let mut hard_equality_rows = Vec::with_capacity(hard_equalities.len());
    for (canonical_index, equality) in hard_equalities {
        let backend_row = hard_equality_rows.len();
        hard_equality_rows.push(QpCanonicalRow {
            canonical_index,
            provenance: equality.provenance().clone(),
            derived_block: DerivedBlockId::from_residual(equality.provenance().residual()),
            derived_row: DerivedRowId::from_residual(equality.provenance().residual()),
            derived_column: equality
                .field()
                .map(|_| DerivedColumnId::from_residual(equality.provenance().residual())),
            coefficients: reduced_affine_row(
                representation,
                problem,
                equality.field(),
                equality.latent_coefficients(),
            )?,
            rhs: equality.target(),
            violation_variable: None,
            provenance_edges: vec![QpProvenanceEdge {
                provenance: equality.provenance().clone(),
                derived_block: DerivedBlockId::from_residual(equality.provenance().residual()),
                derived_row: DerivedRowId::from_residual(equality.provenance().residual()),
                derived_column: equality
                    .field()
                    .map(|_| DerivedColumnId::from_residual(equality.provenance().residual())),
                backend_row,
                backend_column: None,
                cone: QpConeRole::Zero,
            }],
        });
    }
    let mut affine_inequality_rows = Vec::with_capacity(problem.affine_inequalities.len());
    let base_variables = variables - soft_violation_variables;
    let mut next_violation_variable = base_variables;
    for (canonical_index, inequality) in problem.affine_inequalities.iter().enumerate() {
        let multiplier = inequality.sense().upper_form_multiplier();
        let mut coefficients = reduced_affine_row(
            representation,
            problem,
            inequality.field(),
            inequality.latent_coefficients(),
        )?;
        coefficients.resize(variables, 0.0);
        coefficients
            .iter_mut()
            .for_each(|entry| *entry *= multiplier);
        let violation_variable = inequality.violation_channel().map(|_| {
            let variable = next_violation_variable;
            next_violation_variable += 1;
            coefficients[variable] = -1.0;
            variable
        });
        affine_inequality_rows.push(QpCanonicalRow {
            canonical_index,
            provenance: inequality.provenance().clone(),
            derived_block: DerivedBlockId::from_residual(inequality.provenance().residual()),
            derived_row: DerivedRowId::from_residual(inequality.provenance().residual()),
            derived_column: inequality
                .field()
                .map(|_| DerivedColumnId::from_residual(inequality.provenance().residual())),
            coefficients,
            rhs: inequality.upper_form_bound(),
            violation_variable,
            provenance_edges: inequality
                .source_provenances()
                .iter()
                .map(|provenance| QpProvenanceEdge {
                    provenance: provenance.clone(),
                    derived_block: DerivedBlockId::from_residual(provenance.residual()),
                    derived_row: DerivedRowId::from_residual(provenance.residual()),
                    derived_column: inequality
                        .field()
                        .map(|_| DerivedColumnId::from_residual(provenance.residual())),
                    backend_row: hard_equality_rows.len() + canonical_index,
                    backend_column: violation_variable,
                    cone: QpConeRole::Nonnegative,
                })
                .collect(),
        });
    }
    let violation_nonnegative_rows = affine_inequality_rows
        .iter()
        .filter_map(|row| row.violation_variable.map(|variable| (row, variable)))
        .enumerate()
        .map(|(soft_index, (relation_row, variable))| {
            let mut coefficients = vec![0.0; variables];
            coefficients[variable] = -1.0;
            let provenance = relation_row
                .provenance_edges
                .first()
                .expect("a soft inequality retains source provenance")
                .provenance
                .clone();
            QpViolationNonnegativeRow {
                coefficients,
                rhs: 0.0,
                violation_variable: variable,
                provenance_edge: QpProvenanceEdge {
                    derived_block: DerivedBlockId::from_residual(provenance.residual()),
                    derived_row: DerivedRowId::from_residual(provenance.residual()),
                    derived_column: Some(DerivedColumnId::from_residual(provenance.residual())),
                    provenance,
                    backend_row: hard_equality_rows.len()
                        + affine_inequality_rows.len()
                        + soft_index,
                    backend_column: Some(variable),
                    cone: QpConeRole::Nonnegative,
                },
            }
        })
        .collect::<Vec<_>>();
    let soft_objective_blocks = assemble_qp_soft_objectives(representation, problem)?;
    let field_energy_scale = representation.field_energy_normalization().factor()
        / representation.solve_coordinate_transform().length().powi(3);
    let mut hessian = vec![0.0; variables * variables];
    for row in 0..reduced_field_variables {
        for column in 0..reduced_field_variables {
            hessian[row * variables + column] =
                field_energy_scale * representation.reduced_kernel_pairing().get(row, column);
        }
    }
    let mut linear_objective = vec![0.0; variables];
    for objective in &soft_objective_blocks {
        let dimension = objective.rows.len();
        for row in 0..variables {
            for column in 0..variables {
                hessian[row * variables + column] += (0..dimension)
                    .flat_map(|left| {
                        (0..dimension).map(move |right| {
                            objective.rows[left][row]
                                * objective.precision[left * dimension + right]
                                * objective.rows[right][column]
                        })
                    })
                    .sum::<f64>();
            }
            linear_objective[row] -= (0..dimension)
                .flat_map(|left| {
                    (0..dimension).map(move |right| {
                        objective.rows[left][row]
                            * objective.precision[left * dimension + right]
                            * objective.targets[right]
                    })
                })
                .sum::<f64>();
        }
    }
    for (relation, row) in problem
        .affine_inequalities
        .iter()
        .zip(&affine_inequality_rows)
    {
        let (Some(channel), Some(variable)) =
            (relation.violation_channel(), row.violation_variable)
        else {
            continue;
        };
        match channel.loss() {
            CanonicalViolationLoss::QuadraticPenalty { weight } => {
                hessian[variable * variables + variable] += weight;
            }
            CanonicalViolationLoss::LinearViolationPenalty { weight } => {
                linear_objective[variable] += weight;
            }
        }
    }
    let rows = hard_equality_rows
        .iter()
        .chain(&affine_inequality_rows)
        .map(|row| row.coefficients.as_slice())
        .chain(
            violation_nonnegative_rows
                .iter()
                .map(|row| row.coefficients.as_slice()),
        )
        .collect::<Vec<_>>();
    let constraints = rows
        .iter()
        .flat_map(|row| row.iter().copied())
        .collect::<Vec<_>>();
    let constraint_rhs = hard_equality_rows
        .iter()
        .map(|row| row.rhs)
        .chain(affine_inequality_rows.iter().map(|row| row.rhs))
        .chain(violation_nonnegative_rows.iter().map(|row| row.rhs))
        .collect::<Vec<_>>();
    Ok(ConvexQpForm {
        variables,
        reduced_field_variables,
        polynomial_variables: POLYNOMIAL_DIMENSION,
        semantic_latents: problem.semantic_latents.len(),
        soft_violation_variables,
        hessian,
        linear_objective,
        constraints,
        constraint_rhs,
        equality_constraints: hard_equality_rows.len(),
        inequality_constraints: affine_inequality_rows.len() + violation_nonnegative_rows.len(),
        hard_equality_rows,
        affine_inequality_rows,
        violation_nonnegative_rows,
        soft_objective_blocks,
        capacity,
    })
}

fn assemble_qp_soft_objectives(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
) -> Result<Vec<QpSoftObjectiveBlock>, CubicExecutionFailure> {
    let mut canonical_offset = 0;
    let mut blocks = Vec::with_capacity(problem.soft_objectives.len());
    for (objective_index, objective) in problem.soft_objectives.iter().enumerate() {
        let dimension = objective.residuals().len();
        let canonical_indices =
            (canonical_offset..canonical_offset + dimension).collect::<Vec<_>>();
        canonical_offset += dimension;
        let relations = canonical_indices
            .iter()
            .map(|index| &problem.soft_equalities[*index])
            .collect::<Vec<_>>();
        let mut rows = Vec::with_capacity(dimension);
        for relation in &relations {
            rows.push(reduced_affine_row(
                representation,
                problem,
                Some(relation.field()),
                &[],
            )?);
        }
        blocks.push(QpSoftObjectiveBlock {
            objective_index,
            canonical_indices,
            provenances: relations
                .iter()
                .map(|relation| relation.provenance().clone())
                .collect(),
            residuals: objective.residuals().to_vec(),
            derived_blocks: relations
                .iter()
                .map(|relation| DerivedBlockId::from_residual(relation.provenance().residual()))
                .collect(),
            derived_rows: relations
                .iter()
                .map(|relation| DerivedRowId::from_residual(relation.provenance().residual()))
                .collect(),
            derived_columns: relations
                .iter()
                .map(|relation| DerivedColumnId::from_residual(relation.provenance().residual()))
                .collect(),
            rows,
            targets: relations.iter().map(|relation| relation.target()).collect(),
            precision: objective.loss().precision_matrix(dimension),
            whitening: objective.loss().whitening_matrix(dimension),
            inverse_whitening: objective.loss().inverse_whitening_matrix(dimension),
            covariance_group: objective.covariance_group().cloned(),
            block_kind: objective.block_kind().clone(),
        });
    }
    Ok(blocks)
}

fn reduced_affine_row(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    field: Option<&crate::functional::FunctionalUse>,
    latent_coefficients: &[SemanticLatentCoefficient],
) -> Result<Vec<f64>, CubicExecutionFailure> {
    let reduced_dimension = representation.null_space().reduced_dimension();
    let mut row = vec![0.0; reduced_dimension + POLYNOMIAL_DIMENSION];
    if let Some(field) = field {
        let functional_index = representation
            .fitting_uses()
            .iter()
            .position(|usage| usage.functional() == field.functional())
            .ok_or(QpAssemblyFailureReason::InvalidCanonicalProblem)
            .map_err(CubicExecutionFailure::Assembly)?;
        let ambient = (0..representation.kernel_pairing().shape().1)
            .map(|column| {
                representation
                    .kernel_pairing()
                    .get(functional_index, column)
            })
            .collect::<Vec<_>>();
        let reduced = representation
            .null_space()
            .project(&ambient)
            .map_err(|failure| CubicExecutionFailure::Representation(Box::new(failure)))?;
        row[..reduced_dimension].copy_from_slice(&reduced);
        for polynomial in 0..POLYNOMIAL_DIMENSION {
            row[reduced_dimension + polynomial] = representation
                .polynomial_pairing()
                .get(functional_index, polynomial);
        }
    }
    row.extend(std::iter::repeat_n(0.0, problem.semantic_latents.len()));
    for term in latent_coefficients {
        row[reduced_dimension + POLYNOMIAL_DIMENSION + term.latent] = term.coefficient;
    }
    row.extend(std::iter::repeat_n(
        0.0,
        problem
            .affine_inequalities
            .iter()
            .filter(|relation| relation.violation_channel().is_some())
            .count(),
    ));
    Ok(row)
}

fn scale_qp_form(form: &ConvexQpForm) -> Result<ScaledQpForm, CubicExecutionFailure> {
    let constraints_count = form.equality_constraints + form.inequality_constraints;
    if form.hessian.len() != form.variables.saturating_mul(form.variables)
        || form.linear_objective.len() != form.variables
        || form.constraints.len() != constraints_count.saturating_mul(form.variables)
        || form.constraint_rhs.len() != constraints_count
    {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::InvalidShape,
        ));
    }
    let original_hessian = form.hessian.clone();
    let original_objective = form.linear_objective.clone();
    let original_constraints = form.constraints.clone();
    let original_rhs = form.constraint_rhs.clone();
    let mut hessian = original_hessian.clone();
    let mut linear_objective = original_objective.clone();
    let mut constraints = original_constraints.clone();
    let mut constraint_rhs = original_rhs.clone();
    let mut cumulative_variable_exponents = vec![0_i32; form.variables];
    let mut cumulative_constraint_exponents = vec![0_i32; constraints_count];
    let mut rounds = Vec::with_capacity(EQUALITY_KKT_POLICY_V1.ruiz_rounds);
    for _ in 0..EQUALITY_KKT_POLICY_V1.ruiz_rounds {
        let variable_exponents = (0..form.variables)
            .map(|column| {
                let hessian_norm = (0..form.variables)
                    .map(|row| hessian[row * form.variables + column].abs())
                    .fold(0.0_f64, f64::max);
                let constraint_norm = (0..constraints_count)
                    .map(|row| constraints[row * form.variables + column].abs())
                    .fold(0.0_f64, f64::max);
                bounded_ruiz_exponent(
                    hessian_norm.max(constraint_norm),
                    cumulative_variable_exponents[column],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let variable_factors = variable_exponents
            .iter()
            .map(|exponent| 2.0_f64.powi(*exponent))
            .collect::<Vec<_>>();
        for row in 0..form.variables {
            for column in 0..form.variables {
                hessian[row * form.variables + column] *=
                    variable_factors[row] * variable_factors[column];
            }
        }
        for column in 0..form.variables {
            linear_objective[column] *= variable_factors[column];
            cumulative_variable_exponents[column] += variable_exponents[column];
        }
        for row in 0..constraints_count {
            for column in 0..form.variables {
                constraints[row * form.variables + column] *= variable_factors[column];
            }
        }
        let constraint_exponents = (0..constraints_count)
            .map(|row| {
                let norm = (0..form.variables)
                    .map(|column| constraints[row * form.variables + column].abs())
                    .fold(0.0_f64, f64::max);
                bounded_ruiz_exponent(norm, cumulative_constraint_exponents[row])
            })
            .collect::<Result<Vec<_>, _>>()?;
        for row in 0..constraints_count {
            let factor = 2.0_f64.powi(constraint_exponents[row]);
            cumulative_constraint_exponents[row] += constraint_exponents[row];
            constraint_rhs[row] *= factor;
            for column in 0..form.variables {
                constraints[row * form.variables + column] *= factor;
            }
        }
        rounds.push(QpScalingRoundEvidence {
            variable_exponents,
            constraint_exponents,
        });
    }
    if cumulative_variable_exponents
        .iter()
        .chain(&cumulative_constraint_exponents)
        .any(|exponent| exponent.abs() > EQUALITY_KKT_POLICY_V1.ruiz_cumulative_exponent_limit)
    {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::ScalingLimitExceeded,
        ));
    }
    let scaling = QpScalingEvidence {
        saturated_outside_target: count_saturated_qp_rows(
            form.variables,
            constraints_count,
            &hessian,
            &constraints,
            &cumulative_variable_exponents,
            &cumulative_constraint_exponents,
        ),
        rounds,
        cumulative_variable_exponents,
        cumulative_constraint_exponents,
    };
    let variable_factors = scaling.variable_factors();
    let constraint_factors = scaling.constraint_factors();
    let round_trip_error = recover_scaled_form_error(
        form.variables,
        constraints_count,
        &hessian,
        &linear_objective,
        &constraints,
        &constraint_rhs,
        &variable_factors,
        &constraint_factors,
        &original_hessian,
        &original_objective,
        &original_constraints,
        &original_rhs,
    );
    if !round_trip_error.is_finite()
        || round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit
    {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::ScalingLimitExceeded,
        ));
    }
    Ok(ScaledQpForm {
        hessian,
        linear_objective,
        constraints,
        constraint_rhs,
        equality_constraints: form.equality_constraints,
        inequality_constraints: form.inequality_constraints,
        scaling,
        round_trip_error,
    })
}

fn bounded_ruiz_exponent(norm: f64, cumulative: i32) -> Result<i32, CubicExecutionFailure> {
    if norm == 0.0 {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::ZeroScalingNorm,
        ));
    }
    if !norm.is_finite() {
        return Err(CubicExecutionFailure::Assembly(
            QpAssemblyFailureReason::NonFiniteScalingNorm,
        ));
    }
    let desired = (-0.5 * norm.log2()).round();
    let desired = desired.clamp(i32::MIN as f64, i32::MAX as f64) as i32;
    let round_limit = EQUALITY_KKT_POLICY_V1.ruiz_single_round_exponent_limit;
    let cumulative_limit = EQUALITY_KKT_POLICY_V1.ruiz_cumulative_exponent_limit;
    Ok(desired.clamp(-round_limit, round_limit).clamp(
        -cumulative_limit - cumulative,
        cumulative_limit - cumulative,
    ))
}

fn count_saturated_qp_rows(
    variables: usize,
    constraints_count: usize,
    hessian: &[f64],
    constraints: &[f64],
    variable_exponents: &[i32],
    constraint_exponents: &[i32],
) -> usize {
    let limit = EQUALITY_KKT_POLICY_V1.ruiz_cumulative_exponent_limit;
    let variable_count = (0..variables)
        .filter(|column| {
            let norm = (0..variables)
                .map(|row| hessian[row * variables + *column].abs())
                .chain(
                    (0..constraints_count).map(|row| constraints[row * variables + *column].abs()),
                )
                .fold(0.0_f64, f64::max);
            !(0.5..=2.0).contains(&norm) && variable_exponents[*column].abs() == limit
        })
        .count();
    variable_count
        + (0..constraints_count)
            .filter(|row| {
                let norm = (0..variables)
                    .map(|column| constraints[*row * variables + column].abs())
                    .fold(0.0_f64, f64::max);
                !(0.5..=2.0).contains(&norm) && constraint_exponents[*row].abs() == limit
            })
            .count()
}

#[allow(clippy::too_many_arguments)]
fn recover_scaled_form_error(
    variables: usize,
    constraints_count: usize,
    hessian: &[f64],
    linear_objective: &[f64],
    constraints: &[f64],
    rhs: &[f64],
    variable_factors: &[f64],
    constraint_factors: &[f64],
    original_hessian: &[f64],
    original_objective: &[f64],
    original_constraints: &[f64],
    original_rhs: &[f64],
) -> f64 {
    let hessian_error = (0..variables * variables)
        .map(|index| {
            let row = index / variables;
            let column = index % variables;
            relative_error(
                hessian[index] / (variable_factors[row] * variable_factors[column]),
                original_hessian[index],
            )
        })
        .fold(0.0_f64, f64::max);
    let objective_error = (0..variables)
        .map(|column| {
            relative_error(
                linear_objective[column] / variable_factors[column],
                original_objective[column],
            )
        })
        .fold(0.0_f64, f64::max);
    let constraint_error = (0..constraints_count * variables)
        .map(|index| {
            let row = index / variables;
            let column = index % variables;
            relative_error(
                constraints[index] / (constraint_factors[row] * variable_factors[column]),
                original_constraints[index],
            )
        })
        .fold(0.0_f64, f64::max);
    let rhs_error = (0..constraints_count)
        .map(|row| relative_error(rhs[row] / constraint_factors[row], original_rhs[row]))
        .fold(0.0_f64, f64::max);
    hessian_error
        .max(objective_error)
        .max(constraint_error)
        .max(rhs_error)
}

fn execute_qp_attempts(
    form: &ConvexQpForm,
    scaled: &ScaledQpForm,
    fault: Option<QpFaultInjection>,
) -> Result<
    (
        ClarabelCandidateEnvelope,
        ConvexResidualEvidence,
        Vec<QpAttemptRecord>,
        usize,
        f64,
    ),
    CubicExecutionFailure,
> {
    let variables = scaled.linear_objective.len();
    let constraints_count = scaled.constraint_rhs.len();
    let mut attempts = Vec::with_capacity(2);
    for (sequence, profile) in [
        ClarabelAttemptProfile::Standard,
        ClarabelAttemptProfile::Robust,
    ]
    .into_iter()
    .enumerate()
    {
        let mut candidate = clarabel_backend::solve_qp(
            ClarabelQpInput {
                variables,
                constraints: constraints_count,
                hessian: &scaled.hessian,
                linear_objective: &scaled.linear_objective,
                constraint_matrix: &scaled.constraints,
                constraint_rhs: &scaled.constraint_rhs,
                equality_constraints: scaled.equality_constraints,
                inequality_constraints: scaled.inequality_constraints,
            },
            profile,
            sequence,
        )
        .map_err(|failure| CubicExecutionFailure::BackendAdapter(Box::new(failure)))?;
        if matches!(fault, Some(QpFaultInjection::InfeasibilityCertificate))
            && matches!(
                candidate.attempt.termination,
                ClarabelTermination::PrimalInfeasible | ClarabelTermination::AlmostPrimalInfeasible
            )
        {
            candidate.dual.fill(0.0);
        }
        let finite = candidate
            .primal
            .iter()
            .chain(&candidate.dual)
            .chain(&candidate.slack)
            .all(|value| value.is_finite());
        if sequence == 0 && matches!(fault, Some(QpFaultInjection::StandardRetry)) {
            candidate.attempt.termination = ClarabelTermination::AlmostSolved;
        }
        let candidate_termination = matches!(
            candidate.attempt.termination,
            ClarabelTermination::Solved | ClarabelTermination::AlmostSolved
        );
        let infeasibility_termination = matches!(
            candidate.attempt.termination,
            ClarabelTermination::PrimalInfeasible | ClarabelTermination::AlmostPrimalInfeasible
        );
        let mut residuals =
            (finite && candidate_termination).then(|| convex_residuals(scaled, &candidate));
        if sequence == 0 && matches!(fault, Some(QpFaultInjection::StandardRetry)) {
            residuals = residuals.map(|mut residuals| {
                residuals.primal = 1.0;
                residuals
            });
        }
        if sequence == 0 && matches!(fault, Some(QpFaultInjection::BackendResidual)) {
            residuals = residuals.map(|mut residuals| {
                residuals.primal = 1.0;
                residuals
            });
        }
        let internal_scaling_error = backend_internal_scaling_round_trip_error(&candidate);
        let fingerprint_verified = candidate.attempt.backend.crate_name
            == clarabel_backend::CRATE_NAME
            && candidate.attempt.backend.crate_version == clarabel_backend::CRATE_VERSION
            && candidate.attempt.backend.features == clarabel_backend::FEATURES
            && candidate.attempt.backend.direct_solver == clarabel_backend::DIRECT_SOLVER
            && candidate.attempt.settings == clarabel_backend::expected_settings(profile)
            && candidate.attempt.requested_threads == clarabel_backend::REQUESTED_THREADS
            && candidate
                .attempt
                .linear_solver
                .to_ascii_lowercase()
                .contains(clarabel_backend::DIRECT_SOLVER)
            && internal_scaling_error.is_finite()
            && internal_scaling_error <= EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit;
        let threads_verified =
            candidate.attempt.actual_threads == clarabel_backend::REQUESTED_THREADS;
        let residual_verified = residuals.is_some_and(|residuals| {
            residuals.is_finite()
                && residuals.maximum() <= EQUALITY_KKT_POLICY_V1.convex_backend_residual_limit
        });
        let infeasibility_certificate =
            (finite && threads_verified && fingerprint_verified && infeasibility_termination)
                .then(|| validate_primal_infeasibility_certificate(form, scaled, &candidate))
                .flatten();
        let failure_reason = if !finite {
            Some(QpAttemptFailureReason::NonFiniteCandidate)
        } else if !threads_verified {
            Some(QpAttemptFailureReason::ThreadContractViolation)
        } else if !fingerprint_verified {
            Some(QpAttemptFailureReason::BackendFingerprintMismatch)
        } else if candidate_termination && !residual_verified {
            Some(QpAttemptFailureReason::BackendResidualExceeded)
        } else if infeasibility_termination && infeasibility_certificate.is_none() {
            Some(QpAttemptFailureReason::InvalidInfeasibilityCertificate)
        } else if !candidate_termination {
            if infeasibility_certificate.is_some() {
                None
            } else {
                Some(QpAttemptFailureReason::UnverifiedTermination)
            }
        } else {
            None
        };
        attempts.push(QpAttemptRecord {
            backend: candidate.attempt.clone(),
            georbf_scaling: scaled.scaling.clone(),
            georbf_scaling_round_trip_error: scaled.round_trip_error,
            residuals,
            infeasibility_certificate,
            failure_reason,
        });
        if let Some(evidence) = infeasibility_certificate {
            return Err(CubicExecutionFailure::ValidatedInfeasible { evidence, attempts });
        }
        if failure_reason.is_none() {
            return Ok((
                candidate,
                residuals.expect("an accepted candidate has residual evidence"),
                attempts,
                sequence,
                internal_scaling_error,
            ));
        }
        if matches!(candidate.attempt.termination, ClarabelTermination::Solved)
            || matches!(
                failure_reason,
                Some(
                    QpAttemptFailureReason::ThreadContractViolation
                        | QpAttemptFailureReason::BackendFingerprintMismatch
                )
            )
        {
            let observed = residuals
                .map(ConvexResidualEvidence::maximum)
                .unwrap_or(f64::INFINITY);
            return Err(CubicExecutionFailure::BackendContract {
                attempts,
                observed,
                limit: EQUALITY_KKT_POLICY_V1.convex_backend_residual_limit,
            });
        }
    }
    Err(CubicExecutionFailure::AttemptsExhausted { attempts })
}

fn validate_primal_infeasibility_certificate(
    form: &ConvexQpForm,
    scaled: &ScaledQpForm,
    candidate: &ClarabelCandidateEnvelope,
) -> Option<ValidatedInfeasibilityEvidence> {
    let constraint_factors = scaled.scaling.constraint_factors();
    if constraint_factors.len() != candidate.dual.len() {
        return None;
    }
    let recovered_ray = candidate
        .dual
        .iter()
        .zip(constraint_factors)
        .map(|(value, factor)| value * factor)
        .collect::<Vec<_>>();
    validate_primal_infeasibility_ray(
        form.variables,
        &form.constraints,
        &form.constraint_rhs,
        form.equality_constraints,
        &recovered_ray,
    )
}

fn validate_primal_infeasibility_ray(
    variables: usize,
    constraints: &[f64],
    constraint_rhs: &[f64],
    equality_constraints: usize,
    candidate_ray: &[f64],
) -> Option<ValidatedInfeasibilityEvidence> {
    let constraints_count = constraint_rhs.len();
    if candidate_ray.len() != constraints_count
        || constraints.len() != constraints_count.saturating_mul(variables)
        || equality_constraints > constraints_count
    {
        return None;
    }
    let ray_norm = candidate_ray
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if !ray_norm.is_finite() || ray_norm <= 0.0 {
        return None;
    }
    let ray = candidate_ray
        .iter()
        .map(|value| value / ray_norm)
        .collect::<Vec<_>>();
    let normalized_ray_norm = ray.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    let matrix_scale = constraints
        .iter()
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let stationarity_residual = (0..variables)
        .map(|column| {
            (0..constraints_count)
                .map(|row| constraints[row * variables + column] * ray[row])
                .sum::<f64>()
                .abs()
        })
        .fold(0.0_f64, f64::max)
        / matrix_scale;
    let dual_cone_violation = ray[equality_constraints..]
        .iter()
        .map(|value| (-value).max(0.0))
        .fold(0.0_f64, f64::max);
    let rhs_scale = constraint_rhs
        .iter()
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let separation_margin = -dot_product(constraint_rhs, &ray) / rhs_scale;
    let evidence = ValidatedInfeasibilityEvidence {
        finite: normalized_ray_norm.is_finite()
            && stationarity_residual.is_finite()
            && dual_cone_violation.is_finite()
            && separation_margin.is_finite(),
        normalized_ray_norm,
        stationarity_residual,
        dual_cone_violation,
        separation_margin,
        residual_limit: EQUALITY_KKT_POLICY_V1.convex_certificate_residual_limit,
        separation_limit: EQUALITY_KKT_POLICY_V1.convex_certificate_separation_limit,
    };
    (evidence.finite
        && (evidence.normalized_ray_norm - 1.0).abs()
            <= EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit
        && evidence.stationarity_residual <= evidence.residual_limit
        && evidence.dual_cone_violation <= evidence.residual_limit
        && evidence.separation_margin >= evidence.separation_limit)
        .then_some(evidence)
}

fn convex_residuals(
    form: &ScaledQpForm,
    candidate: &ClarabelCandidateEnvelope,
) -> ConvexResidualEvidence {
    let variables = form.linear_objective.len();
    let constraints_count = form.constraint_rhs.len();
    if candidate.primal.len() != variables
        || candidate.dual.len() != constraints_count
        || candidate.slack.len() != constraints_count
    {
        return ConvexResidualEvidence {
            primal: f64::INFINITY,
            dual: f64::INFINITY,
            stationarity: f64::INFINITY,
            complementarity: f64::INFINITY,
            relative_gap: f64::INFINITY,
        };
    }
    let affine = dense_matrix_vector_product(
        &form.constraints,
        constraints_count,
        variables,
        &candidate.primal,
    );
    let equation_residual = affine
        .iter()
        .zip(&candidate.slack)
        .zip(&form.constraint_rhs)
        .map(|((affine, slack), rhs)| (affine + slack - rhs).abs())
        .fold(0.0_f64, f64::max);
    let primal_scale = affine
        .iter()
        .chain(&candidate.slack)
        .chain(&form.constraint_rhs)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let zero_slack_violation = candidate.slack[..form.equality_constraints]
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let nonnegative_slack_violation = candidate.slack[form.equality_constraints..]
        .iter()
        .map(|value| (-value).max(0.0))
        .fold(0.0_f64, f64::max);
    let primal = (equation_residual / primal_scale)
        .max(zero_slack_violation / primal_scale)
        .max(nonnegative_slack_violation / primal_scale);
    let hessian_product =
        dense_matrix_vector_product(&form.hessian, variables, variables, &candidate.primal);
    let stationarity_vector = (0..variables)
        .map(|column| {
            hessian_product[column]
                + form.linear_objective[column]
                + (0..constraints_count)
                    .map(|row| form.constraints[row * variables + column] * candidate.dual[row])
                    .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let stationarity_scale = hessian_product
        .iter()
        .chain(&form.linear_objective)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let stationarity = stationarity_vector
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max)
        / stationarity_scale;
    let dual_scale = candidate
        .dual
        .iter()
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let dual = candidate.dual[form.equality_constraints..]
        .iter()
        .map(|value| (-value).max(0.0))
        .fold(0.0_f64, f64::max)
        / dual_scale;
    let primal_objective = 0.5 * dot_product(&candidate.primal, &hessian_product)
        + dot_product(&form.linear_objective, &candidate.primal);
    let dual_objective = -0.5 * dot_product(&candidate.primal, &hessian_product)
        - dot_product(&form.constraint_rhs, &candidate.dual);
    let objective_scale = 1.0 + primal_objective.abs().max(dual_objective.abs());
    ConvexResidualEvidence {
        primal,
        dual,
        stationarity,
        complementarity: dot_product(&candidate.slack, &candidate.dual).abs() / objective_scale,
        relative_gap: (primal_objective - dual_objective).abs() / objective_scale,
    }
}

fn backend_internal_scaling_round_trip_error(candidate: &ClarabelCandidateEnvelope) -> f64 {
    let scaling = &candidate.attempt.internal_scaling;
    let variable_error = if scaling.variable.len() == scaling.inverse_variable.len() {
        scaling
            .variable
            .iter()
            .zip(&scaling.inverse_variable)
            .map(|(forward, inverse)| (forward * inverse - 1.0).abs())
            .fold(0.0_f64, f64::max)
    } else {
        f64::INFINITY
    };
    let constraint_error = if scaling.constraint.len() == scaling.inverse_constraint.len() {
        scaling
            .constraint
            .iter()
            .zip(&scaling.inverse_constraint)
            .map(|(forward, inverse)| (forward * inverse - 1.0).abs())
            .fold(0.0_f64, f64::max)
    } else {
        f64::INFINITY
    };
    let finite = scaling
        .variable
        .iter()
        .chain(&scaling.inverse_variable)
        .chain(&scaling.constraint)
        .chain(&scaling.inverse_constraint)
        .chain([&scaling.objective])
        .all(|value| value.is_finite() && *value > 0.0);
    if finite {
        variable_error.max(constraint_error)
    } else {
        f64::INFINITY
    }
}

#[allow(clippy::too_many_arguments)]
fn recover_and_verify_qp(
    plan: CubicAlgebraicPlan,
    representation: CubicRepresentation,
    problem: CubicCanonicalProblem,
    form: ConvexQpForm,
    scaled: ScaledQpForm,
    candidate: ClarabelCandidateEnvelope,
    backend_residuals: ConvexResidualEvidence,
    attempts: Vec<QpAttemptRecord>,
    accepted_attempt: usize,
    backend_internal_scaling_round_trip_error: f64,
    fault: Option<QpFaultInjection>,
) -> Result<CubicExecutionSolution, CubicExecutionFailure> {
    if !representation
        .solve_coordinate_transform()
        .is_valid_recovery_map()
    {
        return recovery_failure(
            vec![QpRecoveryFailureReason::InvalidRecoveryMap],
            f64::INFINITY,
            f64::INFINITY,
            scaled.round_trip_error,
            attempts,
        );
    }
    let provenance_verified = verifies_qp_form_provenance(&representation, &problem, &form)?;
    if !provenance_verified {
        return recovery_failure(
            vec![QpRecoveryFailureReason::ProvenanceMismatch],
            f64::INFINITY,
            f64::INFINITY,
            scaled.round_trip_error,
            attempts,
        );
    }
    let variable_factors = scaled.scaling.variable_factors();
    let constraint_factors = scaled.scaling.constraint_factors();
    if candidate.primal.len() != form.variables
        || candidate.slack.len() != form.constraint_rhs.len()
        || variable_factors.len() != form.variables
        || constraint_factors.len() != form.constraint_rhs.len()
    {
        return recovery_failure(
            vec![QpRecoveryFailureReason::NonFiniteRecoveredQuantity],
            f64::INFINITY,
            f64::INFINITY,
            scaled.round_trip_error,
            attempts,
        );
    }
    let unscaled_candidate = candidate
        .primal
        .iter()
        .zip(&variable_factors)
        .map(|(value, factor)| value * factor)
        .collect::<Vec<_>>();
    let mut unscaled_slacks = candidate
        .slack
        .iter()
        .zip(&constraint_factors)
        .map(|(value, factor)| value / factor)
        .collect::<Vec<_>>();
    let candidate_scaling_round_trip = candidate
        .primal
        .iter()
        .zip(&unscaled_candidate)
        .zip(&variable_factors)
        .map(|((scaled, unscaled), factor)| relative_error(unscaled / factor, *scaled))
        .chain(
            candidate
                .slack
                .iter()
                .zip(&unscaled_slacks)
                .zip(&constraint_factors)
                .map(|((scaled, unscaled), factor)| relative_error(unscaled * factor, *scaled)),
        )
        .fold(scaled.round_trip_error, f64::max);
    if matches!(fault, Some(QpFaultInjection::BackendStandardForm)) {
        unscaled_slacks[0] += 1.0;
    }
    let reduced = &unscaled_candidate[..form.reduced_field_variables];
    let standard_coefficients = representation
        .null_space()
        .expand(reduced)
        .map_err(|failure| CubicExecutionFailure::Representation(Box::new(failure)))?;
    let mut recovered_reduced = representation
        .null_space()
        .project(&standard_coefficients)
        .map_err(|failure| CubicExecutionFailure::Representation(Box::new(failure)))?;
    if matches!(fault, Some(QpFaultInjection::RecoveryMap)) {
        recovered_reduced[0] += 1.0e-2;
    }
    let reduction_round_trip_error = relative_slice_error(&recovered_reduced, reduced);
    let polynomial_offset = form.reduced_field_variables;
    let standard_polynomial: [f64; POLYNOMIAL_DIMENSION] = unscaled_candidate
        [polynomial_offset..polynomial_offset + POLYNOMIAL_DIMENSION]
        .try_into()
        .expect("the QP Algebraic Plan retains complete Cubic Pi1 coefficients");
    let field = RecoveredCubicField::from_standard_candidate(
        &representation,
        &standard_coefficients,
        standard_polynomial,
    );
    let latent_offset = polynomial_offset + POLYNOMIAL_DIMENSION;
    let latent_values =
        unscaled_candidate[latent_offset..latent_offset + form.semantic_latents].to_vec();
    let semantic_latents = problem
        .semantic_latents
        .iter()
        .zip(&latent_values)
        .map(|(definition, value)| RecoveredSemanticLatent {
            group_id: definition.group_id.clone(),
            field_unit: definition.field_unit.clone(),
            member_source_ids: definition.member_source_ids.clone(),
            value: *value,
        })
        .collect::<Vec<_>>();
    let hard_equalities = problem
        .equalities
        .iter()
        .map(|equality| {
            let value = equality.evaluate(&field, &latent_values);
            RecoveredHardEquality {
                provenance: equality.provenance().clone(),
                dimension: equality.dimension(),
                target: equality.target(),
                value,
                residual: value - equality.target(),
            }
        })
        .collect::<Vec<_>>();
    let soft_equalities = problem
        .soft_equalities
        .iter()
        .map(|relation| {
            let value = relation.evaluate(&field);
            RecoveredSoftEquality {
                provenance: relation.provenance().clone(),
                dimension: relation.dimension(),
                target: relation.target(),
                value,
                residual: value - relation.target(),
            }
        })
        .collect::<Vec<_>>();
    let soft_objectives = problem
        .soft_objectives
        .iter()
        .zip(&form.soft_objective_blocks)
        .map(|(objective, block)| {
            let residual = block
                .canonical_indices
                .iter()
                .map(|index| soft_equalities[*index].residual)
                .collect::<Vec<_>>();
            let dimension = residual.len();
            let whitened_residual =
                dense_matrix_vector_product(&block.whitening, dimension, dimension, &residual);
            let recovered_residual = dense_matrix_vector_product(
                &block.inverse_whitening,
                dimension,
                dimension,
                &whitened_residual,
            );
            RecoveredSoftObjective {
                canonical_indices: block.canonical_indices.clone(),
                loss: objective.loss().clone(),
                covariance_group: objective.covariance_group().cloned(),
                block_kind: objective.block_kind().clone(),
                objective_contribution: 0.5
                    * whitened_residual
                        .iter()
                        .map(|component| component * component)
                        .sum::<f64>(),
                whitening_round_trip_error: relative_slice_error(&recovered_residual, &residual),
                whitened_residual,
            }
        })
        .collect::<Vec<_>>();

    let standard_side_components = std::array::from_fn(|column| {
        (0..standard_coefficients.len())
            .map(|row| {
                representation.polynomial_pairing().get(row, column) * standard_coefficients[row]
            })
            .sum::<f64>()
    });
    let mapped_physical_side = representation
        .solve_coordinate_transform()
        .to_physical_side_condition(standard_side_components);
    let physical_side_components = std::array::from_fn(|column| {
        representation
            .fitting_uses()
            .iter()
            .zip(field.coefficients())
            .map(|(usage, coefficient)| {
                usage.functional().evaluate_affine(
                    if column == 0 { 1.0 } else { 0.0 },
                    std::array::from_fn(|axis| if axis + 1 == column { 1.0 } else { 0.0 }),
                ) * coefficient
            })
            .sum::<f64>()
    });
    let recovered_standard_side = representation
        .solve_coordinate_transform()
        .to_standard_side_condition(physical_side_components);
    let side_condition = PhysicalSideConditionEvidence {
        components: physical_side_components,
        physical_tolerances: representation
            .solve_coordinate_transform()
            .to_physical_side_condition_tolerances(
                [EQUALITY_KKT_POLICY_V1.side_condition_limit; POLYNOMIAL_DIMENSION],
            ),
        standard_components: standard_side_components,
        recovered_standard_components: recovered_standard_side,
        round_trip_error: relative_slice_error(&mapped_physical_side, &physical_side_components)
            .max(relative_slice_error(
                &recovered_standard_side,
                &standard_side_components,
            )),
    };

    let native_cubic_energy = field.native_cubic_energy();
    let field_energy = representation.field_energy_normalization().factor() * native_cubic_energy;
    let standard_energy = dot_product(
        &standard_coefficients,
        &representation
            .kernel_pairing()
            .multiply_vector(&standard_coefficients),
    );
    let recovered_energy = representation.field_energy_normalization().factor() * standard_energy
        / representation.solve_coordinate_transform().length().powi(3);
    let field_energy_round_trip_error = relative_error(field_energy, recovered_energy);
    let field_scale = canonical_field_scale(&representation, &problem, field_energy);
    let hard_relation_tolerances = problem
        .equalities
        .iter()
        .map(|relation| {
            relation_tolerance(
                &representation,
                relation.dimension(),
                relation.target(),
                relation.constant_shift_response(),
                field_scale,
                canonical_gauge_offset(&problem, relation.dimension()),
            )
        })
        .collect::<Vec<_>>();
    let affine_relation_tolerances = problem
        .affine_inequalities
        .iter()
        .map(|relation| {
            relation_tolerance(
                &representation,
                relation.dimension(),
                relation.bound(),
                relation.constant_shift_response(),
                field_scale,
                canonical_gauge_offset(&problem, relation.dimension()),
            )
        })
        .collect::<Vec<_>>();
    let affine_inequalities = problem
        .affine_inequalities
        .iter()
        .zip(&form.affine_inequality_rows)
        .zip(&affine_relation_tolerances)
        .map(|((relation, row), tolerance)| {
            let value = relation.evaluate(&field, &latent_values);
            let margin = relation.physical_margin(value);
            let recovered_violation_channel = row
                .violation_variable
                .map(|variable| unscaled_candidate[variable]);
            let field_violation = (-margin).max(0.0);
            let violation_loss = relation.violation_channel().map(|channel| channel.loss());
            RecoveredAffineInequality {
                provenance: relation.provenance().clone(),
                dimension: relation.dimension(),
                sense: relation.sense(),
                bound: relation.bound(),
                value,
                slack: margin.max(0.0),
                tolerance: tolerance.physical_tolerance,
                violation: recovered_violation_channel.unwrap_or(field_violation),
                recovered_violation_channel,
                violation_loss,
                objective_contribution: violation_loss
                    .zip(recovered_violation_channel)
                    .map(|(loss, recovered)| loss.objective_contribution(recovered)),
                backend_slack: margin + recovered_violation_channel.unwrap_or(0.0),
            }
        })
        .collect::<Vec<_>>();
    let hard_equality_violations = FunctionalViolationEnvelope::from_dimensioned_residuals(
        hard_equalities
            .iter()
            .map(|equality| (equality.dimension, equality.residual)),
    );

    let physical_standard_form_violation =
        standard_form_violation(&form, &unscaled_candidate, &unscaled_slacks);
    if physical_standard_form_violation > EQUALITY_KKT_POLICY_V1.convex_backend_residual_limit {
        return Err(CubicExecutionFailure::BackendContract {
            attempts,
            observed: physical_standard_form_violation,
            limit: EQUALITY_KKT_POLICY_V1.convex_backend_residual_limit,
        });
    }
    let backend_slack_mismatch = affine_inequalities
        .iter()
        .enumerate()
        .map(|(index, inequality)| {
            let backend_row = form.equality_constraints + index;
            (unscaled_slacks[backend_row] - inequality.backend_slack).abs()
                / inequality.backend_slack.abs().max(1.0)
        })
        .chain(
            affine_inequalities
                .iter()
                .filter_map(|relation| relation.recovered_violation_channel)
                .enumerate()
                .map(|(index, violation)| {
                    let backend_row =
                        form.equality_constraints + form.affine_inequality_rows.len() + index;
                    (unscaled_slacks[backend_row] - violation).abs() / violation.abs().max(1.0)
                }),
        )
        .fold(0.0_f64, f64::max);
    let recovered_standard_polynomial = representation
        .solve_coordinate_transform()
        .to_standard(field.physical_polynomial());
    let polynomial_round_trip_error =
        relative_slice_error(&recovered_standard_polynomial, &standard_polynomial);
    let recovered_standard_coefficients = representation
        .solve_coordinate_transform()
        .to_standard_field_coefficients(field.coefficients());
    let field_coefficient_round_trip_error =
        relative_slice_error(&recovered_standard_coefficients, &standard_coefficients);
    let whitening_round_trip_error = soft_objectives
        .iter()
        .map(|objective| objective.whitening_round_trip_error)
        .fold(0.0_f64, f64::max);
    let soft_loss = soft_objectives
        .iter()
        .map(|objective| objective.objective_contribution)
        .sum::<f64>();
    let violation_loss = affine_inequalities
        .iter()
        .filter_map(|relation| relation.objective_contribution)
        .sum::<f64>();
    let total_objective = 0.5 * field_energy + soft_loss + violation_loss;
    let standard_soft_loss = form
        .soft_objective_blocks
        .iter()
        .map(|objective| {
            let residual = objective
                .rows
                .iter()
                .zip(&objective.targets)
                .map(|(row, target)| dot_product(row, &unscaled_candidate) - target)
                .collect::<Vec<_>>();
            let weighted = dense_matrix_vector_product(
                &objective.precision,
                residual.len(),
                residual.len(),
                &residual,
            );
            0.5 * dot_product(&residual, &weighted)
        })
        .sum::<f64>();
    let standard_violation_loss = problem
        .affine_inequalities
        .iter()
        .zip(&form.affine_inequality_rows)
        .filter_map(|(relation, row)| {
            relation
                .violation_channel()
                .zip(row.violation_variable)
                .map(|(channel, variable)| {
                    channel
                        .loss()
                        .objective_contribution(unscaled_candidate[variable])
                })
        })
        .sum::<f64>();
    let mut standard_total_objective =
        0.5 * recovered_energy + standard_soft_loss + standard_violation_loss;
    if matches!(fault, Some(QpFaultInjection::Objective)) {
        standard_total_objective += 1.0;
    }
    let objective_round_trip_error = relative_error(total_objective, standard_total_objective);

    let recovery_finite = standard_coefficients
        .iter()
        .copied()
        .chain(field.coefficients().iter().copied())
        .chain(field.physical_polynomial())
        .chain(latent_values.iter().copied())
        .chain(
            hard_equalities
                .iter()
                .flat_map(|equality| [equality.value, equality.residual]),
        )
        .chain(affine_inequalities.iter().flat_map(|relation| {
            [
                relation.value,
                relation.slack,
                relation.violation,
                relation.recovered_violation_channel.unwrap_or(0.0),
                relation.objective_contribution.unwrap_or(0.0),
            ]
        }))
        .chain(
            soft_equalities
                .iter()
                .flat_map(|relation| [relation.value, relation.residual]),
        )
        .chain([
            field_energy,
            total_objective,
            physical_standard_form_violation,
            backend_slack_mismatch,
            reduction_round_trip_error,
            candidate_scaling_round_trip,
            polynomial_round_trip_error,
            field_coefficient_round_trip_error,
            field_energy_round_trip_error,
            whitening_round_trip_error,
            objective_round_trip_error,
        ])
        .all(f64::is_finite);
    let mut reasons = Vec::new();
    if !recovery_finite {
        reasons.push(QpRecoveryFailureReason::NonFiniteRecoveredQuantity);
    }
    if !side_condition.is_within_policy() {
        reasons.push(QpRecoveryFailureReason::SideConditionViolation);
    }
    if side_condition.round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::SideConditionRoundTripViolation);
    }
    if hard_equalities
        .iter()
        .zip(&hard_relation_tolerances)
        .any(|(relation, tolerance)| relation.residual.abs() > tolerance.physical_tolerance)
    {
        reasons.push(QpRecoveryFailureReason::HardEqualityViolation);
    }
    if affine_inequalities.iter().any(|relation| {
        relation.violation_loss.is_none() && relation.violation > relation.tolerance
            || relation
                .recovered_violation_channel
                .is_some_and(|recovered| {
                    recovered < -relation.tolerance || relation.backend_slack < -relation.tolerance
                })
    }) {
        reasons.push(QpRecoveryFailureReason::AffineInequalityViolation);
    }
    // The slack identity is a backend feasibility equation rather than a
    // coordinate-recovery round trip.  Judge its solve residual with the
    // backend residual policy; the independent physical relation check above
    // still uses the tighter relation tolerance.
    if backend_slack_mismatch > EQUALITY_KKT_POLICY_V1.convex_backend_residual_limit {
        reasons.push(QpRecoveryFailureReason::BackendSlackMismatch);
    }
    if reduction_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::ReductionRoundTripViolation);
    }
    if candidate_scaling_round_trip > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::ScalingRoundTripViolation);
    }
    if polynomial_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::PolynomialRoundTripViolation);
    }
    if field_coefficient_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::FieldCoefficientRoundTripViolation);
    }
    if field_energy_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::FieldEnergyRoundTripViolation);
    }
    if whitening_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::WhiteningRoundTripViolation);
    }
    if objective_round_trip_error > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit {
        reasons.push(QpRecoveryFailureReason::ObjectiveRoundTripViolation);
    }
    if !reasons.is_empty() {
        return recovery_failure(
            reasons,
            physical_standard_form_violation,
            reduction_round_trip_error,
            candidate_scaling_round_trip,
            attempts,
        );
    }

    Ok(CubicExecutionSolution {
        plan,
        representation: representation.evidence().clone(),
        field,
        semantic_latents,
        hard_equalities,
        affine_inequalities,
        soft_equalities,
        soft_objectives,
        hard_equality_violations,
        side_condition,
        field_energy,
        total_objective,
        backend_standard_form_verified: true,
        canonical_acceptance_verified: true,
        qp: Some(CubicQpEvidence {
            capacity: form.capacity,
            scaling: scaled.scaling,
            scaling_round_trip_error: candidate_scaling_round_trip,
            backend_residuals,
            backend_internal_scaling_round_trip_error,
            attempts,
            attempt_plan: QpAttemptPlan::georbf_v1(),
            accepted_attempt,
            physical_standard_form_violation,
            reduction_round_trip_error,
            hard_relation_tolerances,
            affine_relation_tolerances,
            polynomial_round_trip_error,
            field_coefficient_round_trip_error,
            field_energy_round_trip_error,
            whitening_round_trip_error,
            objective_round_trip_error,
            provenance_verified,
        }),
    })
}

fn verifies_qp_form_provenance(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    form: &ConvexQpForm,
) -> Result<bool, CubicExecutionFailure> {
    if form.polynomial_variables != POLYNOMIAL_DIMENSION
        || form.semantic_latents != problem.semantic_latents.len()
        || form.hard_equality_rows.len() != form.equality_constraints
        || form.affine_inequality_rows.len() != problem.affine_inequalities.len()
        || form.affine_inequality_rows.len() + form.violation_nonnegative_rows.len()
            != form.inequality_constraints
        || form.violation_nonnegative_rows.len() != form.soft_violation_variables
        || form.soft_objective_blocks.len() != problem.soft_objectives.len()
    {
        return Ok(false);
    }
    let solver_equalities = problem
        .equalities
        .iter()
        .enumerate()
        .filter(|(_, equality)| {
            equality.participation() == CanonicalEqualityParticipation::SolverConstraint
        })
        .collect::<Vec<_>>();
    for (backend_row, (row, (canonical_index, equality))) in form
        .hard_equality_rows
        .iter()
        .zip(solver_equalities)
        .enumerate()
    {
        let expected = reduced_affine_row(
            representation,
            problem,
            equality.field(),
            equality.latent_coefficients(),
        )?;
        if row.canonical_index != canonical_index
            || row.provenance != *equality.provenance()
            || row.derived_block != DerivedBlockId::from_residual(equality.provenance().residual())
            || row.derived_row != DerivedRowId::from_residual(equality.provenance().residual())
            || row.derived_column
                != equality
                    .field()
                    .map(|_| DerivedColumnId::from_residual(equality.provenance().residual()))
            || row.coefficients != expected
            || row.rhs != equality.target()
            || row.violation_variable.is_some()
            || row.provenance_edges
                != vec![QpProvenanceEdge {
                    provenance: equality.provenance().clone(),
                    derived_block: DerivedBlockId::from_residual(equality.provenance().residual()),
                    derived_row: DerivedRowId::from_residual(equality.provenance().residual()),
                    derived_column: equality
                        .field()
                        .map(|_| DerivedColumnId::from_residual(equality.provenance().residual())),
                    backend_row,
                    backend_column: None,
                    cone: QpConeRole::Zero,
                }]
        {
            return Ok(false);
        }
    }
    for (row, (canonical_index, inequality)) in form
        .affine_inequality_rows
        .iter()
        .zip(problem.affine_inequalities.iter().enumerate())
    {
        let mut expected = reduced_affine_row(
            representation,
            problem,
            inequality.field(),
            inequality.latent_coefficients(),
        )?;
        let multiplier = inequality.sense().upper_form_multiplier();
        expected.iter_mut().for_each(|entry| *entry *= multiplier);
        let expected_violation_variable = inequality.violation_channel().map(|_| {
            form.variables - form.soft_violation_variables
                + problem.affine_inequalities[..canonical_index]
                    .iter()
                    .filter(|relation| relation.violation_channel().is_some())
                    .count()
        });
        if let Some(variable) = expected_violation_variable {
            expected[variable] = -1.0;
        }
        let expected_provenance_edges = inequality
            .source_provenances()
            .iter()
            .map(|provenance| QpProvenanceEdge {
                provenance: provenance.clone(),
                derived_block: DerivedBlockId::from_residual(provenance.residual()),
                derived_row: DerivedRowId::from_residual(provenance.residual()),
                derived_column: inequality
                    .field()
                    .map(|_| DerivedColumnId::from_residual(provenance.residual())),
                backend_row: form.hard_equality_rows.len() + canonical_index,
                backend_column: expected_violation_variable,
                cone: QpConeRole::Nonnegative,
            })
            .collect::<Vec<_>>();
        if row.canonical_index != canonical_index
            || row.provenance != *inequality.provenance()
            || row.derived_block
                != DerivedBlockId::from_residual(inequality.provenance().residual())
            || row.derived_row != DerivedRowId::from_residual(inequality.provenance().residual())
            || row.derived_column
                != inequality
                    .field()
                    .map(|_| DerivedColumnId::from_residual(inequality.provenance().residual()))
            || row.coefficients != expected
            || row.rhs != inequality.upper_form_bound()
            || row.violation_variable != expected_violation_variable
            || row.provenance_edges != expected_provenance_edges
        {
            return Ok(false);
        }
    }
    for (soft_index, (row, (canonical_index, inequality))) in form
        .violation_nonnegative_rows
        .iter()
        .zip(
            problem
                .affine_inequalities
                .iter()
                .enumerate()
                .filter(|(_, inequality)| inequality.violation_channel().is_some()),
        )
        .enumerate()
    {
        let expected_variable = form.variables - form.soft_violation_variables + soft_index;
        let provenance = inequality.provenance();
        let mut expected_coefficients = vec![0.0; form.variables];
        expected_coefficients[expected_variable] = -1.0;
        if inequality.source_provenances().len() != 1
            || inequality.source_provenances().first() != Some(provenance)
            || row.coefficients != expected_coefficients
            || row.rhs != 0.0
            || row.violation_variable != expected_variable
            || row.provenance_edge
                != (QpProvenanceEdge {
                    provenance: provenance.clone(),
                    derived_block: DerivedBlockId::from_residual(provenance.residual()),
                    derived_row: DerivedRowId::from_residual(provenance.residual()),
                    derived_column: Some(DerivedColumnId::from_residual(provenance.residual())),
                    backend_row: form.hard_equality_rows.len()
                        + form.affine_inequality_rows.len()
                        + soft_index,
                    backend_column: Some(expected_variable),
                    cone: QpConeRole::Nonnegative,
                })
            || canonical_index >= problem.affine_inequalities.len()
        {
            return Ok(false);
        }
    }
    for (block, objective) in form
        .soft_objective_blocks
        .iter()
        .zip(&problem.soft_objectives)
    {
        let dimension = objective.residuals().len();
        if block.objective_index >= problem.soft_objectives.len()
            || block.canonical_indices.len() != dimension
            || block.provenances.len() != dimension
            || block.derived_blocks.len() != dimension
            || block.derived_rows.len() != dimension
            || block.derived_columns.len() != dimension
            || block.residuals != objective.residuals()
            || block.targets.len() != dimension
            || block.rows.len() != dimension
            || block.precision != objective.loss().precision_matrix(dimension)
            || block.whitening != objective.loss().whitening_matrix(dimension)
            || block.inverse_whitening != objective.loss().inverse_whitening_matrix(dimension)
            || block.covariance_group.as_ref() != objective.covariance_group()
            || block.block_kind != *objective.block_kind()
        {
            return Ok(false);
        }
        for (component, canonical_index) in block.canonical_indices.iter().enumerate() {
            let Some(relation) = problem.soft_equalities.get(*canonical_index) else {
                return Ok(false);
            };
            let expected =
                reduced_affine_row(representation, problem, Some(relation.field()), &[])?;
            if block.provenances[component] != *relation.provenance()
                || block.residuals[component] != *relation.provenance().residual()
                || block.derived_blocks[component]
                    != DerivedBlockId::from_residual(relation.provenance().residual())
                || block.derived_rows[component]
                    != DerivedRowId::from_residual(relation.provenance().residual())
                || block.derived_columns[component]
                    != DerivedColumnId::from_residual(relation.provenance().residual())
                || block.targets[component] != relation.target()
                || block.rows[component] != expected
            {
                return Ok(false);
            }
        }
    }
    let expected_constraints = form
        .hard_equality_rows
        .iter()
        .chain(&form.affine_inequality_rows)
        .flat_map(|row| row.coefficients.iter().copied())
        .chain(
            form.violation_nonnegative_rows
                .iter()
                .flat_map(|row| row.coefficients.iter().copied()),
        )
        .collect::<Vec<_>>();
    let expected_rhs = form
        .hard_equality_rows
        .iter()
        .chain(&form.affine_inequality_rows)
        .map(|row| row.rhs)
        .chain(form.violation_nonnegative_rows.iter().map(|row| row.rhs))
        .collect::<Vec<_>>();
    Ok(form.constraints == expected_constraints && form.constraint_rhs == expected_rhs)
}

fn canonical_field_scale(
    representation: &CubicRepresentation,
    problem: &CubicCanonicalProblem,
    field_energy: f64,
) -> f64 {
    let length = representation.solve_coordinate_transform().length();
    canonical_characteristic_field_scale(problem, length, field_energy)
}

fn relation_tolerance(
    representation: &CubicRepresentation,
    dimension: FunctionalDimension,
    target: f64,
    constant_shift_response: f64,
    field_scale: f64,
    gauge_offset: f64,
) -> CanonicalRelationToleranceEvidence {
    let characteristic_scale = match dimension {
        FunctionalDimension::FieldValue => field_scale,
        FunctionalDimension::FieldValuePerLength => {
            field_scale / representation.solve_coordinate_transform().length()
        }
    };
    let relation_reference_scale = (target - constant_shift_response * gauge_offset).abs();
    let physical_tolerance = EQUALITY_KKT_POLICY_V1.canonical_characteristic_tolerance_multiplier
        * characteristic_scale
        + EQUALITY_KKT_POLICY_V1.canonical_relation_reference_tolerance_multiplier
            * relation_reference_scale;
    CanonicalRelationToleranceEvidence {
        dimension,
        characteristic_scale,
        relation_reference_scale,
        physical_tolerance,
        standard_tolerance: physical_tolerance,
        scaled_kkt_tolerance: None,
        recovered_physical_tolerance: physical_tolerance,
        round_trip_error: 0.0,
    }
}

fn standard_form_violation(form: &ConvexQpForm, primal: &[f64], slacks: &[f64]) -> f64 {
    let rows = form.constraint_rhs.len();
    let affine = dense_matrix_vector_product(&form.constraints, rows, form.variables, primal);
    let scale = affine
        .iter()
        .chain(slacks)
        .chain(&form.constraint_rhs)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    affine
        .iter()
        .zip(slacks)
        .zip(&form.constraint_rhs)
        .map(|((affine, slack), rhs)| (affine + slack - rhs).abs() / scale)
        .fold(0.0_f64, f64::max)
}

fn recovery_failure<T>(
    reasons: Vec<QpRecoveryFailureReason>,
    backend_standard_form_violation: f64,
    reduction_round_trip_error: f64,
    scaling_round_trip_error: f64,
    attempts: Vec<QpAttemptRecord>,
) -> Result<T, CubicExecutionFailure> {
    Err(CubicExecutionFailure::RecoveryVerification {
        evidence: Box::new(QpRecoveryFailureEvidence {
            reasons,
            backend_standard_form_violation,
            reduction_round_trip_error,
            scaling_round_trip_error,
            no_model_produced: true,
        }),
        attempts,
    })
}

fn relative_error(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic::GlobalAnisotropyMetric;
    use crate::cubic_equality::{
        CanonicalAffineInequality, CanonicalEqualityParticipation, CanonicalHardEquality,
        CanonicalSoftEquality, CanonicalSoftLoss, CanonicalSoftObjective,
        CanonicalViolationChannel, CanonicalViolationLoss, CubicCanonicalProblem,
        SemanticLatentDefinition,
    };
    use crate::functional::{
        CanonicalFunctional, FunctionalDimension, FunctionalTerm, FunctionalUse, RelationId,
        ResidualId, SemanticRolePath, SourceId, UsageProvenance,
    };
    use crate::geometry::FieldUnitLabel;
    use crate::kernel::FieldEnergyNormalization;

    const MANUFACTURED_TARGETS: [f64; 10] = [
        0.204_930_731_175_453_18,
        0.204_930_731_175_453_18,
        -1.670_953_367_752_023_9,
        -2.060_930_563_063_515_7,
        3.156_001_299_339_634_2,
        1.033_400_763_358_782_7,
        0.613_739_534_748_635_2,
        0.764_288_970_221_123_3,
        0.840_004_985_206_390_4,
        0.891_736_553_208_148_7,
    ];

    fn usage(name: &str, support: [f64; 3]) -> FunctionalUse {
        FunctionalUse::new(
            CanonicalFunctional::new(
                FunctionalDimension::FieldValue,
                vec![FunctionalTerm::new(support, 1.0, [0.0; 3])],
            )
            .unwrap(),
            UsageProvenance::new(
                SourceId::new(name),
                None,
                RelationId::new(format!("{name}-relation")),
                ResidualId::new(format!("{name}-residual")),
                SemanticRolePath::new(format!("{name}/value")),
            ),
        )
    }

    fn equality(usage: FunctionalUse, target: f64) -> CanonicalHardEquality {
        CanonicalHardEquality::new(
            Some(usage.clone()),
            Vec::new(),
            usage.provenance().clone(),
            FunctionalDimension::FieldValue,
            target,
            CanonicalEqualityParticipation::SolverConstraint,
        )
    }

    fn equality_problem() -> CubicCanonicalProblem {
        CubicCanonicalProblem {
            equalities: vec![equality(usage("value", [0.0; 3]), 0.0)],
            affine_inequalities: Vec::new(),
            soft_equalities: Vec::new(),
            soft_objectives: Vec::new(),
            semantic_latents: Vec::new(),
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
        }
    }

    fn manufactured_usage(
        name: &str,
        support: [f64; 3],
        value_coefficient: f64,
        gradient_coefficient: [f64; 3],
    ) -> FunctionalUse {
        FunctionalUse::new(
            CanonicalFunctional::new(
                FunctionalDimension::FieldValue,
                vec![FunctionalTerm::new(
                    support,
                    value_coefficient,
                    gradient_coefficient,
                )],
            )
            .unwrap(),
            UsageProvenance::new(
                SourceId::new(name),
                None,
                RelationId::new(format!("{name}-relation")),
                ResidualId::new(format!("{name}-residual")),
                SemanticRolePath::new(format!("{name}/generalized-value")),
            ),
        )
    }

    fn manufactured_problem() -> CubicCanonicalProblem {
        let uses = [
            manufactured_usage("m0", [-1.0, -1.0, -1.0], 1.0, [0.0; 3]),
            manufactured_usage("m1", [1.0, -1.0, -1.0], 1.0, [0.0; 3]),
            manufactured_usage("m2", [-1.0, 1.0, -1.0], 1.0, [0.0; 3]),
            manufactured_usage("m3", [-1.0, -1.0, 1.0], 1.0, [0.0; 3]),
            manufactured_usage("m4", [1.0, 1.0, 1.0], 1.0, [0.0; 3]),
            manufactured_usage("m5", [0.25, -0.5, 0.75], 0.0, [1.0, 0.0, 0.0]),
            manufactured_usage("m6", [-0.75, 0.25, 0.5], 0.0, [0.0, 1.0, 0.0]),
            manufactured_usage("m7", [0.5, 0.75, -0.25], 0.0, [0.0, 0.0, 1.0]),
            manufactured_usage("m8", [0.0, 0.0, 0.0], 1.0, [0.5, -0.25, 0.125]),
            manufactured_usage("m9", [-0.5, 0.625, -0.75], 0.0, [1.0, 1.0, 1.0]),
        ];
        CubicCanonicalProblem {
            equalities: uses
                .into_iter()
                .zip(MANUFACTURED_TARGETS)
                .map(|(usage, target)| equality(usage, target))
                .collect(),
            affine_inequalities: Vec::new(),
            soft_equalities: Vec::new(),
            soft_objectives: Vec::new(),
            semantic_latents: Vec::new(),
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
        }
    }

    fn bounded_manufactured_problem() -> CubicCanonicalProblem {
        let mut problem = manufactured_problem();
        let equality_use = problem.equalities[8].field().unwrap();
        let bound_use = FunctionalUse::new(
            equality_use.functional().clone(),
            UsageProvenance::new(
                SourceId::new("fault-upper"),
                None,
                RelationId::new("fault-upper-relation"),
                ResidualId::new("fault-upper-residual"),
                SemanticRolePath::new("affine-upper/value"),
            ),
        );
        problem
            .affine_inequalities
            .push(CanonicalAffineInequality::upper_bound(
                Some(bound_use.clone()),
                Vec::new(),
                bound_use.provenance().clone(),
                FunctionalDimension::FieldValue,
                MANUFACTURED_TARGETS[8] + 1.0,
            ));
        problem
    }

    fn soft_bounded_manufactured_problem() -> CubicCanonicalProblem {
        let mut problem = manufactured_problem();
        let equality_use = problem.equalities[8].field().unwrap();
        let bound_use = FunctionalUse::new(
            equality_use.functional().clone(),
            UsageProvenance::new(
                SourceId::new("soft-upper"),
                None,
                RelationId::new("soft-upper-relation"),
                ResidualId::new("soft-upper-residual"),
                SemanticRolePath::new("affine-upper/soft-value"),
            ),
        );
        problem
            .affine_inequalities
            .push(CanonicalAffineInequality::new(
                Some(bound_use.clone()),
                Vec::new(),
                bound_use.provenance().clone(),
                FunctionalDimension::FieldValue,
                CanonicalInequalitySense::Upper,
                MANUFACTURED_TARGETS[8] - 0.25,
                Some(CanonicalViolationChannel::new(
                    ResidualId::new("soft-upper-residual"),
                    CanonicalViolationLoss::QuadraticPenalty { weight: 2.0 },
                )),
            ));
        problem.field_energy_normalization = FieldEnergyNormalization::try_new(1.0)
            .expect("a soft bound requires positive field-energy normalization");
        problem
    }

    fn add_shared_latent_and_soft_objective(problem: &mut CubicCanonicalProblem) {
        let latent_provenance = UsageProvenance::new(
            SourceId::new("latent-gauge"),
            Some(GroupId::new("latent-group")),
            RelationId::new("latent-gauge-relation"),
            ResidualId::new("latent-gauge-residual"),
            SemanticRolePath::new("latent-group/gauge"),
        );
        problem.semantic_latents.push(SemanticLatentDefinition {
            group_id: GroupId::new("latent-group"),
            field_unit: FieldUnitLabel::new("manufactured-unit"),
            member_source_ids: vec![SourceId::new("latent-member")],
        });
        problem.equalities.push(CanonicalHardEquality::new(
            None,
            vec![SemanticLatentCoefficient {
                latent: 0,
                coefficient: 1.0,
            }],
            latent_provenance,
            FunctionalDimension::FieldValue,
            2.5,
            CanonicalEqualityParticipation::SolverConstraint,
        ));

        let hard_use = problem.equalities[9]
            .field()
            .expect("the manufactured field equality has a functional")
            .clone();
        let soft_use = FunctionalUse::new(
            hard_use.functional().clone(),
            UsageProvenance::new(
                SourceId::new("soft-check"),
                None,
                RelationId::new("soft-check-relation"),
                ResidualId::new("soft-check-residual"),
                SemanticRolePath::new("soft/check"),
            ),
        );
        let soft = CanonicalSoftEquality::new(soft_use, MANUFACTURED_TARGETS[9] + 0.125);
        problem.soft_objectives.push(CanonicalSoftObjective::new(
            soft.provenance().residual().clone(),
            CanonicalSoftLoss::QuadraticPenalty { weight: 2.0 },
        ));
        problem.soft_equalities.push(soft);
        problem.field_energy_normalization = FieldEnergyNormalization::try_new(3.0)
            .expect("the manufactured objective normalization is positive");
    }

    fn flattened_bounded_problem() -> CubicCanonicalProblem {
        let original = manufactured_problem();
        let flattened_equalities = original
            .equalities
            .into_iter()
            .map(|relation| {
                let usage = relation
                    .field()
                    .expect("the manufactured equality has a functional");
                let term = usage.functional().terms()[0];
                let mut support = term.support();
                let mut gradient = term.gradient_coefficient();
                support[2] = 0.0;
                gradient[2] = 0.0;
                if term.value_coefficient() == 0.0 && gradient == [0.0; 3] {
                    gradient[0] = 1.0;
                }
                let flattened = FunctionalUse::new(
                    CanonicalFunctional::new(
                        usage.functional().dimension(),
                        vec![FunctionalTerm::new(
                            support,
                            term.value_coefficient(),
                            gradient,
                        )],
                    )
                    .expect("the flattened functional remains valid"),
                    usage.provenance().clone(),
                );
                equality(flattened, relation.target())
            })
            .collect::<Vec<_>>();
        let bound_use = flattened_equalities[0]
            .field()
            .expect("the flattened equality has a functional")
            .clone();
        CubicCanonicalProblem {
            equalities: flattened_equalities,
            affine_inequalities: vec![CanonicalAffineInequality::upper_bound(
                Some(bound_use.clone()),
                Vec::new(),
                bound_use.provenance().clone(),
                FunctionalDimension::FieldValue,
                MANUFACTURED_TARGETS[0] + 1.0,
            )],
            soft_equalities: Vec::new(),
            soft_objectives: Vec::new(),
            semantic_latents: Vec::new(),
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
        }
    }

    #[test]
    fn canonical_execution_route_depends_only_on_affine_inequality_capability() {
        let equality_only = equality_problem();
        assert_eq!(
            CubicExecutionCore::plan(&equality_only).form,
            CubicFormKind::SymmetricKkt
        );

        let mut bounded = equality_only.clone();
        let bound_use = usage("upper", [1.0, 0.0, 0.0]);
        bounded
            .affine_inequalities
            .push(CanonicalAffineInequality::upper_bound(
                Some(bound_use.clone()),
                Vec::new(),
                bound_use.provenance().clone(),
                FunctionalDimension::FieldValue,
                1.0,
            ));
        let plan = CubicExecutionCore::plan(&bounded);
        assert_eq!(plan.form, CubicFormKind::ConvexQp);
        assert_eq!(plan.hard_equalities, equality_only.equalities.len());
        assert_eq!(plan.affine_inequalities, 1);
        assert_eq!(
            CubicEqualityCore::solve_canonical(bounded, GlobalAnisotropyMetric::identity(),)
                .expect_err("the Equality entry point must not ignore an affine inequality"),
            CubicEqualityFailure::AffineInequalityRequiresConvexQp
        );
    }

    #[test]
    fn soft_bound_capacity_includes_violation_variable_and_nonnegative_row() {
        let problem = soft_bounded_manufactured_problem();
        let base_variables = canonical_fitting_uses(
            &problem.equalities,
            &problem.soft_equalities,
            &problem.affine_inequalities,
        )
        .len()
            + problem.semantic_latents.len();
        let solver_equalities = problem
            .equalities
            .iter()
            .filter(|equality| {
                equality.participation() == CanonicalEqualityParticipation::SolverConstraint
            })
            .count();

        let solution = CubicExecutionCore::solve(problem, GlobalAnisotropyMetric::identity())
            .expect("the soft-bound QP should fit within its exact capacity plan");
        let qp = solution.qp.expect("a soft bound selects the QP route");

        assert_eq!(qp.capacity.variables, base_variables + 1);
        assert_eq!(qp.capacity.constraints, solver_equalities + 2);
        assert_eq!(solution.affine_inequalities.len(), 1);
        assert!(
            solution.affine_inequalities[0]
                .recovered_violation_channel
                .is_some()
        );
    }

    #[test]
    fn redundant_affine_bound_uses_qp_and_recovers_kkt_canonical_observables() {
        let mut equality_problem = manufactured_problem();
        add_shared_latent_and_soft_objective(&mut equality_problem);
        let equality =
            CubicExecutionCore::solve(equality_problem.clone(), GlobalAnisotropyMetric::identity())
                .expect("the equality capability should execute through KKT");

        let mut bounded_problem = equality_problem;
        let equality_use = bounded_problem.equalities[8].field().unwrap();
        let bound_use = FunctionalUse::new(
            equality_use.functional().clone(),
            UsageProvenance::new(
                SourceId::new("redundant-upper"),
                None,
                RelationId::new("redundant-upper-relation"),
                ResidualId::new("redundant-upper-residual"),
                SemanticRolePath::new("affine-upper/value"),
            ),
        );
        bounded_problem
            .affine_inequalities
            .push(CanonicalAffineInequality::upper_bound(
                Some(bound_use.clone()),
                Vec::new(),
                bound_use.provenance().clone(),
                FunctionalDimension::FieldValue,
                MANUFACTURED_TARGETS[8] + 1.0,
            ));
        let qp = CubicExecutionCore::solve(bounded_problem, GlobalAnisotropyMetric::identity())
            .expect("the affine inequality capability should execute through Clarabel QP");

        assert_eq!(equality.plan.form, CubicFormKind::SymmetricKkt);
        assert_eq!(qp.plan.form, CubicFormKind::ConvexQp);
        assert!(qp.backend_standard_form_verified);
        assert!(qp.canonical_acceptance_verified);
        assert_eq!(qp.affine_inequalities.len(), 1);
        assert!(qp.affine_inequalities[0].slack >= 0.0);
        for (left, right) in equality
            .field
            .coefficients()
            .iter()
            .zip(qp.field.coefficients())
        {
            assert!((left - right).abs() <= 1.0e-8);
        }
        for (left, right) in equality
            .field
            .physical_polynomial()
            .into_iter()
            .zip(qp.field.physical_polynomial())
        {
            assert!((left - right).abs() <= 1.0e-8);
        }
        assert!((equality.field_energy - qp.field_energy).abs() <= 1.0e-8);
        assert!((equality.total_objective - qp.total_objective).abs() <= 1.0e-8);
        assert_eq!(equality.semantic_latents.len(), 1);
        assert_eq!(qp.semantic_latents.len(), 1);
        assert!(
            (equality.semantic_latents[0].value - qp.semantic_latents[0].value).abs() <= 1.0e-8
        );
        assert_eq!(equality.soft_equalities.len(), 1);
        assert_eq!(qp.soft_equalities.len(), 1);
        assert_eq!(equality.soft_objectives.len(), 1);
        assert_eq!(qp.soft_objectives.len(), 1);
        assert!(
            (equality.soft_objectives[0].objective_contribution
                - qp.soft_objectives[0].objective_contribution)
                .abs()
                <= 1.0e-8
        );
        let evidence = qp
            .qp
            .as_ref()
            .expect("the affine bound selects the QP route");
        assert_eq!(evidence.scaling.rounds.len(), 8);
        assert!(evidence.scaling.rounds.iter().all(|round| {
            round
                .variable_exponents
                .iter()
                .chain(&round.constraint_exponents)
                .all(|exponent| (-8..=8).contains(exponent))
        }));
        assert!(
            evidence
                .scaling
                .cumulative_variable_exponents
                .iter()
                .chain(&evidence.scaling.cumulative_constraint_exponents)
                .all(|exponent| (-32..=32).contains(exponent))
        );
        assert_eq!(evidence.attempts.len(), 1);
        assert_eq!(evidence.attempts[0].backend.requested_threads, 1);
        assert_eq!(evidence.attempts[0].backend.actual_threads, 1);
        assert_eq!(evidence.attempts[0].backend.linear_solver, "qdldl");
        assert_eq!(
            evidence.attempts[0].backend.backend,
            clarabel_backend::ClarabelBackendFingerprint {
                crate_name: "clarabel",
                crate_version: "0.11.1",
                features: ["serde"],
                direct_solver: "qdldl",
            }
        );
        assert!(
            evidence.attempts[0]
                .backend
                .settings
                .all_settings
                .contains("tol_feas")
        );
        assert_eq!(evidence.attempts[0].georbf_scaling, evidence.scaling);
        assert_eq!(
            evidence.attempts[0].georbf_scaling_round_trip_error,
            evidence.scaling_round_trip_error
        );
        let certificate = evidence.attempts[0].backend.certificate;
        assert!(certificate.primal_infeasibility_residual.is_finite());
        assert!(certificate.dual_infeasibility_residual.is_finite());
        assert!(certificate.kappa_tau_ratio.is_finite());
        assert_eq!(clarabel_backend::CRATE_NAME, "clarabel");
        assert_eq!(clarabel_backend::CRATE_VERSION, "0.11.1");
        assert_eq!(clarabel_backend::FEATURES, ["serde"]);
        assert!(evidence.backend_internal_scaling_round_trip_error <= 1.0e-11);
        for point in [[0.2, -0.3, 0.4], [-0.7, 0.1, 0.9]] {
            let left = equality.field.sample(point);
            let right = qp.field.sample(point);
            assert!((left.value - right.value).abs() <= 1.0e-8);
            for (left, right) in left.gradient.into_iter().zip(right.gradient) {
                assert!((left - right).abs() <= 1.0e-8);
            }
        }
    }

    #[test]
    fn qp_polynomial_rank_defect_is_structured_before_clarabel_entry() {
        let failure = CubicExecutionCore::solve(
            flattened_bounded_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a missing affine mode must stop before the QP backend");

        match failure {
            CubicExecutionFailure::Representation(failure) => match *failure {
                RepresentationFailure::PolynomialRankDeficient { mode, .. } => {
                    assert_eq!(mode.residual, 0.0);
                    assert!(!mode.execution.solver_invoked);
                    assert!(!mode.execution.hidden_regularization_applied);
                }
                other => panic!("unexpected representation failure: {other:?}"),
            },
            other => panic!("unexpected execution failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_qp_provenance_is_a_recovery_failure_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::Provenance);
        let failure = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a reassociated QP row must fail canonical provenance recovery");

        match failure {
            CubicExecutionFailure::RecoveryVerification { evidence, attempts } => {
                assert_eq!(
                    evidence.reasons,
                    vec![QpRecoveryFailureReason::ProvenanceMismatch]
                );
                assert!(evidence.no_model_produced);
                assert!(!attempts.is_empty());
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_qp_scaling_map_is_rejected_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::ScalingMap);
        let failure = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a corrupted GeoRBF scaling map must fail recovery");

        match failure {
            CubicExecutionFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![QpRecoveryFailureReason::ScalingRoundTripViolation]
                );
                assert!(
                    evidence.scaling_round_trip_error
                        > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_qp_reduction_recovery_map_is_rejected_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::RecoveryMap);
        let failure = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a corrupted Householder recovery map must fail recovery");

        match failure {
            CubicExecutionFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![QpRecoveryFailureReason::ReductionRoundTripViolation]
                );
                assert!(
                    evidence.reduction_round_trip_error
                        > EQUALITY_KKT_POLICY_V1.recovery_round_trip_limit
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_qp_backend_residual_is_a_backend_contract_failure() {
        inject_qp_fault_once(QpFaultInjection::BackendResidual);
        let failure = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a corrupted Clarabel residual envelope must fail before recovery");

        match failure {
            CubicExecutionFailure::BackendContract {
                attempts,
                observed,
                limit,
            } => {
                assert_eq!(attempts.len(), 1);
                assert_eq!(
                    attempts[0].failure_reason,
                    Some(QpAttemptFailureReason::BackendResidualExceeded)
                );
                assert!(observed > limit);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn damaged_recovered_standard_form_is_a_backend_contract_failure() {
        inject_qp_fault_once(QpFaultInjection::BackendStandardForm);
        let failure = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a damaged recovered backend equation is not a canonical recovery failure");

        match failure {
            CubicExecutionFailure::BackendContract {
                attempts,
                observed,
                limit,
            } => {
                assert_eq!(attempts.len(), 1);
                assert!(observed > limit);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn derivative_soft_loss_scale_is_converted_to_field_value_units() {
        let mut problem = manufactured_problem();
        let soft_use = FunctionalUse::new(
            CanonicalFunctional::new(
                FunctionalDimension::FieldValuePerLength,
                vec![FunctionalTerm::new([0.0, 0.0, 0.0], 0.0, [1.0, 0.0, 0.0])],
            )
            .expect("the derivative functional is nonzero"),
            UsageProvenance::new(
                SourceId::new("derivative-scale"),
                None,
                RelationId::new("derivative-scale-relation"),
                ResidualId::new("derivative-scale-residual"),
                SemanticRolePath::new("derivative/scale"),
            ),
        );
        let soft = CanonicalSoftEquality::new(soft_use, 0.0);
        problem.soft_objectives.push(CanonicalSoftObjective::new(
            soft.provenance().residual().clone(),
            CanonicalSoftLoss::QuadraticPenalty { weight: 1.0e-12 },
        ));
        problem.soft_equalities.push(soft);
        let representation = CubicRepresentation::new(
            canonical_fitting_uses(
                &problem.equalities,
                &problem.soft_equalities,
                &problem.affine_inequalities,
            ),
            GlobalAnisotropyMetric::identity(),
        )
        .expect("the derivative-scale representation is valid");

        let expected = representation.solve_coordinate_transform().length() * 1.0e6;
        let actual = canonical_field_scale(&representation, &problem, 0.0);
        assert!(
            relative_error(actual, expected) <= 1.0e-12,
            "actual={actual:e}, expected={expected:e}"
        );
    }

    #[test]
    fn derivative_violation_loss_scale_is_converted_to_field_value_units() {
        let mut problem = manufactured_problem();
        let derivative_use = FunctionalUse::new(
            CanonicalFunctional::new(
                FunctionalDimension::FieldValuePerLength,
                vec![FunctionalTerm::new(
                    [0.125, -0.25, 0.5],
                    0.0,
                    [1.0, 0.0, 0.0],
                )],
            )
            .expect("the derivative functional is nonzero"),
            UsageProvenance::new(
                SourceId::new("derivative-violation-scale"),
                None,
                RelationId::new("derivative-violation-scale-relation"),
                ResidualId::new("derivative-violation-scale-residual"),
                SemanticRolePath::new("derivative-violation/scale"),
            ),
        );
        problem
            .affine_inequalities
            .push(CanonicalAffineInequality::new(
                Some(derivative_use.clone()),
                Vec::new(),
                derivative_use.provenance().clone(),
                FunctionalDimension::FieldValuePerLength,
                CanonicalInequalitySense::Upper,
                0.0,
                Some(CanonicalViolationChannel::new(
                    ResidualId::new("derivative-violation-scale-residual"),
                    CanonicalViolationLoss::QuadraticPenalty { weight: 1.0e-12 },
                )),
            ));
        problem.field_energy_normalization = FieldEnergyNormalization::try_new(1.0)
            .expect("a soft violation requires positive field-energy normalization");

        let length = 3.0;
        let expected = length * 1.0e6;
        let actual = canonical_characteristic_field_scale(&problem, length, 0.0);
        assert!(
            relative_error(actual, expected) <= 1.0e-12,
            "actual={actual:e}, expected={expected:e}"
        );
    }

    #[test]
    fn damaged_qp_objective_round_trip_is_rejected_without_a_model() {
        inject_qp_fault_once(QpFaultInjection::Objective);
        let failure = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect_err("a corrupted standard-form objective must fail canonical recovery");

        match failure {
            CubicExecutionFailure::RecoveryVerification { evidence, .. } => {
                assert_eq!(
                    evidence.reasons,
                    vec![QpRecoveryFailureReason::ObjectiveRoundTripViolation]
                );
                assert!(evidence.no_model_produced);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn unverified_standard_attempt_uses_the_single_robust_retry() {
        inject_qp_fault_once(QpFaultInjection::StandardRetry);
        let solution = CubicExecutionCore::solve(
            bounded_manufactured_problem(),
            GlobalAnisotropyMetric::identity(),
        )
        .expect("the deterministic robust attempt should recover the valid QP");
        let qp = solution.qp.expect("the affine bound selects the QP route");

        assert_eq!(
            qp.attempt_plan.profiles,
            [
                ClarabelAttemptProfile::Standard,
                ClarabelAttemptProfile::Robust,
            ]
        );
        assert_eq!(qp.accepted_attempt, 1);
        assert_eq!(qp.attempts.len(), 2);
        assert_eq!(
            qp.attempts[0].failure_reason,
            Some(QpAttemptFailureReason::BackendResidualExceeded)
        );
        assert_eq!(qp.attempts[1].failure_reason, None);
        assert_eq!(
            qp.attempts[0].backend.settings.feasibility_tolerance,
            1.0e-9
        );
        assert_eq!(
            qp.attempts[1].backend.settings.feasibility_tolerance,
            1.0e-10
        );
        assert!(solution.canonical_acceptance_verified);
    }

    #[test]
    fn farkas_validation_rejects_unseparated_nonstationary_and_wrong_cone_rays() {
        // x <= 0 and x >= 1, written in Clarabel's A x + s = b form.
        let form = ScaledQpForm {
            hessian: vec![1.0],
            linear_objective: vec![0.0],
            constraints: vec![1.0, -1.0],
            constraint_rhs: vec![0.0, -1.0],
            equality_constraints: 0,
            inequality_constraints: 2,
            scaling: QpScalingEvidence {
                rounds: Vec::new(),
                cumulative_variable_exponents: vec![0],
                cumulative_constraint_exponents: vec![0, 0],
                saturated_outside_target: 0,
            },
            round_trip_error: 0.0,
        };

        let valid = validate_primal_infeasibility_ray(
            form.linear_objective.len(),
            &form.constraints,
            &form.constraint_rhs,
            form.equality_constraints,
            &[1.0, 1.0],
        )
        .expect("the normalized positive ray proves the two bounds incompatible");
        assert_eq!(valid.normalized_ray_norm, 1.0);
        assert_eq!(valid.stationarity_residual, 0.0);
        assert_eq!(valid.dual_cone_violation, 0.0);
        assert_eq!(valid.separation_margin, 1.0);
        for invalid in [[1.0, 0.0], [-1.0, -1.0], [0.0, 0.0]] {
            assert!(
                validate_primal_infeasibility_ray(
                    form.linear_objective.len(),
                    &form.constraints,
                    &form.constraint_rhs,
                    form.equality_constraints,
                    &invalid,
                )
                .is_none()
            );
        }
    }
}
