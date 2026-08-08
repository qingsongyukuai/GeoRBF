use faer::diag::{DiagMut, DiagRef};
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::lblt::{factor, solve};
use faer::{Conj, MatMut, MatRef};

use crate::capacity::{
    CapacityExceededEvidence, EqualityCapacityPlan, FaerWorkspaceEvidence, plan_equality_capacity,
};
use crate::faer_backend;
use crate::numerical::{
    EQUALITY_KKT_POLICY_V2, NumericalPolicyId, SpectralAnalysisFailure, SpectralRankDecision,
    analyze_spectral_rank,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RuizRoundEvidence {
    pub(crate) exponents: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RuizScalingEvidence {
    pub(crate) rounds: Vec<RuizRoundEvidence>,
    pub(crate) cumulative_exponents: Vec<i32>,
    pub(crate) saturated_outside_target: Vec<usize>,
}

impl RuizScalingEvidence {
    fn factors(&self) -> Vec<f64> {
        self.cumulative_exponents
            .iter()
            .map(|exponent| 2.0_f64.powi(*exponent))
            .collect()
    }

    pub(crate) fn recover_matrix(&self, scaled: &[f64]) -> Vec<f64> {
        let dimension = self.cumulative_exponents.len();
        let factors = self.factors();
        (0..dimension * dimension)
            .map(|index| {
                let row = index % dimension;
                let column = index / dimension;
                scaled[index] / (factors[row] * factors[column])
            })
            .collect()
    }

    pub(crate) fn recover_rhs(&self, scaled: &[f64]) -> Vec<f64> {
        scaled
            .iter()
            .zip(self.factors())
            .map(|(value, factor)| value / factor)
            .collect()
    }

    pub(crate) fn scale_residual_or_tolerance(&self, physical: &[f64]) -> Vec<f64> {
        physical
            .iter()
            .zip(self.factors())
            .map(|(value, factor)| value * factor)
            .collect()
    }

    pub(crate) fn recover_residual_or_tolerance(&self, scaled: &[f64]) -> Vec<f64> {
        self.recover_rhs(scaled)
    }

    fn recover_solution(&self, scaled: &[f64]) -> Vec<f64> {
        scaled
            .iter()
            .zip(self.factors())
            .map(|(value, factor)| value * factor)
            .collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RuizScalingFailure {
    InvalidShape,
    ZeroNorm { index: usize },
    NonFiniteNorm { index: usize },
}

#[derive(Debug, Clone, PartialEq)]
struct ScaledKkt {
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    evidence: RuizScalingEvidence,
}

fn equilibrate_symmetric_kkt(
    matrix: &[f64],
    rhs: &[f64],
    dimension: usize,
) -> Result<ScaledKkt, RuizScalingFailure> {
    if matrix.len() != dimension.saturating_mul(dimension) || rhs.len() != dimension {
        return Err(RuizScalingFailure::InvalidShape);
    }
    let mut matrix = matrix.to_vec();
    let mut rhs = rhs.to_vec();
    let mut cumulative_exponents = vec![0_i32; dimension];
    let mut rounds = Vec::with_capacity(EQUALITY_KKT_POLICY_V2.ruiz_rounds);
    for _ in 0..EQUALITY_KKT_POLICY_V2.ruiz_rounds {
        let mut exponents = Vec::with_capacity(dimension);
        for row in 0..dimension {
            let norm = (0..dimension)
                .map(|column| matrix[row + column * dimension].abs())
                .fold(0.0_f64, f64::max);
            if norm == 0.0 {
                return Err(RuizScalingFailure::ZeroNorm { index: row });
            }
            if !norm.is_finite() {
                return Err(RuizScalingFailure::NonFiniteNorm { index: row });
            }
            let desired = (-0.5 * norm.log2()).round();
            let desired = desired.clamp(i32::MIN as f64, i32::MAX as f64) as i32;
            let round_limit = EQUALITY_KKT_POLICY_V2.ruiz_single_round_exponent_limit;
            let cumulative_limit = EQUALITY_KKT_POLICY_V2.ruiz_cumulative_exponent_limit;
            let exponent = desired.clamp(-round_limit, round_limit).clamp(
                -cumulative_limit - cumulative_exponents[row],
                cumulative_limit - cumulative_exponents[row],
            );
            exponents.push(exponent);
        }
        let factors = exponents
            .iter()
            .map(|exponent| 2.0_f64.powi(*exponent))
            .collect::<Vec<_>>();
        for column in 0..dimension {
            for row in 0..dimension {
                matrix[row + column * dimension] *= factors[row] * factors[column];
            }
        }
        for row in 0..dimension {
            rhs[row] *= factors[row];
            cumulative_exponents[row] += exponents[row];
        }
        rounds.push(RuizRoundEvidence { exponents });
    }
    let saturated_outside_target = (0..dimension)
        .filter(|row| {
            let norm = (0..dimension)
                .map(|column| matrix[*row + column * dimension].abs())
                .fold(0.0_f64, f64::max);
            !(0.5..=2.0).contains(&norm)
                && cumulative_exponents[*row].abs()
                    == EQUALITY_KKT_POLICY_V2.ruiz_cumulative_exponent_limit
        })
        .collect();
    Ok(ScaledKkt {
        matrix,
        rhs,
        evidence: RuizScalingEvidence {
            rounds,
            cumulative_exponents,
            saturated_outside_target,
        },
    })
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct EqualityKktSystem<'a> {
    pub(crate) primal_variables: usize,
    pub(crate) equality_constraints: usize,
    pub(crate) hessian: &'a [f64],
    pub(crate) equality_jacobian: &'a [f64],
    pub(crate) stationarity_rhs: &'a [f64],
    pub(crate) equality_rhs: &'a [f64],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BackendFingerprint {
    pub(crate) schema_version: u32,
    pub(crate) crate_name: &'static str,
    pub(crate) crate_version: &'static str,
    pub(crate) features: [&'static str; 2],
    pub(crate) algorithm: &'static str,
    pub(crate) target_arch: &'static str,
    pub(crate) target_os: &'static str,
    pub(crate) requested_threads: usize,
    pub(crate) actual_threads: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BackendAttemptSettings {
    Lblt {
        pivoting: &'static str,
        block_size: usize,
        parallelism_threshold: usize,
        factor_workspace_source: &'static str,
        maximum_refinement_steps: usize,
    },
    FullSvd {
        settings_id: &'static str,
        left_vectors: &'static str,
        right_vectors: &'static str,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AlgebraicAnalysisSettings {
    pub(crate) llt_settings_id: &'static str,
    pub(crate) householder_qr_settings_id: &'static str,
    pub(crate) rrqr_algorithm: &'static str,
    pub(crate) rrqr_settings_id: &'static str,
    pub(crate) svd_algorithm: &'static str,
    pub(crate) svd_settings_id: &'static str,
    pub(crate) evd_algorithm: &'static str,
    pub(crate) evd_settings_id: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum VerificationScope {
    BackendStandardForm,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RankClassification {
    FullRank,
    RankDeficient,
    NumericalDecisionGrayZone,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AlgebraicRankEvidence {
    pub(crate) exact_zero_index: Option<usize>,
    pub(crate) rrqr_ratio: f64,
    pub(crate) singular_values: Vec<f64>,
    pub(crate) svd_ratio: f64,
    pub(crate) reject_ratio: f64,
    pub(crate) accept_ratio: f64,
    pub(crate) classification: RankClassification,
    pub(crate) backend_invoked: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Inertia {
    pub(crate) positive: usize,
    pub(crate) negative: usize,
    pub(crate) zero: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum KktAttemptKind {
    BunchKaufmanRefinement,
    SvdRescue,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct KktAttemptPlan {
    pub(crate) numerical_policy: NumericalPolicyId,
    pub(crate) attempts: [KktAttemptKind; 2],
    pub(crate) svd_rescue_requires_confirmed_full_rank: bool,
}

impl KktAttemptPlan {
    fn equality_v1() -> Self {
        Self {
            numerical_policy: EQUALITY_KKT_POLICY_V2.id,
            attempts: [
                KktAttemptKind::BunchKaufmanRefinement,
                KktAttemptKind::SvdRescue,
            ],
            svd_rescue_requires_confirmed_full_rank: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SolveAttemptTermination {
    CandidateProduced,
    NumericalError,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum BackendContractViolationReason {
    NonFiniteCandidate,
    BackwardErrorExceeded { observed: f64, limit: f64 },
    ScalingRoundTripExceeded { observed: f64, limit: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NumericalFailureReason {
    BackendDecompositionFailure,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum KktAttemptFailureReason {
    BackendContract(BackendContractViolationReason),
    Numerical(NumericalFailureReason),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LinearResidualEvidence {
    pub(crate) infinity_norm: f64,
    pub(crate) matrix_infinity_norm: f64,
    pub(crate) solution_infinity_norm: f64,
    pub(crate) rhs_infinity_norm: f64,
    pub(crate) normalized_backward_error: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ScalingSummary {
    pub(crate) method: &'static str,
    pub(crate) rounds: usize,
    pub(crate) saturated_outside_target: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct KktAttemptRecord {
    pub(crate) sequence: usize,
    pub(crate) kind: KktAttemptKind,
    pub(crate) backend: BackendFingerprint,
    pub(crate) settings: BackendAttemptSettings,
    pub(crate) scaling: ScalingSummary,
    pub(crate) requested_threads: usize,
    pub(crate) actual_threads: usize,
    pub(crate) refinement_steps: usize,
    pub(crate) termination: SolveAttemptTermination,
    pub(crate) residual: Option<LinearResidualEvidence>,
    pub(crate) certificate_present: bool,
    pub(crate) failure_reason: Option<KktAttemptFailureReason>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct KktSolveEvidence {
    pub(crate) candidate: Vec<f64>,
    pub(crate) equality_multipliers: Vec<f64>,
    pub(crate) normalized_backward_error: f64,
    pub(crate) numerical_policy: NumericalPolicyId,
    pub(crate) verification_scope: VerificationScope,
    pub(crate) capacity: EqualityCapacityPlan,
    pub(crate) workspace: FaerWorkspaceEvidence,
    pub(crate) backend: BackendFingerprint,
    pub(crate) analysis_settings: AlgebraicAnalysisSettings,
    pub(crate) scaling: RuizScalingEvidence,
    pub(crate) scaling_round_trip_error: f64,
    pub(crate) rank: AlgebraicRankEvidence,
    pub(crate) expected_inertia: Inertia,
    pub(crate) observed_inertia: Inertia,
    pub(crate) attempts: Vec<KktAttemptRecord>,
    pub(crate) attempt_plan: KktAttemptPlan,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CompletedKktAnalysis {
    pub(crate) rank: AlgebraicRankEvidence,
    pub(crate) expected_inertia: Inertia,
    pub(crate) observed_inertia: Inertia,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum KktInputField {
    Hessian,
    EqualityJacobian,
    StationarityRightHandSide,
    EqualityRightHandSide,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum WorkspacePhase {
    RankAnalysis,
    InertiaAnalysis,
    Factor,
    Solve,
    SvdRescue,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlgebraicAnalysisPhase {
    RankConfirmation,
    Inertia,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum KktFailure {
    Capacity(CapacityExceededEvidence),
    InvalidLength {
        field: KktInputField,
        expected: usize,
        actual: usize,
    },
    NonFiniteInput {
        field: KktInputField,
        index: usize,
    },
    WorkspaceAllocation {
        phase: WorkspacePhase,
        bytes: u64,
        alignment: usize,
    },
    Scaling(RuizScalingFailure),
    RankDeficient {
        evidence: AlgebraicRankEvidence,
    },
    NumericalDecisionGrayZone {
        evidence: AlgebraicRankEvidence,
    },
    AlgebraicAnalysisFailure {
        phase: AlgebraicAnalysisPhase,
    },
    UnexpectedInertia {
        expected: Inertia,
        observed: Inertia,
        backend_invoked: bool,
    },
    BackendContractViolation {
        attempt_plan: KktAttemptPlan,
        attempts: Vec<KktAttemptRecord>,
        reason: BackendContractViolationReason,
        analysis: Option<Box<CompletedKktAnalysis>>,
    },
    NumericalFailure {
        attempt_plan: KktAttemptPlan,
        attempts: Vec<KktAttemptRecord>,
        reason: NumericalFailureReason,
        analysis: Option<Box<CompletedKktAnalysis>>,
    },
}

pub(crate) fn solve_equality_kkt(
    system: &EqualityKktSystem<'_>,
) -> Result<KktSolveEvidence, KktFailure> {
    let capacity = plan_equality_capacity(system.primal_variables, system.equality_constraints)
        .map_err(KktFailure::Capacity)?;
    validate_system(system)?;
    let attempt_plan = KktAttemptPlan::equality_v1();

    let primal_variables = system.primal_variables;
    let equality_constraints = system.equality_constraints;
    let dimension = capacity.kkt_dimension;
    let matrix_elements = dimension * dimension;
    let mut kkt = vec![0.0; matrix_elements];
    for row in 0..primal_variables {
        for column in 0..primal_variables {
            kkt[row + column * dimension] = system.hessian[row * primal_variables + column];
        }
    }
    for equality in 0..equality_constraints {
        for variable in 0..primal_variables {
            let value = system.equality_jacobian[equality * primal_variables + variable];
            kkt[primal_variables + equality + variable * dimension] = value;
            kkt[variable + (primal_variables + equality) * dimension] = value;
        }
    }

    let original_rhs = system
        .stationarity_rhs
        .iter()
        .chain(system.equality_rhs)
        .copied()
        .collect::<Vec<_>>();
    if let Some(index) = exact_zero_row(&kkt, dimension) {
        return Err(KktFailure::RankDeficient {
            evidence: exact_rank_deficiency(index, dimension),
        });
    }

    let scaled =
        equilibrate_symmetric_kkt(&kkt, &original_rhs, dimension).map_err(KktFailure::Scaling)?;
    let rank = analyze_rank(&scaled.matrix, dimension)?;
    let expected_inertia = Inertia {
        positive: primal_variables,
        negative: equality_constraints,
        zero: 0,
    };
    let observed_inertia = analyze_inertia(&scaled.matrix, dimension, &rank)?;
    if observed_inertia != expected_inertia {
        return Err(KktFailure::UnexpectedInertia {
            expected: expected_inertia,
            observed: observed_inertia,
            backend_invoked: false,
        });
    }
    let completed_analysis = CompletedKktAnalysis {
        rank: rank.clone(),
        expected_inertia,
        observed_inertia,
    };

    let scaling_summary = ScalingSummary {
        method: "block-aware Ruiz max-norm diagonal congruence",
        rounds: scaled.evidence.rounds.len(),
        saturated_outside_target: scaled.evidence.saturated_outside_target.len(),
    };
    let mut attempts = Vec::with_capacity(2);
    let (first_candidate, refinement_steps) =
        solve_lblt_with_refinement(&scaled.matrix, &scaled.rhs, dimension, &capacity)?;
    let first_residual =
        linear_residual_evidence(&scaled.matrix, dimension, &first_candidate, &scaled.rhs);
    let first_reason = candidate_rejection_reason(&first_candidate, first_residual);
    attempts.push(attempt_record(
        0,
        KktAttemptKind::BunchKaufmanRefinement,
        faer_backend::ALGORITHM,
        scaling_summary,
        refinement_steps,
        Some(first_residual),
        first_reason,
    ));

    let (scaled_solution, selected_backend, normalized_backward_error) = if first_reason.is_none() {
        (
            first_candidate,
            attempts[0].backend.clone(),
            first_residual.normalized_backward_error,
        )
    } else {
        let rescue = faer_backend::solve_with_full_svd(
            MatRef::from_column_major_slice(&scaled.matrix, dimension, dimension),
            &scaled.rhs,
        );
        match rescue {
            Ok(candidate) => {
                let residual =
                    linear_residual_evidence(&scaled.matrix, dimension, &candidate, &scaled.rhs);
                let reason = candidate_rejection_reason(&candidate, residual);
                attempts.push(attempt_record(
                    1,
                    KktAttemptKind::SvdRescue,
                    faer_backend::SVD_ALGORITHM,
                    scaling_summary,
                    0,
                    Some(residual),
                    reason,
                ));
                if let Some(reason) = reason {
                    return Err(exhausted_attempt_failure(
                        attempt_plan,
                        attempts,
                        reason,
                        Some(&completed_analysis),
                    ));
                }
                (
                    candidate,
                    attempts[1].backend.clone(),
                    residual.normalized_backward_error,
                )
            }
            Err(failure) => {
                if let faer_backend::DecompositionFailure::WorkspaceAllocation(failure) = failure {
                    return Err(KktFailure::WorkspaceAllocation {
                        phase: WorkspacePhase::SvdRescue,
                        bytes: failure.bytes,
                        alignment: failure.alignment,
                    });
                }
                let reason = NumericalFailureReason::BackendDecompositionFailure;
                attempts.push(attempt_record(
                    1,
                    KktAttemptKind::SvdRescue,
                    faer_backend::SVD_ALGORITHM,
                    scaling_summary,
                    0,
                    None,
                    Some(KktAttemptFailureReason::Numerical(reason)),
                ));
                return Err(exhausted_attempt_failure(
                    attempt_plan,
                    attempts,
                    KktAttemptFailureReason::Numerical(reason),
                    Some(&completed_analysis),
                ));
            }
        }
    };
    let solution = scaled.evidence.recover_solution(&scaled_solution);
    let scaling_round_trip_error = scaling_round_trip_error(
        &kkt,
        &original_rhs,
        &scaled.matrix,
        &scaled.rhs,
        &scaled.evidence,
    );
    if !scaling_round_trip_error.is_finite()
        || scaling_round_trip_error > EQUALITY_KKT_POLICY_V2.recovery_round_trip_limit
    {
        let reason = BackendContractViolationReason::ScalingRoundTripExceeded {
            observed: scaling_round_trip_error,
            limit: EQUALITY_KKT_POLICY_V2.recovery_round_trip_limit,
        };
        return Err(KktFailure::BackendContractViolation {
            attempt_plan,
            attempts,
            reason,
            analysis: Some(Box::new(completed_analysis)),
        });
    }

    Ok(KktSolveEvidence {
        candidate: solution[..primal_variables].to_vec(),
        equality_multipliers: solution[primal_variables..].to_vec(),
        normalized_backward_error,
        numerical_policy: EQUALITY_KKT_POLICY_V2.id,
        verification_scope: VerificationScope::BackendStandardForm,
        workspace: capacity.faer_workspace.clone(),
        backend: selected_backend,
        analysis_settings: algebraic_analysis_settings(),
        scaling: scaled.evidence,
        scaling_round_trip_error,
        rank,
        expected_inertia,
        observed_inertia,
        attempts,
        attempt_plan,
        capacity,
    })
}

fn exact_zero_row(matrix: &[f64], dimension: usize) -> Option<usize> {
    (0..dimension).find(|row| (0..dimension).all(|column| matrix[*row + column * dimension] == 0.0))
}

fn exact_rank_deficiency(index: usize, dimension: usize) -> AlgebraicRankEvidence {
    let (reject_ratio, accept_ratio) = EQUALITY_KKT_POLICY_V2.spectral_ratio_thresholds(dimension);
    AlgebraicRankEvidence {
        exact_zero_index: Some(index),
        rrqr_ratio: 0.0,
        singular_values: Vec::new(),
        svd_ratio: 0.0,
        reject_ratio,
        accept_ratio,
        classification: RankClassification::RankDeficient,
        backend_invoked: false,
    }
}

fn analyze_rank(matrix: &[f64], dimension: usize) -> Result<AlgebraicRankEvidence, KktFailure> {
    let matrix = MatRef::from_column_major_slice(matrix, dimension, dimension);
    let analysis = analyze_spectral_rank(matrix).map_err(|failure| match failure {
        SpectralAnalysisFailure::WorkspaceAllocation(failure) => KktFailure::WorkspaceAllocation {
            phase: WorkspacePhase::RankAnalysis,
            bytes: failure.bytes,
            alignment: failure.alignment,
        },
        SpectralAnalysisFailure::NumericalError => KktFailure::AlgebraicAnalysisFailure {
            phase: AlgebraicAnalysisPhase::RankConfirmation,
        },
    })?;
    let classification = match analysis.decision {
        SpectralRankDecision::Reject => RankClassification::RankDeficient,
        SpectralRankDecision::GrayZone => RankClassification::NumericalDecisionGrayZone,
        SpectralRankDecision::Accept => RankClassification::FullRank,
    };
    let evidence = AlgebraicRankEvidence {
        exact_zero_index: None,
        rrqr_ratio: analysis.rrqr_ratio,
        singular_values: analysis.singular_values,
        svd_ratio: analysis.svd_ratio,
        reject_ratio: analysis.reject_ratio,
        accept_ratio: analysis.accept_ratio,
        classification,
        backend_invoked: false,
    };
    match classification {
        RankClassification::FullRank => Ok(evidence),
        RankClassification::RankDeficient => Err(KktFailure::RankDeficient { evidence }),
        RankClassification::NumericalDecisionGrayZone => {
            Err(KktFailure::NumericalDecisionGrayZone { evidence })
        }
    }
}

fn analyze_inertia(
    matrix: &[f64],
    dimension: usize,
    rank: &AlgebraicRankEvidence,
) -> Result<Inertia, KktFailure> {
    let eigenvalues = faer_backend::self_adjoint_eigenvalues(MatRef::from_column_major_slice(
        matrix, dimension, dimension,
    ))
    .map_err(|failure| match failure {
        faer_backend::DecompositionFailure::WorkspaceAllocation(failure) => {
            KktFailure::WorkspaceAllocation {
                phase: WorkspacePhase::InertiaAnalysis,
                bytes: failure.bytes,
                alignment: failure.alignment,
            }
        }
        faer_backend::DecompositionFailure::NumericalError => {
            KktFailure::AlgebraicAnalysisFailure {
                phase: AlgebraicAnalysisPhase::Inertia,
            }
        }
    })?;
    let scale = eigenvalues
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let reject = rank.reject_ratio * scale;
    let accept = rank.accept_ratio * scale;
    let mut inertia = Inertia {
        positive: 0,
        negative: 0,
        zero: 0,
    };
    for value in eigenvalues {
        if value.abs() <= reject {
            inertia.zero += 1;
        } else if value.abs() < accept {
            return Err(KktFailure::NumericalDecisionGrayZone {
                evidence: AlgebraicRankEvidence {
                    classification: RankClassification::NumericalDecisionGrayZone,
                    ..rank.clone()
                },
            });
        } else if value > 0.0 {
            inertia.positive += 1;
        } else {
            inertia.negative += 1;
        }
    }
    Ok(inertia)
}

fn solve_lblt_with_refinement(
    matrix: &[f64],
    rhs: &[f64],
    dimension: usize,
    capacity: &EqualityCapacityPlan,
) -> Result<(Vec<f64>, usize), KktFailure> {
    let mut solution = rhs.to_vec();
    let mut factors = matrix.to_vec();
    let mut subdiagonal = vec![0.0; dimension];
    let mut permutation = vec![0usize; dimension];
    let mut inverse_permutation = vec![0usize; dimension];

    let factor_requirement = faer_backend::factor_workspace_requirement(dimension);
    let mut factor_memory =
        MemBuffer::try_new(factor_requirement).map_err(|_| KktFailure::WorkspaceAllocation {
            phase: WorkspacePhase::Factor,
            bytes: capacity.faer_workspace.factor.bytes,
            alignment: capacity.faer_workspace.factor.alignment,
        })?;
    let (_, permutation) = factor::cholesky_in_place(
        MatMut::from_column_major_slice_mut(&mut factors, dimension, dimension),
        DiagMut::from_slice_mut(&mut subdiagonal),
        &mut permutation,
        &mut inverse_permutation,
        faer_backend::parallelism(),
        MemStack::new(&mut factor_memory),
        faer_backend::lblt_params(),
    );
    drop(factor_memory);

    let solve_requirement = faer_backend::solve_workspace_requirement(dimension);
    let mut solve_memory =
        MemBuffer::try_new(solve_requirement).map_err(|_| KktFailure::WorkspaceAllocation {
            phase: WorkspacePhase::Solve,
            bytes: capacity.faer_workspace.solve.bytes,
            alignment: capacity.faer_workspace.solve.alignment,
        })?;
    solve::solve_in_place_with_conj(
        MatRef::from_column_major_slice(&factors, dimension, dimension),
        MatRef::from_column_major_slice(&factors, dimension, dimension).diagonal(),
        DiagRef::from_slice(&subdiagonal),
        Conj::No,
        permutation,
        MatMut::from_column_major_slice_mut(&mut solution, dimension, 1),
        faer_backend::parallelism(),
        MemStack::new(&mut solve_memory),
    );

    let mut refinement_steps = 0;
    while refinement_steps < EQUALITY_KKT_POLICY_V2.kkt_max_refinement_steps {
        let error =
            linear_residual_evidence(matrix, dimension, &solution, rhs).normalized_backward_error;
        if !error.is_finite()
            || error <= EQUALITY_KKT_POLICY_V2.backend_standard_form_backward_error_limit
            || solution.iter().any(|value| !value.is_finite())
        {
            break;
        }
        let mut correction = residual(matrix, dimension, &solution, rhs);
        solve::solve_in_place_with_conj(
            MatRef::from_column_major_slice(&factors, dimension, dimension),
            MatRef::from_column_major_slice(&factors, dimension, dimension).diagonal(),
            DiagRef::from_slice(&subdiagonal),
            Conj::No,
            permutation,
            MatMut::from_column_major_slice_mut(&mut correction, dimension, 1),
            faer_backend::parallelism(),
            MemStack::new(&mut solve_memory),
        );
        for (value, correction) in solution.iter_mut().zip(correction) {
            *value += correction;
        }
        refinement_steps += 1;
    }
    Ok((solution, refinement_steps))
}

fn residual(matrix: &[f64], dimension: usize, solution: &[f64], rhs: &[f64]) -> Vec<f64> {
    (0..dimension)
        .map(|row| {
            rhs[row]
                - (0..dimension)
                    .map(|column| matrix[row + column * dimension] * solution[column])
                    .sum::<f64>()
        })
        .collect()
}

fn candidate_rejection_reason(
    candidate: &[f64],
    residual: LinearResidualEvidence,
) -> Option<KktAttemptFailureReason> {
    if candidate.iter().any(|value| !value.is_finite()) {
        Some(KktAttemptFailureReason::BackendContract(
            BackendContractViolationReason::NonFiniteCandidate,
        ))
    } else if !residual.normalized_backward_error.is_finite()
        || residual.normalized_backward_error
            > EQUALITY_KKT_POLICY_V2.backend_standard_form_backward_error_limit
    {
        Some(KktAttemptFailureReason::BackendContract(
            BackendContractViolationReason::BackwardErrorExceeded {
                observed: residual.normalized_backward_error,
                limit: EQUALITY_KKT_POLICY_V2.backend_standard_form_backward_error_limit,
            },
        ))
    } else {
        None
    }
}

fn attempt_record(
    sequence: usize,
    kind: KktAttemptKind,
    algorithm: &'static str,
    scaling: ScalingSummary,
    refinement_steps: usize,
    residual: Option<LinearResidualEvidence>,
    failure_reason: Option<KktAttemptFailureReason>,
) -> KktAttemptRecord {
    KktAttemptRecord {
        sequence,
        kind,
        backend: backend_fingerprint(algorithm),
        settings: match kind {
            KktAttemptKind::BunchKaufmanRefinement => BackendAttemptSettings::Lblt {
                pivoting: faer_backend::PIVOTING,
                block_size: faer_backend::BLOCK_SIZE,
                parallelism_threshold: faer_backend::PARALLELISM_THRESHOLD,
                factor_workspace_source: faer_backend::FACTOR_WORKSPACE_SOURCE,
                maximum_refinement_steps: EQUALITY_KKT_POLICY_V2.kkt_max_refinement_steps,
            },
            KktAttemptKind::SvdRescue => BackendAttemptSettings::FullSvd {
                settings_id: faer_backend::SVD_SETTINGS_ID,
                left_vectors: "full",
                right_vectors: "full",
            },
        },
        scaling,
        requested_threads: faer_backend::REQUESTED_THREADS,
        actual_threads: faer_backend::parallelism().degree(),
        refinement_steps,
        termination: match failure_reason {
            Some(KktAttemptFailureReason::Numerical(_)) => SolveAttemptTermination::NumericalError,
            Some(KktAttemptFailureReason::BackendContract(_)) | None => {
                SolveAttemptTermination::CandidateProduced
            }
        },
        residual,
        certificate_present: false,
        failure_reason,
    }
}

fn exhausted_attempt_failure(
    attempt_plan: KktAttemptPlan,
    attempts: Vec<KktAttemptRecord>,
    reason: KktAttemptFailureReason,
    analysis: Option<&CompletedKktAnalysis>,
) -> KktFailure {
    match reason {
        KktAttemptFailureReason::BackendContract(reason) => KktFailure::BackendContractViolation {
            attempt_plan,
            attempts,
            reason,
            analysis: analysis.cloned().map(Box::new),
        },
        KktAttemptFailureReason::Numerical(reason) => KktFailure::NumericalFailure {
            attempt_plan,
            attempts,
            reason,
            analysis: analysis.cloned().map(Box::new),
        },
    }
}

fn scaling_round_trip_error(
    matrix: &[f64],
    rhs: &[f64],
    scaled_matrix: &[f64],
    scaled_rhs: &[f64],
    scaling: &RuizScalingEvidence,
) -> f64 {
    scaling
        .recover_matrix(scaled_matrix)
        .into_iter()
        .zip(matrix)
        .chain(scaling.recover_rhs(scaled_rhs).into_iter().zip(rhs))
        .map(|(actual, expected)| (actual - expected).abs() / expected.abs().max(1.0))
        .fold(0.0_f64, f64::max)
}

fn validate_system(system: &EqualityKktSystem<'_>) -> Result<(), KktFailure> {
    let expected_hessian = system
        .primal_variables
        .checked_mul(system.primal_variables)
        .unwrap_or(usize::MAX);
    let expected_jacobian = system
        .equality_constraints
        .checked_mul(system.primal_variables)
        .unwrap_or(usize::MAX);
    validate_slice(KktInputField::Hessian, system.hessian, expected_hessian)?;
    validate_slice(
        KktInputField::EqualityJacobian,
        system.equality_jacobian,
        expected_jacobian,
    )?;
    validate_slice(
        KktInputField::StationarityRightHandSide,
        system.stationarity_rhs,
        system.primal_variables,
    )?;
    validate_slice(
        KktInputField::EqualityRightHandSide,
        system.equality_rhs,
        system.equality_constraints,
    )
}

fn validate_slice(field: KktInputField, values: &[f64], expected: usize) -> Result<(), KktFailure> {
    if values.len() != expected {
        return Err(KktFailure::InvalidLength {
            field,
            expected,
            actual: values.len(),
        });
    }
    if let Some((index, _)) = values
        .iter()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(KktFailure::NonFiniteInput { field, index });
    }
    Ok(())
}

fn linear_residual_evidence(
    matrix: &[f64],
    dimension: usize,
    solution: &[f64],
    rhs: &[f64],
) -> LinearResidualEvidence {
    let residual_norm = (0..dimension)
        .map(|row| {
            let product = (0..dimension)
                .map(|column| matrix[row + column * dimension] * solution[column])
                .sum::<f64>();
            (product - rhs[row]).abs()
        })
        .fold(0.0_f64, f64::max);
    let matrix_norm = (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix[row + column * dimension].abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let solution_norm = solution
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let rhs_norm = rhs.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    let denominator = matrix_norm * solution_norm + rhs_norm;
    let normalized_backward_error = if denominator == 0.0 {
        residual_norm
    } else {
        residual_norm / denominator
    };
    LinearResidualEvidence {
        infinity_norm: residual_norm,
        matrix_infinity_norm: matrix_norm,
        solution_infinity_norm: solution_norm,
        rhs_infinity_norm: rhs_norm,
        normalized_backward_error,
    }
}

fn backend_fingerprint(algorithm: &'static str) -> BackendFingerprint {
    BackendFingerprint {
        schema_version: 1,
        crate_name: faer_backend::CRATE_NAME,
        crate_version: faer_backend::CRATE_VERSION,
        features: faer_backend::FEATURES,
        algorithm,
        target_arch: std::env::consts::ARCH,
        target_os: std::env::consts::OS,
        requested_threads: faer_backend::REQUESTED_THREADS,
        actual_threads: faer_backend::parallelism().degree(),
    }
}

fn algebraic_analysis_settings() -> AlgebraicAnalysisSettings {
    AlgebraicAnalysisSettings {
        llt_settings_id: faer_backend::LLT_SETTINGS_ID,
        householder_qr_settings_id: faer_backend::HOUSEHOLDER_QR_SETTINGS_ID,
        rrqr_algorithm: faer_backend::RRQR_ALGORITHM,
        rrqr_settings_id: faer_backend::RRQR_SETTINGS_ID,
        svd_algorithm: faer_backend::SVD_ALGORITHM,
        svd_settings_id: faer_backend::SVD_SETTINGS_ID,
        evd_algorithm: faer_backend::EVD_ALGORITHM,
        evd_settings_id: faer_backend::EVD_SETTINGS_ID,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::capacity::REPORT_FIXED_BYTES;

    #[test]
    fn symmetric_equality_kkt_returns_verified_candidate_and_backend_evidence() {
        let system = EqualityKktSystem {
            primal_variables: 2,
            equality_constraints: 1,
            hessian: &[2.0, 0.0, 0.0, 2.0],
            equality_jacobian: &[1.0, 1.0],
            stationarity_rhs: &[2.0, 0.0],
            equality_rhs: &[0.0],
        };

        let evidence = solve_equality_kkt(&system).expect("manufactured KKT should solve");

        assert_eq!(evidence.candidate, vec![0.5, -0.5]);
        assert_eq!(evidence.equality_multipliers, vec![1.0]);
        assert!(evidence.normalized_backward_error <= 1.0e-11);
        assert_eq!(evidence.capacity.kkt_dimension, 3);
        assert_eq!(evidence.workspace.factor.bytes, 64);
        assert_eq!(evidence.workspace.solve.bytes, 64);
        assert_eq!(evidence.backend.schema_version, 1);
        assert_eq!(evidence.backend.crate_name, "faer");
        assert_eq!(evidence.backend.crate_version, "0.24.4");
        assert_eq!(evidence.backend.features, ["linalg", "std"]);
        assert_eq!(evidence.backend.algorithm, "LBLT Bunch-Kaufman");
        assert_eq!(
            evidence.attempts[0].settings,
            BackendAttemptSettings::Lblt {
                pivoting: "PartialDiag",
                block_size: 64,
                parallelism_threshold: 16_384,
                factor_workspace_source: "faer::linalg::cholesky::lblt::factor::cholesky_in_place_scratch",
                maximum_refinement_steps: 2,
            }
        );
        assert_eq!(evidence.backend.requested_threads, 1);
        assert_eq!(evidence.backend.actual_threads, 1);
        assert_eq!(
            evidence.analysis_settings.rrqr_settings_id,
            "georbf-faer-rrqr-v1:block-size=256,blocking=2304,parallel=49152"
        );
        assert_eq!(
            evidence.analysis_settings.svd_settings_id,
            faer_backend::SVD_SETTINGS_ID
        );
        assert_eq!(
            evidence.analysis_settings.evd_settings_id,
            faer_backend::EVD_SETTINGS_ID
        );
        assert_eq!(evidence.numerical_policy.as_str(), "georbf-v2");
        assert_eq!(evidence.scaling.rounds.len(), 8);
        assert_eq!(evidence.rank.classification, RankClassification::FullRank);
        assert!(evidence.rank.rrqr_ratio > evidence.rank.accept_ratio);
        assert!(evidence.rank.svd_ratio > evidence.rank.accept_ratio);
        assert_eq!(
            evidence.expected_inertia,
            Inertia {
                positive: 2,
                negative: 1,
                zero: 0,
            }
        );
        assert_eq!(evidence.observed_inertia, evidence.expected_inertia);
        assert_eq!(evidence.attempts.len(), 1);
        assert_eq!(
            evidence.attempt_plan.attempts,
            [
                KktAttemptKind::BunchKaufmanRefinement,
                KktAttemptKind::SvdRescue,
            ]
        );
        assert!(
            evidence
                .attempt_plan
                .svd_rescue_requires_confirmed_full_rank
        );
        assert_eq!(
            evidence.attempts[0].kind,
            KktAttemptKind::BunchKaufmanRefinement
        );
        assert_eq!(
            evidence.attempts[0].termination,
            SolveAttemptTermination::CandidateProduced
        );
        let residual = evidence.attempts[0]
            .residual
            .expect("an accepted candidate records complete residual evidence");
        assert!(residual.infinity_norm.is_finite());
        assert!(residual.matrix_infinity_norm.is_finite());
        assert!(residual.solution_infinity_norm.is_finite());
        assert!(residual.rhs_infinity_norm.is_finite());
        assert_eq!(
            residual.normalized_backward_error,
            evidence.normalized_backward_error
        );
        assert!(evidence.scaling_round_trip_error <= 1.0e-11);
        assert_eq!(
            evidence.verification_scope,
            VerificationScope::BackendStandardForm
        );
    }

    #[test]
    fn exact_zero_structure_is_rejected_before_scaling_or_backend_invocation() {
        let failure = solve_equality_kkt(&EqualityKktSystem {
            primal_variables: 1,
            equality_constraints: 0,
            hessian: &[0.0],
            equality_jacobian: &[],
            stationarity_rhs: &[0.0],
            equality_rhs: &[],
        })
        .expect_err("an exact zero KKT row is structurally rank deficient");

        match failure {
            KktFailure::RankDeficient { evidence } => {
                assert_eq!(evidence.exact_zero_index, Some(0));
                assert_eq!(evidence.classification, RankClassification::RankDeficient);
                assert!(!evidence.backend_invoked);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn wrong_kkt_inertia_is_rejected_before_a_candidate_is_produced() {
        let failure = solve_equality_kkt(&EqualityKktSystem {
            primal_variables: 2,
            equality_constraints: 1,
            hessian: &[-2.0, 0.0, 0.0, -2.0],
            equality_jacobian: &[1.0, 1.0],
            stationarity_rhs: &[0.0, 0.0],
            equality_rhs: &[0.0],
        })
        .expect_err("a convex Equality KKT must have the policy's expected inertia");

        match failure {
            KktFailure::UnexpectedInertia {
                expected,
                observed,
                backend_invoked,
            } => {
                assert_eq!(expected.positive, 2);
                assert_eq!(expected.negative, 1);
                assert_ne!(observed, expected);
                assert!(!backend_invoked);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn nonfinite_candidates_exhaust_the_single_permitted_svd_rescue_and_fail_closed() {
        let failure = solve_equality_kkt(&EqualityKktSystem {
            primal_variables: 1,
            equality_constraints: 0,
            hessian: &[1.0e-320],
            equality_jacobian: &[],
            stationarity_rhs: &[1.0],
            equality_rhs: &[],
        })
        .expect_err("neither backend attempt may return a nonfinite candidate");

        match failure {
            KktFailure::BackendContractViolation {
                attempt_plan,
                attempts,
                reason,
                ..
            } => {
                assert_eq!(reason, BackendContractViolationReason::NonFiniteCandidate);
                assert_eq!(attempt_plan.numerical_policy.as_str(), "georbf-v2");
                assert_eq!(attempts.len(), 2);
                assert_eq!(attempts[0].kind, KktAttemptKind::BunchKaufmanRefinement);
                assert_eq!(attempts[1].kind, KktAttemptKind::SvdRescue);
                assert_eq!(
                    attempts[1].settings,
                    BackendAttemptSettings::FullSvd {
                        settings_id: faer_backend::SVD_SETTINGS_ID,
                        left_vectors: "full",
                        right_vectors: "full",
                    }
                );
                assert!(attempts[0].scaling.saturated_outside_target > 0);
                assert!(attempts.iter().all(|attempt| {
                    attempt.termination == SolveAttemptTermination::CandidateProduced
                        && attempt.residual.is_some()
                        && attempt.backend.crate_version == "0.24.4"
                        && attempt.requested_threads == 1
                        && attempt.actual_threads == 1
                }));
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn finite_backward_error_damage_exhausts_attempts_without_a_candidate() {
        let residual = LinearResidualEvidence {
            infinity_norm: 1.0e-6,
            matrix_infinity_norm: 1.0,
            solution_infinity_norm: 1.0,
            rhs_infinity_norm: 1.0,
            normalized_backward_error: 1.0e-6 / 3.0,
        };
        let reason = candidate_rejection_reason(&[1.0], residual)
            .expect("a finite candidate above the backward-error limit is rejected");
        let scaling = ScalingSummary {
            method: "block-aware Ruiz max-norm diagonal congruence",
            rounds: 8,
            saturated_outside_target: 0,
        };
        let attempts = vec![
            attempt_record(
                0,
                KktAttemptKind::BunchKaufmanRefinement,
                faer_backend::ALGORITHM,
                scaling,
                2,
                Some(residual),
                Some(reason),
            ),
            attempt_record(
                1,
                KktAttemptKind::SvdRescue,
                faer_backend::SVD_ALGORITHM,
                scaling,
                0,
                Some(residual),
                Some(reason),
            ),
        ];
        let failure =
            exhausted_attempt_failure(KktAttemptPlan::equality_v1(), attempts, reason, None);

        match failure {
            KktFailure::BackendContractViolation {
                attempts, reason, ..
            } => {
                assert_eq!(attempts.len(), 2);
                assert!(attempts.iter().all(|attempt| {
                    attempt.termination == SolveAttemptTermination::CandidateProduced
                        && attempt
                            .residual
                            .is_some_and(|evidence| evidence.normalized_backward_error.is_finite())
                }));
                assert_eq!(
                    reason,
                    BackendContractViolationReason::BackwardErrorExceeded {
                        observed: residual.normalized_backward_error,
                        limit: EQUALITY_KKT_POLICY_V2.backend_standard_form_backward_error_limit,
                    }
                );
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn kkt_rank_band_preserves_a_numerical_decision_gray_zone() {
        let coupling = 1.0 - 2.0e-13;
        let failure = solve_equality_kkt(&EqualityKktSystem {
            primal_variables: 2,
            equality_constraints: 0,
            hessian: &[1.0, coupling, coupling, 1.0],
            equality_jacobian: &[],
            stationarity_rhs: &[1.0, 1.0],
            equality_rhs: &[],
        })
        .expect_err("a full-rank guess is forbidden between the rank bands");

        match failure {
            KktFailure::NumericalDecisionGrayZone { evidence } => {
                assert_eq!(
                    evidence.classification,
                    RankClassification::NumericalDecisionGrayZone
                );
                assert!(evidence.svd_ratio > evidence.reject_ratio);
                assert!(evidence.svd_ratio < evidence.accept_ratio);
                assert!(!evidence.backend_invoked);
            }
            other => panic!("unexpected failure: {other:?}"),
        }
    }

    #[test]
    fn capacity_report_envelope_covers_the_fixed_kkt_evidence() {
        assert!(std::mem::size_of::<KktSolveEvidence>() as u64 <= REPORT_FIXED_BYTES);
    }

    #[test]
    fn ruiz_congruence_runs_eight_quantized_rounds_and_is_reversible() {
        let original = vec![2.0_f64.powi(-40), 0.0, 0.0, 2.0_f64.powi(20)];
        let rhs = vec![3.0, -5.0];

        let scaled = equilibrate_symmetric_kkt(&original, &rhs, 2)
            .expect("a finite KKT with no zero row is scalable");

        assert_eq!(scaled.evidence.rounds.len(), 8);
        assert!(scaled.evidence.rounds.iter().all(|round| {
            round
                .exponents
                .iter()
                .all(|exponent| (-8..=8).contains(exponent))
        }));
        assert!(
            scaled
                .evidence
                .cumulative_exponents
                .iter()
                .all(|exponent| (-32..=32).contains(exponent))
        );
        assert_eq!(scaled.matrix[1], scaled.matrix[2]);

        let recovered_matrix = scaled.evidence.recover_matrix(&scaled.matrix);
        let recovered_rhs = scaled.evidence.recover_rhs(&scaled.rhs);
        let physical_tolerances = vec![1.0e-10, 2.0e-8];
        let scaled_tolerances = scaled
            .evidence
            .scale_residual_or_tolerance(&physical_tolerances);
        assert_eq!(
            scaled
                .evidence
                .recover_residual_or_tolerance(&scaled_tolerances),
            physical_tolerances
        );
        for (actual, expected) in recovered_matrix.into_iter().zip(original) {
            assert_eq!(actual, expected);
        }
        for (actual, expected) in recovered_rhs.into_iter().zip(rhs) {
            assert_eq!(actual, expected);
        }
    }
}
