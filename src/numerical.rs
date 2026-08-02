use faer::MatRef;

use crate::faer_backend::{self, DecompositionFailure, WorkspaceAllocationFailure};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NumericalPolicyId(pub(crate) &'static str);

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EqualityKktNumericalPolicy {
    pub(crate) id: NumericalPolicyId,
    pub(crate) backend_standard_form_backward_error_limit: f64,
    pub(crate) spectral_reject_multiplier: f64,
    pub(crate) spectral_accept_multiplier: f64,
    pub(crate) reduced_symmetry_multiplier: f64,
    pub(crate) null_space_defect_limit: f64,
    pub(crate) affine_reproduction_limit: f64,
    pub(crate) side_condition_limit: f64,
    pub(crate) canonical_characteristic_tolerance_multiplier: f64,
    pub(crate) canonical_relation_reference_tolerance_multiplier: f64,
    pub(crate) recovery_round_trip_limit: f64,
    pub(crate) metric_determinant_one_multiplier: f64,
    pub(crate) ruiz_rounds: usize,
    pub(crate) ruiz_single_round_exponent_limit: i32,
    pub(crate) ruiz_cumulative_exponent_limit: i32,
    pub(crate) kkt_max_refinement_steps: usize,
}

pub(crate) const EQUALITY_KKT_POLICY_V1: EqualityKktNumericalPolicy = EqualityKktNumericalPolicy {
    id: NumericalPolicyId("georbf-v1"),
    backend_standard_form_backward_error_limit: 1.0e-11,
    spectral_reject_multiplier: 64.0,
    spectral_accept_multiplier: 4096.0,
    reduced_symmetry_multiplier: 256.0,
    null_space_defect_limit: 1.0e-12,
    affine_reproduction_limit: 1.0e-11,
    side_condition_limit: 1.0e-10,
    canonical_characteristic_tolerance_multiplier: 1.0e-10,
    canonical_relation_reference_tolerance_multiplier: 1.0e-8,
    recovery_round_trip_limit: 1.0e-11,
    metric_determinant_one_multiplier: 64.0,
    ruiz_rounds: 8,
    ruiz_single_round_exponent_limit: 8,
    ruiz_cumulative_exponent_limit: 32,
    kkt_max_refinement_steps: 2,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SpectralRankDecision {
    Reject,
    GrayZone,
    Accept,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SpectralRankEvidence {
    pub(crate) rrqr_ratio: f64,
    pub(crate) singular_values: Vec<f64>,
    pub(crate) svd_ratio: f64,
    pub(crate) reject_ratio: f64,
    pub(crate) accept_ratio: f64,
    pub(crate) rank: usize,
    pub(crate) decision: SpectralRankDecision,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SpectralAnalysisFailure {
    WorkspaceAllocation(WorkspaceAllocationFailure),
    NumericalError,
}

pub(crate) fn analyze_spectral_rank(
    matrix: MatRef<'_, f64>,
) -> Result<SpectralRankEvidence, SpectralAnalysisFailure> {
    let diagonal = faer_backend::rrqr_diagonal(matrix)
        .map_err(SpectralAnalysisFailure::WorkspaceAllocation)?;
    let rrqr_largest = diagonal.iter().copied().fold(0.0_f64, f64::max);
    let rrqr_smallest = diagonal.iter().copied().fold(f64::INFINITY, f64::min);
    let rrqr_ratio = if rrqr_largest > 0.0 {
        rrqr_smallest / rrqr_largest
    } else {
        0.0
    };
    let singular_values =
        faer_backend::singular_values(matrix).map_err(|failure| match failure {
            DecompositionFailure::WorkspaceAllocation(failure) => {
                SpectralAnalysisFailure::WorkspaceAllocation(failure)
            }
            DecompositionFailure::NumericalError => SpectralAnalysisFailure::NumericalError,
        })?;
    let largest = singular_values.first().copied().unwrap_or(0.0);
    let smallest = singular_values.last().copied().unwrap_or(0.0);
    let dimension = matrix.nrows().max(matrix.ncols());
    let (reject_ratio, accept_ratio) = EQUALITY_KKT_POLICY_V1.spectral_ratio_thresholds(dimension);
    let svd_ratio = if largest > 0.0 {
        smallest / largest
    } else {
        0.0
    };
    let rank = singular_values
        .iter()
        .filter(|singular_value| **singular_value > reject_ratio * largest)
        .count();
    Ok(SpectralRankEvidence {
        rrqr_ratio,
        singular_values,
        svd_ratio,
        reject_ratio,
        accept_ratio,
        rank,
        decision: EQUALITY_KKT_POLICY_V1.classify_spectral_ratio(dimension, svd_ratio),
    })
}

impl EqualityKktNumericalPolicy {
    pub(crate) fn spectral_thresholds(self, dimension: usize, scale: f64) -> (f64, f64) {
        let (reject_ratio, accept_ratio) = self.spectral_ratio_thresholds(dimension);
        (reject_ratio * scale, accept_ratio * scale)
    }

    pub(crate) fn spectral_ratio_thresholds(self, dimension: usize) -> (f64, f64) {
        let unit_scale = f64::EPSILON * dimension as f64;
        (
            self.spectral_reject_multiplier * unit_scale,
            self.spectral_accept_multiplier * unit_scale,
        )
    }

    pub(crate) fn classify_spectral_ratio(
        self,
        dimension: usize,
        ratio: f64,
    ) -> SpectralRankDecision {
        let (reject, accept) = self.spectral_ratio_thresholds(dimension);
        if ratio <= reject {
            SpectralRankDecision::Reject
        } else if ratio < accept {
            SpectralRankDecision::GrayZone
        } else {
            SpectralRankDecision::Accept
        }
    }
}
