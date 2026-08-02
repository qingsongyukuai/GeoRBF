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
    pub(crate) field_value_recovery_limit: f64,
    pub(crate) field_derivative_recovery_limit: f64,
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
    field_value_recovery_limit: 1.0e-8,
    field_derivative_recovery_limit: 1.0e-8,
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
