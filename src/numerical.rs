#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NumericalPolicyId(pub(crate) &'static str);

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EqualityKktNumericalPolicy {
    pub(crate) id: NumericalPolicyId,
    pub(crate) backend_standard_form_backward_error_limit: f64,
}

pub(crate) const EQUALITY_KKT_POLICY_V1: EqualityKktNumericalPolicy = EqualityKktNumericalPolicy {
    id: NumericalPolicyId("georbf-v1"),
    backend_standard_form_backward_error_limit: 1.0e-11,
};
