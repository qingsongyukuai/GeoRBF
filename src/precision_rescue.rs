//! Bounded precision primitives for canonical Cubic rescue calculations.
//!
//! This module deliberately owns its arithmetic instead of delegating to a
//! native multiprecision runtime. Its conformance seam is crate-private; issue
//! #42 can reuse the values without exposing matrices or precision choices in
//! GeoRBF's public API.

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::cubic::GlobalAnisotropyMetric;
use crate::functional::CanonicalFunctional;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct DoubleDouble {
    high: f64,
    low: f64,
}

impl DoubleDouble {
    pub(crate) fn from_components(high: f64, low: f64) -> Self {
        let (high, low) = two_sum(high, low);
        Self { high, low }
    }

    pub(crate) fn high(self) -> f64 {
        self.high
    }

    pub(crate) fn low(self) -> f64 {
        self.low
    }

    pub(crate) fn is_zero(self) -> bool {
        self.high == 0.0 && self.low == 0.0
    }

    pub(crate) fn is_positive(self) -> bool {
        self.high > 0.0 || (self.high == 0.0 && self.low > 0.0)
    }

    pub(crate) fn is_negative(self) -> bool {
        self.high < 0.0 || (self.high == 0.0 && self.low < 0.0)
    }

    pub(crate) fn sqrt(self) -> Self {
        if self.is_zero() {
            return Self::from(0.0);
        }
        if self.is_negative() {
            return Self {
                high: f64::NAN,
                low: f64::NAN,
            };
        }

        let mut root = Self::from(self.high.sqrt());
        // Two Newton corrections retain the second word even when the first
        // correction rounds into the leading word.
        for _ in 0..2 {
            root += (self / root - root) * Self::from(0.5);
        }
        root
    }
}

impl From<f64> for DoubleDouble {
    fn from(value: f64) -> Self {
        Self {
            high: value,
            low: 0.0,
        }
    }
}

impl Add for DoubleDouble {
    type Output = Self;

    fn add(self, right: Self) -> Self::Output {
        let (sum, high_error) = two_sum(self.high, right.high);
        let (low_sum, low_error) = two_sum(self.low, right.low);
        let (sum, carry) = two_sum(sum, high_error + low_sum);
        Self::from_components(sum, carry + low_error)
    }
}

impl AddAssign for DoubleDouble {
    fn add_assign(&mut self, right: Self) {
        *self = *self + right;
    }
}

impl Neg for DoubleDouble {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            high: -self.high,
            low: -self.low,
        }
    }
}

impl Sub for DoubleDouble {
    type Output = Self;

    fn sub(self, right: Self) -> Self::Output {
        self + -right
    }
}

impl SubAssign for DoubleDouble {
    fn sub_assign(&mut self, right: Self) {
        *self = *self - right;
    }
}

impl Mul for DoubleDouble {
    type Output = Self;

    fn mul(self, right: Self) -> Self::Output {
        let (product, error) = two_product(self.high, right.high);
        let cross = self.high * right.low + self.low * right.high;
        let (high, carry) = two_sum(product, error + cross);
        Self::from_components(high, carry + self.low * right.low)
    }
}

impl Div for DoubleDouble {
    type Output = Self;

    fn div(self, right: Self) -> Self::Output {
        let quotient_high = self.high / right.high;
        let remainder = self - right * Self::from(quotient_high);
        let quotient_middle = remainder.high / right.high;
        let remainder = remainder - right * Self::from(quotient_middle);
        let quotient_low = remainder.high / right.high;
        Self::from(quotient_high) + Self::from(quotient_middle) + Self::from(quotient_low)
    }
}

fn two_sum(left: f64, right: f64) -> (f64, f64) {
    let sum = left + right;
    let right_virtual = sum - left;
    let error = (left - (sum - right_virtual)) + (right - right_virtual);
    (sum, error)
}

fn two_product(left: f64, right: f64) -> (f64, f64) {
    let product = left * right;
    (product, left.mul_add(right, -product))
}

#[derive(Debug, Clone, Copy)]
struct DoubleDoubleCubicJet {
    value: DoubleDouble,
    gradient_x: [DoubleDouble; 3],
    gradient_y: [DoubleDouble; 3],
    mixed_xy: [[DoubleDouble; 3]; 3],
}

pub(crate) fn cubic_pairing_dd(
    left: &CanonicalFunctional,
    right: &CanonicalFunctional,
    metric: &GlobalAnisotropyMetric,
) -> DoubleDouble {
    let mut pairing = DoubleDouble::from(0.0);
    for left_term in left.terms() {
        for right_term in right.terms() {
            let jet = cubic_jet_dd(left_term.support(), right_term.support(), metric.matrix());
            let left_value = DoubleDouble::from(left_term.value_coefficient());
            let right_value = DoubleDouble::from(right_term.value_coefficient());
            let left_gradient = left_term.gradient_coefficient().map(DoubleDouble::from);
            let right_gradient = right_term.gradient_coefficient().map(DoubleDouble::from);
            pairing += left_value * right_value * jet.value;
            pairing += left_value * dot_dd(right_gradient, jet.gradient_y);
            pairing += right_value * dot_dd(left_gradient, jet.gradient_x);
            pairing += dot_dd(
                left_gradient,
                matrix_vector_dd(jet.mixed_xy, right_gradient),
            );
        }
    }
    pairing
}

fn cubic_jet_dd(x: [f64; 3], y: [f64; 3], metric: [[f64; 3]; 3]) -> DoubleDoubleCubicJet {
    let delta =
        std::array::from_fn(|axis| DoubleDouble::from(x[axis]) - DoubleDouble::from(y[axis]));
    if delta.into_iter().all(DoubleDouble::is_zero) {
        return DoubleDoubleCubicJet {
            value: DoubleDouble::from(0.0),
            gradient_x: [DoubleDouble::from(0.0); 3],
            gradient_y: [DoubleDouble::from(0.0); 3],
            mixed_xy: [[DoubleDouble::from(0.0); 3]; 3],
        };
    }

    let displacement_scale = delta
        .iter()
        .map(|component| component.high.abs().max(component.low.abs()))
        .fold(0.0_f64, f64::max);
    let scale = DoubleDouble::from(displacement_scale);
    let scaled_delta = delta.map(|component| component / scale);
    let metric_dd = metric.map(|row| row.map(DoubleDouble::from));
    let metric_delta = matrix_vector_dd(metric_dd, scaled_delta);
    let scaled_radius = dot_dd(scaled_delta, metric_delta).sqrt();
    let radius = scale * scaled_radius;
    let radius_cubed = radius * radius * radius;
    let gradient_scale = DoubleDouble::from(3.0) * scale * scale * scaled_radius;
    let gradient_x = metric_delta.map(|component| gradient_scale * component);
    let gradient_y = gradient_x.map(Neg::neg);
    let mixed_xy = std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            -DoubleDouble::from(3.0)
                * scale
                * (scaled_radius * metric_dd[row][column]
                    + metric_delta[row] * metric_delta[column] / scaled_radius)
        })
    });
    DoubleDoubleCubicJet {
        value: radius_cubed,
        gradient_x,
        gradient_y,
        mixed_xy,
    }
}

fn dot_dd(left: [DoubleDouble; 3], right: [DoubleDouble; 3]) -> DoubleDouble {
    left.into_iter()
        .zip(right)
        .fold(DoubleDouble::from(0.0), |sum, (left, right)| {
            sum + left * right
        })
}

fn matrix_vector_dd(
    matrix: [[DoubleDouble; 3]; 3],
    vector: [DoubleDouble; 3],
) -> [DoubleDouble; 3] {
    matrix.map(|row| dot_dd(row, vector))
}

pub(crate) fn symmetric_schur_entry(
    pairing: DoubleDouble,
    left_factors: &[DoubleDouble],
    right_factors: &[DoubleDouble],
) -> DoubleDouble {
    assert_eq!(left_factors.len(), right_factors.len());
    left_factors
        .iter()
        .copied()
        .zip(right_factors.iter().copied())
        .fold(pairing, |entry, (left, right)| entry - left * right)
}

#[cfg(test)]
mod tests {
    use crate::cubic::GlobalAnisotropyMetric;
    use crate::functional::{CanonicalFunctional, FunctionalDimension, FunctionalTerm};
    use crate::oracle_fixture::{hex_values, verify_artifact_identity};

    use super::{DoubleDouble, cubic_pairing_dd, symmetric_schur_entry};

    const CASES: &str =
        include_str!("../validation/oracle/precision-rescue-v1/cases/precision-rescue.json");
    const FIXTURE: &str =
        include_str!("../validation/oracle/precision-rescue-v1/fixtures/precision-rescue.json");
    const GENERATOR: &str = include_str!("../validation/oracle/precision-rescue-v1/generate.py");
    const MANIFEST: &str =
        include_str!("../validation/oracle/precision-rescue-v1/source-manifest.json");

    fn fixture_section(marker: &str) -> &'static str {
        let start = FIXTURE
            .find(marker)
            .unwrap_or_else(|| panic!("missing oracle fixture marker {marker}"));
        &FIXTURE[start..]
    }

    fn expected(marker: &str) -> (f64, f64) {
        let section = fixture_section(marker);
        (
            hex_values(section, "high", 1)[0],
            hex_values(section, "low", 1)[0],
        )
    }

    fn independent_two_sum(left: f64, right: f64) -> (f64, f64) {
        let sum = left + right;
        let right_virtual = sum - left;
        let error = (left - (sum - right_virtual)) + (right - right_virtual);
        (sum, error)
    }

    fn assert_oracle_dd(actual: DoubleDouble, expected: (f64, f64)) {
        let (high_delta, high_roundoff) = independent_two_sum(actual.high(), -expected.0);
        let (low_delta, low_roundoff) = independent_two_sum(actual.low(), -expected.1);
        let (combined, combined_roundoff) = independent_two_sum(high_delta, low_delta);
        let absolute_error =
            combined.abs() + combined_roundoff.abs() + high_roundoff.abs() + low_roundoff.abs();
        let scale = expected.0.abs().max(expected.1.abs());
        let limit = 8.0 * f64::EPSILON * f64::EPSILON * scale;
        assert!(
            absolute_error <= limit,
            "actual={actual:?}, expected={expected:?}, error={absolute_error:e}, limit={limit:e}"
        );
    }

    #[test]
    fn precision_rescue_corpus_has_a_fixed_independent_160_digit_identity() {
        verify_artifact_identity(
            CASES,
            "81275ae5686b9af66a7124999d8ff034a273087d9f07835bac046a09d2f0c2cd",
        );
        verify_artifact_identity(
            GENERATOR,
            "0bbf50406b7cecc12d8445b46e75ced13196606b3c2b75e201b2333fbf8b7c19",
        );
        verify_artifact_identity(
            FIXTURE,
            "d5c169e58455b09bd0b3fc78f1b968ca85f0406f1642dcf1c12b8801e9735c9e",
        );
        verify_artifact_identity(
            MANIFEST,
            "6db4bb77b50ee35fdfc1022b0c4555b89649f2c45c3fbde65bcb040fae7dd789",
        );
        assert!(MANIFEST.contains("\"working_decimal_digits\": 160"));
        assert!(MANIFEST.contains("\"packages\": []"));
        assert!(MANIFEST.contains("GeoRBF issue #41"));
    }

    #[test]
    fn double_double_basic_arithmetic_matches_the_independent_oracle() {
        let add = DoubleDouble::from_components(1.0e16, 1.0)
            + DoubleDouble::from_components(-1.0e16, 0.5);
        assert_oracle_dd(add, expected("\"operation\": \"add\""));

        let left = DoubleDouble::from_components(1.0, 2.0_f64.powi(-53));
        let right = DoubleDouble::from_components(1.0, -2.0_f64.powi(-54));
        assert_oracle_dd(left - right, expected("\"operation\": \"subtract\""));
        assert_oracle_dd(left * right, expected("\"operation\": \"multiply\""));

        let divisor = DoubleDouble::from_components(3.0, -2.0_f64.powi(-53));
        assert_oracle_dd(left / divisor, expected("\"operation\": \"divide\""));
        assert_oracle_dd(
            DoubleDouble::from_components(2.0, 2.0_f64.powi(-53)).sqrt(),
            expected("\"operation\": \"sqrt\""),
        );
    }

    #[test]
    fn cubic_canonical_pairing_matches_the_independent_oracle() {
        let left = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![
                FunctionalTerm::new([0.0; 3], 2.0, [1.0, -0.5, 0.25]),
                FunctionalTerm::new([0.25, -0.5, 1.25], -0.75, [-0.125, 0.375, 0.5]),
            ],
        )
        .expect("oracle left functional is finite");
        let right = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![
                FunctionalTerm::new([1.0, -2.0, 0.5], -1.5, [0.75, 1.0, -0.25]),
                FunctionalTerm::new([-0.75, 0.5, -1.0], 0.625, [0.5, -0.25, 0.125]),
            ],
        )
        .expect("oracle right functional is finite");
        let metric =
            GlobalAnisotropyMetric::new([[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]])
                .expect("oracle metric is determinant-one SPD");

        assert_oracle_dd(
            cubic_pairing_dd(&left, &right, &metric),
            expected("\"cubic_pairing\""),
        );
    }

    #[test]
    fn cubic_canonical_pairing_preserves_a_nonzero_close_support_mixed_derivative() {
        let left = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![FunctionalTerm::new(
                [1.0e-200, 0.0, 0.0],
                0.0,
                [1.0, 0.0, 0.0],
            )],
        )
        .expect("close-support left functional is finite");
        let right = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![FunctionalTerm::new([0.0; 3], 0.0, [1.0, 0.0, 0.0])],
        )
        .expect("close-support right functional is finite");

        let pairing = cubic_pairing_dd(&left, &right, &GlobalAnisotropyMetric::identity());

        assert!((pairing.high() + 6.0e-200).abs() <= 1.0e-214);
        assert!(pairing.low().is_finite());
    }

    #[test]
    fn symmetric_schur_oracle_distinguishes_small_positive_zero_and_negative() {
        let one = [DoubleDouble::from(1.0)];
        let small = DoubleDouble::from_components(1.0, 2.0_f64.powi(-100));
        let positive = symmetric_schur_entry(small, &one, &one);
        assert_oracle_dd(
            positive,
            expected("\"classification\": \"strictly_positive\""),
        );
        assert!(positive.is_positive());

        let zero = symmetric_schur_entry(DoubleDouble::from(1.0), &one, &one);
        assert_oracle_dd(zero, expected("\"classification\": \"algebraic_zero\""));
        assert!(zero.is_zero());

        let negative = symmetric_schur_entry(
            DoubleDouble::from_components(1.0, -2.0_f64.powi(-100)),
            &one,
            &one,
        );
        assert_oracle_dd(
            negative,
            expected("\"classification\": \"negative_curvature\""),
        );
        assert!(negative.is_negative());
    }

    #[test]
    fn symmetric_schur_accumulation_matches_the_general_oracle_in_both_orders() {
        let diagonal = DoubleDouble::from_components(2.0, 2.0_f64.powi(-53));
        let left = [
            DoubleDouble::from_components(1.0, 2.0_f64.powi(-53)),
            DoubleDouble::from_components(0.5, -2.0_f64.powi(-54)),
        ];
        let right = [
            DoubleDouble::from_components(1.0, -2.0_f64.powi(-54)),
            DoubleDouble::from_components(0.25, 2.0_f64.powi(-55)),
        ];
        let oracle = expected("\"classification\": \"general_symmetric_entry\"");

        assert_oracle_dd(symmetric_schur_entry(diagonal, &left, &right), oracle);
        assert_oracle_dd(symmetric_schur_entry(diagonal, &right, &left), oracle);
    }
}
