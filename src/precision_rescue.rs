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

    pub(crate) fn to_f64(self) -> f64 {
        self.high + self.low
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
pub(crate) struct DoubleDoubleCubicJet {
    value: DoubleDouble,
    gradient_x: [DoubleDouble; 3],
    gradient_y: [DoubleDouble; 3],
    mixed_xy: [[DoubleDouble; 3]; 3],
}

impl DoubleDoubleCubicJet {
    pub(crate) fn value(self) -> DoubleDouble {
        self.value
    }

    pub(crate) fn gradient_x(self) -> [DoubleDouble; 3] {
        self.gradient_x
    }

    pub(crate) fn gradient_y(self) -> [DoubleDouble; 3] {
        self.gradient_y
    }

    pub(crate) fn mixed_xy(self) -> [[DoubleDouble; 3]; 3] {
        self.mixed_xy
    }
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

pub(crate) fn cubic_pairing_dd_certified(
    left: &CanonicalFunctional,
    right: &CanonicalFunctional,
    metric: &GlobalAnisotropyMetric,
) -> CertifiedDoubleDouble {
    let mut pairing = DoubleDouble::from(0.0);
    let mut absolute_scale = 0.0;
    let mut operations = 0;
    for left_term in left.terms() {
        for right_term in right.terms() {
            let jet = cubic_jet_dd(left_term.support(), right_term.support(), metric.matrix());
            let left_value = DoubleDouble::from(left_term.value_coefficient());
            let right_value = DoubleDouble::from(right_term.value_coefficient());
            let left_gradient = left_term.gradient_coefficient().map(DoubleDouble::from);
            let right_gradient = right_term.gradient_coefficient().map(DoubleDouble::from);
            let contributions = [
                left_value * right_value * jet.value,
                left_value * dot_dd(right_gradient, jet.gradient_y),
                right_value * dot_dd(left_gradient, jet.gradient_x),
                dot_dd(
                    left_gradient,
                    matrix_vector_dd(jet.mixed_xy, right_gradient),
                ),
            ];
            for contribution in contributions {
                pairing += contribution;
                absolute_scale += contribution.to_f64().abs();
                operations += 8;
            }
        }
    }
    CertifiedDoubleDouble::new(pairing, absolute_scale, operations)
}

pub(crate) fn cubic_jet_dd(
    x: [f64; 3],
    y: [f64; 3],
    metric: [[f64; 3]; 3],
) -> DoubleDoubleCubicJet {
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

pub(crate) const MAX_RESCUED_MODES: usize = 64;
pub(crate) const DOUBLE_DOUBLE_PRECISION_BITS: u32 = 106;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PrecisionRescueConclusion {
    Positive,
    AlgebraicZero,
    NegativeCurvature,
    GrayZone,
    CapacityExceeded,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SymmetricRescueResult {
    pub(crate) rescued_modes: usize,
    pub(crate) conclusion: PrecisionRescueConclusion,
    pub(crate) permutation: Vec<usize>,
    pub(crate) lower: Vec<DoubleDouble>,
    pub(crate) pivot_lower_bounds: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CertifiedDoubleDouble {
    pub(crate) value: DoubleDouble,
    pub(crate) error: f64,
}

impl CertifiedDoubleDouble {
    pub(crate) fn new(value: DoubleDouble, absolute_scale: f64, operations: usize) -> Self {
        Self {
            value,
            error: 8.0 * operations.max(1) as f64 * f64::EPSILON * f64::EPSILON * absolute_scale,
        }
    }

    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        let magnitude = value.to_f64().abs();
        Self {
            value,
            error: self.error / (2.0 * magnitude) + 8.0 * f64::EPSILON * f64::EPSILON * magnitude,
        }
    }

    fn div(self, right: Self) -> Self {
        let value = self.value / right.value;
        let denominator = right.value.to_f64().abs();
        let magnitude = value.to_f64().abs();
        Self {
            value,
            error: (self.error + magnitude * right.error) / denominator
                + 8.0 * f64::EPSILON * f64::EPSILON * magnitude,
        }
    }

    fn subtract_product(self, left: Self, right: Self) -> Self {
        let value = self.value - left.value * right.value;
        let product_scale = left.value.to_f64().abs() * right.value.to_f64().abs();
        Self {
            value,
            error: self.error
                + left.value.to_f64().abs() * right.error
                + right.value.to_f64().abs() * left.error
                + left.error * right.error
                + 8.0 * f64::EPSILON * f64::EPSILON * product_scale,
        }
    }
}

/// Classifies and factors a small symmetric Schur block with deterministic
/// diagonal pivoting. The matrix is row-major and is expected to have been
/// rebuilt from canonical inputs rather than promoted from rounded f64 Gram
/// entries.
pub(crate) fn classify_symmetric_schur(
    matrix: &[CertifiedDoubleDouble],
    dimension: usize,
    algebraic_zero: impl Fn(&[DoubleDouble]) -> bool,
) -> SymmetricRescueResult {
    assert_eq!(matrix.len(), dimension * dimension);
    if dimension > MAX_RESCUED_MODES {
        return SymmetricRescueResult {
            rescued_modes: dimension,
            conclusion: PrecisionRescueConclusion::CapacityExceeded,
            permutation: (0..dimension).collect(),
            lower: Vec::new(),
            pivot_lower_bounds: Vec::new(),
        };
    }

    let mut schur = matrix.to_vec();
    let mut permutation = (0..dimension).collect::<Vec<_>>();
    let mut lower = vec![
        CertifiedDoubleDouble {
            value: DoubleDouble::from(0.0),
            error: 0.0,
        };
        dimension * dimension
    ];
    let mut pivot_lower_bounds = Vec::with_capacity(dimension);
    for pivot in 0..dimension {
        let selected = (pivot..dimension)
            .max_by(|left, right| {
                schur[*left * dimension + *left]
                    .value
                    .to_f64()
                    .total_cmp(&schur[*right * dimension + *right].value.to_f64())
            })
            .expect("the remaining symmetric block is nonempty");
        if selected != pivot {
            for column in 0..dimension {
                schur.swap(pivot * dimension + column, selected * dimension + column);
            }
            for row in 0..dimension {
                schur.swap(row * dimension + pivot, row * dimension + selected);
            }
            for column in 0..pivot {
                lower.swap(pivot * dimension + column, selected * dimension + column);
            }
            permutation.swap(pivot, selected);
        }

        let diagonal = schur[pivot * dimension + pivot];
        let lower_bound = next_down(diagonal.value.to_f64() - diagonal.error);
        let upper_bound = next_up(diagonal.value.to_f64() + diagonal.error);
        if upper_bound < 0.0 {
            return SymmetricRescueResult {
                rescued_modes: dimension,
                conclusion: PrecisionRescueConclusion::NegativeCurvature,
                permutation,
                lower: lower.into_iter().map(|entry| entry.value).collect(),
                pivot_lower_bounds,
            };
        }
        if lower_bound <= 0.0 {
            let null_mode = schur_null_mode(&lower, &permutation, dimension, pivot);
            let conclusion = if diagonal.value.is_zero() && algebraic_zero(&null_mode) {
                PrecisionRescueConclusion::AlgebraicZero
            } else {
                PrecisionRescueConclusion::GrayZone
            };
            return SymmetricRescueResult {
                rescued_modes: dimension,
                conclusion,
                permutation,
                lower: lower.into_iter().map(|entry| entry.value).collect(),
                pivot_lower_bounds,
            };
        }

        pivot_lower_bounds.push(lower_bound);
        let root = diagonal.sqrt();
        lower[pivot * dimension + pivot] = root;
        for row in (pivot + 1)..dimension {
            let coordinate = schur[row * dimension + pivot].div(root);
            lower[row * dimension + pivot] = coordinate;
        }
        for row in (pivot + 1)..dimension {
            for column in row..dimension {
                let updated = schur[column * dimension + row].subtract_product(
                    lower[column * dimension + pivot],
                    lower[row * dimension + pivot],
                );
                schur[column * dimension + row] = updated;
                schur[row * dimension + column] = updated;
            }
        }
    }
    SymmetricRescueResult {
        rescued_modes: dimension,
        conclusion: PrecisionRescueConclusion::Positive,
        permutation,
        lower: lower.into_iter().map(|entry| entry.value).collect(),
        pivot_lower_bounds,
    }
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        value
    } else if value == 0.0 {
        f64::from_bits(1)
    } else if value > 0.0 {
        f64::from_bits(value.to_bits() + 1)
    } else {
        f64::from_bits(value.to_bits() - 1)
    }
}

fn next_down(value: f64) -> f64 {
    -next_up(-value)
}

fn schur_null_mode(
    lower: &[CertifiedDoubleDouble],
    permutation: &[usize],
    dimension: usize,
    pivot: usize,
) -> Vec<DoubleDouble> {
    let mut factor_mode = vec![DoubleDouble::from(0.0); dimension];
    factor_mode[pivot] = DoubleDouble::from(1.0);
    for row in (0..pivot).rev() {
        let tail = ((row + 1)..=pivot).fold(DoubleDouble::from(0.0), |sum, column| {
            sum + lower[column * dimension + row].value * factor_mode[column]
        });
        factor_mode[row] = -tail / lower[row * dimension + row].value;
    }
    let mut original_mode = vec![DoubleDouble::from(0.0); dimension];
    for (factor_index, original_index) in permutation.iter().copied().enumerate() {
        original_mode[original_index] = factor_mode[factor_index];
    }
    original_mode
}

#[cfg(test)]
mod tests {
    use crate::cubic::GlobalAnisotropyMetric;
    use crate::functional::{CanonicalFunctional, FunctionalDimension, FunctionalTerm};
    use crate::oracle_fixture::{hex_values, verify_artifact_identity};

    use super::{DoubleDouble, cubic_pairing_dd, symmetric_schur_entry};

    use super::{
        CertifiedDoubleDouble, MAX_RESCUED_MODES, PrecisionRescueConclusion,
        classify_symmetric_schur,
    };

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

    #[test]
    fn bounded_symmetric_rescue_matches_the_independent_schur_conclusions() {
        let one = DoubleDouble::from(1.0);
        let small_positive = DoubleDouble::from_components(1.0, 2.0_f64.powi(-100));
        let small_negative = DoubleDouble::from_components(1.0, -2.0_f64.powi(-100));
        let factors = [one];
        let positive = symmetric_schur_entry(small_positive, &factors, &factors);
        let zero = symmetric_schur_entry(one, &factors, &factors);
        let negative = symmetric_schur_entry(small_negative, &factors, &factors);

        assert_eq!(
            classify_symmetric_schur(&[CertifiedDoubleDouble::new(positive, 1.0, 1)], 1, |_| {
                false
            })
            .conclusion,
            PrecisionRescueConclusion::Positive,
        );
        assert_eq!(
            classify_symmetric_schur(&[CertifiedDoubleDouble::new(zero, 1.0, 1)], 1, |_| true)
                .conclusion,
            PrecisionRescueConclusion::AlgebraicZero,
        );
        assert_eq!(
            classify_symmetric_schur(&[CertifiedDoubleDouble::new(negative, 1.0, 1)], 1, |_| {
                false
            })
            .conclusion,
            PrecisionRescueConclusion::NegativeCurvature,
        );
    }

    #[test]
    fn bounded_symmetric_rescue_accepts_64_modes_and_never_truncates_65() {
        let diagonal = |dimension: usize| {
            (0..dimension)
                .flat_map(|row| {
                    (0..dimension).map(move |column| {
                        CertifiedDoubleDouble::new(
                            DoubleDouble::from(f64::from(row == column)),
                            f64::from(row == column),
                            1,
                        )
                    })
                })
                .collect::<Vec<_>>()
        };

        let accepted =
            classify_symmetric_schur(&diagonal(MAX_RESCUED_MODES), MAX_RESCUED_MODES, |_| false);
        assert_eq!(accepted.rescued_modes, MAX_RESCUED_MODES);
        assert_eq!(accepted.conclusion, PrecisionRescueConclusion::Positive);

        let rejected = classify_symmetric_schur(
            &diagonal(MAX_RESCUED_MODES + 1),
            MAX_RESCUED_MODES + 1,
            |_| false,
        );
        assert_eq!(rejected.rescued_modes, MAX_RESCUED_MODES + 1);
        assert_eq!(
            rejected.conclusion,
            PrecisionRescueConclusion::CapacityExceeded,
        );
    }

    #[test]
    fn certified_rescue_returns_gray_when_the_propagated_interval_spans_zero() {
        let residual = DoubleDouble::from(2.0_f64.powi(-100));
        let unresolved = CertifiedDoubleDouble::new(residual, 1.0, 3);

        assert_eq!(
            classify_symmetric_schur(&[unresolved], 1, |_| false).conclusion,
            PrecisionRescueConclusion::GrayZone,
        );
    }
}
